from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .annotation import annotate_with_all_versions
from .config import merge_cli_overrides
from .discovery import ROOT, resolve_data_root
from .doublet import run_doublet_detection
from .embedding import run_embedding
from .outputs import write_sample_outputs, write_tuning_selection_artifacts
from .qc import apply_qc_filters, compute_qc_metrics
from .read_inputs import read_sample_input
from .models import RunConfig
from .tuning_metrics import (
    score_annotation_metrics,
    score_embedding_metrics,
    score_qc_metrics,
    summarize_candidate_score,
)
from .tuning_presets import default_tuning_presets


@dataclass(frozen=True)
class CandidateSpec:
    qc_preset_id: str
    azimuth_preset_id: str
    embedding_preset_id: str

    @property
    def candidate_id(self) -> str:
        return (
            f"{self.qc_preset_id}__{self.azimuth_preset_id}__{self.embedding_preset_id}"
        )


@dataclass(frozen=True)
class CandidateEvaluation:
    candidate_id: str
    total_score: float
    reason_code: str


def _resolve_candidate_config(candidate_id: str, base_config: RunConfig) -> RunConfig:
    qc_preset_id, azimuth_preset_id, embedding_preset_id = candidate_id.split("__", 2)
    presets = default_tuning_presets()
    return merge_cli_overrides(
        base_config,
        qc=presets.qc[qc_preset_id],
        azimuth=presets.azimuth[azimuth_preset_id],
        embedding=presets.embedding[embedding_preset_id],
    )


def _candidate_output_root(output_dir: Path, candidate_id: str) -> Path:
    return Path(output_dir) / "tuning" / candidate_id


def _mean_confidence(series: object, pass_qc_mask) -> float:
    if series is None:
        return 0.0
    values = series.loc[pass_qc_mask]  # type: ignore[index]
    numeric = values.dropna().astype(float)
    if len(numeric) == 0:
        return 0.0
    return float(numeric.mean())


def _low_confidence_fraction(series: object, pass_qc_mask) -> float:
    if series is None:
        return 0.0
    values = series.loc[pass_qc_mask].fillna(False).astype(bool)  # type: ignore[index]
    if len(values) == 0:
        return 0.0
    return float(values.mean())


def _embedding_metrics(adata) -> tuple[float, float]:
    pass_qc_mask = adata.obs["pass_qc"].fillna(False).astype(bool)
    pass_qc = int(pass_qc_mask.sum())
    if pass_qc <= 0:
        return 0.0, 1.0

    clustered = (
        adata.obs.loc[pass_qc_mask, "cluster"] if "cluster" in adata.obs else None
    )
    if clustered is None:
        return 0.0, 1.0

    counts = clustered.dropna().astype(str).value_counts()
    if len(counts) == 0:
        return 0.0, 1.0

    separation_score = 1.0 if len(counts) > 1 else 0.5
    fragmentation_penalty = float((counts < 5).mean()) if len(counts) > 0 else 1.0
    return separation_score, fragmentation_penalty


@dataclass(frozen=True)
class BoundedTuningResult:
    best_candidate_id: str
    best_total_score: float
    reason_code: str
    tuning_dir: Path


def _candidate_priority_order() -> list[CandidateSpec]:
    return [
        CandidateSpec("baseline", "baseline", "baseline"),
        CandidateSpec("strict", "baseline", "stable"),
        CandidateSpec("lenient", "baseline", "separated"),
        CandidateSpec("baseline", "conservative", "baseline"),
        CandidateSpec("baseline", "smooth", "baseline"),
        CandidateSpec("strict", "conservative", "stable"),
        CandidateSpec("strict", "smooth", "stable"),
        CandidateSpec("lenient", "conservative", "separated"),
        CandidateSpec("lenient", "smooth", "separated"),
    ]


def _enumerate_candidates(config: RunConfig) -> list[CandidateSpec]:
    presets = default_tuning_presets()
    valid_qc = set(presets.qc)
    valid_azimuth = set(presets.azimuth)
    valid_embedding = set(presets.embedding)
    max_candidates = max(1, int(config.tuning.max_candidates))

    ordered = [
        candidate
        for candidate in _candidate_priority_order()
        if candidate.qc_preset_id in valid_qc
        and candidate.azimuth_preset_id in valid_azimuth
        and candidate.embedding_preset_id in valid_embedding
    ]
    return ordered[:max_candidates]


def evaluate_candidate(
    *,
    candidate_id: str,
    sample_id: str,
    gse: str,
    input_sample: object,
    config: RunConfig,
    output_dir: Path,
) -> CandidateEvaluation:
    candidate_config = _resolve_candidate_config(candidate_id, config)
    reference_dir = resolve_data_root(ROOT) / "reference"
    adata = read_sample_input(input_sample)
    adata = compute_qc_metrics(adata, candidate_config)
    adata = run_doublet_detection(adata, candidate_config)
    adata = apply_qc_filters(adata, candidate_config)
    adata = run_embedding(adata, candidate_config)
    adata = annotate_with_all_versions(
        adata,
        reference_dir=reference_dir,
        azimuth_model_dir=None,
        celltypist_model_path=None,
        singler_model_path=None,
        scanvi_model_path=None,
        methods=["cima", "azimuth"],
        azimuth_config=candidate_config.azimuth,
    )

    candidate_output_root = _candidate_output_root(output_dir, candidate_id)
    write_sample_outputs(
        adata,
        output_root=candidate_output_root,
        gse=gse,
        sample_id=sample_id,
        config=candidate_config,
    )

    pass_qc_mask = adata.obs["pass_qc"].fillna(False).astype(bool)
    n_cells_total = int(adata.n_obs)
    n_cells_pass_qc = int(pass_qc_mask.sum())
    qc_score = score_qc_metrics(
        n_cells_total=n_cells_total,
        n_cells_pass_qc=n_cells_pass_qc,
    )

    annotation_status = dict(adata.uns.get("annotation_method_status", {})).get(
        "azimuth", {}
    )
    annotation_score = score_annotation_metrics(
        method_status=str(annotation_status.get("status", "")),
        confidence_mean=_mean_confidence(adata.obs.get("azimuth_score"), pass_qc_mask),
        low_confidence_fraction=_low_confidence_fraction(
            adata.obs.get("azimuth_low_confidence"), pass_qc_mask
        ),
    )

    separation_score, fragmentation_penalty = _embedding_metrics(adata)
    embedding_score = score_embedding_metrics(
        separation_score=separation_score,
        fragmentation_penalty=fragmentation_penalty,
    )

    summary = summarize_candidate_score(
        qc_score=qc_score,
        annotation_score=annotation_score,
        embedding_score=embedding_score,
        reason_code=(
            str(annotation_status.get("status", "annotation_unavailable"))
            if annotation_status
            else "annotation_unavailable"
        ),
    )
    return CandidateEvaluation(
        candidate_id=candidate_id,
        total_score=float(summary.total_score),
        reason_code=summary.reason_code,
    )


def run_bounded_tuning(
    *,
    sample_id: str,
    gse: str,
    input_sample: object,
    output_dir: Path,
    config: RunConfig,
) -> BoundedTuningResult:
    candidate_specs = _enumerate_candidates(config)
    evaluations: list[CandidateEvaluation] = []

    for candidate in candidate_specs:
        try:
            evaluation = evaluate_candidate(
                candidate_id=candidate.candidate_id,
                sample_id=sample_id,
                gse=gse,
                input_sample=input_sample,
                config=config,
                output_dir=output_dir,
            )
        except Exception as exc:
            evaluation = CandidateEvaluation(
                candidate_id=candidate.candidate_id,
                total_score=0.0,
                reason_code=f"evaluation_error:{exc}",
            )
        evaluations.append(evaluation)

    if not evaluations:
        raise ValueError("No tuning candidates were generated")

    successful_evaluations = [
        evaluation
        for evaluation in evaluations
        if not str(evaluation.reason_code).startswith("evaluation_error:")
    ]

    if not successful_evaluations:
        write_tuning_selection_artifacts(
            output_dir=Path(output_dir),
            sample_id=sample_id,
            gse=gse,
            candidate_specs=candidate_specs,
            evaluations=evaluations,
            best_candidate_id="",
        )
        raise ValueError("All tuning candidates failed evaluation")

    best = max(successful_evaluations, key=lambda evaluation: evaluation.total_score)
    tuning_dir = write_tuning_selection_artifacts(
        output_dir=Path(output_dir),
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id=best.candidate_id,
    )
    return BoundedTuningResult(
        best_candidate_id=best.candidate_id,
        best_total_score=float(best.total_score),
        reason_code=best.reason_code,
        tuning_dir=tuning_dir,
    )


__all__ = [
    "BoundedTuningResult",
    "CandidateEvaluation",
    "CandidateSpec",
    "evaluate_candidate",
    "run_bounded_tuning",
]
