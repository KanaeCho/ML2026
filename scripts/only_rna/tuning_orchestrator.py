from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, cast

import anndata as ad
import pandas as pd

from .annotation import annotate_with_all_versions
from .config import merge_cli_overrides
from .discovery import ROOT, DiscoveredSample, resolve_data_root
from .doublet import run_doublet_detection
from .embedding import run_embedding
from .outputs import write_sample_outputs
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
    del candidate_id
    return Path(output_dir)


def _prune_stale_candidate_outputs(output_dir: Path, candidate_specs: list[CandidateSpec]) -> None:
    del candidate_specs
    tuning_dir = Path(output_dir) / "tuning"
    if tuning_dir.exists():
        shutil.rmtree(tuning_dir)


def _mean_confidence(series: pd.Series | None, pass_qc_mask: pd.Series) -> float:
    if series is None:
        return 0.0
    values = series.loc[pass_qc_mask]
    numeric = values.dropna().astype(float)
    if len(numeric) == 0:
        return 0.0
    return float(numeric.mean())


def _low_confidence_fraction(series: pd.Series | None, pass_qc_mask: pd.Series) -> float:
    if series is None:
        return 0.0
    values = series.loc[pass_qc_mask].fillna(False).astype(bool)
    if len(values) == 0:
        return 0.0
    return float(values.mean())


def _embedding_metrics(adata: ad.AnnData) -> tuple[float, float]:
    obs = cast(pd.DataFrame, adata.obs)
    pass_qc_mask = cast(pd.Series, obs["pass_qc"].fillna(False).astype(bool))
    pass_qc = int(pass_qc_mask.to_numpy(dtype=bool).sum())
    if pass_qc <= 0:
        return 0.0, 1.0

    clustered = (
        obs.loc[pass_qc_mask, "cluster"] if "cluster" in obs.columns else None
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
    ]


def _enumerate_candidates(config: RunConfig) -> list[CandidateSpec]:
    del config
    return _candidate_priority_order()


def evaluate_candidate(
    *,
    candidate_id: str,
    sample_id: str,
    gse: str,
    input_sample: DiscoveredSample,
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
        methods=["azimuth"],
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

    obs = cast(pd.DataFrame, adata.obs)
    pass_qc_mask = cast(pd.Series, obs["pass_qc"].fillna(False).astype(bool))
    n_cells_total = int(adata.n_obs)
    n_cells_pass_qc = int(pass_qc_mask.to_numpy(dtype=bool).sum())
    qc_score = score_qc_metrics(
        n_cells_total=n_cells_total,
        n_cells_pass_qc=n_cells_pass_qc,
    )

    annotation_status = dict(cast(Any, adata.uns.get("annotation_method_status", {}))).get(
        "azimuth", {}
    )
    annotation_score = score_annotation_metrics(
        method_status=str(annotation_status.get("status", "")),
        confidence_mean=_mean_confidence(
            cast(pd.Series | None, obs["azimuth_score"])
            if "azimuth_score" in obs.columns
            else None,
            pass_qc_mask,
        ),
        low_confidence_fraction=_low_confidence_fraction(
            cast(pd.Series | None, obs["azimuth_low_confidence"])
            if "azimuth_low_confidence" in obs.columns
            else None,
            pass_qc_mask,
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
    input_sample: DiscoveredSample,
    output_dir: Path,
    config: RunConfig,
) -> BoundedTuningResult:
    candidate_specs = _enumerate_candidates(config)
    _prune_stale_candidate_outputs(Path(output_dir), candidate_specs)
    if not candidate_specs:
        raise ValueError("No tuning candidates were generated")

    candidate = candidate_specs[0]
    try:
        best = evaluate_candidate(
            candidate_id=candidate.candidate_id,
            sample_id=sample_id,
            gse=gse,
            input_sample=input_sample,
            config=config,
            output_dir=output_dir,
        )
    except Exception as exc:
        raise ValueError(f"Baseline tuning candidate failed evaluation: {exc}") from exc

    return BoundedTuningResult(
        best_candidate_id=best.candidate_id,
        best_total_score=float(best.total_score),
        reason_code=best.reason_code,
        tuning_dir=Path(output_dir),
    )


__all__ = [
    "BoundedTuningResult",
    "CandidateEvaluation",
    "CandidateSpec",
    "evaluate_candidate",
    "run_bounded_tuning",
]
