from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from scripts.only_rna.discovery import DiscoveredSample


def test_default_tuning_presets_are_bounded_and_named() -> None:
    from scripts.only_rna.tuning_presets import default_tuning_presets

    presets = default_tuning_presets()

    assert set(presets.qc.keys()) == {"baseline", "strict", "lenient"}
    assert set(presets.azimuth.keys()) == {"baseline", "conservative", "smooth"}
    assert set(presets.embedding.keys()) == {"baseline", "separated", "stable"}


def test_candidate_score_includes_qc_annotation_embedding_and_reason_code() -> None:
    from scripts.only_rna.tuning_metrics import summarize_candidate_score

    score = summarize_candidate_score(
        qc_score=0.8,
        annotation_score=0.7,
        embedding_score=0.6,
        reason_code="balanced_default",
    )

    assert score.total_score == pytest.approx(0.7)
    assert score.reason_code == "balanced_default"


def test_score_qc_metrics_rewards_reasonable_retention() -> None:
    from scripts.only_rna.tuning_metrics import score_qc_metrics

    score = score_qc_metrics(n_cells_total=1000, n_cells_pass_qc=700)

    assert score == pytest.approx(0.7)


def test_score_annotation_metrics_penalizes_low_confidence_and_blocked_status() -> None:
    from scripts.only_rna.tuning_metrics import score_annotation_metrics

    ok_score = score_annotation_metrics(
        method_status="ok",
        confidence_mean=0.8,
        low_confidence_fraction=0.1,
    )
    blocked_score = score_annotation_metrics(
        method_status="error",
        confidence_mean=0.8,
        low_confidence_fraction=0.1,
    )

    assert ok_score == pytest.approx(0.72)
    assert blocked_score == 0.0


def test_score_embedding_metrics_penalizes_fragmentation() -> None:
    from scripts.only_rna.tuning_metrics import score_embedding_metrics

    compact_score = score_embedding_metrics(
        separation_score=0.8,
        fragmentation_penalty=0.1,
    )
    fragmented_score = score_embedding_metrics(
        separation_score=0.8,
        fragmentation_penalty=0.4,
    )

    assert compact_score == pytest.approx(0.7)
    assert fragmented_score == pytest.approx(0.4)


def test_run_tuning_selects_best_candidate_and_writes_selection_artifacts(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.tuning_orchestrator import (
        CandidateEvaluation,
        run_bounded_tuning,
    )

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.evaluate_candidate",
        lambda *args, **kwargs: CandidateEvaluation(
            candidate_id=kwargs["candidate_id"],
            total_score={
                "baseline__baseline__baseline": 0.6,
                "strict__baseline__stable": 0.9,
            }[kwargs["candidate_id"]],
            reason_code="test",
        ),
    )

    result = run_bounded_tuning(
        sample_id="GSM1",
        gse="GSE1",
        input_sample=object(),
        output_dir=tmp_path,
        config=config,
    )

    assert result.best_candidate_id == "strict__baseline__stable"
    assert (tmp_path / "tuning" / "selection_summary.json").exists()
    assert (tmp_path / "tuning" / "selected_params.json").exists()
    assert (tmp_path / "tuning" / "candidates.csv").exists()


def test_enumerate_candidates_respects_priority_and_max_candidates() -> None:
    from scripts.only_rna.config import load_run_config, merge_cli_overrides
    from scripts.only_rna.tuning_orchestrator import _enumerate_candidates

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    config = merge_cli_overrides(config, tuning__max_candidates=2)

    candidates = _enumerate_candidates(config)

    assert [candidate.candidate_id for candidate in candidates] == [
        "baseline__baseline__baseline",
        "strict__baseline__stable",
    ]


def test_run_tuning_raises_when_all_candidates_fail_but_keeps_candidate_audit(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config, merge_cli_overrides
    from scripts.only_rna.tuning_orchestrator import run_bounded_tuning

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    config = merge_cli_overrides(config, tuning__max_candidates=2)

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.evaluate_candidate",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("boom")),
    )

    with pytest.raises(ValueError, match="All tuning candidates failed evaluation"):
        run_bounded_tuning(
            sample_id="GSM1",
            gse="GSE1",
            input_sample=object(),
            output_dir=tmp_path,
            config=config,
        )

    assert (tmp_path / "tuning" / "candidates.csv").exists()
    assert not (tmp_path / "tuning" / "selection_summary.json").exists()
    assert not (tmp_path / "tuning" / "selected_params.json").exists()


def _make_sample(
    *, gse: str = "GSE123456", sample_id: str = "GSM123456"
) -> DiscoveredSample:
    return DiscoveredSample(
        gse=gse,
        sample_id=sample_id,
        input_type="triplet",
        sample_kind="gsm",
        supported=True,
        note="fixture",
        source_name=f"{sample_id}_matrix.mtx.gz",
        matrix_path=Path(f"/tmp/{sample_id}_matrix.mtx.gz"),
        barcodes_path=Path(f"/tmp/{sample_id}_barcodes.tsv.gz"),
        features_path=Path(f"/tmp/{sample_id}_features.tsv.gz"),
    )


def test_evaluate_candidate_runs_pipeline_with_candidate_config_and_writes_candidate_outputs(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.tuning_orchestrator import evaluate_candidate

    base_config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    sample = _make_sample(gse="GSE167363", sample_id="GSM5102900")

    stage_calls: list[tuple[str, object]] = []
    captured_output_roots: list[Path] = []

    def fake_read_sample_input(input_sample):
        stage_calls.append(("read", input_sample.sample_id))
        adata = ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 2.0]], dtype=float),
            obs=pd.DataFrame(index=["cell-1", "cell-2"]),
            var=pd.DataFrame(index=["GeneA", "GeneB"]),
        )
        adata.obs["gse"] = input_sample.gse
        adata.obs["sample_id"] = input_sample.sample_id
        adata.obs["input_type"] = input_sample.input_type
        return adata

    def fake_compute_qc_metrics(adata, config):
        stage_calls.append(("qc_metrics", config.qc.min_counts))
        out = adata.copy()
        out.obs["n_counts"] = [1000.0, 600.0]
        out.obs["n_genes"] = [400, 250]
        out.obs["pct_mt"] = [5.0, 8.0]
        out.obs["pct_ribo"] = [20.0, 18.0]
        return out

    def fake_run_doublet_detection(adata, config):
        stage_calls.append(("doublet", config.qc.max_pct_mt))
        out = adata.copy()
        out.obs["doublet_score"] = [0.01, 0.02]
        out.obs["is_doublet"] = [False, False]
        return out

    def fake_apply_qc_filters(adata, config):
        stage_calls.append(("qc_filters", config.qc.min_genes))
        out = adata.copy()
        out.obs["fails_count_floor"] = [False, False]
        out.obs["fails_gene_floor"] = [False, True]
        out.obs["fails_mt_ceiling"] = [False, False]
        out.obs["fails_ribo_ceiling"] = [False, False]
        out.obs["fails_doublet"] = [False, False]
        out.obs["pass_qc"] = [True, False]
        return out

    def fake_run_embedding(adata, config):
        stage_calls.append(("embedding", config.embedding.n_neighbors))
        out = adata.copy()
        out.obs["cluster"] = pd.Series(
            ["0", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["umap_1"] = [0.1, np.nan]
        out.obs["umap_2"] = [0.2, np.nan]
        return out

    def fake_annotate_with_all_versions(adata, reference_dir, **kwargs):
        stage_calls.append(("annotation", kwargs["azimuth_config"].k_weight))
        out = adata.copy()
        out.obs["cima_l1"] = pd.Series(
            ["CD4_T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l2"] = pd.Series(
            ["CD4_T_CM", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l1_masked"] = pd.Series(
            ["CD4_T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l1_score"] = [0.9, np.nan]
        out.obs["cima_l1_score_margin"] = [0.2, np.nan]
        out.obs["cima_l2_score"] = [0.85, np.nan]
        out.obs["cima_l2_score_margin"] = [0.15, np.nan]
        out.obs["cima_l1_low_confidence"] = [False, pd.NA]
        out.obs["azimuth_cell_type"] = pd.Series(
            ["CD4 T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["azimuth_score"] = [0.8, np.nan]
        out.obs["azimuth_score_margin"] = [0.1, np.nan]
        out.obs["azimuth_low_confidence"] = [False, pd.NA]
        out.uns["annotation_method_status"] = {
            "azimuth": {"status": "ok", "detail": kwargs["azimuth_config"].reference}
        }
        return out

    def fake_write_sample_outputs(adata, output_root, gse, sample_id, config):
        output_root = Path(output_root)
        captured_output_roots.append(output_root)
        stage_calls.append(("outputs", output_root))
        sample_dir = output_root / gse / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "metadata.csv").write_text("fixture\n", encoding="utf-8")
        return sample_dir

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.read_sample_input",
        fake_read_sample_input,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.compute_qc_metrics",
        fake_compute_qc_metrics,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.run_doublet_detection",
        fake_run_doublet_detection,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.apply_qc_filters",
        fake_apply_qc_filters,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.run_embedding",
        fake_run_embedding,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.annotate_with_all_versions",
        fake_annotate_with_all_versions,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.write_sample_outputs",
        fake_write_sample_outputs,
    )

    evaluation = evaluate_candidate(
        candidate_id="strict__conservative__stable",
        sample_id=sample.sample_id,
        gse=sample.gse,
        input_sample=sample,
        config=base_config,
        output_dir=tmp_path,
    )

    assert evaluation.candidate_id == "strict__conservative__stable"
    assert evaluation.reason_code != "not_evaluated_yet"
    assert evaluation.total_score > 0.0
    assert [name for name, _ in stage_calls] == [
        "read",
        "qc_metrics",
        "doublet",
        "qc_filters",
        "embedding",
        "annotation",
        "outputs",
    ]
    assert stage_calls[1] == ("qc_metrics", 800)
    assert stage_calls[3] == ("qc_filters", 500)
    assert stage_calls[4] == ("embedding", 20)
    assert stage_calls[5] == ("annotation", 30)

    expected_output_root = tmp_path / "tuning" / "strict__conservative__stable"
    assert captured_output_roots == [expected_output_root]
    assert (
        expected_output_root / sample.gse / sample.sample_id / "metadata.csv"
    ).exists()


def test_evaluate_candidate_uses_resolved_reference_dir(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.tuning_orchestrator import evaluate_candidate

    base_config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    sample = _make_sample(gse="GSE167363", sample_id="GSM5102900")
    expected_reference_dir = tmp_path / "external_data_root" / "reference"

    def fake_read_sample_input(input_sample):
        adata = ad.AnnData(
            X=np.array([[1.0, 0.0], [0.0, 2.0]], dtype=float),
            obs=pd.DataFrame(index=["cell-1", "cell-2"]),
            var=pd.DataFrame(index=["GeneA", "GeneB"]),
        )
        adata.obs["gse"] = input_sample.gse
        adata.obs["sample_id"] = input_sample.sample_id
        adata.obs["input_type"] = input_sample.input_type
        return adata

    def fake_compute_qc_metrics(adata, config):
        out = adata.copy()
        out.obs["n_counts"] = [1000.0, 600.0]
        out.obs["n_genes"] = [400, 250]
        out.obs["pct_mt"] = [5.0, 8.0]
        out.obs["pct_ribo"] = [20.0, 18.0]
        return out

    def fake_run_doublet_detection(adata, config):
        out = adata.copy()
        out.obs["doublet_score"] = [0.01, 0.02]
        out.obs["is_doublet"] = [False, False]
        return out

    def fake_apply_qc_filters(adata, config):
        out = adata.copy()
        out.obs["fails_count_floor"] = [False, False]
        out.obs["fails_gene_floor"] = [False, True]
        out.obs["fails_mt_ceiling"] = [False, False]
        out.obs["fails_ribo_ceiling"] = [False, False]
        out.obs["fails_doublet"] = [False, False]
        out.obs["pass_qc"] = [True, False]
        return out

    def fake_run_embedding(adata, config):
        out = adata.copy()
        out.obs["cluster"] = pd.Series(
            ["0", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["umap_1"] = [0.1, np.nan]
        out.obs["umap_2"] = [0.2, np.nan]
        return out

    captured_reference_dirs: list[Path] = []

    def fake_annotate_with_all_versions(adata, reference_dir, **kwargs):
        captured_reference_dirs.append(Path(reference_dir))
        out = adata.copy()
        out.obs["cima_l1"] = pd.Series(
            ["CD4_T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l2"] = pd.Series(
            ["CD4_T_CM", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l1_masked"] = pd.Series(
            ["CD4_T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["cima_l1_score"] = [0.9, np.nan]
        out.obs["cima_l1_score_margin"] = [0.2, np.nan]
        out.obs["cima_l2_score"] = [0.85, np.nan]
        out.obs["cima_l2_score_margin"] = [0.15, np.nan]
        out.obs["cima_l1_low_confidence"] = [False, pd.NA]
        out.obs["azimuth_cell_type"] = pd.Series(
            ["CD4 T", pd.NA], index=out.obs_names, dtype="string"
        )
        out.obs["azimuth_score"] = [0.8, np.nan]
        out.obs["azimuth_score_margin"] = [0.1, np.nan]
        out.obs["azimuth_low_confidence"] = [False, pd.NA]
        out.uns["annotation_method_status"] = {
            "azimuth": {"status": "ok", "detail": kwargs["azimuth_config"].reference}
        }
        return out

    def fake_write_sample_outputs(adata, output_root, gse, sample_id, config):
        sample_dir = Path(output_root) / gse / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "metadata.csv").write_text("fixture\n", encoding="utf-8")
        return sample_dir

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.resolve_data_root",
        lambda cwd=None: expected_reference_dir.parent,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.read_sample_input",
        fake_read_sample_input,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.compute_qc_metrics",
        fake_compute_qc_metrics,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.run_doublet_detection",
        fake_run_doublet_detection,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.apply_qc_filters",
        fake_apply_qc_filters,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.run_embedding",
        fake_run_embedding,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.annotate_with_all_versions",
        fake_annotate_with_all_versions,
    )
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.write_sample_outputs",
        fake_write_sample_outputs,
    )

    evaluate_candidate(
        candidate_id="baseline__baseline__baseline",
        sample_id=sample.sample_id,
        gse=sample.gse,
        input_sample=sample,
        config=base_config,
        output_dir=tmp_path,
    )

    assert captured_reference_dirs == [expected_reference_dir]
