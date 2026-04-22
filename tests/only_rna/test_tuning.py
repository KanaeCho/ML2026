from __future__ import annotations

from pathlib import Path
from typing import cast

import anndata as ad
import matplotlib.axes
import numpy as np
import pandas as pd
import pytest

from scripts.only_rna.discovery import DiscoveredSample


def _make_overview_candidate_adata(offset: float = 0.0) -> ad.AnnData:
    obs_index = pd.Index(["cell-1", "cell-2"], dtype=object)
    var_index = pd.Index(["GeneA", "GeneB"], dtype=object)
    return ad.AnnData(
        X=np.array([[1.0, 0.0], [0.0, 1.0]], dtype=float),
        obs=pd.DataFrame(
            {
                "umap_1": [0.1 + offset, 1.1 + offset],
                "umap_2": [1.0, 2.0],
                "azimuth_cell_type": pd.Series(
                    ["CD4 T", "B cell"],
                    index=obs_index,
                    dtype=object,
                ),
                "azimuth_cell_type_l1_raw": pd.Series(
                    ["CD4 T", "B"],
                    index=obs_index,
                    dtype=object,
                ),
                "azimuth_cell_type_l2_raw": pd.Series(
                    ["CD4 TCM", "B naive"],
                    index=obs_index,
                    dtype=object,
                ),
                "azimuth_cima_l1": pd.Series(
                    ["CD4_T", "B"],
                    index=obs_index,
                    dtype=object,
                ),
            },
            index=obs_index,
        ),
        var=pd.DataFrame(index=var_index),
    )


def test_default_tuning_presets_are_bounded_and_named() -> None:
    from scripts.only_rna.tuning_presets import default_tuning_presets

    presets = default_tuning_presets()

    assert set(presets.qc.keys()) == {"baseline"}
    assert set(presets.azimuth.keys()) == {"baseline"}
    assert set(presets.embedding.keys()) == {"baseline"}


def test_qc_calibration_profiles_match_stricter_mainline_defaults() -> None:
    from scripts.only_rna.qc_calibration import _profile_configs
    from scripts.only_rna.config import load_run_config

    profiles = _profile_configs()
    default_config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )

    baseline = profiles["baseline"].qc
    stricter = profiles["stricter_v1"].qc

    assert baseline.counts_lower_nmads == default_config.qc.counts_lower_nmads
    assert baseline.genes_lower_nmads == default_config.qc.genes_lower_nmads
    assert baseline.pct_mt_upper_nmads == default_config.qc.pct_mt_upper_nmads
    assert stricter.counts_lower_nmads == baseline.counts_lower_nmads
    assert stricter.genes_lower_nmads == baseline.genes_lower_nmads
    assert stricter.pct_mt_upper_nmads == baseline.pct_mt_upper_nmads
    assert stricter.pct_ribo_upper_nmads == baseline.pct_ribo_upper_nmads


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


def test_run_tuning_executes_single_baseline_candidate_and_writes_selection_artifacts(
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

    seen_candidate_ids: list[str] = []

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.evaluate_candidate",
        lambda *args, **kwargs: seen_candidate_ids.append(kwargs["candidate_id"])
        or CandidateEvaluation(
            candidate_id=kwargs["candidate_id"],
            total_score=0.6,
            reason_code="test",
        ),
    )

    result = run_bounded_tuning(
        sample_id="GSM1",
        gse="GSE1",
        input_sample=_make_sample(gse="GSE1", sample_id="GSM1"),
        output_dir=tmp_path,
        config=config,
    )

    assert seen_candidate_ids == ["baseline__baseline__baseline"]
    assert result.best_candidate_id == "baseline__baseline__baseline"
    assert result.tuning_dir == tmp_path


def test_write_tuning_selection_artifacts_emits_azimuth_candidate_overview_image(
    tmp_path: Path,
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [
        CandidateSpec("baseline", "baseline", "baseline"),
        CandidateSpec("lenient", "baseline", "separated"),
    ]
    evaluations = [
        CandidateEvaluation(
            candidate_id="baseline__baseline__baseline",
            total_score=0.71,
            reason_code="ok",
        ),
        CandidateEvaluation(
            candidate_id="lenient__baseline__separated",
            total_score=0.82,
            reason_code="ok",
        ),
    ]

    for candidate_id, offset in [
        ("baseline__baseline__baseline", 0.0),
        ("lenient__baseline__separated", 1.0),
    ]:
        candidate_dir = tmp_path / "tuning" / candidate_id / gse / sample_id
        candidate_dir.mkdir(parents=True, exist_ok=True)
        _make_overview_candidate_adata(offset=offset).write_h5ad(
            candidate_dir / f"{sample_id}.h5ad"
        )

    tuning_dir = write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="lenient__baseline__separated",
    )

    overview_cima_l1_path = tuning_dir / "umap_rna_candidates_overview_cima_l1.png"
    overview_pbmcref_path = tuning_dir / "umap_rna_candidates_overview_pbmcref.png"
    assert overview_cima_l1_path.exists()
    assert overview_pbmcref_path.exists()


def test_write_tuning_selection_artifacts_sorts_overview_candidates_by_parameter_id(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [
        CandidateSpec("strict", "smooth", "stable"),
        CandidateSpec("baseline", "baseline", "baseline"),
        CandidateSpec("lenient", "conservative", "separated"),
    ]
    evaluations = [
        CandidateEvaluation(
            candidate_id="strict__smooth__stable",
            total_score=0.41,
            reason_code="ok",
        ),
        CandidateEvaluation(
            candidate_id="baseline__baseline__baseline",
            total_score=0.81,
            reason_code="ok",
        ),
        CandidateEvaluation(
            candidate_id="lenient__conservative__separated",
            total_score=0.72,
            reason_code="ok",
        ),
    ]

    for candidate_id, offset in [
        ("strict__smooth__stable", 0.0),
        ("baseline__baseline__baseline", 1.0),
        ("lenient__conservative__separated", 2.0),
    ]:
        candidate_dir = tmp_path / "tuning" / candidate_id / gse / sample_id
        candidate_dir.mkdir(parents=True, exist_ok=True)
        _make_overview_candidate_adata(offset=offset).write_h5ad(
            candidate_dir / f"{sample_id}.h5ad"
        )

    captured: dict[str, list[tuple[str, float, ad.AnnData]]] = {}

    def fake_save_azimuth_candidate_overview(**kwargs):
        key = Path(kwargs["output_path"]).stem
        captured[key] = kwargs["candidates"]
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_azimuth_candidate_overview",
        fake_save_azimuth_candidate_overview,
    )

    write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="baseline__baseline__baseline",
    )

    expected = [
        "baseline__baseline__baseline",
        "lenient__conservative__separated",
        "strict__smooth__stable",
    ]
    assert [
        candidate_id
        for candidate_id, _, _ in cast(
            list[tuple[str, float, ad.AnnData]],
            captured["umap_rna_candidates_overview_cima_l1"],
        )
    ] == expected
    assert [
        candidate_id
        for candidate_id, _, _ in cast(
            list[tuple[str, float, ad.AnnData]],
            captured["umap_rna_candidates_overview_pbmcref"],
        )
    ] == expected


def test_save_azimuth_candidate_overview_uses_smaller_points_for_dense_panels(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.plotting import save_azimuth_candidate_overview

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    candidates = [
        ("baseline__baseline__baseline", 0.81, _make_overview_candidate_adata()),
        (
            "lenient__conservative__separated",
            0.72,
            _make_overview_candidate_adata(offset=1.0),
        ),
    ]

    captured_sizes: list[float] = []
    original_scatter = matplotlib.axes.Axes.scatter

    def capture_scatter(self, *args, **kwargs):
        size = kwargs.get("s")
        if size is not None:
            captured_sizes.append(float(size))
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

    output_path = tmp_path / "overview.png"
    save_azimuth_candidate_overview(
        candidates=candidates,
        output_path=output_path,
        title="Overview",
        config=config,
    )

    assert output_path.exists()
    assert captured_sizes
    assert max(captured_sizes) <= float(config.plotting.point_size) * 0.35


def test_save_azimuth_candidate_overview_uses_aligned_cima_l1_field(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.plotting import save_azimuth_candidate_overview

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    candidates = [
        ("baseline__baseline__baseline", 0.81, _make_overview_candidate_adata()),
    ]

    captured: dict[str, object] = {}

    def fake_plot_frame(adata, color_key):
        captured["color_key"] = color_key
        return pd.DataFrame(
            {
                "umap_1": [0.1, 1.1],
                "umap_2": [1.0, 2.0],
                color_key: ["CD4_T", "B"],
            }
        )

    def fake_palette(color_key, categories):
        captured["palette_key"] = color_key
        return {category: f"C{idx}" for idx, category in enumerate(categories)}

    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        captured["legend_title"] = kwargs.get("title")
        return original_legend(self, *args, **kwargs)

    monkeypatch.setattr("scripts.only_rna.plotting._plot_frame", fake_plot_frame)
    monkeypatch.setattr(
        "scripts.only_rna.plotting.get_hierarchical_palette", fake_palette
    )
    monkeypatch.setattr(matplotlib.axes.Axes, "legend", capture_legend)

    output_path = tmp_path / "overview.png"
    save_azimuth_candidate_overview(
        candidates=candidates,
        output_path=output_path,
        title="Overview",
        config=config,
    )

    assert output_path.exists()
    assert captured["color_key"] == "azimuth_cima_l1"
    assert captured["palette_key"] == "azimuth_cima_l1"
    assert captured["legend_title"] == "CIMA L1"


def test_write_tuning_selection_artifacts_uses_compact_point_size_for_overview(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [
        CandidateSpec("baseline", "baseline", "baseline"),
        CandidateSpec("lenient", "baseline", "separated"),
    ]
    evaluations = [
        CandidateEvaluation(
            candidate_id="baseline__baseline__baseline",
            total_score=0.71,
            reason_code="ok",
        ),
        CandidateEvaluation(
            candidate_id="lenient__baseline__separated",
            total_score=0.82,
            reason_code="ok",
        ),
    ]

    for candidate_id, offset in [
        ("baseline__baseline__baseline", 0.0),
        ("lenient__baseline__separated", 1.0),
    ]:
        candidate_dir = tmp_path / "tuning" / candidate_id / gse / sample_id
        candidate_dir.mkdir(parents=True, exist_ok=True)
        _make_overview_candidate_adata(offset=offset).write_h5ad(
            candidate_dir / f"{sample_id}.h5ad"
        )

    captured: dict[str, float] = {}

    def fake_save_azimuth_candidate_overview(**kwargs):
        captured[Path(kwargs["output_path"]).stem] = float(
            kwargs["config"].plotting.point_size
        )
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_azimuth_candidate_overview",
        fake_save_azimuth_candidate_overview,
    )

    write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="lenient__baseline__separated",
    )

    assert captured["umap_rna_candidates_overview_cima_l1"] <= 8.0
    assert captured["umap_rna_candidates_overview_pbmcref"] <= 8.0


def test_write_tuning_selection_artifacts_uses_three_by_three_square_layout_for_overviews(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [
        CandidateSpec("baseline", "baseline", "baseline"),
        CandidateSpec("baseline", "conservative", "baseline"),
        CandidateSpec("baseline", "smooth", "baseline"),
        CandidateSpec("lenient", "baseline", "separated"),
        CandidateSpec("lenient", "conservative", "separated"),
        CandidateSpec("lenient", "smooth", "separated"),
        CandidateSpec("strict", "baseline", "stable"),
        CandidateSpec("strict", "conservative", "stable"),
        CandidateSpec("strict", "smooth", "stable"),
    ]
    evaluations = [
        CandidateEvaluation(
            candidate_id=spec.candidate_id,
            total_score=0.9 - idx * 0.01,
            reason_code="ok",
        )
        for idx, spec in enumerate(candidate_specs)
    ]

    for idx, spec in enumerate(candidate_specs):
        candidate_dir = tmp_path / "tuning" / spec.candidate_id / gse / sample_id
        candidate_dir.mkdir(parents=True, exist_ok=True)
        _make_overview_candidate_adata(offset=float(idx)).write_h5ad(
            candidate_dir / f"{sample_id}.h5ad"
        )

    captured: dict[str, tuple[float, float]] = {}

    def fake_save_azimuth_candidate_overview(**kwargs):
        plotting = kwargs["config"].plotting
        captured[Path(kwargs["output_path"]).stem] = (
            float(plotting.umap_width),
            float(plotting.umap_height),
        )
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_azimuth_candidate_overview",
        fake_save_azimuth_candidate_overview,
    )

    write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="lenient__baseline__separated",
    )

    assert captured["umap_rna_candidates_overview_cima_l1"][0] == pytest.approx(
        captured["umap_rna_candidates_overview_cima_l1"][1]
    )
    assert captured["umap_rna_candidates_overview_pbmcref"][0] == pytest.approx(
        captured["umap_rna_candidates_overview_pbmcref"][1]
    )


def test_write_tuning_selection_artifacts_uses_fine_pbmcref_labels_for_raw_overview(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [
        CandidateSpec("baseline", "baseline", "baseline"),
    ]
    evaluations = [
        CandidateEvaluation(
            candidate_id="baseline__baseline__baseline",
            total_score=0.71,
            reason_code="ok",
        ),
    ]

    candidate_dir = (
        tmp_path / "tuning" / "baseline__baseline__baseline" / gse / sample_id
    )
    candidate_dir.mkdir(parents=True, exist_ok=True)
    _make_overview_candidate_adata(offset=0.0).write_h5ad(
        candidate_dir / f"{sample_id}.h5ad"
    )

    captured: dict[str, tuple[str | None, str | None]] = {}

    def fake_save_azimuth_candidate_overview(**kwargs):
        captured[Path(kwargs["output_path"]).stem] = (
            kwargs.get("color_key"),
            kwargs.get("legend_title"),
        )
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_azimuth_candidate_overview",
        fake_save_azimuth_candidate_overview,
    )

    write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="baseline__baseline__baseline",
    )

    assert captured["umap_rna_candidates_overview_pbmcref"] == (
        "azimuth_cell_type_l2_raw",
        "pbmcref",
    )


def test_save_azimuth_candidate_overview_adds_text_labels_to_pbmcref_panels(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.plotting import save_azimuth_candidate_overview

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    candidates = [
        ("baseline__baseline__baseline", 0.81, _make_overview_candidate_adata()),
    ]

    captured: list[str] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label in {"CD4 TCM", "B naive"}:
            captured.append(str(label))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "pbmcref_overview.png"
    save_azimuth_candidate_overview(
        candidates=candidates,
        output_path=output_path,
        title="pbmcref overview",
        config=config,
        color_key="azimuth_cell_type_l2_raw",
        legend_title="pbmcref",
    )

    assert output_path.exists()
    assert "CD4 TCM" in captured
    assert "B naive" in captured


def test_save_azimuth_candidate_overview_skips_text_labels_for_pbmcref_panels(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.plotting import save_azimuth_candidate_overview

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    candidates = [
        ("baseline__baseline__baseline", 0.81, _make_overview_candidate_adata()),
    ]

    captured: list[str] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label in {"CD4 TCM", "B naive"}:
            captured.append(str(label))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "pbmcref_overview.png"
    save_azimuth_candidate_overview(
        candidates=candidates,
        output_path=output_path,
        title="pbmcref overview",
        config=config,
        color_key="azimuth_cell_type_l2_raw",
        legend_title="pbmcref",
    )

    assert output_path.exists()
    assert captured == []


def test_write_tuning_selection_artifacts_emits_pbmcref_highlight_overview(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.outputs import write_tuning_selection_artifacts
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, CandidateSpec

    gse = "GSE157007"
    sample_id = "GSM4750304"
    candidate_specs = [CandidateSpec("baseline", "baseline", "baseline")]
    evaluations = [
        CandidateEvaluation(
            candidate_id="baseline__baseline__baseline",
            total_score=0.71,
            reason_code="ok",
        ),
    ]

    candidate_dir = (
        tmp_path / "tuning" / "baseline__baseline__baseline" / gse / sample_id
    )
    candidate_dir.mkdir(parents=True, exist_ok=True)
    _make_overview_candidate_adata(offset=0.0).write_h5ad(
        candidate_dir / f"{sample_id}.h5ad"
    )

    captured_outputs: list[str] = []

    def fake_save_azimuth_candidate_overview(**kwargs):
        captured_outputs.append(Path(kwargs["output_path"]).name)
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    def fake_save_highlight_category_overview(**kwargs):
        captured_outputs.append(Path(kwargs["output_path"]).name)
        kwargs["output_path"].write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_azimuth_candidate_overview",
        fake_save_azimuth_candidate_overview,
    )
    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_highlight_category_overview",
        fake_save_highlight_category_overview,
    )

    write_tuning_selection_artifacts(
        output_dir=tmp_path,
        sample_id=sample_id,
        gse=gse,
        candidate_specs=candidate_specs,
        evaluations=evaluations,
        best_candidate_id="baseline__baseline__baseline",
    )

    assert "umap_rna_candidates_overview_pbmcref.png" in captured_outputs
    assert "umap_rna_candidates_overview_pbmcref_highlight.png" in captured_outputs


def test_enumerate_candidates_always_returns_single_baseline_candidate() -> None:
    from scripts.only_rna.config import load_run_config, merge_cli_overrides
    from scripts.only_rna.tuning_orchestrator import _enumerate_candidates

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    config = merge_cli_overrides(config, tuning__max_candidates=2)

    candidates = _enumerate_candidates(config)

    assert [candidate.candidate_id for candidate in candidates] == [
        "baseline__baseline__baseline",
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

    with pytest.raises(ValueError, match="Baseline tuning candidate failed evaluation"):
        run_bounded_tuning(
            sample_id="GSM1",
            gse="GSE1",
            input_sample=_make_sample(gse="GSE1", sample_id="GSM1"),
            output_dir=tmp_path,
            config=config,
        )


def test_run_tuning_prunes_stale_non_baseline_candidate_directories(
    tmp_path: Path, monkeypatch
) -> None:
    from scripts.only_rna.config import load_run_config
    from scripts.only_rna.tuning_orchestrator import CandidateEvaluation, run_bounded_tuning

    config = load_run_config(
        Path("/mnt/f/ydd/ML2026/scripts/only_rna/default_config.yaml")
    )
    tuning_dir = tmp_path / "tuning"
    stale_dir = tuning_dir / "lenient__baseline__separated"
    keep_dir = tuning_dir / "baseline__baseline__baseline"
    stale_dir.mkdir(parents=True, exist_ok=True)
    keep_dir.mkdir(parents=True, exist_ok=True)

    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.evaluate_candidate",
        lambda *args, **kwargs: CandidateEvaluation(
            candidate_id=kwargs["candidate_id"],
            total_score=0.6,
            reason_code="ok",
        ),
    )

    run_bounded_tuning(
        sample_id="GSM1",
        gse="GSE1",
        input_sample=_make_sample(gse="GSE1", sample_id="GSM1"),
        output_dir=tmp_path,
        config=config,
    )

    assert not tuning_dir.exists()
    assert not stale_dir.exists()


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
        stage_calls.append(("qc_metrics", config.qc.count_floor_min))
        out = adata.copy()
        out.obs["n_counts"] = [1000.0, 600.0]
        out.obs["n_genes"] = [400, 250]
        out.obs["pct_mt"] = [5.0, 8.0]
        out.obs["pct_ribo"] = [20.0, 18.0]
        return out

    def fake_run_doublet_detection(adata, config):
        stage_calls.append(("doublet", config.qc.pct_mt_ceiling_max))
        out = adata.copy()
        out.obs["doublet_score"] = [0.01, 0.02]
        out.obs["is_doublet"] = [False, False]
        return out

    def fake_apply_qc_filters(adata, config):
        stage_calls.append(("qc_filters", config.qc.gene_floor_min))
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
        assert kwargs["methods"] == ["azimuth"]
        stage_calls.append(("annotation", kwargs["azimuth_config"].k_weight))
        out = adata.copy()
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
        candidate_id="baseline__baseline__baseline",
        sample_id=sample.sample_id,
        gse=sample.gse,
        input_sample=sample,
        config=base_config,
        output_dir=tmp_path,
    )

    assert evaluation.candidate_id == "baseline__baseline__baseline"
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
    assert stage_calls[1] == ("qc_metrics", 100)
    assert stage_calls[3] == ("qc_filters", 100)
    assert stage_calls[4] == ("embedding", 15)
    assert stage_calls[5] == ("annotation", 50)

    expected_output_root = tmp_path
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
        assert kwargs["methods"] == ["azimuth"]
        captured_reference_dirs.append(Path(reference_dir))
        out = adata.copy()
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
