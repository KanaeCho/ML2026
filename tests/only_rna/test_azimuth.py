from __future__ import annotations

import importlib
from types import SimpleNamespace

import anndata as ad
import numpy as np
import pandas as pd

from scripts.only_rna.models import AzimuthConfig
from scripts.only_rna.annotation import annotate_with_all_versions
from scripts.only_rna.label_alignment import align_pbmcref_to_cima_l1


def _make_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [10.0, 0.0, 1.0],
                [0.0, 8.0, 2.0],
                [1.0, 1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"pass_qc": [True, True, False]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )


def test_run_azimuth_annotation_returns_real_labels_for_pass_qc_cells(monkeypatch):
    azimuth_module = importlib.import_module("scripts.only_rna.azimuth")

    fake = pd.Series(["CD4 T", "B"], index=["cell-1", "cell-2"], dtype="string")
    monkeypatch.setattr(
        azimuth_module,
        "_run_azimuth_r",
        lambda adata, reference="pbmcref", annotation_level="l1", max_cells=None, k_weight=50, n_trees=20, mapping_score_k=100: (
            fake
        ),
    )

    result = azimuth_module.run_azimuth_annotation(
        _make_adata(),
        config=AzimuthConfig(enabled=True),
        annotation_level="l1",
    )

    assert result.status == "ok"
    assert result.detail == "pbmcref"
    assert result.labels is not None
    assert list(result.labels.index) == ["cell-1", "cell-2"]
    assert result.labels.tolist() == ["CD4 T", "B"]


def test_run_azimuth_annotation_reports_disabled_without_running(monkeypatch):
    azimuth_module = importlib.import_module("scripts.only_rna.azimuth")

    def _should_not_run(**_kwargs):
        raise AssertionError(
            "_run_azimuth_r should not be called when Azimuth is disabled"
        )

    monkeypatch.setattr(azimuth_module, "_run_azimuth_r", _should_not_run)

    result = azimuth_module.run_azimuth_annotation(
        _make_adata(),
        config=AzimuthConfig(enabled=False),
        annotation_level="l1",
    )

    assert result.status == "disabled"
    assert "disabled" in result.detail.lower()
    assert result.labels is None


def test_run_azimuth_annotation_reports_errors_without_fallback_labels(monkeypatch):
    azimuth_module = importlib.import_module("scripts.only_rna.azimuth")

    def _raise(*args, **kwargs):
        del args, kwargs
        raise RuntimeError("Azimuth crashed")

    monkeypatch.setattr(azimuth_module, "_run_azimuth_r", _raise)

    result = azimuth_module.run_azimuth_annotation(
        _make_adata(),
        config=AzimuthConfig(enabled=True),
        annotation_level="l1",
    )

    assert result.status == "error"
    assert "Azimuth crashed" in result.detail
    assert result.labels is None


def test_annotate_with_all_versions_uses_first_configured_azimuth_annotation_level(
    tmp_path, monkeypatch
):
    adata = ad.AnnData(
        X=np.array(
            [
                [10.0, 0.0, 1.0],
                [0.0, 8.0, 2.0],
                [1.0, 1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "pass_qc": [True, False, True],
                "umap_1": [0.0, np.nan, 2.0],
                "umap_2": [1.0, np.nan, 3.0],
            },
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(
            {"feature_id": ["GeneA", "GeneB", "GeneC"]}, index=["g1", "g2", "g3"]
        ),
    )

    def _fake_run_azimuth_annotation(
        _adata,
        *,
        config: AzimuthConfig,
        annotation_level: str = "l1",
        max_cells: int | None = None,
    ) -> azimuth_module.AzimuthAnnotationResult:
        assert config.enabled is True
        assert annotation_level == "l2"
        assert max_cells is None
        return azimuth_module.AzimuthAnnotationResult(
            labels=pd.Series(
                ["Memory B", "Naive B"], index=["cell-1", "cell-3"], dtype="string"
            ),
            labels_by_level={
                "l2": pd.Series(
                    ["Memory B", "Naive B"],
                    index=["cell-1", "cell-3"],
                    dtype="string",
                )
            },
            status="ok",
            detail="pbmcref",
        )

    azimuth_module = importlib.import_module("scripts.only_rna.azimuth")
    monkeypatch.setattr(
        "scripts.only_rna.annotation.run_azimuth_annotation",
        _fake_run_azimuth_annotation,
    )

    out = annotate_with_all_versions(
        adata,
        reference_dir=tmp_path,
        methods=["azimuth"],
        azimuth_config=AzimuthConfig(enabled=True, annotation_levels=("l2",)),
    )

    assert out.obs.loc["cell-1", "azimuth_cell_type"] == "Memory B"
    assert out.obs.loc["cell-3", "azimuth_cell_type"] == "Naive B"


def test_annotate_with_all_versions_preserves_raw_azimuth_and_adds_cima_l1_aligned_fields(
    tmp_path, monkeypatch
):
    adata = ad.AnnData(
        X=np.array(
            [
                [10.0, 0.0, 1.0],
                [0.0, 8.0, 2.0],
                [1.0, 1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {
                "pass_qc": [True, False, True],
                "umap_1": [0.0, np.nan, 2.0],
                "umap_2": [1.0, np.nan, 3.0],
            },
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(
            {"feature_id": ["GeneA", "GeneB", "GeneC"]}, index=["g1", "g2", "g3"]
        ),
    )

    def _fake_run_azimuth_annotation(
        _adata,
        *,
        config: AzimuthConfig,
        annotation_level: str = "l1",
        max_cells: int | None = None,
    ):
        assert config.enabled is True
        assert annotation_level == "l1"
        assert max_cells is None
        raw_l1 = pd.Series(
            ["CD4 T", "B"],
            index=["cell-1", "cell-3"],
            dtype="string",
        )
        raw_l2 = pd.Series(
            ["Memory CD4 T", "Naive B"],
            index=["cell-1", "cell-3"],
            dtype="string",
        )
        return SimpleNamespace(
            labels=raw_l1,
            labels_by_level={"l1": raw_l1, "l2": raw_l2},
            status="ok",
            detail="pbmcref",
        )

    monkeypatch.setattr(
        "scripts.only_rna.annotation.run_azimuth_annotation",
        _fake_run_azimuth_annotation,
    )

    out = annotate_with_all_versions(
        adata,
        reference_dir=tmp_path,
        methods=["azimuth"],
        azimuth_config=AzimuthConfig(enabled=True, annotation_levels=("l1", "l2")),
    )

    assert out.obs.loc["cell-1", "azimuth_cell_type"] == "CD4 T"
    assert out.obs.loc["cell-1", "azimuth_cell_type_l1_raw"] == "CD4 T"
    assert out.obs.loc["cell-1", "azimuth_cell_type_l2_raw"] == "Memory CD4 T"
    assert out.obs.loc["cell-1", "azimuth_cima_l1"] == "CD4_T"
    assert bool(out.obs.loc["cell-1", "azimuth_cima_l1_unmapped"]) is False

    assert out.obs.loc["cell-3", "azimuth_cell_type"] == "B"
    assert out.obs.loc["cell-3", "azimuth_cell_type_l1_raw"] == "B"
    assert out.obs.loc["cell-3", "azimuth_cell_type_l2_raw"] == "Naive B"
    assert out.obs.loc["cell-3", "azimuth_cima_l1"] == "B"
    assert bool(out.obs.loc["cell-3", "azimuth_cima_l1_unmapped"]) is False

    assert pd.isna(out.obs.loc["cell-2", "azimuth_cell_type_l1_raw"])
    assert pd.isna(out.obs.loc["cell-2", "azimuth_cell_type_l2_raw"])
    assert pd.isna(out.obs.loc["cell-2", "azimuth_cima_l1"])


def test_align_pbmcref_to_cima_l1_prefers_specific_raw_l2_over_coarse_raw_l1():
    raw_l1 = pd.Series(
        ["CD8 T", "CD8 T", "CD8 T", "other"],
        index=["cell-1", "cell-2", "cell-3", "cell-4"],
        dtype="string",
    )
    raw_l2 = pd.Series(
        ["CD4 TCM", "ILC", "HSPC", "CD8 TEM"],
        index=["cell-1", "cell-2", "cell-3", "cell-4"],
        dtype="string",
    )

    aligned, unmapped = align_pbmcref_to_cima_l1(raw_l1, raw_l2)

    assert aligned.loc["cell-1"] == "CD4_T"
    assert bool(unmapped.loc["cell-1"]) is False
    assert aligned.loc["cell-2"] == "ILC"
    assert bool(unmapped.loc["cell-2"]) is False
    assert aligned.loc["cell-3"] == "Unknown"
    assert bool(unmapped.loc["cell-3"]) is True
    assert aligned.loc["cell-4"] == "CD8_T"
    assert bool(unmapped.loc["cell-4"]) is False
