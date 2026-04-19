from __future__ import annotations

import importlib

import anndata as ad
import numpy as np
import pandas as pd

from scripts.only_rna.models import AzimuthConfig


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
