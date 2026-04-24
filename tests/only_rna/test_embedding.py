from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd

from scripts.only_rna.embedding import run_embedding
from scripts.only_rna.models import (
    EmbeddingConfig,
    PlottingConfig,
    QcThresholds,
    RunConfig,
)


def _make_run_config(**embedding_overrides) -> RunConfig:
    return RunConfig(
        qc=QcThresholds(
            min_counts=500,
            min_genes=300,
            max_pct_mt=20.0,
            max_pct_ribo=60.0,
        ),
        plotting=PlottingConfig(
            umap_width=8.0,
            umap_height=6.0,
            dpi=150,
            point_size=3.0,
            legend_fontsize=10.0,
            legend_title_fontsize=11.0,
        ),
        embedding=EmbeddingConfig(**embedding_overrides),
    )


def test_run_embedding_uses_explicit_embedding_config(monkeypatch):
    adata = ad.AnnData(
        X=np.arange(1600, dtype=float).reshape(8, 200) + 1.0,
        obs=pd.DataFrame(
            {"pass_qc": [True] * 8}, index=[f"cell-{i}" for i in range(8)]
        ),
        var=pd.DataFrame(index=[f"Gene{i}" for i in range(200)]),
    )
    config = _make_run_config(
        n_top_genes=123,
        n_pcs=7,
        n_neighbors=6,
        resolution=2.5,
        min_dist=0.12,
        spread=1.7,
        random_state=11,
    )

    captured: dict[str, dict[str, object]] = {}

    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.pp.normalize_total", lambda *a, **k: None
    )
    monkeypatch.setattr("scripts.only_rna.embedding.sc.pp.log1p", lambda *a, **k: None)

    def fake_hvg(current, **kwargs):
        del current
        captured["hvg"] = kwargs

    def fake_pca(current, **kwargs):
        captured["pca"] = kwargs
        current.obsm["X_pca"] = np.ones((current.n_obs, kwargs["n_comps"]), dtype=float)

    def fake_neighbors(current, **kwargs):
        del current
        captured["neighbors"] = kwargs

    def fake_leiden(current, **kwargs):
        captured["leiden"] = kwargs
        current.obs["cluster"] = pd.Series(
            ["0"] * current.n_obs,
            index=current.obs_names,
            dtype="string",
        )

    def fake_umap(current, **kwargs):
        captured["umap"] = kwargs
        current.obsm["X_umap"] = np.column_stack(
            [
                np.arange(current.n_obs, dtype=float),
                np.zeros(current.n_obs, dtype=float),
            ]
        )

    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.pp.highly_variable_genes", fake_hvg
    )
    monkeypatch.setattr("scripts.only_rna.embedding.sc.pp.scale", lambda *a, **k: None)
    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.pca", fake_pca)
    monkeypatch.setattr("scripts.only_rna.embedding.sc.pp.neighbors", fake_neighbors)
    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.leiden", fake_leiden)
    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.umap", fake_umap)

    out = run_embedding(adata, config)

    assert captured["hvg"]["n_top_genes"] == 123
    assert captured["pca"]["n_comps"] == 7
    assert captured["neighbors"]["n_neighbors"] == 6
    assert captured["neighbors"]["n_pcs"] == 7
    assert captured["leiden"]["resolution"] == 2.5
    assert captured["umap"]["min_dist"] == 0.12
    assert captured["umap"]["spread"] == 1.7
    assert captured["umap"]["random_state"] == 11
    assert out.obs["cluster"].notna().all()
