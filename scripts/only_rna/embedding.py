from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .models import RunConfig


def _pass_qc_mask(adata: ad.AnnData) -> pd.Series:
    if "pass_qc" not in adata.obs:
        raise KeyError("adata.obs must contain 'pass_qc'")
    return adata.obs["pass_qc"].fillna(False).astype(bool)


def _initialize_embedding_columns(obs: pd.DataFrame) -> None:
    obs["cluster"] = pd.Series(pd.NA, index=obs.index, dtype="string")
    obs["umap_1"] = pd.Series(np.nan, index=obs.index, dtype=float)
    obs["umap_2"] = pd.Series(np.nan, index=obs.index, dtype=float)


def _run_scanpy_embedding(pass_qc_adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    del config

    out = pass_qc_adata.copy()
    if out.n_obs == 0:
        return out

    if out.n_obs == 1:
        out.obs["cluster"] = pd.Series(["0"], index=out.obs_names, dtype="string")
        out.obsm["X_umap"] = np.array([[0.0, 0.0]], dtype=float)
        return out

    sc.pp.normalize_total(out, target_sum=1e4)
    sc.pp.log1p(out)

    n_top_genes = max(1, min(2000, out.n_vars))
    sc.pp.highly_variable_genes(
        out,
        n_top_genes=n_top_genes,
        flavor="seurat",
        subset=True,
    )

    n_comps = max(1, min(20, out.n_obs - 1, out.n_vars))
    sc.pp.scale(out, max_value=10)
    sc.tl.pca(out, n_comps=n_comps)

    n_neighbors = max(1, min(15, out.n_obs - 1))
    sc.pp.neighbors(out, n_neighbors=n_neighbors, n_pcs=n_comps)

    if out.n_obs >= 3:
        try:
            sc.tl.leiden(out, key_added="cluster")
        except ImportError:
            out.obs["cluster"] = pd.Series(
                ["0"] * out.n_obs, index=out.obs_names, dtype="string"
            )
        try:
            sc.tl.umap(out)
        except ImportError:
            out.obsm["X_umap"] = np.column_stack(
                [np.arange(out.n_obs, dtype=float), np.zeros(out.n_obs, dtype=float)]
            )
    else:
        out.obs["cluster"] = pd.Series(
            ["0"] * out.n_obs, index=out.obs_names, dtype="string"
        )
        out.obsm["X_umap"] = np.column_stack(
            [np.arange(out.n_obs, dtype=float), np.zeros(out.n_obs, dtype=float)]
        )

    if "cluster" not in out.obs:
        out.obs["cluster"] = pd.Series(
            ["0"] * out.n_obs, index=out.obs_names, dtype="string"
        )
    else:
        out.obs["cluster"] = out.obs["cluster"].astype("string")

    if "X_umap" not in out.obsm:
        out.obsm["X_umap"] = np.column_stack(
            [np.arange(out.n_obs, dtype=float), np.zeros(out.n_obs, dtype=float)]
        )

    return out


def run_embedding(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    out = adata.copy()
    _initialize_embedding_columns(out.obs)

    pass_qc_mask = _pass_qc_mask(out)
    if not pass_qc_mask.any():
        return out

    embedded = _run_scanpy_embedding(out[pass_qc_mask].copy(), config)
    umap = np.asarray(embedded.obsm["X_umap"], dtype=float)

    out.obs.loc[embedded.obs_names, "cluster"] = embedded.obs["cluster"].astype(
        "string"
    )
    out.obs.loc[embedded.obs_names, "umap_1"] = umap[:, 0]
    out.obs.loc[embedded.obs_names, "umap_2"] = umap[:, 1]
    return out


__all__ = ["run_embedding"]
