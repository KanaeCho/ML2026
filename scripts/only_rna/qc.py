from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .models import RunConfig


def _sum_axis(matrix, axis: int) -> np.ndarray:
    summed = matrix.sum(axis=axis)
    return np.asarray(summed).ravel().astype(float)


def _nnz_axis(matrix, axis: int) -> np.ndarray:
    if sparse.issparse(matrix):
        counts = matrix.getnnz(axis=axis)
    else:
        counts = np.count_nonzero(np.asarray(matrix), axis=axis)
    return np.asarray(counts).ravel().astype(int)


def _fractional_percent(counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        percent = np.divide(
            counts,
            totals,
            out=np.zeros_like(totals, dtype=float),
            where=totals > 0,
        )
    return percent * 100.0


def compute_qc_metrics(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    del config

    out = adata.copy()
    matrix = out.X
    var_names = pd.Index(out.var_names.astype(str))

    total_counts = _sum_axis(matrix, axis=1)
    n_genes = _nnz_axis(matrix, axis=1)

    mt_mask = np.asarray(var_names.str.startswith("MT-"), dtype=bool)
    ribo_mask = np.asarray(
        var_names.str.startswith("RPS") | var_names.str.startswith("RPL"),
        dtype=bool,
    )

    mt_counts = (
        _sum_axis(matrix[:, mt_mask], axis=1) if mt_mask.any() else np.zeros(out.n_obs)
    )
    ribo_counts = (
        _sum_axis(matrix[:, ribo_mask], axis=1)
        if ribo_mask.any()
        else np.zeros(out.n_obs)
    )

    out.obs["n_counts"] = total_counts
    out.obs["n_genes"] = n_genes
    out.obs["pct_mt"] = _fractional_percent(mt_counts, total_counts)
    out.obs["pct_ribo"] = _fractional_percent(ribo_counts, total_counts)
    return out


def apply_qc_filters(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    out = adata.copy()

    if "is_doublet" in out.obs:
        is_doublet = out.obs["is_doublet"].fillna(False).astype(bool)
    else:
        is_doublet = pd.Series(False, index=out.obs_names, dtype=bool)

    fails_count_floor = out.obs["n_counts"] < config.qc.min_counts
    fails_gene_floor = out.obs["n_genes"] < config.qc.min_genes
    fails_mt_ceiling = out.obs["pct_mt"] > config.qc.max_pct_mt
    fails_ribo_ceiling = out.obs["pct_ribo"] > config.qc.max_pct_ribo
    fails_doublet = is_doublet

    out.obs["fails_count_floor"] = fails_count_floor.astype(bool)
    out.obs["fails_gene_floor"] = fails_gene_floor.astype(bool)
    out.obs["fails_mt_ceiling"] = fails_mt_ceiling.astype(bool)
    out.obs["fails_ribo_ceiling"] = fails_ribo_ceiling.astype(bool)
    out.obs["fails_doublet"] = fails_doublet.astype(bool)

    pass_qc = ~(
        out.obs["fails_count_floor"]
        | out.obs["fails_gene_floor"]
        | out.obs["fails_mt_ceiling"]
        | out.obs["fails_ribo_ceiling"]
        | out.obs["fails_doublet"]
    )
    out.obs["pass_qc"] = pass_qc.astype(bool)
    return out


__all__ = ["compute_qc_metrics", "apply_qc_filters"]
