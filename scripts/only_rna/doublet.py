from __future__ import annotations

import anndata as ad
import pandas as pd
import scanpy as sc

from .models import RunConfig


def _normalize_doublet_columns(adata: ad.AnnData) -> ad.AnnData:
    out = adata.copy()

    if "doublet_score" in out.obs:
        out.obs["doublet_score"] = (
            pd.to_numeric(out.obs["doublet_score"], errors="coerce")
            .fillna(0.0)
            .astype(float)
        )
    else:
        out.obs["doublet_score"] = 0.0

    if "is_doublet" in out.obs:
        out.obs["is_doublet"] = out.obs["is_doublet"].fillna(False).astype(bool)
    elif "predicted_doublet" in out.obs:
        out.obs["is_doublet"] = out.obs["predicted_doublet"].fillna(False).astype(bool)
    else:
        out.obs["is_doublet"] = False

    return out


def run_doublet_detection(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    del config

    out = adata.copy()
    if "doublet_score" in out.obs and "is_doublet" in out.obs:
        return _normalize_doublet_columns(out)

    try:
        sc.pp.scrublet(out)
    except Exception:
        return _normalize_doublet_columns(out)

    if "doublet_score" not in out.obs and "doublet_scores" in out.obs:
        out.obs["doublet_score"] = out.obs["doublet_scores"]

    return _normalize_doublet_columns(out)


__all__ = ["run_doublet_detection"]
