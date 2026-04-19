from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import subprocess
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .models import AzimuthConfig


@dataclass(frozen=True)
class AzimuthAnnotationResult:
    labels: pd.Series | None
    status: str
    detail: str


def _pass_qc_subset(adata: ad.AnnData) -> ad.AnnData:
    pass_qc = adata.obs.get("pass_qc", pd.Series(False, index=adata.obs.index)).fillna(
        False
    )
    pass_qc = pass_qc.astype(bool)
    return adata[pass_qc].copy()


def _run_azimuth_r(
    adata: ad.AnnData,
    *,
    reference: str = "pbmcref",
    annotation_level: str = "l1",
    max_cells: int | None = None,
    k_weight: int = 50,
    n_trees: int = 20,
    mapping_score_k: int = 100,
) -> pd.Series:
    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return pd.Series(dtype="string")

    if max_cells is not None:
        query = query[: min(max_cells, query.n_obs)].copy()
    else:
        query = query.copy()

    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        mtx_path = temp_path / "matrix.mtx"
        barcodes_path = temp_path / "barcodes.tsv"
        features_path = temp_path / "features.tsv"
        output_path = temp_path / "labels.csv"

        matrix = (
            query.X.T.tocsr()
            if sparse.issparse(query.X)
            else sparse.csr_matrix(np.asarray(query.X).T)
        )
        mmwrite(mtx_path, matrix)
        pd.Series(query.obs_names, dtype="string").to_csv(
            barcodes_path,
            sep="\t",
            index=False,
            header=False,
        )

        feature_ids = (
            query.var["feature_id"].astype("string")
            if "feature_id" in query.var.columns
            else pd.Series(query.var_names, index=query.var_names, dtype="string")
        )
        feature_names = (
            query.var["feature_name"].astype("string")
            if "feature_name" in query.var.columns
            else pd.Series(query.var_names, index=query.var_names, dtype="string")
        )
        pd.DataFrame({0: feature_ids.to_numpy(), 1: feature_names.to_numpy()}).to_csv(
            features_path,
            sep="\t",
            index=False,
            header=False,
        )

        r_code = f"""
        suppressPackageStartupMessages(library(Azimuth))
        suppressPackageStartupMessages(library(Seurat))
        suppressPackageStartupMessages(library(future))
        options(future.globals.maxSize = 2 * 1024^3)
        future::plan(\"sequential\")
        counts <- ReadMtx(
          mtx = '{mtx_path.as_posix()}',
          cells = '{barcodes_path.as_posix()}',
          features = '{features_path.as_posix()}',
          cell.column = 1,
          feature.column = 2
        )
        query <- CreateSeuratObject(counts = counts)
        mapped <- RunAzimuth(
          query = query,
          reference = '{reference}',
          verbose = FALSE,
          k.weight = {int(k_weight)},
          n.trees = {int(n_trees)},
          mapping.score.k = {int(mapping_score_k)}
        )
        md <- mapped[[]]
        pred_col <- paste0('predicted.celltype.{annotation_level}')
        out <- data.frame(cell_id = rownames(md), label = md[[pred_col]], row.names = NULL)
        write.csv(out, '{output_path.as_posix()}', row.names = FALSE)
        """
        subprocess.run(
            ["Rscript", "-e", r_code],
            check=True,
            capture_output=True,
            text=True,
        )
        labels = pd.read_csv(output_path)
        return pd.Series(
            labels["label"].astype("string").to_numpy(),
            index=labels["cell_id"].tolist(),
            dtype="string",
        )


def run_azimuth_annotation(
    adata: ad.AnnData,
    *,
    config: AzimuthConfig,
    annotation_level: str = "l1",
    max_cells: int | None = None,
) -> AzimuthAnnotationResult:
    if not config.enabled:
        return AzimuthAnnotationResult(
            labels=None,
            status="disabled",
            detail="Azimuth disabled in config",
        )

    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return AzimuthAnnotationResult(
            labels=None,
            status="no_pass_qc",
            detail="No pass_qc cells available for Azimuth",
        )

    try:
        labels = _run_azimuth_r(
            adata,
            reference=config.reference,
            annotation_level=annotation_level,
            max_cells=max_cells,
            k_weight=config.k_weight,
            n_trees=config.n_trees,
            mapping_score_k=config.mapping_score_k,
        )
    except Exception as exc:
        return AzimuthAnnotationResult(
            labels=None,
            status="error",
            detail=str(exc),
        )

    return AzimuthAnnotationResult(
        labels=labels,
        status="ok",
        detail=config.reference,
    )


__all__ = [
    "AzimuthAnnotationResult",
    "_pass_qc_subset",
    "_run_azimuth_r",
    "run_azimuth_annotation",
]
