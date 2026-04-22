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
    labels_by_level: dict[str, pd.Series]
    status: str
    detail: str
    scores: pd.Series | None = None
    score_margins: pd.Series | None = None
    low_confidence: pd.Series | None = None


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
) -> pd.DataFrame:
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

        # Large shared/GSE-level RNA queries can exceed the default future globals cap
        # inside RunAzimuth even in sequential mode.
        r_code = f"""
        suppressPackageStartupMessages(library(Azimuth))
        suppressPackageStartupMessages(library(Seurat))
        suppressPackageStartupMessages(library(future))
        options(future.globals.maxSize = 8 * 1024^3)
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
        out <- data.frame(
          cell_id = rownames(md),
          label = md[[paste0('predicted.celltype.{annotation_level}')]],
          label_l1 = if ('predicted.celltype.l1' %in% colnames(md)) md[['predicted.celltype.l1']] else NA,
          label_l2 = if ('predicted.celltype.l2' %in% colnames(md)) md[['predicted.celltype.l2']] else NA,
          label_l1_score = if ('predicted.celltype.l1.score' %in% colnames(md)) md[['predicted.celltype.l1.score']] else NA,
          label_l2_score = if ('predicted.celltype.l2.score' %in% colnames(md)) md[['predicted.celltype.l2.score']] else NA,
          mapping_score = if ('mapping.score' %in% colnames(md)) md[['mapping.score']] else NA,
          row.names = NULL
        )
        write.csv(out, '{output_path.as_posix()}', row.names = FALSE)
        """
        subprocess.run(
            ["Rscript", "-e", r_code],
            check=True,
            capture_output=True,
            text=True,
        )
        labels = pd.read_csv(output_path)
        return labels


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
            labels_by_level={},
            status="disabled",
            detail="Azimuth disabled in config",
        )

    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return AzimuthAnnotationResult(
            labels=None,
            labels_by_level={},
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
            labels_by_level={},
            status="error",
            detail=str(exc),
        )

    if isinstance(labels, pd.Series):
        label_series = labels.astype("string")
        labels_by_level = {annotation_level: label_series}
        score_series = None
        score_margin_series = None
        low_confidence_series = None
    else:
        label_series = pd.Series(
            labels["label"].astype("string").to_numpy(),
            index=labels["cell_id"].tolist(),
            dtype="string",
        )
        labels_by_level: dict[str, pd.Series] = {}
        for level_name, column in (("l1", "label_l1"), ("l2", "label_l2")):
            if column in labels.columns:
                labels_by_level[level_name] = pd.Series(
                    labels[column].astype("string").to_numpy(),
                    index=labels["cell_id"].tolist(),
                    dtype="string",
                )
        labels_by_level.setdefault(annotation_level, label_series)

        score_column = f"label_{annotation_level}_score"
        score_series = None
        if score_column in labels.columns:
            score_series = pd.Series(
                pd.to_numeric(labels[score_column], errors="coerce").to_numpy(),
                index=labels["cell_id"].tolist(),
                dtype=float,
            )
        elif "mapping_score" in labels.columns:
            score_series = pd.Series(
                pd.to_numeric(labels["mapping_score"], errors="coerce").to_numpy(),
                index=labels["cell_id"].tolist(),
                dtype=float,
            )

        if score_series is not None:
            # Current Azimuth export exposes one confidence score column, so reuse it
            # as a monotonic confidence margin surrogate until top2 probabilities are exported.
            score_margin_series = score_series.copy()
            low_confidence_series = pd.Series(
                score_series.fillna(0.0) < 0.5,
                index=score_series.index,
                dtype="boolean",
            )
        else:
            score_margin_series = None
            low_confidence_series = None

    return AzimuthAnnotationResult(
        labels=label_series,
        labels_by_level=labels_by_level,
        scores=score_series,
        score_margins=score_margin_series,
        low_confidence=low_confidence_series,
        status="ok",
        detail=config.reference,
    )


__all__ = [
    "AzimuthAnnotationResult",
    "_pass_qc_subset",
    "_run_azimuth_r",
    "run_azimuth_annotation",
]
