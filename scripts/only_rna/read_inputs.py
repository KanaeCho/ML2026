from __future__ import annotations

import gzip
import tarfile
import tempfile
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.io import mmread
from scipy.sparse import csc_matrix

from .discovery import DiscoveredSample


def _open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


def _read_lines(path: Path) -> list[str]:
    with _open_maybe_gzip(path) as handle:
        return [line.rstrip("\n") for line in handle]


def _read_feature_table(path: Path) -> pd.DataFrame:
    with _open_maybe_gzip(path) as handle:
        return pd.read_csv(handle, sep="\t", header=None, dtype=str)


def _normalize_adata(adata: ad.AnnData, sample: DiscoveredSample) -> ad.AnnData:
    normalized = adata.copy()
    normalized.obs_names = normalized.obs_names.astype(str)
    normalized.var_names = normalized.var_names.astype(str)
    normalized.obs_names_make_unique()
    normalized.var_names_make_unique()
    normalized.obs["gse"] = sample.gse
    normalized.obs["sample_id"] = sample.sample_id
    normalized.obs["sample"] = sample.sample_id
    normalized.obs["dataset"] = sample.gse
    normalized.obs["input_type"] = sample.input_type
    normalized.obs["individual_id"] = sample.individual_id or ""
    normalized.obs["donor"] = sample.donor or sample.individual_id or ""
    normalized.obs["age"] = sample.age or ""
    normalized.obs["health"] = sample.health or ""
    return normalized


def _read_triplet(
    matrix_path: Path, barcodes_path: Path, features_path: Path
) -> ad.AnnData:
    with _open_maybe_gzip(matrix_path) as handle:
        counts = csc_matrix(mmread(handle))

    barcodes = _read_lines(barcodes_path)
    features = _read_feature_table(features_path)

    feature_ids = features.iloc[:, 0].astype(str)
    if features.shape[1] >= 2:
        feature_names = features.iloc[:, 1].astype(str)
    else:
        feature_names = feature_ids.copy()

    if features.shape[1] >= 3:
        feature_types = features.iloc[:, 2].astype(str)
        if (feature_types == "Gene Expression").any():
            keep = feature_types == "Gene Expression"
            counts = counts[keep.to_numpy(), :]
            feature_ids = feature_ids.loc[keep].reset_index(drop=True)
            feature_names = feature_names.loc[keep].reset_index(drop=True)

    adata = ad.AnnData(
        X=counts.transpose().tocsr(),
        obs=pd.DataFrame(index=pd.Index(barcodes, dtype=str)),
        var=pd.DataFrame(
            {
                "feature_id": feature_ids.to_list(),
                "feature_name": feature_names.to_list(),
            },
            index=pd.Index(feature_names.to_list(), dtype=str),
        ),
    )
    return adata


def _find_archived_triplet(extract_dir: Path) -> tuple[Path, Path, Path]:
    files = [path for path in extract_dir.rglob("*") if path.is_file()]

    def _pick(candidates: list[Path]) -> Path:
        if not candidates:
            raise FileNotFoundError(
                "Archive is missing matrix/barcodes/features triplet"
            )
        return sorted(candidates)[0]

    matrix_path = _pick(
        [path for path in files if path.name in {"matrix.mtx", "matrix.mtx.gz"}]
    )
    barcodes_path = _pick(
        [path for path in files if path.name in {"barcodes.tsv", "barcodes.tsv.gz"}]
    )
    features_path = _pick(
        [
            path
            for path in files
            if path.name
            in {"features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"}
        ]
    )
    return matrix_path, barcodes_path, features_path


def _read_archive(archive_path: Path) -> ad.AnnData:
    with tempfile.TemporaryDirectory(prefix="ml2026_only_rna_archive_") as tmp_dir:
        extract_dir = Path(tmp_dir)
        with tarfile.open(archive_path, "r:*") as handle:
            handle.extractall(extract_dir)
        matrix_path, barcodes_path, features_path = _find_archived_triplet(extract_dir)
        return _read_triplet(matrix_path, barcodes_path, features_path)


def _read_h5(h5_path: Path) -> ad.AnnData:
    adata = sc.read_10x_h5(h5_path, gex_only=False)
    if "feature_types" in adata.var.columns:
        feature_types = adata.var["feature_types"].astype(str)
        if (feature_types == "Gene Expression").any():
            adata = adata[:, feature_types == "Gene Expression"].copy()
    return adata


def read_sample_input(sample: DiscoveredSample) -> ad.AnnData:
    if sample.input_type == "triplet":
        if (
            not sample.matrix_path
            or not sample.barcodes_path
            or not sample.features_path
        ):
            raise ValueError("Triplet sample is missing required paths")
        adata = _read_triplet(
            sample.matrix_path, sample.barcodes_path, sample.features_path
        )
        return _normalize_adata(adata, sample)

    if sample.input_type == "archive":
        if not sample.archive_path:
            raise ValueError("Archive sample is missing archive_path")
        adata = _read_archive(sample.archive_path)
        return _normalize_adata(adata, sample)

    if sample.input_type == "h5":
        if not sample.h5_path:
            raise ValueError("H5 sample is missing h5_path")
        adata = _read_h5(sample.h5_path)
        return _normalize_adata(adata, sample)

    raise ValueError(f"Unsupported RNA input_type: {sample.input_type}")


__all__ = ["read_sample_input"]
