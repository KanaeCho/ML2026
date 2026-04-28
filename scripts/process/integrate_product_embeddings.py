#!/usr/bin/env python3

from __future__ import annotations

import gzip
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import scipy.io
import scipy.sparse as sp
import harmonypy as hm
from sklearn.metrics import pairwise_distances


RNA_PRODUCTS = {"only_rna", "co_rna"}
ATAC_PRODUCTS = {"only_atac", "co_atac"}


@dataclass(frozen=True)
class IntegrationConfig:
    product: str
    product_dir: Path
    data_root: Path
    n_components: int = 30
    max_umap_fit_cells: int = 100_000
    n_clusters: int = 30
    random_state: int = 42
    batch_key: str = "sample_id"
    integration_method: str = "bbknn"
    neighbors_within_batch: int = 1
    bbknn_trim: int = 60
    leiden_resolution: float = 1.0
    rna_min_cima_l1_score: float = 0.0


def read_lines(path: Path) -> list[str]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n").split("\t")[0] for line in handle]


def read_feature_model(reference_dir: Path, product: str) -> pd.DataFrame:
    cima_dir = reference_dir / "cima"
    if product in RNA_PRODUCTS:
        return pd.read_csv(cima_dir / "cima_rna_reference_pca_features.tsv.gz", sep="\t")
    return pd.read_csv(cima_dir / "cima_atac_reference_lsi_features.tsv.gz", sep="\t")


def normalize_dense_rows(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return matrix / norms


def normalize_rows(matrix: sp.spmatrix) -> sp.csr_matrix:
    csr = sp.csr_matrix(matrix, dtype=np.float32, copy=True)
    row_sums = np.asarray(csr.sum(axis=1)).ravel().astype(np.float32)
    row_sums[row_sums == 0] = 1.0
    csr.data /= np.repeat(row_sums, np.diff(csr.indptr))
    return csr


def project_rna_sample(sample_dir: Path, sample_id: str, feature_model: pd.DataFrame, n_components: int) -> np.ndarray:
    h5ad_path = sample_dir / f"{sample_id}.h5ad"
    if not h5ad_path.exists():
        raise FileNotFoundError(h5ad_path)

    adata = ad.read_h5ad(h5ad_path)
    pass_mask = adata.obs["pass_qc"].fillna(False).astype(bool).to_numpy()
    pass_adata = adata[pass_mask]
    reference_features = feature_model["feature_id"].astype(str).tolist()
    dim_cols = [c for c in feature_model.columns if c.startswith("pc_dim_")][:n_components]
    loadings = feature_model[dim_cols].to_numpy(dtype=np.float32)
    gene_mean = feature_model["gene_mean"].to_numpy(dtype=np.float32)
    gene_std = feature_model["gene_std"].to_numpy(dtype=np.float32)
    gene_std[gene_std == 0] = 1.0

    candidates = []
    if "feature_name" in pass_adata.var.columns:
        candidates.append(pd.Index(pass_adata.var["feature_name"].astype(str)))
    candidates.append(pd.Index(pass_adata.var_names.astype(str)))
    if "feature_id" in pass_adata.var.columns:
        candidates.append(pd.Index(pass_adata.var["feature_id"].astype(str)))
    ref_set = set(reference_features)
    best = max(candidates, key=lambda idx: len(set(idx.astype(str)) & ref_set))
    lookup = pd.Series(np.arange(pass_adata.n_vars), index=best)
    lookup = lookup[~lookup.index.duplicated(keep="first")]

    ref_positions: list[int] = []
    query_positions: list[int] = []
    for ref_idx, feature in enumerate(reference_features):
        query_idx = lookup.get(feature)
        if query_idx is None or pd.isna(query_idx):
            continue
        ref_positions.append(ref_idx)
        query_positions.append(int(query_idx))

    query = np.zeros((pass_adata.n_obs, len(reference_features)), dtype=np.float32)
    if query_positions:
        selected = pass_adata.X[:, query_positions]
        if sp.issparse(pass_adata.X):
            cell_totals = np.asarray(pass_adata.X.sum(axis=1)).ravel().astype(np.float32)
        else:
            cell_totals = np.asarray(pass_adata.X, dtype=np.float32).sum(axis=1)
        cell_totals[cell_totals <= 0] = 1.0
        scale = (10000.0 / cell_totals).astype(np.float32)
        if sp.issparse(selected):
            normalized = sp.csr_matrix(selected, dtype=np.float32, copy=True)
            normalized.data *= np.repeat(scale, np.diff(normalized.indptr))
            normalized.data = np.log1p(normalized.data)
            selected_array = normalized.toarray().astype(np.float32, copy=False)
        else:
            selected_array = np.asarray(selected, dtype=np.float32) * scale[:, None]
            selected_array = np.log1p(selected_array).astype(np.float32, copy=False)
        query[:, ref_positions] = selected_array

    scaled = (query - gene_mean) / gene_std
    return np.asarray(scaled @ loadings, dtype=np.float32)


def assign_rna_cima_labels_from_embedding(
    metadata: pd.DataFrame,
    embedding: np.ndarray,
    reference_dir: Path,
    n_components: int,
) -> pd.DataFrame:
    cima_dir = reference_dir / "cima"
    l1 = pd.read_csv(cima_dir / "cima_rna_reference_l1_centroids.tsv", sep="\t")
    l2 = pd.read_csv(cima_dir / "cima_rna_reference_l2_centroids.tsv", sep="\t")
    hierarchy = pd.read_csv(cima_dir / "cima_rna_celltype_hierarchy.csv")
    dims = [c for c in l1.columns if c.startswith("pc_dim_")][:n_components]
    query = normalize_dense_rows(embedding[:, : len(dims)])

    l1_matrix = normalize_dense_rows(l1[dims].to_numpy(dtype=np.float32))
    l1_scores = query @ l1_matrix.T
    l1_idx = np.argmax(l1_scores, axis=1)
    l1_labels = l1["cell_type_l1"].astype(str).to_numpy()[l1_idx]
    metadata["integrated_cima_l1"] = l1_labels
    metadata["integrated_cima_l1_score"] = l1_scores[np.arange(len(query)), l1_idx]

    l2_matrix_all = normalize_dense_rows(l2[dims].to_numpy(dtype=np.float32))
    l2_labels_all = l2["cell_type_l2"].astype(str).to_numpy()
    l2_by_l1 = {
        str(k): set(v["cell_type_l2"].astype(str).tolist())
        for k, v in hierarchy.groupby("cell_type_l1", sort=False)
    }
    integrated_l2: list[str] = []
    integrated_l2_score: list[float] = []
    for i, l1_label in enumerate(l1_labels):
        allowed = l2_by_l1.get(str(l1_label), set())
        mask = np.asarray([label in allowed for label in l2_labels_all], dtype=bool)
        if not mask.any():
            integrated_l2.append("")
            integrated_l2_score.append(float("nan"))
            continue
        scores = query[i : i + 1] @ l2_matrix_all[mask].T
        idx = int(np.argmax(scores[0]))
        integrated_l2.append(l2_labels_all[mask][idx])
        integrated_l2_score.append(float(scores[0, idx]))
    metadata["integrated_cima_l2"] = integrated_l2
    metadata["integrated_cima_l2_score"] = integrated_l2_score
    return metadata


def project_atac_sample(sample_dir: Path, metadata: pd.DataFrame, feature_model: pd.DataFrame, n_components: int) -> np.ndarray:
    matrix_dir = sample_dir / "matrix"
    matrix_path = matrix_dir / "matrix.mtx"
    features_path = matrix_dir / "features.tsv.gz"
    barcodes_path = matrix_dir / "barcodes.tsv.gz"
    if not features_path.exists():
        features_path = matrix_dir / "features.tsv"
    if not barcodes_path.exists():
        barcodes_path = matrix_dir / "barcodes.tsv"

    feature_index = feature_model["feature_index"].astype(int).to_numpy() - 1
    feature_id = feature_model["feature_id"].astype(str).to_numpy()
    dim_cols = [c for c in feature_model.columns if c.startswith("dim_")][:n_components]
    idf = feature_model["idf"].to_numpy(dtype=np.float32)
    loadings = feature_model[dim_cols].to_numpy(dtype=np.float32)

    observed_features = np.asarray(read_lines(features_path), dtype=str)[feature_index]
    if not np.array_equal(observed_features, feature_id):
        raise ValueError(f"ATAC feature mismatch in {sample_dir}")

    barcodes = read_lines(barcodes_path)
    barcode_lookup = {barcode: i for i, barcode in enumerate(barcodes)}
    cell_barcodes = metadata["source_cell_id"].astype(str).tolist()
    cell_indices = [barcode_lookup[barcode] for barcode in cell_barcodes]

    matrix = sp.csr_matrix(scipy.io.mmread(matrix_path), dtype=np.float32)
    selected = matrix[feature_index, :][:, cell_indices]
    selected = sp.csc_matrix(selected, dtype=np.float32, copy=True)
    if selected.nnz:
        selected.data[:] = 1.0
        col_totals = np.asarray(selected.sum(axis=0)).ravel().astype(np.float32)
        col_totals[col_totals == 0] = 1.0
        selected.data /= np.repeat(col_totals, np.diff(selected.indptr))
    tfidf = sp.diags(idf, dtype=np.float32) @ selected
    return np.asarray(tfidf.T @ loadings, dtype=np.float32)


def sample_paths_from_metadata(cells: pd.DataFrame) -> list[tuple[str, str, Path, pd.DataFrame]]:
    groups = []
    for (gse, sample_id, source_dir), group in cells.groupby(["gse", "sample_id", "source_output_dir"], sort=True):
        groups.append((str(gse), str(sample_id), Path(str(source_dir)), group.copy()))
    return groups


def run_scanpy_tool_integration(
    embedding: np.ndarray, metadata: pd.DataFrame, config: IntegrationConfig
) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    obs = pd.DataFrame(
        {config.batch_key: metadata[config.batch_key].astype(str).to_numpy()},
        index=pd.Index(metadata.index.astype(str), dtype=object),
    )
    adata = ad.AnnData(obs=obs)
    adata.obsm["X_product"] = np.asarray(embedding, dtype=np.float32)
    method_detail: dict[str, Any] = {
        "requested_integration_method": config.integration_method,
        "neighbor_graph_tool": config.integration_method,
        "cluster_method": "scanpy.tl.leiden",
        "fallback_reason": "",
    }

    if config.integration_method == "bbknn":
        sce.pp.bbknn(
            adata,
            batch_key=config.batch_key,
            use_rep="X_product",
            n_pcs=adata.obsm["X_product"].shape[1],
            neighbors_within_batch=config.neighbors_within_batch,
            trim=config.bbknn_trim,
            approx=True,
            use_annoy=True,
            annoy_n_trees=20,
            metric="euclidean",
            computation="annoy",
        )
        neighbor_rep = "X_product"
    elif config.integration_method == "harmony":
        harmony_out = hm.run_harmony(
            adata.obsm["X_product"],
            adata.obs,
            [config.batch_key],
            random_state=config.random_state,
        )
        harmony_embedding = np.asarray(harmony_out.Z_corr, dtype=np.float32)
        if harmony_embedding.shape[0] != adata.n_obs and harmony_embedding.shape[1] == adata.n_obs:
            harmony_embedding = harmony_embedding.T
        adata.obsm["X_harmony"] = harmony_embedding
        sc.pp.neighbors(
            adata,
            n_neighbors=30,
            use_rep="X_harmony",
            metric="cosine",
            random_state=config.random_state,
        )
        method_detail["neighbor_graph_tool"] = "scanpy.pp.neighbors"
        method_detail["harmony_adjusted_basis"] = "X_harmony"
        neighbor_rep = "X_harmony"
    elif config.integration_method == "scanpy_neighbors":
        sc.pp.neighbors(
            adata,
            n_neighbors=30,
            use_rep="X_product",
            metric="cosine",
            random_state=config.random_state,
        )
        neighbor_rep = "X_product"
    else:
        raise ValueError(f"Unsupported integration method: {config.integration_method}")
    method_detail["neighbor_representation"] = neighbor_rep

    sc.tl.umap(adata, random_state=config.random_state, min_dist=0.25)
    sc.tl.leiden(
        adata,
        resolution=config.leiden_resolution,
        key_added="integrated_cluster",
        random_state=config.random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    coords = np.asarray(adata.obsm["X_umap"], dtype=np.float32)
    clusters = adata.obs["integrated_cluster"].astype(str).to_numpy()
    graph = adata.obsp.get("connectivities")
    if graph is not None:
        method_detail["neighbor_graph_nnz"] = int(graph.nnz)
        method_detail["neighbor_graph_bytes_estimated"] = int(
            graph.data.nbytes + graph.indices.nbytes + graph.indptr.nbytes
        )
    return coords, clusters, method_detail


def label_purity(metadata: pd.DataFrame, cluster_col: str, label_col: str) -> float:
    if label_col not in metadata.columns:
        return float("nan")
    total = 0
    dominant = 0
    for _, group in metadata.groupby(cluster_col, sort=False):
        labels = group[label_col].dropna().astype(str)
        if labels.empty:
            continue
        counts = labels.value_counts()
        total += int(counts.sum())
        dominant += int(counts.max())
    return float(dominant / total) if total else float("nan")


def mixing_metrics(integrated: np.ndarray, metadata: pd.DataFrame, random_state: int) -> dict[str, float]:
    n = len(metadata)
    if n < 2 or metadata["sample_id"].nunique() < 2:
        return {"fraction_cells_with_mixed_sample_neighbors": float("nan"), "mean_neighbor_sample_entropy": float("nan")}
    rng = np.random.default_rng(random_state)
    take = min(n, 5000)
    idx = rng.choice(np.arange(n), size=take, replace=False)
    distances = pairwise_distances(integrated[idx], integrated, metric="cosine", n_jobs=1)
    k = min(11, n)
    neighbor_idx = np.argpartition(distances, kth=k - 1, axis=1)[:, :k]
    samples = metadata["sample_id"].astype(str).to_numpy()
    mixed = []
    entropy = []
    for row_i, neighbors in zip(idx, neighbor_idx, strict=False):
        neighbors = neighbors[neighbors != row_i][:10]
        labels = samples[neighbors]
        counts = pd.Series(labels).value_counts(normalize=True).to_numpy(dtype=float)
        mixed.append(bool((labels != samples[row_i]).any()))
        entropy.append(float(-(counts * np.log2(counts + 1e-12)).sum()))
    return {
        "fraction_cells_with_mixed_sample_neighbors": float(np.mean(mixed)),
        "mean_neighbor_sample_entropy": float(np.mean(entropy)),
    }


def compute_metrics(metadata: pd.DataFrame, integrated: np.ndarray, original_x: str | None, original_y: str | None) -> dict[str, Any]:
    coord_values = metadata[["integrated_umap_1", "integrated_umap_2"]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
    coord_mask = np.isfinite(coord_values).all(axis=1)
    metadata = metadata.loc[coord_mask].copy()
    integrated = integrated[coord_mask]
    if metadata.empty:
        return {"n_integrated_cells_for_metrics": 0}
    clusters = metadata["integrated_cluster"].astype(str)
    counts = clusters.value_counts()
    l1_col = "integrated_cima_l1" if "integrated_cima_l1" in metadata.columns else None
    if l1_col is None:
        l1_col = "cima_cell_type_l1" if "cima_cell_type_l1" in metadata.columns else "azimuth_cima_l1"
    l2_col = "integrated_cima_l2" if "integrated_cima_l2" in metadata.columns else None
    if l2_col is None:
        l2_col = "cima_cell_type_l2" if "cima_cell_type_l2" in metadata.columns else "cima_l2"
    metrics: dict[str, Any] = {
        "umap_1_std": float(pd.to_numeric(metadata["integrated_umap_1"], errors="coerce").std()),
        "umap_2_std": float(pd.to_numeric(metadata["integrated_umap_2"], errors="coerce").std()),
        "umap_unique_coordinate_fraction": float(metadata[["integrated_umap_1", "integrated_umap_2"]].drop_duplicates().shape[0] / len(metadata)),
        "n_integrated_cells_for_metrics": int(len(metadata)),
        "n_integrated_clusters": int(counts.size),
        "largest_cluster_fraction": float(counts.max() / counts.sum()),
        "cluster_sample_id_purity": label_purity(metadata, "integrated_cluster", "sample_id"),
        "cluster_cima_l1_purity": label_purity(metadata, "integrated_cluster", l1_col),
        "cluster_cima_l2_purity": label_purity(metadata, "integrated_cluster", l2_col),
        "cluster_cima_l1_column": l1_col,
        "cluster_cima_l2_column": l2_col,
    }
    metrics.update(mixing_metrics(integrated, metadata, random_state=42))
    if original_x and original_y:
        old = metadata[[original_x, original_y]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
        new = metadata[["integrated_umap_1", "integrated_umap_2"]].to_numpy(dtype=np.float32)
        metrics["integrated_equals_original_umap"] = bool(np.allclose(old, new, equal_nan=True))
    return metrics


def original_umap_columns(metadata: pd.DataFrame) -> tuple[str | None, str | None]:
    for x_col, y_col in [("umap_1", "umap_2"), ("umap_atac_1", "umap_atac_2"), ("cima_ref_umap_1", "cima_ref_umap_2")]:
        if x_col in metadata.columns and y_col in metadata.columns:
            return x_col, y_col
    return None, None


def integrate_product_embeddings(config: IntegrationConfig) -> dict[str, Any]:
    cells_path = config.product_dir / "manifests" / "cells_metadata.csv"
    cells = pd.read_csv(cells_path, low_memory=False)
    reference_dir = config.data_root / "reference"
    feature_model = read_feature_model(reference_dir, config.product)
    embeddings: list[np.ndarray] = []
    metadata_parts: list[pd.DataFrame] = []

    for _, sample_id, source_dir, group in sample_paths_from_metadata(cells):
        if config.product in RNA_PRODUCTS:
            emb = project_rna_sample(source_dir, sample_id, feature_model, config.n_components)
        else:
            emb = project_atac_sample(source_dir, group, feature_model, config.n_components)
        if len(emb) != len(group):
            raise ValueError(f"Embedding row count mismatch for {source_dir}: {len(emb)} != {len(group)}")
        embeddings.append(emb)
        metadata_parts.append(group)

    metadata = pd.concat(metadata_parts, ignore_index=True, sort=False)
    embedding = np.vstack(embeddings).astype(np.float32)
    if config.product in RNA_PRODUCTS:
        metadata = assign_rna_cima_labels_from_embedding(
            metadata, embedding, reference_dir, embedding.shape[1]
        )
    integration_mask = np.ones(len(metadata), dtype=bool)
    exclusion_reason = np.full(len(metadata), "", dtype=object)
    if config.product == "only_rna" and config.rna_min_cima_l1_score > 0:
        scores = pd.to_numeric(metadata["integrated_cima_l1_score"], errors="coerce").to_numpy(dtype=float)
        integration_mask = np.isfinite(scores) & (scores >= config.rna_min_cima_l1_score)
        exclusion_reason[~integration_mask] = "low_integrated_cima_l1_score"
        if int(integration_mask.sum()) < 1000:
            raise ValueError(
                f"RNA integration confidence mask retained only {int(integration_mask.sum())} cells; "
                f"threshold={config.rna_min_cima_l1_score}"
            )
    metadata["integration_included"] = integration_mask
    metadata["integration_exclusion_reason"] = exclusion_reason
    umap_coords, clusters, method_detail = run_scanpy_tool_integration(
        embedding[integration_mask], metadata.loc[integration_mask].copy(), config
    )

    metadata["integrated_umap_1"] = np.nan
    metadata["integrated_umap_2"] = np.nan
    metadata["integrated_cluster"] = "not_integrated"
    metadata.loc[integration_mask, "integrated_umap_1"] = umap_coords[:, 0]
    metadata.loc[integration_mask, "integrated_umap_2"] = umap_coords[:, 1]
    metadata.loc[integration_mask, "integrated_cluster"] = clusters
    metadata["integration_method"] = config.integration_method
    metadata["integration_feature_space"] = "cima_rna_pca" if config.product in RNA_PRODUCTS else "cima_atac_lsi"
    metadata.to_csv(cells_path, index=False)

    integration_dir = config.product_dir / "integration"
    integration_dir.mkdir(parents=True, exist_ok=True)
    original_x, original_y = original_umap_columns(metadata)
    metrics = compute_metrics(metadata, embedding, original_x, original_y)
    pd.DataFrame([metrics]).to_csv(integration_dir / "integration_metrics.csv", index=False)
    metadata.loc[metadata["integration_included"].astype(bool)].groupby("sample_id")["integrated_cluster"].value_counts(normalize=True).rename("fraction").reset_index().to_csv(
        integration_dir / "sample_mixing_summary.csv", index=False
    )
    summary: dict[str, Any] = {
        "integration_status": "success",
        "product": config.product,
        "n_cells": int(len(metadata)),
        "n_cells_integrated": int(integration_mask.sum()),
        "n_cells_excluded_from_integration": int((~integration_mask).sum()),
        "rna_min_cima_l1_score": config.rna_min_cima_l1_score if config.product == "only_rna" else None,
        "n_samples": int(metadata["sample_id"].nunique()),
        "n_components": int(embedding.shape[1]),
        "embedding_bytes_estimated": int(embedding.nbytes),
        "integration_method": config.integration_method,
        "batch_key": config.batch_key,
        "neighbor_graph_tool": method_detail.get("neighbor_graph_tool", config.integration_method),
        "cluster_method": method_detail.get("cluster_method", "scanpy.tl.leiden"),
        "neighbors_within_batch": config.neighbors_within_batch if config.integration_method == "bbknn" else None,
        "bbknn_trim": config.bbknn_trim if config.integration_method == "bbknn" else None,
        "harmony_adjusted_basis": method_detail.get("harmony_adjusted_basis"),
        "leiden_resolution": config.leiden_resolution,
        "peak_memory_strategy": "per-sample projection + product-level low-dimensional BBKNN/Scanpy graph",
        "coordinate_source": "integrated_umap",
        "original_coordinate_columns": [original_x, original_y],
        **method_detail,
        "metrics": metrics,
    }
    (integration_dir / "integration_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    return summary
