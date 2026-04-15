from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(frozen=True)
class CimaReference:
    feature_ids: list[str]
    gene_mean: np.ndarray
    gene_std: np.ndarray
    loadings: np.ndarray
    pc_columns: list[str]
    l1_centroids: pd.DataFrame
    l2_centroids: pd.DataFrame
    l2_by_l1: dict[str, list[str]]


ANNOTATION_STRING_COLUMNS = ["cima_l1", "cima_l2", "cima_l1_masked"]
ANNOTATION_FLOAT_COLUMNS = [
    "cima_l1_score",
    "cima_l1_score_margin",
    "cima_l2_score",
    "cima_l2_score_margin",
]
ANNOTATION_BOOL_COLUMNS = ["cima_l1_low_confidence"]


def _read_feature_model(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="gzip")


def _normalize_rows(matrix: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(matrix, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return matrix / norms


def load_cima_reference(reference_dir: Path) -> CimaReference:
    cima_dir = Path(reference_dir) / "cima"
    feature_model = _read_feature_model(
        cima_dir / "cima_rna_reference_pca_features.tsv.gz"
    )
    l1_centroids = pd.read_csv(
        cima_dir / "cima_rna_reference_l1_centroids.tsv", sep="\t"
    )
    l2_centroids = pd.read_csv(
        cima_dir / "cima_rna_reference_l2_centroids.tsv", sep="\t"
    )
    hierarchy = pd.read_csv(cima_dir / "cima_rna_celltype_hierarchy.csv")

    pc_columns = [
        column for column in feature_model.columns if column.startswith("pc_dim_")
    ]
    l2_by_l1 = {
        str(l1_label): sorted(group["cell_type_l2"].astype(str).unique().tolist())
        for l1_label, group in hierarchy.groupby("cell_type_l1", sort=False)
    }

    return CimaReference(
        feature_ids=feature_model["feature_id"].astype(str).tolist(),
        gene_mean=feature_model["gene_mean"].to_numpy(dtype=float),
        gene_std=np.where(
            feature_model["gene_std"].to_numpy(dtype=float) == 0,
            1.0,
            feature_model["gene_std"].to_numpy(dtype=float),
        ),
        loadings=feature_model[pc_columns].to_numpy(dtype=float),
        pc_columns=pc_columns,
        l1_centroids=l1_centroids.set_index("cell_type_l1", drop=True),
        l2_centroids=l2_centroids.set_index("cell_type_l2", drop=True),
        l2_by_l1=l2_by_l1,
    )


def _pass_qc_mask(adata: ad.AnnData) -> pd.Series:
    if "pass_qc" not in adata.obs:
        raise KeyError("adata.obs must contain 'pass_qc'")
    return adata.obs["pass_qc"].fillna(False).astype(bool)


def _initialize_annotation_columns(obs: pd.DataFrame) -> None:
    for column in ANNOTATION_STRING_COLUMNS:
        obs[column] = pd.Series(pd.NA, index=obs.index, dtype="string")
    for column in ANNOTATION_FLOAT_COLUMNS:
        obs[column] = pd.Series(np.nan, index=obs.index, dtype=float)
    for column in ANNOTATION_BOOL_COLUMNS:
        obs[column] = pd.Series(pd.NA, index=obs.index, dtype="boolean")


def _candidate_feature_keys(adata: ad.AnnData) -> list[pd.Index]:
    candidates: list[pd.Index] = []
    if "feature_name" in adata.var.columns:
        candidates.append(pd.Index(adata.var["feature_name"].astype(str), dtype=str))
    candidates.append(pd.Index(adata.var_names.astype(str), dtype=str))
    if "feature_id" in adata.var.columns:
        candidates.append(pd.Index(adata.var["feature_id"].astype(str), dtype=str))
    return candidates


def _build_query_matrix(adata: ad.AnnData, reference: CimaReference) -> np.ndarray:
    matrix = (
        adata.X.toarray().astype(float, copy=False)
        if sparse.issparse(adata.X)
        else np.asarray(adata.X, dtype=float)
    )
    query = np.zeros((adata.n_obs, len(reference.feature_ids)), dtype=float)
    reference_ids = set(reference.feature_ids)
    best_index = max(
        _candidate_feature_keys(adata),
        key=lambda idx: len(set(idx.astype(str)) & reference_ids),
    )
    feature_lookup = pd.Series(np.arange(adata.n_vars), index=best_index)
    feature_lookup = feature_lookup[~feature_lookup.index.duplicated(keep="first")]

    matched_reference_positions = []
    matched_query_positions = []
    for ref_idx, feature_id in enumerate(reference.feature_ids):
        query_idx = feature_lookup.get(feature_id)
        if query_idx is None or pd.isna(query_idx):
            continue
        matched_reference_positions.append(ref_idx)
        matched_query_positions.append(int(query_idx))

    if matched_reference_positions:
        query[:, matched_reference_positions] = matrix[:, matched_query_positions]
    return query


def _project_query(adata: ad.AnnData, reference: CimaReference) -> np.ndarray:
    query = _build_query_matrix(adata, reference)
    scaled = (query - reference.gene_mean) / reference.gene_std
    return scaled @ reference.loadings


def _top_two_scores(
    similarity: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    order = np.argsort(-similarity, axis=1, kind="mergesort")
    top_index = order[:, 0]
    second_index = order[:, 1] if similarity.shape[1] >= 2 else top_index
    top_score = similarity[np.arange(similarity.shape[0]), top_index]
    margin = top_score - similarity[np.arange(similarity.shape[0]), second_index]
    return top_index, top_score, margin


def annotate_with_cima(adata: ad.AnnData, reference_dir: Path) -> ad.AnnData:
    out = adata.copy()
    _initialize_annotation_columns(out.obs)

    pass_qc_mask = _pass_qc_mask(out)
    if not pass_qc_mask.any():
        return out

    reference = load_cima_reference(reference_dir)
    pass_qc = out[pass_qc_mask].copy()
    query_embeddings = _project_query(pass_qc, reference)
    query_norm = _normalize_rows(query_embeddings)

    l1_matrix = reference.l1_centroids[reference.pc_columns].to_numpy(dtype=float)
    l1_norm = _normalize_rows(l1_matrix)
    l1_similarity = query_norm @ l1_norm.T
    l1_top_index, l1_scores, l1_margins = _top_two_scores(l1_similarity)
    l1_labels = reference.l1_centroids.index.to_numpy()[l1_top_index]

    l2_matrix = reference.l2_centroids[reference.pc_columns].to_numpy(dtype=float)
    l2_norm = _normalize_rows(l2_matrix)
    l2_labels = []
    l2_scores = []
    l2_margins = []

    for row_idx, l1_label in enumerate(l1_labels):
        allowed = reference.l2_by_l1.get(str(l1_label), [])
        if not allowed:
            l2_labels.append(pd.NA)
            l2_scores.append(np.nan)
            l2_margins.append(np.nan)
            continue

        keep = reference.l2_centroids.index.isin(allowed)
        allowed_labels = reference.l2_centroids.index[keep].to_numpy()
        allowed_norm = l2_norm[keep]
        similarities = query_norm[row_idx : row_idx + 1] @ allowed_norm.T
        top_idx, top_score, margin = _top_two_scores(similarities)
        l2_labels.append(allowed_labels[top_idx[0]])
        l2_scores.append(float(top_score[0]))
        l2_margins.append(float(margin[0]))

    low_confidence = (l1_scores < 0.35) | (l1_margins < 0.05)
    masked_l1 = np.where(low_confidence, "Unknown", l1_labels)

    pass_qc.obs["cima_l1"] = pd.Series(
        l1_labels, index=pass_qc.obs_names, dtype="string"
    )
    pass_qc.obs["cima_l2"] = pd.Series(
        l2_labels, index=pass_qc.obs_names, dtype="string"
    )
    pass_qc.obs["cima_l1_score"] = l1_scores
    pass_qc.obs["cima_l1_score_margin"] = l1_margins
    pass_qc.obs["cima_l2_score"] = np.asarray(l2_scores, dtype=float)
    pass_qc.obs["cima_l2_score_margin"] = np.asarray(l2_margins, dtype=float)
    pass_qc.obs["cima_l1_low_confidence"] = pd.Series(
        low_confidence, index=pass_qc.obs_names, dtype="boolean"
    )
    pass_qc.obs["cima_l1_masked"] = pd.Series(
        masked_l1, index=pass_qc.obs_names, dtype="string"
    )

    columns = (
        ANNOTATION_STRING_COLUMNS + ANNOTATION_FLOAT_COLUMNS + ANNOTATION_BOOL_COLUMNS
    )
    out.obs.loc[pass_qc.obs_names, columns] = pass_qc.obs[columns]
    return out


__all__ = ["CimaReference", "load_cima_reference", "annotate_with_cima"]
