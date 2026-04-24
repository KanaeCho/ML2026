from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread


PALETTE = {
    "B": "#2C7BB6",
    "CD4_T": "#D7191C",
    "CD8_T&unconvensional_T": "#FDAE61",
    "ILC": "#762A83",
    "Myeloid": "#1A9641",
}


def normalize_rows(values: np.ndarray) -> np.ndarray:
    norms = np.linalg.norm(values, axis=1, keepdims=True)
    norms[norms == 0] = 1.0
    return values / norms


def load_barcodes(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.strip() for line in handle]


def load_feature_model(reference_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    cima_dir = reference_dir
    feature_model = pd.read_csv(
        cima_dir / "cima_atac_reference_lsi_features.tsv.gz", sep="\t"
    )
    l1_centroids = pd.read_csv(
        cima_dir / "cima_atac_reference_l1_centroids.tsv", sep="\t"
    )
    return feature_model, l1_centroids


def compute_cluster_labels(
    sample_dir: Path, reference_dir: Path
) -> tuple[dict[str, str], list[dict[str, object]], pd.DataFrame]:
    metadata_qc = pd.read_csv(sample_dir / "metadata_qc.csv")
    barcodes = load_barcodes(sample_dir / "matrix" / "barcodes.tsv.gz")
    counts = mmread(sample_dir / "matrix" / "matrix.mtx").tocsr().astype(np.float64)
    feature_model, l1_centroids = load_feature_model(reference_dir)

    selected_idx = feature_model["feature_index"].to_numpy(dtype=int) - 1
    counts = counts[selected_idx, :]
    if counts.nnz > 0:
        counts.data[:] = 1.0

    cell_peak_totals = np.asarray(counts.sum(axis=0)).ravel()
    cell_peak_totals[cell_peak_totals == 0] = 1.0
    tf = counts @ sparse.diags(1.0 / cell_peak_totals)
    tfidf = sparse.diags(feature_model["idf"].to_numpy(dtype=float)) @ tf

    loading_cols = [col for col in feature_model.columns if col.startswith("dim_")]
    query_embeddings = tfidf.T @ feature_model[loading_cols].to_numpy(dtype=float)
    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
    cell_idx = np.array(
        [barcode_to_idx[barcode] for barcode in metadata_qc["cell_barcode"]]
    )
    query_embeddings = query_embeddings[cell_idx, :]

    metadata_qc = metadata_qc.copy()
    metadata_qc["seurat_clusters"] = metadata_qc["seurat_clusters"].astype(str)

    cluster_embeddings = []
    cluster_names: list[str] = []
    for cluster, sub in metadata_qc.groupby("seurat_clusters", sort=True):
        cluster_names.append(cluster)
        cluster_embeddings.append(
            np.asarray(query_embeddings[sub.index, :]).mean(axis=0)
        )
    cluster_embeddings = np.vstack(cluster_embeddings)

    centroid_embeddings = l1_centroids[
        [col for col in l1_centroids.columns if col.startswith("dim_")]
    ].to_numpy(dtype=float)
    cluster_similarity = (
        normalize_rows(cluster_embeddings[:, 1:])
        @ normalize_rows(centroid_embeddings[:, 1:]).T
    )
    top_index = cluster_similarity.argmax(axis=1)
    top_scores = cluster_similarity[np.arange(cluster_similarity.shape[0]), top_index]
    margins = []
    for row in cluster_similarity:
        sorted_row = np.sort(row)[::-1]
        margins.append(sorted_row[0] - sorted_row[1] if len(sorted_row) > 1 else np.nan)

    labels = l1_centroids["cell_type_l1"].tolist()
    cluster_labels = {
        cluster_names[i]: labels[top_index[i]] for i in range(len(cluster_names))
    }
    summaries = []
    for i, cluster in enumerate(cluster_names):
        summaries.append(
            {
                "seurat_cluster": cluster,
                "cima_cluster_centroid_l1": cluster_labels[cluster],
                "n_cells": int((metadata_qc["seurat_clusters"] == cluster).sum()),
                "similarity_score": float(top_scores[i]),
                "similarity_margin": float(margins[i]),
            }
        )
    return cluster_labels, summaries, metadata_qc


def write_cluster_annotation(
    metadata_qc: pd.DataFrame, cluster_labels: dict[str, str], output_path: Path
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "celltype", "seurat_cluster"])
        for _, row in metadata_qc.iterrows():
            writer.writerow(
                [
                    row["cell_barcode"],
                    cluster_labels[row["seurat_clusters"]],
                    row["seurat_clusters"],
                ]
            )


def write_cluster_summary(
    summary_rows: list[dict[str, object]], output_path: Path
) -> None:
    fieldnames = [
        "seurat_cluster",
        "cima_cluster_centroid_l1",
        "n_cells",
        "similarity_score",
        "similarity_margin",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


def render_umap(
    metadata_qc: pd.DataFrame,
    cluster_labels: dict[str, str],
    output_path: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))
    by_label: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for _, row in metadata_qc.iterrows():
        label = cluster_labels[row["seurat_clusters"]]
        by_label[label].append((float(row["umap_atac_1"]), float(row["umap_atac_2"])))

    for label in ["B", "CD4_T", "CD8_T&unconvensional_T", "ILC", "Myeloid"]:
        pts = by_label.get(label)
        if not pts:
            continue
        xs, ys = zip(*pts)
        ax.scatter(xs, ys, s=4, c=PALETTE[label], label=label, linewidths=0, alpha=0.85)

    ax.set_title(title, fontsize=18, weight="bold")
    ax.set_xlabel("UMAP_1")
    ax.set_ylabel("UMAP_2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Render TEA-seq CIMA cluster-centroid broad labels on ATAC UMAP"
    )
    parser.add_argument("--gse", required=True)
    parser.add_argument("--gsm", required=True)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--sample-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--reference-dir", default=None)
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    sample_dir = (
        Path(args.sample_dir).resolve()
        if args.sample_dir
        else project_root / "output" / args.gse / args.gsm
    )
    output_dir = Path(args.output_dir).resolve() if args.output_dir else sample_dir
    reference_dir = (
        Path(args.reference_dir).resolve()
        if args.reference_dir
        else project_root / "data" / "reference" / "cima"
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    cluster_labels, summaries, metadata_qc = compute_cluster_labels(
        sample_dir, reference_dir
    )

    annotation_path = output_dir / "cima_cluster_centroid_labels.csv"
    summary_path = output_dir / "cima_cluster_centroid_label_summary.csv"
    umap_path = output_dir / "umap_atac_cima_cell_type_l1_cluster_centroid.png"

    write_cluster_annotation(metadata_qc, cluster_labels, annotation_path)
    write_cluster_summary(summaries, summary_path)
    render_umap(
        metadata_qc,
        cluster_labels,
        umap_path,
        f"{args.gsm} query-native UMAP by CIMA L1 (cluster centroid)",
    )

    print(f"Annotation: {annotation_path}")
    print(f"Summary: {summary_path}")
    print(f"UMAP: {umap_path}")


if __name__ == "__main__":
    main()
