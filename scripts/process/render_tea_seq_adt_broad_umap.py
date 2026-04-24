from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, cast

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csc_matrix


PALETTE = {
    "B": "#2C7BB6",
    "CD4_T": "#D7191C",
    "CD8_T&unconvensional_T": "#FDAE61",
    "ILC": "#762A83",
    "Myeloid": "#1A9641",
    "Unknown": "#7F7F7F",
}

MARKERS = [
    "TCRab",
    "CD4",
    "CD8a",
    "CD56",
    "CD161",
    "CD11b",
    "CD11c",
    "HLA.DR",
    "TCRgd",
    "CD21",
    "CD24",
]


def normalize_barcode(value: str) -> str:
    value = str(value)
    return value[:-2] if value.endswith("-1") else value


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def find_h5(raw_dir: Path, gsm: str) -> Path:
    matches = sorted(raw_dir.glob(f"{gsm}_*.h5"))
    if not matches:
        raise FileNotFoundError(f"No h5 file found for {gsm} in {raw_dir}")
    if len(matches) > 1:
        raise RuntimeError(f"Multiple h5 files found for {gsm}: {matches}")
    return matches[0]


def assign_cluster_label(marker_medians: dict[str, float]) -> str:
    tcrab = marker_medians.get("TCRab", math.nan)
    cd4 = marker_medians.get("CD4", math.nan)
    cd8 = marker_medians.get("CD8a", math.nan)
    cd56 = marker_medians.get("CD56", math.nan)
    cd161 = marker_medians.get("CD161", math.nan)
    cd11b = marker_medians.get("CD11b", math.nan)
    cd11c = marker_medians.get("CD11c", math.nan)
    hla_dr = marker_medians.get("HLA.DR", math.nan)
    tcrgd = marker_medians.get("TCRgd", math.nan)
    cd21 = marker_medians.get("CD21", math.nan)
    cd24 = marker_medians.get("CD24", math.nan)

    def safe(value: float) -> float:
        return 0.0 if math.isnan(value) else value

    if tcrab >= 3 and cd4 - cd8 >= 1.5 and cd56 < 1.2 and cd11b < 1.2:
        return "CD4_T"
    if tcrab >= 3 and cd8 - cd4 >= 1.5 and cd56 < 1.5:
        return "CD8_T&unconvensional_T"
    if cd56 >= 1.5 or cd161 >= 2.5:
        return "ILC"
    if cd11b >= 1.5 or cd11c >= 1.5 or hla_dr >= 1.5:
        return "Myeloid"
    if cd21 >= 1.5 or cd24 >= 1.5:
        return "B"

    fallback_scores = {
        "CD4_T": safe(tcrab) + safe(cd4) + max(safe(cd4) - safe(cd8), 0.0),
        "CD8_T&unconvensional_T": safe(tcrab)
        + safe(cd8)
        + max(safe(cd8) - safe(cd4), 0.0)
        + 0.5 * safe(tcrgd),
        "ILC": safe(cd56) + safe(cd161) + 0.5 * safe(tcrgd),
        "Myeloid": safe(cd11b) + safe(cd11c) + safe(hla_dr),
        "B": safe(cd21) + safe(cd24),
    }
    return max(fallback_scores.items(), key=lambda item: item[1])[0]


def compute_cluster_labels(
    metadata_rows: list[dict[str, str]], h5_path: Path
) -> tuple[dict[str, str], list[dict[str, object]]]:
    orig_to_cluster = {
        normalize_barcode(row["filtered_metadata_original_barcodes"]): row[
            "seurat_clusters"
        ]
        for row in metadata_rows
    }

    with h5py.File(h5_path, "r") as handle:
        matrix_group = cast(Any, handle["matrix"])
        observations_group = cast(Any, matrix_group["observations"])
        orig_barcodes = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in observations_group["original_barcodes"][()]
        ]
        adt = cast(Any, handle["ADT"])
        matrix = csc_matrix(
            (adt["data"][()], adt["indices"][()], adt["indptr"][()]),
            shape=tuple(int(v) for v in adt["shape"][()]),
        )
        features_group = cast(Any, adt["features"])
        features = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in features_group["id"][()]
        ]

    values = np.log1p(matrix.toarray().astype(float))
    feature_idx = {name: idx for idx, name in enumerate(features)}
    cluster_to_cols: dict[str, list[int]] = defaultdict(list)
    for col_idx, original in enumerate(orig_barcodes):
        cluster = orig_to_cluster.get(normalize_barcode(original))
        if cluster is not None:
            cluster_to_cols[cluster].append(col_idx)

    cluster_labels: dict[str, str] = {}
    summaries: list[dict[str, object]] = []
    for cluster in sorted(cluster_to_cols, key=lambda value: int(value)):
        cols = cluster_to_cols[cluster]
        medians = {
            marker: float(np.median(values[feature_idx[marker], cols]))
            for marker in MARKERS
            if marker in feature_idx
        }
        label = assign_cluster_label(medians)
        cluster_labels[cluster] = label
        summaries.append(
            {
                "seurat_cluster": cluster,
                "adt_broad_label": label,
                "n_cells": len(cols),
                **{
                    f"median_{marker}": medians.get(marker, math.nan)
                    for marker in MARKERS
                },
            }
        )
    return cluster_labels, summaries


def write_cluster_annotation(
    metadata_rows: list[dict[str, str]],
    cluster_labels: dict[str, str],
    output_path: Path,
) -> None:
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "celltype", "seurat_cluster"])
        for row in metadata_rows:
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
        "adt_broad_label",
        "n_cells",
        *[f"median_{marker}" for marker in MARKERS],
    ]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)


def render_umap(
    metadata_rows: list[dict[str, str]],
    cluster_labels: dict[str, str],
    output_path: Path,
    title: str,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 8))
    points_by_label: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in metadata_rows:
        label = cluster_labels[row["seurat_clusters"]]
        points_by_label[label].append(
            (float(row["umap_atac_1"]), float(row["umap_atac_2"]))
        )

    order = ["B", "CD4_T", "CD8_T&unconvensional_T", "ILC", "Myeloid", "Unknown"]
    for label in order:
        pts = points_by_label.get(label)
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
        description="Render TEA-seq ADT-derived broad labels on ATAC UMAP"
    )
    parser.add_argument("--gse", required=True)
    parser.add_argument("--gsm", required=True)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--sample-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    project_root = Path(args.project_root).resolve()
    sample_dir = (
        Path(args.sample_dir).resolve()
        if args.sample_dir
        else project_root / "output" / args.gse / args.gsm
    )
    output_dir = Path(args.output_dir).resolve() if args.output_dir else sample_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata_qc_path = sample_dir / "metadata_qc.csv"
    if not metadata_qc_path.exists():
        raise FileNotFoundError(f"metadata_qc.csv not found: {metadata_qc_path}")
    metadata_rows = read_csv_rows(metadata_qc_path)
    required_cols = {
        "cell_barcode",
        "seurat_clusters",
        "filtered_metadata_original_barcodes",
        "umap_atac_1",
        "umap_atac_2",
    }
    missing = required_cols.difference(metadata_rows[0].keys())
    if missing:
        raise RuntimeError(
            f"metadata_qc.csv missing required columns: {sorted(missing)}"
        )

    h5_path = find_h5(project_root / "data" / "raw" / args.gse, args.gsm)
    cluster_labels, summary_rows = compute_cluster_labels(metadata_rows, h5_path)

    annotation_path = output_dir / "adt_cluster_broad_labels.csv"
    summary_path = output_dir / "adt_cluster_broad_label_summary.csv"
    umap_path = output_dir / "umap_atac_adt_cluster_broad_celltype.png"

    write_cluster_annotation(metadata_rows, cluster_labels, annotation_path)
    write_cluster_summary(summary_rows, summary_path)
    render_umap(
        metadata_rows,
        cluster_labels,
        umap_path,
        f"{args.gsm} ATAC UMAP by ADT broad label",
    )

    print(f"Annotation: {annotation_path}")
    print(f"Summary: {summary_path}")
    print(f"UMAP: {umap_path}")


if __name__ == "__main__":
    main()
