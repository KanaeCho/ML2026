from __future__ import annotations

# pyright: reportMissingImports=false, reportMissingModuleSource=false

import argparse
import csv
import gzip
import json
import random
import re
from pathlib import Path
from typing import Any, cast

import h5py
import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA


def decode_categories(group: Any) -> tuple[list[str], list[str]]:
    categories = [
        value.decode() if isinstance(value, bytes) else str(value)
        for value in group["categories"][:]
    ]
    codes = group["codes"][:]
    values = [categories[int(code)] if code >= 0 else "" for code in codes]
    return categories, values


def decode_vector(node: Any) -> list[str]:
    if isinstance(node, h5py.Group) and {"categories", "codes"} <= set(node.keys()):
        return decode_categories(node)[1]
    values = node[:]
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def parse_csv_arg(value: str | None) -> list[str] | None:
    if value is None:
        return None
    items = [item.strip() for item in value.split(",") if item.strip()]
    return items or None


def select_rows_balanced_by_l1(
    by_l1_l4: dict[str, dict[str, list[tuple[int, str, str]]]],
    per_l1_total: int | None,
    per_l4: int,
) -> tuple[list[int], list[str], list[str], list[str], list[str]]:
    selected_rows: list[int] = []
    selected_l1: list[str] = []
    selected_l2: list[str] = []
    selected_l3: list[str] = []
    selected_l4: list[str] = []

    for l1_label in sorted(by_l1_l4):
        l4_groups = {key: list(value) for key, value in by_l1_l4[l1_label].items() if value}
        if not l4_groups:
            continue

        chosen: dict[str, list[tuple[int, str, str]]] = {label: [] for label in l4_groups}
        if per_l1_total is None:
            for l4_label, rows in l4_groups.items():
                take = min(len(rows), per_l4)
                chosen[l4_label] = sorted(rows[:take], key=lambda item: item[0])
        else:
            pointers = {label: 0 for label in l4_groups}
            cycle_labels = list(sorted(l4_groups))
            while sum(len(rows) for rows in chosen.values()) < per_l1_total:
                progressed = False
                for l4_label in cycle_labels:
                    pointer = pointers[l4_label]
                    rows = l4_groups[l4_label]
                    if pointer >= len(rows):
                        continue
                    chosen[l4_label].append(rows[pointer])
                    pointers[l4_label] += 1
                    progressed = True
                    if sum(len(rows) for rows in chosen.values()) >= per_l1_total:
                        break
                if not progressed:
                    break

        for l4_label in sorted(chosen):
            rows = sorted(chosen[l4_label], key=lambda item: item[0])
            selected_rows.extend([row_idx for row_idx, _, _ in rows])
            selected_l1.extend([l1_label] * len(rows))
            selected_l2.extend([l2 for _, l2, _ in rows])
            selected_l3.extend([l3 for _, _, l3 in rows])
            selected_l4.extend([l4_label] * len(rows))

    return selected_rows, selected_l1, selected_l2, selected_l3, selected_l4


def choose_hvgs(var_group: Any, n_hvg: int) -> np.ndarray:
    if "variances_norm" not in var_group:
        raise KeyError("RNA reference var is missing variances_norm")
    scores = np.asarray(var_group["variances_norm"][:], dtype=np.float64)
    finite = np.isfinite(scores)
    if not finite.any():
        raise ValueError("No finite variances_norm values found in RNA reference")
    scores[~finite] = -np.inf
    order = np.argsort(scores)[::-1]
    keep = order[: min(n_hvg, len(order))]
    return np.sort(keep.astype(np.int32))


def write_feature_model(
    path: Path,
    keep: np.ndarray,
    feature_ids: list[str],
    gene_mean: np.ndarray,
    gene_std: np.ndarray,
    components: np.ndarray,
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["feature_index", "feature_id", "gene_mean", "gene_std"]
            + [f"pc_dim_{index + 1}" for index in range(components.shape[0])]
        )
        for position, feature_idx in enumerate(keep):
            writer.writerow(
                [
                    int(feature_idx + 1),
                    feature_ids[int(feature_idx)],
                    float(gene_mean[position]),
                    float(gene_std[position]),
                ]
                + [float(value) for value in components[:, position]]
            )


def write_hierarchy_csv(
    path: Path,
    level_1: list[str],
    level_2: list[str],
    level_3: list[str],
    level_4: list[str],
) -> None:
    seen: set[tuple[str, str, str, str]] = set()
    rows: list[tuple[str, str, str, str]] = []
    for values in zip(level_1, level_2, level_3, level_4):
        if "" in values:
            continue
        if values in seen:
            continue
        seen.add(values)
        rows.append(values)

    rows.sort()
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell_type_l1", "cell_type_l2", "cell_type_l3", "cell_type_l4"])
        writer.writerows(rows)


def write_centroids(
    path: Path, level_name: str, labels: list[str], embeddings: np.ndarray
) -> None:
    unique_labels: list[str] = []
    seen: set[str] = set()
    for label in labels:
        if label and label not in seen:
            unique_labels.append(label)
            seen.add(label)

    groups: dict[str, list[int]] = {label: [] for label in unique_labels}
    for index, label in enumerate(labels):
        if label:
            groups[label].append(index)

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [f"cell_type_{level_name}"]
            + [f"pc_dim_{index + 1}" for index in range(embeddings.shape[1])]
        )
        for label in unique_labels:
            writer.writerow(
                [label] + [float(value) for value in embeddings[groups[label]].mean(axis=0)]
            )


def build_selected_matrix(handle: h5py.File, selected_rows: list[int]) -> sparse.csr_matrix:
    x_group = cast(Any, handle["X"])
    indptr = cast(Any, x_group["indptr"])[:]
    indices_ds = cast(Any, x_group["indices"])
    data_ds = cast(Any, x_group["data"])
    n_features = int(cast(Any, x_group.attrs["shape"])[1])

    index_chunks: list[np.ndarray] = []
    data_chunks: list[np.ndarray] = []
    row_indptr = [0]
    for row in selected_rows:
        start = int(indptr[row])
        end = int(indptr[row + 1])
        row_indices = np.asarray(indices_ds[start:end], dtype=np.int32)
        row_data = np.asarray(data_ds[start:end], dtype=np.float32)
        index_chunks.append(row_indices)
        data_chunks.append(row_data)
        row_indptr.append(row_indptr[-1] + len(row_indices))

    all_indices = np.concatenate(index_chunks) if index_chunks else np.array([], dtype=np.int32)
    all_data = np.concatenate(data_chunks) if data_chunks else np.array([], dtype=np.float32)
    return sparse.csr_matrix(
        (all_data, all_indices, np.asarray(row_indptr, dtype=np.int64)),
        shape=(len(selected_rows), n_features),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--per-l4", type=int, default=200)
    parser.add_argument("--per-l1-total", type=int, default=None)
    parser.add_argument("--include-l1", type=str, default=None)
    parser.add_argument("--exclude-l4-regex", type=str, default=None)
    parser.add_argument("--n-hvg", type=int, default=3000)
    parser.add_argument("--components", type=int, default=50)
    parser.add_argument("--target-sum", type=float, default=10000.0)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    include_l1 = parse_csv_arg(args.include_l1)
    exclude_l4_regex = re.compile(args.exclude_l4_regex) if args.exclude_l4_regex else None

    with h5py.File(args.h5ad, "r") as handle:
        obs = cast(Any, handle["obs"])
        cell_type_l1_values = decode_vector(cast(Any, obs["cell_type_l1"]))
        cell_type_l2_values = decode_vector(cast(Any, obs["cell_type_l2"]))
        cell_type_l3_values = decode_vector(cast(Any, obs["cell_type_l3"]))
        cell_type_l4_values = decode_vector(cast(Any, obs["cell_type_l4"]))

        by_l1_l4: dict[str, dict[str, list[tuple[int, str, str]]]] = {}
        kept_l4: set[str] = set()
        for index, l4_label in enumerate(cell_type_l4_values):
            if not l4_label:
                continue
            l1_label = cell_type_l1_values[index]
            if include_l1 is not None and l1_label not in include_l1:
                continue
            if exclude_l4_regex is not None and exclude_l4_regex.search(l4_label):
                continue
            kept_l4.add(l4_label)
            by_l1_l4.setdefault(l1_label, {}).setdefault(l4_label, []).append(
                (index, cell_type_l2_values[index], cell_type_l3_values[index])
            )

        selected_rows, selected_l1, selected_l2, selected_l3, selected_l4 = (
            select_rows_balanced_by_l1(by_l1_l4, args.per_l1_total, args.per_l4)
        )
        if not selected_rows:
            raise ValueError("No RNA reference cells selected after applying filters")

        selected_matrix = build_selected_matrix(handle, selected_rows)
        var = cast(Any, handle["var"])
        feature_ids = decode_vector(cast(Any, var["_index"]))
        keep = choose_hvgs(var, args.n_hvg)

    selected_matrix = selected_matrix[:, keep]
    dense_reference = np.asarray(selected_matrix.toarray(), dtype=np.float32)

    gene_mean = dense_reference.mean(axis=0)
    gene_std = dense_reference.std(axis=0)
    gene_std[gene_std == 0] = 1.0
    scaled = (dense_reference - gene_mean) / gene_std

    n_components = min(args.components, scaled.shape[0] - 1, scaled.shape[1])
    if n_components < 2:
        raise ValueError("Need at least 2 PCA components to build RNA reference")

    pca = PCA(n_components=n_components, svd_solver="randomized", random_state=args.seed)
    embeddings = pca.fit_transform(scaled)
    components = np.asarray(pca.components_, dtype=np.float32)

    write_feature_model(
        output_dir / "cima_rna_reference_pca_features.tsv.gz",
        keep=keep,
        feature_ids=feature_ids,
        gene_mean=gene_mean,
        gene_std=gene_std,
        components=components,
    )
    write_hierarchy_csv(
        output_dir / "cima_rna_celltype_hierarchy.csv",
        selected_l1,
        selected_l2,
        selected_l3,
        selected_l4,
    )
    write_centroids(output_dir / "cima_rna_reference_l1_centroids.tsv", "l1", selected_l1, embeddings)
    write_centroids(output_dir / "cima_rna_reference_l2_centroids.tsv", "l2", selected_l2, embeddings)
    write_centroids(output_dir / "cima_rna_reference_l3_centroids.tsv", "l3", selected_l3, embeddings)
    write_centroids(output_dir / "cima_rna_reference_l4_centroids.tsv", "l4", selected_l4, embeddings)

    model_metadata = {
        "source_h5ad": str(Path(args.h5ad).resolve()),
        "modality": "RNA",
        "per_l4": args.per_l4,
        "per_l1_total": args.per_l1_total,
        "include_l1": include_l1,
        "exclude_l4_regex": args.exclude_l4_regex,
        "kept_l4_labels": len(kept_l4),
        "n_reference_cells": len(selected_rows),
        "n_hvg": args.n_hvg,
        "actual_n_hvg": int(len(keep)),
        "n_pca_components": int(n_components),
        "target_sum": args.target_sum,
        "seed": args.seed,
        "explained_variance_total": float(np.sum(pca.explained_variance_ratio_)),
    }
    (output_dir / "cima_rna_reference_model.json").write_text(
        json.dumps(model_metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    print(
        f"Built RNA reference model with {len(selected_rows)} cells, "
        f"{len(keep)} genes, {n_components} PCs"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
