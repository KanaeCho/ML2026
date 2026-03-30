from __future__ import annotations

# pyright: reportMissingImports=false, reportMissingModuleSource=false

import argparse
import csv
import gzip
import json
import random
from pathlib import Path
from typing import Any, cast

import h5py
import numpy as np
from scipy import sparse
from sklearn.decomposition import TruncatedSVD


def decode_categories(group: Any) -> tuple[list[str], list[str]]:
    categories = [
        x.decode() if isinstance(x, bytes) else str(x) for x in group["categories"][:]
    ]
    codes = group["codes"][:]
    values = [categories[int(code)] if code >= 0 else "" for code in codes]
    return categories, values


def write_feature_model(
    path: Path,
    keep: np.ndarray,
    feature_ids: list[str],
    idf: np.ndarray,
    components: np.ndarray,
) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            ["feature_index", "feature_id", "idf"]
            + [f"dim_{i + 1}" for i in range(components.shape[0])]
        )
        for position, feature_idx in enumerate(keep):
            writer.writerow(
                [
                    int(feature_idx + 1),
                    feature_ids[int(feature_idx)],
                    float(idf[feature_idx]),
                ]
                + [float(x) for x in components[:, position]]
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
        writer.writerow(
            ["cell_type_l1", "cell_type_l2", "cell_type_l3", "cell_type_l4"]
        )
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
    for idx, label in enumerate(labels):
        if label:
            groups[label].append(idx)

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [f"cell_type_{level_name}"]
            + [f"dim_{i + 1}" for i in range(embeddings.shape[1])]
        )
        for label in unique_labels:
            writer.writerow(
                [label] + [float(x) for x in embeddings[groups[label]].mean(axis=0)]
            )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--per-l4", type=int, default=100)
    parser.add_argument("--top-features", type=int, default=12000)
    parser.add_argument("--components", type=int, default=31)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    random.seed(args.seed)
    np.random.seed(args.seed)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(args.h5ad, "r") as handle:
        obs = cast(Any, handle["obs"])
        cell_type_l1_categories, cell_type_l1_values = decode_categories(
            cast(Any, obs["cell_type_l1"])
        )
        cell_type_l2_categories, cell_type_l2_values = decode_categories(
            cast(Any, obs["cell_type_l2"])
        )
        cell_type_l3_categories, cell_type_l3_values = decode_categories(
            cast(Any, obs["cell_type_l3"])
        )
        cell_type_l4_categories, cell_type_l4_values = decode_categories(
            cast(Any, obs["cell_type_l4"])
        )

        del cell_type_l1_categories, cell_type_l2_categories, cell_type_l3_categories

        by_l4: dict[str, list[int]] = {label: [] for label in cell_type_l4_categories}
        for idx, label in enumerate(cell_type_l4_values):
            if label:
                by_l4[label].append(idx)

        selected_rows: list[int] = []
        selected_l1: list[str] = []
        selected_l2: list[str] = []
        selected_l3: list[str] = []
        selected_l4: list[str] = []
        for label in cell_type_l4_categories:
            rows = by_l4[label]
            if len(rows) > args.per_l4:
                rows = random.sample(rows, args.per_l4)
                rows.sort()
            selected_rows.extend(rows)
            selected_l1.extend([cell_type_l1_values[row] for row in rows])
            selected_l2.extend([cell_type_l2_values[row] for row in rows])
            selected_l3.extend([cell_type_l3_values[row] for row in rows])
            selected_l4.extend([cell_type_l4_values[row] for row in rows])

        x_group = cast(Any, handle["X"])
        indptr = cast(Any, x_group["indptr"])[:]
        indices_ds = cast(Any, x_group["indices"])
        n_features = int(cast(Any, x_group.attrs["shape"])[1])
        feature_ids = [
            value.decode() if isinstance(value, bytes) else str(value)
            for value in cast(Any, cast(Any, handle["var"])["index"])[:]
        ]

        index_chunks: list[np.ndarray] = []
        row_indptr = [0]
        for row in selected_rows:
            start = int(indptr[row])
            end = int(indptr[row + 1])
            row_indices = np.asarray(indices_ds[start:end], dtype=np.int32)
            index_chunks.append(row_indices)
            row_indptr.append(row_indptr[-1] + len(row_indices))

    all_indices = (
        np.concatenate(index_chunks) if index_chunks else np.array([], dtype=np.int32)
    )
    all_data = np.ones(len(all_indices), dtype=np.float32)
    reference_matrix = sparse.csr_matrix(
        (all_data, all_indices, np.asarray(row_indptr, dtype=np.int64)),
        shape=(len(selected_rows), n_features),
    )

    row_sums = cast(
        Any, np.maximum(np.asarray(reference_matrix.sum(axis=1)).ravel(), 1.0)
    )
    row_scale = cast(Any, np.expand_dims(row_sums, axis=1))
    tf = cast(Any, reference_matrix.multiply(1.0 / row_scale).tocsr())
    document_frequency = cast(Any, np.asarray(reference_matrix.getnnz(axis=0)).ravel())
    idf = cast(
        Any, np.log1p(reference_matrix.shape[0] / np.maximum(document_frequency, 1))
    )
    tfidf = cast(Any, tf.multiply(idf).tocsr())
    feature_scores = np.asarray(tfidf.sum(axis=0)).ravel()
    keep = np.argsort(feature_scores)[::-1][: args.top_features]
    keep.sort()
    tfidf_subset = tfidf[:, keep].tocsr()

    svd = TruncatedSVD(
        n_components=args.components, algorithm="randomized", random_state=args.seed
    )
    embeddings = svd.fit_transform(tfidf_subset)

    write_feature_model(
        output_dir / "cima_atac_reference_lsi_features.tsv.gz",
        keep,
        feature_ids,
        idf,
        svd.components_,
    )
    write_hierarchy_csv(
        output_dir / "cima_atac_celltype_hierarchy.csv",
        cell_type_l1_values,
        cell_type_l2_values,
        cell_type_l3_values,
        cell_type_l4_values,
    )
    write_centroids(
        output_dir / "cima_atac_reference_l1_centroids.tsv",
        "l1",
        selected_l1,
        embeddings,
    )
    write_centroids(
        output_dir / "cima_atac_reference_l2_centroids.tsv",
        "l2",
        selected_l2,
        embeddings,
    )
    write_centroids(
        output_dir / "cima_atac_reference_l3_centroids.tsv",
        "l3",
        selected_l3,
        embeddings,
    )
    write_centroids(
        output_dir / "cima_atac_reference_l4_centroids.tsv",
        "l4",
        selected_l4,
        embeddings,
    )

    metadata = {
        "source_h5ad": str(Path(args.h5ad).resolve()),
        "per_l4": args.per_l4,
        "n_reference_cells": len(selected_rows),
        "top_features": args.top_features,
        "components": args.components,
        "seed": args.seed,
    }
    (output_dir / "cima_atac_reference_model.json").write_text(
        json.dumps(metadata, indent=2)
    )
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
