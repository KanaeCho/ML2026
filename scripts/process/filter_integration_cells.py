#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"


def read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return [line.rstrip("\n") for line in handle]


def write_lines_gz(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for value in values:
            handle.write(f"{value}\n")


def parse_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def parse_comma_list(value: str | None) -> set[str]:
    if not value:
        return set()
    return {item.strip() for item in value.split(",") if item.strip()}


def parse_cluster_score_rules(values: list[str]) -> dict[str, float]:
    rules: dict[str, float] = {}
    for raw in values:
        for piece in raw.split(","):
            piece = piece.strip()
            if not piece:
                continue
            cluster, sep, threshold = piece.partition(":")
            if not sep:
                raise SystemExit(
                    f"Invalid --exclude-cluster-doublet-score-ge rule: {piece!r}; expected CLUSTER:THRESHOLD"
                )
            cluster = cluster.strip()
            threshold = threshold.strip()
            if not cluster or not threshold:
                raise SystemExit(
                    f"Invalid --exclude-cluster-doublet-score-ge rule: {piece!r}; expected CLUSTER:THRESHOLD"
                )
            rules[cluster] = float(threshold)
    return rules


def safe_float(value: str | None) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def iter_matrix_market_entries(path: Path):
    with path.open("r", encoding="utf-8") as handle:
        header_seen = False
        dims_seen = False
        for line in handle:
            if not header_seen:
                header_seen = True
                continue
            if line.startswith("%"):
                continue
            if not dims_seen:
                dims_seen = True
                continue
            yield line.rstrip("\n")


def parse_matrix_market_shape(path: Path) -> tuple[int, int, int]:
    with path.open("r", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n")
        if not header.startswith("%%MatrixMarket matrix coordinate"):
            raise ValueError(f"Unexpected MatrixMarket header in {path}")
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if len(parts) != 3:
                raise ValueError(f"Invalid matrix shape line in {path}: {line!r}")
            return int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError(f"Missing matrix shape in {path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Filter a merged integration input by annotated metadata columns."
    )
    parser.add_argument(
        "--input-dir",
        default=str(OUTPUT_ROOT / "integration_merged"),
        help="Merged integration input directory with matrix.mtx, features.tsv.gz, barcodes.tsv.gz, merged_metadata.csv",
    )
    parser.add_argument(
        "--annotated-metadata",
        default=str(
            OUTPUT_ROOT
            / "integration_merged_analysis"
            / "integrated_metadata_celltyped.csv.gz"
        ),
        help="Annotated cell-level metadata used for filtering",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Output directory for filtered integration input",
    )
    parser.add_argument(
        "--filter-column",
        default="celltype",
        help="Metadata column in --annotated-metadata to filter on",
    )
    parser.add_argument(
        "--exclude-values",
        default="Unknown / non-PBMC-like",
        help="Comma-separated values in --filter-column to exclude",
    )
    parser.add_argument(
        "--exclude-clusters",
        default="",
        help="Comma-separated seurat_clusters values to exclude",
    )
    parser.add_argument(
        "--exclude-doublet-class",
        action="store_true",
        default=False,
        help="Exclude cells whose scDblFinder.class is explicit doublet",
    )
    parser.add_argument(
        "--exclude-doublet-score-ge",
        type=float,
        default=None,
        help="Exclude cells with scDblFinder.score >= this threshold",
    )
    parser.add_argument(
        "--exclude-cluster-doublet-score-ge",
        action="append",
        default=[],
        help="Repeatable or comma-separated rules of the form CLUSTER:THRESHOLD",
    )
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    annotated_metadata_path = Path(args.annotated_metadata)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    matrix_path = input_dir / "matrix.mtx"
    features_path = input_dir / "features.tsv.gz"
    barcodes_path = input_dir / "barcodes.tsv.gz"
    merged_metadata_path = input_dir / "merged_metadata.csv"
    required = [
        matrix_path,
        features_path,
        barcodes_path,
        merged_metadata_path,
        annotated_metadata_path,
    ]
    if not all(path.exists() for path in required):
        raise SystemExit("Missing required merged matrix or annotated metadata inputs.")

    exclude_values = parse_comma_list(args.exclude_values)
    exclude_clusters = parse_comma_list(args.exclude_clusters)
    cluster_score_rules = parse_cluster_score_rules(
        args.exclude_cluster_doublet_score_ge
    )
    has_any_rule = any(
        [
            bool(exclude_values),
            bool(exclude_clusters),
            args.exclude_doublet_class,
            args.exclude_doublet_score_ge is not None,
            bool(cluster_score_rules),
        ]
    )
    if not has_any_rule:
        raise SystemExit("At least one filtering rule is required.")

    annotated_fields, annotated_rows = parse_csv(annotated_metadata_path)
    if exclude_values and args.filter_column not in annotated_fields:
        raise SystemExit(
            f"Missing filter column in annotated metadata: {args.filter_column}"
        )
    if "global_cell_id" not in annotated_fields:
        raise SystemExit("Annotated metadata must contain global_cell_id")
    if exclude_clusters and "seurat_clusters" not in annotated_fields:
        raise SystemExit(
            "Annotated metadata must contain seurat_clusters for --exclude-clusters"
        )
    if args.exclude_doublet_class and "scDblFinder.class" not in annotated_fields:
        raise SystemExit(
            "Annotated metadata must contain scDblFinder.class for --exclude-doublet-class"
        )
    if (
        args.exclude_doublet_score_ge is not None or cluster_score_rules
    ) and "scDblFinder.score" not in annotated_fields:
        raise SystemExit(
            "Annotated metadata must contain scDblFinder.score for score-based filters"
        )

    keep_ids: set[str] = set()
    excluded_counter: Counter[str] = Counter()
    kept_counter: Counter[str] = Counter()
    removal_reason_counts: Counter[str] = Counter()
    for row in annotated_rows:
        reasons: list[str] = []
        value = row.get(args.filter_column, "")
        cluster_value = row.get("seurat_clusters", "")
        score_value = safe_float(row.get("scDblFinder.score"))

        if exclude_values and value in exclude_values:
            reasons.append(f"{args.filter_column}:{value}")
            excluded_counter[value] += 1
        if exclude_clusters and cluster_value in exclude_clusters:
            reasons.append(f"cluster:{cluster_value}")
        if (
            args.exclude_doublet_class
            and row.get("scDblFinder.class", "").strip().lower() == "doublet"
        ):
            reasons.append("scDblFinder.class:doublet")
        if (
            args.exclude_doublet_score_ge is not None
            and score_value is not None
            and score_value >= args.exclude_doublet_score_ge
        ):
            reasons.append(f"scDblFinder.score>={args.exclude_doublet_score_ge:g}")
        cluster_threshold = cluster_score_rules.get(cluster_value)
        if (
            cluster_threshold is not None
            and score_value is not None
            and score_value >= cluster_threshold
        ):
            reasons.append(
                f"cluster:{cluster_value}:scDblFinder.score>={cluster_threshold:g}"
            )

        if reasons:
            for reason in reasons:
                removal_reason_counts[reason] += 1
            continue

        keep_ids.add(row["global_cell_id"])
        kept_counter[value] += 1

    if not keep_ids:
        raise SystemExit("No cells remain after filtering.")

    merged_fields, merged_rows = parse_csv(merged_metadata_path)
    if "global_cell_id" not in merged_fields:
        raise SystemExit("merged_metadata.csv must contain global_cell_id")

    merged_by_id = {row["global_cell_id"]: row for row in merged_rows}
    merged_barcodes = read_lines_gz(barcodes_path)

    selected_rows: list[dict[str, str]] = []
    selected_ids: list[str] = []
    position_map: dict[int, int] = {}
    for original_pos, global_cell_id in enumerate(merged_barcodes, start=1):
        if global_cell_id not in keep_ids:
            continue
        merged_row = merged_by_id.get(global_cell_id)
        if merged_row is None:
            raise SystemExit(f"Missing merged metadata row for {global_cell_id}")
        position_map[original_pos] = len(position_map) + 1
        selected_ids.append(global_cell_id)
        selected_rows.append(merged_row)

    if not selected_rows:
        raise SystemExit("No selected rows remain after matching back to merged input.")

    write_lines_gz(output_dir / "barcodes.tsv.gz", selected_ids)
    write_lines_gz(output_dir / "features.tsv.gz", read_lines_gz(features_path))

    with (output_dir / "merged_metadata.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=merged_fields)
        writer.writeheader()
        writer.writerows(selected_rows)

    entries_path = output_dir / "matrix.entries.tmp"
    total_nnz = 0
    with entries_path.open("w", encoding="utf-8") as entries_handle:
        for entry in iter_matrix_market_entries(matrix_path):
            row_str, col_str, value_str = entry.split()
            new_col = position_map.get(int(col_str))
            if new_col is None:
                continue
            entries_handle.write(f"{row_str} {new_col} {value_str}\n")
            total_nnz += 1

    n_rows, _, _ = parse_matrix_market_shape(matrix_path)
    with (output_dir / "matrix.mtx").open("w", encoding="utf-8") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% filtered integration matrix\n")
        handle.write(f"{n_rows} {len(selected_rows)} {total_nnz}\n")
        with entries_path.open("r", encoding="utf-8") as entries_handle:
            shutil.copyfileobj(entries_handle, handle)
    entries_path.unlink()

    summary = {
        "input_dir": str(input_dir),
        "annotated_metadata": str(annotated_metadata_path),
        "output_dir": str(output_dir),
        "filter_column": args.filter_column,
        "exclude_values": sorted(exclude_values),
        "exclude_clusters": sorted(exclude_clusters),
        "exclude_doublet_class": args.exclude_doublet_class,
        "exclude_doublet_score_ge": args.exclude_doublet_score_ge,
        "exclude_cluster_doublet_score_ge": {
            key: cluster_score_rules[key] for key in sorted(cluster_score_rules)
        },
        "input_cells": len(merged_barcodes),
        "selected_cells": len(selected_rows),
        "removed_cells": len(merged_barcodes) - len(selected_rows),
        "kept_value_counts": dict(sorted(kept_counter.items())),
        "excluded_value_counts": dict(sorted(excluded_counter.items())),
        "removal_reason_counts": dict(sorted(removal_reason_counts.items())),
        "matrix_rows": n_rows,
        "matrix_cols": len(selected_rows),
        "matrix_nnz": total_nnz,
    }
    (output_dir / "filter_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    print(f"[filter] kept cells: {len(selected_rows)}")
    print(f"[filter] removed cells: {len(merged_barcodes) - len(selected_rows)}")
    print(f"[filter] output: {output_dir}")


if __name__ == "__main__":
    main()
