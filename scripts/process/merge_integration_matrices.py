#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


ROOT = Path("/Users/kanae/Works/ML2026")
OUTPUT_ROOT = ROOT / "output"


@dataclass(frozen=True)
class SampleInput:
    gse: str
    gsm: str
    sample_dir: Path
    integration_dir: Path
    matrix_file: Path
    barcodes_file: Path
    features_file: Path
    metadata_file: Path


def read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def write_lines_gz(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for value in values:
            handle.write(f"{value}\n")


def discover_inputs(output_root: Path) -> list[SampleInput]:
    samples: list[SampleInput] = []
    for gse_dir in sorted(path for path in output_root.iterdir() if path.is_dir() and path.name.startswith("GSE")):
        for sample_dir in sorted(path for path in gse_dir.iterdir() if path.is_dir() and path.name.startswith("GSM")):
            integration_dir = sample_dir / "integration_qc"
            matrix_dir = integration_dir / "matrix"
            metadata_file = integration_dir / "metadata_integration_qc.csv"
            matrix_file = matrix_dir / "matrix.mtx"
            barcodes_file = matrix_dir / "barcodes.tsv.gz"
            features_file = matrix_dir / "features.tsv.gz"
            if not all(path.exists() for path in [metadata_file, matrix_file, barcodes_file, features_file]):
                continue
            samples.append(
                SampleInput(
                    gse=gse_dir.name,
                    gsm=sample_dir.name,
                    sample_dir=sample_dir,
                    integration_dir=integration_dir,
                    matrix_file=matrix_file,
                    barcodes_file=barcodes_file,
                    features_file=features_file,
                    metadata_file=metadata_file,
                )
            )
    return samples


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
            return tuple(int(x) for x in parts)
    raise ValueError(f"Missing matrix shape in {path}")


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


def verify_feature_alignment(samples: list[SampleInput]) -> list[str]:
    if not samples:
        raise ValueError("No integration_qc inputs found")
    reference_features = read_lines_gz(samples[0].features_file)
    reference_count = len(reference_features)
    for sample in samples[1:]:
        current_features = read_lines_gz(sample.features_file)
        if len(current_features) != reference_count:
            raise ValueError(
                f"Feature count mismatch: {sample.gse}/{sample.gsm} has {len(current_features)} "
                f"vs {reference_count} in {samples[0].gse}/{samples[0].gsm}"
            )
        if current_features != reference_features:
            for idx, (expected, observed) in enumerate(zip(reference_features, current_features), start=1):
                if expected != observed:
                    raise ValueError(
                        f"Feature order/content mismatch at line {idx}: "
                        f"{sample.gse}/{sample.gsm} has {observed}, expected {expected}"
                    )
            raise ValueError(f"Feature mismatch detected for {sample.gse}/{sample.gsm}")
    return reference_features


def merge_metadata(samples: list[SampleInput], out_dir: Path) -> tuple[pd.DataFrame, list[dict]]:
    merged_frames: list[pd.DataFrame] = []
    sample_manifest: list[dict] = []
    column_offset = 0
    for sample in samples:
        barcodes = read_lines_gz(sample.barcodes_file)
        metadata = pd.read_csv(sample.metadata_file)
        if "cell_barcode" not in metadata.columns:
            raise ValueError(f"cell_barcode column missing in {sample.metadata_file}")
        barcode_order = {barcode: idx for idx, barcode in enumerate(barcodes)}
        metadata = metadata[metadata["cell_barcode"].isin(barcode_order)].copy()
        metadata["__barcode_order"] = metadata["cell_barcode"].map(barcode_order)
        metadata = metadata.sort_values("__barcode_order").drop(columns="__barcode_order")
        if len(metadata) != len(barcodes):
            raise ValueError(
                f"Barcode/metadata mismatch in {sample.gse}/{sample.gsm}: "
                f"{len(barcodes)} matrix barcodes vs {len(metadata)} metadata rows"
            )

        metadata["source_gse"] = sample.gse
        metadata["source_gsm"] = sample.gsm
        metadata["global_cell_id"] = [f"{sample.gsm}:{barcode}" for barcode in barcodes]
        merged_frames.append(metadata)

        sample_manifest.append(
            {
                "gse": sample.gse,
                "gsm": sample.gsm,
                "sample_dir": str(sample.sample_dir),
                "integration_dir": str(sample.integration_dir),
                "matrix_file": str(sample.matrix_file),
                "metadata_file": str(sample.metadata_file),
                "cell_count": len(barcodes),
                "column_start_1based": column_offset + 1,
                "column_end_1based": column_offset + len(barcodes),
            }
        )
        column_offset += len(barcodes)

    merged_metadata = pd.concat(merged_frames, ignore_index=True)
    merged_metadata.to_csv(out_dir / "merged_metadata.csv", index=False)
    with (out_dir / "sample_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_manifest[0].keys()))
        writer.writeheader()
        writer.writerows(sample_manifest)

    write_lines_gz(out_dir / "barcodes.tsv.gz", merged_metadata["global_cell_id"].tolist())
    return merged_metadata, sample_manifest


def merge_matrix_market(samples: list[SampleInput], features: list[str], out_dir: Path) -> dict:
    shapes = [parse_matrix_market_shape(sample.matrix_file) for sample in samples]
    n_rows = shapes[0][0]
    if any(shape[0] != n_rows for shape in shapes):
        raise ValueError("Matrix row counts are not identical across samples")

    total_cols = sum(shape[1] for shape in shapes)
    total_nnz = sum(shape[2] for shape in shapes)
    matrix_path = out_dir / "matrix.mtx"
    tmp_matrix_path = out_dir / "matrix.mtx.tmp"
    tmp_features_path = out_dir / "features.tsv.gz.tmp"

    with tmp_matrix_path.open("w", encoding="utf-8") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% merged integration_qc peak-by-cell matrix\n")
        handle.write(f"{n_rows} {total_cols} {total_nnz}\n")

        col_offset = 0
        for idx, (sample, shape) in enumerate(zip(samples, shapes), start=1):
            expected_cols = len(read_lines_gz(sample.barcodes_file))
            if shape[1] != expected_cols:
                raise ValueError(
                    f"Matrix/barcode count mismatch for {sample.gse}/{sample.gsm}: "
                    f"matrix has {shape[1]} columns, barcodes file has {expected_cols}"
                )
            print(
                f"[merge] matrix {idx}/{len(samples)} {sample.gse}/{sample.gsm} "
                f"cols={shape[1]} offset={col_offset}",
                flush=True,
            )
            for entry in iter_matrix_market_entries(sample.matrix_file):
                row_str, col_str, value_str = entry.split()
                merged_col = int(col_str) + col_offset
                handle.write(f"{row_str} {merged_col} {value_str}\n")
            col_offset += shape[1]

    write_lines_gz(tmp_features_path, features)
    tmp_matrix_path.replace(matrix_path)
    tmp_features_path.replace(out_dir / "features.tsv.gz")
    return {"n_rows": n_rows, "n_cols": total_cols, "n_nnz": total_nnz, "matrix_file": str(matrix_path)}


def main() -> None:
    parser = argparse.ArgumentParser(description="Strictly verify and merge integration_qc peak-by-cell matrices.")
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--merged-dir", default=str(OUTPUT_ROOT / "integration_merged"))
    args = parser.parse_args()

    output_root = Path(args.output_root)
    merged_dir = Path(args.merged_dir)
    merged_dir.mkdir(parents=True, exist_ok=True)

    samples = discover_inputs(output_root)
    if not samples:
        raise SystemExit("No integration_qc inputs found under output/")

    features = verify_feature_alignment(samples)
    merged_metadata, sample_manifest = merge_metadata(samples, merged_dir)
    matrix_info = merge_matrix_market(samples, features, merged_dir)

    summary = {
        "samples": len(samples),
        "features": len(features),
        "cells": int(len(merged_metadata)),
        "matrix_rows": matrix_info["n_rows"],
        "matrix_cols": matrix_info["n_cols"],
        "matrix_nnz": matrix_info["n_nnz"],
        "merged_dir": str(merged_dir),
    }
    (merged_dir / "merge_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    print(f"[merge] samples: {summary['samples']}")
    print(f"[merge] features: {summary['features']}")
    print(f"[merge] cells: {summary['cells']}")
    print(f"[merge] matrix: {matrix_info['matrix_file']}")
    print(f"[merge] metadata: {merged_dir / 'merged_metadata.csv'}")
    print(f"[merge] manifest: {merged_dir / 'sample_manifest.csv'}")


if __name__ == "__main__":
    main()
