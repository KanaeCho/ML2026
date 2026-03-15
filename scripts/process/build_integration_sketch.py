#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
import shutil
from dataclasses import dataclass
from pathlib import Path


ROOT = Path("/Users/kanae/Works/ML2026")
OUTPUT_ROOT = ROOT / "output"


@dataclass(frozen=True)
class SampleInput:
    gse: str
    gsm: str
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
    for gse_dir in sorted(
        path
        for path in output_root.iterdir()
        if path.is_dir() and path.name.startswith("GSE")
    ):
        for sample_dir in sorted(
            path
            for path in gse_dir.iterdir()
            if path.is_dir() and path.name.startswith("GSM")
        ):
            integration_dir = sample_dir / "integration_qc"
            matrix_dir = integration_dir / "matrix"
            metadata_file = integration_dir / "metadata_integration_qc.csv"
            matrix_file = matrix_dir / "matrix.mtx"
            barcodes_file = matrix_dir / "barcodes.tsv.gz"
            features_file = matrix_dir / "features.tsv.gz"
            if not all(
                path.exists()
                for path in [metadata_file, matrix_file, barcodes_file, features_file]
            ):
                continue
            samples.append(
                SampleInput(
                    gse=gse_dir.name,
                    gsm=sample_dir.name,
                    matrix_file=matrix_file,
                    barcodes_file=barcodes_file,
                    features_file=features_file,
                    metadata_file=metadata_file,
                )
            )
    return samples


def verify_features(samples: list[SampleInput]) -> list[str]:
    reference = read_lines_gz(samples[0].features_file)
    for sample in samples[1:]:
        current = read_lines_gz(sample.features_file)
        if current != reference:
            raise ValueError(f"Feature mismatch in {sample.gse}/{sample.gsm}")
    return reference


def parse_shape(path: Path) -> tuple[int, int, int]:
    with path.open("r", encoding="utf-8") as handle:
        _ = handle.readline()
        for line in handle:
            if line.startswith("%"):
                continue
            n_rows, n_cols, n_nnz = (int(x) for x in line.strip().split())
            return n_rows, n_cols, n_nnz
    raise ValueError(f"Invalid matrix file: {path}")


def iter_entries(path: Path):
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


def read_metadata_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = [dict(row) for row in reader]
    return fieldnames, rows


def build_sample_selection(
    sample: SampleInput,
    cells_per_sample: int,
    rng: random.Random,
) -> tuple[list[str], list[dict[str, str]], dict[int, int], dict]:
    barcodes = read_lines_gz(sample.barcodes_file)
    metadata_fields, metadata_rows = read_metadata_rows(sample.metadata_file)
    order = {barcode: idx for idx, barcode in enumerate(barcodes)}
    filtered_rows = [
        row for row in metadata_rows if row.get("cell_barcode", "") in order
    ]
    filtered_rows.sort(key=lambda row: order[row["cell_barcode"]])
    if len(filtered_rows) != len(barcodes):
        raise ValueError(f"Metadata/barcode mismatch in {sample.gse}/{sample.gsm}")

    sample_size = min(cells_per_sample, len(barcodes))
    selected_positions = sorted(rng.sample(range(len(barcodes)), sample_size))
    selected_barcodes = [barcodes[idx] for idx in selected_positions]
    extra_fields = ["source_gse", "source_gsm", "global_cell_id"]
    selected_metadata_rows: list[dict[str, str]] = []
    for barcode, position in zip(selected_barcodes, selected_positions):
        row = dict(filtered_rows[position])
        row["source_gse"] = sample.gse
        row["source_gsm"] = sample.gsm
        row["global_cell_id"] = f"{sample.gsm}:{barcode}"
        selected_metadata_rows.append(row)

    position_map = {
        orig_idx + 1: new_idx
        for new_idx, orig_idx in enumerate(selected_positions, start=1)
    }
    manifest = {
        "gse": sample.gse,
        "gsm": sample.gsm,
        "available_cells": len(barcodes),
        "selected_cells": sample_size,
    }
    return (
        metadata_fields + extra_fields,
        selected_metadata_rows,
        position_map,
        manifest,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build a balanced integration sketch matrix from per-sample integration_qc outputs."
    )
    parser.add_argument("--output-root", default=str(OUTPUT_ROOT))
    parser.add_argument("--sketch-dir", default=str(OUTPUT_ROOT / "integration_sketch"))
    parser.add_argument("--cells-per-sample", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1234)
    args = parser.parse_args()

    output_root = Path(args.output_root)
    sketch_dir = Path(args.sketch_dir)
    sketch_dir.mkdir(parents=True, exist_ok=True)

    samples = discover_inputs(output_root)
    if not samples:
        raise SystemExit("No integration_qc inputs found.")

    features = verify_features(samples)
    rng = random.Random(args.seed)

    metadata_fieldnames: list[str] = []
    selected_rows: list[dict[str, str]] = []
    sample_manifest = []
    position_maps: list[dict[int, int]] = []
    total_cols = 0

    for sample in samples:
        selected_fieldnames, sample_rows, position_map, manifest = (
            build_sample_selection(sample, args.cells_per_sample, rng)
        )
        for field in selected_fieldnames:
            if field not in metadata_fieldnames:
                metadata_fieldnames.append(field)
        selected_rows.extend(sample_rows)
        sample_manifest.append(manifest)
        position_maps.append(position_map)
        total_cols += len(position_map)
        print(
            f"[sketch] {sample.gse}/{sample.gsm}: {manifest['selected_cells']}/{manifest['available_cells']} cells",
            flush=True,
        )

    with (sketch_dir / "merged_metadata.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=metadata_fieldnames)
        writer.writeheader()
        writer.writerows(selected_rows)
    write_lines_gz(
        sketch_dir / "barcodes.tsv.gz", [row["global_cell_id"] for row in selected_rows]
    )
    write_lines_gz(sketch_dir / "features.tsv.gz", features)

    entries_path = sketch_dir / "matrix.entries.tmp"
    total_nnz = 0
    with entries_path.open("w", encoding="utf-8") as entries_handle:
        global_offset = 0
        for idx, (sample, position_map) in enumerate(
            zip(samples, position_maps), start=1
        ):
            sample_nnz = 0
            print(
                f"[sketch] matrix {idx}/{len(samples)} {sample.gse}/{sample.gsm}",
                flush=True,
            )
            for entry in iter_entries(sample.matrix_file):
                row_str, col_str, value_str = entry.split()
                original_col = int(col_str)
                new_local_col = position_map.get(original_col)
                if new_local_col is None:
                    continue
                merged_col = global_offset + new_local_col
                entries_handle.write(f"{row_str} {merged_col} {value_str}\n")
                sample_nnz += 1
            sample_manifest[idx - 1]["selected_nnz"] = sample_nnz
            total_nnz += sample_nnz
            global_offset += len(position_map)

    matrix_path = sketch_dir / "matrix.mtx"
    with matrix_path.open("w", encoding="utf-8") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% integration sketch matrix\n")
        handle.write(f"{len(features)} {total_cols} {total_nnz}\n")
        with entries_path.open("r", encoding="utf-8") as entries_handle:
            shutil.copyfileobj(entries_handle, handle)
    entries_path.unlink()

    with (sketch_dir / "sample_manifest.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_manifest[0].keys()))
        writer.writeheader()
        writer.writerows(sample_manifest)

    summary = {
        "samples": len(samples),
        "features": len(features),
        "cells_per_sample": args.cells_per_sample,
        "total_cells": total_cols,
        "total_nnz": total_nnz,
        "sketch_dir": str(sketch_dir),
    }
    (sketch_dir / "sketch_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )

    print(f"[sketch] total cells: {total_cols}")
    print(f"[sketch] matrix: {matrix_path}")


if __name__ == "__main__":
    main()
