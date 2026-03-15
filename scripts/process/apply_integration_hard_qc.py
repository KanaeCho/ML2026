#!/usr/bin/env python3

import argparse
import csv
import gzip
import shutil
from pathlib import Path

import pandas as pd
from scipy.io import mmread, mmwrite


def read_lines_gz(path: Path) -> list[str]:
    with gzip.open(path, "rt") as handle:
        return [line.rstrip("\n") for line in handle]


def write_lines_gz(path: Path, values: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        for value in values:
            handle.write(f"{value}\n")


def apply_hard_filters(
    metadata_qc: pd.DataFrame,
    min_count: int,
    max_count: int,
    min_tss: float,
    min_frip: float,
    max_blacklist: float,
    max_nucleosome: float,
    min_unique_ratio: float | None,
) -> pd.Series:
    keep = (
        metadata_qc["nCount_ATAC"].between(min_count, max_count)
        & (metadata_qc["TSS.enrichment"] >= min_tss)
        & (metadata_qc["FRiP"] >= min_frip)
        & (metadata_qc["blacklist_fraction"] <= max_blacklist)
        & (metadata_qc["nucleosome_signal"] <= max_nucleosome)
    )
    if min_unique_ratio is not None and "unique_ratio" in metadata_qc.columns:
        unique_ratio = pd.to_numeric(metadata_qc["unique_ratio"], errors="coerce")
        keep &= unique_ratio.isna() | (unique_ratio >= min_unique_ratio)
    return keep


def process_sample(sample_dir: Path, args: argparse.Namespace) -> dict:
    metadata_qc_file = sample_dir / "metadata_qc.csv"
    matrix_dir = sample_dir / "matrix"
    matrix_file = matrix_dir / "matrix.mtx"
    barcodes_file = matrix_dir / "barcodes.tsv.gz"
    features_file = matrix_dir / "features.tsv.gz"

    if not metadata_qc_file.exists():
        raise FileNotFoundError(f"Missing metadata_qc.csv under {sample_dir}")
    if not matrix_file.exists():
        raise FileNotFoundError(f"Missing matrix.mtx under {sample_dir}")
    if not barcodes_file.exists():
        raise FileNotFoundError(f"Missing barcodes.tsv.gz under {sample_dir}")
    if not features_file.exists():
        raise FileNotFoundError(f"Missing features.tsv.gz under {sample_dir}")

    metadata_qc = pd.read_csv(metadata_qc_file)
    metadata_qc["nCount_ATAC"] = pd.to_numeric(metadata_qc["nCount_ATAC"], errors="coerce")
    metadata_qc["TSS.enrichment"] = pd.to_numeric(metadata_qc["TSS.enrichment"], errors="coerce")
    metadata_qc["FRiP"] = pd.to_numeric(metadata_qc["FRiP"], errors="coerce")
    metadata_qc["blacklist_fraction"] = pd.to_numeric(metadata_qc["blacklist_fraction"], errors="coerce")
    metadata_qc["nucleosome_signal"] = pd.to_numeric(metadata_qc["nucleosome_signal"], errors="coerce")
    if "unique_ratio" in metadata_qc.columns:
        metadata_qc["unique_ratio"] = pd.to_numeric(metadata_qc["unique_ratio"], errors="coerce")

    keep_mask = apply_hard_filters(
        metadata_qc=metadata_qc,
        min_count=args.min_count,
        max_count=args.max_count,
        min_tss=args.min_tss,
        min_frip=args.min_frip,
        max_blacklist=args.max_blacklist,
        max_nucleosome=args.max_nucleosome,
        min_unique_ratio=args.min_unique_ratio,
    )
    kept_metadata = metadata_qc.loc[keep_mask].copy()
    keep_barcodes = set(kept_metadata["cell_barcode"].tolist())

    matrix_barcodes = read_lines_gz(barcodes_file)
    keep_positions = [idx for idx, barcode in enumerate(matrix_barcodes) if barcode in keep_barcodes]
    kept_barcodes = [matrix_barcodes[idx] for idx in keep_positions]

    barcode_order = {barcode: idx for idx, barcode in enumerate(kept_barcodes)}
    kept_metadata = kept_metadata[kept_metadata["cell_barcode"].isin(barcode_order)].copy()
    kept_metadata["__barcode_order"] = kept_metadata["cell_barcode"].map(barcode_order)
    kept_metadata = kept_metadata.sort_values("__barcode_order").drop(columns="__barcode_order")

    sparse_matrix = mmread(matrix_file).tocsc()
    filtered_matrix = sparse_matrix[:, keep_positions]

    integration_dir = sample_dir / "integration_qc"
    integration_matrix_dir = integration_dir / "matrix"
    integration_matrix_dir.mkdir(parents=True, exist_ok=True)

    integration_metadata_file = integration_dir / "metadata_integration_qc.csv"
    integration_summary_file = integration_dir / "integration_hard_qc_summary.csv"
    integration_matrix_file = integration_matrix_dir / "matrix.mtx"
    integration_barcodes_file = integration_matrix_dir / "barcodes.tsv.gz"
    integration_features_file = integration_matrix_dir / "features.tsv.gz"

    kept_metadata.to_csv(integration_metadata_file, index=False)
    mmwrite(str(integration_matrix_file), filtered_matrix)
    write_lines_gz(integration_barcodes_file, kept_barcodes)
    shutil.copyfile(features_file, integration_features_file)

    summary_rows = [
        ["gse", sample_dir.parent.name],
        ["gsm", sample_dir.name],
        ["source_metadata_qc_cells", len(metadata_qc)],
        ["integration_qc_cells", len(kept_metadata)],
        ["integration_qc_rate_vs_metadata_qc", f"{(len(kept_metadata) / len(metadata_qc) * 100) if len(metadata_qc) else 0:.2f}%"],
        ["min_nCount_ATAC", args.min_count],
        ["max_nCount_ATAC", args.max_count],
        ["min_TSS_enrichment", args.min_tss],
        ["min_FRiP", args.min_frip],
        ["max_blacklist_fraction", args.max_blacklist],
        ["max_nucleosome_signal", args.max_nucleosome],
        ["min_unique_ratio", "" if args.min_unique_ratio is None else args.min_unique_ratio],
    ]
    with open(integration_summary_file, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(summary_rows)

    return {
        "gse": sample_dir.parent.name,
        "gsm": sample_dir.name,
        "metadata_qc_cells": len(metadata_qc),
        "integration_qc_cells": len(kept_metadata),
        "integration_qc_rate": (len(kept_metadata) / len(metadata_qc) * 100) if len(metadata_qc) else 0.0,
    }


def discover_samples(output_root: Path, gse: str | None) -> list[Path]:
    pattern = str(output_root / (gse if gse else "*") / "GSM*")
    return sorted(Path(p) for p in output_root.glob(f"{gse if gse else '*'}" + "/GSM*"))


def main() -> None:
    parser = argparse.ArgumentParser(description="Apply integration hard QC on existing metadata_qc + matrix outputs.")
    parser.add_argument("--output-root", default="/Users/kanae/Works/ML2026/output")
    parser.add_argument("--gse")
    parser.add_argument("--min-count", type=int, default=1000)
    parser.add_argument("--max-count", type=int, default=100000)
    parser.add_argument("--min-tss", type=float, default=4.0)
    parser.add_argument("--min-frip", type=float, default=0.35)
    parser.add_argument("--max-blacklist", type=float, default=0.05)
    parser.add_argument("--max-nucleosome", type=float, default=4.0)
    parser.add_argument("--min-unique-ratio", type=float, default=None)
    args = parser.parse_args()

    output_root = Path(args.output_root)
    sample_dirs = discover_samples(output_root, args.gse)
    if not sample_dirs:
        raise SystemExit("No sample directories found.")

    results = []
    for sample_dir in sample_dirs:
        result = process_sample(sample_dir, args)
        results.append(result)
        print(
            f"[hard-qc] {result['gse']}/{result['gsm']}: "
            f"{result['integration_qc_cells']}/{result['metadata_qc_cells']} "
            f"({result['integration_qc_rate']:.2f}%)"
        )

    summary_file = output_root / "integration_hard_qc_summary.csv"
    pd.DataFrame(results).to_csv(summary_file, index=False)
    print(f"[hard-qc] summary: {summary_file}")


if __name__ == "__main__":
    main()
