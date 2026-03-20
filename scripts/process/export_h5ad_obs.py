#!/usr/bin/env python3

from __future__ import annotations

import argparse
import gzip
import importlib
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Export .h5ad obs metadata to CSV/CSV.GZ for downstream label mapping."
    )
    parser.add_argument("--input-h5ad", required=True, help="Input .h5ad path")
    parser.add_argument(
        "--output-csv", required=True, help="Output .csv or .csv.gz path"
    )
    parser.add_argument(
        "--obs-name-col",
        default="obs_name",
        help="Column name used for AnnData obs_names",
    )
    parser.add_argument(
        "--columns",
        default="",
        help="Optional comma-separated subset of obs columns to export",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    try:
        ad = importlib.import_module("anndata")
    except Exception as exc:  # pragma: no cover - dependency guard
        sys.stderr.write(
            "Missing optional dependency 'anndata'. Install it before exporting H5AD metadata.\n"
        )
        sys.stderr.write(f"Import error: {exc}\n")
        return 1

    input_path = Path(args.input_h5ad).resolve()
    output_path = Path(args.output_csv).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(input_path, backed="r")
    obs = adata.obs.copy()
    obs.insert(0, args.obs_name_col, adata.obs_names.astype(str))

    requested = [col.strip() for col in args.columns.split(",") if col.strip()]
    if requested:
        if args.obs_name_col not in requested:
            requested = [args.obs_name_col, *requested]
        missing = [col for col in requested if col not in obs.columns]
        if missing:
            raise SystemExit(f"Missing requested obs columns: {', '.join(missing)}")
        obs = obs[requested]

    if output_path.suffix == ".gz":
        with gzip.open(output_path, "wt", encoding="utf-8", newline="") as handle:
            obs.to_csv(handle, index=False)
    else:
        obs.to_csv(output_path, index=False)

    print(f"Exported {len(obs):,} obs rows to {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
