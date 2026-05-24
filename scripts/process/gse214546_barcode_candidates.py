#!/usr/bin/env python3
"""Create barcode candidate files from GSE214546 filtered metadata."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RAW_ROOT = ROOT / "data" / "raw" / "GSE214546"
DEFAULT_OUTPUT_ROOT = ROOT / "data" / "reference" / "GSE214546" / "atac_barcodes_param_contrast"


def _clean_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "singlet"})


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(0.0, index=frame.index)
    return pd.to_numeric(frame[column], errors="coerce").fillna(0.0)


def _rank_score(frame: pd.DataFrame) -> pd.Series:
    n_unique = _numeric(frame, "n_unique").rank(pct=True)
    tss = _numeric(frame, "TSSEnrichment").rank(pct=True)
    peaks = _numeric(frame, "peaks_frac").rank(pct=True)
    altius = _numeric(frame, "altius_frac").rank(pct=True)
    return 0.35 * tss + 0.30 * peaks + 0.25 * n_unique + 0.10 * altius


def _singlets(frame: pd.DataFrame) -> pd.DataFrame:
    if "singlet" not in frame.columns:
        return frame.copy()
    return frame.loc[_clean_bool(frame["singlet"])].copy()


def select_candidate(frame: pd.DataFrame, candidate: str) -> tuple[pd.DataFrame, dict[str, object]]:
    selectors: dict[str, Callable[[pd.DataFrame], pd.DataFrame]] = {
        "baseline_filtered_metadata": lambda df: df.copy(),
        "singlet_only": _singlets,
        "singlet_top20000_nunique": lambda df: _singlets(df).assign(
            _score=_numeric(_singlets(df), "n_unique")
        ).sort_values("_score", ascending=False).head(20000),
        "singlet_top16000_nunique": lambda df: _singlets(df).assign(
            _score=_numeric(_singlets(df), "n_unique")
        ).sort_values("_score", ascending=False).head(16000),
        "singlet_top12000_nunique": lambda df: _singlets(df).assign(
            _score=_numeric(_singlets(df), "n_unique")
        ).sort_values("_score", ascending=False).head(12000),
        "singlet_top8000_nunique": lambda df: _singlets(df).assign(
            _score=_numeric(_singlets(df), "n_unique")
        ).sort_values("_score", ascending=False).head(8000),
        "singlet_tss_peaks_top12000": lambda df: _singlets(df).assign(
            _score=_rank_score(_singlets(df))
        ).sort_values("_score", ascending=False).head(12000),
        "singlet_min_unique5000": lambda df: _singlets(df).loc[_numeric(_singlets(df), "n_unique") >= 5000].copy(),
    }
    if candidate not in selectors:
        raise ValueError(f"Unknown candidate: {candidate}")
    selected = selectors[candidate](frame)
    selected = selected.drop(columns=["_score"], errors="ignore")
    summary = {
        "candidate": candidate,
        "input_barcodes": int(len(frame)),
        "selected_barcodes": int(len(selected)),
        "selection_fraction": float(len(selected) / max(len(frame), 1)),
    }
    for column in ["n_unique", "n_fragments", "TSSEnrichment", "peaks_frac", "altius_frac"]:
        if column in selected.columns and len(selected):
            summary[f"median_{column}"] = float(pd.to_numeric(selected[column], errors="coerce").median())
    return selected, summary


def find_filtered_metadata(raw_root: Path, sample_id: str) -> Path:
    matches = sorted(raw_root.glob(f"{sample_id}*filtered_metadata.csv.gz"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one filtered_metadata for {sample_id}, found {len(matches)}")
    return matches[0]


def write_candidate(raw_root: Path, output_root: Path, sample_id: str, candidate: str) -> Path:
    metadata_path = find_filtered_metadata(raw_root, sample_id)
    frame = pd.read_csv(metadata_path, dtype=str, keep_default_na=False)
    if "barcodes" not in frame.columns:
        raise KeyError(f"{metadata_path} missing required 'barcodes' column")
    selected, summary = select_candidate(frame, candidate)
    out_dir = output_root / candidate / sample_id
    out_dir.mkdir(parents=True, exist_ok=True)
    barcode_path = out_dir / "filtered_barcodes.tsv.gz"
    with gzip.open(barcode_path, "wt", encoding="utf-8") as handle:
        for barcode in selected["barcodes"].astype(str):
            if barcode:
                handle.write(f"{barcode}\n")
    summary.update(
        {
            "sample_id": sample_id,
            "source_filtered_metadata": str(metadata_path),
            "barcode_file": str(barcode_path),
        }
    )
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return barcode_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--raw-root", type=Path, default=DEFAULT_RAW_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()
    barcode_path = write_candidate(args.raw_root, args.output_root, args.sample_id, args.candidate)
    print(f"barcode_file={barcode_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
