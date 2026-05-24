#!/usr/bin/env python3
"""Create GSE282769 barcode candidates from existing ATAC metadata."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path
from typing import Callable

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_METADATA_ROOT = ROOT / "output" / "atac" / "GSE282769"
DEFAULT_OUTPUT_ROOT = ROOT / "data" / "reference" / "GSE282769" / "atac_barcodes_param_contrast"


def _numeric(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series(0.0, index=frame.index)
    return pd.to_numeric(frame[column], errors="coerce").fillna(0.0)


def _high_quality(frame: pd.DataFrame, min_frip: float, min_tss: float) -> pd.DataFrame:
    frip = _numeric(frame, "FRiP")
    tss = _numeric(frame, "TSS.enrichment")
    return frame.loc[(frip >= min_frip) & (tss >= min_tss)].copy()


def _ranked(frame: pd.DataFrame) -> pd.DataFrame:
    scored = frame.copy()
    scored["_score"] = (
        0.45 * _numeric(scored, "FRiP").rank(pct=True)
        + 0.35 * _numeric(scored, "TSS.enrichment").rank(pct=True)
        + 0.20 * _numeric(scored, "fragments").rank(pct=True)
    )
    return scored.sort_values("_score", ascending=False)


def _top_quality(frame: pd.DataFrame, min_frip: float, min_tss: float, cap: int) -> pd.DataFrame:
    return _ranked(_high_quality(frame, min_frip=min_frip, min_tss=min_tss)).head(cap)


def _top_ranked(frame: pd.DataFrame, cap: int) -> pd.DataFrame:
    return _ranked(frame).head(cap)


def select_candidate(frame: pd.DataFrame, candidate: str) -> tuple[pd.DataFrame, dict[str, object]]:
    selectors: dict[str, Callable[[pd.DataFrame], pd.DataFrame]] = {
        "frip06_tss2": lambda df: _high_quality(df, min_frip=0.6, min_tss=2.0),
        "frip06_tss2_top8000": lambda df: _top_quality(df, min_frip=0.6, min_tss=2.0, cap=8000),
        "frip05_tss2_top8000": lambda df: _top_quality(df, min_frip=0.5, min_tss=2.0, cap=8000),
        "frip04_tss2_top8000": lambda df: _top_quality(df, min_frip=0.4, min_tss=2.0, cap=8000),
        "ranked_top8000": lambda df: _top_ranked(df, cap=8000),
    }
    if candidate not in selectors:
        raise ValueError(f"Unknown candidate: {candidate}")

    selected = selectors[candidate](frame).drop(columns=["_score"], errors="ignore")
    summary: dict[str, object] = {
        "candidate": candidate,
        "input_barcodes": int(len(frame)),
        "selected_barcodes": int(len(selected)),
        "selection_fraction": float(len(selected) / max(len(frame), 1)),
    }
    for column in ["FRiP", "TSS.enrichment", "fragments", "nCount_ATAC", "nFeature_ATAC", "scDblFinder.score"]:
        if column in selected.columns and len(selected):
            summary[f"median_{column}"] = float(pd.to_numeric(selected[column], errors="coerce").median())
    if "scDblFinder.class" in selected.columns and len(selected):
        summary["baseline_singlet_fraction"] = float(
            selected["scDblFinder.class"].astype(str).str.lower().eq("singlet").mean()
        )
    return selected, summary


def metadata_path(metadata_root: Path, sample_id: str) -> Path:
    path = metadata_root / sample_id / "metadata.csv"
    if not path.exists():
        raise FileNotFoundError(f"Missing metadata.csv for {sample_id}: {path}")
    return path


def write_candidate(metadata_root: Path, output_root: Path, sample_id: str, candidate: str) -> Path:
    source_path = metadata_path(metadata_root, sample_id)
    frame = pd.read_csv(source_path, dtype=str, keep_default_na=False)
    if "cell_barcode" not in frame.columns:
        raise KeyError(f"{source_path} missing required 'cell_barcode' column")
    selected, summary = select_candidate(frame, candidate)

    out_dir = output_root / candidate / sample_id
    out_dir.mkdir(parents=True, exist_ok=True)
    barcode_path = out_dir / "filtered_barcodes.tsv.gz"
    with gzip.open(barcode_path, "wt", encoding="utf-8") as handle:
        for barcode in selected["cell_barcode"].astype(str):
            if barcode:
                handle.write(f"{barcode}\n")

    summary.update(
        {
            "sample_id": sample_id,
            "source_metadata": str(source_path),
            "barcode_file": str(barcode_path),
        }
    )
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return barcode_path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--metadata-root", type=Path, default=DEFAULT_METADATA_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    args = parser.parse_args()

    barcode_path = write_candidate(args.metadata_root, args.output_root, args.sample_id, args.candidate)
    print(f"barcode_file={barcode_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
