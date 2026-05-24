#!/usr/bin/env python3
"""Write top fragment-count barcodes from a fragments.tsv.gz file."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def count_barcodes(fragment_file: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    with gzip.open(fragment_file, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 4 and parts[3]:
                counts[parts[3]] += 1
    return counts


def write_top_barcodes(fragment_file: Path, output_dir: Path, top_n: int) -> Path:
    counts = count_barcodes(fragment_file)
    selected = counts.most_common(top_n)
    output_dir.mkdir(parents=True, exist_ok=True)
    barcode_file = output_dir / "filtered_barcodes.tsv.gz"
    with gzip.open(barcode_file, "wt", encoding="utf-8") as handle:
        for barcode, _ in selected:
            handle.write(f"{barcode}\n")
    summary = {
        "fragment_file": str(fragment_file),
        "barcode_file": str(barcode_file),
        "unique_barcodes": len(counts),
        "selected_barcodes": len(selected),
        "top_n": top_n,
        "min_selected_fragments": selected[-1][1] if selected else 0,
        "max_selected_fragments": selected[0][1] if selected else 0,
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return barcode_file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fragment-file", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--top-n", type=int, default=8000)
    args = parser.parse_args()
    barcode_file = write_top_barcodes(args.fragment_file, args.output_dir, args.top_n)
    print(f"barcode_file={barcode_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
