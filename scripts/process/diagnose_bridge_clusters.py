#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Mapping


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_csv_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        return list(reader.fieldnames or []), [dict(row) for row in reader]


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_csv_gz(
    path: Path, fieldnames: list[str], rows: list[dict[str, object]]
) -> None:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def safe_float(value: str | float | int | None) -> float | None:
    if value is None:
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(parsed):
        return None
    return parsed


def parse_float_list(value: str) -> list[float]:
    out: list[float] = []
    for piece in value.split(","):
        piece = piece.strip()
        if not piece:
            continue
        out.append(float(piece))
    return out


def quantile(values: list[float], q: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    idx = (len(ordered) - 1) * q
    lo = int(math.floor(idx))
    hi = int(math.ceil(idx))
    if lo == hi:
        return ordered[lo]
    frac = idx - lo
    return ordered[lo] * (1 - frac) + ordered[hi] * frac


def format_float(value: float | None, ndigits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{ndigits}f}"


def format_pct(value: float | None) -> str:
    if value is None:
        return ""
    return f"{value * 100:.2f}%"


def build_validation_lookup(
    validation_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    return {
        row["seurat_clusters"]: row
        for row in validation_rows
        if row.get("seurat_clusters")
    }


def build_celltype_lookup(
    celltyped_rows: list[dict[str, str]],
) -> dict[str, dict[str, str]]:
    return {
        row["global_cell_id"]: row
        for row in celltyped_rows
        if row.get("global_cell_id")
    }


def suspicion_flags(
    cluster_row: Mapping[str, Any], candidate_threshold: float
) -> list[str]:
    flags: list[str] = []
    frac_ge_candidate = float(
        cluster_row[f"doublet_score_frac_ge_{candidate_threshold:g}"]
    )
    top_gse_fraction = safe_float(cluster_row.get("top_gse_fraction")) or 0.0
    top_gsm_fraction = safe_float(cluster_row.get("top_gsm_fraction")) or 0.0
    review_priority = str(cluster_row.get("review_priority", ""))
    annotation_confidence = str(cluster_row.get("annotation_confidence", ""))
    review_reasons = str(cluster_row.get("review_reasons", ""))
    q95 = safe_float(cluster_row.get("doublet_score_q95")) or 0.0

    if frac_ge_candidate >= 0.20:
        flags.append("doublet_enriched_cluster")
    if q95 >= max(candidate_threshold + 0.1, 0.6):
        flags.append("high_doublet_tail")
    if top_gsm_fraction >= 0.35:
        flags.append("sample_dominant")
    if top_gse_fraction >= 0.80:
        flags.append("gse_dominant")
    if annotation_confidence == "low":
        flags.append("low_annotation_confidence")
    if review_priority == "high":
        flags.append("high_review_priority")
    if "broad_label_mixed" in review_reasons:
        flags.append("broad_label_mixed")
    if "subtype_mixed" in review_reasons:
        flags.append("subtype_mixed")
    return flags


def suspicion_score(flags: list[str]) -> int:
    weights = {
        "doublet_enriched_cluster": 4,
        "high_doublet_tail": 2,
        "sample_dominant": 3,
        "gse_dominant": 2,
        "low_annotation_confidence": 2,
        "high_review_priority": 2,
        "broad_label_mixed": 2,
        "subtype_mixed": 1,
    }
    return sum(weights.get(flag, 1) for flag in flags)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Diagnose bridge-like scATAC clusters using QC, residual doublet score, and sample dominance."
    )
    parser.add_argument(
        "--analysis-dir",
        default=str(OUTPUT_ROOT / "integration_merged_without_unknown_analysis"),
    )
    parser.add_argument(
        "--validation-summary",
        default=None,
        help="Optional cluster_annotation_validation_summary.csv path",
    )
    parser.add_argument(
        "--celltyped-metadata",
        default=None,
        help="Optional integrated_metadata_celltyped*.csv.gz path for celltype labels in candidate-cell output",
    )
    parser.add_argument(
        "--doublet-score-thresholds",
        default="0.4,0.5,0.6",
        help="Comma-separated doublet-score thresholds to summarize",
    )
    parser.add_argument(
        "--candidate-score-threshold",
        type=float,
        default=0.5,
        help="Threshold for candidate doublet-enriched cells output",
    )
    parser.add_argument(
        "--top-clusters",
        type=int,
        default=10,
        help="How many top suspicious clusters to highlight in the markdown report",
    )
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    metadata_path = analysis_dir / "integrated_metadata.csv.gz"
    if not metadata_path.exists():
        raise SystemExit(f"Missing integrated metadata: {metadata_path}")

    validation_path = (
        Path(args.validation_summary)
        if args.validation_summary
        else analysis_dir / "cluster_annotation_validation_summary.csv"
    )
    celltyped_path = (
        Path(args.celltyped_metadata)
        if args.celltyped_metadata
        else analysis_dir / "integrated_metadata_celltyped.csv.gz"
    )

    _, metadata_rows = read_csv_rows(metadata_path)
    validation_lookup: dict[str, dict[str, str]] = {}
    if validation_path.exists():
        _, validation_rows = read_csv_rows(validation_path)
        validation_lookup = build_validation_lookup(validation_rows)

    celltype_lookup: dict[str, dict[str, str]] = {}
    if celltyped_path.exists():
        _, celltyped_rows = read_csv_rows(celltyped_path)
        celltype_lookup = build_celltype_lookup(celltyped_rows)

    thresholds = sorted(set(parse_float_list(args.doublet_score_thresholds)))
    if args.candidate_score_threshold not in thresholds:
        thresholds.append(args.candidate_score_threshold)
        thresholds = sorted(thresholds)

    cluster_acc: defaultdict[str, dict[str, Any]] = defaultdict(
        lambda: {
            "rows": 0,
            "doublet_scores": [],
            "nCount_ATAC": [],
            "FRiP": [],
            "TSS.enrichment": [],
            "nucleosome_signal": [],
            "mixing_source_gse": [],
            "mixing_source_gsm": [],
            "source_gse": Counter(),
            "source_gsm": Counter(),
            "threshold_counts": Counter(),
        }
    )
    candidate_cells: list[dict[str, Any]] = []

    for row in metadata_rows:
        cluster = row["seurat_clusters"]
        acc = cluster_acc[cluster]
        acc["rows"] += 1
        score = safe_float(row.get("scDblFinder.score"))
        if score is not None:
            acc["doublet_scores"].append(score)
            for threshold in thresholds:
                if score >= threshold:
                    acc["threshold_counts"][threshold] += 1
        for key in [
            "nCount_ATAC",
            "FRiP",
            "TSS.enrichment",
            "nucleosome_signal",
            "mixing_source_gse",
            "mixing_source_gsm",
        ]:
            value = safe_float(row.get(key))
            if value is not None:
                acc[key].append(value)
        acc["source_gse"][row.get("source_gse", "")] += 1
        acc["source_gsm"][row.get("source_gsm", "")] += 1

        if score is not None and score >= args.candidate_score_threshold:
            celltyped = celltype_lookup.get(row.get("global_cell_id", ""), {})
            candidate_cells.append(
                {
                    "global_cell_id": row.get("global_cell_id", ""),
                    "seurat_clusters": cluster,
                    "scDblFinder.score": format_float(score, 6),
                    "source_gse": row.get("source_gse", ""),
                    "source_gsm": row.get("source_gsm", ""),
                    "nCount_ATAC": row.get("nCount_ATAC", ""),
                    "FRiP": row.get("FRiP", ""),
                    "TSS.enrichment": row.get("TSS.enrichment", ""),
                    "nucleosome_signal": row.get("nucleosome_signal", ""),
                    "celltype": celltyped.get("celltype", ""),
                    "celltype_subtype": celltyped.get("celltype_subtype", ""),
                }
            )

    cluster_rows: list[dict[str, object]] = []
    for cluster, acc in sorted(cluster_acc.items(), key=lambda item: int(item[0])):
        total = acc["rows"]
        top_gse, top_gse_count = acc["source_gse"].most_common(1)[0]
        top_gsm, top_gsm_count = acc["source_gsm"].most_common(1)[0]
        validation = validation_lookup.get(cluster, {})
        row: dict[str, Any] = {
            "seurat_clusters": cluster,
            "total_cells": total,
            "celltype": validation.get("celltype", ""),
            "celltype_subtype": validation.get("celltype_subtype", ""),
            "annotation_confidence": validation.get("annotation_confidence", ""),
            "review_priority": validation.get("review_priority", ""),
            "review_reasons": validation.get("review_reasons", ""),
            "top_marker_genes": validation.get("top_marker_genes", ""),
            "top_reference_celltype": validation.get("top_reference_celltype", ""),
            "top_reference_celltype_fraction": validation.get(
                "top_reference_celltype_fraction", ""
            ),
            "doublet_score_mean": format_float(
                sum(acc["doublet_scores"]) / len(acc["doublet_scores"])
                if acc["doublet_scores"]
                else None,
                6,
            ),
            "doublet_score_median": format_float(
                median(acc["doublet_scores"]) if acc["doublet_scores"] else None, 6
            ),
            "doublet_score_q95": format_float(quantile(acc["doublet_scores"], 0.95), 6),
            "nCount_ATAC_median": format_float(
                median(acc["nCount_ATAC"]) if acc["nCount_ATAC"] else None, 2
            ),
            "FRiP_median": format_float(
                median(acc["FRiP"]) if acc["FRiP"] else None, 6
            ),
            "TSS_enrichment_median": format_float(
                median(acc["TSS.enrichment"]) if acc["TSS.enrichment"] else None, 6
            ),
            "nucleosome_signal_median": format_float(
                median(acc["nucleosome_signal"]) if acc["nucleosome_signal"] else None,
                6,
            ),
            "mixing_source_gse_median": format_float(
                median(acc["mixing_source_gse"]) if acc["mixing_source_gse"] else None,
                6,
            ),
            "mixing_source_gsm_median": format_float(
                median(acc["mixing_source_gsm"]) if acc["mixing_source_gsm"] else None,
                6,
            ),
            "top_gse": top_gse,
            "top_gse_fraction": format_float(top_gse_count / total, 4),
            "top_gsm": top_gsm,
            "top_gsm_fraction": format_float(top_gsm_count / total, 4),
        }
        for threshold in thresholds:
            row[f"doublet_score_count_ge_{threshold:g}"] = acc["threshold_counts"][
                threshold
            ]
            row[f"doublet_score_frac_ge_{threshold:g}"] = format_float(
                acc["threshold_counts"][threshold] / total,
                4,
            )
        flags = suspicion_flags(row, args.candidate_score_threshold)
        row["suspicion_flags"] = "; ".join(flags)
        row["suspicion_score"] = suspicion_score(flags)
        cluster_rows.append(row)

    cluster_rows.sort(
        key=lambda row: (
            int(str(row["suspicion_score"])),
            safe_float(
                str(
                    row.get(
                        f"doublet_score_frac_ge_{args.candidate_score_threshold:g}", ""
                    )
                )
            )
            or 0.0,
            safe_float(str(row.get("top_gsm_fraction", ""))) or 0.0,
            safe_float(str(row.get("top_gse_fraction", ""))) or 0.0,
            -int(str(row["seurat_clusters"])),
        ),
        reverse=True,
    )

    cluster_summary_path = analysis_dir / "bridge_cluster_diagnostic_summary.csv"
    candidate_cells_path = analysis_dir / "bridge_candidate_doublet_cells.csv.gz"
    report_path = analysis_dir / "bridge_diagnostic_report.md"
    summary_json_path = analysis_dir / "bridge_diagnostic_summary.json"

    cluster_fieldnames = list(cluster_rows[0].keys()) if cluster_rows else []
    if cluster_rows:
        write_csv(cluster_summary_path, cluster_fieldnames, cluster_rows)

    candidate_cells.sort(
        key=lambda row: (
            int(row["seurat_clusters"]),
            -float(row["scDblFinder.score"]),
            row["global_cell_id"],
        )
    )
    candidate_fieldnames = (
        list(candidate_cells[0].keys())
        if candidate_cells
        else [
            "global_cell_id",
            "seurat_clusters",
            "scDblFinder.score",
            "source_gse",
            "source_gsm",
            "nCount_ATAC",
            "FRiP",
            "TSS.enrichment",
            "nucleosome_signal",
            "celltype",
            "celltype_subtype",
        ]
    )
    write_csv_gz(candidate_cells_path, candidate_fieldnames, candidate_cells)

    top_rows = cluster_rows[: args.top_clusters]
    report_lines = [
        "# Bridge cluster diagnostic report",
        "",
        f"- Analysis directory: `{analysis_dir}`",
        f"- Cells scanned: {len(metadata_rows):,}",
        f"- Candidate score threshold: `{args.candidate_score_threshold}`",
        f"- Doublet summary thresholds: `{', '.join(str(x) for x in thresholds)}`",
        f"- Candidate cells meeting score threshold: {len(candidate_cells):,}",
        "",
        "## Top suspicious clusters",
        "",
    ]
    for row in top_rows:
        report_lines.extend(
            [
                f"### Cluster {row['seurat_clusters']} — suspicion score {row['suspicion_score']}",
                f"- Broad label: `{row.get('celltype', '')}` / `{row.get('celltype_subtype', '')}`",
                f"- Review priority: `{row.get('review_priority', '')}`",
                f"- Doublet burden >= {args.candidate_score_threshold}: `{format_pct(safe_float(row.get(f'doublet_score_frac_ge_{args.candidate_score_threshold:g}')))} ` ({row.get(f'doublet_score_count_ge_{args.candidate_score_threshold:g}', 0)} cells)",
                f"- Doublet score mean / q95: `{row.get('doublet_score_mean', '')}` / `{row.get('doublet_score_q95', '')}`",
                f"- Median QC: `nCount={row.get('nCount_ATAC_median', '')}`, `FRiP={row.get('FRiP_median', '')}`, `TSS={row.get('TSS_enrichment_median', '')}`, `nucleosome={row.get('nucleosome_signal_median', '')}`",
                f"- Top batch: `GSE={row.get('top_gse', '')}` ({format_pct(safe_float(row.get('top_gse_fraction')))}), `GSM={row.get('top_gsm', '')}` ({format_pct(safe_float(row.get('top_gsm_fraction')))} )",
                f"- Flags: `{row.get('suspicion_flags', '')}`",
                f"- Review reasons: `{row.get('review_reasons', '')}`",
                f"- Top markers: `{row.get('top_marker_genes', '')}`",
                "",
            ]
        )
    report_lines.extend(
        [
            "## Suggested next actions",
            "",
            "1. Start with cluster-directed sensitivity reruns rather than global QC retuning.",
            "2. Use `bridge_candidate_doublet_cells.csv.gz` to remove high-score cells from the top suspicious cluster(s).",
            "3. Compare current dims against LSI-cleaned reruns after excluding QC-loaded dimensions.",
            "4. Treat low-confidence mixed clusters as unresolved until structural cleanup is complete.",
            "",
            "## Output files",
            "",
            f"- `{cluster_summary_path.name}`",
            f"- `{candidate_cells_path.name}`",
            f"- `{summary_json_path.name}`",
        ]
    )
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    summary_payload = {
        "analysis_dir": str(analysis_dir),
        "cells_scanned": len(metadata_rows),
        "candidate_score_threshold": args.candidate_score_threshold,
        "doublet_score_thresholds": thresholds,
        "candidate_cell_count": len(candidate_cells),
        "top_suspicious_clusters": [row["seurat_clusters"] for row in top_rows],
    }
    summary_json_path.write_text(
        json.dumps(summary_payload, indent=2) + "\n", encoding="utf-8"
    )

    print(f"[bridge-diagnose] wrote {cluster_summary_path}")
    print(f"[bridge-diagnose] wrote {candidate_cells_path}")
    print(f"[bridge-diagnose] wrote {report_path}")


if __name__ == "__main__":
    main()
