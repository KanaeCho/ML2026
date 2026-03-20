#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle))


def normalize_cluster(value: str) -> str:
    return value[1:] if value.startswith("g") else value


def parse_gene_list(value: str) -> list[str]:
    return [item.strip() for item in value.split(";") if item.strip()]


def label_keywords(celltype: str, subtype: str) -> list[str]:
    text = f"{celltype}_{subtype}".replace("/", "_")
    tokens = re.split(r"[^A-Za-z0-9]+", text)
    ignore = {
        "T",
        "NK",
        "B",
        "CD",
        "CELL",
        "LINEAGE",
        "LIKE",
        "UNKNOWN",
        "PBMC",
        "NON",
        "INNATE",
        "MYELOID",
        "ACTIVATED",
        "PLATELET",
        "MEGAKARYOCYTE",
        "DENDRITIC",
    }
    out = []
    for token in tokens:
        if len(token) < 3:
            continue
        upper = token.upper()
        if upper in ignore:
            continue
        if re.search(r"[A-Z]", token):
            out.append(upper)
    return out


def safe_float(value: str | None) -> float | None:
    if value in (None, "", "NaN", "nan"):
        return None
    return float(value)


def top_counter(counter: Counter[str]) -> tuple[str, int, float]:
    if not counter:
        return "", 0, 0.0
    label, count = counter.most_common(1)[0]
    total = sum(counter.values())
    return label, count, count / total if total else 0.0


def marker_lookup(rows: list[dict[str, str]], top_n: int = 5) -> dict[str, list[str]]:
    grouped: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        cluster = normalize_cluster(row["cluster"])
        gene = (row.get("gene_name") or "").strip()
        if gene and gene not in grouped[cluster]:
            grouped[cluster].append(gene)
    return {cluster: genes[:top_n] for cluster, genes in grouped.items()}


def batch_dominance(
    rows: list[dict[str, str]],
) -> tuple[dict[str, Counter[str]], dict[str, Counter[str]]]:
    gse_counts: dict[str, Counter[str]] = defaultdict(Counter)
    gsm_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        cluster = row["seurat_clusters"]
        gse_counts[cluster][row["source_gse"]] += 1
        gsm_counts[cluster][row["source_gsm"]] += 1
    return gse_counts, gsm_counts


def reference_overlap(
    rows: list[dict[str, str]],
) -> tuple[dict[str, Counter[str]], dict[str, Counter[str]]]:
    broad: dict[str, Counter[str]] = defaultdict(Counter)
    subtype: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        cluster = row["seurat_clusters"]
        broad[cluster][row["reference_celltype"]] += int(row["N"])
        subtype[cluster][row["reference_subtype"]] += int(row["N"])
    return broad, subtype


def priority_rank(priority: str) -> int:
    order = {"high": 0, "medium": 1, "low": 2}
    return order.get(priority, 3)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Review annotation consistency for integrated scATAC outputs."
    )
    parser.add_argument(
        "--analysis-dir",
        default=str(OUTPUT_ROOT / "integration_merged_without_unknown_analysis"),
        help="Analysis directory containing integration and annotation outputs",
    )
    args = parser.parse_args()

    analysis_dir = Path(args.analysis_dir)
    annotation_map_path = analysis_dir / "cluster_celltype_annotation_map.csv"
    marker_path = analysis_dir / "cluster_top_accessible_peaks.csv"
    metadata_path = analysis_dir / "integrated_metadata.csv.gz"
    summary_path = analysis_dir / "integration_summary.json"
    overlap_path = analysis_dir / "cluster_annotation_reference_overlap.csv"

    required = [annotation_map_path, marker_path, metadata_path, summary_path]
    if not all(path.exists() for path in required):
        raise SystemExit("Missing required analysis outputs for annotation review")

    annotation_rows = read_csv_rows(annotation_map_path)
    marker_rows = read_csv_rows(marker_path)
    metadata_rows = read_csv_rows(metadata_path)
    marker_genes = marker_lookup(marker_rows)
    gse_counts, gsm_counts = batch_dominance(metadata_rows)
    summary_payload = json.loads(summary_path.read_text(encoding="utf-8"))

    broad_overlap: dict[str, Counter[str]] = defaultdict(Counter)
    subtype_overlap: dict[str, Counter[str]] = defaultdict(Counter)
    if overlap_path.exists():
        broad_overlap, subtype_overlap = reference_overlap(read_csv_rows(overlap_path))

    validation_rows: list[dict[str, str | int | float]] = []
    review_blocks: list[str] = []
    flagged_clusters: list[dict[str, str | int | float]] = []

    for row in annotation_rows:
        cluster = row["seurat_clusters"]
        current_markers = marker_genes.get(cluster, [])
        expected_tokens = label_keywords(row["celltype"], row["celltype_subtype"])
        overlap_tokens = sorted(
            {
                token
                for token in expected_tokens
                if token in {gene.upper() for gene in current_markers}
            }
        )

        top_gse, top_gse_count, top_gse_fraction = top_counter(
            gse_counts.get(cluster, Counter())
        )
        top_gsm, top_gsm_count, top_gsm_fraction = top_counter(
            gsm_counts.get(cluster, Counter())
        )
        ref_broad, ref_broad_count, ref_broad_fraction = top_counter(
            broad_overlap.get(cluster, Counter())
        )
        ref_sub, ref_sub_count, ref_sub_fraction = top_counter(
            subtype_overlap.get(cluster, Counter())
        )

        reasons = []
        confidence = row.get("annotation_confidence", "")
        if confidence == "low":
            reasons.append("low_annotation_confidence")
        elif confidence == "medium":
            reasons.append("medium_annotation_confidence")
        if top_gsm_fraction >= 0.35:
            reasons.append("sample_dominant")
        if top_gse_fraction >= 0.8:
            reasons.append("gse_dominant")
        if ref_broad_fraction and ref_broad_fraction < 0.7:
            reasons.append("broad_label_mixed")
        if ref_sub_fraction and ref_sub_fraction < 0.5:
            reasons.append("subtype_mixed")
        if ref_broad and ref_broad != row["celltype"]:
            reasons.append("reference_broad_mismatch")
        if ref_sub and ref_sub != row["celltype_subtype"]:
            reasons.append("reference_subtype_mismatch")
        if expected_tokens and not overlap_tokens:
            reasons.append("label_marker_keyword_not_seen")

        high_priority_reasons = {
            "low_annotation_confidence",
            "sample_dominant",
            "reference_broad_mismatch",
            "reference_subtype_mismatch",
            "broad_label_mixed",
            "subtype_mixed",
        }

        review_priority = "low"
        if any(reason in high_priority_reasons for reason in reasons):
            review_priority = "high"
        elif confidence == "medium":
            review_priority = "medium"

        validation = {
            "seurat_clusters": cluster,
            "celltype": row["celltype"],
            "celltype_subtype": row["celltype_subtype"],
            "annotation_confidence": confidence,
            "total_cells": int(row.get("total_cells") or 0),
            "celltype_fraction": safe_float(row.get("celltype_fraction")),
            "subtype_fraction": safe_float(row.get("subtype_fraction")),
            "top_marker_genes": "; ".join(current_markers),
            "expected_label_keywords": "; ".join(expected_tokens),
            "matched_label_keywords": "; ".join(overlap_tokens),
            "top_gse": top_gse,
            "top_gse_count": top_gse_count,
            "top_gse_fraction": round(top_gse_fraction, 4),
            "top_gsm": top_gsm,
            "top_gsm_count": top_gsm_count,
            "top_gsm_fraction": round(top_gsm_fraction, 4),
            "top_reference_celltype": ref_broad,
            "top_reference_celltype_count": ref_broad_count,
            "top_reference_celltype_fraction": round(ref_broad_fraction, 4),
            "top_reference_subtype": ref_sub,
            "top_reference_subtype_count": ref_sub_count,
            "top_reference_subtype_fraction": round(ref_sub_fraction, 4),
            "review_priority": review_priority,
            "review_reasons": "; ".join(reasons),
        }
        validation_rows.append(validation)
        if review_priority != "low":
            flagged_clusters.append(validation)

    validation_rows.sort(
        key=lambda row: (
            priority_rank(str(row["review_priority"])),
            int(str(row["seurat_clusters"])),
        )
    )
    flagged_clusters.sort(
        key=lambda row: (
            priority_rank(str(row["review_priority"])),
            int(str(row["seurat_clusters"])),
        )
    )

    summary_csv_path = analysis_dir / "cluster_annotation_validation_summary.csv"
    with summary_csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(validation_rows[0].keys()))
        writer.writeheader()
        writer.writerows(validation_rows)

    review_lines = [
        "# Annotation validation report",
        "",
        f"- Analysis directory: `{analysis_dir}`",
        f"- Cells integrated: {summary_payload.get('total_cells', 'NA')}",
        f"- Clusters reviewed: {summary_payload.get('clusters', 'NA')}",
        f"- Flagged clusters: {len(flagged_clusters)}",
        "- Validation combines current top marker genes, dominant GSE/GSM fractions, and prior-annotation overlap when available.",
        "",
        "## Review priorities",
        "",
    ]
    for row in flagged_clusters:
        review_lines.extend(
            [
                f"### Cluster {row['seurat_clusters']} — {row['celltype']} / {row['celltype_subtype']}",
                f"- Priority: `{row['review_priority']}`",
                f"- Confidence: `{row['annotation_confidence']}`",
                f"- Current top marker genes: `{row['top_marker_genes']}`",
                f"- Dominant GSE: `{row['top_gse']}` ({row['top_gse_fraction']:.2%})",
                f"- Dominant GSM: `{row['top_gsm']}` ({row['top_gsm_fraction']:.2%})",
                f"- Dominant prior broad label: `{row['top_reference_celltype']}` ({row['top_reference_celltype_fraction']:.2%})",
                f"- Dominant prior subtype: `{row['top_reference_subtype']}` ({row['top_reference_subtype_fraction']:.2%})",
                f"- Review reasons: `{row['review_reasons']}`",
                "",
            ]
        )
    if not flagged_clusters:
        review_lines.extend(["- No clusters triggered the current review rules.", ""])

    review_md_path = analysis_dir / "annotation_validation_report.md"
    review_md_path.write_text("\n".join(review_lines), encoding="utf-8")

    print(f"[review] wrote {summary_csv_path}")
    print(f"[review] wrote {review_md_path}")


if __name__ == "__main__":
    main()
