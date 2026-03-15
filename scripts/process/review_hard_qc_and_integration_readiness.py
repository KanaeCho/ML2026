#!/usr/bin/env python3

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"
REVIEW_DIR = OUTPUT_ROOT / "hard_qc_review"

THRESHOLDS = {
    "nCount_ATAC_min": 1000.0,
    "nCount_ATAC_max": 100000.0,
    "TSS.enrichment_min": 4.0,
    "FRiP_min": 0.35,
    "blacklist_fraction_max": 0.05,
    "nucleosome_signal_max": 4.0,
}

METRICS = [
    "nCount_ATAC",
    "TSS.enrichment",
    "FRiP",
    "blacklist_fraction",
    "nucleosome_signal",
]

METRIC_PLOT_CONFIG = {
    "nCount_ATAC": {
        "label": "nCount_ATAC",
        "threshold": 1000.0,
        "mode": "min",
        "log10": True,
    },
    "TSS.enrichment": {
        "label": "TSS enrichment",
        "threshold": 4.0,
        "mode": "min",
        "log10": False,
    },
    "FRiP": {"label": "FRiP", "threshold": 0.35, "mode": "min", "log10": False},
    "blacklist_fraction": {
        "label": "Blacklist fraction",
        "threshold": 0.05,
        "mode": "max",
        "log10": False,
    },
    "nucleosome_signal": {
        "label": "Nucleosome signal",
        "threshold": 4.0,
        "mode": "max",
        "log10": False,
    },
}

DATASET_COLORS = {
    "GSE190992": "#0f766e",
    "GSE283744": "#b45309",
}


@dataclass(frozen=True)
class SamplePaths:
    gse: str
    gsm: str
    metadata_qc: Path
    metadata_integration_qc: Path


def to_float(value: str) -> float | None:
    value = value.strip()
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def safe_percent(part: int, whole: int) -> float:
    return (part / whole * 100.0) if whole else 0.0


def percentile(sorted_values: list[float], q: float) -> float | None:
    if not sorted_values:
        return None
    if len(sorted_values) == 1:
        return sorted_values[0]
    pos = (len(sorted_values) - 1) * q
    low = math.floor(pos)
    high = math.ceil(pos)
    if low == high:
        return sorted_values[low]
    frac = pos - low
    return sorted_values[low] * (1 - frac) + sorted_values[high] * frac


def fmt_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def discover_samples(output_root: Path) -> list[SamplePaths]:
    samples: list[SamplePaths] = []
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
            metadata_qc = sample_dir / "metadata_qc.csv"
            metadata_integration_qc = (
                sample_dir / "integration_qc" / "metadata_integration_qc.csv"
            )
            if metadata_qc.exists() and metadata_integration_qc.exists():
                samples.append(
                    SamplePaths(
                        gse=gse_dir.name,
                        gsm=sample_dir.name,
                        metadata_qc=metadata_qc,
                        metadata_integration_qc=metadata_integration_qc,
                    )
                )
    return samples


def evaluate_failures(row: dict[str, str]) -> list[str]:
    failures: list[str] = []
    ncount = to_float(row.get("nCount_ATAC", ""))
    tss = to_float(row.get("TSS.enrichment", ""))
    frip = to_float(row.get("FRiP", ""))
    blacklist = to_float(row.get("blacklist_fraction", ""))
    nucleosome = to_float(row.get("nucleosome_signal", ""))

    if ncount is None or ncount < THRESHOLDS["nCount_ATAC_min"]:
        failures.append("nCount_ATAC_low")
    if ncount is None or ncount > THRESHOLDS["nCount_ATAC_max"]:
        failures.append("nCount_ATAC_high")
    if tss is None or tss < THRESHOLDS["TSS.enrichment_min"]:
        failures.append("TSS.enrichment_low")
    if frip is None or frip < THRESHOLDS["FRiP_min"]:
        failures.append("FRiP_low")
    if blacklist is None or blacklist > THRESHOLDS["blacklist_fraction_max"]:
        failures.append("blacklist_fraction_high")
    if nucleosome is None or nucleosome > THRESHOLDS["nucleosome_signal_max"]:
        failures.append("nucleosome_signal_high")
    return failures


def read_post_barcodes(path: Path) -> set[str]:
    barcodes: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            barcode = row.get("cell_barcode", "")
            if barcode:
                barcodes.add(barcode)
    return barcodes


def collect_review(samples: list[SamplePaths]) -> dict:
    overall_before: dict[str, list[float]] = {metric: [] for metric in METRICS}
    overall_after: dict[str, list[float]] = {metric: [] for metric in METRICS}
    removed_metrics: dict[str, list[float]] = {metric: [] for metric in METRICS}
    dataset_before: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {metric: [] for metric in METRICS}
    )
    dataset_after: dict[str, dict[str, list[float]]] = defaultdict(
        lambda: {metric: [] for metric in METRICS}
    )
    failure_counts: Counter[str] = Counter()
    failure_counts_by_dataset: dict[str, Counter[str]] = defaultdict(Counter)
    combo_counts: Counter[str] = Counter()
    sample_rows: list[dict[str, object]] = []
    validation_rows: list[dict[str, object]] = []

    total_before = 0
    total_after = 0

    for sample in samples:
        post_barcodes = read_post_barcodes(sample.metadata_integration_qc)
        before_rows = 0
        recomputed_after = 0
        metric_before = {metric: [] for metric in METRICS}
        metric_after = {metric: [] for metric in METRICS}
        sample_failure_counts: Counter[str] = Counter()

        with sample.metadata_qc.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                before_rows += 1
                failures = evaluate_failures(row)
                keep = not failures
                barcode = row.get("cell_barcode", "")
                if keep:
                    recomputed_after += 1
                if keep != (barcode in post_barcodes):
                    raise ValueError(
                        f"Hard-QC validation mismatch for {sample.gse}/{sample.gsm} barcode {barcode}: "
                        f"recomputed keep={keep}, post file contains={barcode in post_barcodes}"
                    )

                for metric in METRICS:
                    value = to_float(row.get(metric, ""))
                    if value is None:
                        continue
                    overall_before[metric].append(value)
                    dataset_before[sample.gse][metric].append(value)
                    metric_before[metric].append(value)
                    if keep:
                        overall_after[metric].append(value)
                        dataset_after[sample.gse][metric].append(value)
                        metric_after[metric].append(value)
                    else:
                        removed_metrics[metric].append(value)

                if failures:
                    unique_failures = sorted(set(failures))
                    combo_counts[" + ".join(unique_failures)] += 1
                    for failure in unique_failures:
                        failure_counts[failure] += 1
                        failure_counts_by_dataset[sample.gse][failure] += 1
                        sample_failure_counts[failure] += 1

        total_before += before_rows
        total_after += recomputed_after
        validation_rows.append(
            {
                "gse": sample.gse,
                "gsm": sample.gsm,
                "metadata_qc_cells": before_rows,
                "recomputed_integration_qc_cells": recomputed_after,
                "file_integration_qc_cells": len(post_barcodes),
            }
        )
        sample_rows.append(
            {
                "gse": sample.gse,
                "gsm": sample.gsm,
                "before_cells": before_rows,
                "after_cells": recomputed_after,
                "retention_rate": safe_percent(recomputed_after, before_rows),
                "top_failure_reason": sample_failure_counts.most_common(1)[0][0]
                if sample_failure_counts
                else "",
                "failed_nCount_ATAC_low": sample_failure_counts["nCount_ATAC_low"],
                "failed_nCount_ATAC_high": sample_failure_counts["nCount_ATAC_high"],
                "failed_TSS_enrichment_low": sample_failure_counts[
                    "TSS.enrichment_low"
                ],
                "failed_FRiP_low": sample_failure_counts["FRiP_low"],
                "failed_blacklist_fraction_high": sample_failure_counts[
                    "blacklist_fraction_high"
                ],
                "failed_nucleosome_signal_high": sample_failure_counts[
                    "nucleosome_signal_high"
                ],
            }
        )

    dataset_rows = []
    for dataset in sorted(dataset_before):
        before_count = len(dataset_before[dataset]["nCount_ATAC"])
        after_count = len(dataset_after[dataset]["nCount_ATAC"])
        dataset_rows.append(
            {
                "gse": dataset,
                "before_cells": before_count,
                "after_cells": after_count,
                "retention_rate": safe_percent(after_count, before_count),
                "failed_nCount_ATAC_low": failure_counts_by_dataset[dataset][
                    "nCount_ATAC_low"
                ],
                "failed_nCount_ATAC_high": failure_counts_by_dataset[dataset][
                    "nCount_ATAC_high"
                ],
                "failed_TSS_enrichment_low": failure_counts_by_dataset[dataset][
                    "TSS.enrichment_low"
                ],
                "failed_FRiP_low": failure_counts_by_dataset[dataset]["FRiP_low"],
                "failed_blacklist_fraction_high": failure_counts_by_dataset[dataset][
                    "blacklist_fraction_high"
                ],
                "failed_nucleosome_signal_high": failure_counts_by_dataset[dataset][
                    "nucleosome_signal_high"
                ],
            }
        )

    return {
        "sample_rows": sample_rows,
        "dataset_rows": dataset_rows,
        "validation_rows": validation_rows,
        "overall_before": overall_before,
        "overall_after": overall_after,
        "dataset_before": dataset_before,
        "dataset_after": dataset_after,
        "removed_metrics": removed_metrics,
        "failure_counts": failure_counts,
        "failure_counts_by_dataset": failure_counts_by_dataset,
        "combo_counts": combo_counts,
        "total_before": total_before,
        "total_after": total_after,
    }


def summarize_metric(values: list[float]) -> dict[str, float | None]:
    ordered = sorted(values)
    if not ordered:
        return {
            "n": 0,
            "min": None,
            "p05": None,
            "p25": None,
            "median": None,
            "p75": None,
            "p95": None,
            "max": None,
        }
    return {
        "n": len(ordered),
        "min": ordered[0],
        "p05": percentile(ordered, 0.05),
        "p25": percentile(ordered, 0.25),
        "median": median(ordered),
        "p75": percentile(ordered, 0.75),
        "p95": percentile(ordered, 0.95),
        "max": ordered[-1],
    }


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_metric_summary(review: dict) -> None:
    rows: list[dict[str, object]] = []
    for metric in METRICS:
        for stage, values in (
            ("before_hard_qc", review["overall_before"][metric]),
            ("after_hard_qc", review["overall_after"][metric]),
        ):
            summary = summarize_metric(values)
            rows.append(
                {
                    "scope": "overall",
                    "dataset": "all",
                    "metric": metric,
                    "stage": stage,
                    **summary,
                }
            )
    for dataset in sorted(review["dataset_before"]):
        dataset_before = review["dataset_before"][dataset]
        dataset_after = review["dataset_after"][dataset]
        for metric in METRICS:
            rows.append(
                {
                    "scope": "dataset",
                    "dataset": dataset,
                    "metric": metric,
                    "stage": "before_hard_qc",
                    **summarize_metric(dataset_before[metric]),
                }
            )
            rows.append(
                {
                    "scope": "dataset",
                    "dataset": dataset,
                    "metric": metric,
                    "stage": "after_hard_qc",
                    **summarize_metric(dataset_after[metric]),
                }
            )
    write_csv(
        REVIEW_DIR / "hard_qc_metric_summary.csv",
        rows,
        [
            "scope",
            "dataset",
            "metric",
            "stage",
            "n",
            "min",
            "p05",
            "p25",
            "median",
            "p75",
            "p95",
            "max",
        ],
    )


def write_failure_summaries(review: dict) -> None:
    reason_rows = []
    for reason, count in sorted(review["failure_counts"].items()):
        reason_rows.append(
            {"scope": "overall", "dataset": "all", "reason": reason, "cells": count}
        )
    for dataset, counts in sorted(review["failure_counts_by_dataset"].items()):
        for reason, count in sorted(counts.items()):
            reason_rows.append(
                {
                    "scope": "dataset",
                    "dataset": dataset,
                    "reason": reason,
                    "cells": count,
                }
            )
    write_csv(
        REVIEW_DIR / "hard_qc_failure_reasons.csv",
        reason_rows,
        ["scope", "dataset", "reason", "cells"],
    )

    combo_rows = [
        {"failure_combo": combo, "cells": count}
        for combo, count in review["combo_counts"].most_common(20)
    ]
    write_csv(
        REVIEW_DIR / "hard_qc_failure_combinations_top20.csv",
        combo_rows,
        ["failure_combo", "cells"],
    )


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#fffdf8"/>',
        "<style>text{font-family:Arial,Helvetica,sans-serif;fill:#1f2937} .small{font-size:12px} .title{font-size:22px;font-weight:700} .panel{font-size:15px;font-weight:700} .axis{font-size:11px;fill:#4b5563} .legend{font-size:12px}</style>",
    ]


def histogram(values: list[float], bins: int, vmin: float, vmax: float) -> list[int]:
    counts = [0] * bins
    if vmax <= vmin:
        counts[-1] = len(values)
        return counts
    for value in values:
        if value < vmin:
            idx = 0
        elif value >= vmax:
            idx = bins - 1
        else:
            idx = int((value - vmin) / (vmax - vmin) * bins)
            if idx == bins:
                idx = bins - 1
        counts[idx] += 1
    return counts


def format_axis_value(value: float, log10_mode: bool) -> str:
    if log10_mode:
        raw = 10**value
        if raw >= 1000:
            return f"{int(round(raw / 1000))}k"
        return str(int(round(raw)))
    return f"{value:.2f}" if abs(value) < 10 else f"{value:.1f}"


def draw_hist_panel(
    lines: list[str],
    x: int,
    y: int,
    width: int,
    height: int,
    title: str,
    before_values: list[float],
    after_values: list[float],
    threshold: float,
    mode: str,
    log10_mode: bool,
) -> None:
    panel_pad = 36
    plot_x = x + panel_pad
    plot_y = y + 20
    plot_w = width - panel_pad - 18
    plot_h = height - 60
    lines.append(
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="10" fill="#fffaf0" stroke="#e5e7eb"/>'
    )
    lines.append(f'<text x="{x + 12}" y="{y + 18}" class="panel">{title}</text>')

    transform = (
        (lambda values: [math.log10(v + 1.0) for v in values])
        if log10_mode
        else (lambda values: list(values))
    )
    before_t = transform(before_values)
    after_t = transform(after_values)
    all_values = before_t + after_t
    if not all_values:
        return
    vmin = min(all_values)
    vmax = max(all_values)
    if vmax == vmin:
        vmax = vmin + 1.0

    bins = 28
    before_counts = histogram(before_t, bins, vmin, vmax)
    after_counts = histogram(after_t, bins, vmin, vmax)
    max_count = max(max(before_counts), max(after_counts), 1)
    bin_w = plot_w / bins

    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y + plot_h}" x2="{plot_x + plot_w}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )
    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )

    for idx, count in enumerate(before_counts):
        bar_h = plot_h * count / max_count
        lines.append(
            f'<rect x="{plot_x + idx * bin_w + 1:.2f}" y="{plot_y + plot_h - bar_h:.2f}" width="{max(bin_w - 2, 1):.2f}" height="{bar_h:.2f}" fill="#f59e0b" fill-opacity="0.35"/>'
        )
    for idx, count in enumerate(after_counts):
        bar_h = plot_h * count / max_count
        lines.append(
            f'<rect x="{plot_x + idx * bin_w + 1:.2f}" y="{plot_y + plot_h - bar_h:.2f}" width="{max(bin_w - 2, 1):.2f}" height="{bar_h:.2f}" fill="#0f766e" fill-opacity="0.45"/>'
        )

    threshold_t = math.log10(threshold + 1.0) if log10_mode else threshold
    threshold_x = plot_x + (threshold_t - vmin) / (vmax - vmin) * plot_w
    threshold_x = max(plot_x, min(plot_x + plot_w, threshold_x))
    lines.append(
        f'<line x1="{threshold_x:.2f}" y1="{plot_y}" x2="{threshold_x:.2f}" y2="{plot_y + plot_h}" stroke="#dc2626" stroke-dasharray="5,4"/>'
    )
    anchor = threshold_x + 4 if mode == "min" else threshold_x - 4
    text_anchor = "start" if mode == "min" else "end"
    label = f"{'>=' if mode == 'min' else '<='} {threshold:g}"
    lines.append(
        f'<text x="{anchor:.2f}" y="{plot_y + 12}" class="axis" text-anchor="{text_anchor}">{label}</text>'
    )

    for frac in (0.0, 0.5, 1.0):
        xv = plot_x + plot_w * frac
        value = vmin + (vmax - vmin) * frac
        lines.append(
            f'<text x="{xv:.2f}" y="{plot_y + plot_h + 16}" class="axis" text-anchor="middle">{format_axis_value(value, log10_mode)}</text>'
        )
    for frac in (0.0, 0.5, 1.0):
        yv = plot_y + plot_h - plot_h * frac
        lines.append(
            f'<text x="{plot_x - 6}" y="{yv + 4:.2f}" class="axis" text-anchor="end">{int(max_count * frac)}</text>'
        )


def draw_retention_panel(
    lines: list[str],
    x: int,
    y: int,
    width: int,
    height: int,
    sample_rows: list[dict[str, object]],
) -> None:
    lines.append(
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="10" fill="#fffaf0" stroke="#e5e7eb"/>'
    )
    lines.append(
        f'<text x="{x + 12}" y="{y + 18}" class="panel">Per-sample hard-QC retention</text>'
    )
    plot_x = x + 42
    plot_y = y + 26
    plot_w = width - 56
    plot_h = height - 54
    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y + plot_h}" x2="{plot_x + plot_w}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )
    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )
    ordered = sorted(sample_rows, key=lambda row: float(row["retention_rate"]))
    bar_w = max(plot_w / max(len(ordered), 1) - 2, 1)
    for idx, row in enumerate(ordered):
        rate = float(row["retention_rate"])
        bar_h = plot_h * rate / 100.0
        color = DATASET_COLORS.get(str(row["gse"]), "#64748b")
        bx = plot_x + idx * (bar_w + 2)
        by = plot_y + plot_h - bar_h
        lines.append(
            f'<rect x="{bx:.2f}" y="{by:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" fill="{color}" fill-opacity="0.8"/>'
        )
        if idx % 3 == 0:
            lines.append(
                f'<text x="{bx + bar_w / 2:.2f}" y="{plot_y + plot_h + 14}" class="axis" text-anchor="end" transform="rotate(-60 {bx + bar_w / 2:.2f},{plot_y + plot_h + 14})">{row["gsm"]}</text>'
            )
    for tick in range(0, 101, 25):
        ty = plot_y + plot_h - plot_h * tick / 100.0
        lines.append(
            f'<line x1="{plot_x}" y1="{ty:.2f}" x2="{plot_x + plot_w}" y2="{ty:.2f}" stroke="#e5e7eb"/>'
        )
        lines.append(
            f'<text x="{plot_x - 8}" y="{ty + 4:.2f}" class="axis" text-anchor="end">{tick}%</text>'
        )
    legend_x = x + width - 170
    legend_y = y + 22
    for idx, dataset in enumerate(sorted(DATASET_COLORS)):
        color = DATASET_COLORS[dataset]
        lines.append(
            f'<rect x="{legend_x}" y="{legend_y + idx * 18}" width="12" height="12" fill="{color}"/>'
        )
        lines.append(
            f'<text x="{legend_x + 18}" y="{legend_y + 10 + idx * 18}" class="legend">{dataset}</text>'
        )


def draw_failure_panel(
    lines: list[str],
    x: int,
    y: int,
    width: int,
    height: int,
    dataset_rows: list[dict[str, object]],
) -> None:
    lines.append(
        f'<rect x="{x}" y="{y}" width="{width}" height="{height}" rx="10" fill="#fffaf0" stroke="#e5e7eb"/>'
    )
    lines.append(
        f'<text x="{x + 12}" y="{y + 18}" class="panel">Cells failing each hard-QC rule</text>'
    )
    reasons = [
        ("failed_nCount_ATAC_low", "low count"),
        ("failed_nCount_ATAC_high", "high count"),
        ("failed_TSS_enrichment_low", "low TSS"),
        ("failed_FRiP_low", "low FRiP"),
        ("failed_blacklist_fraction_high", "high blacklist"),
        ("failed_nucleosome_signal_high", "high nucleosome"),
    ]
    plot_x = x + 54
    plot_y = y + 34
    plot_w = width - 72
    plot_h = height - 60
    max_count = 1
    for row in dataset_rows:
        for key, _ in reasons:
            max_count = max(max_count, int(row[key]))
    group_w = plot_w / len(reasons)
    bar_w = max(group_w / (len(dataset_rows) + 1), 10)
    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y + plot_h}" x2="{plot_x + plot_w}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )
    lines.append(
        f'<line x1="{plot_x}" y1="{plot_y}" x2="{plot_x}" y2="{plot_y + plot_h}" stroke="#374151"/>'
    )
    for tick in range(0, 5):
        frac = tick / 4
        ty = plot_y + plot_h - plot_h * frac
        lines.append(
            f'<line x1="{plot_x}" y1="{ty:.2f}" x2="{plot_x + plot_w}" y2="{ty:.2f}" stroke="#e5e7eb"/>'
        )
        lines.append(
            f'<text x="{plot_x - 8}" y="{ty + 4:.2f}" class="axis" text-anchor="end">{int(max_count * frac)}</text>'
        )
    for ridx, (key, label) in enumerate(reasons):
        gx = plot_x + ridx * group_w
        lines.append(
            f'<text x="{gx + group_w / 2:.2f}" y="{plot_y + plot_h + 16}" class="axis" text-anchor="middle">{label}</text>'
        )
        for didx, row in enumerate(dataset_rows):
            count = int(row[key])
            bar_h = plot_h * count / max_count
            color = DATASET_COLORS.get(str(row["gse"]), "#64748b")
            bx = gx + 8 + didx * (bar_w + 4)
            by = plot_y + plot_h - bar_h
            lines.append(
                f'<rect x="{bx:.2f}" y="{by:.2f}" width="{bar_w:.2f}" height="{bar_h:.2f}" fill="{color}" fill-opacity="0.82"/>'
            )
    legend_x = x + width - 170
    legend_y = y + 22
    for idx, dataset in enumerate(sorted(DATASET_COLORS)):
        color = DATASET_COLORS[dataset]
        lines.append(
            f'<rect x="{legend_x}" y="{legend_y + idx * 18}" width="12" height="12" fill="{color}"/>'
        )
        lines.append(
            f'<text x="{legend_x + 18}" y="{legend_y + 10 + idx * 18}" class="legend">{dataset}</text>'
        )


def write_metric_svg(review: dict) -> None:
    width = 1500
    height = 1140
    lines = svg_header(width, height)
    lines.append('<text x="36" y="42" class="title">Hard-QC before/after review</text>')
    lines.append(
        f'<text x="36" y="66" class="small">Before: {review["total_before"]:,} cells | After: {review["total_after"]:,} cells | Retained: {safe_percent(review["total_after"], review["total_before"]):.2f}% | Thresholds from scripts/process/apply_integration_hard_qc.py</text>'
    )
    panel_w = 460
    panel_h = 260
    positions = [
        (36, 96, "nCount_ATAC"),
        (518, 96, "TSS.enrichment"),
        (1000, 96, "FRiP"),
        (36, 380, "blacklist_fraction"),
        (518, 380, "nucleosome_signal"),
    ]
    for x, y, metric in positions:
        cfg = METRIC_PLOT_CONFIG[metric]
        draw_hist_panel(
            lines,
            x,
            y,
            panel_w,
            panel_h,
            cfg["label"],
            review["overall_before"][metric],
            review["overall_after"][metric],
            cfg["threshold"],
            cfg["mode"],
            cfg["log10"],
        )
    draw_retention_panel(lines, 1000, 380, 460, 300, review["sample_rows"])
    draw_failure_panel(lines, 36, 704, 1424, 360, review["dataset_rows"])
    lines.append(
        '<rect x="1190" y="70" width="12" height="12" fill="#f59e0b" fill-opacity="0.35"/>'
    )
    lines.append('<text x="1208" y="80" class="legend">before hard-QC</text>')
    lines.append(
        '<rect x="1320" y="70" width="12" height="12" fill="#0f766e" fill-opacity="0.45"/>'
    )
    lines.append('<text x="1338" y="80" class="legend">after hard-QC</text>')
    lines.append(
        '<line x1="1450" y1="76" x2="1478" y2="76" stroke="#dc2626" stroke-dasharray="5,4"/>'
    )
    lines.append('<text x="1482" y="80" class="legend">threshold</text>')
    lines.append("</svg>")
    (REVIEW_DIR / "hard_qc_before_after_metrics.svg").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def estimate_sparse_memory_gb(
    matrix_nnz: int, matrix_cols: int, matrix_rows: int
) -> dict[str, float]:
    index_bytes = 4.0
    data_bytes = 8.0
    indptr_bytes = 4.0
    csc_bytes = (
        matrix_nnz * (index_bytes + data_bytes) + (matrix_cols + 1) * indptr_bytes
    )
    double_copy_bytes = csc_bytes * 2.0
    lsi_embedding_bytes = matrix_cols * 30 * 8.0
    feature_stats_bytes = matrix_rows * 8.0 * 6.0
    return {
        "csc_matrix_gb": csc_bytes / 1024**3,
        "double_copy_gb": double_copy_bytes / 1024**3,
        "lsi_embedding_30d_gb": lsi_embedding_bytes / 1024**3,
        "feature_stats_overhead_gb": feature_stats_bytes / 1024**3,
    }


def write_report(
    review: dict, merge_summary: dict, mem_info: dict[str, object]
) -> None:
    sample_rows = sorted(
        review["sample_rows"], key=lambda row: float(row["retention_rate"])
    )
    worst_samples = sample_rows[:5]
    best_samples = sample_rows[-5:]
    top_combos = review["combo_counts"].most_common(5)
    memory_est = estimate_sparse_memory_gb(
        int(merge_summary.get("matrix_nnz", 0)),
        int(merge_summary.get("matrix_cols", 0)),
        int(merge_summary.get("matrix_rows", 0)),
    )
    available_gb = float(mem_info.get("available_gb", 0.0))
    lines = []
    lines.append("# Hard-QC review and integration readiness")
    lines.append("")
    lines.append("## Current state")
    lines.append("")
    lines.append(
        f"- Single-sample QC complete for {len(review['sample_rows'])} samples."
    )
    lines.append(
        f"- Hard-QC complete for all samples: {review['total_after']:,} / {review['total_before']:,} cells retained ({safe_percent(review['total_after'], review['total_before']):.2f}%)."
    )
    lines.append(
        f"- Merged integration matrix exists with {merge_summary.get('features', 0):,} peaks, {merge_summary.get('cells', 0):,} cells, and {merge_summary.get('matrix_nnz', 0):,} non-zero entries."
    )
    lines.append("")
    lines.append("## How the hard thresholds are currently determined")
    lines.append("")
    lines.append(
        "- The thresholds are fixed constants in `scripts/process/apply_integration_hard_qc.py`, not learned from the 28-sample cohort."
    )
    lines.append(
        "- Active rules: `nCount_ATAC >= 1000`, `nCount_ATAC <= 100000`, `TSS.enrichment >= 4`, `FRiP >= 0.35`, `blacklist_fraction <= 0.05`, `nucleosome_signal <= 4`."
    )
    lines.append(
        "- `unique_ratio` is intentionally not enforced by default because its availability is inconsistent across datasets."
    )
    lines.append("")
    lines.append("## What the filter is doing")
    lines.append("")
    for row in sorted(review["dataset_rows"], key=lambda item: str(item["gse"])):
        lines.append(
            f"- {row['gse']}: {int(row['after_cells']):,} / {int(row['before_cells']):,} cells retained ({float(row['retention_rate']):.2f}%). "
            f"Most common losses are TSS/FRiP related rather than blacklist or nucleosome cutoffs."
        )
    lines.append("")
    lines.append("Top failure combinations:")
    for combo, count in top_combos:
        lines.append(f"- {combo}: {count:,} cells")
    lines.append("")
    lines.append("Worst-retained samples:")
    for row in worst_samples:
        lines.append(
            f"- {row['gse']}/{row['gsm']}: {float(row['retention_rate']):.2f}% retained, top failure `{row['top_failure_reason']}`."
        )
    lines.append("")
    lines.append("Best-retained samples:")
    for row in reversed(best_samples):
        lines.append(
            f"- {row['gse']}/{row['gsm']}: {float(row['retention_rate']):.2f}% retained."
        )
    lines.append("")
    lines.append("## Where the project is blocked")
    lines.append("")
    lines.append(
        "- The repo already has hard-QC outputs and a merged matrix, but it does not yet have a committed downstream integration script that runs LSI/Harmony/UMAP on the merged object."
    )
    lines.append(
        "- The hard-threshold step was implemented before its cohort-level evaluation, so the blocker is now justification and scaling strategy, not the filtering code itself."
    )
    lines.append(
        "- `AGENTS.md` was stale and still said only 3/28 samples had completed QC, so the docs needed to catch up with the real outputs."
    )
    lines.append("")
    lines.append("## Merge-first vs batch-correct-while-merging")
    lines.append("")
    lines.append(
        "- For scATAC, merge-first on a common peak set is the standard starting point; batch correction is then applied on low-dimensional embeddings (typically LSI), not on the raw peak matrix itself."
    )
    lines.append(
        "- The current merged matrix is therefore the correct canonical input, but the risky step is trying to hold the full sparse count matrix plus all Seurat/Signac intermediates in memory at once."
    )
    lines.append(
        "- A safer workflow here is two-stage: use the merged matrix as truth, but start with a sketch/subsample for method tuning and QC diagnostics, then run the full merged object only once the LSI/Harmony settings are fixed."
    )
    lines.append("")
    lines.append("## Memory assessment")
    lines.append("")
    lines.append(
        f"- Current machine snapshot: {mem_info.get('cpu_count', 'NA')} CPU threads, {mem_info.get('total_gb', 'NA')} GiB RAM total, {mem_info.get('available_gb', 'NA')} GiB currently available, {mem_info.get('disk_free_gb', 'NA')} GiB disk free."
    )
    lines.append(
        f"- Sparse CSC storage estimate for the merged count matrix alone: ~{memory_est['csc_matrix_gb']:.2f} GiB."
    )
    lines.append(
        f"- If the workflow creates one extra matrix-sized copy during TF-IDF/SVD prep, budget rises to ~{memory_est['double_copy_gb']:.2f} GiB before Seurat object overhead."
    )
    lines.append(
        f"- The 30D cell embedding itself is small (~{memory_est['lsi_embedding_30d_gb']:.2f} GiB), so Harmony is not the true bottleneck; matrix loading and LSI prep are."
    )
    if available_gb and memory_est["double_copy_gb"] > available_gb * 0.6:
        lines.append(
            "- Conclusion: full in-memory integration may fit on this machine only if the implementation is careful, but it is not comfortable enough to be the first experiment. A sketch-first pass is the safer choice."
        )
    else:
        lines.append(
            "- Conclusion: the machine likely can hold the sparse matrix, but the first full run still deserves caution because Seurat/Signac often creates additional transient objects during preprocessing."
        )
    lines.append("")
    lines.append("## Recommended next step")
    lines.append("")
    lines.append(
        "1. Use the merged hard-QC matrix as the source of truth, but build a balanced sketch for exploratory LSI/Harmony/UMAP parameter tuning."
    )
    lines.append(
        "2. Confirm batch variables (`GSE`, `GSM`) and plot pre-Harmony LSI colored by both dataset and QC metrics."
    )
    lines.append(
        "3. Run Harmony on the sketch first; if mixing and marker structure look reasonable, scale the same settings to the full merged object."
    )
    lines.append(
        "4. Keep the hard-QC review outputs in this directory as the written justification for the current filter before moving to integration."
    )
    (REVIEW_DIR / "integration_readiness_report.md").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def current_memory_info() -> dict[str, object]:
    import shutil

    total_kib = 0
    available_kib = 0
    with Path("/proc/meminfo").open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("MemTotal:"):
                total_kib = int(line.split()[1])
            elif line.startswith("MemAvailable:"):
                available_kib = int(line.split()[1])
    usage = shutil.disk_usage(ROOT)
    return {
        "cpu_count": os_cpu_count(),
        "total_gb": round(total_kib / 1024**2, 2),
        "available_gb": round(available_kib / 1024**2, 2),
        "disk_free_gb": round(usage.free / 1024**3, 2),
    }


def os_cpu_count() -> int | None:
    try:
        import os

        return os.cpu_count()
    except Exception:
        return None


def main() -> None:
    REVIEW_DIR.mkdir(parents=True, exist_ok=True)
    samples = discover_samples(OUTPUT_ROOT)
    if not samples:
        raise SystemExit(
            "No samples with both metadata_qc.csv and integration_qc outputs found."
        )

    review = collect_review(samples)

    sample_rows = sorted(
        review["sample_rows"], key=lambda row: (str(row["gse"]), str(row["gsm"]))
    )
    write_csv(
        REVIEW_DIR / "hard_qc_sample_summary.csv",
        sample_rows,
        [
            "gse",
            "gsm",
            "before_cells",
            "after_cells",
            "retention_rate",
            "top_failure_reason",
            "failed_nCount_ATAC_low",
            "failed_nCount_ATAC_high",
            "failed_TSS_enrichment_low",
            "failed_FRiP_low",
            "failed_blacklist_fraction_high",
            "failed_nucleosome_signal_high",
        ],
    )
    write_csv(
        REVIEW_DIR / "hard_qc_validation.csv",
        review["validation_rows"],
        [
            "gse",
            "gsm",
            "metadata_qc_cells",
            "recomputed_integration_qc_cells",
            "file_integration_qc_cells",
        ],
    )
    write_csv(
        REVIEW_DIR / "hard_qc_dataset_summary.csv",
        sorted(review["dataset_rows"], key=lambda row: str(row["gse"])),
        [
            "gse",
            "before_cells",
            "after_cells",
            "retention_rate",
            "failed_nCount_ATAC_low",
            "failed_nCount_ATAC_high",
            "failed_TSS_enrichment_low",
            "failed_FRiP_low",
            "failed_blacklist_fraction_high",
            "failed_nucleosome_signal_high",
        ],
    )
    write_failure_summaries(review)
    write_metric_summary(review)
    write_metric_svg(review)

    merge_summary_path = OUTPUT_ROOT / "integration_merged" / "merge_summary.json"
    merge_summary = (
        json.loads(merge_summary_path.read_text(encoding="utf-8"))
        if merge_summary_path.exists()
        else {}
    )
    mem_info = current_memory_info()
    write_report(review, merge_summary, mem_info)

    machine_summary = {
        "thresholds": THRESHOLDS,
        "samples": len(samples),
        "total_before": review["total_before"],
        "total_after": review["total_after"],
        "retention_rate": safe_percent(review["total_after"], review["total_before"]),
        "memory_snapshot": mem_info,
        "merged_matrix": merge_summary,
    }
    (REVIEW_DIR / "hard_qc_review_summary.json").write_text(
        json.dumps(machine_summary, indent=2) + "\n", encoding="utf-8"
    )
    print(f"[review] wrote outputs to {REVIEW_DIR}")


if __name__ == "__main__":
    main()
