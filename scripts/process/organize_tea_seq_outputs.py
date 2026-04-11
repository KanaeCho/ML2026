#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import shutil
import statistics
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"
PROCESS_DIR = ROOT / "scripts" / "process"
CLUSTER_RENDER_SCRIPT = PROCESS_DIR / "render_tea_seq_cima_cluster_umap.py"
ADT_RENDER_SCRIPT = PROCESS_DIR / "render_tea_seq_adt_broad_umap.py"

CORE_ARTIFACTS = (
    "GSM_PLACEHOLDER_seurat_qc.rds",
    "metadata.csv",
    "metadata_qc.csv",
    "qc_overview.png",
    "qc_summary.csv",
    "run_status.json",
    "validation_result.csv",
    "umap_atac_clusters.png",
    "matrix/matrix.mtx",
    "matrix/barcodes.tsv.gz",
    "matrix/features.tsv.gz",
)

ACCEPTED_ARTIFACTS = (
    "cima_cluster_centroid_labels.csv",
    "cima_cluster_centroid_label_summary.csv",
    "umap_atac_cima_cell_type_l1_cluster_centroid.png",
    "adt_cluster_broad_labels.csv",
    "adt_cluster_broad_label_summary.csv",
    "umap_atac_adt_cluster_broad_celltype.png",
)

QC_METRICS = (
    "input_cells",
    "pass_qc",
    "qc_rate",
    "median_TSS_enrichment",
    "median_FRiP",
    "median_fragments",
    "query_cluster_count",
    "cima_unique_l4_labels",
    "median_cima_l4_score",
)


def parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_float(value: str | None) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    if text.endswith("%"):
        text = text[:-1]
    try:
        return float(text)
    except ValueError:
        return None


def format_float(value: float | None, digits: int = 4) -> str:
    if value is None:
        return ""
    return f"{value:.{digits}f}"


def format_percent(value: float | None, digits: int = 2) -> str:
    if value is None:
        return ""
    return f"{value * 100:.{digits}f}%"


def read_metric_table(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        return {row["metric"]: row["value"] for row in csv.DictReader(handle)}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def discover_sample_dirs(dataset_dir: Path) -> list[Path]:
    return sorted(
        path
        for path in dataset_dir.iterdir()
        if path.is_dir() and path.name.startswith("GSM")
    )


def core_artifact_names(sample_id: str) -> list[str]:
    return [
        name.replace("GSM_PLACEHOLDER", sample_id) for name in CORE_ARTIFACTS
    ]


def relpath(path: Path, project_root: Path) -> str:
    try:
        return str(path.relative_to(project_root))
    except ValueError:
        return str(path)


def compute_cluster_purity_from_labels(rows: list[dict[str, str]]) -> tuple[float | None, float | None]:
    if not rows:
        return None, None

    label_col = None
    for candidate in ("cima_cell_type_l1", "cima_cell_type_l1_masked"):
        if candidate in rows[0]:
            label_col = candidate
            break
    if label_col is None or "seurat_clusters" not in rows[0]:
        return None, None

    cluster_purities: list[float] = []
    cluster_sizes: dict[str, int] = {}
    for cluster, cluster_rows in group_rows(rows, "seurat_clusters").items():
        labels = [
            row.get(label_col, "").strip()
            for row in cluster_rows
            if row.get(label_col, "").strip()
        ]
        if not labels:
            continue
        counts = Counter(labels)
        purity = counts.most_common(1)[0][1] / len(cluster_rows)
        cluster_purities.append(purity)
        cluster_sizes[cluster] = len(cluster_rows)

    if not cluster_purities:
        return None, None

    median_purity = statistics.median(cluster_purities)
    return median_purity, None


def group_rows(rows: list[dict[str, str]], key: str) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        grouped[str(row.get(key, ""))].append(row)
    return grouped


def summarize_validation_result(
    validation_path: Path,
    score_threshold: float,
    margin_threshold: float,
    low_purity_threshold: float = 0.8,
) -> dict[str, float | int | None]:
    rows = read_csv_rows(validation_path)
    if not rows:
        return {
            "cell_count": 0,
            "low_conf_cell_count": 0,
            "low_conf_cell_frac": None,
            "cluster_majority_purity_median": None,
            "low_purity_cluster_cell_frac": None,
        }

    has_explicit_low_conf = any(
        row.get("cima_l1_low_confidence", "").strip() for row in rows
    )
    low_conf_flags: list[bool] = []
    for row in rows:
        if has_explicit_low_conf:
            low_conf_flags.append(parse_bool(row.get("cima_l1_low_confidence")))
            continue
        score = parse_float(row.get("cima_l4_score"))
        margin = parse_float(row.get("cima_l4_score_margin"))
        low_score = score is not None and score < score_threshold
        low_margin = margin is not None and margin < margin_threshold
        low_conf_flags.append(low_score or low_margin)

    has_explicit_purity = any(
        row.get("cima_l1_cluster_purity", "").strip() for row in rows
    )
    cluster_purity_values: dict[str, float] = {}
    if has_explicit_purity:
        for row in rows:
            cluster = str(row.get("seurat_clusters", ""))
            purity = parse_float(row.get("cima_l1_cluster_purity"))
            if cluster and purity is not None:
                cluster_purity_values[cluster] = purity
    else:
        label_column = (
            "cima_cell_type_l1_masked"
            if "cima_cell_type_l1_masked" in rows[0]
            else "cima_cell_type_l1"
        )
        for cluster, cluster_rows in group_rows(rows, "seurat_clusters").items():
            labels = [
                row.get(label_column, "").strip()
                for row in cluster_rows
                if row.get(label_column, "").strip()
            ]
            if not labels:
                continue
            counts = Counter(labels)
            cluster_purity_values[cluster] = counts.most_common(1)[0][1] / len(
                cluster_rows
            )

    purity_values = list(cluster_purity_values.values())
    low_purity_cell_count = 0
    if cluster_purity_values:
        for row in rows:
            cluster = str(row.get("seurat_clusters", ""))
            purity = cluster_purity_values.get(cluster)
            if purity is not None and purity < low_purity_threshold:
                low_purity_cell_count += 1

    return {
        "cell_count": len(rows),
        "low_conf_cell_count": sum(low_conf_flags),
        "low_conf_cell_frac": sum(low_conf_flags) / len(rows),
        "cluster_majority_purity_median": (
            statistics.median(purity_values) if purity_values else None
        ),
        "low_purity_cluster_cell_frac": (
            low_purity_cell_count / len(rows) if rows else None
        ),
    }


def render_cluster_centroid_outputs(
    project_root: Path, gse: str, gsm: str, sample_dir: Path
) -> None:
    command = [
        sys.executable,
        str(CLUSTER_RENDER_SCRIPT),
        "--project-root",
        str(project_root),
        "--gse",
        gse,
        "--gsm",
        gsm,
        "--sample-dir",
        str(sample_dir),
        "--output-dir",
        str(sample_dir),
    ]
    subprocess.run(command, check=True, cwd=project_root)


def render_adt_broad_outputs(
    project_root: Path, gse: str, gsm: str, sample_dir: Path
) -> None:
    command = [
        sys.executable,
        str(ADT_RENDER_SCRIPT),
        "--project-root",
        str(project_root),
        "--gse",
        gse,
        "--gsm",
        gsm,
        "--sample-dir",
        str(sample_dir),
        "--output-dir",
        str(sample_dir),
    ]
    subprocess.run(command, check=True, cwd=project_root)


def record_manifest_row(
    rows: list[dict[str, str]],
    project_root: Path,
    path: Path,
    action: str,
    sample_id: str = "",
    destination: Path | None = None,
) -> None:
    rows.append(
        {
            "sample_id": sample_id,
            "relative_path": relpath(path, project_root),
            "action": action,
            "destination": relpath(destination, project_root) if destination else "",
        }
    )


def maybe_move_file(
    source: Path,
    destination: Path,
    *,
    dry_run: bool,
    manifest_rows: list[dict[str, str]],
    project_root: Path,
    sample_id: str,
) -> None:
    if destination.exists():
        record_manifest_row(
            manifest_rows,
            project_root,
            source,
            "already_present",
            sample_id=sample_id,
            destination=destination,
        )
        return
    if dry_run:
        record_manifest_row(
            manifest_rows,
            project_root,
            source,
            "would_move",
            sample_id=sample_id,
            destination=destination,
        )
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(source), str(destination))
    record_manifest_row(
        manifest_rows,
        project_root,
        destination,
        "moved",
        sample_id=sample_id,
        destination=destination,
    )


def ensure_sample_outputs(
    *,
    project_root: Path,
    gse: str,
    sample_dir: Path,
    dry_run: bool,
    manifest_rows: list[dict[str, str]],
) -> None:
    gsm = sample_dir.name
    legacy_l1_png = (
        sample_dir.parent
        / "l1"
        / f"{gsm}_umap_atac_cima_cell_type_l1_cluster_centroid.png"
    )
    cluster_png = sample_dir / "umap_atac_cima_cell_type_l1_cluster_centroid.png"
    if legacy_l1_png.exists():
        maybe_move_file(
            legacy_l1_png,
            cluster_png,
            dry_run=dry_run,
            manifest_rows=manifest_rows,
            project_root=project_root,
            sample_id=gsm,
        )

    cluster_artifacts = (
        sample_dir / "cima_cluster_centroid_labels.csv",
        sample_dir / "cima_cluster_centroid_label_summary.csv",
        cluster_png,
    )
    if not all(path.exists() for path in cluster_artifacts):
        if dry_run:
            record_manifest_row(
                manifest_rows,
                project_root,
                sample_dir,
                "would_render_cluster_centroid",
                sample_id=gsm,
            )
        else:
            render_cluster_centroid_outputs(project_root, gse, gsm, sample_dir)

    adt_artifacts = (
        sample_dir / "adt_cluster_broad_labels.csv",
        sample_dir / "adt_cluster_broad_label_summary.csv",
        sample_dir / "umap_atac_adt_cluster_broad_celltype.png",
    )
    if not all(path.exists() for path in adt_artifacts):
        if dry_run:
            record_manifest_row(
                manifest_rows,
                project_root,
                sample_dir,
                "would_render_adt_broad",
                sample_id=gsm,
            )
        else:
            render_adt_broad_outputs(project_root, gse, gsm, sample_dir)


def collect_sample_row(
    *,
    project_root: Path,
    sample_dir: Path,
    score_threshold: float,
    margin_threshold: float,
    low_purity_threshold: float,
) -> tuple[dict[str, str], dict[str, str], bool]:
    gsm = sample_dir.name
    qc_summary = read_metric_table(sample_dir / "qc_summary.csv")
    validation_summary = summarize_validation_result(
        sample_dir / "validation_result.csv",
        score_threshold=score_threshold,
        margin_threshold=margin_threshold,
        low_purity_threshold=low_purity_threshold,
    )

    core_names = core_artifact_names(gsm)
    missing_core = [name for name in core_names if not (sample_dir / name).exists()]
    missing_accepted = [
        name for name in ACCEPTED_ARTIFACTS if not (sample_dir / name).exists()
    ]

    sample_row = {
        "sample_id": gsm,
        "core_complete": str(not missing_core),
        "accepted_complete": str(not missing_accepted),
        "missing_core_count": str(len(missing_core)),
        "missing_accepted_count": str(len(missing_accepted)),
        "missing_core": ";".join(missing_core),
        "missing_accepted": ";".join(missing_accepted),
    }
    for metric in QC_METRICS:
        sample_row[metric] = qc_summary.get(metric, "")

    validation_row = {
        "sample_id": gsm,
        "validation_cells": str(validation_summary["cell_count"]),
        "low_conf_cell_count": str(validation_summary["low_conf_cell_count"]),
        "low_conf_cell_frac": format_percent(
            validation_summary["low_conf_cell_frac"]
        ),
        "cluster_majority_purity_median": format_float(
            validation_summary["cluster_majority_purity_median"]
        ),
        "low_purity_cluster_cell_frac": format_percent(
            validation_summary["low_purity_cluster_cell_frac"]
        ),
    }
    return sample_row, validation_row, not missing_accepted


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def sort_fraction_rows(rows: list[dict[str, str]], key: str) -> list[dict[str, str]]:
    def rank(row: dict[str, str]) -> float:
        value = parse_float(row.get(key))
        return -1.0 if value is None else value

    return sorted(rows, key=rank, reverse=True)


def write_markdown_summary(
    *,
    summary_path: Path,
    gse: str,
    sample_rows: list[dict[str, str]],
    validation_rows: list[dict[str, str]],
) -> None:
    sample_count = len(sample_rows)
    accepted_complete = sum(
        row["accepted_complete"].lower() == "true" for row in sample_rows
    )
    qc_rates = [parse_float(row.get("qc_rate")) for row in sample_rows]
    qc_rates = [value for value in qc_rates if value is not None]
    tss_values = [parse_float(row.get("median_TSS_enrichment")) for row in sample_rows]
    tss_values = [value for value in tss_values if value is not None]
    frip_values = [parse_float(row.get("median_FRiP")) for row in sample_rows]
    frip_values = [value for value in frip_values if value is not None]

    lines = [
        f"# {gse} TEA-seq QC Audit",
        "",
        f"- Samples audited: {sample_count}",
        f"- Accepted TEA-seq outputs complete: {accepted_complete}/{sample_count}",
    ]
    if qc_rates:
        lines.append(
            f"- QC rate range: {min(qc_rates):.2f}% to {max(qc_rates):.2f}%"
        )
    if tss_values:
        lines.append(
            f"- Median TSS enrichment range: {min(tss_values):.2f} to {max(tss_values):.2f}"
        )
    if frip_values:
        lines.append(
            f"- Median FRiP range: {min(frip_values):.4f} to {max(frip_values):.4f}"
        )

    low_conf_sorted = sort_fraction_rows(validation_rows, "low_conf_cell_frac")
    lines.extend(["", "## Priority Review Samples", ""])
    if not low_conf_sorted:
        lines.append("- No validation rows available.")
    else:
        for row in low_conf_sorted[:5]:
            lines.append(
                "- "
                + f"{row['sample_id']}: low_conf={row['low_conf_cell_frac'] or 'NA'}, "
                + f"cluster_purity_median={row['cluster_majority_purity_median'] or 'NA'}, "
                + f"low_purity_cluster_cells={row['low_purity_cluster_cell_frac'] or 'NA'}"
            )

    missing_sorted = sorted(
        sample_rows, key=lambda row: int(row["missing_accepted_count"]), reverse=True
    )
    lines.extend(["", "## Missing Accepted Artifacts", ""])
    for row in missing_sorted[:5]:
        if row["missing_accepted_count"] == "0":
            continue
        lines.append(
            f"- {row['sample_id']}: {row['missing_accepted_count']} missing -> {row['missing_accepted']}"
        )
    if lines[-1] == "## Missing Accepted Artifacts":
        lines.append("- None.")

    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def enumerate_legacy_files(dataset_dir: Path) -> list[Path]:
    legacy_files: list[Path] = []
    for relative_dir in ("l1", "gex_refined_cima_compare"):
        legacy_root = dataset_dir / relative_dir
        if not legacy_root.exists():
            continue
        legacy_files.extend(sorted(path for path in legacy_root.rglob("*") if path.is_file()))
    return legacy_files


def mark_legacy_retention(
    *,
    project_root: Path,
    dataset_dir: Path,
    manifest_rows: list[dict[str, str]],
    action: str,
) -> None:
    for path in enumerate_legacy_files(dataset_dir):
        record_manifest_row(manifest_rows, project_root, path, action)


def maybe_delete_legacy_dirs(
    *,
    project_root: Path,
    dataset_dir: Path,
    delete_legacy: bool,
    dry_run: bool,
    accepted_complete: bool,
    manifest_rows: list[dict[str, str]],
) -> None:
    if not delete_legacy:
        mark_legacy_retention(
            project_root=project_root,
            dataset_dir=dataset_dir,
            manifest_rows=manifest_rows,
            action="kept_legacy",
        )
        return

    if not accepted_complete:
        mark_legacy_retention(
            project_root=project_root,
            dataset_dir=dataset_dir,
            manifest_rows=manifest_rows,
            action="kept_incomplete",
        )
        return

    for relative_dir in ("gex_refined_cima_compare", "l1"):
        legacy_root = dataset_dir / relative_dir
        if not legacy_root.exists():
            continue
        for path in sorted(path for path in legacy_root.rglob("*") if path.is_file()):
            record_manifest_row(
                manifest_rows,
                project_root,
                path,
                "would_delete" if dry_run else "deleted",
            )
        if not dry_run:
            shutil.rmtree(legacy_root)


def organize_dataset(
    *,
    project_root: Path,
    gse: str,
    dry_run: bool,
    delete_legacy: bool,
    score_threshold: float,
    margin_threshold: float,
    low_purity_threshold: float = 0.8,
) -> dict[str, object]:
    dataset_dir = project_root / "output" / gse
    if not dataset_dir.exists():
        raise FileNotFoundError(f"Dataset output directory not found: {dataset_dir}")

    manifest_rows: list[dict[str, str]] = []
    sample_rows: list[dict[str, str]] = []
    validation_rows: list[dict[str, str]] = []

    for sample_dir in discover_sample_dirs(dataset_dir):
        ensure_sample_outputs(
            project_root=project_root,
            gse=gse,
            sample_dir=sample_dir,
            dry_run=dry_run,
            manifest_rows=manifest_rows,
        )

    accepted_complete = True
    for sample_dir in discover_sample_dirs(dataset_dir):
        sample_row, validation_row, sample_complete = collect_sample_row(
            project_root=project_root,
            sample_dir=sample_dir,
            score_threshold=score_threshold,
            margin_threshold=margin_threshold,
            low_purity_threshold=low_purity_threshold,
        )
        sample_rows.append(sample_row)
        validation_rows.append(validation_row)
        accepted_complete = accepted_complete and sample_complete

    maybe_delete_legacy_dirs(
        project_root=project_root,
        dataset_dir=dataset_dir,
        delete_legacy=delete_legacy,
        dry_run=dry_run,
        accepted_complete=accepted_complete,
        manifest_rows=manifest_rows,
    )

    audit_dir = dataset_dir / "qc_audit"
    sample_table_path = audit_dir / "tea_seq_sample_audit.csv"
    validation_table_path = audit_dir / "tea_seq_validation_audit.csv"
    manifest_path = audit_dir / "tea_seq_legacy_manifest.csv"
    summary_path = audit_dir / "tea_seq_qc_summary.md"

    write_csv(sample_table_path, sample_rows)
    write_csv(validation_table_path, validation_rows)
    write_csv(manifest_path, manifest_rows)
    write_markdown_summary(
        summary_path=summary_path,
        gse=gse,
        sample_rows=sample_rows,
        validation_rows=validation_rows,
    )

    return {
        "dataset_dir": dataset_dir,
        "audit_dir": audit_dir,
        "sample_rows": sample_rows,
        "validation_rows": validation_rows,
        "legacy_manifest_rows": manifest_rows,
        "accepted_complete": accepted_complete,
        "sample_table_path": sample_table_path,
        "validation_table_path": validation_table_path,
        "manifest_path": manifest_path,
        "summary_path": summary_path,
    }


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Organize TEA-seq accepted outputs and produce a dataset-level QC audit"
    )
    parser.add_argument("--gse", required=True)
    parser.add_argument("--project-root", default=str(ROOT))
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--delete-legacy", action="store_true")
    parser.add_argument("--score-threshold", type=float, default=0.4)
    parser.add_argument("--margin-threshold", type=float, default=0.1)
    parser.add_argument("--low-purity-threshold", type=float, default=0.8)
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    result = organize_dataset(
        project_root=Path(args.project_root).resolve(),
        gse=args.gse,
        dry_run=args.dry_run,
        delete_legacy=args.delete_legacy,
        score_threshold=args.score_threshold,
        margin_threshold=args.margin_threshold,
        low_purity_threshold=args.low_purity_threshold,
    )

    print(f"Audit dir: {result['audit_dir']}")
    print(f"Sample audit: {result['sample_table_path']}")
    print(f"Validation audit: {result['validation_table_path']}")
    print(f"Legacy manifest: {result['manifest_path']}")
    print(f"Summary: {result['summary_path']}")
    print(f"Accepted complete: {result['accepted_complete']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
