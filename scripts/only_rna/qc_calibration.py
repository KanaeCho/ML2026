from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
from typing import Any, cast

import pandas as pd

from .annotation import annotate_with_all_versions
from .config import load_default_config, merge_cli_overrides
from .discovery import DiscoveredSample, ROOT, discover_rna_samples, resolve_data_root, selected_rna_gses
from .doublet import run_doublet_detection
from .embedding import run_embedding
from .outputs import write_sample_outputs
from .qc import apply_qc_filters, compute_qc_metrics
from .read_inputs import read_sample_input
from .tuning_metrics import score_annotation_metrics


DEFAULT_CONFIG_PATH = ROOT / "scripts" / "only_rna" / "default_config.yaml"
DEFAULT_OUTPUT_ROOT = ROOT / "output" / "rna_qc_calibration"


def _parse_sample_token(token: str) -> tuple[str, str]:
    normalized = token.strip()
    for sep in (":", "/"):
        if sep in normalized:
            gse, sample_id = normalized.split(sep, 1)
            return gse.strip(), sample_id.strip()
    raise ValueError(f"Invalid sample token '{token}'. Use GSE:SAMPLE or GSE/SAMPLE")


def _discover_supported_sample(gse: str, sample_id: str) -> DiscoveredSample:
    data_root = resolve_data_root(ROOT)
    selected = selected_rna_gses(data_root / "reference")
    for sample in discover_rna_samples(data_root / "raw", selected):
        if sample.gse == gse and sample.sample_id == sample_id:
            if not sample.supported:
                raise ValueError(f"Unsupported sample: {gse}/{sample_id}: {sample.note}")
            if sample.sample_kind == "gse_shared":
                raise ValueError(
                    f"Calibration currently targets explicit GSM samples only: {gse}/{sample_id}"
                )
            return sample
    raise FileNotFoundError(f"RNA sample not found: {gse}/{sample_id}")


def _profile_configs():
    mainline = load_default_config(DEFAULT_CONFIG_PATH)
    baseline = merge_cli_overrides(
        mainline,
        qc__counts_lower_nmads=3.0,
        qc__genes_lower_nmads=3.0,
        qc__pct_mt_upper_nmads=3.0,
        qc__pct_ribo_upper_nmads=3.5,
    )
    stricter_v1 = merge_cli_overrides(
        mainline,
        qc__counts_lower_nmads=2.5,
        qc__genes_lower_nmads=2.5,
        qc__pct_mt_upper_nmads=2.5,
    )
    return {
        "baseline": baseline,
        "stricter_v1": stricter_v1,
    }


def _run_profile(
    *,
    sample: DiscoveredSample,
    profile_name: str,
    output_root: Path,
):
    profile_config = _profile_configs()[profile_name]
    adata = read_sample_input(sample)
    adata = compute_qc_metrics(adata, profile_config)
    adata = run_doublet_detection(adata, profile_config)
    adata = apply_qc_filters(adata, profile_config)
    adata = run_embedding(adata, profile_config)
    adata = annotate_with_all_versions(
        adata,
        reference_dir=resolve_data_root(ROOT) / "reference",
        azimuth_model_dir=None,
        celltypist_model_path=None,
        singler_model_path=None,
        scanvi_model_path=None,
        methods=(profile_config.annotation.methods if profile_config.annotation else ["azimuth"]),
        azimuth_config=profile_config.azimuth,
    )

    profile_output_root = output_root / profile_name
    sample_output_dir = write_sample_outputs(
        adata,
        output_root=profile_output_root,
        gse=sample.gse,
        sample_id=sample.sample_id,
        config=profile_config,
    )
    return adata, profile_config, sample_output_dir


def _numeric_mean(frame: pd.DataFrame, column: str) -> float:
    if column not in frame.columns or len(frame) == 0:
        return 0.0
    numeric = pd.Series(pd.to_numeric(frame[column], errors="coerce"), dtype=float)
    values = numeric.dropna()
    if len(values) == 0:
        return 0.0
    return float(values.mean())


def _cluster_label_purity(frame: pd.DataFrame, label_key: str) -> float:
    if len(frame) == 0 or "cluster" not in frame.columns or label_key not in frame.columns:
        return 0.0

    working = pd.DataFrame(
        {
            "cluster": frame["cluster"].astype("string"),
            label_key: frame[label_key].astype("string"),
        }
    )
    mask = working["cluster"].notna() & working[label_key].notna()
    working = working.loc[mask].copy()
    if len(working) == 0:
        return 0.0

    dominant_total = 0
    total = 0
    for _, group in working.groupby("cluster", sort=False):
        counts = group[label_key].value_counts(dropna=True)
        if len(counts) == 0:
            continue
        dominant_total += int(counts.iloc[0])
        total += int(len(group))
    if total <= 0:
        return 0.0
    return float(dominant_total) / float(total)


def _summarize_profile(
    *,
    sample: DiscoveredSample,
    profile_name: str,
    sample_output_dir: Path,
    adata,
    profile_config,
) -> dict[str, object]:
    obs = adata.obs.copy()
    pass_qc_mask = obs["pass_qc"].fillna(False).astype(bool)
    pass_qc_obs = obs.loc[pass_qc_mask]
    qc_thresholds = dict(adata.uns.get("qc_thresholds", {}))
    annotation_status = dict(adata.uns.get("annotation_method_status", {})).get("azimuth", {})
    azimuth_status = str(annotation_status.get("status", ""))

    low_conf_fraction = 0.0
    if "azimuth_low_confidence" in pass_qc_obs.columns and len(pass_qc_obs) > 0:
        low_conf_fraction = float(
            pass_qc_obs["azimuth_low_confidence"].fillna(False).astype(bool).mean()
        )

    azimuth_score_mean = _numeric_mean(pass_qc_obs, "azimuth_score")
    azimuth_score_margin_mean = _numeric_mean(pass_qc_obs, "azimuth_score_margin")
    annotation_score = score_annotation_metrics(
        method_status=azimuth_status,
        confidence_mean=azimuth_score_mean,
        low_confidence_fraction=low_conf_fraction,
    )

    cluster_count = 0
    if "cluster" in pass_qc_obs.columns:
        cluster_count = int(pass_qc_obs["cluster"].dropna().astype(str).nunique())

    azimuth_label_count = 0
    if "azimuth_cell_type" in pass_qc_obs.columns:
        azimuth_label_count = int(
            pass_qc_obs["azimuth_cell_type"].dropna().astype(str).nunique()
        )

    cluster_purity_pbmcref = _cluster_label_purity(pass_qc_obs, "azimuth_cell_type")
    cluster_purity_cima_l1 = _cluster_label_purity(pass_qc_obs, "azimuth_cima_l1")

    cima_l1_unknown_fraction = 0.0
    if "azimuth_cima_l1" in pass_qc_obs.columns and len(pass_qc_obs) > 0:
        cima_l1_unknown_fraction = float(
            pass_qc_obs["azimuth_cima_l1"]
            .fillna("Unknown")
            .astype(str)
            .str.strip()
            .str.lower()
            .eq("unknown")
            .mean()
        )

    return {
        "gse": sample.gse,
        "sample_id": sample.sample_id,
        "profile": profile_name,
        "output_dir": str(sample_output_dir),
        "n_cells_total": int(adata.n_obs),
        "n_cells_pass_qc": int(pass_qc_mask.sum()),
        "pass_qc_fraction": float(pass_qc_mask.mean()) if len(pass_qc_mask) else 0.0,
        "azimuth_status": azimuth_status,
        "azimuth_detail": str(annotation_status.get("detail", "")),
        "azimuth_score_mean": azimuth_score_mean,
        "azimuth_score_margin_mean": azimuth_score_margin_mean,
        "azimuth_low_confidence_fraction": low_conf_fraction,
        "annotation_score": annotation_score,
        "cluster_count": cluster_count,
        "azimuth_label_count": azimuth_label_count,
        "cluster_purity_pbmcref": cluster_purity_pbmcref,
        "cluster_purity_cima_l1": cluster_purity_cima_l1,
        "azimuth_cima_l1_unknown_fraction": cima_l1_unknown_fraction,
        "counts_lower_nmads": float(profile_config.qc.counts_lower_nmads),
        "genes_lower_nmads": float(profile_config.qc.genes_lower_nmads),
        "pct_mt_upper_nmads": float(profile_config.qc.pct_mt_upper_nmads),
        "pct_ribo_upper_nmads": float(profile_config.qc.pct_ribo_upper_nmads),
        "final_min_counts": qc_thresholds.get("min_counts", ""),
        "final_min_genes": qc_thresholds.get("min_genes", ""),
        "final_max_pct_mt": qc_thresholds.get("max_pct_mt", ""),
        "final_max_pct_ribo": qc_thresholds.get("max_pct_ribo", ""),
        "small_sample_rule_used": bool(qc_thresholds.get("small_sample_rule_used", False)),
        "retention_guard_triggered": bool(
            qc_thresholds.get("retention_guard_triggered", False)
        ),
        "zero_mad_metrics": ";".join(qc_thresholds.get("zero_mad_metrics", [])),
    }


def _build_comparison_rows(summary_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    frame = pd.DataFrame(summary_rows)
    rows: list[dict[str, object]] = []
    grouped = frame.groupby(["gse", "sample_id"], sort=True)
    for group_key, group in grouped:
        gse, sample_id = cast(tuple[str, str], group_key)
        by_profile = {
            str(row["profile"]): cast(dict[str, Any], row.to_dict())
            for _, row in group.iterrows()
        }
        baseline = by_profile.get("baseline")
        stricter = by_profile.get("stricter_v1")
        if baseline is None or stricter is None:
            continue
        rows.append(
            {
                "gse": gse,
                "sample_id": sample_id,
                "baseline_pass_qc_fraction": float(baseline["pass_qc_fraction"]),
                "stricter_v1_pass_qc_fraction": float(stricter["pass_qc_fraction"]),
                "delta_pass_qc_fraction": float(stricter["pass_qc_fraction"])
                - float(baseline["pass_qc_fraction"]),
                "baseline_annotation_score": float(baseline["annotation_score"]),
                "stricter_v1_annotation_score": float(stricter["annotation_score"]),
                "delta_annotation_score": float(stricter["annotation_score"])
                - float(baseline["annotation_score"]),
                "baseline_azimuth_score_mean": float(baseline["azimuth_score_mean"]),
                "stricter_v1_azimuth_score_mean": float(stricter["azimuth_score_mean"]),
                "delta_azimuth_score_mean": float(stricter["azimuth_score_mean"])
                - float(baseline["azimuth_score_mean"]),
                "baseline_azimuth_score_margin_mean": float(
                    baseline["azimuth_score_margin_mean"]
                ),
                "stricter_v1_azimuth_score_margin_mean": float(
                    stricter["azimuth_score_margin_mean"]
                ),
                "delta_azimuth_score_margin_mean": float(
                    stricter["azimuth_score_margin_mean"]
                )
                - float(baseline["azimuth_score_margin_mean"]),
                "baseline_azimuth_low_confidence_fraction": float(
                    baseline["azimuth_low_confidence_fraction"]
                ),
                "stricter_v1_azimuth_low_confidence_fraction": float(
                    stricter["azimuth_low_confidence_fraction"]
                ),
                "delta_azimuth_low_confidence_fraction": float(
                    stricter["azimuth_low_confidence_fraction"]
                )
                - float(baseline["azimuth_low_confidence_fraction"]),
                "baseline_cluster_purity_pbmcref": float(baseline["cluster_purity_pbmcref"]),
                "stricter_v1_cluster_purity_pbmcref": float(
                    stricter["cluster_purity_pbmcref"]
                ),
                "delta_cluster_purity_pbmcref": float(
                    stricter["cluster_purity_pbmcref"]
                )
                - float(baseline["cluster_purity_pbmcref"]),
                "baseline_cluster_purity_cima_l1": float(baseline["cluster_purity_cima_l1"]),
                "stricter_v1_cluster_purity_cima_l1": float(
                    stricter["cluster_purity_cima_l1"]
                ),
                "delta_cluster_purity_cima_l1": float(
                    stricter["cluster_purity_cima_l1"]
                )
                - float(baseline["cluster_purity_cima_l1"]),
                "baseline_cluster_count": int(baseline["cluster_count"]),
                "stricter_v1_cluster_count": int(stricter["cluster_count"]),
                "delta_cluster_count": int(stricter["cluster_count"]) - int(baseline["cluster_count"]),
                "baseline_final_min_counts": baseline["final_min_counts"],
                "stricter_v1_final_min_counts": stricter["final_min_counts"],
                "baseline_final_min_genes": baseline["final_min_genes"],
                "stricter_v1_final_min_genes": stricter["final_min_genes"],
                "baseline_final_max_pct_mt": baseline["final_max_pct_mt"],
                "stricter_v1_final_max_pct_mt": stricter["final_max_pct_mt"],
                "baseline_dir": baseline["output_dir"],
                "stricter_v1_dir": stricter["output_dir"],
            }
        )
    return rows


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run baseline vs stricter QC calibration for selected RNA samples"
    )
    parser.add_argument(
        "--sample",
        dest="samples",
        action="append",
        required=True,
        help="Sample token as GSE:SAMPLE or GSE/SAMPLE. Repeat for multiple samples.",
    )
    parser.add_argument(
        "--output-root",
        default=str(DEFAULT_OUTPUT_ROOT),
        help="Calibration output root. Profiles are written under this directory.",
    )
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    output_root = Path(args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    requested = [_parse_sample_token(token) for token in args.samples]
    summary_rows: list[dict[str, object]] = []

    for gse, sample_id in requested:
        sample = _discover_supported_sample(gse, sample_id)
        for profile_name in ("baseline", "stricter_v1"):
            print(f"[qc-calibration] {profile_name} {gse}/{sample_id}")
            adata, profile_config, sample_output_dir = _run_profile(
                sample=sample,
                profile_name=profile_name,
                output_root=output_root,
            )
            summary_rows.append(
                _summarize_profile(
                    sample=sample,
                    profile_name=profile_name,
                    sample_output_dir=sample_output_dir,
                    adata=adata,
                    profile_config=profile_config,
                )
            )

    summary_path = output_root / "qc_profile_summary.csv"
    pd.DataFrame(summary_rows).sort_values(["gse", "sample_id", "profile"]).to_csv(
        summary_path, index=False
    )

    comparison_rows = _build_comparison_rows(summary_rows)
    comparison_path = output_root / "qc_profile_comparison.csv"
    pd.DataFrame(comparison_rows).sort_values(["gse", "sample_id"]).to_csv(
        comparison_path, index=False
    )

    manifest_path = output_root / "qc_profile_manifest.json"
    manifest_path.write_text(
        json.dumps(
            {
                "profiles": {
                    "baseline": asdict(_profile_configs()["baseline"].qc),
                    "stricter_v1": asdict(_profile_configs()["stricter_v1"].qc),
                },
                "samples": [f"{gse}/{sample_id}" for gse, sample_id in requested],
                "summary_csv": str(summary_path),
                "comparison_csv": str(comparison_path),
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"[qc-calibration] summary={summary_path}")
    print(f"[qc-calibration] comparison={comparison_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
