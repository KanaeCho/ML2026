from __future__ import annotations

from dataclasses import asdict
from typing import cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .models import ComputedQcThresholds, QcMetricThresholdAudit, RunConfig


def _sum_axis(matrix, axis: int) -> np.ndarray:
    summed = matrix.sum(axis=axis)
    return np.asarray(summed).ravel().astype(float)


def _nnz_axis(matrix, axis: int) -> np.ndarray:
    if sparse.issparse(matrix):
        counts = matrix.getnnz(axis=axis)
    else:
        counts = np.count_nonzero(np.asarray(matrix), axis=axis)
    return np.asarray(counts).ravel().astype(int)


def _fractional_percent(counts: np.ndarray, totals: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        percent = np.divide(
            counts,
            totals,
            out=np.zeros_like(totals, dtype=float),
            where=totals > 0,
        )
    return percent * 100.0


def _median_and_mad(values: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return 0.0, 0.0
    center = float(np.median(finite))
    mad = float(np.median(np.abs(finite - center)))
    return center, mad


def _compute_lower_tail_mad_threshold_log10(
    values: np.ndarray,
    *,
    nmads: float,
    floor_min: int,
    floor_max: int,
    metric_name: str,
) -> tuple[int, QcMetricThresholdAudit, bool]:
    transformed = np.log10(np.asarray(values, dtype=float) + 1.0)
    center, mad = _median_and_mad(transformed)
    zero_mad = mad <= 0.0
    raw_threshold = center - nmads * mad if not zero_mad else center
    original_scale = max(float((10**raw_threshold) - 1.0), 0.0)
    final_threshold = int(np.ceil(np.clip(original_scale, floor_min, floor_max)))

    guardrails: list[str] = []
    if final_threshold == floor_min and original_scale < floor_min:
        guardrails.append("floor_min")
    if final_threshold == floor_max and original_scale > floor_max:
        guardrails.append("floor_max")
    if zero_mad:
        guardrails.append("zero_mad")

    audit = QcMetricThresholdAudit(
        transform="log10p1",
        direction="lower",
        center=center,
        mad=mad,
        nmads=nmads,
        raw_threshold=raw_threshold,
        final_threshold_original_scale=float(final_threshold),
        n_cells_used=int(np.isfinite(transformed).sum()),
        guardrails_applied=tuple(guardrails),
    )
    return final_threshold, audit, zero_mad


def _compute_upper_tail_mad_threshold(
    values: np.ndarray,
    *,
    nmads: float,
    min_bound: float,
    max_bound: float,
    metric_name: str,
) -> tuple[float, QcMetricThresholdAudit, bool]:
    finite = np.asarray(values, dtype=float)
    center, mad = _median_and_mad(finite)
    zero_mad = mad <= 0.0
    raw_threshold = center + nmads * mad if not zero_mad else center
    final_threshold = float(np.clip(raw_threshold, min_bound, max_bound))

    guardrails: list[str] = []
    if final_threshold == min_bound and raw_threshold < min_bound:
        guardrails.append("ceiling_min")
    if final_threshold == max_bound and raw_threshold > max_bound:
        guardrails.append("ceiling_max")
    if zero_mad:
        guardrails.append("zero_mad")

    audit = QcMetricThresholdAudit(
        transform="identity",
        direction="upper",
        center=center,
        mad=mad,
        nmads=nmads,
        raw_threshold=raw_threshold,
        final_threshold_original_scale=final_threshold,
        n_cells_used=int(np.isfinite(finite).sum()),
        guardrails_applied=tuple(guardrails),
    )
    return final_threshold, audit, zero_mad


def _compute_sample_qc_thresholds(
    obs: pd.DataFrame,
    config: RunConfig,
) -> ComputedQcThresholds:
    sample_id = ""
    if "sample_id" in obs.columns and len(obs) > 0:
        sample_id_series = obs["sample_id"]
        sample_id = str(sample_id_series.iloc[0])

    gse = ""
    if "gse" in obs.columns and len(obs) > 0:
        gse_series = obs["gse"]
        gse = str(gse_series.iloc[0])
    n_cells_total = int(len(obs))
    small_sample_rule_used = n_cells_total < config.qc.min_cells_for_dynamic

    counts_nmads = config.qc.counts_lower_nmads + (0.5 if small_sample_rule_used else 0.0)
    genes_nmads = config.qc.genes_lower_nmads + (0.5 if small_sample_rule_used else 0.0)
    mt_nmads = config.qc.pct_mt_upper_nmads + (0.5 if small_sample_rule_used else 0.0)
    ribo_nmads = config.qc.pct_ribo_upper_nmads + (0.5 if small_sample_rule_used else 0.0)

    min_counts, counts_audit, counts_zero_mad = _compute_lower_tail_mad_threshold_log10(
        obs["n_counts"].to_numpy(dtype=float),
        nmads=counts_nmads,
        floor_min=config.qc.count_floor_min,
        floor_max=config.qc.count_floor_max,
        metric_name="n_counts",
    )
    min_genes, genes_audit, genes_zero_mad = _compute_lower_tail_mad_threshold_log10(
        obs["n_genes"].to_numpy(dtype=float),
        nmads=genes_nmads,
        floor_min=config.qc.gene_floor_min,
        floor_max=config.qc.gene_floor_max,
        metric_name="n_genes",
    )
    max_pct_mt, mt_audit, mt_zero_mad = _compute_upper_tail_mad_threshold(
        obs["pct_mt"].to_numpy(dtype=float),
        nmads=mt_nmads,
        min_bound=config.qc.pct_mt_ceiling_min,
        max_bound=config.qc.pct_mt_ceiling_max,
        metric_name="pct_mt",
    )
    max_pct_ribo, ribo_audit, ribo_zero_mad = _compute_upper_tail_mad_threshold(
        obs["pct_ribo"].to_numpy(dtype=float),
        nmads=ribo_nmads,
        min_bound=config.qc.pct_ribo_ceiling_min,
        max_bound=config.qc.pct_ribo_ceiling_max,
        metric_name="pct_ribo",
    )

    zero_mad_metrics = tuple(
        metric_name
        for metric_name, is_zero in (
            ("n_counts", counts_zero_mad),
            ("n_genes", genes_zero_mad),
            ("pct_mt", mt_zero_mad),
            ("pct_ribo", ribo_zero_mad),
        )
        if is_zero
    )

    return ComputedQcThresholds(
        sample_id=sample_id,
        gse=gse,
        method=config.qc.method,
        n_cells_total=n_cells_total,
        min_counts=min_counts,
        min_genes=min_genes,
        max_pct_mt=max_pct_mt,
        max_pct_ribo=max_pct_ribo,
        n_counts_audit=counts_audit,
        n_genes_audit=genes_audit,
        pct_mt_audit=mt_audit,
        pct_ribo_audit=ribo_audit,
        small_sample_rule_used=small_sample_rule_used,
        zero_mad_metrics=zero_mad_metrics,
        retention_guard_triggered=False,
    )


def _json_safe_qc_thresholds(thresholds: ComputedQcThresholds) -> dict[str, object]:
    payload = asdict(thresholds)
    payload["zero_mad_metrics"] = list(payload.get("zero_mad_metrics", []))
    for metric_key in ("n_counts_audit", "n_genes_audit", "pct_mt_audit", "pct_ribo_audit"):
        metric_payload = dict(cast(dict[str, object], payload.get(metric_key, {})))
        metric_payload["guardrails_applied"] = list(metric_payload.get("guardrails_applied", []))
        payload[metric_key] = metric_payload
    return payload


def compute_qc_metrics(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    del config

    out = adata.copy()
    matrix = out.X
    var_names = pd.Index(out.var_names.astype(str))

    total_counts = _sum_axis(matrix, axis=1)
    n_genes = _nnz_axis(matrix, axis=1)

    mt_mask = np.asarray(var_names.str.startswith("MT-"), dtype=bool)
    ribo_mask = np.asarray(
        var_names.str.startswith("RPS") | var_names.str.startswith("RPL"),
        dtype=bool,
    )

    mt_counts = (
        _sum_axis(matrix[:, mt_mask], axis=1) if mt_mask.any() else np.zeros(out.n_obs)
    )
    ribo_counts = (
        _sum_axis(matrix[:, ribo_mask], axis=1)
        if ribo_mask.any()
        else np.zeros(out.n_obs)
    )

    out.obs["n_counts"] = total_counts
    out.obs["n_genes"] = n_genes
    out.obs["pct_mt"] = _fractional_percent(mt_counts, total_counts)
    out.obs["pct_ribo"] = _fractional_percent(ribo_counts, total_counts)
    return out


def apply_qc_filters(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    out = adata.copy()

    if "is_doublet" in out.obs:
        is_doublet = out.obs["is_doublet"].fillna(False).astype(bool)
    else:
        is_doublet = pd.Series(False, index=out.obs_names, dtype=bool)

    thresholds = _compute_sample_qc_thresholds(cast(pd.DataFrame, out.obs), config)
    out.uns["qc_thresholds"] = _json_safe_qc_thresholds(thresholds)

    fails_count_floor = out.obs["n_counts"] < thresholds.min_counts
    fails_gene_floor = out.obs["n_genes"] < thresholds.min_genes
    fails_mt_ceiling = out.obs["pct_mt"] > thresholds.max_pct_mt
    fails_ribo_ceiling = out.obs["pct_ribo"] > thresholds.max_pct_ribo
    fails_doublet = is_doublet

    out.obs["fails_count_floor"] = fails_count_floor.astype(bool)
    out.obs["fails_gene_floor"] = fails_gene_floor.astype(bool)
    out.obs["fails_mt_ceiling"] = fails_mt_ceiling.astype(bool)
    out.obs["fails_ribo_ceiling"] = fails_ribo_ceiling.astype(bool)
    out.obs["fails_doublet"] = fails_doublet.astype(bool)

    pass_qc = ~(
        out.obs["fails_count_floor"]
        | out.obs["fails_gene_floor"]
        | out.obs["fails_mt_ceiling"]
        | out.obs["fails_ribo_ceiling"]
        | out.obs["fails_doublet"]
    )
    out.obs["pass_qc"] = pass_qc.astype(bool)
    return out


__all__ = ["compute_qc_metrics", "apply_qc_filters"]
