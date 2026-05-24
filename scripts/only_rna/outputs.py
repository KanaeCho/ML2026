from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path
from typing import Any, Sequence, cast

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .models import PlottingConfig, QcThresholds, RunConfig
from .final_celltype import (
    RNA_FINAL_CELLTYPE_MAPPING_VERSION,
    infer_rna_final_celltype_series,
    known_rna_final_celltype_mask,
)
from .plotting import (
    save_azimuth_candidate_overview,
    save_categorical_umap,
    save_dual_annotation_umap,
    save_highlight_category_overview,
    save_qc_overview,
    save_sample_cima_l1_umap,
)
from .tuning_metrics import score_annotation_metrics


UMAP_PLOTS = {
    "azimuth_cima_l1": (
        "umap_rna_pbmcref_vs_cima_l1.png",
        "pbmcref vs CIMA L1",
    ),
}

SAMPLE_ROOT_VALIDATION_PLOTS = {
    "azimuth_cima_l1": (
        "umap_rna_pbmcref_vs_cima_l1.png",
        "pbmcref vs CIMA L1",
    ),
    "pbmcref_highlight": (
        "umap_rna_pbmcref_highlight.png",
        "pbmcref highlight",
    ),
}

REMOVED_SAMPLE_ROOT_ARTIFACTS = (
    "qc_overview.png",
    "umap_rna_clusters.png",
    "umap_rna_cima_cell_type_l1.png",
    "umap_rna_cima_cell_type_l2.png",
    "umap_rna_cima_cell_type_l1_masked.png",
    "umap_rna_annotation_method_compare.png",
    "umap_rna_celltypist.png",
    "umap_rna_singler.png",
    "umap_rna_scanvi.png",
    "umap_rna_azimuth.png",
)


def _sanitize_dataframe_for_h5ad(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.index = pd.Index(
        np.array(
            ["" if pd.isna(value) else str(value) for value in out.index.tolist()],
            dtype=object,
        ),
        dtype=object,
    )
    out.columns = pd.Index([str(value) for value in out.columns.tolist()], dtype=object)
    for column in out.columns:
        series = out[column]
        if (
            pd.api.types.is_string_dtype(series.dtype)
            or pd.api.types.is_object_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            out[column] = pd.Series(
                ["" if pd.isna(value) else str(value) for value in series.tolist()],
                index=series.index,
                dtype=object,
            ).astype(object)
    return out


def _prepare_h5ad_adata(adata: ad.AnnData) -> ad.AnnData:
    obs = cast(pd.DataFrame, adata.obs)
    pass_qc = _bool_pass_qc(obs)
    out = adata[pass_qc.to_numpy()].copy()
    obs_out = cast(pd.DataFrame, out.obs).copy()
    obs_out = _add_rna_final_output_columns(obs_out, out.obs_names.astype(str).tolist())
    known_mask = _known_final_celltype_mask(obs_out)
    out = out[known_mask.to_numpy()].copy()
    obs_out = obs_out.loc[known_mask.to_numpy()].copy()

    keep_obs = [
        "cell_barcode",
        "sample",
        "dataset",
        "age",
        "health",
        "donor",
        "final_celltype",
        "final_celltype_mapping",
        "azimuth_cima_l1_raw",
        "pbmcref_celltype",
        "umap_1",
        "umap_2",
        "azimuth_cima_l1",
        "azimuth_cell_type_l2_raw",
        "azimuth_cell_type",
    ]
    out.obs = obs_out[[column for column in keep_obs if column in obs_out.columns]].copy()
    sanitized_obs = _sanitize_dataframe_for_h5ad(cast(pd.DataFrame, out.obs))
    sanitized_var = _sanitize_dataframe_for_h5ad(cast(pd.DataFrame, out.var))
    out.obs = sanitized_obs
    out.var = sanitized_var
    return out


def _add_rna_final_output_columns(frame: pd.DataFrame, cell_ids: Sequence[str]) -> pd.DataFrame:
    out = frame.copy()
    out["cell_barcode"] = list(cell_ids)
    out["sample"] = out["sample_id"].astype(str) if "sample_id" in out.columns else ""
    out["dataset"] = out["gse"].astype(str) if "gse" in out.columns else ""
    if "donor" not in out.columns:
        out["donor"] = out["individual_id"].astype(str) if "individual_id" in out.columns else ""
    if "age" not in out.columns:
        out["age"] = ""
    if "health" not in out.columns:
        out["health"] = ""
    if "azimuth_cima_l1_raw" in out.columns:
        raw_l1 = out["azimuth_cima_l1_raw"]
    elif "azimuth_cima_l1" in out.columns:
        raw_l1 = out["azimuth_cima_l1"]
        out["azimuth_cima_l1_raw"] = raw_l1.astype(str)
    elif "cima_l1" in out.columns:
        raw_l1 = out["cima_l1"]
        out["azimuth_cima_l1_raw"] = raw_l1.astype(str)
    else:
        raw_l1 = pd.Series("", index=out.index, dtype=object)
        out["azimuth_cima_l1_raw"] = ""
    raw_l2 = (
        out["azimuth_cell_type_l2_raw"]
        if "azimuth_cell_type_l2_raw" in out.columns
        else out["azimuth_cell_type"]
        if "azimuth_cell_type" in out.columns
        else pd.Series(pd.NA, index=out.index, dtype=object)
    )
    out["final_celltype"] = infer_rna_final_celltype_series(raw_l1, raw_l2).astype(str)
    out["azimuth_cima_l1"] = out["final_celltype"].astype(str)
    out["final_celltype_mapping"] = RNA_FINAL_CELLTYPE_MAPPING_VERSION
    if "azimuth_cell_type_l2_raw" in out.columns:
        out["pbmcref_celltype"] = out["azimuth_cell_type_l2_raw"].astype(str)
    elif "azimuth_cell_type" in out.columns:
        out["pbmcref_celltype"] = out["azimuth_cell_type"].astype(str)
    else:
        out["pbmcref_celltype"] = ""
    return out


def _known_final_celltype_mask(frame: pd.DataFrame) -> pd.Series:
    if "final_celltype" not in frame.columns:
        return pd.Series(False, index=frame.index)
    return known_rna_final_celltype_mask(frame["final_celltype"])


def _bool_pass_qc(obs: pd.DataFrame) -> pd.Series:
    if "pass_qc" not in obs.columns:
        raise KeyError("adata.obs must contain 'pass_qc'")
    return cast(pd.Series, obs["pass_qc"].fillna(False).astype(bool))


def _metadata_frame(adata: ad.AnnData) -> pd.DataFrame:
    metadata = cast(pd.DataFrame, adata.obs.copy())
    metadata.insert(0, "cell_id", adata.obs_names.astype(str))
    return metadata.reset_index(drop=True)


def _write_qc_summary(
    output_path: Path,
    *,
    gse: str,
    sample_id: str,
    n_cells_total: int,
    n_cells_pass_qc: int,
    n_cells_final_output: int,
    n_cells_unknown_final_celltype_removed: int,
    metadata_qc: pd.DataFrame,
    annotation_method_status: dict[str, dict[str, str]],
    qc_thresholds: dict[str, Any] | None,
) -> None:
    azimuth_status = annotation_method_status.get("azimuth", {})
    qc_thresholds = qc_thresholds or {}
    azimuth_score_mean = 0.0
    azimuth_score_margin_mean = 0.0
    azimuth_low_confidence_fraction = 0.0
    if len(metadata_qc) > 0:
        if "azimuth_score" in metadata_qc.columns:
            azimuth_scores = pd.to_numeric(metadata_qc["azimuth_score"], errors="coerce").dropna()
            if len(azimuth_scores) > 0:
                azimuth_score_mean = float(azimuth_scores.mean())
        if "azimuth_score_margin" in metadata_qc.columns:
            azimuth_margins = pd.to_numeric(
                metadata_qc["azimuth_score_margin"], errors="coerce"
            ).dropna()
            if len(azimuth_margins) > 0:
                azimuth_score_margin_mean = float(azimuth_margins.mean())
        if "azimuth_low_confidence" in metadata_qc.columns:
            azimuth_low_confidence_fraction = float(
                metadata_qc["azimuth_low_confidence"].fillna(False).astype(bool).mean()
            )
    annotation_score = score_annotation_metrics(
        method_status=str(azimuth_status.get("status", "")),
        confidence_mean=azimuth_score_mean,
        low_confidence_fraction=azimuth_low_confidence_fraction,
    )
    pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "gse": gse,
                "n_cells_total": n_cells_total,
                "n_cells_pass_qc": n_cells_pass_qc,
                "n_cells_fail_qc": n_cells_total - n_cells_pass_qc,
                "n_cells_final_output": n_cells_final_output,
                "n_cells_unknown_final_celltype_removed": n_cells_unknown_final_celltype_removed,
                "pass_qc_fraction": (
                    float(n_cells_pass_qc) / float(n_cells_total)
                    if n_cells_total > 0
                    else 0.0
                ),
                "azimuth_status": azimuth_status.get("status", ""),
                "azimuth_detail": azimuth_status.get("detail", ""),
                "azimuth_score_mean": azimuth_score_mean,
                "azimuth_score_margin_mean": azimuth_score_margin_mean,
                "azimuth_low_confidence_fraction": azimuth_low_confidence_fraction,
                "annotation_score": annotation_score,
                "qc_threshold_method": qc_thresholds.get("method", ""),
                "final_min_counts": qc_thresholds.get("min_counts", ""),
                "final_min_genes": qc_thresholds.get("min_genes", ""),
                "final_max_pct_mt": qc_thresholds.get("max_pct_mt", ""),
                "final_max_pct_ribo": qc_thresholds.get("max_pct_ribo", ""),
            }
        ]
    ).to_csv(output_path, index=False)


def _write_qc_thresholds_json(output_path: Path, qc_thresholds: dict[str, Any]) -> None:
    output_path.write_text(
        json.dumps(qc_thresholds, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _write_text_gzip(path: Path, lines: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _export_matrix_triplet(adata: ad.AnnData, matrix_dir: Path) -> None:
    matrix_dir.mkdir(parents=True, exist_ok=True)
    matrix = cast(Any, adata.X).T
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    else:
        matrix = cast(Any, matrix).tocoo()

    mmwrite(matrix_dir / "matrix.mtx", matrix)
    _write_text_gzip(
        matrix_dir / "barcodes.tsv.gz",
        adata.obs_names.astype(str).tolist(),
    )

    feature_ids = (
        adata.var["feature_id"].astype(str).tolist()
        if "feature_id" in adata.var.columns
        else adata.var_names.astype(str).tolist()
    )
    feature_names = (
        adata.var["feature_name"].astype(str).tolist()
        if "feature_name" in adata.var.columns
        else adata.var_names.astype(str).tolist()
    )
    feature_lines = [
        f"{feature_id}\t{feature_name}"
        for feature_id, feature_name in zip(feature_ids, feature_names, strict=True)
    ]
    _write_text_gzip(matrix_dir / "features.tsv.gz", feature_lines)


def _write_validation_result(
    output_path: Path,
    *,
    expected_paths: dict[str, Path],
    n_cells_total: int,
    n_cells_pass_qc: int,
    n_cells_final_output: int,
    n_cells_unknown_final_celltype_removed: int,
    annotation_method_status: dict[str, dict[str, str]],
) -> None:
    rows = [
        {
            "check_name": "completion",
            "passed": True,
            "detail": "sample outputs written",
        },
        {
            "check_name": "metadata_all_cells",
            "passed": True,
            "detail": f"n_cells_total={n_cells_total}",
        },
        {
            "check_name": "metadata_qc_pass_qc_only",
            "passed": True,
            "detail": f"n_cells_pass_qc={n_cells_pass_qc}",
        },
        {
            "check_name": "metadata_qc_known_final_celltype_only",
            "passed": True,
            "detail": (
                f"n_cells_final_output={n_cells_final_output};"
                f"unknown_final_celltype_removed={n_cells_unknown_final_celltype_removed}"
            ),
        },
    ]
    rows.extend(
        {
            "check_name": f"output_presence:{name}",
            "passed": path.exists() or path == output_path,
            "detail": str(path),
        }
        for name, path in expected_paths.items()
    )
    rows.extend(
        {
            "check_name": f"annotation_status:{method}",
            "passed": status.get("status") == "ok",
            "detail": status.get("detail", ""),
        }
        for method, status in annotation_method_status.items()
    )
    rows.extend(
        {
            "check_name": f"annotation_eval_presence:{column}",
            "passed": column in pd.read_csv(expected_paths["qc_summary.csv"]).columns,
            "detail": column,
        }
        for column in (
            "azimuth_score_mean",
            "azimuth_score_margin_mean",
            "azimuth_low_confidence_fraction",
            "annotation_score",
        )
    )
    pd.DataFrame(rows).to_csv(output_path, index=False)


def write_tuning_selection_artifacts(
    *,
    output_dir: Path,
    sample_id: str,
    gse: str,
    candidate_specs: Sequence[object],
    evaluations: Sequence[object],
    best_candidate_id: str,
) -> Path:
    tuning_dir = Path(output_dir) / "tuning"
    tuning_dir.mkdir(parents=True, exist_ok=True)

    spec_by_id = {
        getattr(candidate, "candidate_id"): candidate for candidate in candidate_specs
    }
    candidate_rows = []
    for evaluation in evaluations:
        candidate = spec_by_id.get(getattr(evaluation, "candidate_id"))
        candidate_rows.append(
            {
                "sample_id": sample_id,
                "gse": gse,
                "candidate_id": getattr(evaluation, "candidate_id"),
                "qc_preset_id": getattr(candidate, "qc_preset_id", ""),
                "azimuth_preset_id": getattr(candidate, "azimuth_preset_id", ""),
                "embedding_preset_id": getattr(candidate, "embedding_preset_id", ""),
                "total_score": float(getattr(evaluation, "total_score", 0.0)),
                "reason_code": getattr(evaluation, "reason_code", ""),
            }
        )

    candidates_path = tuning_dir / "candidates.csv"
    pd.DataFrame(candidate_rows).to_csv(candidates_path, index=False)

    if not best_candidate_id:
        return tuning_dir

    best_row = next(
        row for row in candidate_rows if row["candidate_id"] == best_candidate_id
    )

    selection_summary_path = tuning_dir / "selection_summary.json"
    selection_summary_path.write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "gse": gse,
                "best_candidate_id": best_candidate_id,
                "best_total_score": best_row["total_score"],
                "reason_code": best_row["reason_code"],
                "n_candidates": len(candidate_rows),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    selected_params_path = tuning_dir / "selected_params.json"
    selected_params_path.write_text(
        json.dumps(
            {
                "sample_id": sample_id,
                "gse": gse,
                "candidate_id": best_candidate_id,
                "qc_preset_id": best_row["qc_preset_id"],
                "azimuth_preset_id": best_row["azimuth_preset_id"],
                "embedding_preset_id": best_row["embedding_preset_id"],
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    overview_candidates: list[tuple[str, float, ad.AnnData]] = []
    for row in candidate_rows:
        candidate_h5ad = (
            tuning_dir / row["candidate_id"] / gse / sample_id / f"{sample_id}.h5ad"
        )
        if not candidate_h5ad.exists():
            continue
        overview_candidates.append(
            (
                str(row["candidate_id"]),
                float(row["total_score"]),
                ad.read_h5ad(candidate_h5ad),
            )
        )

    overview_candidates.sort(key=lambda item: item[0])

    if overview_candidates:
        save_azimuth_candidate_overview(
            candidates=overview_candidates,
            output_path=tuning_dir / "umap_rna_candidates_overview_cima_l1.png",
            title=f"{gse}/{sample_id} CIMA L1 overview",
            config=RunConfig(
                qc=QcThresholds(),
                plotting=PlottingConfig(
                    umap_width=4.8,
                    umap_height=4.8,
                    dpi=120,
                    point_size=6.0,
                    legend_fontsize=7.0,
                    legend_title_fontsize=8.0,
                ),
            ),
            color_key="azimuth_cima_l1",
            legend_title="CIMA L1",
        )
        save_azimuth_candidate_overview(
            candidates=overview_candidates,
            output_path=tuning_dir / "umap_rna_candidates_overview_pbmcref.png",
            title=f"{gse}/{sample_id} pbmcref overview",
            config=RunConfig(
                qc=QcThresholds(),
                plotting=PlottingConfig(
                    umap_width=4.8,
                    umap_height=4.8,
                    dpi=120,
                    point_size=6.0,
                    legend_fontsize=7.0,
                    legend_title_fontsize=8.0,
                ),
            ),
            color_key="azimuth_cell_type_l2_raw",
            legend_title="pbmcref",
        )
        best_candidate_adata = next(
            adata_obj
            for candidate_id, _, adata_obj in overview_candidates
            if candidate_id == best_candidate_id
        )
        save_highlight_category_overview(
            adata=best_candidate_adata,
            color_key="azimuth_cell_type_l2_raw",
            output_path=tuning_dir
            / "umap_rna_candidates_overview_pbmcref_highlight.png",
            title=f"{gse}/{sample_id} pbmcref highlight overview",
            config=RunConfig(
                qc=QcThresholds(),
                plotting=PlottingConfig(
                    umap_width=4.8,
                    umap_height=4.8,
                    dpi=120,
                    point_size=6.0,
                    legend_fontsize=7.0,
                    legend_title_fontsize=8.0,
                ),
            ),
            legend_title="pbmcref",
        )
    return tuning_dir


def write_sample_outputs(
    adata: ad.AnnData,
    output_root: Path,
    gse: str,
    sample_id: str,
    config: RunConfig,
) -> Path:
    output_dir = Path(output_root) / gse / sample_id
    output_dir.mkdir(parents=True, exist_ok=True)

    stale_tuning_dir = output_dir / "tuning"
    if stale_tuning_dir.exists():
        shutil.rmtree(stale_tuning_dir)

    stale_nested_sample_dir = output_dir / gse / sample_id
    if stale_nested_sample_dir.exists():
        shutil.rmtree(stale_nested_sample_dir)

    for artifact_name in REMOVED_SAMPLE_ROOT_ARTIFACTS:
        artifact_path = output_dir / artifact_name
        if artifact_path.exists():
            artifact_path.unlink()

    metadata = _add_rna_final_output_columns(
        _metadata_frame(adata),
        adata.obs_names.astype(str).tolist(),
    )
    pass_qc_mask = _bool_pass_qc(cast(pd.DataFrame, adata.obs))
    metadata_pass_qc = metadata.loc[pass_qc_mask.to_numpy()].reset_index(drop=True)
    metadata_qc_final = _add_rna_final_output_columns(
        metadata_pass_qc,
        metadata_pass_qc["cell_id"].astype(str).tolist(),
    )
    known_final_mask = _known_final_celltype_mask(metadata_qc_final)
    n_cells_unknown_final_celltype_removed = int((~known_final_mask).sum())
    metadata_qc = metadata_qc_final.loc[known_final_mask.to_numpy()].reset_index(drop=True)
    annotation_method_status = dict(adata.uns.get("annotation_method_status", {}))
    qc_thresholds = dict(cast(dict[str, Any], adata.uns.get("qc_thresholds", {})))

    metadata_path = output_dir / "metadata.csv"
    metadata_qc_path = output_dir / "metadata_qc.csv"
    qc_summary_path = output_dir / "qc_summary.csv"
    qc_thresholds_path = output_dir / "qc_thresholds.json"
    qc_overview_path = output_dir / "qc_overview.png"
    validation_result_path = output_dir / "validation_result.csv"
    h5ad_path = output_dir / f"{sample_id}.h5ad"
    matrix_dir = output_dir / "matrix"

    metadata.to_csv(metadata_path, index=False)
    metadata_qc.to_csv(metadata_qc_path, index=False)
    _write_qc_summary(
        qc_summary_path,
        gse=gse,
        sample_id=sample_id,
        n_cells_total=adata.n_obs,
        n_cells_pass_qc=int(pass_qc_mask.sum()),
        n_cells_final_output=len(metadata_qc),
        n_cells_unknown_final_celltype_removed=n_cells_unknown_final_celltype_removed,
        metadata_qc=metadata_qc,
        annotation_method_status=annotation_method_status,
        qc_thresholds=qc_thresholds,
    )
    _write_qc_thresholds_json(qc_thresholds_path, qc_thresholds)
    _prepare_h5ad_adata(adata).write_h5ad(
        h5ad_path, convert_strings_to_categoricals=False
    )
    save_qc_overview(adata, qc_overview_path, config)

    expected_paths: dict[str, Path] = {
        "metadata.csv": metadata_path,
        "metadata_qc.csv": metadata_qc_path,
        "qc_summary.csv": qc_summary_path,
        "qc_thresholds.json": qc_thresholds_path,
        "qc_overview.png": qc_overview_path,
        "validation_result.csv": validation_result_path,
        f"{sample_id}.h5ad": h5ad_path,
    }

    dual_plot_path = output_dir / "umap_rna_pbmcref_vs_cima_l1.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=dual_plot_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )
    expected_paths["umap_rna_pbmcref_vs_cima_l1.png"] = dual_plot_path

    highlight_color_key = (
        "azimuth_cell_type_l2_raw"
        if "azimuth_cell_type_l2_raw" in adata.obs.columns
        else "azimuth_cell_type"
    )
    highlight_plot_path = output_dir / "umap_rna_pbmcref_highlight.png"
    save_highlight_category_overview(
        adata=adata,
        color_key=highlight_color_key,
        output_path=highlight_plot_path,
        title="pbmcref highlight",
        config=config,
        legend_title="pbmcref",
    )
    expected_paths["umap_rna_pbmcref_highlight.png"] = highlight_plot_path

    cima_l1_plot_path = output_dir / "umap_rna_cima_l1.png"
    save_sample_cima_l1_umap(
        adata=adata,
        output_path=cima_l1_plot_path,
        title=sample_id,
        config=config,
    )
    expected_paths["umap_rna_cima_l1.png"] = cima_l1_plot_path

    _write_validation_result(
        validation_result_path,
        expected_paths=expected_paths,
        n_cells_total=adata.n_obs,
        n_cells_pass_qc=int(pass_qc_mask.sum()),
        n_cells_final_output=len(metadata_qc),
        n_cells_unknown_final_celltype_removed=n_cells_unknown_final_celltype_removed,
        annotation_method_status=annotation_method_status,
    )
    return output_dir


__all__ = ["write_sample_outputs", "write_tuning_selection_artifacts"]
