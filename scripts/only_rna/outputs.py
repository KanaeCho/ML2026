from __future__ import annotations

import gzip
import json
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .models import RunConfig
from .plotting import save_annotation_method_comparison_umap, save_categorical_umap


UMAP_PLOTS = {
    "cluster": ("umap_rna_clusters.png", "RNA clusters"),
    "cima_l1": ("umap_rna_cima_cell_type_l1.png", "CIMA L1"),
    "cima_l2": ("umap_rna_cima_cell_type_l2.png", "CIMA L2"),
    "cima_l1_masked": ("umap_rna_cima_cell_type_l1_masked.png", "CIMA L1 masked"),
    # Optional alternative annotation overlays (if produced by updated annotation flow)
    "azimuth": ("umap_rna_azimuth.png", "Azimuth"),
    "celltypist": ("umap_rna_celltypist.png", "CellTypist"),
    "singler": ("umap_rna_singler.png", "SingleR"),
    "scanvi": ("umap_rna_scanvi.png", "scANVI"),
}

ANNOTATION_METHOD_COMPARISON_PLOT = (
    "umap_rna_annotation_method_compare.png",
    "Annotation method comparison",
)


def _sanitize_dataframe_for_h5ad(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.index = out.index.astype(str)
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
            )
    return out


def _prepare_h5ad_adata(adata: ad.AnnData) -> ad.AnnData:
    out = adata.copy()
    out.obs = _sanitize_dataframe_for_h5ad(out.obs)
    out.var = _sanitize_dataframe_for_h5ad(out.var)
    out.obs_names = pd.Index(out.obs_names.to_numpy(dtype=str), dtype=object)
    out.var_names = pd.Index(out.var_names.to_numpy(dtype=str), dtype=object)
    return out


def _bool_pass_qc(obs: pd.DataFrame) -> pd.Series:
    if "pass_qc" not in obs.columns:
        raise KeyError("adata.obs must contain 'pass_qc'")
    return obs["pass_qc"].fillna(False).astype(bool)


def _metadata_frame(adata: ad.AnnData) -> pd.DataFrame:
    metadata = adata.obs.copy()
    metadata.insert(0, "cell_id", adata.obs_names.astype(str))
    return metadata.reset_index(drop=True)


def _write_qc_summary(
    output_path: Path,
    *,
    gse: str,
    sample_id: str,
    n_cells_total: int,
    n_cells_pass_qc: int,
    cima_l1_low_confidence_fraction: float,
    annotation_method_status: dict[str, dict[str, str]],
) -> None:
    azimuth_status = annotation_method_status.get("azimuth", {})
    pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "gse": gse,
                "n_cells_total": n_cells_total,
                "n_cells_pass_qc": n_cells_pass_qc,
                "n_cells_fail_qc": n_cells_total - n_cells_pass_qc,
                "pass_qc_fraction": (
                    float(n_cells_pass_qc) / float(n_cells_total)
                    if n_cells_total > 0
                    else 0.0
                ),
                "cima_l1_low_confidence_fraction": float(
                    cima_l1_low_confidence_fraction
                ),
                "azimuth_status": azimuth_status.get("status", ""),
                "azimuth_detail": azimuth_status.get("detail", ""),
            }
        ]
    ).to_csv(output_path, index=False)


def _write_text_gzip(path: Path, lines: list[str]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def _export_matrix_triplet(adata: ad.AnnData, matrix_dir: Path) -> None:
    matrix_dir.mkdir(parents=True, exist_ok=True)
    matrix = adata.X.T
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    else:
        matrix = matrix.tocoo()

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
    pd.DataFrame(rows).to_csv(output_path, index=False)


def write_tuning_selection_artifacts(
    *,
    output_dir: Path,
    sample_id: str,
    gse: str,
    candidate_specs: list[object],
    evaluations: list[object],
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

    metadata = _metadata_frame(adata)
    pass_qc_mask = _bool_pass_qc(adata.obs)
    metadata_qc = metadata.loc[pass_qc_mask.to_numpy()].reset_index(drop=True)
    annotation_method_status = dict(adata.uns.get("annotation_method_status", {}))

    if "cima_l1_low_confidence" in adata.obs.columns:
        low_conf_mask = (
            adata.obs.loc[pass_qc_mask, "cima_l1_low_confidence"]
            .fillna(False)
            .astype(bool)
        )
        cima_l1_low_confidence_fraction = (
            float(low_conf_mask.mean()) if len(low_conf_mask) > 0 else 0.0
        )
    else:
        cima_l1_low_confidence_fraction = 0.0

    metadata_path = output_dir / "metadata.csv"
    metadata_qc_path = output_dir / "metadata_qc.csv"
    qc_summary_path = output_dir / "qc_summary.csv"
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
        cima_l1_low_confidence_fraction=cima_l1_low_confidence_fraction,
        annotation_method_status=annotation_method_status,
    )
    _prepare_h5ad_adata(adata).write_h5ad(
        h5ad_path, convert_strings_to_categoricals=False
    )
    _export_matrix_triplet(adata, matrix_dir)

    expected_paths: dict[str, Path] = {
        "metadata.csv": metadata_path,
        "metadata_qc.csv": metadata_qc_path,
        "qc_summary.csv": qc_summary_path,
        "validation_result.csv": validation_result_path,
        f"{sample_id}.h5ad": h5ad_path,
        "matrix/matrix.mtx": matrix_dir / "matrix.mtx",
        "matrix/barcodes.tsv.gz": matrix_dir / "barcodes.tsv.gz",
        "matrix/features.tsv.gz": matrix_dir / "features.tsv.gz",
    }

    for column, (filename, title) in UMAP_PLOTS.items():
        plot_path = output_dir / filename
        if column in adata.obs.columns:
            save_categorical_umap(
                adata,
                color_key=column,
                output_path=plot_path,
                title=title,
                config=config,
            )
        expected_paths[filename] = plot_path

    comparison_filename, comparison_title = ANNOTATION_METHOD_COMPARISON_PLOT
    comparison_path = output_dir / comparison_filename
    comparison_columns = [
        "azimuth_cell_type",
        "celltypist_cell_type",
        "singler_cell_type",
        "scanvi_cell_type",
    ]
    if all(column in adata.obs.columns for column in comparison_columns):
        save_annotation_method_comparison_umap(
            adata,
            output_path=comparison_path,
            title=comparison_title,
            config=config,
        )
    expected_paths[comparison_filename] = comparison_path

    _write_validation_result(
        validation_result_path,
        expected_paths=expected_paths,
        n_cells_total=adata.n_obs,
        n_cells_pass_qc=int(pass_qc_mask.sum()),
        annotation_method_status=annotation_method_status,
    )
    return output_dir


__all__ = ["write_sample_outputs", "write_tuning_selection_artifacts"]
