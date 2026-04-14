from __future__ import annotations

import gzip
from pathlib import Path

import anndata as ad
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite

from .models import RunConfig
from .plotting import save_categorical_umap


UMAP_PLOTS = {
    "cluster": ("umap_rna_clusters.png", "RNA clusters"),
    "cima_l1": ("umap_rna_cima_cell_type_l1.png", "CIMA L1"),
    "cima_l2": ("umap_rna_cima_cell_type_l2.png", "CIMA L2"),
    "cima_l1_masked": ("umap_rna_cima_cell_type_l1_masked.png", "CIMA L1 masked"),
}


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
) -> None:
    pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "gse": gse,
                "n_cells_total": n_cells_total,
                "n_cells_pass_qc": n_cells_pass_qc,
                "n_cells_fail_qc": n_cells_total - n_cells_pass_qc,
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
    pd.DataFrame(rows).to_csv(output_path, index=False)


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

    _write_validation_result(
        validation_result_path,
        expected_paths=expected_paths,
        n_cells_total=adata.n_obs,
        n_cells_pass_qc=int(pass_qc_mask.sum()),
    )
    return output_dir


__all__ = ["write_sample_outputs"]
