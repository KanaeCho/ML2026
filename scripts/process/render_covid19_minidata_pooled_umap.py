#!/usr/bin/env python3
"""Build pooled exploratory UMAPs for the COVID19-minidata ATAC outputs."""

from __future__ import annotations

import argparse
import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import anndata as ad
import harmonypy as hm
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from sklearn.decomposition import TruncatedSVD
from sklearn.preprocessing import normalize


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_ROOT = ROOT / "output" / "atac" / "COVID19-minidata"
DEFAULT_OUTPUT_DIR = DEFAULT_INPUT_ROOT / "integrated"
DEFAULT_METADATA = ROOT / "data" / "reference" / "COVID19-minidata" / "sample_metadata.csv"

BASE_L1_PALETTE = {
    "B": "#2C7BB6",
    "CD4_T": "#D7191C",
    "CD8_T&unconvensional_T": "#FDAE61",
    "Myeloid": "#1A9641",
    "ILC": "#762A83",
    "Unknown": "#8C8C8C",
}
GROUP_PALETTE = {
    "healthy_control": "#2C7BB6",
    "COVID-19_recovered": "#D95F02",
    "healthy": "#2C7BB6",
    "COVID-19_recovered": "#D95F02",
}


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--sample-metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--n-components", type=int, default=30)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--min-dist", type=float, default=0.35)
    parser.add_argument("--leiden-resolution", type=float, default=0.5)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--point-size", type=float, default=12.0)
    return parser.parse_args()


def sanitize_dataframe(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out.index = pd.Index(["" if pd.isna(v) else str(v) for v in out.index], dtype=object)
    out.columns = pd.Index([str(v) for v in out.columns], dtype=object)
    for column in out.columns:
        series = out[column]
        if (
            pd.api.types.is_object_dtype(series.dtype)
            or pd.api.types.is_string_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            out[column] = pd.Series(
                ["" if pd.isna(v) else str(v) for v in series.tolist()],
                index=series.index,
                dtype=object,
            )
    return out


def discover_h5ad(input_root: Path) -> list[Path]:
    paths = []
    for path in sorted(input_root.glob("*/*.h5ad")):
        if path.parent.name == "integrated":
            continue
        paths.append(path)
    if not paths:
        raise FileNotFoundError(f"No sample h5ad files found under {input_root}")
    return paths


def load_sample_metadata(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    frame = pd.read_csv(path, dtype=str).fillna("")
    if "sample_id" not in frame.columns:
        return {}
    return {str(row["sample_id"]): {str(k): str(v) for k, v in row.items()} for _, row in frame.iterrows()}


def load_pooled(paths: list[Path], sample_metadata: dict[str, dict[str, str]]) -> ad.AnnData:
    adatas: list[ad.AnnData] = []
    var_names: pd.Index | None = None
    for path in paths:
        sample_id = path.parent.name
        sample = ad.read_h5ad(path)
        if var_names is None:
            var_names = pd.Index(sample.var_names.astype(str))
        elif not var_names.equals(pd.Index(sample.var_names.astype(str))):
            raise ValueError(f"Feature order mismatch in {path}")
        sample.obs = sanitize_dataframe(sample.obs)
        sample.var = sanitize_dataframe(sample.var)
        sample.obs["sample_id"] = sample_id
        sample.obs["source_h5ad"] = str(path)
        for key, value in sample_metadata.get(sample_id, {}).items():
            if key == "sample_id":
                continue
            sample.obs[key] = value
        sample.obs_names = pd.Index([f"{sample_id}:{barcode}" for barcode in sample.obs_names.astype(str)], dtype=object)
        adatas.append(sample)

    pooled = ad.concat(adatas, join="inner", merge="first", index_unique=None)
    pooled.obs = sanitize_dataframe(pooled.obs)
    pooled.var = sanitize_dataframe(pooled.var)
    pooled.uns["source_h5ad_files"] = [str(path) for path in paths]
    pooled.uns["dataset"] = "COVID19-minidata"
    pooled.uns["created_at"] = utc_now()
    return pooled


def binary_tfidf_lsi(matrix: sp.spmatrix, n_components: int, random_state: int) -> np.ndarray:
    binary = sp.csr_matrix(matrix, dtype=np.float32, copy=True)
    if binary.nnz:
        binary.data[:] = 1.0
    tf = normalize(binary, norm="l1", axis=1, copy=True)
    peak_counts = np.asarray(binary.sum(axis=0)).ravel().astype(np.float32)
    idf = np.log1p(binary.shape[0] / np.maximum(peak_counts, 1.0)).astype(np.float32)
    tfidf = tf @ sp.diags(idf, dtype=np.float32)
    n_components = min(n_components, max(2, min(tfidf.shape) - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=random_state)
    lsi = svd.fit_transform(tfidf).astype(np.float32)
    if lsi.shape[1] > 1:
        lsi = lsi[:, 1:]
    return lsi


def run_raw_umap(adata: ad.AnnData, n_neighbors: int, min_dist: float, resolution: float, random_state: int) -> None:
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_lsi", metric="cosine", random_state=random_state)
    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="raw_cluster",
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    adata.obs["raw_umap_1"] = adata.obsm["X_umap"][:, 0]
    adata.obs["raw_umap_2"] = adata.obsm["X_umap"][:, 1]


def run_harmony_umap(adata: ad.AnnData, n_neighbors: int, min_dist: float, resolution: float, random_state: int) -> dict[str, Any]:
    harmony = hm.run_harmony(adata.obsm["X_lsi"], adata.obs, ["sample_id"], random_state=random_state)
    corrected = np.asarray(harmony.Z_corr, dtype=np.float32)
    if corrected.shape[0] != adata.n_obs and corrected.shape[1] == adata.n_obs:
        corrected = corrected.T
    adata.obsm["X_harmony"] = corrected
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X_harmony", metric="cosine", random_state=random_state)
    sc.tl.umap(adata, min_dist=min_dist, random_state=random_state)
    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added="integrated_cluster",
        random_state=random_state,
        flavor="igraph",
        n_iterations=2,
        directed=False,
    )
    adata.obs["integrated_umap_1"] = adata.obsm["X_umap"][:, 0]
    adata.obs["integrated_umap_2"] = adata.obsm["X_umap"][:, 1]
    return {
        "method": "harmony",
        "batch_key": "sample_id",
        "n_iterations": int(getattr(harmony, "converged", False) is False) if hasattr(harmony, "converged") else None,
    }


def palette_for(values: list[str], column: str) -> dict[str, str]:
    labels = [label for label in dict.fromkeys(str(v) for v in values if str(v) and str(v) != "nan")]
    if column in {"final_celltype", "cima_cell_type_l1", "cima_cell_type_l1_masked"}:
        base = {label: BASE_L1_PALETTE[label] for label in labels if label in BASE_L1_PALETTE}
    elif column == "group":
        base = {label: GROUP_PALETTE[label] for label in labels if label in GROUP_PALETTE}
    else:
        base = {}
    missing = [label for label in labels if label not in base]
    cmap = plt.get_cmap("tab20", max(len(missing), 1))
    for i, label in enumerate(missing):
        base[label] = matplotlib.colors.to_hex(cmap(i))
    return base


def plot_umap(
    metadata: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    out_path: Path,
    title: str,
    point_size: float,
) -> None:
    plot_df = metadata[[x_col, y_col, color_col]].copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df[color_col] = plot_df[color_col].fillna("Unknown").astype(str).replace({"": "Unknown", "nan": "Unknown"})
    plot_df = plot_df.dropna(subset=[x_col, y_col])
    palette = palette_for(plot_df[color_col].tolist(), color_col)

    fig, ax = plt.subplots(figsize=(8.5, 7.5), constrained_layout=True)
    for label, group in plot_df.groupby(color_col, sort=True):
        ax.scatter(
            group[x_col],
            group[y_col],
            s=point_size,
            alpha=0.82,
            linewidths=0,
            c=palette.get(str(label), "#666666"),
            label=f"{label} (n={len(group):,})",
            rasterized=True,
        )
    ax.set_title(f"{title}\nn={len(plot_df):,} QC cells", fontsize=14, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=1.4, fontsize=9)
    fig.savefig(out_path, dpi=220)
    plt.close(fig)


def write_summaries(adata: ad.AnnData, output_dir: Path, config: dict[str, Any], harmony_detail: dict[str, Any]) -> None:
    metadata = adata.obs.copy()
    metadata.to_csv(output_dir / "metadata.csv", index=True, index_label="cell_id")

    sample_qc = (
        metadata.groupby(["sample_id", "group", "health", "donor", "sex"], dropna=False)
        .size()
        .reset_index(name="n_qc_cells")
        .sort_values("sample_id")
    )
    sample_qc.to_csv(output_dir / "sample_qc_summary.csv", index=False)

    celltype_by_group = (
        metadata.groupby(["group", "final_celltype"], dropna=False)
        .size()
        .reset_index(name="n_cells")
    )
    celltype_by_group["fraction_within_group"] = celltype_by_group["n_cells"] / celltype_by_group.groupby("group")["n_cells"].transform("sum")
    celltype_by_group.to_csv(output_dir / "celltype_by_group.csv", index=False)

    summary = {
        "dataset": "COVID19-minidata",
        "created_at": utc_now(),
        "n_cells": int(adata.n_obs),
        "n_peaks": int(adata.n_vars),
        "n_samples": int(metadata["sample_id"].nunique()),
        "raw_umap": {
            "basis": "X_lsi",
            "cluster_column": "raw_cluster",
        },
        "integrated_umap": {
            "basis": "X_harmony",
            "cluster_column": "integrated_cluster",
            **harmony_detail,
        },
        "config": config,
        "interpretation_warning": "Exploratory pooled visualization only; COVID19-minidata has low FRiP/TSS and high CIMA low-confidence fraction.",
    }
    (output_dir / "integration_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = discover_h5ad(args.input_root)
    sample_metadata = load_sample_metadata(args.sample_metadata)
    pooled = load_pooled(paths, sample_metadata)

    pooled.obsm["X_lsi"] = binary_tfidf_lsi(pooled.X, args.n_components, args.random_state)
    run_raw_umap(pooled, args.n_neighbors, args.min_dist, args.leiden_resolution, args.random_state)
    harmony_detail = run_harmony_umap(pooled, args.n_neighbors, args.min_dist, args.leiden_resolution, args.random_state)

    config = {
        "n_components_requested": args.n_components,
        "n_components_used_after_dropping_first_lsi": int(pooled.obsm["X_lsi"].shape[1]),
        "n_neighbors": args.n_neighbors,
        "min_dist": args.min_dist,
        "leiden_resolution": args.leiden_resolution,
        "random_state": args.random_state,
    }
    write_summaries(pooled, args.output_dir, config, harmony_detail)

    for prefix, x_col, y_col in [
        ("raw", "raw_umap_1", "raw_umap_2"),
        ("integrated", "integrated_umap_1", "integrated_umap_2"),
    ]:
        for color_col in ["sample_id", "group", "final_celltype", "cima_cell_type_l1_masked", "sex", "comorbidity"]:
            if color_col not in pooled.obs.columns:
                continue
            plot_umap(
                pooled.obs,
                x_col,
                y_col,
                color_col,
                args.output_dir / f"umap_{prefix}_by_{color_col}.png",
                f"COVID19-minidata {prefix} pooled UMAP by {color_col}",
                args.point_size,
            )

    pooled.obs = sanitize_dataframe(pooled.obs)
    pooled.var = sanitize_dataframe(pooled.var)
    pooled.write_h5ad(args.output_dir / "covid19_minidata_pooled.h5ad", convert_strings_to_categoricals=False)
    print(f"pooled_cells={pooled.n_obs}")
    print(f"pooled_peaks={pooled.n_vars}")
    print(f"output_dir={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
