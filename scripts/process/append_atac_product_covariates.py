from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
import pandas as pd
from scipy import sparse


ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = ROOT / "output"
DATA_ROOT = ROOT / "data"

PRODUCT_DIRS = {
    "only_atac": OUTPUT_ROOT / "1.only_atac",
    "co_atac": OUTPUT_ROOT / "3.co_atac",
}

GSE190992_DONORS = {
    "GSM5737281": "PTID2",
    "GSM5737282": "PTID2",
    "GSM5737283": "PTID2",
    "GSM5737284": "PTID2",
    "GSM5737285": "PTID2",
    "GSM5737286": "PTID2",
    "GSM5737287": "PTID4",
    "GSM5737288": "PTID4",
    "GSM5737289": "PTID4",
    "GSM5737290": "PTID4",
    "GSM5737291": "PTID4",
    "GSM5737292": "PTID4",
    "GSM5737293": "PTID5",
    "GSM5737294": "PTID5",
    "GSM5737295": "PTID5",
    "GSM5737296": "PTID6",
    "GSM5737297": "PTID6",
    "GSM5737298": "PTID6",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Append audit covariates to ATAC product metadata.")
    parser.add_argument(
        "--products",
        nargs="+",
        choices=sorted(PRODUCT_DIRS),
        default=sorted(PRODUCT_DIRS),
        help="ATAC products to update.",
    )
    return parser.parse_args()


def clean_value(value: Any) -> str:
    if pd.isna(value):
        return ""
    text = str(value).strip()
    if text.lower() in {"nan", "none", "<na>"}:
        return ""
    return text


def normalize_health(value: Any) -> str:
    text = clean_value(value)
    lower = text.lower()
    if not text:
        return ""
    if text in {"健康", "Healthy", "healthy control", "Normal"} or lower in {"healthy", "healthy control", "normal"}:
        return "healthy"
    if "covid" in lower or "ncov" in lower:
        return "COVID-19"
    if "rsv" in lower:
        return "RSV"
    if "flu" in lower or "influenza" in lower:
        return "influenza"
    if text == "发病":
        return "symptomatic"
    if text == "暴露":
        return "exposed"
    return text


def load_dataset_covariates(reference_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = reference_dir / "datasets.xlsx"
    df = pd.read_excel(path)
    sample_col = "样本名(GSM*)"
    gse_col = "数据集(GSE*)"
    assay_col = "测序数据(scATAC/scRNA)"
    covariates: dict[tuple[str, str], dict[str, str]] = {}
    atac = df[df[assay_col].astype(str).str.contains("scATAC", na=False)].copy()
    for _, row in atac.iterrows():
        gse = clean_value(row[gse_col])
        sample_id = clean_value(row[sample_col])
        if not gse or not sample_id:
            continue
        covariates[(gse, sample_id)] = {
            "x_dataset": gse,
            "x_sample": sample_id,
            "x_donor": "",
            "x_age": clean_value(row.get("年龄", "")),
            "x_health": normalize_health(row.get("健康状态", "")),
        }
    for (gse, sample_id), values in covariates.items():
        if gse == "GSE190992":
            values["x_donor"] = GSE190992_DONORS.get(sample_id, "")
            values["x_age"] = "25-38"
            values["x_health"] = "healthy"
    return covariates


def load_gse283744_donors_from_titles(reference_dir: Path) -> dict[str, str]:
    path = reference_dir / "gse283744_atac_sample_titles.csv"
    if not path.exists():
        return {}
    df = pd.read_csv(path)
    return {clean_value(row["sample_id"]): clean_value(row["donor"]) for _, row in df.iterrows()}


def load_only_atac_covariates(reference_dir: Path) -> dict[tuple[str, str], dict[str, str]]:
    covariates = load_dataset_covariates(reference_dir)
    gse283744_donors = load_gse283744_donors_from_titles(reference_dir)
    for key, values in covariates.items():
        if key[0] == "GSE283744":
            values["x_donor"] = gse283744_donors.get(key[1], "")
    return covariates


def load_co_atac_covariates(data_root: Path) -> dict[tuple[str, str], dict[str, str]]:
    path = data_root / "raw" / "7555405" / "sample_layout.tsv"
    df = pd.read_csv(path, sep="\t")
    covariates: dict[tuple[str, str], dict[str, str]] = {}
    for _, row in df.iterrows():
        sample_id = clean_value(row["sample_id"])
        covariates[("7555405", sample_id)] = {
            "x_dataset": "7555405",
            "x_sample": sample_id,
            "x_donor": clean_value(row["donor"]),
            "x_age": clean_value(row["age"]),
            "x_health": "healthy",
        }
    return covariates


def append_covariates(df: pd.DataFrame, covariates: dict[tuple[str, str], dict[str, str]]) -> pd.DataFrame:
    out = df.copy()
    keys = list(zip(out["gse"].astype(str), out["sample_id"].astype(str), strict=False))
    for column in ["x_dataset", "x_sample", "x_donor", "x_age", "x_health"]:
        out[column] = [covariates.get(key, {}).get(column, "") for key in keys]
    return out


def find_h5ad(product_dir: Path, product: str) -> Path | None:
    preferred = product_dir / f"{product}.h5ad"
    if preferred.exists():
        return preferred
    candidates = sorted(product_dir.glob("*.h5ad"))
    if len(candidates) == 1:
        return candidates[0]
    return None


def summarize_coverage(product: str, samples: pd.DataFrame, cells: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {"product": product, "n_samples": int(len(samples)), "n_cells": int(len(cells))}
    for column in ["x_dataset", "x_sample", "x_donor", "x_age", "x_health"]:
        sample_mask = samples[column].astype(str).str.len() > 0
        cell_mask = cells[column].astype(str).str.len() > 0
        summary[f"samples_with_{column}"] = int(sample_mask.sum())
        summary[f"cells_with_{column}"] = int(cell_mask.sum())
    unmapped = samples[(samples["x_dataset"].astype(str).str.len() == 0) | (samples["x_age"].astype(str).str.len() == 0)]
    summary["unmapped_samples"] = [
        {"gse": clean_value(row["gse"]), "sample_id": clean_value(row["sample_id"])} for _, row in unmapped.iterrows()
    ]
    return summary


def update_h5ad(product_dir: Path, product: str, cells: pd.DataFrame) -> str:
    h5ad_path = find_h5ad(product_dir, product)
    if h5ad_path is None:
        return "skipped:no_unique_h5ad"
    existing = ad.read_h5ad(h5ad_path)
    if len(existing.obs) != len(cells):
        return f"skipped:obs_length_mismatch:{len(existing.obs)}!={len(cells)}"
    index_col = "global_cell_id" if "global_cell_id" in cells.columns else None
    obs = cells.copy()
    if index_col:
        obs.index = pd.Index([str(value) for value in obs[index_col].tolist()], dtype=object, name=None)
    else:
        obs.index = pd.Index([str(value) for value in range(len(obs))], dtype=object)
    for column in obs.columns:
        if not pd.api.types.is_numeric_dtype(obs[column]) and not pd.api.types.is_bool_dtype(obs[column]):
            obs[column] = pd.Series([clean_value(value) for value in obs[column].tolist()], index=obs.index, dtype=object)
    var = pd.DataFrame(index=pd.Index([], dtype=object))
    adata = ad.AnnData(X=sparse.csr_matrix((len(obs), 0)), obs=obs, var=var)
    if {"integrated_umap_1", "integrated_umap_2"}.issubset(obs.columns):
        adata.obsm["X_integrated_umap"] = obs[["integrated_umap_1", "integrated_umap_2"]].to_numpy(dtype="float32")
    if {"cima_ref_umap_1", "cima_ref_umap_2"}.issubset(obs.columns):
        adata.obsm["X_cima_ref_umap"] = obs[["cima_ref_umap_1", "cima_ref_umap_2"]].to_numpy(dtype="float32")
    if {"umap_atac_1", "umap_atac_2"}.issubset(obs.columns):
        adata.obsm["X_umap_atac"] = obs[["umap_atac_1", "umap_atac_2"]].to_numpy(dtype="float32")
    adata.write_h5ad(h5ad_path, convert_strings_to_categoricals=False)
    return str(h5ad_path)


def update_product(product: str) -> dict[str, Any]:
    product_dir = PRODUCT_DIRS[product]
    samples_path = product_dir / "manifests" / "samples.csv"
    cells_path = product_dir / "manifests" / "cells_metadata.csv"
    samples = pd.read_csv(samples_path, low_memory=False)
    cells = pd.read_csv(cells_path, low_memory=False)
    if product == "only_atac":
        covariates = load_only_atac_covariates(DATA_ROOT / "reference")
    elif product == "co_atac":
        covariates = load_co_atac_covariates(DATA_ROOT)
    else:  # pragma: no cover - protected by argparse choices
        raise ValueError(product)
    samples = append_covariates(samples, covariates)
    cells = append_covariates(cells, covariates)
    samples.to_csv(samples_path, index=False)
    cells.to_csv(cells_path, index=False)
    summary = summarize_coverage(product, samples, cells)
    summary["h5ad"] = update_h5ad(product_dir, product, cells)
    summary_path = product_dir / "manifests" / "covariate_coverage.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    args = parse_args()
    summaries = [update_product(product) for product in args.products]
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
