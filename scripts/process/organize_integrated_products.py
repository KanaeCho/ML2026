#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.process.render_product_umap_panels import render_panels
from scripts.process.integrate_product_embeddings import (
    IntegrationConfig,
    integrate_product_embeddings,
)

PRODUCT_NAMES = ["only_atac", "only_rna", "co_atac", "co_rna"]
PRODUCT_DIRS = {
    "only_atac": "1.only_atac",
    "only_rna": "2.only_rna",
    "co_atac": "3.co_atac",
    "co_rna": "4.co_rna",
}
PRODUCT_BRANCH = {
    "only_atac": "only_atac",
    "only_rna": "only_rna",
    "co_atac": "co",
    "co_rna": "co",
}
PRODUCT_MODALITY = {
    "only_atac": "ATAC",
    "only_rna": "RNA",
    "co_atac": "ATAC",
    "co_rna": "RNA",
}
DEFAULT_INTEGRATION_BATCH_KEY = "sample_id"
DEFAULT_BBKNN_NEIGHBORS_WITHIN_BATCH = 1
DEFAULT_BBKNN_TRIM = 60
DEFAULT_INTEGRATION_METHOD = "bbknn"
RNA_ONLY_INTEGRATION_METHOD = "harmony"
RNA_ONLY_INTEGRATION_BATCH_KEY = "gse"
RNA_ONLY_BBKNN_NEIGHBORS_WITHIN_BATCH = 3
RNA_ONLY_BBKNN_TRIM = 30
RNA_ONLY_MIN_CIMA_L1_SCORE = 0.5
EXCLUDED_OUTPUT_NAMES = {
    "1.only_atac",
    "2.only_rna",
    "3.co_atac",
    "4.co_rna",
    "co",
    "reference",
    "rna",
    "atac",
}
EXCLUDED_ONLY_ATAC_GSES = {
    # GSE206284 is the removed co2 dataset, not part of the accepted legacy
    # only_atac product. Including it creates isolated
    # sample-driven structures in product-level integration.
    "GSE206284",
    # GSE282769 was downloaded and partially processed from atac.xlsx, but the
    # CIMA ATAC labels did not resolve readable cell types for product review.
    "GSE282769",
}
EXCLUDED_ONLY_RNA_GSES = {
    # GSE206284 is the removed co2 dataset and should not enter only_rna
    # product integration if stale outputs are present on disk.
    "GSE206284",
}
COMMON_REQUIRED = [
    "metadata.csv",
    "metadata_qc.csv",
    "qc_summary.csv",
    "validation_result.csv",
    "run_status.json",
]
RNA_REQUIRED = [
    *COMMON_REQUIRED,
    "qc_overview.png",
]
ATAC_REQUIRED = [
    *COMMON_REQUIRED,
    "qc_overview.png",
    "umap_cima_cell_type_l1.png",
    "umap_cima_cell_type_l2.png",
]
LEGACY_ATAC_REQUIRED = [
    "qc_summary.csv",
    "validation_result.csv",
    "run_status.json",
    "umap_cima_cell_type_l1.png",
]


@dataclass(frozen=True)
class SampleProduct:
    product: str
    gse: str
    sample_id: str
    source_dir: Path
    required_files: tuple[str, ...]
    branch: str
    modality: str
    is_co_sample: bool = False
    co_dataset: str = ""
    co_dataset_id: str = ""


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Organize annotated branch outputs into integrated product directories."
    )
    parser.add_argument(
        "--products",
        default="all",
        help="Comma-separated products: only_atac,only_rna,co_atac,co_rna,all",
    )
    parser.add_argument("--output-root", default=str(ROOT / "output"))
    parser.add_argument("--copy-mode", choices=["symlink", "copy"], default="symlink")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--skip-integration", action="store_true")
    parser.add_argument("--integration-n-components", type=int, default=30)
    parser.add_argument("--integration-max-umap-fit-cells", type=int, default=100_000)
    parser.add_argument("--integration-clusters", type=int, default=30)
    parser.add_argument("--integration-batch-key", default=None)
    parser.add_argument(
        "--integration-method",
        choices=["bbknn", "scanpy_neighbors", "harmony"],
        default=None,
    )
    parser.add_argument("--bbknn-neighbors-within-batch", type=int, default=None)
    parser.add_argument("--bbknn-trim", type=int, default=None)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--rna-min-cima-l1-score", type=float, default=None)
    parser.add_argument("--include-incomplete", action="store_true")
    return parser.parse_args()


def selected_products(raw: str) -> list[str]:
    if raw.strip() == "all":
        return list(PRODUCT_NAMES)
    values = [value.strip() for value in raw.split(",") if value.strip()]
    unknown = sorted(set(values) - set(PRODUCT_NAMES))
    if unknown:
        raise ValueError(f"Unknown product(s): {', '.join(unknown)}")
    return values


def default_integration_batch_key(product: str) -> str:
    if product == "only_rna":
        return RNA_ONLY_INTEGRATION_BATCH_KEY
    return DEFAULT_INTEGRATION_BATCH_KEY


def default_integration_method(product: str) -> str:
    if product in {"only_rna", "co_rna"}:
        return RNA_ONLY_INTEGRATION_METHOD
    return DEFAULT_INTEGRATION_METHOD


def default_bbknn_neighbors_within_batch(product: str) -> int:
    if product == "only_rna":
        return RNA_ONLY_BBKNN_NEIGHBORS_WITHIN_BATCH
    return DEFAULT_BBKNN_NEIGHBORS_WITHIN_BATCH


def default_bbknn_trim(product: str) -> int:
    if product == "only_rna":
        return RNA_ONLY_BBKNN_TRIM
    return DEFAULT_BBKNN_TRIM


def default_rna_min_cima_l1_score(product: str) -> float:
    if product == "only_rna":
        return RNA_ONLY_MIN_CIMA_L1_SCORE
    return 0.0


def read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}


def run_status_for(sample_dir: Path) -> tuple[str, bool | None]:
    payload = read_json(sample_dir / "run_status.json")
    status = str(payload.get("status", "")) if payload else ""
    outputs_complete = payload.get("outputs_complete") if payload else None
    if isinstance(outputs_complete, bool):
        return status, outputs_complete
    return status, None


def required_with_h5(sample: SampleProduct) -> list[str]:
    required = list(sample.required_files)
    if sample.modality == "RNA":
        required.append(f"{sample.sample_id}.h5ad")
    if sample.modality == "ATAC":
        required.append(f"{sample.sample_id}.h5ad")
    return required


def missing_required_files(sample: SampleProduct) -> list[str]:
    missing = []
    for relative in required_with_h5(sample):
        path = sample.source_dir / relative
        if not path.exists():
            missing.append(relative)
    if sample.modality == "RNA" and not any(
        (sample.source_dir / name).exists()
        for name in ["umap_rna_cima_l1.png", "umap_rna_pbmcref_vs_cima_l1.png", "umap_rna_pbmcref_highlight.png"]
    ):
        missing.append("rna_umap_png")
    return missing


def discover_only_rna(output_root: Path) -> list[SampleProduct]:
    source_root = output_root / "rna"
    if not source_root.exists():
        return []
    samples = []
    for meta_path in sorted(source_root.glob("*/*/metadata_qc.csv")):
        if "tuning" in meta_path.parts:
            continue
        sample_dir = meta_path.parent
        gse = sample_dir.parent.name
        if gse in EXCLUDED_ONLY_RNA_GSES:
            continue
        sample_id = sample_dir.name
        samples.append(
            SampleProduct(
                product="only_rna",
                gse=gse,
                sample_id=sample_id,
                source_dir=sample_dir,
                required_files=tuple(RNA_REQUIRED),
                branch=PRODUCT_BRANCH["only_rna"],
                modality=PRODUCT_MODALITY["only_rna"],
            )
        )
    return samples


def discover_only_atac(output_root: Path, atac_workbook: Path | None = None) -> list[SampleProduct]:
    if atac_workbook is None:
        atac_workbook = ROOT / "data" / "reference" / "atac.xlsx"
    selected: set[tuple[str, str]] | None = None
    if atac_workbook.exists():
        rules = pd.read_excel(atac_workbook, dtype=str).fillna("")
        if {"dataset", "sample"}.issubset(rules.columns):
            rules["dataset"] = rules["dataset"].astype(str).str.strip()
            rules["sample"] = rules["sample"].astype(str).str.strip()
            selected = set(zip(rules["dataset"], rules["sample"], strict=False))
    samples = []
    source_root = output_root / "atac"
    if not source_root.exists():
        return []
    for gse_dir in sorted(path for path in source_root.iterdir() if path.is_dir()):
        if gse_dir.name in EXCLUDED_OUTPUT_NAMES or not gse_dir.name.startswith("GSE"):
            continue
        if gse_dir.name in EXCLUDED_ONLY_ATAC_GSES:
            continue
        if gse_dir.name.endswith("_tea_seq"):
            continue
        for sample_dir in sorted(path for path in gse_dir.iterdir() if path.is_dir()):
            if selected is not None and (gse_dir.name, sample_dir.name) not in selected:
                continue
            if sample_dir.name == "qc_audit" or not (sample_dir / "validation_result.csv").exists():
                continue
            samples.append(
                SampleProduct(
                    product="only_atac",
                    gse=gse_dir.name,
                    sample_id=sample_dir.name,
                    source_dir=sample_dir,
                    required_files=tuple(LEGACY_ATAC_REQUIRED),
                    branch=PRODUCT_BRANCH["only_atac"],
                    modality=PRODUCT_MODALITY["only_atac"],
                )
            )
    return samples


def discover_co_atac(output_root: Path) -> list[SampleProduct]:
    source_root = output_root / "co" / "atac" / "7555405"
    return discover_co_modality(source_root, "co_atac", tuple(ATAC_REQUIRED))


def discover_co_rna(output_root: Path) -> list[SampleProduct]:
    source_root = output_root / "co" / "rna" / "7555405"
    return discover_co_modality(source_root, "co_rna", tuple(RNA_REQUIRED))


def discover_co_modality(source_root: Path, product: str, required: tuple[str, ...]) -> list[SampleProduct]:
    if not source_root.exists():
        return []
    samples = []
    for sample_dir in sorted(path for path in source_root.iterdir() if path.is_dir()):
        if not (sample_dir / "metadata_qc.csv").exists():
            continue
        samples.append(
            SampleProduct(
                product=product,
                gse="7555405",
                sample_id=sample_dir.name,
                source_dir=sample_dir,
                required_files=required,
                branch=PRODUCT_BRANCH[product],
                modality=PRODUCT_MODALITY[product],
                is_co_sample=True,
                co_dataset="co1",
                co_dataset_id="7555405",
            )
        )
    return samples


def discover_product(
    output_root: Path, product: str, atac_workbook: Path | None = None
) -> list[SampleProduct]:
    if product == "only_rna":
        return discover_only_rna(output_root)
    if product == "only_atac":
        return discover_only_atac(output_root, atac_workbook=atac_workbook)
    if product == "co_atac":
        return discover_co_atac(output_root)
    if product == "co_rna":
        return discover_co_rna(output_root)
    raise ValueError(f"Unknown product: {product}")


def ensure_clean_product_dir(product_dir: Path, force: bool) -> None:
    if product_dir.exists():
        if not force:
            return
        shutil.rmtree(product_dir)
    product_dir.mkdir(parents=True, exist_ok=True)


def materialize_sample(sample: SampleProduct, product_dir: Path, copy_mode: str, force: bool) -> Path:
    dest = product_dir / "samples" / sample.gse / sample.sample_id
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if force:
            if dest.is_symlink() or dest.is_file():
                dest.unlink()
            else:
                shutil.rmtree(dest)
        else:
            return dest
    if copy_mode == "symlink":
        try:
            dest.symlink_to(sample.source_dir, target_is_directory=True)
        except OSError:
            dest.mkdir(parents=True, exist_ok=True)
            (dest / "SOURCE_OUTPUT_DIR.txt").write_text(
                str(sample.source_dir), encoding="utf-8"
            )
    else:
        shutil.copytree(sample.source_dir, dest)
    return dest


def sample_row(sample: SampleProduct, product_sample_dir: Path) -> dict[str, object]:
    missing = missing_required_files(sample)
    status, outputs_complete = run_status_for(sample.source_dir)
    sample_complete = not missing and (status in {"", "success"})
    return {
        "product": sample.product,
        "branch": sample.branch,
        "modality": sample.modality,
        "gse": sample.gse,
        "sample_id": sample.sample_id,
        "individual_id": "",
        "source_output_dir": str(sample.source_dir),
        "product_sample_dir": str(product_sample_dir),
        "sample_complete": sample_complete,
        "missing_files": ";".join(missing),
        "run_status": status,
        "outputs_complete": outputs_complete if outputs_complete is not None else "",
        "is_co_sample": sample.is_co_sample,
        "co_dataset": sample.co_dataset,
        "co_dataset_id": sample.co_dataset_id,
    }


def add_traceability(meta: pd.DataFrame, row: dict[str, object]) -> pd.DataFrame:
    meta = meta.copy()
    if "cell_id" not in meta.columns:
        if "cell_barcode" in meta.columns:
            meta.insert(0, "cell_id", meta["cell_barcode"].astype(str))
        else:
            meta.insert(0, "cell_id", [f"cell_{i}" for i in range(len(meta))])
    source_cell_id = pd.Series(meta["cell_id"], index=meta.index).astype(str)
    for column in [
        "product",
        "branch",
        "modality",
        "gse",
        "sample_id",
        "is_co_sample",
        "co_dataset",
        "co_dataset_id",
        "source_output_dir",
    ]:
        meta[column] = row.get(column, "")
    has_individual = bool(meta["individual_id"].notna().any()) if "individual_id" in meta.columns else False
    if not has_individual:
        meta["individual_id"] = row["sample_id"] if row.get("is_co_sample") else ""
    meta["source_cell_id"] = source_cell_id
    meta["global_cell_id"] = (
        meta["product"].astype(str)
        + "__"
        + meta["gse"].astype(str)
        + "__"
        + meta["sample_id"].astype(str)
        + "__"
        + source_cell_id
    )
    return meta


def read_sample_metadata(sample: SampleProduct, row: dict[str, object]) -> pd.DataFrame:
    h5ad_path = sample.source_dir / f"{sample.sample_id}.h5ad"
    if h5ad_path.exists():
        h5ad = ad.read_h5ad(h5ad_path)
        meta = pd.DataFrame(h5ad.obs.copy()).reset_index(drop=True)
    else:
        metadata_path = sample.source_dir / "metadata_qc.csv"
        if not metadata_path.exists() and sample.modality == "ATAC":
            metadata_path = sample.source_dir / "validation_result.csv"
        meta = pd.read_csv(metadata_path, low_memory=False)
    return add_traceability(meta, row)


def write_metadata_object(cells: pd.DataFrame, product_dir: Path, product: str) -> None:
    if cells.empty:
        return
    obs = cells.copy()
    obs.index = pd.Index(pd.Series(obs["global_cell_id"], index=obs.index).astype(str).tolist(), dtype=object)
    for column in obs.columns:
        if pd.api.types.is_bool_dtype(obs[column]):
            obs[column] = obs[column].astype(bool)
        elif pd.api.types.is_numeric_dtype(obs[column]):
            obs[column] = pd.to_numeric(obs[column], errors="coerce")
        else:
            obs[column] = obs[column].map(lambda value: "" if pd.isna(value) else str(value)).astype(object)
    var = pd.DataFrame(index=pd.Index([], dtype=object))
    adata = ad.AnnData(X=sparse.csr_matrix((len(obs), 0), dtype=np.float32), obs=obs, var=var)
    for x_col, y_col, key in [
        ("integrated_umap_1", "integrated_umap_2", "X_integrated_umap"),
        ("cima_ref_umap_1", "cima_ref_umap_2", "X_cima_ref_umap"),
        ("umap_atac_1", "umap_atac_2", "X_umap_atac"),
        ("umap_1", "umap_2", "X_umap"),
    ]:
        if x_col in obs.columns and y_col in obs.columns:
            coords = obs[[x_col, y_col]].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=np.float32)
            if np.isfinite(coords).any():
                adata.obsm[key] = coords
    adata.uns["product"] = product
    adata.uns["created_at"] = utc_now()
    adata.write_h5ad(product_dir / f"{product}.h5ad", convert_strings_to_categoricals=False)


def write_output_inventory(product_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted((product_dir / "samples").glob("**/*")):
        if path.is_dir() and not path.is_symlink():
            continue
        rows.append(
            {
                "relative_path": str(path.relative_to(product_dir)),
                "path": str(path),
                "is_symlink": path.is_symlink(),
                "target": str(path.resolve()) if path.is_symlink() else "",
            }
        )
    inventory = pd.DataFrame(rows)
    inventory.to_csv(product_dir / "manifests" / "output_files.csv", index=False)
    return inventory


def concat_existing_csv(samples: list[SampleProduct], filename: str, rows: list[dict[str, object]]) -> pd.DataFrame:
    parts = []
    row_by_source = {str(row["source_output_dir"]): row for row in rows}
    for sample in samples:
        path = sample.source_dir / filename
        if not path.exists():
            continue
        df = pd.read_csv(path, low_memory=False)
        row = row_by_source[str(sample.source_dir)]
        for column in ["product", "branch", "modality", "gse", "sample_id", "is_co_sample", "co_dataset", "co_dataset_id"]:
            df[column] = row.get(column, "")
        parts.append(df)
    return pd.concat(parts, ignore_index=True, sort=False) if parts else pd.DataFrame()


def organize_product(
    output_root: Path,
    product: str,
    copy_mode: str,
    force: bool,
    skip_figures: bool,
    skip_integration: bool,
    integration_n_components: int,
    integration_max_umap_fit_cells: int,
    integration_clusters: int,
    integration_batch_key: str,
    integration_method: str,
    bbknn_neighbors_within_batch: int,
    bbknn_trim: int,
    leiden_resolution: float,
    rna_min_cima_l1_score: float,
    include_incomplete: bool,
    atac_workbook: Path | None = None,
) -> dict[str, object]:
    product_dir = output_root / PRODUCT_DIRS[product]
    ensure_clean_product_dir(product_dir, force)
    for name in ["manifests", "qc", "figures", "samples"]:
        (product_dir / name).mkdir(parents=True, exist_ok=True)

    samples = discover_product(output_root, product, atac_workbook=atac_workbook)
    rows = []
    cell_parts = []
    complete_samples = []
    for sample in samples:
        dest = materialize_sample(sample, product_dir, copy_mode, force)
        row = sample_row(sample, dest)
        rows.append(row)
        if row["sample_complete"] or include_incomplete:
            try:
                cell_parts.append(read_sample_metadata(sample, row))
                complete_samples.append(sample)
            except Exception as exc:
                row["sample_complete"] = False
                row["missing_files"] = f"{row['missing_files']};metadata_read_error:{exc}"

    samples_df = pd.DataFrame(rows)
    samples_df.to_csv(product_dir / "manifests" / "samples.csv", index=False)
    cells = pd.concat(cell_parts, ignore_index=True, sort=False) if cell_parts else pd.DataFrame()
    cells_path = product_dir / "manifests" / "cells_metadata.csv"
    cells.to_csv(cells_path, index=False)

    integration_status: dict[str, Any] = {
        "integration_status": "skipped_by_user" if skip_integration else "skipped_no_cells",
        "coordinate_source": "metadata_fallback",
    }
    if not skip_integration and not cells.empty:
        integration_status = integrate_product_embeddings(
            IntegrationConfig(
                product=product,
                product_dir=product_dir,
                data_root=ROOT / "data",
                n_components=integration_n_components,
                max_umap_fit_cells=integration_max_umap_fit_cells,
                n_clusters=integration_clusters,
                batch_key=integration_batch_key,
                integration_method=integration_method,
                neighbors_within_batch=bbknn_neighbors_within_batch,
                bbknn_trim=bbknn_trim,
                leiden_resolution=leiden_resolution,
                rna_min_cima_l1_score=rna_min_cima_l1_score,
            )
        )
        cells = pd.read_csv(cells_path, low_memory=False)

    write_metadata_object(cells, product_dir, product)
    write_output_inventory(product_dir)

    qc_summary = concat_existing_csv(complete_samples, "qc_summary.csv", rows)
    validation = concat_existing_csv(complete_samples, "validation_result.csv", rows)
    qc_summary.to_csv(product_dir / "qc" / "sample_qc_summary.csv", index=False)
    validation.to_csv(product_dir / "qc" / "validation_summary.csv", index=False)

    figure_status: dict[str, str] = {"status": "skipped", "detail": "--skip-figures"}
    if not skip_figures and not cells.empty:
        figure_status = render_panels(
            metadata_csv=product_dir / "manifests" / "cells_metadata.csv",
            output_dir=product_dir / "figures",
            product=product,
        )

    n_complete = 0
    n_incomplete = 0
    if not samples_df.empty:
        complete_values = pd.Series(samples_df["sample_complete"]).astype(bool).to_numpy()
        n_complete = int(np.count_nonzero(complete_values))
        n_incomplete = int(len(complete_values) - n_complete)

    status: dict[str, Any] = {
        "product": product,
        "product_dir": str(product_dir),
        "created_at": utc_now(),
        "copy_mode": copy_mode,
        "source_roots": sorted({str(sample.source_dir.parent) for sample in samples}),
        "n_samples_discovered": len(samples),
        "n_samples_complete": n_complete,
        "n_samples_incomplete": n_incomplete,
        "n_cells_metadata": len(cells),
        "figures": figure_status,
        "cross_modal_integration": False,
        "integration_status": integration_status.get("integration_status", "unknown"),
        "coordinate_source": integration_status.get("coordinate_source", "unknown"),
        "integration": integration_status,
    }
    (product_dir / "product_status.json").write_text(
        json.dumps(status, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return status


def cmd_organize_products(args: argparse.Namespace) -> int:
    output_root = Path(args.output_root)
    statuses = []
    for product in selected_products(args.products):
        integration_batch_key = args.integration_batch_key or default_integration_batch_key(product)
        integration_method = args.integration_method or default_integration_method(product)
        bbknn_neighbors_within_batch = (
            args.bbknn_neighbors_within_batch
            if args.bbknn_neighbors_within_batch is not None
            else default_bbknn_neighbors_within_batch(product)
        )
        bbknn_trim = args.bbknn_trim if args.bbknn_trim is not None else default_bbknn_trim(product)
        rna_min_cima_l1_score = (
            args.rna_min_cima_l1_score
            if args.rna_min_cima_l1_score is not None
            else default_rna_min_cima_l1_score(product)
        )
        status = organize_product(
            output_root=output_root,
            product=product,
            copy_mode=args.copy_mode,
            force=args.force,
            skip_figures=args.skip_figures,
            skip_integration=args.skip_integration,
            integration_n_components=args.integration_n_components,
            integration_max_umap_fit_cells=args.integration_max_umap_fit_cells,
            integration_clusters=args.integration_clusters,
            integration_batch_key=integration_batch_key,
            integration_method=integration_method,
            bbknn_neighbors_within_batch=bbknn_neighbors_within_batch,
            bbknn_trim=bbknn_trim,
            leiden_resolution=args.leiden_resolution,
            rna_min_cima_l1_score=rna_min_cima_l1_score,
            include_incomplete=args.include_incomplete,
        )
        statuses.append(status)
        print(
            f"[{product}] samples={status['n_samples_complete']}/{status['n_samples_discovered']} "
            f"cells={status['n_cells_metadata']} dir={status['product_dir']}"
        )
    return 0


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args() if argv is None else parse_args_from(argv)
    return cmd_organize_products(args)


def parse_args_from(argv: Iterable[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Organize annotated branch outputs into integrated product directories."
    )
    parser.add_argument("--products", default="all")
    parser.add_argument("--output-root", default=str(ROOT / "output"))
    parser.add_argument("--copy-mode", choices=["symlink", "copy"], default="symlink")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--skip-integration", action="store_true")
    parser.add_argument("--integration-n-components", type=int, default=30)
    parser.add_argument("--integration-max-umap-fit-cells", type=int, default=100_000)
    parser.add_argument("--integration-clusters", type=int, default=30)
    parser.add_argument("--integration-batch-key", default=None)
    parser.add_argument(
        "--integration-method",
        choices=["bbknn", "scanpy_neighbors", "harmony"],
        default=None,
    )
    parser.add_argument("--bbknn-neighbors-within-batch", type=int, default=None)
    parser.add_argument("--bbknn-trim", type=int, default=None)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--rna-min-cima-l1-score", type=float, default=None)
    parser.add_argument("--include-incomplete", action="store_true")
    return parser.parse_args(list(argv))


if __name__ == "__main__":
    raise SystemExit(cmd_organize_products(parse_args()))
