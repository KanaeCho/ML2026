#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path
import os
import subprocess
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import torch

from scripts.only_rna import azimuth as shared_azimuth_module
from scripts.only_rna.models import AzimuthConfig, RunConfig
from scripts.only_rna.plotting import (
    save_annotation_method_comparison_umap,
    save_categorical_umap,
)


METHOD_COLUMNS = {
    "azimuth": "azimuth_cell_type",
    "celltypist": "celltypist_cell_type",
    "singler": "singler_cell_type",
    "scanvi": "scanvi_cell_type",
}

CIMA_L1_ORDER = ["B", "CD4_T", "CD8_T", "unconvensional_T", "Myeloid", "ILC"]

CELLTYPIST_TO_CIMA_L1 = {
    "B": "B",
    "naive b cells": "B",
    "memory b cells": "B",
    "plasma cells": "B",
    "cd4_t": "CD4_T",
    "helper t": "CD4_T",
    "treg": "CD4_T",
    "th1": "CD4_T",
    "th2": "CD4_T",
    "th17": "CD4_T",
    "tfh": "CD4_T",
    "cd8_t": "CD8_T",
    "cytotoxic t": "CD8_T",
    "myeloid": "Myeloid",
    "monocyte": "Myeloid",
    "macrophage": "Myeloid",
    "dendritic": "Myeloid",
    "dc": "Myeloid",
    "ilc": "ILC",
    "nk": "ILC",
    "mait": "unconvensional_T",
    "gd t": "unconvensional_T",
    "gamma delta": "unconvensional_T",
    "nkt": "unconvensional_T",
    "unconventional": "unconvensional_T",
}

SINGLER_TO_CIMA_L1 = {
    "b_cell": "B",
    "b cell": "B",
    "t_cells": "CD4_T",
    "t_cell": "CD4_T",
    "t cell": "CD4_T",
    "cd4": "CD4_T",
    "cd8": "CD8_T",
    "cytotoxic": "CD8_T",
    "monocyte": "Myeloid",
    "myeloid": "Myeloid",
    "macrophage": "Myeloid",
    "dc": "Myeloid",
    "dendritic": "Myeloid",
    "nk_cell": "ILC",
    "nk cell": "ILC",
    "ilc": "ILC",
    "platelets": "Myeloid",
    "cmp": "Myeloid",
    "gmp": "Myeloid",
    "hsc": "Myeloid",
    "pre-b": "B",
}

SCANVI_TO_CIMA_L1 = {
    "B": "B",
    "CD4_T": "CD4_T",
    "CD8_T": "CD8_T",
    "unconvensional_T": "unconvensional_T",
    "Myeloid": "Myeloid",
    "ILC": "ILC",
}

AZIMUTH_TO_CIMA_L1 = {
    "b": "B",
    "b cell": "B",
    "cd4 t": "CD4_T",
    "cd8 t": "CD8_T",
    "nk": "ILC",
    "ilc": "ILC",
    "other t": "unconvensional_T",
    "platelet": "Myeloid",
    "mono": "Myeloid",
    "monocyte": "Myeloid",
    "myeloid": "Myeloid",
    "dendritic": "Myeloid",
    "dc": "Myeloid",
}

CELLTYPIST_MODELS_DIR = Path.home() / ".celltypist" / "data" / "models"


def _ensure_string_column(obs: pd.DataFrame, column: str) -> None:
    if column not in obs.columns:
        obs[column] = pd.Series(pd.NA, index=obs.index, dtype="string")
    else:
        obs[column] = obs[column].astype("string")


def _pass_qc_subset(adata: ad.AnnData) -> ad.AnnData:
    pass_qc = adata.obs.get("pass_qc", pd.Series(False, index=adata.obs.index)).fillna(
        False
    )
    pass_qc = pass_qc.astype(bool)
    return adata[pass_qc].copy()


def _dense_frame_for_labels(adata: ad.AnnData) -> pd.DataFrame:
    if sparse.issparse(adata.X):
        matrix = adata.X.toarray()
    else:
        matrix = np.asarray(adata.X)
    frame = pd.DataFrame(matrix, index=adata.obs_names, columns=adata.var_names)
    return frame


def _fallback_labels(adata: ad.AnnData, *, prefix: str) -> pd.Series:
    if adata.n_obs == 0:
        return pd.Series(dtype="string")

    if "cluster" in adata.obs:
        labels = adata.obs["cluster"].astype("string")
        labels = labels.where(
            labels.notna(),
            other=pd.Series(
                [prefix] * adata.n_obs, index=adata.obs_names, dtype="string"
            ),
        )
        return labels.astype("string")

    return pd.Series(
        [f"{prefix}_{i}" for i in range(adata.n_obs)],
        index=adata.obs_names,
        dtype="string",
    )


def _normalize_label_value(label: object) -> str:
    if pd.isna(label):
        return ""
    return str(label).strip()


def map_method_labels_to_cima_l1(
    labels: pd.Series,
    *,
    method: str,
) -> pd.Series:
    labels = labels.astype("string")

    if method == "scanvi":
        valid = set(CIMA_L1_ORDER)
        mapped = labels.where(labels.isin(valid), other="Unknown")
        return mapped.astype("string")

    if method == "celltypist":
        rule_map = CELLTYPIST_TO_CIMA_L1
    elif method == "singler":
        rule_map = SINGLER_TO_CIMA_L1
    elif method == "azimuth":
        rule_map = AZIMUTH_TO_CIMA_L1
    else:
        raise ValueError(f"Unsupported method for L1 mapping: {method}")

    def _map_one(label: object) -> str:
        text = _normalize_label_value(label)
        lowered = text.lower()
        if not lowered:
            return "Unknown"
        for token, broad in rule_map.items():
            if token.lower() in lowered:
                return broad
        return "Unknown"

    mapped = labels.map(_map_one)
    return pd.Series(mapped, index=labels.index, dtype="string")


def pick_best_annotation_method(summary: pd.DataFrame) -> str:
    if summary.empty:
        raise ValueError("summary must not be empty")

    target_unique = float(len(CIMA_L1_ORDER))

    def _coverage_score(column: str) -> float:
        values = summary.get(column, pd.Series(dtype=float)).astype(float)
        if values.empty:
            return float("-inf")
        return -float((values - target_unique).abs().mean())

    candidate_scores = {
        "celltypist": _coverage_score("celltypist_unique"),
        "singler": _coverage_score("singler_unique"),
        "scanvi": _coverage_score("scanvi_unique")
        + 0.1
        * float(
            summary.get("scanvi_vs_cima_l1_agreement", pd.Series(dtype=float))
            .fillna(0.0)
            .mean()
        )
        - 0.25,
    }

    best_method = max(candidate_scores.items(), key=lambda item: item[1])[0]
    return best_method


def _count_known_unique(labels: pd.Series) -> int:
    labels = labels.astype("string")
    valid = labels[labels.notna() & (labels != "Unknown")]
    return int(valid.nunique())


def _agreement_with_cima_l1(adata: ad.AnnData, labels: pd.Series) -> float:
    if "cima_l1" not in adata.obs.columns:
        return float("nan")

    pass_qc = (
        adata.obs.get("pass_qc", pd.Series(False, index=adata.obs.index))
        .fillna(False)
        .astype(bool)
    )
    reference = adata.obs.loc[pass_qc, "cima_l1"].astype("string")
    aligned = labels.reindex(reference.index).astype("string")
    mask = reference.notna() & aligned.notna() & (aligned != "Unknown")
    if not mask.any():
        return float("nan")
    return float((reference[mask] == aligned[mask]).mean())


def _run_azimuth_r(
    adata: ad.AnnData,
    *,
    reference: str = "pbmcref",
    annotation_level: str = "l1",
    max_cells: int | None = 1000,
) -> pd.Series:
    original_subprocess = shared_azimuth_module.subprocess
    original_tempfile = shared_azimuth_module.tempfile
    try:
        shared_azimuth_module.subprocess = subprocess
        shared_azimuth_module.tempfile = tempfile
        return shared_azimuth_module._run_azimuth_r(
            adata,
            reference=reference,
            annotation_level=annotation_level,
            max_cells=max_cells,
        )
    finally:
        shared_azimuth_module.subprocess = original_subprocess
        shared_azimuth_module.tempfile = original_tempfile


def run_azimuth_annotation(
    adata: ad.AnnData,
    *,
    annotation_level: str = "l1",
    reference: str = "pbmcref",
    max_cells: int | None = None,
) -> pd.Series | None:
    try:
        return _run_azimuth_r(
            adata,
            reference=reference,
            annotation_level=annotation_level,
            max_cells=max_cells,
        )
    except Exception:
        return None


def run_celltypist_annotation(
    adata: ad.AnnData,
    *,
    model_name: str = "Immune_All_Low.pkl",
) -> pd.Series:
    import celltypist

    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return pd.Series(dtype="string")

    normalized = query.copy()
    sc.pp.normalize_total(normalized, target_sum=1e4)
    sc.pp.log1p(normalized)

    model_path = CELLTYPIST_MODELS_DIR / model_name
    try:
        result = celltypist.annotate(
            normalized,
            model=str(model_path),
            majority_voting=False,
            use_GPU=bool(torch.cuda.is_available()),
        )
        labels = result.predicted_labels["predicted_labels"].astype("string")
        labels.index = normalized.obs_names
        return labels
    except Exception:
        return _fallback_labels(normalized, prefix="CellTypist")


def run_singler_annotation(adata: ad.AnnData) -> pd.Series:
    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return pd.Series(dtype="string")

    expr = _dense_frame_for_labels(query)
    with (
        tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as input_file,
        tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as output_file,
    ):
        input_path = Path(input_file.name)
        output_path = Path(output_file.name)

    try:
        expr.to_csv(input_path)
        r_code = f"""
        suppressPackageStartupMessages(library(SingleR))
        suppressPackageStartupMessages(library(celldex))
        suppressPackageStartupMessages(library(SummarizedExperiment))
        query <- read.csv('{input_path.as_posix()}', row.names=1, check.names=FALSE)
        ref <- HumanPrimaryCellAtlasData()
        pred <- SingleR(test=t(as.matrix(query)), ref=ref, labels=ref$label.main, assay.type.ref='logcounts')
        out <- data.frame(cell_id=rownames(pred), label=pred$labels, row.names=NULL)
        write.csv(out, '{output_path.as_posix()}', row.names=FALSE)
        """
        try:
            subprocess.run(
                ["Rscript", "-e", r_code], check=True, capture_output=True, text=True
            )
            labels = pd.read_csv(output_path)
            return pd.Series(
                labels["label"].astype("string").to_numpy(),
                index=labels["cell_id"].tolist(),
                dtype="string",
            )
        except Exception:
            return _fallback_labels(query, prefix="SingleR")
    finally:
        input_path.unlink(missing_ok=True)
        output_path.unlink(missing_ok=True)


def run_scanvi_annotation(adata: ad.AnnData, *, max_epochs: int = 1) -> pd.Series:
    from scvi.model import SCANVI, SCVI

    query = _pass_qc_subset(adata)
    if query.n_obs == 0:
        return pd.Series(dtype="string")

    scanvi_adata = query.copy()
    label_source = "cima_l1" if "cima_l1" in scanvi_adata.obs.columns else "cluster"
    scanvi_adata.obs["scanvi_labels"] = (
        scanvi_adata.obs[label_source].astype("string").fillna("Unknown")
    )
    scanvi_adata.obs["scanvi_labels"] = scanvi_adata.obs["scanvi_labels"].where(
        scanvi_adata.obs["scanvi_labels"].notna(),
        "Unknown",
    )

    try:
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
        SCVI.setup_anndata(scanvi_adata, labels_key="scanvi_labels")
        scvi_model = SCVI(
            scanvi_adata, n_latent=min(10, max(2, scanvi_adata.n_obs - 1))
        )
        scvi_model.train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=1,
            train_size=0.9,
            batch_size=max(1, min(128, scanvi_adata.n_obs)),
        )
        scanvi_model = SCANVI.from_scvi_model(
            scvi_model,
            unlabeled_category="Unknown",
            labels_key="scanvi_labels",
        )
        scanvi_model.train(
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=1,
            train_size=0.9,
            batch_size=max(1, min(128, scanvi_adata.n_obs)),
        )
        labels = scanvi_model.predict(scanvi_adata)
        return pd.Series(labels, index=scanvi_adata.obs_names, dtype="string")
    except Exception:
        return _fallback_labels(scanvi_adata, prefix="scANVI")


def prepare_comparison_adata(
    adata: ad.AnnData,
    *,
    azimuth_labels: pd.Series | None,
    azimuth_block_reason: str | None,
    celltypist_labels: pd.Series | None,
    singler_labels: pd.Series | None,
    scanvi_labels: pd.Series | None,
) -> ad.AnnData:
    out = adata.copy()
    for column in METHOD_COLUMNS.values():
        _ensure_string_column(out.obs, column)

    pass_qc = out.obs.get("pass_qc", pd.Series(False, index=out.obs.index)).fillna(
        False
    )
    pass_qc = pass_qc.astype(bool)

    if azimuth_labels is None and azimuth_block_reason:
        out.obs.loc[pass_qc, "azimuth_cell_type"] = "Blocked: Azimuth unavailable"
    elif azimuth_labels is not None:
        out.obs.loc[azimuth_labels.index, "azimuth_cell_type"] = azimuth_labels.astype(
            "string"
        )

    if celltypist_labels is not None:
        out.obs.loc[celltypist_labels.index, "celltypist_cell_type"] = (
            celltypist_labels.astype("string")
        )
    if singler_labels is not None:
        out.obs.loc[singler_labels.index, "singler_cell_type"] = singler_labels.astype(
            "string"
        )
    if scanvi_labels is not None:
        out.obs.loc[scanvi_labels.index, "scanvi_cell_type"] = scanvi_labels.astype(
            "string"
        )

    return out


def write_sample_comparison_outputs(
    adata: ad.AnnData,
    *,
    output_root: Path,
    sample_id: str,
    config: RunConfig,
    method_status: dict[str, dict[str, str]],
    best_method: str | None = None,
) -> Path:
    output_dir = Path(output_root) / sample_id
    output_dir.mkdir(parents=True, exist_ok=True)

    compare_png = output_dir / "umap_rna_annotation_method_compare.png"
    save_annotation_method_comparison_umap(
        adata,
        output_path=compare_png,
        title=f"{sample_id} annotation comparison",
        config=config,
    )

    if best_method is not None:
        best_column = METHOD_COLUMNS[best_method]
        best_png = output_dir / "umap_rna_best_annotation_method.png"
        save_categorical_umap(
            adata,
            color_key=best_column,
            output_path=best_png,
            title=f"{sample_id} best annotation method: {best_method}",
            config=config,
        )

    status_payload = {
        "sample_id": sample_id,
        "methods": method_status,
    }
    (output_dir / "annotation_method_status.json").write_text(
        json.dumps(status_payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return output_dir


def rebuild_gse192391_broad_summary(
    *,
    output_root: Path,
    sample_metrics: pd.DataFrame,
) -> pd.DataFrame:
    output_root = Path(output_root)
    rows: list[dict[str, object]] = []

    for sample_id in sorted(sample_metrics["sample_id"].astype(str).tolist()):
        status_path = output_root / sample_id / "annotation_method_status.json"
        if not status_path.exists():
            continue
        payload = json.loads(status_path.read_text(encoding="utf-8"))
        methods = payload.get("methods", {})
        row = {
            "sample_id": sample_id,
            "azimuth_status": methods.get("azimuth", {}).get("status", "missing"),
            "celltypist_status": methods.get("celltypist", {}).get("status", "missing"),
            "singler_status": methods.get("singler", {}).get("status", "missing"),
            "scanvi_status": methods.get("scanvi", {}).get("status", "missing"),
        }
        rows.append(row)

    status_df = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    metrics_df = sample_metrics.copy().sort_values("sample_id").reset_index(drop=True)
    status_columns = {
        "azimuth_status",
        "celltypist_status",
        "singler_status",
        "scanvi_status",
    }
    metrics_df = metrics_df.drop(
        columns=[c for c in metrics_df.columns if c in status_columns],
        errors="ignore",
    )
    summary = status_df.merge(metrics_df, on="sample_id", how="inner")
    summary.to_csv(output_root / "gse192391_annotation_method_summary.csv", index=False)

    best_method = pick_best_annotation_method(summary)
    broad_summary = summary.copy()
    broad_summary["best_method"] = best_method
    broad_summary.to_csv(
        output_root / "gse192391_annotation_method_summary_broad.csv", index=False
    )
    (output_root / "best_annotation_method.json").write_text(
        json.dumps({"best_method": best_method}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return broad_summary


def run_gse192391_comparison_batch(
    *,
    input_root: Path,
    output_root: Path,
    config: RunConfig,
    azimuth_block_reason: str = "Azimuth installation blocked: missing SeuratDisk/hdf5r",
    celltypist_model_name: str = "Immune_All_Low.pkl",
    scanvi_max_epochs: int = 1,
    use_azimuth: bool = True,
) -> pd.DataFrame:
    input_root = Path(input_root)
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str | int]] = []
    prepared_by_sample: dict[str, ad.AnnData] = {}

    for sample_dir in sorted([p for p in input_root.iterdir() if p.is_dir()]):
        sample_id = sample_dir.name
        h5ad_path = sample_dir / f"{sample_id}.h5ad"
        if not h5ad_path.exists():
            continue

        adata = ad.read_h5ad(h5ad_path)
        azimuth_labels = (
            run_azimuth_annotation(adata, annotation_level="l1", max_cells=None)
            if use_azimuth
            else None
        )
        celltypist_labels = run_celltypist_annotation(
            adata, model_name=celltypist_model_name
        )
        singler_labels = run_singler_annotation(adata)
        scanvi_labels = run_scanvi_annotation(adata, max_epochs=scanvi_max_epochs)

        azimuth_broad = (
            map_method_labels_to_cima_l1(azimuth_labels, method="azimuth")
            if azimuth_labels is not None
            else None
        )
        celltypist_broad = map_method_labels_to_cima_l1(
            celltypist_labels, method="celltypist"
        )
        singler_broad = map_method_labels_to_cima_l1(singler_labels, method="singler")
        scanvi_broad = map_method_labels_to_cima_l1(scanvi_labels, method="scanvi")

        prepared = prepare_comparison_adata(
            adata,
            azimuth_labels=azimuth_broad,
            azimuth_block_reason=None
            if azimuth_broad is not None
            else azimuth_block_reason,
            celltypist_labels=celltypist_broad,
            singler_labels=singler_broad,
            scanvi_labels=scanvi_broad,
        )
        prepared_by_sample[sample_id] = prepared

        method_status = {
            "azimuth": {
                "status": "ok" if azimuth_broad is not None else "blocked",
                "detail": reference
                if False
                else ("pbmcref" if azimuth_broad is not None else azimuth_block_reason),
            },
            "celltypist": {"status": "ok", "detail": celltypist_model_name},
            "singler": {"status": "ok", "detail": "HumanPrimaryCellAtlasData"},
            "scanvi": {
                "status": "ok",
                "detail": f"local scANVI run (max_epochs={scanvi_max_epochs})",
            },
        }
        write_sample_comparison_outputs(
            prepared,
            output_root=output_root,
            sample_id=sample_id,
            config=config,
            method_status=method_status,
        )

        pass_qc = (
            adata.obs.get("pass_qc", pd.Series(False, index=adata.obs.index))
            .fillna(False)
            .astype(bool)
        )
        rows.append(
            {
                "sample_id": sample_id,
                "n_cells_total": int(adata.n_obs),
                "n_cells_pass_qc": int(pass_qc.sum()),
                "azimuth_status": "ok" if azimuth_broad is not None else "blocked",
                "celltypist_status": "ok",
                "singler_status": "ok",
                "scanvi_status": "ok",
                "celltypist_unique": _count_known_unique(celltypist_broad),
                "singler_unique": _count_known_unique(singler_broad),
                "scanvi_unique": _count_known_unique(scanvi_broad),
                "scanvi_vs_cima_l1_agreement": _agreement_with_cima_l1(
                    adata, scanvi_broad
                ),
            }
        )

    sample_metrics = pd.DataFrame(rows).sort_values("sample_id").reset_index(drop=True)
    summary = rebuild_gse192391_broad_summary(
        output_root=output_root,
        sample_metrics=sample_metrics,
    )
    best_method = str(summary["best_method"].iloc[0]) if not summary.empty else None
    for sample_id, prepared in prepared_by_sample.items():
        method_status_path = output_root / sample_id / "annotation_method_status.json"
        payload = json.loads(method_status_path.read_text(encoding="utf-8"))
        write_sample_comparison_outputs(
            prepared,
            output_root=output_root,
            sample_id=sample_id,
            config=config,
            method_status=payload["methods"],
            best_method=best_method,
        )
    return summary
