from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, cast

import anndata as ad
import pandas as pd

from .final_celltype import (
    RNA_FINAL_CELLTYPE_MAPPING_VERSION,
    infer_rna_final_celltype_series,
    known_rna_final_celltype_mask,
)


def _utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _apply_final_celltype(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "azimuth_cima_l1_raw" in out.columns:
        raw_l1 = out["azimuth_cima_l1_raw"]
    elif "azimuth_cima_l1" in out.columns:
        raw_l1 = out["azimuth_cima_l1"]
        out["azimuth_cima_l1_raw"] = raw_l1.astype(str)
    elif "final_celltype" in out.columns:
        raw_l1 = out["final_celltype"]
        out["azimuth_cima_l1_raw"] = raw_l1.astype(str)
    else:
        raw_l1 = pd.Series("", index=out.index, dtype=object)
        out["azimuth_cima_l1_raw"] = ""

    if "azimuth_cell_type_l2_raw" in out.columns:
        raw_l2 = out["azimuth_cell_type_l2_raw"]
    elif "pbmcref_celltype" in out.columns:
        raw_l2 = out["pbmcref_celltype"]
    elif "azimuth_cell_type" in out.columns:
        raw_l2 = out["azimuth_cell_type"]
    else:
        raw_l2 = pd.Series(pd.NA, index=out.index, dtype=object)

    out["final_celltype"] = infer_rna_final_celltype_series(raw_l1, raw_l2).astype(str)
    out["azimuth_cima_l1"] = out["final_celltype"].astype(str)
    out["final_celltype_mapping"] = RNA_FINAL_CELLTYPE_MAPPING_VERSION
    return out


def _known_mask(frame: pd.DataFrame) -> pd.Series:
    return known_rna_final_celltype_mask(frame["final_celltype"])


def _read_csv(path: Path) -> pd.DataFrame | None:
    if not path.exists():
        return None
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _update_qc_summary(path: Path, n_final: int, n_removed: int) -> None:
    if not path.exists():
        return
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if len(frame) == 0:
        return
    frame.loc[0, "n_cells_final_output"] = str(n_final)
    frame.loc[0, "n_cells_unknown_final_celltype_removed"] = str(n_removed)
    frame.loc[0, "final_celltype_mapping"] = RNA_FINAL_CELLTYPE_MAPPING_VERSION
    frame.to_csv(path, index=False)


def _update_validation(path: Path, n_final: int, n_removed: int) -> None:
    if not path.exists():
        return
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    check_name = frame["check_name"] if "check_name" in frame.columns else pd.Series("", index=frame.index)
    mask = check_name.eq("metadata_qc_known_final_celltype_only")
    if bool(mask.any()):
        frame.loc[mask, "detail"] = (
            f"n_cells_final_output={n_final};"
            f"unknown_final_celltype_removed={n_removed}"
        )
    frame.to_csv(path, index=False)


def _update_status(path: Path, payload: dict[str, Any]) -> None:
    existing: dict[str, Any] = {}
    if path.exists():
        existing = json.loads(path.read_text(encoding="utf-8"))
    existing["final_celltype_backfill"] = payload
    path.write_text(json.dumps(existing, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _backfill_h5ad(path: Path) -> tuple[int, int]:
    adata = ad.read_h5ad(path)
    obs = _apply_final_celltype(cast(pd.DataFrame, adata.obs.copy()))
    known = _known_mask(obs)
    removed = int((~known).sum())
    adata = adata[known.to_numpy()].copy()
    obs_out = obs.loc[known.to_numpy()].copy()
    obs_out.index = pd.Index([str(value) for value in obs_out.index.tolist()], dtype=object)
    obs_out.columns = pd.Index([str(value) for value in obs_out.columns.tolist()], dtype=object)
    for column in obs_out.columns:
        series = obs_out[column]
        if (
            pd.api.types.is_string_dtype(series.dtype)
            or pd.api.types.is_object_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            obs_out[column] = pd.Series(
                ["" if pd.isna(value) else str(value) for value in series.tolist()],
                index=series.index,
                dtype=object,
            )
    adata.obs = obs_out
    var_out = cast(pd.DataFrame, adata.var.copy())
    var_out.index = pd.Index([str(value) for value in var_out.index.tolist()], dtype=object)
    var_out.columns = pd.Index([str(value) for value in var_out.columns.tolist()], dtype=object)
    for column in var_out.columns:
        series = var_out[column]
        if (
            pd.api.types.is_string_dtype(series.dtype)
            or pd.api.types.is_object_dtype(series.dtype)
            or isinstance(series.dtype, pd.CategoricalDtype)
        ):
            var_out[column] = pd.Series(
                ["" if pd.isna(value) else str(value) for value in series.tolist()],
                index=series.index,
                dtype=object,
            )
    adata.var = var_out
    adata.write_h5ad(path, convert_strings_to_categoricals=False)
    return int(adata.n_obs), removed


def _sample_dirs(root: Path) -> list[Path]:
    if not root.exists():
        raise FileNotFoundError(f"Backfill root not found: {root}")
    return sorted(path for path in root.iterdir() if path.is_dir())


def backfill_sample_dir(sample_dir: Path, *, dry_run: bool = False) -> dict[str, Any]:
    metadata_path = sample_dir / "metadata.csv"
    metadata_qc_path = sample_dir / "metadata_qc.csv"
    sample_id = sample_dir.name
    h5ad_path = sample_dir / f"{sample_id}.h5ad"

    metadata = _read_csv(metadata_path)
    metadata_qc = _read_csv(metadata_qc_path)
    if metadata is None or metadata_qc is None or not h5ad_path.exists():
        return {
            "sample_dir": str(sample_dir),
            "status": "skipped",
            "detail": "missing metadata.csv, metadata_qc.csv, or h5ad",
        }

    metadata_new = _apply_final_celltype(metadata)
    metadata_qc_new = _apply_final_celltype(metadata_qc)
    known_qc = _known_mask(metadata_qc_new)
    metadata_qc_new = metadata_qc_new.loc[known_qc.to_numpy()].reset_index(drop=True)
    removed_qc = int((~known_qc).sum())

    result: dict[str, Any] = {
        "sample_dir": str(sample_dir),
        "status": "dry_run" if dry_run else "updated",
        "metadata_cells": int(len(metadata_new)),
        "metadata_qc_cells": int(len(metadata_qc_new)),
        "metadata_qc_removed_unknown_final_celltype": removed_qc,
    }
    if dry_run:
        result["final_celltype_counts"] = (
            metadata_qc_new["final_celltype"].value_counts().sort_index().to_dict()
        )
        return result

    metadata_new.to_csv(metadata_path, index=False)
    metadata_qc_new.to_csv(metadata_qc_path, index=False)
    h5ad_n, h5ad_removed = _backfill_h5ad(h5ad_path)
    _update_qc_summary(sample_dir / "qc_summary.csv", h5ad_n, h5ad_removed)
    _update_validation(sample_dir / "validation_result.csv", h5ad_n, h5ad_removed)
    _update_status(
        sample_dir / "run_status.json",
        {
            "version": RNA_FINAL_CELLTYPE_MAPPING_VERSION,
            "updated_at": _utc_now(),
            "h5ad_cells": h5ad_n,
            "h5ad_removed_unknown_final_celltype": h5ad_removed,
        },
    )
    result["h5ad_cells"] = h5ad_n
    result["h5ad_removed_unknown_final_celltype"] = h5ad_removed
    return result


def cmd_backfill_rna_final_celltype(args) -> int:
    root = Path(args.root)
    rows = [backfill_sample_dir(sample_dir, dry_run=args.dry_run) for sample_dir in _sample_dirs(root)]
    print(pd.DataFrame(rows).to_string(index=False))
    failed = [row for row in rows if row.get("status") == "skipped"]
    return 1 if failed and getattr(args, "fail_on_skip", False) else 0


__all__ = ["cmd_backfill_rna_final_celltype", "backfill_sample_dir"]
