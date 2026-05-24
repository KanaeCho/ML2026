from __future__ import annotations

import pandas as pd


RNA_FINAL_CELLTYPE_MAPPING_VERSION = "rna_l1_to_5class_with_l2_nk_v1"
RNA_FINAL_CELLTYPES = {"CD4_T", "CD8_T", "B", "Myeloid", "NK"}


def _clean_label(value: object) -> str:
    if value is None:
        return ""
    if not isinstance(value, (str, int, float, bool)):
        return str(value).strip()
    if pd.isna(value):
        return ""
    return str(value).strip()


def infer_rna_final_celltype(raw_l1: object, raw_l2: object = None) -> str:
    """Map RNA annotations to the project 5-class final cell type vocabulary."""
    l2 = _clean_label(raw_l2)
    if "nk" in l2.lower():
        return "NK"

    l1 = _clean_label(raw_l1)
    mapping = {
        "B": "B",
        "CD4_T": "CD4_T",
        "CD8_T": "CD8_T",
        "Myeloid": "Myeloid",
        "NK": "NK",
        "ILC": "NK",
        "unconvensional_T": "CD8_T",
    }
    return mapping.get(l1, "Unknown")


def infer_rna_final_celltype_series(
    raw_l1: pd.Series,
    raw_l2: pd.Series | None = None,
) -> pd.Series:
    if raw_l2 is None:
        raw_l2 = pd.Series(pd.NA, index=raw_l1.index, dtype="object")
    else:
        raw_l2 = raw_l2.reindex(raw_l1.index)
    return pd.Series(
        [infer_rna_final_celltype(l1, l2) for l1, l2 in zip(raw_l1, raw_l2, strict=True)],
        index=raw_l1.index,
        dtype="object",
    )


def known_rna_final_celltype_mask(labels: pd.Series) -> pd.Series:
    cleaned = labels.fillna("").astype(str).str.strip()
    return cleaned.isin(RNA_FINAL_CELLTYPES)


__all__ = [
    "RNA_FINAL_CELLTYPE_MAPPING_VERSION",
    "RNA_FINAL_CELLTYPES",
    "infer_rna_final_celltype",
    "infer_rna_final_celltype_series",
    "known_rna_final_celltype_mask",
]
