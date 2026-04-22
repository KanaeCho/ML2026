from __future__ import annotations

import pandas as pd


_L1_DIRECT_MAP: dict[str, str] = {
    "B": "B",
    "CD4 T": "CD4_T",
    "CD8 T": "CD8_T",
    "DC": "Myeloid",
    "Mono": "Myeloid",
    "NK": "ILC",
    "other T": "unconvensional_T",
}

_L2_FALLBACK_RULES: list[tuple[list[str], str]] = [
    (["cd4", "treg", "tfh"], "CD4_T"),
    (["cd8", "cytotoxic"], "CD8_T"),
    (["mono", "dc", "cdc", "pdc", "macrophage"], "Myeloid"),
    (["nk", "ilc"], "ILC"),
    (["mait", "nkt", "gdt", "γδt", "gamma delta", "gamma-delta"], "unconvensional_T"),
    (["b cell", "memory b", "naive b", "plasma", "plasmablast", "b"], "B"),
]


def _align_one(raw_l1: object, raw_l2: object) -> str:
    if not pd.isna(raw_l2):
        lowered = str(raw_l2).strip().lower()
        for keywords, aligned in _L2_FALLBACK_RULES:
            if any(keyword in lowered for keyword in keywords):
                return aligned
        return "Unknown"

    if not pd.isna(raw_l1):
        direct = _L1_DIRECT_MAP.get(str(raw_l1).strip())
        if direct is not None:
            return direct

    return "Unknown"


def align_pbmcref_to_cima_l1(
    raw_l1: pd.Series | None,
    raw_l2: pd.Series | None,
) -> tuple[pd.Series, pd.Series]:
    if raw_l1 is None and raw_l2 is None:
        empty_index = pd.Index([], dtype=object)
        return (
            pd.Series([], index=empty_index, dtype="string"),
            pd.Series([], index=empty_index, dtype="boolean"),
        )

    if raw_l1 is None:
        assert raw_l2 is not None
        index = raw_l2.index
        raw_l1 = pd.Series(pd.NA, index=index, dtype="string")
    elif raw_l2 is None:
        index = raw_l1.index
        raw_l2 = pd.Series(pd.NA, index=index, dtype="string")
    else:
        index = raw_l1.index.union(raw_l2.index)
        raw_l1 = raw_l1.reindex(index)
        raw_l2 = raw_l2.reindex(index)

    aligned = pd.Series(
        [_align_one(raw_l1.loc[idx], raw_l2.loc[idx]) for idx in index],
        index=index,
        dtype="string",
    )
    unmapped = pd.Series(aligned == "Unknown", index=index, dtype="boolean")
    return aligned, unmapped


__all__ = ["align_pbmcref_to_cima_l1"]
