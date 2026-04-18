from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex, to_rgb

from .models import RunConfig

# ---------------------------------------------------------------------------
# Hierarchical CIMA palette
# ---------------------------------------------------------------------------

BASE_L1_PALETTE: dict[str, str] = {
    "B": "#2C7BB6",
    "CD4_T": "#D7191C",
    "CD8_T": "#E8521A",
    "unconvensional_T": "#FDAE61",
    "Myeloid": "#1A9641",
    "ILC": "#762A83",
    "Unknown": "0.74",
}

FALLBACK_COLORS: list[str] = [
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#FEE08B",
    "#F46D43",
    "#A50026",
]


def _mix_color(source: str, target: str, amount: float) -> str:
    amount = min(max(amount, 0.0), 1.0)
    src_rgb = np.asarray(to_rgb(source), dtype=np.float32)
    tgt_rgb = np.asarray(to_rgb(target), dtype=np.float32)
    blended = src_rgb * (1 - amount) + tgt_rgb * amount
    return to_hex(
        (float(blended[0]), float(blended[1]), float(blended[2])),
        keep_alpha=False,
    )


def _make_shade_palette(base_color: str, labels: list[str]) -> dict[str, str]:
    labels = list(dict.fromkeys(labels))
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: base_color}
    colors: list[str] = []
    for frac in np.linspace(0, 1, len(labels)):
        if frac <= 0.5:
            amount = 0.55 * (1 - frac / 0.5)
            colors.append(_mix_color(base_color, "#FFFFFF", amount))
        else:
            amount = 0.18 * ((frac - 0.5) / 0.5)
            colors.append(_mix_color(base_color, "#000000", amount))
    return dict(zip(labels, colors, strict=False))


def build_cima_l1_palette(l1_labels: list[str]) -> dict[str, str]:
    palette: dict[str, str] = dict(BASE_L1_PALETTE)
    missing = [label for label in l1_labels if label not in palette]
    for idx, label in enumerate(missing):
        palette[label] = FALLBACK_COLORS[idx % len(FALLBACK_COLORS)]
    return {label: palette[label] for label in l1_labels}


def build_cima_l2_palette(l2_by_l1: dict[str, list[str]]) -> dict[str, str]:
    l1_labels = sorted(l2_by_l1.keys())
    l1_palette = build_cima_l1_palette(l1_labels)
    palette: dict[str, str] = {}
    for l1 in l1_labels:
        children = sorted(set(l2_by_l1[l1]))
        palette.update(_make_shade_palette(l1_palette[l1], children))
    return palette


def _infer_l2_parent(label: str) -> str:
    prefix = label.split("_")[0]
    if prefix in ("CD4", "CD8"):
        return prefix + "_T"
    if prefix == "ILC":
        return "ILC"
    if label.endswith("_B") or label.endswith(" B"):
        return "B"
    if prefix in ("gd", "MAIT", "NKT", "unconvensional") or label in (
        "gd_T",
        "MAIT",
        "NKT",
    ):
        return "unconvensional_T"
    if prefix in ("Mono", "cDC", "pDC", "Macro", "Myeloid"):
        return "Myeloid"
    return prefix


def get_hierarchical_palette(
    color_key: str, categories: list[str]
) -> dict[str, tuple[float, ...]] | None:
    if color_key == "cluster":
        return None
    if color_key in (
        "cima_l1",
        "cima_l1_masked",
        "azimuth_cell_type",
        "celltypist_cell_type",
        "singler_cell_type",
        "scanvi_cell_type",
    ):
        hex_palette = build_cima_l1_palette(categories)
    elif color_key == "cima_l2":
        l2_by_l1: dict[str, list[str]] = {}
        for label in categories:
            parent = _infer_l2_parent(label)
            l2_by_l1.setdefault(parent, []).append(label)
        hex_palette = build_cima_l2_palette(l2_by_l1)
    else:
        return None
    return {label: to_rgb(hex_palette[label]) for label in categories}


# ---------------------------------------------------------------------------


def _plot_frame(adata: ad.AnnData, color_key: str) -> pd.DataFrame:
    if color_key not in adata.obs.columns:
        raise KeyError(f"adata.obs must contain '{color_key}'")

    if "umap_1" not in adata.obs.columns or "umap_2" not in adata.obs.columns:
        return pd.DataFrame(columns=["umap_1", "umap_2", color_key])

    frame = pd.DataFrame(
        {
            "umap_1": pd.to_numeric(adata.obs.get("umap_1"), errors="coerce"),
            "umap_2": pd.to_numeric(adata.obs.get("umap_2"), errors="coerce"),
            color_key: adata.obs[color_key].astype("string"),
        },
        index=adata.obs_names,
    )
    return frame.dropna(subset=["umap_1", "umap_2", color_key])


def save_categorical_umap(
    adata: ad.AnnData,
    color_key: str,
    output_path: Path,
    title: str,
    config: RunConfig,
) -> None:
    plotting = config.plotting
    display_point_size = max(
        float(plotting.point_size) * 1.1, float(plotting.point_size) + 1.0
    )
    legend_location = getattr(plotting, "legend_location", "center left")
    legend_ncols = int(getattr(plotting, "legend_ncols", 0))
    legend_bbox_to_anchor = getattr(plotting, "legend_bbox_to_anchor", (1.02, 0.5))
    legend_markerscale = float(getattr(plotting, "legend_markerscale", 4.0))
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frame = _plot_frame(adata, color_key)
    categories: list[str] = []
    if not frame.empty:
        categories = sorted(frame[color_key].astype(str).unique().tolist())
        if legend_ncols <= 0:
            if color_key == "cima_l2" and len(categories) >= 18:
                legend_ncols = 2
            elif len(categories) >= 24:
                legend_ncols = 3
            elif len(categories) >= 12:
                legend_ncols = 2
            else:
                legend_ncols = 1

    figure_width: float = float(plotting.umap_width)
    if categories:
        figure_width = figure_width + max(1.5, 0.9 * float(legend_ncols))

    fig, ax = plt.subplots(
        figsize=(figure_width, plotting.umap_height),
        dpi=plotting.dpi,
    )
    ax.set_box_aspect(1)

    if frame.empty:
        ax.set_title(title)
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.text(
            0.5, 0.5, "No valid cells", ha="center", va="center", transform=ax.transAxes
        )
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        hier_palette = get_hierarchical_palette(color_key, categories)
        if hier_palette is not None:
            for category in categories:
                subset = frame.loc[frame[color_key].astype(str) == category]
                ax.scatter(
                    subset["umap_1"],
                    subset["umap_2"],
                    s=display_point_size,
                    c=[hier_palette[category]],
                    label=category,
                    linewidths=0,
                    alpha=0.9,
                )
        else:
            cmap = plt.get_cmap("tab20", max(len(categories), 1))
            for idx, category in enumerate(categories):
                subset = frame.loc[frame[color_key].astype(str) == category]
                ax.scatter(
                    subset["umap_1"],
                    subset["umap_2"],
                    s=display_point_size,
                    c=[cmap(idx)],
                    label=category,
                    linewidths=0,
                    alpha=0.9,
                )

        legend = ax.legend(
            title=color_key,
            loc=legend_location,
            bbox_to_anchor=legend_bbox_to_anchor,
            ncol=legend_ncols,
            frameon=False,
            fontsize=plotting.legend_fontsize,
            title_fontsize=plotting.legend_title_fontsize,
            markerscale=legend_markerscale,
        )
        if legend is not None:
            legend._legend_box.align = "left"

        ax.set_title(title)
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_aspect("auto")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0))
    fig.savefig(output_path, dpi=plotting.dpi, format="png")
    plt.close(fig)


def save_annotation_method_comparison_umap(
    adata: ad.AnnData,
    output_path: Path,
    title: str,
    config: RunConfig,
) -> None:
    method_columns = [
        ("azimuth_cell_type", "Azimuth"),
        ("celltypist_cell_type", "CellTypist"),
        ("singler_cell_type", "SingleR"),
        ("scanvi_cell_type", "scANVI"),
    ]
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(
            float(config.plotting.umap_width) * 2.0,
            float(config.plotting.umap_height) * 2.0,
        ),
        dpi=config.plotting.dpi,
    )
    axes_array = np.atleast_1d(axes).ravel()

    for ax, (color_key, panel_title) in zip(axes_array, method_columns, strict=True):
        frame = (
            _plot_frame(adata, color_key)
            if color_key in adata.obs.columns
            else pd.DataFrame()
        )
        ax.set_title(panel_title)
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_box_aspect(1)

        if frame.empty:
            ax.text(
                0.5,
                0.5,
                "No valid cells",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        categories = sorted(frame[color_key].astype(str).unique().tolist())
        hier_palette = get_hierarchical_palette(color_key, categories)
        if hier_palette is None:
            cmap = plt.get_cmap("tab20", max(len(categories), 1))
            color_lookup = {
                category: cmap(idx) for idx, category in enumerate(categories)
            }
        else:
            color_lookup = hier_palette

        for category in categories:
            subset = frame.loc[frame[color_key].astype(str) == category]
            ax.scatter(
                subset["umap_1"],
                subset["umap_2"],
                s=max(
                    float(config.plotting.point_size) * 1.1,
                    float(config.plotting.point_size) + 1.0,
                ),
                c=[color_lookup[category]],
                label=category,
                linewidths=0,
                alpha=0.9,
            )
        ax.legend(
            title=color_key,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=config.plotting.legend_fontsize,
            title_fontsize=config.plotting.legend_title_fontsize,
            markerscale=3.0,
        )

    fig.suptitle(title)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.98))
    fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
    plt.close(fig)


__all__ = [
    "save_categorical_umap",
    "save_annotation_method_comparison_umap",
    "build_cima_l1_palette",
    "build_cima_l2_palette",
    "get_hierarchical_palette",
]
