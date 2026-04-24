from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex, to_rgb
from typing import Any, cast
from matplotlib.lines import Line2D
from scipy.stats import gaussian_kde

from .label_alignment import align_pbmcref_to_cima_l1
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


def _make_family_variant_palette(base_color: str, labels: list[str]) -> dict[str, str]:
    labels = list(dict.fromkeys(labels))
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: base_color}

    variants: list[str] = []
    offsets = np.linspace(-0.12, 0.12, len(labels))
    for offset in offsets:
        if offset < 0:
            variants.append(_mix_color(base_color, "#FFFFFF", abs(float(offset))))
        elif offset > 0:
            variants.append(_mix_color(base_color, "#000000", float(offset)))
        else:
            variants.append(base_color)
    return dict(zip(labels, variants, strict=False))


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
    if color_key in {
        "azimuth_cell_type",
        "azimuth_cell_type_l1_raw",
        "azimuth_cell_type_l2_raw",
    }:
        raw_series = pd.Series(
            categories, index=pd.Index(categories, dtype=object), dtype="string"
        )
        raw_l1 = pd.Series(pd.NA, index=raw_series.index, dtype="string")
        for category in categories:
            stripped = str(category).strip()
            if stripped in {"B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other T"}:
                raw_l1.loc[category] = stripped
        raw_l2 = raw_series if color_key != "azimuth_cell_type_l1_raw" else None
        aligned, _ = align_pbmcref_to_cima_l1(raw_l1, raw_l2)
        grouped: dict[str, list[str]] = {}
        for category in categories:
            grouped.setdefault(str(aligned.loc[category]), []).append(category)

        family_palette = build_cima_l1_palette(list(grouped.keys()))
        palette: dict[str, tuple[float, ...]] = {}
        for family, family_labels in grouped.items():
            sorted_labels = sorted(family_labels)
            if len(sorted_labels) == 1:
                shades = {sorted_labels[0]: family_palette[family]}
            else:
                shades = _make_family_variant_palette(
                    family_palette[family], sorted_labels
                )
            for label, hex_color in shades.items():
                palette[label] = to_rgb(hex_color)
        return palette
    if color_key in (
        "cima_l1",
        "cima_l1_masked",
        "azimuth_cima_l1",
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
    obs = cast(pd.DataFrame, adata.obs)
    if color_key not in obs.columns:
        raise KeyError(f"adata.obs must contain '{color_key}'")

    if "umap_1" not in obs.columns or "umap_2" not in obs.columns:
        return pd.DataFrame(columns=["umap_1", "umap_2", color_key])

    frame = pd.DataFrame(
        {
            "umap_1": pd.to_numeric(obs["umap_1"], errors="coerce"),
            "umap_2": pd.to_numeric(obs["umap_2"], errors="coerce"),
            color_key: obs[color_key].astype("string"),
        },
        index=adata.obs_names,
    )
    return frame.dropna(subset=["umap_1", "umap_2", color_key])


def _add_category_text_labels(
    ax: plt.Axes,
    frame: pd.DataFrame,
    color_key: str,
) -> None:
    if frame.empty:
        return

    def _family_for_label(label: str) -> str:
        raw_series = pd.Series([label], dtype="string")
        raw_l1 = pd.Series([pd.NA], dtype="string")
        stripped = label.strip()
        if stripped in {"B", "CD4 T", "CD8 T", "DC", "Mono", "NK", "other T"}:
            raw_l1.iloc[0] = stripped
        aligned, _ = align_pbmcref_to_cima_l1(raw_l1, raw_series)
        return str(aligned.iloc[0])

    def _component_groups(points: np.ndarray, radius: float) -> list[np.ndarray]:
        n_points = len(points)
        if n_points == 0:
            return []
        visited = np.zeros(n_points, dtype=bool)
        groups: list[np.ndarray] = []
        for start in range(n_points):
            if visited[start]:
                continue
            stack = [start]
            visited[start] = True
            component: list[int] = []
            while stack:
                idx = stack.pop()
                component.append(idx)
                deltas = points - points[idx]
                distances = np.sqrt((deltas**2).sum(axis=1))
                neighbors = np.where((distances <= radius) & (~visited))[0]
                for neighbor in neighbors:
                    visited[neighbor] = True
                    stack.append(int(neighbor))
            groups.append(np.array(component, dtype=int))
        return groups

    x_span = (
        float(frame["umap_1"].max() - frame["umap_1"].min()) if not frame.empty else 0.0
    )
    y_span = (
        float(frame["umap_2"].max() - frame["umap_2"].min()) if not frame.empty else 0.0
    )
    cluster_radius = max(max(x_span, y_span) * 0.08, 0.5)

    for label, subset in frame.groupby(color_key, dropna=True):
        label = str(label)
        if label.strip().lower() in {"unknown", "na", "<na>"}:
            continue
        if _family_for_label(label) == "Unknown":
            continue

        points = subset[["umap_1", "umap_2"]].to_numpy(dtype=float)
        components = _component_groups(points, cluster_radius)
        min_component_size = max(2, int(np.ceil(len(points) * 0.40)))
        category_colors = get_hierarchical_palette(color_key, [label]) or {}
        label_facecolor = category_colors.get(label, to_rgb("#FFFFFF"))
        for component in components:
            if len(component) < min_component_size:
                continue
            component_points = points[component]
            center = np.median(component_points, axis=0)
            ax.text(
                float(center[0]),
                float(center[1]),
                label,
                ha="center",
                va="center",
                fontsize=8,
                fontweight="bold",
                color="black",
                bbox={
                    "boxstyle": "round,pad=0.18",
                    "facecolor": label_facecolor,
                    "edgecolor": "none",
                    "alpha": 0.72,
                },
                zorder=10,
            )


def save_highlight_category_overview(
    *,
    adata: ad.AnnData,
    color_key: str,
    output_path: Path,
    title: str,
    config: RunConfig,
    legend_title: str,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frame = _plot_frame(adata, color_key)
    if frame.empty:
        fig, ax = plt.subplots(
            figsize=(
                float(config.plotting.umap_width),
                float(config.plotting.umap_height),
            ),
            dpi=config.plotting.dpi,
        )
        ax.set_title(title)
        ax.text(
            0.5, 0.5, "No valid cells", ha="center", va="center", transform=ax.transAxes
        )
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
        fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
        plt.close(fig)
        return

    categories = sorted(frame[color_key].astype(str).unique().tolist())
    n_panels = len(categories)
    ncols = min(4, n_panels) if n_panels > 1 else 1
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(
            float(config.plotting.umap_width) * ncols,
            float(config.plotting.umap_height) * nrows,
        ),
        dpi=config.plotting.dpi,
    )
    axes_array = np.atleast_1d(axes).ravel()

    color_lookup = get_hierarchical_palette(color_key, categories) or {
        category: plt.get_cmap("tab20", max(len(categories), 1))(idx)
        for idx, category in enumerate(categories)
    }

    background = frame.copy()
    highlight_point_size = max(float(config.plotting.point_size) * 0.72, 1.6)
    background_point_size = max(float(config.plotting.point_size) * 0.18, 0.7)

    def _draw_smooth_component_outline(
        ax: plt.Axes, points: np.ndarray, color: tuple[float, ...]
    ) -> None:
        if len(points) < 8:
            return

        x = points[:, 0]
        y = points[:, 1]
        x_pad = max((float(x.max()) - float(x.min())) * 0.22, 0.45)
        y_pad = max((float(y.max()) - float(y.min())) * 0.22, 0.45)
        grid_x, grid_y = np.meshgrid(
            np.linspace(float(x.min()) - x_pad, float(x.max()) + x_pad, 180),
            np.linspace(float(y.min()) - y_pad, float(y.max()) + y_pad, 180),
        )

        try:
            kde = gaussian_kde(points.T, bw_method=0.35)
            density = kde(np.vstack([grid_x.ravel(), grid_y.ravel()])).reshape(grid_x.shape)
        except np.linalg.LinAlgError:
            return
        except ValueError:
            return

        max_density = float(np.nanmax(density))
        if not np.isfinite(max_density) or max_density <= 0.0:
            return

        contour_level = max_density * 0.24
        ax.contour(
            grid_x,
            grid_y,
            density,
            levels=[contour_level],
            colors=[color],
            linewidths=1.6,
            linestyles=[(0, (4.0, 2.4))],
            alpha=0.98,
            zorder=5,
        )

    def _component_groups(points: np.ndarray, radius: float) -> list[np.ndarray]:
        n_points = len(points)
        if n_points == 0:
            return []
        visited = np.zeros(n_points, dtype=bool)
        groups: list[np.ndarray] = []
        for start in range(n_points):
            if visited[start]:
                continue
            stack = [start]
            visited[start] = True
            component: list[int] = []
            while stack:
                idx = stack.pop()
                component.append(idx)
                deltas = points - points[idx]
                distances = np.sqrt((deltas**2).sum(axis=1))
                neighbors = np.where((distances <= radius) & (~visited))[0]
                for neighbor in neighbors:
                    visited[neighbor] = True
                    stack.append(int(neighbor))
            groups.append(np.array(component, dtype=int))
        return groups

    x_span = float(background["umap_1"].max() - background["umap_1"].min())
    y_span = float(background["umap_2"].max() - background["umap_2"].min())
    outline_radius = max(max(x_span, y_span) * 0.08, 0.5)

    for panel_idx, (ax, category) in enumerate(zip(axes_array, categories, strict=False)):
        ax.set_title(
            category,
            fontsize=max(config.plotting.legend_title_fontsize + 5.0, 14.0),
            fontweight="bold",
            pad=4.0,
        )
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_box_aspect(1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=max(config.plotting.legend_fontsize - 0.5, 6.0))

        ax.scatter(
            background["umap_1"],
            background["umap_2"],
            s=background_point_size,
            c=[to_rgb("#B5B5B5")],
            linewidths=0,
            alpha=0.34,
        )

        subset = frame.loc[frame[color_key].astype(str) == category]
        subset_points = subset[["umap_1", "umap_2"]].to_numpy(dtype=float)
        ax.scatter(
            subset["umap_1"],
            subset["umap_2"],
            s=highlight_point_size,
            c=[color_lookup[category]],
            linewidths=0,
            alpha=0.98,
            label=category,
        )

        for component in _component_groups(subset_points, outline_radius):
            if len(component) < max(8, int(np.ceil(len(subset_points) * 0.12))):
                continue
            _draw_smooth_component_outline(ax, subset_points[component], color_lookup[category])

        legend_handle = Line2D(
            [0],
            [0],
            marker="o",
            linestyle=(0, (4.0, 2.4)),
            color=tuple(color_lookup[category][:3]),
            markerfacecolor=tuple(color_lookup[category][:3]),
            markeredgewidth=0.0,
            markersize=5.5,
            linewidth=1.4,
        )
        ax.legend(
            handles=[legend_handle],
            labels=[category],
            title=None,
            loc="lower right",
            frameon=False,
            fontsize=max(config.plotting.legend_fontsize + 0.2, 6.8),
            handletextpad=0.35,
            borderpad=0.15,
        )

        if panel_idx % ncols != 0:
            ax.set_ylabel("")
        if panel_idx // ncols != nrows - 1:
            ax.set_xlabel("")

    for ax in axes_array[n_panels:]:
        ax.axis("off")

    fig.suptitle(title, fontsize=max(config.plotting.legend_title_fontsize + 0.5, 9.0))
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
    plt.close(fig)


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
            legend_box = getattr(legend, "_legend_box", None)
            if legend_box is not None:
                legend_box.align = "left"

        ax.set_title(title)
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_aspect("auto")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0))
    fig.savefig(output_path, dpi=plotting.dpi, format="png")
    plt.close(fig)


def save_sample_cima_l1_umap(
    adata: ad.AnnData,
    output_path: Path,
    title: str,
    config: RunConfig,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cima_l1_key = "azimuth_cima_l1" if "azimuth_cima_l1" in adata.obs.columns else ""
    frame = _plot_frame(adata, cima_l1_key) if cima_l1_key else pd.DataFrame()
    if not frame.empty:
        frame = frame.loc[
            ~frame[cima_l1_key].astype(str).str.strip().str.lower().eq("unknown")
        ]

    plotting = config.plotting
    fig, ax = plt.subplots(
        figsize=(float(plotting.umap_width), float(plotting.umap_height)),
        dpi=plotting.dpi,
    )
    ax.set_box_aspect(1)

    if frame.empty:
        ax.set_title(title, fontsize=max(plotting.legend_title_fontsize + 2.0, 12.0), fontweight="bold")
        ax.text(0.5, 0.5, "No labeled cells", ha="center", va="center", transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
    else:
        categories = sorted(frame[cima_l1_key].astype(str).unique().tolist())
        palette = get_hierarchical_palette("cima_l1", categories) or {
            category: plt.get_cmap("tab20", max(len(categories), 1))(idx)
            for idx, category in enumerate(categories)
        }
        point_size = max(float(plotting.point_size) * 1.15, float(plotting.point_size) + 1.0)

        for category in categories:
            subset = frame.loc[frame[cima_l1_key].astype(str) == category]
            ax.scatter(
                subset["umap_1"],
                subset["umap_2"],
                s=point_size,
                c=[palette[category]],
                label=category,
                linewidths=0,
                alpha=0.92,
            )

        ax.legend(
            title="CIMA L1",
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=plotting.legend_fontsize,
            title_fontsize=plotting.legend_title_fontsize,
            markerscale=2.8,
        )
        ax.set_title(title, fontsize=max(plotting.legend_title_fontsize + 2.0, 12.0), fontweight="bold")
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.tick_params(labelsize=max(plotting.legend_fontsize - 1.0, 7.0))

    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.1)
        spine.set_edgecolor("#4A4A4A")

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0))
    fig.savefig(output_path, dpi=plotting.dpi, format="png")
    plt.close(fig)


def save_qc_overview(
    adata: ad.AnnData,
    output_path: Path,
    config: RunConfig,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    qc_thresholds = dict(cast(dict[str, Any], adata.uns.get("qc_thresholds", {})))
    panel_specs = [
        (
            "n_counts",
            "fails_count_floor",
            float(qc_thresholds.get("min_counts", config.qc.count_floor_min)),
            cast(dict[str, Any], qc_thresholds.get("n_counts_audit", {})),
        ),
        (
            "n_genes",
            "fails_gene_floor",
            float(qc_thresholds.get("min_genes", config.qc.gene_floor_min)),
            cast(dict[str, Any], qc_thresholds.get("n_genes_audit", {})),
        ),
        (
            "pct_mt",
            "fails_mt_ceiling",
            float(qc_thresholds.get("max_pct_mt", config.qc.pct_mt_ceiling_max)),
            cast(dict[str, Any], qc_thresholds.get("pct_mt_audit", {})),
        ),
        (
            "pct_ribo",
            "fails_ribo_ceiling",
            float(qc_thresholds.get("max_pct_ribo", config.qc.pct_ribo_ceiling_max)),
            cast(dict[str, Any], qc_thresholds.get("pct_ribo_audit", {})),
        ),
    ]

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(float(config.plotting.umap_width), float(config.plotting.umap_height)),
        dpi=config.plotting.dpi,
    )
    axes_array = np.atleast_1d(axes).ravel()

    obs = cast(pd.DataFrame, adata.obs)

    for ax, (metric, fail_flag, threshold, audit) in zip(axes_array, panel_specs, strict=True):
        raw_metric = obs[metric] if metric in obs.columns else None
        if raw_metric is None:
            values = pd.Series([], dtype=float)
        else:
            values = cast(pd.Series, pd.to_numeric(raw_metric, errors="coerce")).dropna()

        transform = str(audit.get("transform", "identity"))
        direction = str(audit.get("direction", "upper"))
        center = float(audit.get("center", np.nan))
        mad = float(audit.get("mad", np.nan))
        nmads = float(audit.get("nmads", np.nan))
        raw_threshold = float(audit.get("raw_threshold", np.nan))
        final_threshold = float(audit.get("final_threshold_original_scale", threshold))
        guardrails = [str(item) for item in audit.get("guardrails_applied", [])]

        if transform == "log10p1":
            display_values = np.log10(values.to_numpy(dtype=float) + 1.0)
            display_center = center
            display_raw_threshold = raw_threshold
            display_final_threshold = np.log10(max(final_threshold, 0.0) + 1.0)
            axis_label = f"{metric} [log10(x + 1)]"
            method_label = "lower-tail MAD in log10(x + 1)"
        else:
            display_values = values.to_numpy(dtype=float)
            display_center = center
            display_raw_threshold = raw_threshold
            display_final_threshold = final_threshold
            axis_label = metric
            method_label = "upper-tail MAD in original scale"

        raw_fail_series = obs[fail_flag] if fail_flag in obs.columns else None
        if raw_fail_series is None:
            fail_series = pd.Series(False, index=obs.index, dtype=bool)
        else:
            fail_series = cast(pd.Series, raw_fail_series.fillna(False).astype(bool))
        fail_fraction = float(fail_series.mean()) if len(fail_series) else 0.0

        if display_values.size > 0:
            ax.hist(
                display_values,
                bins=36,
                color="#4C72B0",
                alpha=0.88,
                edgecolor="white",
                linewidth=0.3,
            )
        if np.isfinite(display_center):
            ax.axvline(
                display_center,
                color="#1F77B4",
                linestyle="-",
                linewidth=1.2,
            )
        if np.isfinite(display_raw_threshold):
            ax.axvline(
                display_raw_threshold,
                color="#FF7F0E",
                linestyle=":",
                linewidth=1.4,
            )
        if np.isfinite(display_final_threshold):
            ax.axvline(
                display_final_threshold,
                color="#D62728",
                linestyle="--",
                linewidth=1.6,
            )

        ax.set_title(f"{metric} | {method_label} | fail {fail_fraction:.1%}")
        ax.set_xlabel(axis_label)
        ax.set_ylabel("cells")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=max(config.plotting.legend_fontsize - 1.0, 7.0))
        ax.text(
            0.02,
            0.98,
            (
                f"median={center:.3f} | MAD={mad:.3f} | nMAD={nmads:.2f}\n"
                f"raw_thr={raw_threshold:.3f} | final_thr={final_threshold:.2f}"
            ),
            ha="left",
            va="top",
            transform=ax.transAxes,
            fontsize=max(config.plotting.legend_fontsize - 1.7, 6.2),
            color="#444444",
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "edgecolor": "#D8D8D8",
                "alpha": 0.88,
            },
        )
        ax.text(
            0.98,
            0.98,
            "blue=median\norange=raw MAD cutoff\nred=final applied cutoff",
            ha="right",
            va="top",
            transform=ax.transAxes,
            fontsize=max(config.plotting.legend_fontsize - 1.9, 6.0),
            color="#444444",
        )
        if guardrails:
            ax.text(
                0.98,
                0.08,
                f"guardrail: {', '.join(guardrails)}",
                ha="right",
                va="bottom",
                transform=ax.transAxes,
                fontsize=max(config.plotting.legend_fontsize - 1.9, 6.0),
                color="#8C2D04",
            )

    fig.tight_layout(rect=(0.0, 0.0, 1.0, 1.0))
    fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
    plt.close(fig)


def save_dual_annotation_umap(
    *,
    adata: ad.AnnData,
    output_path: Path,
    title: str,
    config: RunConfig,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    raw_color_key = (
        "azimuth_cell_type_l2_raw"
        if "azimuth_cell_type_l2_raw" in adata.obs.columns
        else "azimuth_cell_type"
    )
    panel_specs = [
        (raw_color_key, "pbmcref"),
        ("azimuth_cima_l1", "CIMA L1"),
    ]

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(
            float(config.plotting.umap_width) * 2.0,
            float(config.plotting.umap_height),
        ),
        dpi=config.plotting.dpi,
    )
    axes_array = np.atleast_1d(axes).ravel()
    display_point_size = max(
        float(config.plotting.point_size) * 1.1,
        float(config.plotting.point_size) + 1.0,
    )

    for ax, (color_key, panel_title) in zip(axes_array, panel_specs, strict=True):
        frame = (
            _plot_frame(adata, color_key)
            if color_key in adata.obs.columns
            else pd.DataFrame()
        )
        ax.set_title(panel_title)
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_box_aspect(1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=max(config.plotting.legend_fontsize - 0.5, 7.0))

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
        color_lookup = hier_palette or {
            category: plt.get_cmap("tab20", max(len(categories), 1))(idx)
            for idx, category in enumerate(categories)
        }

        for category in categories:
            subset = frame.loc[frame[color_key].astype(str) == category]
            ax.scatter(
                subset["umap_1"],
                subset["umap_2"],
                s=display_point_size,
                c=[color_lookup[category]],
                label=category,
                linewidths=0,
                alpha=0.9,
            )

        if panel_title == "pbmcref":
            _add_category_text_labels(ax, frame, color_key)

        ax.legend(
            title=panel_title,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=config.plotting.legend_fontsize,
            title_fontsize=config.plotting.legend_title_fontsize,
            markerscale=2.4,
            ncol=2 if panel_title == "pbmcref" and len(categories) >= 16 else 1,
        )

    fig.suptitle(title, fontsize=max(config.plotting.legend_title_fontsize + 0.5, 10.0))
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
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
    "save_sample_cima_l1_umap",
    "save_dual_annotation_umap",
    "save_qc_overview",
    "save_annotation_method_comparison_umap",
    "save_azimuth_candidate_overview",
    "save_highlight_category_overview",
    "build_cima_l1_palette",
    "build_cima_l2_palette",
    "get_hierarchical_palette",
]


def save_azimuth_candidate_overview(
    *,
    candidates: list[tuple[str, float, ad.AnnData]],
    output_path: Path,
    title: str,
    config: RunConfig,
    color_key: str | None = None,
    legend_title: str | None = None,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if not candidates:
        fig, ax = plt.subplots(
            figsize=(
                float(config.plotting.umap_width),
                float(config.plotting.umap_height),
            ),
            dpi=config.plotting.dpi,
        )
        ax.set_title(title)
        ax.text(
            0.5,
            0.5,
            "No candidate UMAPs",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        ax.set_xticks([])
        ax.set_yticks([])
        fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.95))
        fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
        plt.close(fig)
        return

    n_panels = len(candidates)
    ncols = min(3, n_panels) if n_panels > 1 else 1
    nrows = int(np.ceil(n_panels / ncols))
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(
            float(config.plotting.umap_width) * ncols,
            float(config.plotting.umap_height) * nrows,
        ),
        dpi=config.plotting.dpi,
    )
    axes_array = np.atleast_1d(axes).ravel()
    overview_point_size = max(
        float(config.plotting.point_size) * 0.30,
        0.8,
    )

    resolved_color_key = color_key or (
        "azimuth_cima_l1"
        if any("azimuth_cima_l1" in adata.obs.columns for _, _, adata in candidates)
        else "azimuth_cell_type"
    )
    resolved_legend_title = legend_title or (
        "CIMA L1" if resolved_color_key == "azimuth_cima_l1" else "azimuth_cell_type"
    )

    for ax, (candidate_id, total_score, adata) in zip(
        axes_array, candidates, strict=False
    ):
        frame = (
            _plot_frame(adata, resolved_color_key)
            if resolved_color_key in adata.obs.columns
            else pd.DataFrame()
        )
        ax.set_title(f"{candidate_id}\nscore={total_score:.3f}")
        ax.set_xlabel("umap_1")
        ax.set_ylabel("umap_2")
        ax.set_box_aspect(1)

        if frame.empty:
            ax.text(
                0.5,
                0.5,
                "No valid Azimuth UMAP",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            ax.set_xticks([])
            ax.set_yticks([])
            continue

        categories = sorted(frame[resolved_color_key].astype(str).unique().tolist())
        hier_palette = get_hierarchical_palette(resolved_color_key, categories)
        color_lookup = hier_palette or {
            category: plt.get_cmap("tab20", max(len(categories), 1))(idx)
            for idx, category in enumerate(categories)
        }

        for category in categories:
            subset = frame.loc[frame[resolved_color_key].astype(str) == category]
            ax.scatter(
                subset["umap_1"],
                subset["umap_2"],
                s=overview_point_size,
                c=[color_lookup[category]],
                label=category,
                linewidths=0,
                alpha=0.9,
            )

        ax.legend(
            title=resolved_legend_title,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            frameon=False,
            fontsize=config.plotting.legend_fontsize,
            title_fontsize=config.plotting.legend_title_fontsize,
            markerscale=3.0,
        )

    for ax in axes_array[n_panels:]:
        ax.axis("off")

    fig.suptitle(title)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.96))
    fig.savefig(output_path, dpi=config.plotting.dpi, format="png")
    plt.close(fig)
