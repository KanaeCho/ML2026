from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd

from .models import RunConfig


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


__all__ = ["save_categorical_umap"]
