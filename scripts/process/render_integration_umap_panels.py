#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import to_hex, to_rgb


BASE_L1_PALETTE = {
    "B": "#2C7BB6",
    "CD4_T": "#D7191C",
    "CD8_T&unconvensional_T": "#FDAE61",
    "Myeloid": "#1A9641",
    "ILC": "#762A83",
}
FALLBACK_COLORS = [
    "#5E4FA2",
    "#3288BD",
    "#66C2A5",
    "#ABDDA4",
    "#FEE08B",
    "#F46D43",
    "#A50026",
]
BACKGROUND_COLOR = "#D9D9D9"
REPO_ROOT = Path(__file__).resolve().parents[2]
OUTPUT_ROOT = REPO_ROOT / "output"
ANALYSIS_OUTPUT_ROOT = OUTPUT_ROOT / "1.only_atac"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render panel UMAP figures from an existing integration metadata CSV"
    )
    parser.add_argument(
        "--metadata-csv",
        default=str(
            ANALYSIS_OUTPUT_ROOT
            / "accepted_integration_bbknn_gsm"
            / "accepted_integration_metadata.csv"
        ),
    )
    parser.add_argument(
        "--reference-dir",
        default=str(OUTPUT_ROOT / "reference" / "cima"),
    )
    parser.add_argument(
        "--output-dir",
        default=str(ANALYSIS_OUTPUT_ROOT / "accepted_integration_bbknn_gsm"),
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=1.0,
        help="Point size for panel scatter plots",
    )
    parser.add_argument(
        "--ncols",
        type=int,
        default=4,
        help="Default panel column count for larger category sets",
    )
    return parser.parse_args()


def mix_color(color: str, target: str, amount: float) -> str:
    amount = min(max(amount, 0.0), 1.0)
    source_rgb = np.asarray(to_rgb(color), dtype=np.float32)
    target_rgb = np.asarray(to_rgb(target), dtype=np.float32)
    blended = source_rgb * (1 - amount) + target_rgb * amount
    return to_hex(
        (float(blended[0]), float(blended[1]), float(blended[2])), keep_alpha=False
    )


def make_shade_palette(base_color: str, labels: list[str]) -> dict[str, str]:
    labels = list(dict.fromkeys(labels))
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: base_color}
    colors = []
    for frac in np.linspace(0, 1, len(labels)):
        if frac <= 0.5:
            amount = 0.55 * (1 - frac / 0.5)
            colors.append(mix_color(base_color, "#FFFFFF", amount))
        else:
            amount = 0.18 * ((frac - 0.5) / 0.5)
            colors.append(mix_color(base_color, "#000000", amount))
    return dict(zip(labels, colors, strict=False))


def build_cima_palettes(hierarchy: pd.DataFrame) -> dict[str, dict[str, str]]:
    l1_labels = hierarchy["cell_type_l1"].drop_duplicates().tolist()
    l1_palette = dict(BASE_L1_PALETTE)
    missing = [label for label in l1_labels if label not in l1_palette]
    for label, color in zip(missing, FALLBACK_COLORS, strict=False):
        l1_palette[label] = color

    def level_palette(level_col: str) -> dict[str, str]:
        palette: dict[str, str] = {}
        for l1 in l1_labels:
            labels = sorted(
                hierarchy.loc[hierarchy["cell_type_l1"] == l1, level_col]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            palette.update(make_shade_palette(l1_palette[l1], labels))
        return palette

    return {
        "cima_cell_type_l1": {label: l1_palette[label] for label in l1_labels},
        "cima_cell_type_l2": level_palette("cell_type_l2"),
    }


def choose_layout(n_panels: int, default_ncols: int) -> tuple[int, int]:
    if n_panels <= 3:
        ncols = n_panels
    elif n_panels <= 6:
        ncols = 3
    elif n_panels <= 12:
        ncols = 4
    else:
        ncols = default_ncols
    nrows = math.ceil(n_panels / ncols)
    return nrows, ncols


def plot_panel_categorical(
    df: pd.DataFrame,
    x: str,
    y: str,
    color_col: str,
    palette: dict[str, str],
    out_path: Path,
    title: str,
    point_size: float,
    default_ncols: int,
) -> None:
    values = [
        label for label in palette if label in set(df[color_col].dropna().astype(str))
    ]
    if not values:
        raise ValueError(f"No categories found for {color_col}")

    nrows, ncols = choose_layout(len(values), default_ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(5.2 * ncols, 4.8 * nrows),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes_array = np.atleast_1d(axes).ravel()
    x_values = df[x].to_numpy(dtype=np.float32)
    y_values = df[y].to_numpy(dtype=np.float32)
    color_values = df[color_col].astype(str).to_numpy()

    for ax, label in zip(axes_array, values, strict=False):
        mask = color_values == label
        count = int(mask.sum())
        ax.scatter(
            x_values[~mask],
            y_values[~mask],
            s=point_size,
            alpha=0.12,
            c=BACKGROUND_COLOR,
            linewidths=0,
            rasterized=True,
        )
        ax.scatter(
            x_values[mask],
            y_values[mask],
            s=point_size,
            alpha=0.85,
            c=palette[label],
            linewidths=0,
            rasterized=True,
        )
        ax.set_title(f"{label} (n={count:,})", fontsize=11, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor("white")

    for ax in axes_array[len(values) :]:
        ax.axis("off")

    fig.suptitle(title, fontsize=18, fontweight="bold")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    metadata_csv = Path(args.metadata_csv)
    reference_dir = Path(args.reference_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(metadata_csv)
    hierarchy = pd.read_csv(reference_dir / "cima_atac_celltype_hierarchy.csv")
    palettes = build_cima_palettes(hierarchy)

    meta["health_status_plot"] = (
        meta["health_status"]
        .astype(str)
        .map(
            lambda value: {
                "健康": "Healthy",
                "COVID-19": "COVID-19",
                "RSV": "RSV",
                "未知": "Unknown",
            }.get(value, value)
        )
        .fillna(meta["health_status"].astype(str))
    )

    gse_values = sorted(meta["gse"].dropna().astype(str).unique().tolist())
    gse_palette = {
        g: c for g, c in zip(gse_values, ["#1f78b4", "#e31a1c"], strict=False)
    }
    health_values = sorted(
        meta["health_status_plot"].dropna().astype(str).unique().tolist()
    )
    health_palette_base = {
        "Healthy": "#2C7BB6",
        "COVID-19": "#D7191C",
        "RSV": "#1A9641",
        "Unknown": "#7F7F7F",
    }
    health_palette = {
        health: health_palette_base.get(health, "#7F7F7F") for health in health_values
    }

    plot_panel_categorical(
        meta,
        "integrated_umap_1",
        "integrated_umap_2",
        "cima_cell_type_l1",
        palettes["cima_cell_type_l1"],
        output_dir / "accepted_integration_cima_l1_panels.png",
        "Accepted-sample integration panels by CIMA L1",
        args.point_size,
        args.ncols,
    )
    plot_panel_categorical(
        meta,
        "integrated_umap_1",
        "integrated_umap_2",
        "cima_cell_type_l2",
        palettes["cima_cell_type_l2"],
        output_dir / "accepted_integration_cima_l2_panels.png",
        "Accepted-sample integration panels by CIMA L2",
        args.point_size,
        args.ncols,
    )
    plot_panel_categorical(
        meta,
        "integrated_umap_1",
        "integrated_umap_2",
        "gse",
        gse_palette,
        output_dir / "accepted_integration_gse_panels.png",
        "Accepted-sample integration panels by GSE",
        args.point_size,
        args.ncols,
    )
    plot_panel_categorical(
        meta,
        "integrated_umap_1",
        "integrated_umap_2",
        "health_status_plot",
        health_palette,
        output_dir / "accepted_integration_health_status_panels.png",
        "Accepted-sample integration panels by health status",
        args.point_size,
        args.ncols,
    )

    print(output_dir)


if __name__ == "__main__":
    main()
