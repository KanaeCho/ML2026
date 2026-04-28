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


BACKGROUND_COLOR = "#D9D9D9"
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Render old only_atac-style UMAP panel figures for product metadata."
    )
    parser.add_argument("--metadata-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--product", required=True)
    parser.add_argument("--point-size", type=float, default=1.0)
    parser.add_argument("--ncols", type=int, default=4)
    return parser.parse_args()


def mix_color(color: str, target: str, amount: float) -> str:
    amount = min(max(amount, 0.0), 1.0)
    source_rgb = np.asarray(to_rgb(color), dtype=np.float32)
    target_rgb = np.asarray(to_rgb(target), dtype=np.float32)
    blended = source_rgb * (1 - amount) + target_rgb * amount
    return to_hex(tuple(float(v) for v in blended), keep_alpha=False)


def make_shade_palette(base_color: str, labels: list[str]) -> dict[str, str]:
    labels = list(dict.fromkeys(labels))
    if not labels:
        return {}
    if len(labels) == 1:
        return {labels[0]: base_color}

    colors = []
    for frac in np.linspace(0, 1, len(labels)):
        if frac <= 0.5:
            colors.append(mix_color(base_color, "#FFFFFF", 0.55 * (1 - frac / 0.5)))
        else:
            colors.append(mix_color(base_color, "#000000", 0.18 * ((frac - 0.5) / 0.5)))
    return dict(zip(labels, colors, strict=False))


def default_palette(labels: list[str]) -> dict[str, str]:
    labels = [str(label) for label in labels if str(label) and str(label) != "nan"]
    palette = {label: BASE_L1_PALETTE[label] for label in labels if label in BASE_L1_PALETTE}
    missing = [label for label in labels if label not in palette]
    if not missing:
        return palette
    cmap = plt.get_cmap("tab20", max(len(missing), 1))
    for i, label in enumerate(missing):
        fallback = FALLBACK_COLORS[i % len(FALLBACK_COLORS)] if len(missing) <= len(FALLBACK_COLORS) else to_hex(cmap(i))
        palette[label] = fallback
    return palette


def l2_palette_from_l1(df: pd.DataFrame, l1_col: str, l2_col: str) -> dict[str, str]:
    if l1_col not in df.columns or l2_col not in df.columns:
        return {}
    l1_palette = default_palette(sorted(df[l1_col].dropna().astype(str).unique().tolist()))
    palette: dict[str, str] = {}
    for l1 in sorted(df[l1_col].dropna().astype(str).unique().tolist()):
        labels = sorted(
            df.loc[df[l1_col].astype(str) == l1, l2_col]
            .dropna()
            .astype(str)
            .unique()
            .tolist()
        )
        palette.update(make_shade_palette(l1_palette.get(l1, "#7F7F7F"), labels))
    return palette


def choose_layout(n_panels: int, default_ncols: int) -> tuple[int, int]:
    if n_panels <= 3:
        ncols = max(n_panels, 1)
    elif n_panels <= 6:
        ncols = 3
    elif n_panels <= 12:
        ncols = 4
    else:
        ncols = default_ncols
    return math.ceil(n_panels / ncols), ncols


def detect_umap_columns(df: pd.DataFrame) -> tuple[str, str] | None:
    candidates = [
        ("integrated_umap_1", "integrated_umap_2"),
        ("cima_ref_umap_1", "cima_ref_umap_2"),
        ("umap_atac_1", "umap_atac_2"),
        ("umap_1", "umap_2"),
    ]
    for x_col, y_col in candidates:
        if x_col in df.columns and y_col in df.columns:
            return x_col, y_col
    return None


def first_existing(df: pd.DataFrame, columns: list[str]) -> str | None:
    for column in columns:
        if column in df.columns and df[column].notna().any():
            return column
    return None


def plot_panel_categorical(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    color_col: str,
    palette: dict[str, str],
    out_path: Path,
    title: str,
    point_size: float,
    default_ncols: int,
) -> None:
    plot_df = df[[x_col, y_col, color_col]].copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df = plot_df.dropna(subset=[x_col, y_col, color_col])
    plot_df[color_col] = plot_df[color_col].astype(str)
    values = [label for label in palette if label in set(plot_df[color_col])]
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
    x_values = plot_df[x_col].to_numpy(dtype=np.float32)
    y_values = plot_df[y_col].to_numpy(dtype=np.float32)
    color_values = plot_df[color_col].to_numpy()

    for ax, label in zip(axes_array, values, strict=False):
        mask = color_values == label
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
        ax.set_title(f"{label} (n={int(mask.sum()):,})", fontsize=11, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor("white")

    for ax in axes_array[len(values) :]:
        ax.axis("off")

    fig.suptitle(title, fontsize=18, fontweight="bold")
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def render_panels(
    metadata_csv: Path,
    output_dir: Path,
    product: str,
    point_size: float = 1.0,
    ncols: int = 4,
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(metadata_csv, low_memory=False)
    umap_cols = detect_umap_columns(df)
    if umap_cols is None:
        return {"status": "skipped", "detail": "missing_umap_columns"}

    x_col, y_col = umap_cols
    rendered: dict[str, str] = {"status": "ok", "x_col": x_col, "y_col": y_col}
    l1_col = first_existing(
        df,
        ["integrated_cima_l1", "cima_cell_type_l1", "azimuth_cima_l1", "cima_l1_masked", "cima_l1"],
    )
    l2_col = first_existing(
        df,
        ["integrated_cima_l2", "cima_cell_type_l2", "cima_l2", "azimuth_cell_type_l2_raw"],
    )

    plot_specs: list[tuple[str, str, dict[str, str], str]] = []
    if l1_col:
        labels = sorted(df[l1_col].dropna().astype(str).unique().tolist())
        plot_specs.append((l1_col, "cima_l1", default_palette(labels), f"{product} panels by CIMA L1"))
    if l2_col:
        palette = l2_palette_from_l1(df, l1_col, l2_col) if l1_col else {}
        if not palette:
            palette = default_palette(sorted(df[l2_col].dropna().astype(str).unique().tolist()))
        plot_specs.append((l2_col, "cima_l2", palette, f"{product} panels by CIMA L2"))
    for column, suffix, title in [
        ("integrated_cluster", "integrated_cluster", f"{product} panels by integrated cluster"),
        ("gse", "gse", f"{product} panels by GSE"),
        ("sample_id", "sample", f"{product} panels by sample"),
    ]:
        if column in df.columns:
            labels = sorted(df[column].dropna().astype(str).unique().tolist())
            plot_specs.append((column, suffix, default_palette(labels), title))

    rendered["coordinate_source"] = "integrated_umap" if (x_col, y_col) == ("integrated_umap_1", "integrated_umap_2") else "fallback"
    for color_col, suffix, palette, title in plot_specs:
        out_path = output_dir / f"{product}_{suffix}_panels.png"
        try:
            plot_panel_categorical(df, x_col, y_col, color_col, palette, out_path, title, point_size, ncols)
            rendered[f"figure:{suffix}"] = str(out_path)
            rendered[f"column:{suffix}"] = color_col
        except Exception as exc:  # pragma: no cover - kept as product-level audit detail
            rendered[f"figure:{suffix}"] = f"skipped:{exc}"
    return rendered


def main() -> None:
    args = parse_args()
    result = render_panels(
        metadata_csv=Path(args.metadata_csv),
        output_dir=Path(args.output_dir),
        product=args.product,
        point_size=args.point_size,
        ncols=args.ncols,
    )
    print(result)


if __name__ == "__main__":
    main()
