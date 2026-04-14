from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QcThresholds:
    # Minimal quality control thresholds
    min_counts: int
    min_genes: int
    max_pct_mt: float
    max_pct_ribo: float


@dataclass(frozen=True)
class PlottingConfig:
    # Basic plotting configuration for UMAP visualization
    umap_width: float
    umap_height: float
    dpi: int
    point_size: float
    legend_fontsize: float
    legend_title_fontsize: float


@dataclass(frozen=True)
class RunConfig:
    qc: QcThresholds
    plotting: PlottingConfig
