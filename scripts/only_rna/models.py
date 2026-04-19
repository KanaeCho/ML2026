from __future__ import annotations

from dataclasses import dataclass
from typing import List


@dataclass(frozen=True)
class AnnotationConfig:
    # List of annotation backends to run. E.g. ["cima", "azimuth", "cell_typist"]
    methods: List[str]


@dataclass(frozen=True)
class EmbeddingConfig:
    n_top_genes: int = 2000
    n_pcs: int = 20
    n_neighbors: int = 15
    resolution: float = 1.0
    min_dist: float = 0.5
    spread: float = 1.0
    random_state: int = 0


@dataclass(frozen=True)
class AzimuthConfig:
    enabled: bool = False
    reference: str = "pbmcref"
    annotation_levels: tuple[str, ...] = ("l1", "l2")
    k_weight: int = 50
    n_trees: int = 20
    mapping_score_k: int = 100


@dataclass(frozen=True)
class TuningConfig:
    qc_preset_family: str = "default"
    azimuth_preset_family: str = "default"
    embedding_preset_family: str = "default"
    max_candidates: int = 9


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
    embedding: EmbeddingConfig = EmbeddingConfig()
    azimuth: AzimuthConfig = AzimuthConfig()
    tuning: TuningConfig = TuningConfig()
    annotation: AnnotationConfig | None = None
