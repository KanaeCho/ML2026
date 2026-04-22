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
    max_candidates: int = 1


@dataclass(frozen=True)
class QcThresholds:
    method: str = "dynamic_hybrid_mad"
    counts_lower_nmads: float = 3.0
    genes_lower_nmads: float = 3.0
    pct_mt_upper_nmads: float = 3.0
    pct_ribo_upper_nmads: float = 3.5
    min_cells_for_dynamic: int = 50
    count_floor_min: int = 100
    count_floor_max: int = 1500
    gene_floor_min: int = 100
    gene_floor_max: int = 1200
    pct_mt_ceiling_min: float = 5.0
    pct_mt_ceiling_max: float = 40.0
    pct_ribo_ceiling_min: float = 20.0
    pct_ribo_ceiling_max: float = 80.0


@dataclass(frozen=True)
class QcMetricThresholdAudit:
    transform: str
    direction: str
    center: float
    mad: float
    nmads: float
    raw_threshold: float
    final_threshold_original_scale: float
    n_cells_used: int
    guardrails_applied: tuple[str, ...] = ()


@dataclass(frozen=True)
class ComputedQcThresholds:
    sample_id: str
    gse: str
    method: str
    n_cells_total: int
    min_counts: int
    min_genes: int
    max_pct_mt: float
    max_pct_ribo: float
    n_counts_audit: QcMetricThresholdAudit
    n_genes_audit: QcMetricThresholdAudit
    pct_mt_audit: QcMetricThresholdAudit
    pct_ribo_audit: QcMetricThresholdAudit
    small_sample_rule_used: bool = False
    zero_mad_metrics: tuple[str, ...] = ()
    retention_guard_triggered: bool = False


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
