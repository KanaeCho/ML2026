from __future__ import annotations

from pathlib import Path
from dataclasses import asdict, is_dataclass
from typing import Any, Iterable, cast
import yaml

from .models import (
    AnnotationConfig,
    AzimuthConfig,
    EmbeddingConfig,
    PlottingConfig,
    QcThresholds,
    RunConfig,
    TuningConfig,
)


def load_run_config(path: Path) -> RunConfig:
    """Load the YAML default RunConfig from the given path.

    YAML structure is expected to have two top-level keys:
    - qc: { min_counts, min_genes, max_pct_mt, max_pct_ribo }
    - plotting: { umap_width, umap_height, dpi, point_size, legend_fontsize, legend_title_fontsize }
    """
    path = Path(path)
    data = yaml.safe_load(path.read_text()) or {}
    qc_map = data.get("qc", {})
    plot_map = data.get("plotting", {})
    embedding_map = data.get("embedding", {})
    azimuth_map = data.get("azimuth", {})
    tuning_map = data.get("tuning", {})

    qc = QcThresholds(
        method=str(qc_map.get("method", "dynamic_hybrid_mad")),
        counts_lower_nmads=float(qc_map.get("counts_lower_nmads", 3.0)),
        genes_lower_nmads=float(qc_map.get("genes_lower_nmads", 3.0)),
        pct_mt_upper_nmads=float(qc_map.get("pct_mt_upper_nmads", 3.0)),
        pct_ribo_upper_nmads=float(qc_map.get("pct_ribo_upper_nmads", 3.5)),
        min_cells_for_dynamic=int(qc_map.get("min_cells_for_dynamic", 50)),
        count_floor_min=int(qc_map.get("count_floor_min", 100)),
        count_floor_max=int(qc_map.get("count_floor_max", 1500)),
        gene_floor_min=int(qc_map.get("gene_floor_min", 100)),
        gene_floor_max=int(qc_map.get("gene_floor_max", 1200)),
        pct_mt_ceiling_min=float(qc_map.get("pct_mt_ceiling_min", 5.0)),
        pct_mt_ceiling_max=float(qc_map.get("pct_mt_ceiling_max", 40.0)),
        pct_ribo_ceiling_min=float(qc_map.get("pct_ribo_ceiling_min", 20.0)),
        pct_ribo_ceiling_max=float(qc_map.get("pct_ribo_ceiling_max", 80.0)),
    )

    plotting = PlottingConfig(
        umap_width=float(plot_map.get("umap_width", 0.0)),
        umap_height=float(plot_map.get("umap_height", 0.0)),
        dpi=int(plot_map.get("dpi", 0)),
        point_size=float(plot_map.get("point_size", 0.0)),
        legend_fontsize=float(plot_map.get("legend_fontsize", 0.0)),
        legend_title_fontsize=float(plot_map.get("legend_title_fontsize", 0.0)),
    )

    embedding = EmbeddingConfig(
        n_top_genes=int(embedding_map.get("n_top_genes", 2000)),
        n_pcs=int(embedding_map.get("n_pcs", 20)),
        n_neighbors=int(embedding_map.get("n_neighbors", 15)),
        resolution=float(embedding_map.get("resolution", 1.0)),
        min_dist=float(embedding_map.get("min_dist", 0.5)),
        spread=float(embedding_map.get("spread", 1.0)),
        random_state=int(embedding_map.get("random_state", 0)),
    )

    azimuth = AzimuthConfig(
        enabled=bool(azimuth_map.get("enabled", False)),
        reference=str(azimuth_map.get("reference", "pbmcref")),
        annotation_levels=tuple(azimuth_map.get("annotation_levels", ("l1", "l2"))),
        k_weight=int(azimuth_map.get("k_weight", 50)),
        n_trees=int(azimuth_map.get("n_trees", 20)),
        mapping_score_k=int(azimuth_map.get("mapping_score_k", 100)),
    )

    tuning = TuningConfig(
        qc_preset_family=str(tuning_map.get("qc_preset_family", "default")),
        azimuth_preset_family=str(tuning_map.get("azimuth_preset_family", "default")),
        embedding_preset_family=str(
            tuning_map.get("embedding_preset_family", "default")
        ),
        max_candidates=int(tuning_map.get("max_candidates", 1)),
    )

    anno_map = data.get("annotation", {}) or {}
    annotation_cfg = None
    if isinstance(anno_map, dict):
        annotation_cfg = AnnotationConfig(methods=list(anno_map.get("methods", [])))

    return RunConfig(
        qc=qc,
        plotting=plotting,
        embedding=embedding,
        azimuth=azimuth,
        tuning=tuning,
        annotation=annotation_cfg,
    )


def load_default_config(path: Path) -> RunConfig:
    return load_run_config(path)


def merge_cli_overrides(base: RunConfig, **overrides) -> RunConfig:
    """Merge CLI-style overrides into a RunConfig and return a new RunConfig.

    Overrides are expected to be provided using nested keys in the form:
      qc__min_genes=350 or plotting__umap_width=11.5
    This function only touches the nested fields that are explicitly overridden,
    leaving all other fields unchanged.
    """
    # Start from the base values as dicts to allow selective overrides
    qc = {
        "method": base.qc.method,
        "counts_lower_nmads": base.qc.counts_lower_nmads,
        "genes_lower_nmads": base.qc.genes_lower_nmads,
        "pct_mt_upper_nmads": base.qc.pct_mt_upper_nmads,
        "pct_ribo_upper_nmads": base.qc.pct_ribo_upper_nmads,
        "min_cells_for_dynamic": base.qc.min_cells_for_dynamic,
        "count_floor_min": base.qc.count_floor_min,
        "count_floor_max": base.qc.count_floor_max,
        "gene_floor_min": base.qc.gene_floor_min,
        "gene_floor_max": base.qc.gene_floor_max,
        "pct_mt_ceiling_min": base.qc.pct_mt_ceiling_min,
        "pct_mt_ceiling_max": base.qc.pct_mt_ceiling_max,
        "pct_ribo_ceiling_min": base.qc.pct_ribo_ceiling_min,
        "pct_ribo_ceiling_max": base.qc.pct_ribo_ceiling_max,
    }
    plotting = {
        "umap_width": base.plotting.umap_width,
        "umap_height": base.plotting.umap_height,
        "dpi": base.plotting.dpi,
        "point_size": base.plotting.point_size,
        "legend_fontsize": base.plotting.legend_fontsize,
        "legend_title_fontsize": base.plotting.legend_title_fontsize,
    }
    embedding = {
        "n_top_genes": base.embedding.n_top_genes,
        "n_pcs": base.embedding.n_pcs,
        "n_neighbors": base.embedding.n_neighbors,
        "resolution": base.embedding.resolution,
        "min_dist": base.embedding.min_dist,
        "spread": base.embedding.spread,
        "random_state": base.embedding.random_state,
    }
    azimuth = {
        "enabled": base.azimuth.enabled,
        "reference": base.azimuth.reference,
        "annotation_levels": base.azimuth.annotation_levels,
        "k_weight": base.azimuth.k_weight,
        "n_trees": base.azimuth.n_trees,
        "mapping_score_k": base.azimuth.mapping_score_k,
    }
    tuning = {
        "qc_preset_family": base.tuning.qc_preset_family,
        "azimuth_preset_family": base.tuning.azimuth_preset_family,
        "embedding_preset_family": base.tuning.embedding_preset_family,
        "max_candidates": base.tuning.max_candidates,
    }
    annotation = base.annotation

    def _normalize_section(value: Any) -> Any:
        if is_dataclass(value):
            return asdict(cast(Any, value))
        return value

    for key, value in overrides.items():
        # Support plan-aligned simple top-level overrides first
        if key == "min_genes":
            qc["gene_floor_min"] = int(value)
            continue
        if key == "mt_max":
            qc["pct_mt_ceiling_max"] = float(value)
            continue
        if key in ("min_counts", "min_genes", "max_pct_mt", "max_pct_ribo"):
            # Overlays provided as simple overrides (nested under qc)
            # This keeps compatibility with existing CLI shapes like qc__min_genes=...
            # The actual assignment happens below via the nested handling.
            pass
        if "__" in key:
            section, field = key.split("__", 1)
            if section == "qc" and field in qc:
                qc[field] = value
            elif section == "plotting" and field in plotting:
                plotting[field] = value
            elif section == "embedding" and field in embedding:
                embedding[field] = value
            elif section == "azimuth" and field in azimuth:
                azimuth[field] = value
            elif section == "tuning" and field in tuning:
                tuning[field] = value
        else:
            # If top-level override is provided as a dict/dataclass, replace the whole section
            if key == "qc":
                qc = _normalize_section(value)  # type: ignore[assignment]
            elif key == "plotting":
                plotting = _normalize_section(value)  # type: ignore[assignment]
            elif key == "embedding":
                embedding = _normalize_section(value)  # type: ignore[assignment]
            elif key == "azimuth":
                azimuth = _normalize_section(value)  # type: ignore[assignment]
            elif key == "tuning":
                tuning = _normalize_section(value)  # type: ignore[assignment]

    new_qc = QcThresholds(
        method=str(qc["method"]),
        counts_lower_nmads=float(qc["counts_lower_nmads"]),
        genes_lower_nmads=float(qc["genes_lower_nmads"]),
        pct_mt_upper_nmads=float(qc["pct_mt_upper_nmads"]),
        pct_ribo_upper_nmads=float(qc["pct_ribo_upper_nmads"]),
        min_cells_for_dynamic=int(qc["min_cells_for_dynamic"]),
        count_floor_min=int(qc["count_floor_min"]),
        count_floor_max=int(qc["count_floor_max"]),
        gene_floor_min=int(qc["gene_floor_min"]),
        gene_floor_max=int(qc["gene_floor_max"]),
        pct_mt_ceiling_min=float(qc["pct_mt_ceiling_min"]),
        pct_mt_ceiling_max=float(qc["pct_mt_ceiling_max"]),
        pct_ribo_ceiling_min=float(qc["pct_ribo_ceiling_min"]),
        pct_ribo_ceiling_max=float(qc["pct_ribo_ceiling_max"]),
    )
    new_plotting = PlottingConfig(
        umap_width=float(plotting["umap_width"]),
        umap_height=float(plotting["umap_height"]),
        dpi=int(plotting["dpi"]),
        point_size=float(plotting["point_size"]),
        legend_fontsize=float(plotting["legend_fontsize"]),
        legend_title_fontsize=float(plotting["legend_title_fontsize"]),
    )

    new_embedding = EmbeddingConfig(
        n_top_genes=int(embedding["n_top_genes"]),
        n_pcs=int(embedding["n_pcs"]),
        n_neighbors=int(embedding["n_neighbors"]),
        resolution=float(embedding["resolution"]),
        min_dist=float(embedding["min_dist"]),
        spread=float(embedding["spread"]),
        random_state=int(embedding["random_state"]),
    )

    new_azimuth = AzimuthConfig(
        enabled=bool(azimuth["enabled"]),
        reference=str(azimuth["reference"]),
        annotation_levels=tuple(
            cast(Iterable[str], azimuth["annotation_levels"])
        ),
        k_weight=int(cast(Any, azimuth["k_weight"])),
        n_trees=int(cast(Any, azimuth["n_trees"])),
        mapping_score_k=int(cast(Any, azimuth["mapping_score_k"])),
    )

    new_tuning = TuningConfig(
        qc_preset_family=str(tuning["qc_preset_family"]),
        azimuth_preset_family=str(tuning["azimuth_preset_family"]),
        embedding_preset_family=str(tuning["embedding_preset_family"]),
        max_candidates=int(tuning["max_candidates"]),
    )

    return RunConfig(
        qc=new_qc,
        plotting=new_plotting,
        embedding=new_embedding,
        azimuth=new_azimuth,
        tuning=new_tuning,
        annotation=annotation,
    )
