from __future__ import annotations

from pathlib import Path
from dataclasses import asdict, is_dataclass
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
        min_counts=int(qc_map.get("min_counts", 0)),
        min_genes=int(qc_map.get("min_genes", 0)),
        max_pct_mt=float(qc_map.get("max_pct_mt", 0.0)),
        max_pct_ribo=float(qc_map.get("max_pct_ribo", 0.0)),
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
        max_candidates=int(tuning_map.get("max_candidates", 9)),
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
        "min_counts": base.qc.min_counts,
        "min_genes": base.qc.min_genes,
        "max_pct_mt": base.qc.max_pct_mt,
        "max_pct_ribo": base.qc.max_pct_ribo,
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

    def _normalize_section(value):
        if is_dataclass(value):
            return asdict(value)
        return value

    for key, value in overrides.items():
        # Support plan-aligned simple top-level overrides first
        if key == "min_genes":
            qc["min_genes"] = int(value)
            continue
        if key == "mt_max":
            qc["max_pct_mt"] = float(value)
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
        min_counts=int(qc["min_counts"]),
        min_genes=int(qc["min_genes"]),
        max_pct_mt=float(qc["max_pct_mt"]),
        max_pct_ribo=float(qc["max_pct_ribo"]),
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
        annotation_levels=tuple(azimuth["annotation_levels"]),
        k_weight=int(azimuth["k_weight"]),
        n_trees=int(azimuth["n_trees"]),
        mapping_score_k=int(azimuth["mapping_score_k"]),
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
