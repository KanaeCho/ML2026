from __future__ import annotations

from pathlib import Path
import yaml

from .models import RunConfig, QcThresholds, PlottingConfig, AnnotationConfig


def load_default_config(path: Path) -> RunConfig:
    """Load the YAML default RunConfig from the given path.

    YAML structure is expected to have two top-level keys:
      - qc: { min_counts, min_genes, max_pct_mt, max_pct_ribo }
      - plotting: { umap_width, umap_height, dpi, point_size, legend_fontsize, legend_title_fontsize }
    """
    path = Path(path)
    data = yaml.safe_load(path.read_text()) or {}
    qc_map = data.get("qc", {})
    plot_map = data.get("plotting", {})

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

    anno_map = data.get("annotation", {}) or {}
    annotation_cfg = None
    if isinstance(anno_map, dict):
        annotation_cfg = AnnotationConfig(methods=list(anno_map.get("methods", [])))

    return RunConfig(qc=qc, plotting=plotting, annotation=annotation_cfg)


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
    annotation = base.annotation

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
        else:
            # If top-level override is provided as a dict, replace the whole section
            if key == "qc":
                qc = value  # type: ignore[assignment]
            elif key == "plotting":
                plotting = value  # type: ignore[assignment]

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

    return RunConfig(qc=new_qc, plotting=new_plotting, annotation=annotation)
