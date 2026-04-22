from __future__ import annotations

from dataclasses import dataclass

from .models import AzimuthConfig, EmbeddingConfig, QcThresholds


@dataclass(frozen=True)
class TuningPresetFamilies:
    qc: dict[str, QcThresholds]
    azimuth: dict[str, AzimuthConfig]
    embedding: dict[str, EmbeddingConfig]


def default_tuning_presets() -> TuningPresetFamilies:
    return TuningPresetFamilies(
        qc={
            "baseline": QcThresholds(
                method="dynamic_hybrid_mad",
                counts_lower_nmads=3.0,
                genes_lower_nmads=3.0,
                pct_mt_upper_nmads=3.0,
                pct_ribo_upper_nmads=3.5,
                min_cells_for_dynamic=50,
                count_floor_min=100,
                count_floor_max=1500,
                gene_floor_min=100,
                gene_floor_max=1200,
                pct_mt_ceiling_min=5.0,
                pct_mt_ceiling_max=40.0,
                pct_ribo_ceiling_min=20.0,
                pct_ribo_ceiling_max=80.0,
            ),
        },
        azimuth={
            "baseline": AzimuthConfig(enabled=True),
        },
        embedding={
            "baseline": EmbeddingConfig(),
        },
    )


__all__ = ["TuningPresetFamilies", "default_tuning_presets"]
