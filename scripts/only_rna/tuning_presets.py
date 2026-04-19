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
                min_counts=500,
                min_genes=300,
                max_pct_mt=20.0,
                max_pct_ribo=60.0,
            ),
            "strict": QcThresholds(
                min_counts=800,
                min_genes=500,
                max_pct_mt=15.0,
                max_pct_ribo=50.0,
            ),
            "lenient": QcThresholds(
                min_counts=300,
                min_genes=200,
                max_pct_mt=25.0,
                max_pct_ribo=70.0,
            ),
        },
        azimuth={
            "baseline": AzimuthConfig(enabled=True),
            "conservative": AzimuthConfig(
                enabled=True,
                k_weight=30,
                n_trees=40,
                mapping_score_k=50,
            ),
            "smooth": AzimuthConfig(
                enabled=True,
                k_weight=80,
                n_trees=20,
                mapping_score_k=150,
            ),
        },
        embedding={
            "baseline": EmbeddingConfig(),
            "separated": EmbeddingConfig(
                n_neighbors=10,
                resolution=1.2,
                min_dist=0.15,
                spread=1.2,
            ),
            "stable": EmbeddingConfig(
                n_neighbors=20,
                resolution=0.6,
                min_dist=0.5,
                spread=1.0,
            ),
        },
    )


__all__ = ["TuningPresetFamilies", "default_tuning_presets"]
