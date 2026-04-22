from __future__ import annotations

from pathlib import Path

from scripts.only_rna import config as config_module
from scripts.only_rna.models import AzimuthConfig, EmbeddingConfig, QcThresholds


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = REPO_ROOT / "scripts" / "only_rna" / "default_config.yaml"


def test_default_config_pivots_mainline_to_azimuth() -> None:
    config = config_module.load_run_config(DEFAULT_CONFIG_PATH)

    assert config.annotation is not None
    assert config.annotation.methods == ["azimuth"]
    assert config.azimuth.enabled is True
    assert config.azimuth.reference == "pbmcref"


def test_load_run_config_reads_embedding_azimuth_and_tuning_sections(
    tmp_path: Path,
) -> None:
    assert hasattr(config_module, "load_run_config")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
qc:
  method: dynamic_hybrid_mad
  counts_lower_nmads: 2.5
  genes_lower_nmads: 2.8
  pct_mt_upper_nmads: 3.2
  pct_ribo_upper_nmads: 3.7
  min_cells_for_dynamic: 80
  count_floor_min: 120
  count_floor_max: 1600
  gene_floor_min: 140
  gene_floor_max: 1400
  pct_mt_ceiling_min: 6.0
  pct_mt_ceiling_max: 35.0
  pct_ribo_ceiling_min: 25.0
  pct_ribo_ceiling_max: 75.0
plotting:
  umap_width: 10.0
  umap_height: 8.0
  dpi: 180
  point_size: 3.0
  legend_fontsize: 10.0
  legend_title_fontsize: 11.0
embedding:
  n_top_genes: 1500
  n_pcs: 20
  n_neighbors: 12
  resolution: 0.8
  min_dist: 0.35
azimuth:
  enabled: true
  reference: pbmcref
  annotation_levels:
    - l1
    - l2
  k_weight: 50
  mapping_score_k: 100
tuning:
  qc_preset_family: default
  azimuth_preset_family: default
  embedding_preset_family: default
  max_candidates: 9
""".strip(),
        encoding="utf-8",
    )

    config = config_module.load_run_config(config_path)

    assert config.qc.method == "dynamic_hybrid_mad"
    assert config.qc.counts_lower_nmads == 2.5
    assert config.qc.genes_lower_nmads == 2.8
    assert config.qc.min_cells_for_dynamic == 80
    assert config.qc.count_floor_min == 120
    assert config.qc.pct_mt_ceiling_max == 35.0
    assert config.embedding.n_top_genes == 1500
    assert config.embedding.min_dist == 0.35
    assert config.azimuth.enabled is True
    assert config.azimuth.reference == "pbmcref"
    assert list(config.azimuth.annotation_levels) == ["l1", "l2"]
    assert config.tuning.max_candidates == 9


def test_default_config_limits_tuning_to_single_baseline_candidate() -> None:
    config = config_module.load_run_config(DEFAULT_CONFIG_PATH)

    assert config.tuning.qc_preset_family == "baseline_only"
    assert config.tuning.azimuth_preset_family == "baseline_only"
    assert config.tuning.embedding_preset_family == "baseline_only"
    assert config.tuning.max_candidates == 1


def test_merge_cli_overrides_updates_embedding_and_azimuth_fields() -> None:
    assert hasattr(config_module, "load_run_config")

    config = config_module.load_run_config(DEFAULT_CONFIG_PATH)

    updated = config_module.merge_cli_overrides(
        config,
        embedding__n_neighbors=18,
        embedding__min_dist=0.2,
        azimuth__k_weight=80,
    )

    assert updated.embedding.n_neighbors == 18
    assert updated.embedding.min_dist == 0.2
    assert updated.azimuth.k_weight == 80


def test_merge_cli_overrides_accepts_runtime_config_dataclasses() -> None:
    config = config_module.load_run_config(DEFAULT_CONFIG_PATH)

    updated = config_module.merge_cli_overrides(
        config,
        qc=QcThresholds(
            method="dynamic_hybrid_mad",
            counts_lower_nmads=2.4,
            genes_lower_nmads=2.9,
            pct_mt_upper_nmads=3.1,
            pct_ribo_upper_nmads=3.8,
            min_cells_for_dynamic=75,
            count_floor_min=130,
            count_floor_max=1700,
            gene_floor_min=160,
            gene_floor_max=1500,
            pct_mt_ceiling_min=7.0,
            pct_mt_ceiling_max=30.0,
            pct_ribo_ceiling_min=22.0,
            pct_ribo_ceiling_max=70.0,
        ),
        embedding=EmbeddingConfig(
            n_top_genes=2000,
            n_pcs=20,
            n_neighbors=20,
            resolution=0.6,
            min_dist=0.5,
            spread=1.0,
            random_state=0,
        ),
        azimuth=AzimuthConfig(
            enabled=True,
            reference="pbmcref",
            annotation_levels=("l1", "l2"),
            k_weight=30,
            n_trees=40,
            mapping_score_k=50,
        ),
    )

    assert updated.qc.counts_lower_nmads == 2.4
    assert updated.qc.genes_lower_nmads == 2.9
    assert updated.qc.count_floor_min == 130
    assert updated.qc.pct_mt_ceiling_max == 30.0
    assert updated.embedding.n_neighbors == 20
    assert updated.azimuth.k_weight == 30
    assert updated.azimuth.n_trees == 40
