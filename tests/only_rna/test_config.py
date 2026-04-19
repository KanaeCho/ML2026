from __future__ import annotations

from pathlib import Path

from scripts.only_rna import config as config_module
from scripts.only_rna.models import AzimuthConfig, EmbeddingConfig, QcThresholds


REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_PATH = REPO_ROOT / "scripts" / "only_rna" / "default_config.yaml"


def test_load_run_config_reads_embedding_azimuth_and_tuning_sections(
    tmp_path: Path,
) -> None:
    assert hasattr(config_module, "load_run_config")

    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
qc:
  min_counts: 600
  min_genes: 350
  max_pct_mt: 18.0
  max_pct_ribo: 55.0
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

    assert config.qc.min_counts == 600
    assert config.embedding.n_top_genes == 1500
    assert config.embedding.min_dist == 0.35
    assert config.azimuth.enabled is True
    assert config.azimuth.reference == "pbmcref"
    assert list(config.azimuth.annotation_levels) == ["l1", "l2"]
    assert config.tuning.max_candidates == 9


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
            min_counts=800,
            min_genes=500,
            max_pct_mt=15.0,
            max_pct_ribo=50.0,
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

    assert updated.qc.min_counts == 800
    assert updated.qc.min_genes == 500
    assert updated.embedding.n_neighbors == 20
    assert updated.azimuth.k_weight == 30
    assert updated.azimuth.n_trees == 40
