from __future__ import annotations

import gzip
import tarfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.io import mmwrite
from scipy import sparse

from scripts.only_rna.config import load_default_config, merge_cli_overrides
from scripts.only_rna.annotation import (
    annotate_with_cima,
    load_cima_reference,
)
from scripts.only_rna.discovery import DiscoveredSample
from scripts.only_rna.doublet import run_doublet_detection
from scripts.only_rna.embedding import run_embedding
from scripts.only_rna.models import RunConfig, QcThresholds, PlottingConfig
from scripts.only_rna.qc import apply_qc_filters, compute_qc_metrics
from scripts.only_rna.read_inputs import read_sample_input


def test_load_default_config():
    cfg = load_default_config(Path("scripts/only_rna/default_config.yaml"))
    assert isinstance(cfg, RunConfig)
    assert isinstance(cfg.qc, QcThresholds)
    assert isinstance(cfg.plotting, PlottingConfig)


def test_defaults_values_taken_from_yaml():
    cfg = load_default_config(Path("scripts/only_rna/default_config.yaml"))
    assert cfg.qc.min_counts == 500
    assert cfg.qc.min_genes == 300
    assert cfg.qc.max_pct_mt == 20.0
    assert cfg.qc.max_pct_ribo == 60.0


def test_plotting_dimensions_positive():
    cfg = load_default_config(Path("scripts/only_rna/default_config.yaml"))
    assert cfg.plotting.umap_width > 0
    assert cfg.plotting.umap_height > 0


def test_merge_cli_overrides_only_touch_selected():
    base = load_default_config(Path("scripts/only_rna/default_config.yaml"))
    # Use plan-aligned simple override names
    merged = merge_cli_overrides(base, min_genes=350, mt_max=18.0)

    assert isinstance(merged, RunConfig)
    assert merged.qc.min_genes == 350
    assert merged.qc.max_pct_mt == 18.0
    # Ensure other fields are preserved from base
    assert merged.qc.min_counts == base.qc.min_counts
    assert merged.plotting.umap_width == base.plotting.umap_width
    assert merged.plotting.umap_height == base.plotting.umap_height
    # ensure untouched fields preserved


def _make_run_config() -> RunConfig:
    return RunConfig(
        qc=QcThresholds(
            min_counts=5,
            min_genes=2,
            max_pct_mt=20.0,
            max_pct_ribo=40.0,
        ),
        plotting=PlottingConfig(
            umap_width=8.0,
            umap_height=6.0,
            dpi=100,
            point_size=3.0,
            legend_fontsize=8.0,
            legend_title_fontsize=9.0,
        ),
    )


def _make_qc_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [1, 1, 0, 8],
                [0, 0, 0, 2],
                [5, 4, 1, 0],
                [0, 1, 1, 8],
            ],
            dtype=float,
        ),
        obs=pd.DataFrame(index=["cell-1", "cell-2", "cell-3", "cell-4"]),
        var=pd.DataFrame(index=["MT-CO1", "RPS3", "RPL5", "GeneA"]),
    )


def test_compute_qc_metrics_adds_expected_columns():
    adata = _make_qc_adata()
    config = _make_run_config()

    out = compute_qc_metrics(adata, config)

    for key in ["n_counts", "n_genes", "pct_mt", "pct_ribo"]:
        assert key in out.obs

    assert out.obs["n_counts"].tolist() == [10.0, 2.0, 10.0, 10.0]
    assert out.obs["n_genes"].tolist() == [3, 1, 3, 3]
    assert out.obs["pct_mt"].tolist() == [10.0, 0.0, 50.0, 0.0]
    assert out.obs["pct_ribo"].tolist() == [10.0, 0.0, 50.0, 20.0]


def test_apply_qc_filters_sets_fail_flags_and_doublet_hard_filter():
    adata = _make_qc_adata()
    config = _make_run_config()

    with_metrics = compute_qc_metrics(adata, config)
    with_metrics.obs["is_doublet"] = [False, False, False, True]

    out = apply_qc_filters(with_metrics, config)

    assert out.obs["fails_count_floor"].tolist() == [False, True, False, False]
    assert out.obs["fails_gene_floor"].tolist() == [False, True, False, False]
    assert out.obs["fails_mt_ceiling"].tolist() == [False, False, True, False]
    assert out.obs["fails_ribo_ceiling"].tolist() == [False, False, True, False]
    assert out.obs["fails_doublet"].tolist() == [False, False, False, True]
    assert out.obs["pass_qc"].tolist() == [True, False, False, False]


def test_run_doublet_detection_normalizes_existing_columns_without_scrublet(
    monkeypatch,
):
    adata = _make_qc_adata()
    config = _make_run_config()
    adata.obs["doublet_score"] = ["0.1", None, 0.0, 0.2]
    adata.obs["is_doublet"] = [0, 1, False, True]

    def fail_scrublet(*args, **kwargs):
        raise AssertionError("scrublet should not be called when columns exist")

    monkeypatch.setattr("scripts.only_rna.doublet.sc.pp.scrublet", fail_scrublet)

    out = run_doublet_detection(adata, config)

    assert out is not adata
    assert out.obs["doublet_score"].tolist() == [0.1, 0.0, 0.0, 0.2]
    assert out.obs["is_doublet"].tolist() == [False, True, False, True]


def test_run_doublet_detection_writes_columns_from_scrublet(monkeypatch):
    adata = _make_qc_adata()
    config = _make_run_config()

    def fake_scrublet(target_adata, **kwargs):
        target_adata.obs["doublet_score"] = [0.01, 0.8, 0.2, 0.9]
        target_adata.obs["predicted_doublet"] = [False, True, False, True]

    monkeypatch.setattr("scripts.only_rna.doublet.sc.pp.scrublet", fake_scrublet)

    out = run_doublet_detection(adata, config)

    assert "doublet_score" in out.obs
    assert "is_doublet" in out.obs
    assert out.obs["doublet_score"].tolist() == [0.01, 0.8, 0.2, 0.9]
    assert out.obs["is_doublet"].tolist() == [False, True, False, True]


def _write_text_lines(path: Path, lines: list[str]) -> None:
    content = "\n".join(lines) + "\n"
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(content)
        return

    path.write_text(content, encoding="utf-8")


def _write_triplet(
    base_dir: Path,
    *,
    feature_filename: str = "features.tsv",
    feature_rows: list[list[str]] | None = None,
    counts: np.ndarray | None = None,
) -> tuple[Path, Path, Path]:
    matrix = np.asarray(
        counts
        if counts is not None
        else np.array(
            [
                [1, 0],
                [0, 2],
                [3, 4],
                [5, 6],
            ]
        )
    )
    matrix_path = base_dir / "matrix.mtx"
    mmwrite(matrix_path, sparse.coo_matrix(matrix))

    barcodes_path = base_dir / "barcodes.tsv"
    _write_text_lines(barcodes_path, ["cell-1", "cell-2"])

    features_path = base_dir / feature_filename
    rows = feature_rows or [
        ["ENSG0001", "GeneA", "Gene Expression"],
        ["ENSG0002", "ADT1", "Antibody Capture"],
        ["ENSG0003", "GeneA", "Gene Expression"],
        ["ENSG0004", "GeneC", "Gene Expression"],
    ]
    _write_text_lines(features_path, ["\t".join(row) for row in rows])
    return matrix_path, barcodes_path, features_path


def _make_sample(
    *,
    gse: str = "GSE123456",
    sample_id: str = "GSM123456",
    input_type: str = "triplet",
    sample_kind: str = "gsm",
    matrix_path: Path | None = None,
    barcodes_path: Path | None = None,
    features_path: Path | None = None,
    h5_path: Path | None = None,
    archive_path: Path | None = None,
) -> DiscoveredSample:
    return DiscoveredSample(
        gse=gse,
        sample_id=sample_id,
        input_type=input_type,
        sample_kind=sample_kind,
        supported=True,
        note="fixture",
        source_name=sample_id,
        matrix_path=matrix_path,
        barcodes_path=barcodes_path,
        features_path=features_path,
        h5_path=h5_path,
        archive_path=archive_path,
    )


def test_read_triplet_returns_adata_with_sample_metadata(tmp_path: Path):
    matrix_path, barcodes_path, features_path = _write_triplet(tmp_path)
    sample = _make_sample(
        matrix_path=matrix_path,
        barcodes_path=barcodes_path,
        features_path=features_path,
    )

    adata = read_sample_input(sample)

    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert list(adata.obs_names) == ["cell-1", "cell-2"]
    assert list(adata.var_names) == ["GeneA", "GeneA-1", "GeneC"]
    assert adata.obs["gse"].nunique() == 1
    assert adata.obs["gse"].iloc[0] == sample.gse
    assert adata.obs["sample_id"].nunique() == 1
    assert adata.obs["sample_id"].iloc[0] == sample.sample_id
    assert adata.obs["input_type"].iloc[0] == "triplet"


def test_read_shared_gse_triplet_sets_gse_sample_id(tmp_path: Path):
    matrix_path, barcodes_path, features_path = _write_triplet(tmp_path)
    sample = _make_sample(
        gse="GSE149689",
        sample_id="GSE149689",
        sample_kind="gse_shared",
        matrix_path=matrix_path,
        barcodes_path=barcodes_path,
        features_path=features_path,
    )

    adata = read_sample_input(sample)

    assert adata.n_obs == 2
    assert adata.obs["sample_id"].nunique() == 1
    assert adata.obs["sample_id"].iloc[0] == "GSE149689"
    assert adata.obs["gse"].iloc[0] == "GSE149689"


def test_read_triplet_supports_genes_tsv_gz(tmp_path: Path):
    matrix_path, barcodes_path, features_path = _write_triplet(
        tmp_path,
        feature_filename="genes.tsv.gz",
        feature_rows=[
            ["ENSG0001", "GeneX"],
            ["ENSG0002", "GeneY"],
            ["ENSG0003", "GeneY"],
        ],
        counts=np.array(
            [
                [1, 0],
                [0, 2],
                [3, 4],
            ]
        ),
    )
    sample = _make_sample(
        matrix_path=matrix_path,
        barcodes_path=barcodes_path,
        features_path=features_path,
    )

    adata = read_sample_input(sample)

    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert list(adata.var_names) == ["GeneX", "GeneY", "GeneY-1"]
    assert list(adata.var["feature_id"]) == ["ENSG0001", "ENSG0002", "ENSG0003"]


def test_read_archive_input_unpacks_and_reads_triplet(tmp_path: Path):
    extracted_dir = tmp_path / "nested" / "filtered_feature_bc_matrix"
    extracted_dir.mkdir(parents=True)
    _write_triplet(extracted_dir)

    archive_path = tmp_path / "matrix_bundle.tar.gz"
    with tarfile.open(archive_path, "w:gz") as handle:
        handle.add(extracted_dir.parent, arcname="bundle")

    sample = _make_sample(
        input_type="archive",
        archive_path=archive_path,
        matrix_path=None,
        barcodes_path=None,
        features_path=None,
    )

    adata = read_sample_input(sample)

    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert adata.obs["input_type"].iloc[0] == "archive"


def test_read_h5_prefers_gene_expression_when_present(tmp_path: Path, monkeypatch):
    h5_path = tmp_path / "sample.h5"
    h5_path.write_bytes(b"placeholder")
    sample = _make_sample(
        input_type="h5",
        h5_path=h5_path,
        matrix_path=None,
        barcodes_path=None,
        features_path=None,
    )

    calls: dict[str, object] = {}

    def fake_read_10x_h5(path: Path, *, gex_only: bool = True):
        calls["path"] = path
        calls["gex_only"] = gex_only
        return ad.AnnData(
            X=np.array(
                [
                    [1, 2, 3],
                    [4, 5, 6],
                ]
            ),
            obs=pd.DataFrame(index=["cell-1", "cell-2"]),
            var=pd.DataFrame(
                {
                    "feature_types": [
                        "Gene Expression",
                        "Antibody Capture",
                        "Gene Expression",
                    ]
                },
                index=["GeneA", "ADT1", "GeneA"],
            ),
        )

    monkeypatch.setattr("scripts.only_rna.read_inputs.sc.read_10x_h5", fake_read_10x_h5)

    adata = read_sample_input(sample)

    assert calls == {"path": h5_path, "gex_only": False}
    assert adata.n_obs == 2
    assert adata.n_vars == 2
    assert list(adata.var_names) == ["GeneA", "GeneA-1"]
    assert adata.obs["sample_id"].iloc[0] == sample.sample_id
    assert adata.obs["input_type"].iloc[0] == "h5"


def _make_embedding_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [5.0, 0.0, 1.0],
                [0.0, 2.0, 0.0],
                [4.0, 1.0, 0.0],
            ]
        ),
        obs=pd.DataFrame(
            {"pass_qc": [True, False, True]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )


def test_run_embedding_writes_cluster_and_umap_for_pass_qc_only(monkeypatch):
    adata = _make_embedding_adata()
    config = _make_run_config()

    def fake_embedding(pass_qc_adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
        del config
        out = pass_qc_adata.copy()
        out.obs["cluster"] = pd.Series(["0", "1"], index=out.obs_names, dtype="string")
        out.obsm["X_umap"] = np.array([[1.5, 2.5], [3.5, 4.5]], dtype=float)
        return out

    monkeypatch.setattr(
        "scripts.only_rna.embedding._run_scanpy_embedding", fake_embedding
    )

    out = run_embedding(adata, config)

    assert out is not adata
    assert list(out.obs.loc[["cell-1", "cell-3"], "cluster"]) == ["0", "1"]
    assert out.obs.loc["cell-2", "cluster"] is pd.NA or pd.isna(
        out.obs.loc["cell-2", "cluster"]
    )
    assert out.obs.loc[["cell-1", "cell-3"], "umap_1"].tolist() == [1.5, 3.5]
    assert out.obs.loc[["cell-1", "cell-3"], "umap_2"].tolist() == [2.5, 4.5]
    assert pd.isna(out.obs.loc["cell-2", "umap_1"])
    assert pd.isna(out.obs.loc["cell-2", "umap_2"])


def test_run_embedding_falls_back_when_leiden_dependency_missing(monkeypatch):
    adata = ad.AnnData(
        X=np.array(
            [
                [5.0, 1.0, 0.0],
                [0.0, 4.0, 1.0],
                [1.0, 0.0, 5.0],
                [2.0, 2.0, 2.0],
            ]
        ),
        obs=pd.DataFrame(
            {"pass_qc": [True, True, True, False]},
            index=["cell-1", "cell-2", "cell-3", "cell-4"],
        ),
        var=pd.DataFrame(index=["GeneA", "GeneB", "GeneC"]),
    )
    config = _make_run_config()

    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.pp.normalize_total", lambda *a, **k: None
    )
    monkeypatch.setattr("scripts.only_rna.embedding.sc.pp.log1p", lambda *a, **k: None)
    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.pp.highly_variable_genes", lambda *a, **k: None
    )
    monkeypatch.setattr("scripts.only_rna.embedding.sc.pp.scale", lambda *a, **k: None)
    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.tl.pca",
        lambda current, **kwargs: current.obsm.__setitem__(
            "X_pca", np.asarray(current.X, dtype=float)
        ),
    )
    monkeypatch.setattr(
        "scripts.only_rna.embedding.sc.pp.neighbors", lambda *a, **k: None
    )

    def raise_missing_igraph(*args, **kwargs):
        raise ImportError("Please install the igraph package")

    def fake_umap(current, **kwargs):
        del kwargs
        current.obsm["X_umap"] = np.array(
            [[10.0, 1.0], [11.5, 2.5], [13.0, 4.0]], dtype=float
        )

    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.leiden", raise_missing_igraph)
    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.umap", fake_umap)

    out = run_embedding(adata, config)

    assert list(out.obs.loc[["cell-1", "cell-2", "cell-3"], "cluster"]) == [
        "0",
        "0",
        "0",
    ]
    assert out.obs.loc[["cell-1", "cell-2", "cell-3"], "umap_1"].tolist() == [
        10.0,
        11.5,
        13.0,
    ]
    assert out.obs.loc[["cell-1", "cell-2", "cell-3"], "umap_2"].tolist() == [
        1.0,
        2.5,
        4.0,
    ]
    assert out.obs.loc["cell-4", "cluster"] is pd.NA or pd.isna(
        out.obs.loc["cell-4", "cluster"]
    )
    assert pd.isna(out.obs.loc["cell-4", "umap_1"])
    assert pd.isna(out.obs.loc["cell-4", "umap_2"])


def _write_cima_reference_assets(
    reference_dir: Path,
    *,
    feature_rows: list[dict[str, object]] | None = None,
    l1_rows: list[dict[str, object]] | None = None,
    l2_rows: list[dict[str, object]] | None = None,
    hierarchy_rows: list[dict[str, object]] | None = None,
) -> Path:
    cima_dir = reference_dir / "cima"
    cima_dir.mkdir(parents=True, exist_ok=True)

    feature_table = pd.DataFrame(
        feature_rows
        or [
            {
                "feature_id": "GeneA",
                "gene_mean": 0.0,
                "gene_std": 1.0,
                "pc_dim_1": 1.0,
                "pc_dim_2": 0.0,
            },
            {
                "feature_id": "GeneB",
                "gene_mean": 0.0,
                "gene_std": 1.0,
                "pc_dim_1": 0.0,
                "pc_dim_2": 1.0,
            },
        ]
    )
    with gzip.open(
        cima_dir / "cima_rna_reference_pca_features.tsv.gz", "wt", encoding="utf-8"
    ) as handle:
        feature_table.to_csv(handle, sep="\t", index=False)

    pd.DataFrame(
        l1_rows
        or [
            {"cell_type_l1": "L1_A", "pc_dim_1": 1.0, "pc_dim_2": 0.0},
            {"cell_type_l1": "L1_B", "pc_dim_1": 0.0, "pc_dim_2": 1.0},
        ]
    ).to_csv(cima_dir / "cima_rna_reference_l1_centroids.tsv", sep="\t", index=False)

    pd.DataFrame(
        l2_rows
        or [
            {"cell_type_l2": "L2_A1", "pc_dim_1": 0.8, "pc_dim_2": 0.2},
            {"cell_type_l2": "L2_A2", "pc_dim_1": 0.6, "pc_dim_2": 0.4},
            {"cell_type_l2": "L2_B1", "pc_dim_1": 0.1, "pc_dim_2": 0.95},
        ]
    ).to_csv(cima_dir / "cima_rna_reference_l2_centroids.tsv", sep="\t", index=False)

    pd.DataFrame(
        hierarchy_rows
        or [
            {"cell_type_l1": "L1_A", "cell_type_l2": "L2_A1"},
            {"cell_type_l1": "L1_A", "cell_type_l2": "L2_A2"},
            {"cell_type_l1": "L1_B", "cell_type_l2": "L2_B1"},
        ]
    ).to_csv(cima_dir / "cima_rna_celltype_hierarchy.csv", index=False)

    return cima_dir


def _make_annotation_adata() -> ad.AnnData:
    return ad.AnnData(
        X=np.array(
            [
                [10.0, 0.0],
                [0.0, 10.0],
                [1.0, 1.0],
            ]
        ),
        obs=pd.DataFrame(
            {"pass_qc": [True, False, True]},
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(
            {"feature_id": ["GeneA", "GeneB"]},
            index=["gene-a", "gene-b"],
        ),
    )


def test_load_cima_reference_reads_expected_files_and_hierarchy(tmp_path: Path):
    _write_cima_reference_assets(tmp_path)

    reference = load_cima_reference(tmp_path)

    assert reference.feature_ids == ["GeneA", "GeneB"]
    assert reference.loadings.shape == (2, 2)
    assert list(reference.l1_centroids.index) == ["L1_A", "L1_B"]
    assert list(reference.l2_centroids.index) == ["L2_A1", "L2_A2", "L2_B1"]
    assert reference.l2_by_l1 == {
        "L1_A": ["L2_A1", "L2_A2"],
        "L1_B": ["L2_B1"],
    }


def test_annotate_with_cima_only_writes_pass_qc_and_applies_masking(tmp_path: Path):
    _write_cima_reference_assets(tmp_path)
    adata = _make_annotation_adata()

    out = annotate_with_cima(adata, tmp_path)

    assert out is not adata
    assert out.obs.loc["cell-1", "cima_l1"] == "L1_A"
    assert out.obs.loc["cell-1", "cima_l2"] == "L2_A1"
    assert bool(out.obs.loc["cell-1", "cima_l1_low_confidence"]) is False
    assert out.obs.loc["cell-1", "cima_l1_masked"] == "L1_A"

    assert out.obs.loc["cell-3", "cima_l1"] == "L1_A"
    assert bool(out.obs.loc["cell-3", "cima_l1_low_confidence"]) is True
    assert out.obs.loc["cell-3", "cima_l1_masked"] == "Unknown"

    for column in [
        "cima_l1",
        "cima_l2",
        "cima_l1_score",
        "cima_l1_score_margin",
        "cima_l2_score",
        "cima_l2_score_margin",
        "cima_l1_low_confidence",
        "cima_l1_masked",
    ]:
        assert pd.isna(out.obs.loc["cell-2", column])


def test_annotate_with_cima_constrains_l2_by_l1_hierarchy(tmp_path: Path):
    _write_cima_reference_assets(
        tmp_path,
        l2_rows=[
            {"cell_type_l2": "L2_A1", "pc_dim_1": 0.6, "pc_dim_2": 0.4},
            {"cell_type_l2": "L2_B1", "pc_dim_1": 0.99, "pc_dim_2": 0.01},
        ],
        hierarchy_rows=[
            {"cell_type_l1": "L1_A", "cell_type_l2": "L2_A1"},
            {"cell_type_l1": "L1_B", "cell_type_l2": "L2_B1"},
        ],
    )
    adata = ad.AnnData(
        X=np.array([[10.0, 0.0]]),
        obs=pd.DataFrame({"pass_qc": [True]}, index=["cell-1"]),
        var=pd.DataFrame({"feature_id": ["GeneA", "GeneB"]}, index=["g1", "g2"]),
    )

    out = annotate_with_cima(adata, tmp_path)

    assert out.obs.loc["cell-1", "cima_l1"] == "L1_A"
    assert out.obs.loc["cell-1", "cima_l2"] == "L2_A1"
