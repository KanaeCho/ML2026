from __future__ import annotations

import json
from pathlib import Path
from collections.abc import Mapping, Sequence

import pandas as pd

from scripts.process.organize_integrated_products import (
    default_integration_method,
    discover_co_atac,
    discover_co_rna,
    discover_only_atac,
    discover_only_rna,
    organize_product,
)
from scripts.process.render_product_umap_panels import render_panels
from scripts.process.integrate_product_embeddings import compute_metrics


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def write_status(path: Path) -> None:
    path.write_text(json.dumps({"status": "success", "outputs_complete": True}), encoding="utf-8")


def make_rna_sample(root: Path, prefix: Path, sample_id: str) -> Path:
    sample_dir = root / prefix / sample_id
    rows = [
        {
            "cell_id": f"{sample_id}_cell1",
            "pass_qc": True,
            "umap_1": 1.0,
            "umap_2": 2.0,
            "azimuth_cima_l1": "B",
            "cima_l2": "B naive",
        }
    ]
    write_csv(sample_dir / "metadata.csv", rows)
    write_csv(sample_dir / "metadata_qc.csv", rows)
    write_csv(sample_dir / "qc_summary.csv", [{"n_cells_pass_qc": 1}])
    write_csv(sample_dir / "validation_result.csv", [{"completion": True}])
    write_status(sample_dir / "run_status.json")
    (sample_dir / "qc_overview.png").write_text("png")
    (sample_dir / "umap_rna_cima_l1.png").write_text("png")
    (sample_dir / f"{sample_id}.h5ad").write_text("h5ad")
    (sample_dir / "matrix").mkdir(parents=True)
    (sample_dir / "matrix" / "matrix.mtx").write_text("matrix")
    (sample_dir / "matrix" / "barcodes.tsv.gz").write_text("barcodes")
    (sample_dir / "matrix" / "features.tsv.gz").write_text("features")
    return sample_dir


def make_atac_sample(root: Path, prefix: Path, sample_id: str, legacy: bool = False) -> Path:
    sample_dir = root / prefix / sample_id
    rows = [
        {
            "cell_barcode": f"{sample_id}_cell1",
            "pass_qc": True,
            "cima_ref_umap_1": 1.0,
            "cima_ref_umap_2": 2.0,
            "cima_cell_type_l1": "B",
            "cima_cell_type_l2": "B naive",
        }
    ]
    if not legacy:
        write_csv(sample_dir / "metadata.csv", rows)
        write_csv(sample_dir / "metadata_qc.csv", rows)
        (sample_dir / "qc_overview.png").write_text("png")
        (sample_dir / "umap_cima_cell_type_l2.png").write_text("png")
        (sample_dir / f"{sample_id}_seurat_qc.rds").write_text("rds")
    write_csv(sample_dir / "qc_summary.csv", [{"n_cells_pass_qc": 1}])
    write_csv(sample_dir / "validation_result.csv", rows)
    write_status(sample_dir / "run_status.json")
    (sample_dir / "umap_cima_cell_type_l1.png").write_text("png")
    (sample_dir / "matrix").mkdir(parents=True)
    (sample_dir / "matrix" / "matrix.mtx").write_text("matrix")
    (sample_dir / "matrix" / "barcodes.tsv").write_text("barcodes")
    (sample_dir / "matrix" / "features.tsv").write_text("features")
    if not legacy:
        (sample_dir / "matrix" / "barcodes.tsv.gz").write_text("barcodes")
        (sample_dir / "matrix" / "features.tsv.gz").write_text("features")
    return sample_dir


def test_discovers_four_product_sources(tmp_path: Path) -> None:
    make_rna_sample(tmp_path, Path("rna/GSE1"), "GSM1")
    make_rna_sample(tmp_path, Path("rna/GSE206284"), "GSM6249236")
    make_rna_sample(tmp_path, Path("rna/GSE1/GSM1/tuning/baseline/GSE1"), "GSM1")
    make_atac_sample(tmp_path, Path("GSE2"), "GSM2", legacy=True)
    make_atac_sample(tmp_path, Path("GSE206284"), "GSM6254833", legacy=True)
    make_atac_sample(tmp_path, Path("co/atac/7555405"), "donorA_Day0")
    make_rna_sample(tmp_path, Path("co/rna/7555405"), "donorA_Day0")

    assert [s.sample_id for s in discover_only_rna(tmp_path)] == ["GSM1"]
    assert [s.sample_id for s in discover_only_atac(tmp_path)] == ["GSM2"]
    assert [s.sample_id for s in discover_co_atac(tmp_path)] == ["donorA_Day0"]
    assert [s.sample_id for s in discover_co_rna(tmp_path)] == ["donorA_Day0"]


def test_rna_defaults_to_harmony_and_atac_defaults_to_bbknn() -> None:
    assert default_integration_method("only_rna") == "harmony"
    assert default_integration_method("co_rna") == "harmony"
    assert default_integration_method("only_atac") == "bbknn"
    assert default_integration_method("co_atac") == "bbknn"


def test_organize_product_writes_traceable_manifests(tmp_path: Path) -> None:
    make_rna_sample(tmp_path, Path("rna/GSE1"), "GSM1")

    status = organize_product(
        output_root=tmp_path,
        product="only_rna",
        copy_mode="symlink",
        force=True,
        skip_figures=True,
        skip_integration=True,
        integration_n_components=30,
        integration_max_umap_fit_cells=100_000,
        integration_clusters=30,
        integration_batch_key="sample_id",
        integration_method="bbknn",
        bbknn_neighbors_within_batch=1,
        bbknn_trim=60,
        leiden_resolution=1.0,
        rna_min_cima_l1_score=0.0,
        include_incomplete=False,
    )

    product_dir = tmp_path / "2.only_rna"
    samples = pd.read_csv(product_dir / "manifests" / "samples.csv")
    cells = pd.read_csv(product_dir / "manifests" / "cells_metadata.csv")
    assert status["n_samples_complete"] == 1
    assert samples.loc[0, "product"] == "only_rna"
    assert cells.loc[0, "global_cell_id"].startswith("only_rna__GSE1__GSM1__")
    assert (product_dir / "samples" / "GSE1" / "GSM1").is_symlink()
    assert (product_dir / "only_rna.h5ad").exists()
    assert status["integration_status"] == "skipped_by_user"


def test_organize_only_atac_uses_current_product_dir(tmp_path: Path) -> None:
    make_atac_sample(tmp_path, Path("GSE2"), "GSM2", legacy=True)

    status = organize_product(
        output_root=tmp_path,
        product="only_atac",
        copy_mode="symlink",
        force=True,
        skip_figures=True,
        skip_integration=True,
        integration_n_components=30,
        integration_max_umap_fit_cells=100_000,
        integration_clusters=30,
        integration_batch_key="sample_id",
        integration_method="bbknn",
        bbknn_neighbors_within_batch=1,
        bbknn_trim=60,
        leiden_resolution=1.0,
        rna_min_cima_l1_score=0.0,
        include_incomplete=False,
    )

    product_dir = tmp_path / "1.only_atac"
    samples = pd.read_csv(product_dir / "manifests" / "samples.csv")
    cells = pd.read_csv(product_dir / "manifests" / "cells_metadata.csv")
    assert status["n_samples_complete"] == 1
    assert samples.loc[0, "product"] == "only_atac"
    assert cells.loc[0, "global_cell_id"].startswith("only_atac__GSE2__GSM2__")
    assert (product_dir / "only_atac.h5ad").exists()


def test_renderer_creates_old_style_panel(tmp_path: Path) -> None:
    metadata = tmp_path / "cells.csv"
    write_csv(
        metadata,
        [
            {"umap_1": 0.0, "umap_2": 0.0, "gse": "GSE1", "sample_id": "S1", "azimuth_cima_l1": "B"},
            {"umap_1": 1.0, "umap_2": 1.0, "gse": "GSE1", "sample_id": "S2", "azimuth_cima_l1": "CD4_T"},
        ],
    )

    result = render_panels(metadata, tmp_path / "figures", "only_rna")

    assert result["status"] == "ok"
    assert (tmp_path / "figures" / "only_rna_cima_l1_panels.png").exists()


def test_renderer_prefers_integrated_umap_and_cluster_panel(tmp_path: Path) -> None:
    metadata = tmp_path / "cells.csv"
    write_csv(
        metadata,
        [
            {
                "integrated_umap_1": 0.0,
                "integrated_umap_2": 0.0,
                "umap_1": 99.0,
                "umap_2": 99.0,
                "gse": "GSE1",
                "sample_id": "S1",
                "integrated_cluster": "0",
                "azimuth_cima_l1": "B",
            },
            {
                "integrated_umap_1": 1.0,
                "integrated_umap_2": 1.0,
                "umap_1": 100.0,
                "umap_2": 100.0,
                "gse": "GSE2",
                "sample_id": "S2",
                "integrated_cluster": "1",
                "azimuth_cima_l1": "CD4_T",
            },
        ],
    )

    result = render_panels(metadata, tmp_path / "figures", "only_rna")

    assert result["coordinate_source"] == "integrated_umap"
    assert (tmp_path / "figures" / "only_rna_integrated_cluster_panels.png").exists()


def test_renderer_prefers_integrated_cima_labels(tmp_path: Path) -> None:
    metadata = tmp_path / "cells.csv"
    write_csv(
        metadata,
        [
            {
                "integrated_umap_1": 0.0,
                "integrated_umap_2": 0.0,
                "azimuth_cima_l1": "B",
                "azimuth_cell_type_l2_raw": "B naive",
                "integrated_cima_l1": "CD4_T",
                "integrated_cima_l2": "CD4_TCM",
            },
            {
                "integrated_umap_1": 1.0,
                "integrated_umap_2": 1.0,
                "azimuth_cima_l1": "B",
                "azimuth_cell_type_l2_raw": "B naive",
                "integrated_cima_l1": "CD8_T",
                "integrated_cima_l2": "CD8_TEM",
            },
        ],
    )

    result = render_panels(metadata, tmp_path / "figures", "only_rna")

    assert result["column:cima_l1"] == "integrated_cima_l1"
    assert result["column:cima_l2"] == "integrated_cima_l2"


def test_metrics_ignore_cells_without_integrated_coordinates() -> None:
    metadata = pd.DataFrame(
        {
            "integrated_umap_1": [0.0, 1.0, pd.NA],
            "integrated_umap_2": [0.0, 1.0, pd.NA],
            "integrated_cluster": ["0", "0", ""],
            "integrated_cima_l1": ["B", "B", "Myeloid"],
            "integrated_cima_l2": ["B naive", "B naive", "Mono"],
            "sample_id": ["S1", "S2", "S3"],
        }
    )
    embedding = pd.DataFrame([[0.0, 0.0], [1.0, 1.0], [99.0, 99.0]]).to_numpy()

    metrics = compute_metrics(metadata, embedding, None, None)

    assert metrics["n_integrated_cells_for_metrics"] == 2
    assert metrics["cluster_cima_l1_purity"] == 1.0
    assert metrics["cluster_cima_l1_column"] == "integrated_cima_l1"
