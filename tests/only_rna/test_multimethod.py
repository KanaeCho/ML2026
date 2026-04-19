from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.only_rna.annotation import annotate_with_all_versions
from scripts.only_rna.azimuth import AzimuthAnnotationResult
from scripts.only_rna.models import (
    AzimuthConfig,
    AnnotationConfig,
    PlottingConfig,
    QcThresholds,
    RunConfig,
)
from scripts.only_rna.outputs import write_sample_outputs
from scripts.only_rna.config import merge_cli_overrides
from scripts.only_rna.plotting import save_annotation_method_comparison_umap


def _make_run_config() -> RunConfig:
    return RunConfig(
        qc=QcThresholds(
            min_counts=500,
            min_genes=300,
            max_pct_mt=20.0,
            max_pct_ribo=60.0,
        ),
        plotting=PlottingConfig(
            umap_width=4.0,
            umap_height=3.0,
            dpi=50,
            point_size=20.0,
            legend_fontsize=7.0,
            legend_title_fontsize=8.0,
        ),
        annotation=AnnotationConfig(
            methods=["cima", "azimuth", "cell_typist", "singler", "scanvi"]
        ),
    )


def _write_cima_reference_assets(base_dir: Path) -> Path:
    cima_dir = base_dir / "cima"
    cima_dir.mkdir(parents=True, exist_ok=True)

    pd.DataFrame(
        {
            "feature_id": ["GeneA", "GeneB"],
            "gene_mean": [0.0, 0.0],
            "gene_std": [1.0, 1.0],
            "pc_dim_1": [1.0, 0.0],
            "pc_dim_2": [0.0, 1.0],
        }
    ).to_csv(
        cima_dir / "cima_rna_reference_pca_features.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )

    pd.DataFrame(
        [
            {"cell_type_l1": "L1_A", "pc_dim_1": 1.0, "pc_dim_2": 0.0},
            {"cell_type_l1": "L1_B", "pc_dim_1": 0.0, "pc_dim_2": 1.0},
        ]
    ).to_csv(cima_dir / "cima_rna_reference_l1_centroids.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {"cell_type_l2": "L2_A1", "pc_dim_1": 1.0, "pc_dim_2": 0.0},
            {"cell_type_l2": "L2_B1", "pc_dim_1": 0.0, "pc_dim_2": 1.0},
        ]
    ).to_csv(cima_dir / "cima_rna_reference_l2_centroids.tsv", sep="\t", index=False)

    pd.DataFrame(
        [
            {"cell_type_l1": "L1_A", "cell_type_l2": "L2_A1"},
            {"cell_type_l1": "L1_B", "cell_type_l2": "L2_B1"},
        ]
    ).to_csv(cima_dir / "cima_rna_celltype_hierarchy.csv", index=False)

    return base_dir


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
            {
                "pass_qc": [True, False, True],
                "umap_1": [0.0, np.nan, 2.0],
                "umap_2": [1.0, np.nan, 3.0],
            },
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame({"feature_id": ["GeneA", "GeneB"]}, index=["g1", "g2"]),
    )


def test_annotate_with_all_versions_initializes_expected_method_columns(tmp_path: Path):
    reference_dir = _write_cima_reference_assets(tmp_path)
    adata = _make_annotation_adata()

    out = annotate_with_all_versions(
        adata,
        reference_dir=reference_dir,
        methods=["cima", "azimuth", "cell_typist", "singler", "scanvi"],
    )

    expected_string_columns = [
        "cima_l1",
        "cima_l2",
        "cima_l1_masked",
        "azimuth_cell_type",
        "celltypist_cell_type",
        "singler_cell_type",
        "scanvi_cell_type",
    ]
    expected_float_columns = [
        "cima_l1_score",
        "cima_l1_score_margin",
        "cima_l2_score",
        "cima_l2_score_margin",
        "azimuth_score",
        "azimuth_score_margin",
        "celltypist_score",
        "celltypist_score_margin",
        "singler_score",
        "singler_score_margin",
        "scanvi_score",
        "scanvi_score_margin",
    ]
    expected_bool_columns = [
        "cima_l1_low_confidence",
        "azimuth_low_confidence",
        "celltypist_low_confidence",
        "singler_low_confidence",
        "scanvi_low_confidence",
    ]

    for column in expected_string_columns:
        assert str(out.obs[column].dtype) == "string"
    for column in expected_float_columns:
        assert pd.api.types.is_float_dtype(out.obs[column])
    for column in expected_bool_columns:
        assert str(out.obs[column].dtype) == "boolean"

    assert out.obs.loc["cell-1", "cima_l1"] == "L1_A"
    assert pd.isna(out.obs.loc["cell-2", "azimuth_cell_type"])
    assert pd.isna(out.obs.loc["cell-2", "scanvi_low_confidence"])


def test_annotate_with_all_versions_uses_shared_azimuth_result_and_records_status(
    tmp_path: Path, monkeypatch
):
    reference_dir = _write_cima_reference_assets(tmp_path)
    adata = _make_annotation_adata()

    def _fake_run_azimuth_annotation(
        _adata,
        *,
        config: AzimuthConfig,
        annotation_level: str = "l1",
        max_cells: int | None = None,
    ) -> AzimuthAnnotationResult:
        assert config.enabled is True
        assert annotation_level == "l1"
        assert max_cells is None
        return AzimuthAnnotationResult(
            labels=pd.Series(
                ["CD4 T", "B"], index=["cell-1", "cell-3"], dtype="string"
            ),
            status="ok",
            detail="pbmcref",
        )

    monkeypatch.setattr(
        "scripts.only_rna.annotation.run_azimuth_annotation",
        _fake_run_azimuth_annotation,
    )

    out = annotate_with_all_versions(
        adata,
        reference_dir=reference_dir,
        methods=["cima", "azimuth"],
        azimuth_config=AzimuthConfig(enabled=True),
    )

    assert out.obs.loc["cell-1", "azimuth_cell_type"] == "CD4 T"
    assert out.obs.loc["cell-3", "azimuth_cell_type"] == "B"
    assert pd.isna(out.obs.loc["cell-2", "azimuth_cell_type"])
    assert out.uns["annotation_method_status"]["azimuth"] == {
        "status": "ok",
        "detail": "pbmcref",
    }


def test_annotate_with_all_versions_keeps_na_azimuth_labels_and_records_error_status(
    tmp_path: Path, monkeypatch
):
    reference_dir = _write_cima_reference_assets(tmp_path)
    adata = _make_annotation_adata()

    monkeypatch.setattr(
        "scripts.only_rna.annotation.run_azimuth_annotation",
        lambda *_args, **_kwargs: AzimuthAnnotationResult(
            labels=None,
            status="error",
            detail="Azimuth crashed",
        ),
    )

    out = annotate_with_all_versions(
        adata,
        reference_dir=reference_dir,
        methods=["cima", "azimuth"],
        azimuth_config=AzimuthConfig(enabled=True),
    )

    assert pd.isna(out.obs.loc["cell-1", "azimuth_cell_type"])
    assert pd.isna(out.obs.loc["cell-3", "azimuth_cell_type"])
    assert pd.isna(out.obs.loc["cell-1", "azimuth_score"])
    assert pd.isna(out.obs.loc["cell-3", "azimuth_score_margin"])
    assert pd.isna(out.obs.loc["cell-1", "azimuth_low_confidence"])
    assert out.uns["annotation_method_status"]["azimuth"] == {
        "status": "error",
        "detail": "Azimuth crashed",
    }


def test_write_sample_outputs_emits_four_method_comparison_umap(tmp_path: Path):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(4)),
        obs=pd.DataFrame(
            {
                "pass_qc": [True, True, True, True],
                "cluster": ["0", "0", "1", "1"],
                "umap_1": [0.0, 1.0, 2.0, 3.0],
                "umap_2": [0.0, 1.0, 0.5, 1.5],
                "cima_l1": ["L1_A", "L1_A", "L1_B", "L1_B"],
                "cima_l2": ["L2_A1", "L2_A1", "L2_B1", "L2_B1"],
                "cima_l1_masked": ["L1_A", "L1_A", "L1_B", "L1_B"],
                "azimuth_cell_type": ["T", "T", "B", "B"],
                "celltypist_cell_type": ["T", "T", "B", "B"],
                "singler_cell_type": ["T", "T", "B", "B"],
                "scanvi_cell_type": ["T", "T", "B", "B"],
            },
            index=[f"cell-{i}" for i in range(4)],
        ),
        var=pd.DataFrame(
            {"feature_id": ["ENSG1", "ENSG2", "ENSG3", "ENSG4"]},
            index=["G1", "G2", "G3", "G4"],
        ),
    )

    sample_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE192391",
        sample_id="GSM000001",
        config=config,
    )

    comparison_path = sample_dir / "umap_rna_annotation_method_compare.png"
    assert comparison_path.exists()

    validation = pd.read_csv(sample_dir / "validation_result.csv")
    assert validation.loc[
        validation["check_name"]
        == "output_presence:umap_rna_annotation_method_compare.png",
        "passed",
    ].tolist() == [True]


def test_merge_cli_overrides_preserves_annotation_methods():
    base = _make_run_config()

    merged = merge_cli_overrides(base, qc__min_genes=350)

    assert merged.qc.min_genes == 350
    assert merged.annotation is not None
    assert merged.annotation.methods == [
        "cima",
        "azimuth",
        "cell_typist",
        "singler",
        "scanvi",
    ]


def test_save_annotation_method_comparison_umap_handles_method_specific_columns(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(4)),
        obs=pd.DataFrame(
            {
                "umap_1": [0.0, 1.0, 2.0, 3.0],
                "umap_2": [0.0, 1.0, 0.5, 1.5],
                "azimuth_cell_type": ["T", "T", "B", "B"],
                "celltypist_cell_type": ["T", "T", "B", "B"],
                "singler_cell_type": ["T", "T", "B", "B"],
                "scanvi_cell_type": ["T", "T", "B", "B"],
            },
            index=[f"cell-{i}" for i in range(4)],
        ),
    )

    output_path = tmp_path / "comparison.png"
    save_annotation_method_comparison_umap(
        adata,
        output_path=output_path,
        title="Comparison",
        config=config,
    )

    assert output_path.exists()
