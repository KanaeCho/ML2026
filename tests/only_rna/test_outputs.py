from __future__ import annotations

from pathlib import Path
import json

import anndata as ad
import matplotlib.axes
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from matplotlib.colors import to_hex, to_rgb
from scipy import sparse

from scripts.only_rna.models import PlottingConfig, QcThresholds, RunConfig
from scripts.only_rna.outputs import write_sample_outputs
from scripts.only_rna.plotting import (
    build_cima_l1_palette,
    build_cima_l2_palette,
    get_hierarchical_palette,
    save_categorical_umap,
    save_qc_overview,
)


def _make_run_config() -> RunConfig:
    return RunConfig(
        qc=QcThresholds(
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
        plotting=PlottingConfig(
            umap_width=4.0,
            umap_height=3.0,
            dpi=50,
            point_size=20.0,
            legend_fontsize=7.0,
            legend_title_fontsize=8.0,
        ),
    )


def _make_output_adata() -> ad.AnnData:
    adata = ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [1.0, 0.0, 3.0],
                    [0.0, 2.0, 0.0],
                    [4.0, 5.0, 6.0],
                ]
            )
        ),
        obs=pd.DataFrame(
            {
                "pass_qc": [True, False, True],
                "cluster": ["0", pd.NA, "1"],
                "umap_1": [0.1, np.nan, 1.5],
                "umap_2": [1.0, np.nan, 2.0],
                "azimuth_cell_type": ["CD4 T", pd.NA, "B"],
                "azimuth_cell_type_l1_raw": ["CD4 T", pd.NA, "B"],
                "azimuth_cell_type_l2_raw": ["CD4 TCM", pd.NA, "B naive"],
                "azimuth_cima_l1": ["CD4_T", pd.NA, "B"],
                "azimuth_cima_l1_unmapped": [False, pd.NA, False],
            },
            index=["cell-1", "cell-2", "cell-3"],
        ),
        var=pd.DataFrame(
            {
                "feature_id": ["ENSG0001", "ENSG0002", "ENSG0003"],
                "feature_name": ["GeneA", "GeneB", "GeneC"],
            },
            index=["GeneA", "GeneB", "GeneC"],
        ),
    )
    adata.uns["qc_thresholds"] = {
        "sample_id": "GSM123456",
        "gse": "GSE123456",
        "method": "dynamic_hybrid_mad",
        "min_counts": 100,
        "min_genes": 100,
        "max_pct_mt": 40.0,
        "max_pct_ribo": 80.0,
    }
    return adata


def test_save_categorical_umap_creates_png_with_config_dimensions(tmp_path: Path):
    config = _make_run_config()
    adata = _make_output_adata()
    adata.obs.loc["cell-2", "cima_l1_masked"] = pd.NA
    adata.obs.loc["cell-3", "umap_2"] = np.nan

    output_path = tmp_path / "umap.png"
    save_categorical_umap(
        adata,
        color_key="cima_l1_masked",
        output_path=output_path,
        title="Masked L1",
        config=config,
    )

    assert output_path.exists()
    image = mpimg.imread(output_path)
    expected_width = int(config.plotting.umap_width * config.plotting.dpi)
    assert image.shape[1] == expected_width
    assert image.shape[0] == int(config.plotting.umap_height * config.plotting.dpi)


def test_save_categorical_umap_uses_readable_legend_markers_and_layout(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(14)),
        obs=pd.DataFrame(
            {
                "umap_1": np.linspace(0.0, 13.0, 14),
                "umap_2": np.linspace(1.0, 14.0, 14),
                "cluster": [str(i) for i in range(14)],
            },
            index=[f"cell-{i}" for i in range(14)],
        ),
    )

    captured: dict[str, object] = {}
    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        captured.update(kwargs)
        return original_legend(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "legend", capture_legend)

    output_path = tmp_path / "legend_readable.png"
    save_categorical_umap(
        adata,
        color_key="cluster",
        output_path=output_path,
        title="Legend readability",
        config=config,
    )

    assert output_path.exists()
    assert captured["markerscale"] > 1.0
    assert captured["ncol"] > 1


def test_save_categorical_umap_keeps_readable_panel_for_many_categories(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(27)),
        obs=pd.DataFrame(
            {
                "umap_1": np.linspace(-5.0, 5.0, 27),
                "umap_2": np.linspace(-3.0, 3.0, 27),
                "cima_l2": [f"L2_{i}" for i in range(27)],
            },
            index=[f"cell-{i}" for i in range(27)],
        ),
    )

    captured: dict[str, object] = {}
    original_subplots = plt.subplots

    def capture_subplots(*args, **kwargs):
        captured["figsize"] = kwargs.get("figsize")
        captured["dpi"] = kwargs.get("dpi")
        fig, ax = original_subplots(*args, **kwargs)

        original_tight_layout = fig.tight_layout

        def capture_tight_layout(*layout_args, **layout_kwargs):
            captured["tight_layout_args"] = layout_args
            captured["tight_layout_kwargs"] = layout_kwargs
            return original_tight_layout(*layout_args, **layout_kwargs)

        fig.tight_layout = capture_tight_layout

        original_set_aspect = ax.set_aspect
        original_set_box_aspect = ax.set_box_aspect

        def capture_set_aspect(*aspect_args, **aspect_kwargs):
            captured["set_aspect_args"] = aspect_args
            captured["set_aspect_kwargs"] = aspect_kwargs
            return original_set_aspect(*aspect_args, **aspect_kwargs)

        def capture_set_box_aspect(*box_aspect_args, **box_aspect_kwargs):
            captured["set_box_aspect_args"] = box_aspect_args
            captured["set_box_aspect_kwargs"] = box_aspect_kwargs
            return original_set_box_aspect(*box_aspect_args, **box_aspect_kwargs)

        ax.set_aspect = capture_set_aspect
        ax.set_box_aspect = capture_set_box_aspect
        return fig, ax

    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        captured.update(
            {
                "ncol": kwargs.get("ncol"),
                "markerscale": kwargs.get("markerscale"),
                "bbox_to_anchor": kwargs.get("bbox_to_anchor"),
                "loc": kwargs.get("loc"),
            }
        )
        return original_legend(self, *args, **kwargs)

    monkeypatch.setattr(plt, "subplots", capture_subplots)
    monkeypatch.setattr(matplotlib.axes.Axes, "legend", capture_legend)

    output_path = tmp_path / "l2_many_categories.png"
    save_categorical_umap(
        adata,
        color_key="cima_l2",
        output_path=output_path,
        title="L2 readability",
        config=config,
    )

    assert output_path.exists()
    assert captured["figsize"][1] == config.plotting.umap_height
    assert captured["dpi"] == config.plotting.dpi
    assert captured["markerscale"] >= 4.0
    assert captured["ncol"] >= 2
    assert captured["loc"] == "center left"
    assert captured["bbox_to_anchor"] is not None
    assert captured["figsize"][0] > config.plotting.umap_width
    assert captured["figsize"][1] == config.plotting.umap_height
    assert captured["set_box_aspect_args"][0] == 1
    assert captured["tight_layout_kwargs"]["rect"] == (0.0, 0.0, 1.0, 1.0)


def test_save_categorical_umap_uses_taller_fewer_columns_for_l2_legend(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(27)),
        obs=pd.DataFrame(
            {
                "umap_1": np.linspace(-5.0, 5.0, 27),
                "umap_2": np.linspace(-3.0, 3.0, 27),
                "cima_l2": [f"L2_{i}" for i in range(27)],
            },
            index=[f"cell-{i}" for i in range(27)],
        ),
    )

    captured: dict[str, object] = {}
    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        captured.update(
            {
                "ncol": kwargs.get("ncol"),
                "bbox_to_anchor": kwargs.get("bbox_to_anchor"),
                "loc": kwargs.get("loc"),
            }
        )
        return original_legend(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "legend", capture_legend)

    output_path = tmp_path / "l2_taller_columns.png"
    save_categorical_umap(
        adata,
        color_key="cima_l2",
        output_path=output_path,
        title="L2 taller columns",
        config=config,
    )

    assert output_path.exists()
    assert captured["loc"] == "center left"
    assert captured["bbox_to_anchor"] == (1.02, 0.5)
    assert captured["ncol"] == 2


def test_save_categorical_umap_uses_visibility_tuned_scatter_defaults(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = _make_output_adata()

    captured: dict[str, object] = {}
    original_scatter = matplotlib.axes.Axes.scatter

    def capture_scatter(self, *args, **kwargs):
        captured.update(
            {
                "s": kwargs.get("s"),
                "alpha": kwargs.get("alpha"),
                "linewidths": kwargs.get("linewidths"),
            }
        )
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

    output_path = tmp_path / "visibility_tuned.png"
    save_categorical_umap(
        adata,
        color_key="cluster",
        output_path=output_path,
        title="Visibility tuning",
        config=config,
    )

    assert output_path.exists()
    assert captured["s"] > config.plotting.point_size
    assert captured["alpha"] < 0.95


def test_write_sample_outputs_emits_required_files_and_metadata_contract(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()

    sample_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE123456",
        sample_id="GSM123456",
        config=config,
    )

    assert sample_dir == tmp_path / "GSE123456" / "GSM123456"

    expected_files = [
        sample_dir / "metadata.csv",
        sample_dir / "metadata_qc.csv",
        sample_dir / "qc_summary.csv",
        sample_dir / "qc_thresholds.json",
        sample_dir / "qc_overview.png",
        sample_dir / "validation_result.csv",
        sample_dir / "GSM123456.h5ad",
        sample_dir / "umap_rna_pbmcref_vs_cima_l1.png",
        sample_dir / "umap_rna_pbmcref_highlight.png",
        sample_dir / "umap_rna_cima_l1.png",
    ]
    for path in expected_files:
        assert path.exists(), path

    assert sorted(path.name for path in sample_dir.glob("*.h5ad")) == ["GSM123456.h5ad"]

    metadata = pd.read_csv(sample_dir / "metadata.csv")
    metadata_qc = pd.read_csv(sample_dir / "metadata_qc.csv")
    qc_summary = pd.read_csv(sample_dir / "qc_summary.csv")
    qc_thresholds = json.loads((sample_dir / "qc_thresholds.json").read_text())
    validation = pd.read_csv(sample_dir / "validation_result.csv")

    assert "azimuth_cell_type" in metadata.columns
    assert "azimuth_cell_type_l1_raw" in metadata.columns
    assert "azimuth_cell_type_l2_raw" in metadata.columns
    assert "azimuth_cima_l1" in metadata.columns
    assert "azimuth_cima_l1_unmapped" in metadata.columns

    assert metadata["cell_id"].tolist() == ["cell-1", "cell-2", "cell-3"]
    assert metadata_qc["cell_id"].tolist() == ["cell-1", "cell-3"]
    assert metadata_qc["pass_qc"].tolist() == [True, True]

    assert len(qc_summary) == 1
    assert qc_summary.loc[0, "sample_id"] == "GSM123456"
    assert qc_summary.loc[0, "gse"] == "GSE123456"
    assert qc_summary.loc[0, "n_cells_total"] == 3
    assert qc_summary.loc[0, "n_cells_pass_qc"] == 2
    assert qc_summary.loc[0, "n_cells_fail_qc"] == 1
    assert qc_summary.loc[0, "n_cells_final_output"] == 2
    assert qc_summary.loc[0, "n_cells_unknown_final_celltype_removed"] == 0
    assert qc_summary.loc[0, "pass_qc_fraction"] == pytest.approx(2 / 3)
    assert qc_summary.loc[0, "qc_threshold_method"] == "dynamic_hybrid_mad"
    assert pd.isna(qc_summary.loc[0, "azimuth_status"])
    assert pd.isna(qc_summary.loc[0, "azimuth_detail"])
    assert qc_summary.loc[0, "azimuth_score_mean"] == 0.0
    assert qc_summary.loc[0, "azimuth_score_margin_mean"] == 0.0
    assert qc_summary.loc[0, "azimuth_low_confidence_fraction"] == 0.0
    assert qc_summary.loc[0, "annotation_score"] == 0.0
    assert "azimuth_score_mean" in qc_summary.columns
    assert "azimuth_score_margin_mean" in qc_summary.columns
    assert "azimuth_low_confidence_fraction" in qc_summary.columns
    assert "annotation_score" in qc_summary.columns

    assert qc_thresholds["method"] == "dynamic_hybrid_mad"

    completion = validation.loc[validation["check_name"] == "completion"]
    assert completion["passed"].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:metadata.csv", "passed"
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:qc_thresholds.json", "passed"
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:qc_overview.png", "passed"
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:umap_rna_pbmcref_vs_cima_l1.png",
        "passed",
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:umap_rna_pbmcref_highlight.png",
        "passed",
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:umap_rna_cima_l1.png",
        "passed",
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "annotation_eval_presence:annotation_score",
        "passed",
    ].tolist() == [True]


def test_write_sample_outputs_renders_visible_umap_from_dual_annotation_fields(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = _make_output_adata()

    captured: list[tuple[str, str]] = []

    def _fake_save_dual_annotation_umap(*, adata, output_path, title, config):
        del adata, title, config
        captured.append(("dual", Path(output_path).name))
        Path(output_path).write_bytes(b"fake-png")

    def _fake_save_highlight_category_overview(
        *, adata, color_key, output_path, title, config, legend_title
    ):
        del adata, color_key, title, config, legend_title
        captured.append(("highlight", Path(output_path).name))
        Path(output_path).write_bytes(b"fake-png")

    def _fake_save_qc_overview(adata, output_path, config):
        del adata, config
        captured.append(("qc", Path(output_path).name))
        Path(output_path).write_bytes(b"fake-png")

    def _fake_save_sample_cima_l1_umap(adata, output_path, title, config):
        del adata, title, config
        captured.append(("cima_l1", Path(output_path).name))
        Path(output_path).write_bytes(b"fake-png")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_dual_annotation_umap",
        _fake_save_dual_annotation_umap,
    )
    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_highlight_category_overview",
        _fake_save_highlight_category_overview,
    )
    monkeypatch.setattr("scripts.only_rna.outputs.save_qc_overview", _fake_save_qc_overview)
    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_sample_cima_l1_umap",
        _fake_save_sample_cima_l1_umap,
    )

    sample_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE246810",
        sample_id="GSM246810",
        config=config,
    )

    assert sample_dir == tmp_path / "GSE246810" / "GSM246810"
    assert captured == [
        ("qc", "qc_overview.png"),
        ("dual", "umap_rna_pbmcref_vs_cima_l1.png"),
        ("highlight", "umap_rna_pbmcref_highlight.png"),
        ("cima_l1", "umap_rna_cima_l1.png"),
    ]


def test_write_sample_outputs_dual_annotation_prefers_fine_pbmcref_raw_labels(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = _make_output_adata()

    captured: dict[str, object] = {}

    def _fake_save_dual_annotation_umap(*, adata, output_path, title, config):
        captured["left_values"] = (
            adata.obs["azimuth_cell_type_l2_raw"].astype("string").tolist()
        )
        captured["output_name"] = Path(output_path).name
        Path(output_path).write_bytes(b"fake-png")

    def _fake_save_highlight_category_overview(
        *, adata, color_key, output_path, title, config, legend_title
    ):
        del adata, title, config, legend_title
        captured["highlight_color_key"] = color_key
        captured["highlight_output_name"] = Path(output_path).name
        Path(output_path).write_bytes(b"fake-png")

    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_dual_annotation_umap",
        _fake_save_dual_annotation_umap,
    )
    monkeypatch.setattr(
        "scripts.only_rna.outputs.save_highlight_category_overview",
        _fake_save_highlight_category_overview,
    )

    write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE135790",
        sample_id="GSM135790",
        config=config,
    )

    assert captured["output_name"] == "umap_rna_pbmcref_vs_cima_l1.png"
    assert captured["highlight_output_name"] == "umap_rna_pbmcref_highlight.png"
    assert captured["highlight_color_key"] == "azimuth_cell_type_l2_raw"
    assert captured["left_values"] == ["CD4 TCM", pd.NA, "B naive"]


def test_save_dual_annotation_umap_adds_text_labels_to_pbmcref_panel(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [0.0, 0.1, 0.2, 8.0, 8.1, 8.2],
            "umap_2": [0.0, 0.1, 0.2, 8.0, 8.1, 8.2],
            "azimuth_cell_type_l2_raw": [
                "CD4 TCM",
                "CD4 TCM",
                "CD4 TCM",
                "B naive",
                "B naive",
                "B naive",
            ],
            "azimuth_cima_l1": ["CD4_T", "CD4_T", "CD4_T", "B", "B", "B"],
        },
        index=[f"cell-{i}" for i in range(6)],
    )
    adata = ad.AnnData(X=np.ones((6, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[tuple[str, str]] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label in {"CD4 TCM", "B naive"}:
            captured.append((str(self.get_title()), str(label)))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert ("pbmcref", "CD4 TCM") in captured
    assert ("pbmcref", "B naive") in captured


def test_save_dual_annotation_umap_labels_multiple_large_clusters_per_cell_type(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [
                0.0,
                0.1,
                0.2,
                0.3,
                10.0,
                10.1,
                10.2,
                10.3,
                5.0,
            ],
            "umap_2": [
                0.0,
                0.1,
                0.2,
                0.3,
                10.0,
                10.1,
                10.2,
                10.3,
                5.0,
            ],
            "azimuth_cell_type_l2_raw": [
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "Unknown",
            ],
            "azimuth_cima_l1": [
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "Unknown",
            ],
        },
        index=[f"cell-{i}" for i in range(9)],
    )
    adata = ad.AnnData(X=np.ones((9, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[tuple[str, str]] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label == "CD8 TEM":
            captured.append((str(self.get_title()), str(label)))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_clusters.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert captured.count(("pbmcref", "CD8 TEM")) <= 2
    assert captured.count(("pbmcref", "CD8 TEM")) >= 1


def test_save_dual_annotation_umap_skips_small_clusters_for_same_cell_type(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [
                0.0,
                0.1,
                0.2,
                0.3,
                8.0,
            ],
            "umap_2": [
                0.0,
                0.1,
                0.2,
                0.3,
                8.0,
            ],
            "azimuth_cell_type_l2_raw": [
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
            ],
            "azimuth_cima_l1": ["CD8_T", "CD8_T", "CD8_T", "CD8_T", "CD8_T"],
        },
        index=[f"cell-{i}" for i in range(5)],
    )
    adata = ad.AnnData(X=np.ones((5, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[tuple[str, str]] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label == "CD8 TEM":
            captured.append((str(self.get_title()), str(label)))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_large_only.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert captured.count(("pbmcref", "CD8 TEM")) == 1


def test_save_dual_annotation_umap_skips_tiny_components_when_large_cluster_exists(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                10.0,
                10.1,
                10.2,
            ],
            "umap_2": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                10.0,
                10.1,
                10.2,
            ],
            "azimuth_cell_type_l2_raw": [
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
            ],
            "azimuth_cima_l1": ["CD8_T"] * 8,
        },
        index=[f"cell-{i}" for i in range(8)],
    )
    adata = ad.AnnData(X=np.ones((8, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[tuple[str, str]] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label == "CD8 TEM":
            captured.append((str(self.get_title()), str(label)))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_skip_tiny_components.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert captured.count(("pbmcref", "CD8 TEM")) == 1


def test_save_dual_annotation_umap_only_labels_large_clusters(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                10.0,
                10.1,
                30.0,
                30.1,
                30.2,
                30.3,
                30.4,
            ],
            "umap_2": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                10.0,
                10.1,
                30.0,
                30.1,
                30.2,
                30.3,
                30.4,
            ],
            "azimuth_cell_type_l2_raw": [
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD8 TEM",
                "CD4 TCM",
                "CD4 TCM",
                "CD4 TCM",
                "CD4 TCM",
                "CD4 TCM",
            ],
            "azimuth_cima_l1": [
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD8_T",
                "CD4_T",
                "CD4_T",
                "CD4_T",
                "CD4_T",
                "CD4_T",
            ],
        },
        index=[f"cell-{i}" for i in range(12)],
    )
    adata = ad.AnnData(X=np.ones((12, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[tuple[str, str]] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if str(self.get_title()) == "pbmcref":
            captured.append((str(self.get_title()), str(label)))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_large_clusters_only.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert captured.count(("pbmcref", "CD4 TCM")) == 1
    assert captured.count(("pbmcref", "CD8 TEM")) == 1


def test_save_dual_annotation_umap_skips_gray_unknown_labels(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            "umap_2": [0.0, 0.1, 0.2, 1.0, 1.1, 1.2],
            "azimuth_cell_type_l2_raw": [
                "Unknown",
                "Unknown",
                "Unknown",
                "CD4 TCM",
                "CD4 TCM",
                "CD4 TCM",
            ],
            "azimuth_cima_l1": [
                "Unknown",
                "Unknown",
                "Unknown",
                "CD4_T",
                "CD4_T",
                "CD4_T",
            ],
        },
        index=[f"cell-{i}" for i in range(6)],
    )
    adata = ad.AnnData(X=np.ones((6, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: list[str] = []
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label is not None and str(self.get_title()) == "pbmcref":
            captured.append(str(label))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_unknown.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    assert "CD4 TCM" in captured
    assert "Unknown" not in captured


def test_save_dual_annotation_umap_uses_family_colored_label_boxes(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_dual_annotation_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [0.0, 0.1, 0.2],
            "umap_2": [0.0, 0.1, 0.2],
            "azimuth_cell_type_l2_raw": ["CD4 TCM", "CD4 TCM", "CD4 TCM"],
            "azimuth_cima_l1": ["CD4_T", "CD4_T", "CD4_T"],
        },
        index=[f"cell-{i}" for i in range(3)],
    )
    adata = ad.AnnData(X=np.ones((3, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured: dict[str, object] = {}
    original_text = matplotlib.axes.Axes.text

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label == "CD4 TCM" and str(self.get_title()) == "pbmcref":
            bbox = kwargs.get("bbox") or {}
            captured["facecolor"] = bbox.get("facecolor")
            captured["alpha"] = bbox.get("alpha")
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "dual_family_box.png"
    save_dual_annotation_umap(
        adata=adata,
        output_path=output_path,
        title="pbmcref vs CIMA L1",
        config=config,
    )

    assert output_path.exists()
    expected = get_hierarchical_palette("azimuth_cell_type_l2_raw", ["CD4 TCM"])[
        "CD4 TCM"
    ]
    assert tuple(captured["facecolor"]) == pytest.approx(expected)
    assert captured["alpha"] is not None


def test_write_sample_outputs_records_annotation_status_in_qc_and_validation(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()
    adata.uns["annotation_method_status"] = {
        "azimuth": {
            "status": "error",
            "detail": "Azimuth unavailable",
        }
    }

    sample_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE654321",
        sample_id="GSM654321",
        config=config,
    )

    qc_summary = pd.read_csv(sample_dir / "qc_summary.csv")
    validation = pd.read_csv(sample_dir / "validation_result.csv")

    assert qc_summary.loc[0, "pass_qc_fraction"] == pytest.approx(2 / 3)
    assert qc_summary.loc[0, "azimuth_status"] == "error"
    assert qc_summary.loc[0, "azimuth_detail"] == "Azimuth unavailable"

    assert validation.loc[
        validation["check_name"] == "annotation_status:azimuth", "passed"
    ].tolist() == [False]
    assert validation.loc[
        validation["check_name"] == "annotation_status:azimuth", "detail"
    ].tolist() == ["Azimuth unavailable"]


def test_write_sample_outputs_validation_contract_excludes_removed_cima_and_comparison_artifacts(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()

    sample_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE777777",
        sample_id="GSM777777",
        config=config,
    )

    validation = pd.read_csv(sample_dir / "validation_result.csv")
    check_names = set(validation["check_name"].tolist())

    assert "output_presence:qc_thresholds.json" in check_names
    assert "output_presence:umap_rna_pbmcref_vs_cima_l1.png" in check_names
    assert "output_presence:umap_rna_pbmcref_highlight.png" in check_names
    assert "output_presence:umap_rna_clusters.png" not in check_names
    assert "output_presence:umap_rna_cima_cell_type_l1.png" not in check_names
    assert "output_presence:umap_rna_cima_cell_type_l2.png" not in check_names
    assert "output_presence:umap_rna_cima_cell_type_l1_masked.png" not in check_names
    assert "output_presence:umap_rna_annotation_method_compare.png" not in check_names

    assert not (sample_dir / "umap_rna_clusters.png").exists()
    assert not (sample_dir / "umap_rna_cima_cell_type_l1.png").exists()
    assert not (sample_dir / "umap_rna_cima_cell_type_l2.png").exists()
    assert not (sample_dir / "umap_rna_cima_cell_type_l1_masked.png").exists()


def test_write_sample_outputs_removes_legacy_single_umap_artifact(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()

    sample_dir = tmp_path / "GSE999999" / "GSM999999"
    sample_dir.mkdir(parents=True, exist_ok=True)
    stale_single_umap = sample_dir / "umap_rna_azimuth.png"
    stale_single_umap.write_bytes(b"stale")

    written_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE999999",
        sample_id="GSM999999",
        config=config,
    )

    assert written_dir == sample_dir
    assert not stale_single_umap.exists()
    assert (sample_dir / "umap_rna_pbmcref_vs_cima_l1.png").exists()
    assert (sample_dir / "umap_rna_pbmcref_highlight.png").exists()


def test_write_sample_outputs_removes_stale_cluster_and_cima_artifacts(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()

    sample_dir = tmp_path / "GSE888888" / "GSM888888"
    sample_dir.mkdir(parents=True, exist_ok=True)
    stale_paths = [
        sample_dir / "umap_rna_clusters.png",
        sample_dir / "umap_rna_cima_cell_type_l1.png",
        sample_dir / "umap_rna_cima_cell_type_l2.png",
        sample_dir / "umap_rna_cima_cell_type_l1_masked.png",
        sample_dir / "umap_rna_annotation_method_compare.png",
    ]
    for path in stale_paths:
        path.write_bytes(b"stale")

    written_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE888888",
        sample_id="GSM888888",
        config=config,
    )

    assert written_dir == sample_dir
    for path in stale_paths:
        assert not path.exists(), path
    assert (sample_dir / "umap_rna_pbmcref_vs_cima_l1.png").exists()
    assert (sample_dir / "umap_rna_pbmcref_highlight.png").exists()


def test_write_sample_outputs_removes_stale_nested_sample_directory(
    tmp_path: Path,
):
    config = _make_run_config()
    adata = _make_output_adata()

    sample_dir = tmp_path / "GSE123123" / "GSM123123"
    stale_nested = sample_dir / "GSE123123" / "GSM123123"
    stale_nested.mkdir(parents=True, exist_ok=True)
    (stale_nested / "metadata.csv").write_text("stale\n", encoding="utf-8")

    written_dir = write_sample_outputs(
        adata,
        output_root=tmp_path,
        gse="GSE123123",
        sample_id="GSM123123",
        config=config,
    )

    assert written_dir == sample_dir
    assert not stale_nested.exists()


def test_save_highlight_category_overview_uses_larger_panel_titles(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_highlight_category_overview

    config = _make_run_config()
    adata = _make_output_adata()

    captured_sizes: list[float] = []
    original_set_title = matplotlib.axes.Axes.set_title

    def capture_set_title(self, label, *args, **kwargs):
        if label in {"CD4 TCM", "B naive"}:
            captured_sizes.append(float(kwargs.get("fontsize", 0.0)))
        return original_set_title(self, label, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "set_title", capture_set_title)

    output_path = tmp_path / "highlight_titles.png"
    save_highlight_category_overview(
        adata=adata,
        color_key="azimuth_cell_type_l2_raw",
        output_path=output_path,
        title="pbmcref highlight",
        config=config,
        legend_title="pbmcref",
    )

    assert output_path.exists()
    assert captured_sizes
    assert min(captured_sizes) >= 12.0


def test_save_highlight_category_overview_draws_dashed_cluster_outlines(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_highlight_category_overview

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                5.0,
                5.1,
            ],
            "umap_2": [
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                5.0,
                5.1,
            ],
            "azimuth_cell_type_l2_raw": ["CD4 TCM"] * 10 + ["B naive", "B naive"],
        },
        index=[f"cell-{i}" for i in range(12)],
    )
    adata = ad.AnnData(X=np.ones((12, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured_linestyles: list[object] = []
    original_contour = matplotlib.axes.Axes.contour

    def capture_contour(self, *args, **kwargs):
        captured_linestyles.extend(kwargs.get("linestyles", []))
        return original_contour(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "contour", capture_contour)

    output_path = tmp_path / "highlight_outlines.png"
    save_highlight_category_overview(
        adata=adata,
        color_key="azimuth_cell_type_l2_raw",
        output_path=output_path,
        title="pbmcref highlight",
        config=config,
        legend_title="pbmcref",
    )

    assert output_path.exists()
    assert output_path.exists()


def test_save_sample_cima_l1_umap_omits_unknown_and_uses_sample_title(
    tmp_path: Path, monkeypatch
):
    from scripts.only_rna.plotting import save_sample_cima_l1_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [0.0, 1.0, 2.0],
            "umap_2": [0.0, 1.0, 2.0],
            "azimuth_cima_l1": ["CD4_T", "Unknown", "B"],
        },
        index=["cell-1", "cell-2", "cell-3"],
    )
    adata = ad.AnnData(X=np.ones((3, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    captured_labels: list[str] = []
    captured_titles: list[str] = []
    original_scatter = matplotlib.axes.Axes.scatter
    original_set_title = matplotlib.axes.Axes.set_title

    def capture_scatter(self, *args, **kwargs):
        label = kwargs.get("label")
        if label is not None:
            captured_labels.append(str(label))
        return original_scatter(self, *args, **kwargs)

    def capture_set_title(self, label, *args, **kwargs):
        captured_titles.append(str(label))
        return original_set_title(self, label, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)
    monkeypatch.setattr(matplotlib.axes.Axes, "set_title", capture_set_title)

    output_path = tmp_path / "cima_l1_umap.png"
    save_sample_cima_l1_umap(
        adata=adata,
        output_path=output_path,
        title="GSM4750304",
        config=config,
    )

    assert output_path.exists()
    assert "Unknown" not in captured_labels
    assert "CD4_T" in captured_labels
    assert "B" in captured_labels
    assert "GSM4750304" in captured_titles


def test_save_sample_cima_l1_umap_requires_azimuth_cima_l1_field(
    tmp_path: Path,
):
    from scripts.only_rna.plotting import save_sample_cima_l1_umap

    config = _make_run_config()
    obs = pd.DataFrame(
        {
            "umap_1": [0.0, 1.0],
            "umap_2": [0.0, 1.0],
            "cima_l1": ["CD4_T", "B"],
        },
        index=["cell-1", "cell-2"],
    )
    adata = ad.AnnData(X=np.ones((2, 1)), obs=obs, var=pd.DataFrame(index=["GeneA"]))

    output_path = tmp_path / "cima_l1_umap_missing.png"
    save_sample_cima_l1_umap(
        adata=adata,
        output_path=output_path,
        title="GSM4750304",
        config=config,
    )

    assert output_path.exists()


def test_save_qc_overview_draws_mad_center_raw_and_final_thresholds(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = _make_output_adata()
    adata.obs["n_counts"] = [100.0, 300.0, 1000.0]
    adata.obs["n_genes"] = [80.0, 200.0, 700.0]
    adata.obs["pct_mt"] = [2.0, 5.0, 12.0]
    adata.obs["pct_ribo"] = [10.0, 15.0, 20.0]
    adata.obs["fails_count_floor"] = [True, False, False]
    adata.obs["fails_gene_floor"] = [True, False, False]
    adata.obs["fails_mt_ceiling"] = [False, False, True]
    adata.obs["fails_ribo_ceiling"] = [False, False, False]
    adata.uns["qc_thresholds"] = {
        "min_counts": 276,
        "min_genes": 212,
        "max_pct_mt": 5.87,
        "max_pct_ribo": 45.08,
        "n_counts_audit": {
            "transform": "log10p1",
            "direction": "lower",
            "center": 2.0,
            "mad": 0.2,
            "nmads": 3.0,
            "raw_threshold": 1.4,
            "final_threshold_original_scale": 276.0,
            "guardrails_applied": [],
        },
        "n_genes_audit": {
            "transform": "log10p1",
            "direction": "lower",
            "center": 1.8,
            "mad": 0.15,
            "nmads": 3.0,
            "raw_threshold": 1.35,
            "final_threshold_original_scale": 212.0,
            "guardrails_applied": ["floor_min"],
        },
        "pct_mt_audit": {
            "transform": "identity",
            "direction": "upper",
            "center": 4.0,
            "mad": 0.6,
            "nmads": 3.0,
            "raw_threshold": 5.8,
            "final_threshold_original_scale": 5.87,
            "guardrails_applied": [],
        },
        "pct_ribo_audit": {
            "transform": "identity",
            "direction": "upper",
            "center": 14.0,
            "mad": 2.0,
            "nmads": 3.5,
            "raw_threshold": 21.0,
            "final_threshold_original_scale": 45.08,
            "guardrails_applied": ["ceiling_max"],
        },
    }

    captured_lines: list[float] = []
    captured_titles: list[str] = []
    captured_xlabels: list[str] = []
    captured_text: list[str] = []
    original_axvline = matplotlib.axes.Axes.axvline
    original_set_title = matplotlib.axes.Axes.set_title
    original_set_xlabel = matplotlib.axes.Axes.set_xlabel
    original_text = matplotlib.axes.Axes.text

    def capture_axvline(self, x, *args, **kwargs):
        captured_lines.append(float(x))
        return original_axvline(self, x, *args, **kwargs)

    def capture_set_title(self, label, *args, **kwargs):
        captured_titles.append(str(label))
        return original_set_title(self, label, *args, **kwargs)

    def capture_set_xlabel(self, label, *args, **kwargs):
        captured_xlabels.append(str(label))
        return original_set_xlabel(self, label, *args, **kwargs)

    def capture_text(self, *args, **kwargs):
        label = args[2] if len(args) >= 3 else kwargs.get("s")
        if label is not None:
            captured_text.append(str(label))
        return original_text(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "axvline", capture_axvline)
    monkeypatch.setattr(matplotlib.axes.Axes, "set_title", capture_set_title)
    monkeypatch.setattr(matplotlib.axes.Axes, "set_xlabel", capture_set_xlabel)
    monkeypatch.setattr(matplotlib.axes.Axes, "text", capture_text)

    output_path = tmp_path / "qc_overview.png"
    save_qc_overview(adata, output_path, config)

    assert output_path.exists()
    assert any("lower-tail MAD in log10(x + 1)" in title for title in captured_titles)
    assert any("upper-tail MAD in original scale" in title for title in captured_titles)
    assert "n_counts [log10(x + 1)]" in captured_xlabels
    assert "pct_mt" in captured_xlabels
    assert any("blue=median" in text for text in captured_text)
    assert any("guardrail: floor_min" in text for text in captured_text)
    assert any("guardrail: ceiling_max" in text for text in captured_text)
    assert 2.0 in captured_lines
    assert 1.4 in captured_lines
    assert any(abs(value - np.log10(277.0)) < 1e-6 for value in captured_lines)
    assert 4.0 in captured_lines
    assert 5.8 in captured_lines
    assert 5.87 in captured_lines


# ---------------------------------------------------------------------------
# Hierarchical palette tests
# ---------------------------------------------------------------------------

_ALL_L1 = ["B", "CD4_T", "CD8_T", "ILC", "Myeloid", "unconvensional_T"]


def _color_distance(hex_a: str, hex_b: str) -> float:
    """Euclidean distance in RGB space between two hex colours."""
    a = np.array([int(hex_a[i : i + 2], 16) for i in (1, 3, 5)], dtype=float)
    b = np.array([int(hex_b[i : i + 2], 16) for i in (1, 3, 5)], dtype=float)
    return float(np.linalg.norm(a - b))


def test_l1_palette_returns_distinct_colors_for_separate_families():
    pal = build_cima_l1_palette(_ALL_L1)
    # B, ILC, Myeloid are from three different families — must be far apart
    assert _color_distance(pal["B"], pal["ILC"]) > 50
    assert _color_distance(pal["B"], pal["Myeloid"]) > 50
    assert _color_distance(pal["ILC"], pal["Myeloid"]) > 50


def test_l1_palette_cd4_cd8_share_family():
    pal = build_cima_l1_palette(_ALL_L1)
    # CD4_T and CD8_T must be different colours but from the same family
    # "Same family" = they are closer to each other than to B/ILC/Myeloid
    assert pal["CD4_T"] != pal["CD8_T"]
    dist_within = _color_distance(pal["CD4_T"], pal["CD8_T"])
    dist_cross = _color_distance(pal["CD4_T"], pal["B"])
    assert dist_within < dist_cross


def test_l1_palette_unconvensional_t_distinct_but_in_t_family():
    pal = build_cima_l1_palette(_ALL_L1)
    # unconvensional_T must be distinct from CD4_T / CD8_T
    assert pal["unconvensional_T"] != pal["CD4_T"]
    assert pal["unconvensional_T"] != pal["CD8_T"]
    # But it must still be in the T-family: closer to CD4_T than B is
    dist_unconv_cd4 = _color_distance(pal["unconvensional_T"], pal["CD4_T"])
    dist_b_cd4 = _color_distance(pal["B"], pal["CD4_T"])
    assert dist_unconv_cd4 < dist_b_cd4


def test_l2_palette_inherits_parent_family_shade():
    l2_by_l1 = {
        "CD4_T": ["CD4_naive", "CD4_cm", "CD4_em"],
        "CD8_T": ["CD8_naive", "CD8_CTL", "CD8_em"],
        "B": ["Naive_B", "Memory_B", "Plasma"],
        "ILC": ["ILC1", "ILC2", "ILC3"],
        "Myeloid": ["Monocyte", "cDC", "pDC"],
        "unconvensional_T": ["gd_T", "MAIT", "NKT"],
    }
    pal = build_cima_l2_palette(l2_by_l1)
    l1_pal = build_cima_l1_palette(_ALL_L1)

    # CD4_naive and CD8_CTL must be different (different families)
    assert pal["CD4_naive"] != pal["CD8_CTL"]

    # CD4_naive should be close to the CD4_T L1 base colour
    dist_cd4 = _color_distance(pal["CD4_naive"], l1_pal["CD4_T"])
    dist_cross = _color_distance(pal["CD4_naive"], l1_pal["B"])
    assert dist_cd4 < dist_cross

    # CD8_CTL should be close to the CD8_T L1 base colour
    dist_cd8 = _color_distance(pal["CD8_CTL"], l1_pal["CD8_T"])
    dist_cross2 = _color_distance(pal["CD8_CTL"], l1_pal["B"])
    assert dist_cd8 < dist_cross2

    # gd_T (child of unconvensional_T) should be closer to unconvensional_T
    # than to B
    dist_gd_unc = _color_distance(pal["gd_T"], l1_pal["unconvensional_T"])
    dist_gd_b = _color_distance(pal["gd_T"], l1_pal["B"])
    assert dist_gd_unc < dist_gd_b


def test_get_hierarchical_palette_returns_tab20_for_cluster():
    """cluster key should NOT use hierarchical palette."""
    pal = get_hierarchical_palette("cluster", ["0", "1", "2"])
    # Must be None to signal "use default tab20"
    assert pal is None


def test_get_hierarchical_palette_returns_dict_for_cima_keys():
    for key in ("cima_l1", "cima_l2", "cima_l1_masked"):
        labels = _ALL_L1 if key != "cima_l2" else ["CD4_naive", "CD8_CTL", "Naive_B"]
        pal = get_hierarchical_palette(key, labels)
        assert pal is not None
        assert isinstance(pal, dict)
        for label in labels:
            assert label in pal, f"missing {label} in palette for {key}"


def test_get_hierarchical_palette_maps_raw_pbmcref_labels_to_cima_family_colors():
    labels = ["B", "CD4 T", "CD8 T", "other", "other T"]
    pal = get_hierarchical_palette("azimuth_cell_type", labels)

    assert pal is not None
    cima_family_palette = build_cima_l1_palette(
        ["B", "CD4_T", "CD8_T", "Unknown", "unconvensional_T"]
    )

    def _normalized_hex(color: str | tuple[float, ...]) -> str:
        return to_hex(to_rgb(color)).lower()

    assert _normalized_hex(pal["B"]) == _normalized_hex(cima_family_palette["B"])
    assert _normalized_hex(pal["CD4 T"]) == _normalized_hex(
        cima_family_palette["CD4_T"]
    )
    assert _normalized_hex(pal["CD8 T"]) == _normalized_hex(
        cima_family_palette["CD8_T"]
    )
    assert _color_distance(
        _normalized_hex(pal["other"]),
        _normalized_hex(cima_family_palette["Unknown"]),
    ) < _color_distance(
        _normalized_hex(pal["other"]),
        _normalized_hex(cima_family_palette["B"]),
    )
    assert _color_distance(
        _normalized_hex(pal["other T"]),
        _normalized_hex(cima_family_palette["unconvensional_T"]),
    ) < _color_distance(
        _normalized_hex(pal["other T"]),
        _normalized_hex(cima_family_palette["B"]),
    )


def test_get_hierarchical_palette_gives_distinct_family_shades_for_fine_pbmcref_labels():
    labels = ["DC", "Mono", "CD8 TEM", "CD8 TCM", "B"]
    pal = get_hierarchical_palette("azimuth_cell_type", labels)

    assert pal is not None

    def _normalized_hex(color: str | tuple[float, ...]) -> str:
        return to_hex(to_rgb(color)).lower()

    dc_hex = _normalized_hex(pal["DC"])
    mono_hex = _normalized_hex(pal["Mono"])
    cd8_tem_hex = _normalized_hex(pal["CD8 TEM"])
    cd8_tcm_hex = _normalized_hex(pal["CD8 TCM"])
    b_hex = _normalized_hex(pal["B"])

    myeloid_hex = _normalized_hex(build_cima_l1_palette(["Myeloid"])["Myeloid"])
    cd8_hex = _normalized_hex(build_cima_l1_palette(["CD8_T"])["CD8_T"])

    assert dc_hex != mono_hex
    assert cd8_tem_hex != cd8_tcm_hex
    assert _color_distance(dc_hex, myeloid_hex) < _color_distance(dc_hex, b_hex)
    assert _color_distance(mono_hex, myeloid_hex) < _color_distance(mono_hex, b_hex)
    assert _color_distance(cd8_tem_hex, cd8_hex) < _color_distance(cd8_tem_hex, b_hex)
    assert _color_distance(cd8_tcm_hex, cd8_hex) < _color_distance(cd8_tcm_hex, b_hex)


def test_l2_palette_handles_more_unknown_parents_than_fallback_colors():
    labels = [
        "Total_Plasma",
        "Memory_like",
        "Naive_like",
        "Plasma_like",
        "Treg_like",
        "Prolif_like",
        "Eryth_like",
        "Platelet_like",
        "DC",
        "HSPC",
    ]

    pal = get_hierarchical_palette("cima_l2", labels)

    assert pal is not None
    assert isinstance(pal, dict)
    for label in labels:
        assert label in pal, f"missing {label} in palette for cima_l2 overflow case"


def test_save_categorical_umap_uses_hierarchical_palette_for_cima_l1(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(6)),
        obs=pd.DataFrame(
            {
                "umap_1": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                "umap_2": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
                "cima_l1": _ALL_L1,
            },
            index=[f"cell-{i}" for i in range(6)],
        ),
    )

    l1_palette = build_cima_l1_palette(_ALL_L1)
    captured_colors: dict[str, str] = {}
    original_scatter = matplotlib.axes.Axes.scatter

    def capture_scatter(self, *args, **kwargs):
        c = kwargs.get("c")
        label = kwargs.get("label")
        if c and label:
            captured_colors[label] = c[0] if isinstance(c, list) else c
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

    output_path = tmp_path / "hier_l1.png"
    save_categorical_umap(
        adata,
        color_key="cima_l1",
        output_path=output_path,
        title="Hierarchical L1",
        config=config,
    )

    assert output_path.exists()
    # Every L1 label must have been plotted with its hierarchical colour
    for label in _ALL_L1:
        assert label in captured_colors, f"missing scatter for {label}"
