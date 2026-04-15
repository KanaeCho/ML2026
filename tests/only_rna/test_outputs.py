from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib.axes
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse

from scripts.only_rna.models import PlottingConfig, QcThresholds, RunConfig
from scripts.only_rna.outputs import write_sample_outputs
from scripts.only_rna.plotting import save_categorical_umap


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
    )


def _make_output_adata() -> ad.AnnData:
    return ad.AnnData(
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
                "cima_l1": ["T cell", pd.NA, "B cell"],
                "cima_l2": ["CD4 T", pd.NA, "Naive B"],
                "cima_l1_masked": ["T cell", pd.NA, "Unknown"],
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
    assert image.shape[1] == int(config.plotting.umap_width * config.plotting.dpi)
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
        sample_dir / "validation_result.csv",
        sample_dir / "GSM123456.h5ad",
        sample_dir / "matrix" / "matrix.mtx",
        sample_dir / "matrix" / "barcodes.tsv.gz",
        sample_dir / "matrix" / "features.tsv.gz",
        sample_dir / "umap_rna_clusters.png",
        sample_dir / "umap_rna_cima_cell_type_l1.png",
        sample_dir / "umap_rna_cima_cell_type_l2.png",
        sample_dir / "umap_rna_cima_cell_type_l1_masked.png",
    ]
    for path in expected_files:
        assert path.exists(), path

    assert sorted(path.name for path in sample_dir.glob("*.h5ad")) == ["GSM123456.h5ad"]

    metadata = pd.read_csv(sample_dir / "metadata.csv")
    metadata_qc = pd.read_csv(sample_dir / "metadata_qc.csv")
    qc_summary = pd.read_csv(sample_dir / "qc_summary.csv")
    validation = pd.read_csv(sample_dir / "validation_result.csv")

    assert metadata["cell_id"].tolist() == ["cell-1", "cell-2", "cell-3"]
    assert metadata_qc["cell_id"].tolist() == ["cell-1", "cell-3"]
    assert metadata_qc["pass_qc"].tolist() == [True, True]

    assert qc_summary.to_dict(orient="records") == [
        {
            "sample_id": "GSM123456",
            "gse": "GSE123456",
            "n_cells_total": 3,
            "n_cells_pass_qc": 2,
            "n_cells_fail_qc": 1,
        }
    ]

    completion = validation.loc[validation["check_name"] == "completion"]
    assert completion["passed"].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:metadata.csv", "passed"
    ].tolist() == [True]
    assert validation.loc[
        validation["check_name"] == "output_presence:umap_rna_cima_cell_type_l1.png",
        "passed",
    ].tolist() == [True]
