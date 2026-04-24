from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
import scripts.only_rna.plotting as plotting_module

from scripts.only_rna.plotting import get_hierarchical_palette
from scripts.only_rna.models import (
    AnnotationConfig,
    PlottingConfig,
    QcThresholds,
    RunConfig,
)
from scripts.process.compare_gse192391_annotation_methods import (
    map_method_labels_to_cima_l1,
    pick_best_annotation_method,
    CIMA_L1_ORDER,
    _run_azimuth_r,
    rebuild_gse192391_broad_summary,
    run_azimuth_annotation,
    run_celltypist_annotation,
    run_scanvi_annotation,
    run_singler_annotation,
    run_gse192391_comparison_batch,
    prepare_comparison_adata,
    write_sample_comparison_outputs,
)


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
            methods=["azimuth", "cell_typist", "singler", "scanvi"]
        ),
    )


def _make_sample_adata() -> ad.AnnData:
    obs = pd.DataFrame(
        {
            "pass_qc": [True, True, False],
            "umap_1": [0.0, 1.0, np.nan],
            "umap_2": [1.0, 0.0, np.nan],
        },
        index=["cell-1", "cell-2", "cell-3"],
    )
    obs["cluster"] = pd.Series(
        ["0", "1", None],
        index=obs.index,
        dtype="object",
    )

    return ad.AnnData(
        X=sparse.csr_matrix(
            np.array(
                [
                    [10.0, 0.0, 1.0],
                    [0.0, 8.0, 2.0],
                    [1.0, 1.0, 1.0],
                ]
            )
        ),
        obs=obs,
        var=pd.DataFrame(
            {
                "feature_id": ["GeneA", "GeneB", "GeneC"],
                "feature_name": ["GeneA", "GeneB", "GeneC"],
            },
            index=["GeneA", "GeneB", "GeneC"],
        ),
    )


def test_prepare_comparison_adata_marks_blocked_azimuth_and_preserves_fail_qc_na():
    adata = _make_sample_adata()

    prepared = prepare_comparison_adata(
        adata,
        azimuth_labels=None,
        azimuth_block_reason="Azimuth installation blocked: missing SeuratDisk/hdf5r",
        celltypist_labels=pd.Series(
            ["T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        singler_labels=pd.Series(
            ["T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        scanvi_labels=pd.Series(["T", "B"], index=["cell-1", "cell-2"], dtype="string"),
    )

    assert str(prepared.obs["azimuth_cell_type"].dtype) == "string"
    assert str(prepared.obs["celltypist_cell_type"].dtype) == "string"
    assert str(prepared.obs["singler_cell_type"].dtype) == "string"
    assert str(prepared.obs["scanvi_cell_type"].dtype) == "string"

    assert (
        prepared.obs.loc["cell-1", "azimuth_cell_type"]
        == "Blocked: Azimuth unavailable"
    )
    assert (
        prepared.obs.loc["cell-2", "azimuth_cell_type"]
        == "Blocked: Azimuth unavailable"
    )
    assert prepared.obs.loc["cell-1", "celltypist_cell_type"] == "T"
    assert prepared.obs.loc["cell-2", "singler_cell_type"] == "B"
    assert pd.isna(prepared.obs.loc["cell-3", "azimuth_cell_type"])
    assert pd.isna(prepared.obs.loc["cell-3", "scanvi_cell_type"])


def test_write_sample_comparison_outputs_writes_panel_and_status_json(tmp_path: Path):
    config = _make_run_config()
    prepared = prepare_comparison_adata(
        _make_sample_adata(),
        azimuth_labels=None,
        azimuth_block_reason="Azimuth installation blocked: missing SeuratDisk/hdf5r",
        celltypist_labels=pd.Series(
            ["T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        singler_labels=pd.Series(
            ["T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        scanvi_labels=pd.Series(["T", "B"], index=["cell-1", "cell-2"], dtype="string"),
    )
    method_status = {
        "azimuth": {
            "status": "blocked",
            "detail": "Azimuth installation blocked: missing SeuratDisk/hdf5r",
        },
        "celltypist": {"status": "ok", "detail": "Immune_All_Low.pkl"},
        "singler": {"status": "ok", "detail": "HumanPrimaryCellAtlasData"},
        "scanvi": {"status": "ok", "detail": "local scANVI run"},
    }

    output_dir = write_sample_comparison_outputs(
        prepared,
        output_root=tmp_path,
        sample_id="GSM5746268",
        config=config,
        method_status=method_status,
    )

    compare_png = output_dir / "umap_rna_annotation_method_compare.png"
    status_json = output_dir / "annotation_method_status.json"

    assert compare_png.exists()
    assert status_json.exists()

    payload = json.loads(status_json.read_text(encoding="utf-8"))
    assert payload["sample_id"] == "GSM5746268"
    assert payload["methods"]["azimuth"]["status"] == "blocked"
    assert "SeuratDisk" in payload["methods"]["azimuth"]["detail"]
    assert payload["methods"]["celltypist"]["status"] == "ok"


def test_run_celltypist_annotation_returns_nonempty_labels_for_pass_qc_cells():
    labels = run_celltypist_annotation(
        _make_sample_adata(),
        model_name="Immune_All_Low.pkl",
    )

    assert list(labels.index) == ["cell-1", "cell-2"]
    assert labels.astype("string").notna().all()


def test_run_azimuth_annotation_returns_nonempty_labels_for_pass_qc_cells(monkeypatch):
    fake = pd.Series(["CD4 T", "B"], index=["cell-1", "cell-2"], dtype="string")
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods._run_azimuth_r",
        lambda adata, reference="pbmcref", annotation_level="l1", max_cells=1000: fake,
    )

    labels = run_azimuth_annotation(_make_sample_adata(), annotation_level="l1")

    assert list(labels.index) == ["cell-1", "cell-2"]
    assert labels.astype("string").notna().all()


def test_run_gse192391_comparison_batch_disables_azimuth_cell_cap_by_default(
    tmp_path: Path, monkeypatch
):
    input_root = tmp_path / "rna" / "GSE192391"
    output_root = tmp_path / "compare"
    config = _make_run_config()

    sample_dir = input_root / "GSM5746268"
    sample_dir.mkdir(parents=True, exist_ok=True)
    (sample_dir / "GSM5746268.h5ad").write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.ad.read_h5ad",
        lambda path: _make_sample_adata(),
    )

    captured: dict[str, object] = {}

    def _fake_azimuth(
        adata, annotation_level="l1", reference="pbmcref", max_cells=1000
    ):
        captured["max_cells"] = max_cells
        return pd.Series(["CD4 T", "B"], index=["cell-1", "cell-2"], dtype="string")

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_azimuth_annotation",
        _fake_azimuth,
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_celltypist_annotation",
        lambda adata, model_name="Immune_All_Low.pkl": pd.Series(
            ["Tcm/Naive helper T cells", "Naive B cells"],
            index=["cell-1", "cell-2"],
            dtype="string",
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_singler_annotation",
        lambda adata: pd.Series(
            ["T_cells", "B_cell"], index=["cell-1", "cell-2"], dtype="string"
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_scanvi_annotation",
        lambda adata, max_epochs=1: pd.Series(
            ["CD4_T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
    )

    run_gse192391_comparison_batch(
        input_root=input_root,
        output_root=output_root,
        config=config,
    )

    assert captured["max_cells"] is None


def test_run_azimuth_r_uses_mtx_and_seurat_object_path(tmp_path: Path, monkeypatch):
    commands: dict[str, str] = {}
    labels_csv = tmp_path / "labels.csv"

    def _fake_run(cmd, check, capture_output, text):
        commands["r_code"] = cmd[2]
        labels_csv.write_text(
            "cell_id,label\ncell-1,CD4 T\ncell-2,B\n", encoding="utf-8"
        )

        class _Done:
            returncode = 0

        return _Done()

    original_read_csv = pd.read_csv

    def _spy_read_csv(path, *args, **kwargs):
        path = Path(path)
        if path == labels_csv:
            commands["output_path"] = str(path)
        return original_read_csv(path, *args, **kwargs)

    class _TempDir:
        def __init__(self, path: Path):
            self.path = path

        def __enter__(self):
            return str(self.path)

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.tempfile.TemporaryDirectory",
        lambda: _TempDir(tmp_path),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.subprocess.run",
        _fake_run,
    )
    monkeypatch.setattr(pd, "read_csv", _spy_read_csv)

    labels = _run_azimuth_r(_make_sample_adata(), annotation_level="l1", max_cells=2)

    assert list(labels.index) == ["cell-1", "cell-2"]
    assert "ReadMtx(" in commands["r_code"]
    assert f"mtx = '{(tmp_path / 'matrix.mtx').as_posix()}'" in commands["r_code"]
    assert f"cells = '{(tmp_path / 'barcodes.tsv').as_posix()}'" in commands["r_code"]
    assert (
        f"features = '{(tmp_path / 'features.tsv').as_posix()}'" in commands["r_code"]
    )
    assert "CreateSeuratObject(counts = counts)" in commands["r_code"]
    assert "cells = NULL" not in commands["r_code"]
    assert "features = NULL" not in commands["r_code"]


def test_run_singler_annotation_returns_nonempty_labels_for_pass_qc_cells():
    labels = run_singler_annotation(_make_sample_adata())

    assert list(labels.index) == ["cell-1", "cell-2"]
    assert labels.astype("string").notna().all()


def test_run_scanvi_annotation_returns_nonempty_labels_for_pass_qc_cells():
    labels = run_scanvi_annotation(_make_sample_adata(), max_epochs=1)

    assert list(labels.index) == ["cell-1", "cell-2"]
    assert labels.astype("string").notna().all()


def test_run_gse192391_comparison_batch_writes_per_sample_outputs_and_summary(
    tmp_path: Path, monkeypatch
):
    input_root = tmp_path / "rna" / "GSE192391"
    output_root = tmp_path / "compare"
    config = _make_run_config()

    for sample_id in ["GSM5746268", "GSM5746269"]:
        sample_dir = input_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / f"{sample_id}.h5ad").write_text("fixture\n", encoding="utf-8")

    fake_labels = pd.Series(["T", "B"], index=["cell-1", "cell-2"], dtype="string")

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_azimuth_annotation",
        lambda adata, annotation_level="l1", reference="pbmcref", max_cells=1000: (
            pd.Series(
                ["CD4 T", "B"],
                index=["cell-1", "cell-2"],
                dtype="string",
            )
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_celltypist_annotation",
        lambda adata, model_name="Immune_All_Low.pkl": fake_labels,
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_singler_annotation",
        lambda adata: fake_labels,
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_scanvi_annotation",
        lambda adata, max_epochs=1: fake_labels,
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.ad.read_h5ad",
        lambda path: _make_sample_adata(),
    )

    summary = run_gse192391_comparison_batch(
        input_root=input_root,
        output_root=output_root,
        config=config,
    )

    assert summary["sample_id"].tolist() == ["GSM5746268", "GSM5746269"]
    assert summary["azimuth_status"].tolist() == ["ok", "ok"]
    assert summary["celltypist_status"].tolist() == ["ok", "ok"]

    for sample_id in ["GSM5746268", "GSM5746269"]:
        sample_dir = output_root / sample_id
        assert (sample_dir / "umap_rna_annotation_method_compare.png").exists()
        assert (sample_dir / "annotation_method_status.json").exists()

    summary_csv = output_root / "gse192391_annotation_method_summary.csv"
    assert summary_csv.exists()

    written = pd.read_csv(summary_csv)
    assert written["sample_id"].tolist() == ["GSM5746268", "GSM5746269"]


def test_map_method_labels_to_cima_l1_collapses_fine_labels_to_broad_lineages():
    labels = pd.Series(
        [
            "Tcm/Naive helper T cells",
            "Tem/Temra cytotoxic T cells",
            "Classical monocytes",
            "CD16+ NK cells",
            "Naive B cells",
            "MAIT cells",
        ],
        index=[f"cell-{i}" for i in range(6)],
        dtype="string",
    )

    mapped = map_method_labels_to_cima_l1(labels, method="celltypist")

    assert mapped.tolist() == [
        "CD4_T",
        "CD8_T",
        "Myeloid",
        "ILC",
        "B",
        "unconvensional_T",
    ]


def test_map_method_labels_to_cima_l1_supports_azimuth_l1_labels():
    labels = pd.Series(
        ["CD4 T", "CD8 T", "NK", "B", "Platelet", "other"],
        index=[f"cell-{i}" for i in range(6)],
        dtype="string",
    )

    mapped = map_method_labels_to_cima_l1(labels, method="azimuth")

    assert mapped.tolist() == [
        "CD4_T",
        "CD8_T",
        "ILC",
        "B",
        "Myeloid",
        "Unknown",
    ]


def test_map_method_labels_to_cima_l1_supports_additional_azimuth_broad_labels():
    labels = pd.Series(
        ["Mono", "other T", "DC", "B", "NK"],
        index=[f"cell-{i}" for i in range(5)],
        dtype="string",
    )

    mapped = map_method_labels_to_cima_l1(labels, method="azimuth")

    assert mapped.tolist() == [
        "Myeloid",
        "unconvensional_T",
        "Myeloid",
        "B",
        "ILC",
    ]


def test_pick_best_annotation_method_prefers_more_informative_broad_labels():
    summary = pd.DataFrame(
        [
            {
                "sample_id": "GSM1",
                "celltypist_unique": 6,
                "singler_unique": 4,
                "scanvi_unique": 6,
                "scanvi_vs_cima_l1_agreement": 0.72,
            },
            {
                "sample_id": "GSM2",
                "celltypist_unique": 6,
                "singler_unique": 4,
                "scanvi_unique": 5,
                "scanvi_vs_cima_l1_agreement": 0.70,
            },
        ]
    )

    best = pick_best_annotation_method(summary)

    assert best == "celltypist"


def test_run_gse192391_comparison_batch_writes_broad_summary_and_best_method(
    tmp_path: Path, monkeypatch
):
    input_root = tmp_path / "rna" / "GSE192391"
    output_root = tmp_path / "compare"
    config = _make_run_config()

    for sample_id in ["GSM5746268", "GSM5746269"]:
        sample_dir = input_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / f"{sample_id}.h5ad").write_text("fixture\n", encoding="utf-8")

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.ad.read_h5ad",
        lambda path: _make_sample_adata(),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_azimuth_annotation",
        lambda adata, annotation_level="l1", reference="pbmcref", max_cells=1000: (
            pd.Series(
                ["CD4 T", "B"],
                index=["cell-1", "cell-2"],
                dtype="string",
            )
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_celltypist_annotation",
        lambda adata, model_name="Immune_All_Low.pkl": pd.Series(
            ["Tcm/Naive helper T cells", "Naive B cells"],
            index=["cell-1", "cell-2"],
            dtype="string",
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_singler_annotation",
        lambda adata: pd.Series(
            ["T_cells", "B_cell"],
            index=["cell-1", "cell-2"],
            dtype="string",
        ),
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.run_scanvi_annotation",
        lambda adata, max_epochs=1: pd.Series(
            ["CD4_T", "B"],
            index=["cell-1", "cell-2"],
            dtype="string",
        ),
    )

    summary = run_gse192391_comparison_batch(
        input_root=input_root,
        output_root=output_root,
        config=config,
    )

    assert summary["azimuth_status"].tolist() == ["ok", "ok"]
    assert summary["celltypist_unique"].tolist() == [2, 2]
    assert summary["singler_unique"].tolist() == [2, 2]
    assert summary["scanvi_unique"].tolist() == [2, 2]
    assert (summary["best_method"] == "celltypist").all()

    broad_csv = output_root / "gse192391_annotation_method_summary_broad.csv"
    best_json = output_root / "best_annotation_method.json"
    assert broad_csv.exists()
    assert best_json.exists()

    broad = pd.read_csv(broad_csv)
    assert broad["sample_id"].tolist() == ["GSM5746268", "GSM5746269"]
    assert broad["celltypist_unique"].tolist() == [2, 2]

    payload = json.loads(best_json.read_text(encoding="utf-8"))
    assert payload["best_method"] == "celltypist"

    for sample_id in ["GSM5746268", "GSM5746269"]:
        sample_dir = output_root / sample_id
        assert (sample_dir / "umap_rna_best_annotation_method.png").exists()


def test_rebuild_gse192391_broad_summary_aggregates_existing_sample_outputs(
    tmp_path: Path,
):
    output_root = tmp_path / "compare"
    output_root.mkdir(parents=True, exist_ok=True)

    statuses = {
        "GSM5746268": {
            "sample_id": "GSM5746268",
            "methods": {
                "azimuth": {"status": "ok", "detail": "pbmcref"},
                "celltypist": {"status": "ok", "detail": "Immune_All_Low.pkl"},
                "singler": {"status": "ok", "detail": "HumanPrimaryCellAtlasData"},
                "scanvi": {"status": "ok", "detail": "local scANVI run (max_epochs=1)"},
            },
        },
        "GSM5746269": {
            "sample_id": "GSM5746269",
            "methods": {
                "azimuth": {"status": "ok", "detail": "pbmcref"},
                "celltypist": {"status": "ok", "detail": "Immune_All_Low.pkl"},
                "singler": {"status": "ok", "detail": "HumanPrimaryCellAtlasData"},
                "scanvi": {"status": "ok", "detail": "local scANVI run (max_epochs=1)"},
            },
        },
    }
    for sample_id, payload in statuses.items():
        sample_dir = output_root / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)
        (sample_dir / "annotation_method_status.json").write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
            encoding="utf-8",
        )

    sample_metrics = pd.DataFrame(
        [
            {
                "sample_id": "GSM5746268",
                "n_cells_total": 100,
                "n_cells_pass_qc": 90,
                "celltypist_unique": 6,
                "singler_unique": 4,
                "scanvi_unique": 6,
                "scanvi_vs_cima_l1_agreement": 0.74,
            },
            {
                "sample_id": "GSM5746269",
                "n_cells_total": 120,
                "n_cells_pass_qc": 110,
                "celltypist_unique": 6,
                "singler_unique": 4,
                "scanvi_unique": 5,
                "scanvi_vs_cima_l1_agreement": 0.64,
            },
        ]
    )

    summary = rebuild_gse192391_broad_summary(
        output_root=output_root,
        sample_metrics=sample_metrics,
    )

    assert summary["sample_id"].tolist() == ["GSM5746268", "GSM5746269"]
    assert summary["azimuth_status"].tolist() == ["ok", "ok"]
    assert (summary["best_method"] == "celltypist").all()
    assert (output_root / "gse192391_annotation_method_summary.csv").exists()
    assert (output_root / "gse192391_annotation_method_summary_broad.csv").exists()
    assert (output_root / "best_annotation_method.json").exists()


def test_prepare_comparison_adata_outputs_broad_l1_categories_only():
    adata = _make_sample_adata()

    prepared = prepare_comparison_adata(
        adata,
        azimuth_labels=None,
        azimuth_block_reason="Azimuth installation blocked: missing SeuratDisk/hdf5r",
        celltypist_labels=pd.Series(
            ["CD4_T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        singler_labels=pd.Series(
            ["CD4_T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        scanvi_labels=pd.Series(
            ["CD4_T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
    )

    for column in ["celltypist_cell_type", "singler_cell_type", "scanvi_cell_type"]:
        values = set(prepared.obs[column].dropna().astype(str))
        assert values <= set(CIMA_L1_ORDER)


def test_method_specific_broad_columns_use_cima_l1_palette():
    categories = ["B", "CD4_T", "CD8_T", "Myeloid", "ILC", "unconvensional_T"]

    for color_key in ["celltypist_cell_type", "singler_cell_type", "scanvi_cell_type"]:
        palette = get_hierarchical_palette(color_key, categories)
        assert palette is not None
        assert set(palette.keys()) == set(categories)


def test_method_specific_broad_columns_keep_unknown_as_fixed_gray():
    categories = [
        "B",
        "CD4_T",
        "CD8_T",
        "Myeloid",
        "ILC",
        "unconvensional_T",
        "Unknown",
    ]

    for color_key in ["celltypist_cell_type", "singler_cell_type", "scanvi_cell_type"]:
        palette = get_hierarchical_palette(color_key, categories)
        assert palette is not None
        assert "Unknown" in palette
        assert palette["Unknown"] == (0.74, 0.74, 0.74)


def test_azimuth_broad_column_uses_cima_l1_palette_with_unknown_gray():
    categories = [
        "B",
        "CD4_T",
        "CD8_T",
        "Myeloid",
        "ILC",
        "unconvensional_T",
        "Unknown",
    ]

    palette = get_hierarchical_palette("azimuth_cell_type", categories)

    assert palette is not None
    assert set(palette.keys()) == set(categories)
    assert palette["B"] == plotting_module.to_rgb("#2C7BB6")
    assert palette["CD4_T"] == plotting_module.to_rgb("#D7191C")
    assert palette["Unknown"] == (0.74, 0.74, 0.74)


def test_save_annotation_method_comparison_umap_uses_hierarchical_palette_for_method_columns(
    tmp_path: Path, monkeypatch
):
    adata = prepare_comparison_adata(
        _make_sample_adata(),
        azimuth_labels=None,
        azimuth_block_reason="Azimuth installation blocked: missing SeuratDisk/hdf5r",
        celltypist_labels=pd.Series(
            ["CD4_T", "Unknown"], index=["cell-1", "cell-2"], dtype="string"
        ),
        singler_labels=pd.Series(
            ["CD4_T", "B"], index=["cell-1", "cell-2"], dtype="string"
        ),
        scanvi_labels=pd.Series(
            ["CD8_T", "Myeloid"], index=["cell-1", "cell-2"], dtype="string"
        ),
    )

    def _forbid_tab20(*args, **kwargs):
        raise AssertionError("comparison plot should not fall back to tab20")

    monkeypatch.setattr(plotting_module.plt, "get_cmap", _forbid_tab20)

    plotting_module.save_annotation_method_comparison_umap(
        adata,
        output_path=tmp_path / "compare.png",
        title="compare",
        config=_make_run_config(),
    )

    assert (tmp_path / "compare.png").exists()


def test_run_azimuth_r_writes_feature_by_cell_matrix(tmp_path: Path, monkeypatch):
    from scipy.io import mmread

    labels_csv = tmp_path / "labels.csv"

    def _fake_run(cmd, check, capture_output, text):
        (tmp_path / "captured_r_code.txt").write_text(cmd[2], encoding="utf-8")
        labels_csv.write_text(
            "cell_id,label\ncell-1,CD4 T\ncell-2,B\n", encoding="utf-8"
        )
        return type("Completed", (), {"returncode": 0, "stdout": "", "stderr": ""})()

    class _TempDir:
        def __enter__(self):
            return str(tmp_path)

        def __exit__(self, exc_type, exc, tb):
            return False

    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.subprocess.run",
        _fake_run,
    )
    monkeypatch.setattr(
        "scripts.process.compare_gse192391_annotation_methods.tempfile.TemporaryDirectory",
        lambda: _TempDir(),
    )

    labels = _run_azimuth_r(_make_sample_adata(), annotation_level="l1", max_cells=2)

    assert labels.tolist() == ["CD4 T", "B"]
    matrix = mmread(tmp_path / "matrix.mtx").tocsr()
    assert matrix.shape == (3, 2)
    r_code = (tmp_path / "captured_r_code.txt").read_text(encoding="utf-8")
    assert "options(future.globals.maxSize = 2 * 1024^3)" in r_code
    assert 'future::plan("sequential")' in r_code
