from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO_ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO_ROOT / "scripts" / "process" / "organize_tea_seq_outputs.py"


def load_module():
    spec = importlib.util.spec_from_file_location(
        "organize_tea_seq_outputs", MODULE_PATH
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load module from {MODULE_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class OrganizeTeaSeqOutputsTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.TemporaryDirectory()
        self.project_root = Path(self.temp_dir.name)
        self.dataset_dir = self.project_root / "output" / "GSE214546"
        self.sample_dir = self.dataset_dir / "GSM0000001"
        (self.sample_dir / "matrix").mkdir(parents=True)
        (self.project_root / "data" / "raw" / "GSE214546").mkdir(parents=True)
        self.write_qc_summary(
            self.sample_dir / "qc_summary.csv",
            {
                "input_cells": "100",
                "pass_qc": "80",
                "qc_rate": "80.00%",
                "median_TSS_enrichment": "7.50",
                "median_FRiP": "0.8700",
                "median_fragments": "5000",
                "query_cluster_count": "4",
                "cima_unique_l4_labels": "7",
                "median_cima_l4_score": "0.5200",
            },
        )

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def write_qc_summary(self, path: Path, metrics: dict[str, str]) -> None:
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["metric", "value"])
            for key, value in metrics.items():
                writer.writerow([key, value])

    def write_validation(self, path: Path, rows: list[dict[str, object]]) -> None:
        if not rows:
            raise ValueError("rows must not be empty")
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            for row in rows:
                writer.writerow(row)

    def make_accepted_csvs(self) -> None:
        for name in (
            "cima_cluster_centroid_labels.csv",
            "cima_cluster_centroid_label_summary.csv",
            "adt_cluster_broad_labels.csv",
            "adt_cluster_broad_label_summary.csv",
        ):
            (self.sample_dir / name).write_text("header\n", encoding="utf-8")

    def test_dry_run_reports_l1_move_without_touching_files(self) -> None:
        module = load_module()
        self.make_accepted_csvs()
        self.write_validation(
            self.sample_dir / "validation_result.csv",
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.35",
                    "cima_l4_score_margin": "0.05",
                }
            ],
        )
        legacy_dir = self.dataset_dir / "l1"
        legacy_dir.mkdir(parents=True)
        legacy_png = legacy_dir / "GSM0000001_umap_atac_cima_cell_type_l1_cluster_centroid.png"
        legacy_png.write_text("png", encoding="utf-8")

        result = module.organize_dataset(
            project_root=self.project_root,
            gse="GSE214546",
            dry_run=True,
            delete_legacy=False,
            score_threshold=0.4,
            margin_threshold=0.1,
        )

        self.assertTrue(legacy_png.exists())
        self.assertFalse(
            (
                self.sample_dir / "umap_atac_cima_cell_type_l1_cluster_centroid.png"
            ).exists()
        )
        actions = {(row["relative_path"], row["action"]) for row in result["legacy_manifest_rows"]}
        self.assertIn(
            (
                "output/GSE214546/l1/GSM0000001_umap_atac_cima_cell_type_l1_cluster_centroid.png",
                "would_move",
            ),
            actions,
        )

    def test_real_run_moves_l1_png_and_renders_missing_adt_png(self) -> None:
        module = load_module()
        self.make_accepted_csvs()
        self.write_validation(
            self.sample_dir / "validation_result.csv",
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.55",
                    "cima_l4_score_margin": "0.25",
                }
            ],
        )
        legacy_dir = self.dataset_dir / "l1"
        legacy_dir.mkdir(parents=True)
        legacy_png = legacy_dir / "GSM0000001_umap_atac_cima_cell_type_l1_cluster_centroid.png"
        legacy_png.write_text("png", encoding="utf-8")

        with mock.patch.object(module, "render_adt_broad_outputs") as render_adt:
            def create_adt_png(*_args, **kwargs):
                (
                    self.sample_dir / "umap_atac_adt_cluster_broad_celltype.png"
                ).write_text("png", encoding="utf-8")

            render_adt.side_effect = create_adt_png

            module.organize_dataset(
                project_root=self.project_root,
                gse="GSE214546",
                dry_run=False,
                delete_legacy=False,
                score_threshold=0.4,
                margin_threshold=0.1,
            )

        self.assertFalse(legacy_png.exists())
        self.assertTrue(
            (
                self.sample_dir / "umap_atac_cima_cell_type_l1_cluster_centroid.png"
            ).exists()
        )
        self.assertTrue(
            (self.sample_dir / "umap_atac_adt_cluster_broad_celltype.png").exists()
        )
        render_adt.assert_called_once()

    def test_low_confidence_fallback_uses_score_and_margin_for_old_schema(self) -> None:
        module = load_module()
        validation_path = self.sample_dir / "validation_result.csv"
        self.write_validation(
            validation_path,
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.30",
                    "cima_l4_score_margin": "0.20",
                },
                {
                    "cell_barcode": "cell2",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.70",
                    "cima_l4_score_margin": "0.20",
                },
                {
                    "cell_barcode": "cell3",
                    "seurat_clusters": "1",
                    "cima_cell_type_l1": "B",
                    "cima_l4_score": "0.60",
                    "cima_l4_score_margin": "0.05",
                },
            ],
        )

        summary = module.summarize_validation_result(
            validation_path,
            score_threshold=0.4,
            margin_threshold=0.1,
        )

        self.assertAlmostEqual(summary["low_conf_cell_frac"], 2 / 3, places=6)
        self.assertAlmostEqual(summary["cluster_majority_purity_median"], 1.0, places=6)

    def test_low_confidence_prefers_existing_columns_for_new_schema(self) -> None:
        module = load_module()
        validation_path = self.sample_dir / "validation_result.csv"
        self.write_validation(
            validation_path,
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l1_low_confidence": "TRUE",
                    "cima_l1_cluster_purity": "0.75",
                    "cima_l4_score": "0.90",
                    "cima_l4_score_margin": "0.90",
                },
                {
                    "cell_barcode": "cell2",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l1_low_confidence": "FALSE",
                    "cima_l1_cluster_purity": "0.75",
                    "cima_l4_score": "0.10",
                    "cima_l4_score_margin": "0.01",
                },
            ],
        )

        summary = module.summarize_validation_result(
            validation_path,
            score_threshold=0.4,
            margin_threshold=0.1,
            low_purity_threshold=0.8,
        )

        self.assertAlmostEqual(summary["low_conf_cell_frac"], 0.5, places=6)
        self.assertAlmostEqual(summary["cluster_majority_purity_median"], 0.75, places=6)
        self.assertAlmostEqual(summary["low_purity_cluster_cell_frac"], 1.0, places=6)

    def test_old_schema_does_not_treat_blank_margin_as_low_confidence(self) -> None:
        module = load_module()
        validation_path = self.sample_dir / "validation_result.csv"
        self.write_validation(
            validation_path,
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.70",
                    "cima_l4_score_margin": "",
                },
                {
                    "cell_barcode": "cell2",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.30",
                    "cima_l4_score_margin": "",
                },
            ],
        )

        summary = module.summarize_validation_result(
            validation_path,
            score_threshold=0.4,
            margin_threshold=0.1,
        )

        self.assertAlmostEqual(summary["low_conf_cell_frac"], 0.5, places=6)

    def test_delete_legacy_requires_accepted_outputs_complete(self) -> None:
        module = load_module()
        self.make_accepted_csvs()
        self.write_validation(
            self.sample_dir / "validation_result.csv",
            [
                {
                    "cell_barcode": "cell1",
                    "seurat_clusters": "0",
                    "cima_cell_type_l1": "CD4_T",
                    "cima_l4_score": "0.55",
                    "cima_l4_score_margin": "0.25",
                }
            ],
        )
        l1_dir = self.dataset_dir / "l1"
        l1_dir.mkdir(parents=True)
        (
            l1_dir / "GSM0000001_umap_atac_cima_cell_type_l1_cluster_centroid.png"
        ).write_text("png", encoding="utf-8")
        gex_dir = self.dataset_dir / "gex_refined_cima_compare" / "GSM0000001"
        gex_dir.mkdir(parents=True)
        (gex_dir / "legacy.txt").write_text("legacy", encoding="utf-8")

        with mock.patch.object(module, "render_adt_broad_outputs"):
            result = module.organize_dataset(
                project_root=self.project_root,
                gse="GSE214546",
                dry_run=False,
                delete_legacy=True,
                score_threshold=0.4,
                margin_threshold=0.1,
            )

        self.assertTrue(l1_dir.exists())
        self.assertTrue((self.dataset_dir / "gex_refined_cima_compare").exists())
        legacy_actions = [row["action"] for row in result["legacy_manifest_rows"]]
        self.assertIn("kept_incomplete", legacy_actions)


if __name__ == "__main__":
    unittest.main()
