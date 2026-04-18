from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.process import pipeline


class ResolveDataRootTests(unittest.TestCase):
    def test_prefers_workspace_data_when_present(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp)
            data_dir = project_root / "data"
            data_dir.mkdir()
            fallback_dir = project_root / "fallback_data"
            fallback_dir.mkdir()

            resolved = pipeline.resolve_data_root(
                project_root=project_root, fallback_candidates=[fallback_dir]
            )

            self.assertEqual(resolved, data_dir)

    def test_uses_first_existing_fallback_when_workspace_data_missing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            project_root = Path(tmp) / "project"
            project_root.mkdir()
            fallback_dir = Path(tmp) / "external_data"
            fallback_dir.mkdir()

            resolved = pipeline.resolve_data_root(
                project_root=project_root, fallback_candidates=[fallback_dir]
            )

            self.assertEqual(resolved, fallback_dir)


class DiscoverRnaSamplesTests(unittest.TestCase):
    def test_discovers_per_gsm_triplet_samples(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_root = Path(tmp)
            raw_dir = data_root / "raw" / "GSE167363"
            raw_dir.mkdir(parents=True)
            reference_dir = data_root / "reference"
            reference_dir.mkdir(parents=True)
            (reference_dir / "datasets.xlsx").write_text("placeholder\n")

            self._touch_triplet(raw_dir, "GSM5102900_HC1")
            self._touch_triplet(raw_dir, "GSM5102901_HC2")

            samples = pipeline.discover_rna_samples_from_local_layout(
                raw_dir=data_root / "raw", selected_gses=["GSE167363"]
            )

            self.assertEqual([sample.sample_id for sample in samples], ["GSM5102900", "GSM5102901"])
            self.assertTrue(all(sample.supported for sample in samples))
            self.assertEqual({sample.input_type for sample in samples}, {"triplet"})

    def test_prefers_pbmc_files_for_multi_tissue_gse(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_root = Path(tmp)
            raw_dir = data_root / "raw" / "GSE226039"
            raw_dir.mkdir(parents=True)
            reference_dir = data_root / "reference"
            reference_dir.mkdir(parents=True)
            (reference_dir / "datasets.xlsx").write_text("placeholder\n")

            self._touch_triplet(raw_dir, "GSM7061877_H44_Ileum", features_name="genes")
            self._touch_triplet(raw_dir, "GSM7061878_H44_PBMC", features_name="genes")
            self._touch_triplet(raw_dir, "GSM7061879_H44_Rectum", features_name="genes")

            samples = pipeline.discover_rna_samples_from_local_layout(
                raw_dir=data_root / "raw", selected_gses=["GSE226039"]
            )

            self.assertEqual(len(samples), 1)
            self.assertEqual(samples[0].sample_id, "GSM7061878")
            self.assertIn("PBMC", samples[0].matrix_path.name)

    def test_discovers_shared_gse_triplet_as_one_sample(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_root = Path(tmp)
            raw_dir = data_root / "raw" / "GSE149689"
            raw_dir.mkdir(parents=True)
            reference_dir = data_root / "reference"
            reference_dir.mkdir(parents=True)
            (reference_dir / "datasets.xlsx").write_text("placeholder\n")

            self._touch_shared_triplet(raw_dir, "GSE149689")

            samples = pipeline.discover_rna_samples_from_local_layout(
                raw_dir=data_root / "raw", selected_gses=["GSE149689"]
            )

            self.assertEqual(len(samples), 1)
            self.assertEqual(samples[0].sample_id, "GSE149689")
            self.assertEqual(samples[0].input_type, "triplet")

    def test_discovers_tar_archives(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_root = Path(tmp)
            raw_dir = data_root / "raw" / "GSE231794"
            raw_dir.mkdir(parents=True)
            reference_dir = data_root / "reference"
            reference_dir.mkdir(parents=True)
            (reference_dir / "datasets.xlsx").write_text("placeholder\n")

            (raw_dir / "GSM7300468_HC1_matrix.tar.gz").write_text("placeholder\n")

            samples = pipeline.discover_rna_samples_from_local_layout(
                raw_dir=data_root / "raw", selected_gses=["GSE231794"]
            )

            self.assertEqual(len(samples), 1)
            self.assertEqual(samples[0].sample_id, "GSM7300468")
            self.assertEqual(samples[0].input_type, "archive")

    def test_marks_shared_gene_count_csv_as_unsupported(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            data_root = Path(tmp)
            raw_dir = data_root / "raw" / "GSE198533"
            raw_dir.mkdir(parents=True)
            reference_dir = data_root / "reference"
            reference_dir.mkdir(parents=True)
            (reference_dir / "datasets.xlsx").write_text("placeholder\n")

            (raw_dir / "GSE198533_Raw_gene_counts_matrix.csv.gz").write_text("placeholder\n")

            samples = pipeline.discover_rna_samples_from_local_layout(
                raw_dir=data_root / "raw", selected_gses=["GSE198533"]
            )

            self.assertEqual(len(samples), 1)
            self.assertFalse(samples[0].supported)
            self.assertIn("unsupported", samples[0].note.lower())

    @staticmethod
    def _touch_triplet(directory: Path, stem: str, features_name: str = "features") -> None:
        (directory / f"{stem}_matrix.mtx.gz").write_text("matrix\n")
        (directory / f"{stem}_barcodes.tsv.gz").write_text("barcodes\n")
        (directory / f"{stem}_{features_name}.tsv.gz").write_text("features\n")

    @staticmethod
    def _touch_shared_triplet(directory: Path, stem: str) -> None:
        (directory / f"{stem}_matrix.mtx.gz").write_text("matrix\n")
        (directory / f"{stem}_barcodes.tsv.gz").write_text("barcodes\n")
        (directory / f"{stem}_features.tsv.gz").write_text("features\n")


if __name__ == "__main__":
    unittest.main()
