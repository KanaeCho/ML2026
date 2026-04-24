from __future__ import annotations

from pathlib import Path

from scripts.only_rna.discovery import discover_rna_samples, resolve_data_root


def test_resolve_data_root_prefers_workspace_data(tmp_path: Path, monkeypatch) -> None:
    workspace_data = tmp_path / "data"
    workspace_data.mkdir()
    env_root = tmp_path / "env_data"
    env_root.mkdir()
    monkeypatch.setenv("ML2026_DATA_ROOT", str(env_root))

    resolved = resolve_data_root(tmp_path)

    assert resolved == workspace_data


def test_resolve_data_root_falls_back_to_env(tmp_path: Path, monkeypatch) -> None:
    env_root = tmp_path / "env_data"
    env_root.mkdir()
    monkeypatch.setenv("ML2026_DATA_ROOT", str(env_root))

    resolved = resolve_data_root(tmp_path)

    assert resolved == env_root


def test_discovers_per_gsm_triplet_samples(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw" / "GSE167363"
    raw_dir.mkdir(parents=True)

    _touch_triplet(raw_dir, "GSM5102900_HC1")
    _touch_triplet(raw_dir, "GSM5102901_HC2")

    samples = discover_rna_samples(tmp_path / "raw", ["GSE167363"])

    assert [sample.sample_id for sample in samples] == ["GSM5102900", "GSM5102901"]
    assert {sample.input_type for sample in samples} == {"triplet"}
    assert {sample.sample_kind for sample in samples} == {"gsm"}
    assert all(sample.supported for sample in samples)


def test_prefers_pbmc_files_for_gse226039(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw" / "GSE226039"
    raw_dir.mkdir(parents=True)

    _touch_triplet(raw_dir, "GSM7061877_H44_Ileum", features_name="genes")
    _touch_triplet(raw_dir, "GSM7061878_H44_PBMC", features_name="genes")
    _touch_triplet(raw_dir, "GSM7061879_H44_Rectum", features_name="genes")

    samples = discover_rna_samples(tmp_path / "raw", ["GSE226039"])

    assert len(samples) == 1
    assert samples[0].sample_id == "GSM7061878"
    assert samples[0].sample_kind == "gsm"
    assert "PBMC" in samples[0].source_name.upper()


def test_discovers_shared_gse_triplet_as_supported_gse_level_sample(
    tmp_path: Path,
) -> None:
    raw_dir = tmp_path / "raw" / "GSE149689"
    raw_dir.mkdir(parents=True)

    _touch_shared_triplet(raw_dir, "GSE149689")

    samples = discover_rna_samples(tmp_path / "raw", ["GSE149689"])

    assert len(samples) == 1
    assert samples[0].sample_id == "GSE149689"
    assert samples[0].input_type == "triplet"
    assert samples[0].sample_kind == "gse_shared"
    assert samples[0].supported is True


def test_marks_shared_gene_count_csv_as_unsupported(tmp_path: Path) -> None:
    raw_dir = tmp_path / "raw" / "GSE198533"
    raw_dir.mkdir(parents=True)

    gene_count_path = raw_dir / "GSE198533_Raw_gene_counts_matrix.csv.gz"
    gene_count_path.write_text("placeholder\n")

    samples = discover_rna_samples(tmp_path / "raw", ["GSE198533"])

    assert len(samples) == 1
    assert samples[0].sample_id == "GSE198533"
    assert samples[0].sample_kind == "gse_shared"
    assert samples[0].supported is False
    assert (
        "gene-count" in samples[0].note.lower()
        or "unsupported" in samples[0].note.lower()
    )
    assert samples[0].source_name == gene_count_path.name


def _touch_triplet(directory: Path, stem: str, features_name: str = "features") -> None:
    (directory / f"{stem}_matrix.mtx.gz").write_text("matrix\n")
    (directory / f"{stem}_barcodes.tsv.gz").write_text("barcodes\n")
    (directory / f"{stem}_{features_name}.tsv.gz").write_text("features\n")


def _touch_shared_triplet(directory: Path, stem: str) -> None:
    (directory / f"{stem}_matrix.mtx.gz").write_text("matrix\n")
    (directory / f"{stem}_barcodes.tsv.gz").write_text("barcodes\n")
    (directory / f"{stem}_features.tsv.gz").write_text("features\n")
