from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.longevity import cli


def _write_gzip_lines(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def test_barcode_file_path_uses_sample_subdir(tmp_path: Path) -> None:
    path = cli.barcode_file_path("W202506130001022", tmp_path)
    assert path == tmp_path / "W202506130001022" / "filtered_barcodes.tsv.gz"


def test_barcode_status_missing_empty_incomplete_ready(tmp_path: Path) -> None:
    sample = cli.LongevityAtacSample(
        dataset="longevity",
        sample_id="W1",
        fragment_file=tmp_path / "W1_fragments.tsv.gz",
    )
    assert cli._barcode_status_row(sample, tmp_path)["status"] == "missing"

    _write_gzip_lines(cli.barcode_file_path("W1", tmp_path), [])
    assert cli._barcode_status_row(sample, tmp_path)["status"] == "empty"

    _write_gzip_lines(cli.barcode_file_path("W1", tmp_path), ["AAAC-1"])
    assert cli._barcode_status_row(sample, tmp_path)["status"] == "incomplete"

    cli.barcode_qc_path("W1", tmp_path).write_text("barcode,n_fragments\nAAAC-1,100\n", encoding="utf-8")
    cli.barcode_summary_path("W1", tmp_path).write_text(
        json.dumps({"min_fragments": 200, "max_barcodes": 20000, "min_tss": 2.5}),
        encoding="utf-8",
    )
    row = cli._barcode_status_row(sample, tmp_path)
    assert row["status"] == "ready"
    assert row["n_barcodes"] == 1


def test_override_parsing_applies_sample_values(tmp_path: Path) -> None:
    overrides = tmp_path / "overrides.csv"
    pd.DataFrame(
        [
            {
                "sample_id": "W202506130001022",
                "min_fragments": "200",
                "max_barcodes": "20000",
                "min_tss": "2.5",
                "reason": "known long tail",
            }
        ]
    ).to_csv(overrides, index=False)
    args = SimpleNamespace(
        min_fragments=100,
        max_barcodes=1000,
        min_tss=1.5,
        overrides=overrides,
    )
    thresholds = cli._barcode_thresholds_for_sample("W202506130001022", args)
    assert thresholds["min_fragments"] == 200
    assert thresholds["max_barcodes"] == 20000
    assert thresholds["min_tss"] == 2.5
    assert thresholds["override_reason"] == "known long tail"


def test_run_command_includes_generated_barcode_file(monkeypatch, tmp_path: Path) -> None:
    sample = cli.LongevityAtacSample(
        dataset="longevity",
        sample_id="W1",
        fragment_file=tmp_path / "W1_fragments.tsv.gz",
    )
    _write_gzip_lines(cli.barcode_file_path("W1", tmp_path), ["AAAC-1"])
    monkeypatch.setattr(cli, "barcode_file_path", lambda sample_id, output_root=None: tmp_path / sample_id / "filtered_barcodes.tsv.gz")
    args = SimpleNamespace(
        output_root=tmp_path / "out",
        output_profile="full",
        force=True,
        rscript="Rscript",
        nmads=4,
        min_inferred_fragments=None,
        max_inferred_barcodes=None,
        umap_min_dist=None,
        barcode_output_root=tmp_path,
        barcode_min_fragments=200,
        barcode_max_barcodes=20000,
        barcode_min_tss=2.5,
        barcode_overrides=None,
        dry_run=True,
    )
    assert cli._run_atac_sample(sample, args) == 0
    status = json.loads((tmp_path / "out" / "atac" / "longevity" / "W1" / "run_status.json").read_text(encoding="utf-8"))
    assert "--barcode-file" in status["command"]


def test_run_command_without_generated_barcode_keeps_fallback(monkeypatch, tmp_path: Path) -> None:
    sample = cli.LongevityAtacSample(
        dataset="longevity",
        sample_id="W1",
        fragment_file=tmp_path / "W1_fragments.tsv.gz",
    )
    preprocess_calls = []
    monkeypatch.setattr(cli, "barcode_file_path", lambda sample_id, output_root=None: tmp_path / sample_id / "filtered_barcodes.tsv.gz")
    monkeypatch.setattr(
        cli,
        "_preprocess_atac_barcodes_sample",
        lambda sample, args: preprocess_calls.append((sample.sample_id, args.output_root)) or 0,
    )
    args = SimpleNamespace(
        output_root=tmp_path / "out",
        output_profile="full",
        force=True,
        rscript="Rscript",
        nmads=4,
        min_inferred_fragments=None,
        max_inferred_barcodes=None,
        umap_min_dist=None,
        barcode_output_root=tmp_path,
        barcode_min_fragments=200,
        barcode_max_barcodes=20000,
        barcode_min_tss=2.5,
        barcode_overrides=None,
        dry_run=True,
    )
    assert cli._run_atac_sample(sample, args) == 0
    status = json.loads((tmp_path / "out" / "atac" / "longevity" / "W1" / "run_status.json").read_text(encoding="utf-8"))
    assert "--barcode-file" not in status["command"]
    assert preprocess_calls == [("W1", tmp_path)]


def test_run_command_preprocesses_missing_barcode_before_atac(monkeypatch, tmp_path: Path) -> None:
    sample = cli.LongevityAtacSample(
        dataset="longevity",
        sample_id="W1",
        fragment_file=tmp_path / "W1_fragments.tsv.gz",
    )

    def fake_preprocess(sample, args):
        _write_gzip_lines(cli.barcode_file_path(sample.sample_id, tmp_path), ["AAAC-1"])
        cli.barcode_qc_path(sample.sample_id, tmp_path).write_text("barcode,n_fragments\nAAAC-1,100\n", encoding="utf-8")
        cli.barcode_summary_path(sample.sample_id, tmp_path).write_text(
            json.dumps({"min_fragments": 200, "max_barcodes": 20000, "min_tss": 2.5}),
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(cli, "_preprocess_atac_barcodes_sample", fake_preprocess)
    args = SimpleNamespace(
        output_root=tmp_path / "out",
        output_profile="full",
        force=True,
        rscript="Rscript",
        nmads=4,
        min_inferred_fragments=None,
        max_inferred_barcodes=None,
        umap_min_dist=None,
        barcode_output_root=tmp_path,
        barcode_min_fragments=200,
        barcode_max_barcodes=20000,
        barcode_min_tss=2.5,
        barcode_overrides=None,
        dry_run=True,
    )
    assert cli._run_atac_sample(sample, args) == 0
    status = json.loads((tmp_path / "out" / "atac" / "longevity" / "W1" / "run_status.json").read_text(encoding="utf-8"))
    assert "--barcode-file" in status["command"]


def test_ready_barcode_with_mismatched_thresholds_is_reprocessed(monkeypatch, tmp_path: Path) -> None:
    sample = cli.LongevityAtacSample(
        dataset="longevity",
        sample_id="W1",
        fragment_file=tmp_path / "W1_fragments.tsv.gz",
    )
    _write_gzip_lines(cli.barcode_file_path("W1", tmp_path), ["OLD-1"])
    cli.barcode_qc_path("W1", tmp_path).write_text("barcode,n_fragments\nOLD-1,100\n", encoding="utf-8")
    cli.barcode_summary_path("W1", tmp_path).write_text(
        json.dumps({"min_fragments": 200, "max_barcodes": 50000, "min_tss": 0}),
        encoding="utf-8",
    )
    calls = []

    def fake_preprocess(sample, args):
        calls.append((args.max_barcodes, args.min_tss))
        _write_gzip_lines(cli.barcode_file_path(sample.sample_id, tmp_path), ["NEW-1"])
        cli.barcode_qc_path(sample.sample_id, tmp_path).write_text("barcode,n_fragments\nNEW-1,100\n", encoding="utf-8")
        cli.barcode_summary_path(sample.sample_id, tmp_path).write_text(
            json.dumps({"min_fragments": 200, "max_barcodes": 20000, "min_tss": 2.5}),
            encoding="utf-8",
        )
        return 0

    monkeypatch.setattr(cli, "_preprocess_atac_barcodes_sample", fake_preprocess)
    args = SimpleNamespace(
        output_root=tmp_path / "out",
        output_profile="full",
        force=True,
        rscript="Rscript",
        nmads=4,
        min_inferred_fragments=None,
        max_inferred_barcodes=None,
        umap_min_dist=None,
        barcode_output_root=tmp_path,
        barcode_min_fragments=200,
        barcode_max_barcodes=20000,
        barcode_min_tss=2.5,
        barcode_overrides=None,
        dry_run=True,
    )
    assert cli._run_atac_sample(sample, args) == 0
    assert calls == [(20000, 2.5)]
