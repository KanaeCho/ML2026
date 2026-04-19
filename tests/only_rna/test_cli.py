from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace

from scripts.only_rna import cli
from scripts.only_rna.discovery import DiscoveredSample
from scripts.process import pipeline


def test_cmd_run_rna_sample_rejects_explicit_gse_shared_target(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    shared_sample = _make_sample(
        gse="GSE149689",
        sample_id="GSE149689",
        sample_kind="gse_shared",
    )
    monkeypatch.setattr(cli, "find_rna_sample", lambda gse, sample_id: shared_sample)

    args = SimpleNamespace(
        gse="GSE149689",
        sample_id="GSE149689",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=False,
        dry_run=True,
    )

    returncode = cli.cmd_run_rna_sample(args)

    captured = capsys.readouterr()
    assert returncode != 0
    assert "gse-level shared sample" in captured.err.lower()


def test_cmd_run_rna_gse_includes_supported_shared_triplet_when_present(
    tmp_path: Path, monkeypatch
) -> None:
    seen: list[tuple[str, str]] = []
    samples = [
        _make_sample(gse="GSE149689", sample_id="GSM000001", sample_kind="gsm"),
        _make_sample(
            gse="GSE149689",
            sample_id="GSE149689",
            sample_kind="gse_shared",
        ),
        _make_sample(
            gse="GSE149689",
            sample_id="GSE149689-unsupported",
            sample_kind="gse_shared",
            supported=False,
        ),
    ]
    monkeypatch.setattr(cli, "_discover_selected_rna_samples", lambda gse=None: samples)
    monkeypatch.setattr(
        cli,
        "_route_rna_sample",
        lambda sample, args: seen.append((sample.sample_id, sample.sample_kind)) or 0,
    )

    args = SimpleNamespace(
        gse="GSE149689",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=False,
        dry_run=True,
    )

    returncode = cli.cmd_run_rna_gse(args)

    assert returncode == 0
    assert seen == [("GSM000001", "gsm"), ("GSE149689", "gse_shared")]


def test_cmd_rna_status_uses_new_output_family_with_h5ad(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    complete_sample = _make_sample(gse="GSE123456", sample_id="GSM000001")
    missing_h5ad_sample = _make_sample(gse="GSE123456", sample_id="GSM000002")
    monkeypatch.setattr(
        cli,
        "_discover_selected_rna_samples",
        lambda gse=None: [complete_sample, missing_h5ad_sample],
    )
    monkeypatch.setattr(cli, "DEFAULT_OUTPUT_ROOT", tmp_path)

    _touch_output_family(tmp_path, complete_sample, include_h5ad=True)
    _touch_output_family(tmp_path, missing_h5ad_sample, include_h5ad=False)

    args = SimpleNamespace(gse="GSE123456")

    returncode = cli.cmd_rna_status(args)

    assert returncode == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("gse\tsample_id\tsample_kind\tstatus\toutputs_complete")
    assert any(
        line.startswith("GSE123456\tGSM000001\tgsm\tsuccess\ttrue")
        for line in lines[1:]
    )
    assert any(
        line.startswith("GSE123456\tGSM000002\tgsm\tsuccess\tfalse")
        for line in lines[1:]
    )


def test_pipeline_discover_rna_routes_to_only_rna_cli(monkeypatch) -> None:
    seen: dict[str, str | None] = {}

    def fake_handler(args) -> int:
        seen["command"] = args.command
        seen["gse"] = args.gse
        return 7

    monkeypatch.setattr(pipeline.only_rna_cli, "cmd_discover_rna", fake_handler)
    monkeypatch.setattr(
        sys,
        "argv",
        ["pipeline.py", "discover-rna", "--gse", "GSE167363"],
    )

    returncode = pipeline.main()

    assert returncode == 7
    assert seen == {"command": "discover-rna", "gse": "GSE167363"}


def test_pipeline_tune_rna_sample_routes_to_only_rna_cli(monkeypatch) -> None:
    seen: dict[str, str | None] = {}

    def fake_handler(args) -> int:
        seen["command"] = args.command
        seen["gse"] = args.gse
        seen["sample_id"] = args.sample_id
        return 9

    monkeypatch.setattr(pipeline.only_rna_cli, "cmd_tune_rna_sample", fake_handler)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pipeline.py",
            "tune-rna-sample",
            "--gse",
            "GSE167363",
            "--sample-id",
            "GSM5102900",
        ],
    )

    returncode = pipeline.main()

    assert returncode == 9
    assert seen == {
        "command": "tune-rna-sample",
        "gse": "GSE167363",
        "sample_id": "GSM5102900",
    }


def test_pipeline_tune_rna_gse_routes_to_only_rna_cli(monkeypatch) -> None:
    seen: dict[str, str | None] = {}

    def fake_handler(args) -> int:
        seen["command"] = args.command
        seen["gse"] = args.gse
        return 11

    monkeypatch.setattr(pipeline.only_rna_cli, "cmd_tune_rna_gse", fake_handler)
    monkeypatch.setattr(
        sys,
        "argv",
        ["pipeline.py", "tune-rna-gse", "--gse", "GSE167363"],
    )

    returncode = pipeline.main()

    assert returncode == 11
    assert seen == {"command": "tune-rna-gse", "gse": "GSE167363"}


def test_cmd_discover_rna_prints_discovered_samples(monkeypatch, capsys) -> None:
    discovered = [
        _make_sample(
            gse="GSE167363",
            sample_id="GSM5102900",
            input_type="triplet",
            sample_kind="gsm",
        )
    ]
    monkeypatch.setattr(
        cli, "_discover_selected_rna_samples", lambda gse=None: discovered
    )

    returncode = cli.cmd_discover_rna(SimpleNamespace(gse="GSE167363"))

    assert returncode == 0
    lines = capsys.readouterr().out.strip().splitlines()
    assert lines[0].startswith("gse\tsample_id\tinput_type\tsample_kind")
    assert lines[1].startswith("GSE167363\tGSM5102900\ttriplet\tgsm\ttrue")


def test_cmd_run_rna_sample_executes_processing_chain_when_not_dry_run(
    tmp_path: Path, monkeypatch
) -> None:
    sample = _make_sample(gse="GSE167363", sample_id="GSM5102900")
    calls: list[str] = []

    monkeypatch.setattr(cli, "find_rna_sample", lambda gse, sample_id: sample)

    output_dir = tmp_path / "GSE167363" / "GSM5102900"

    def fake_execute(current_sample, args):
        calls.append(f"execute:{current_sample.sample_id}:{args.output_root}")
        _touch_output_family(tmp_path, current_sample, include_h5ad=True)
        return True

    monkeypatch.setattr(cli, "_execute_rna_sample", fake_execute)

    args = SimpleNamespace(
        gse="GSE167363",
        sample_id="GSM5102900",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=True,
        dry_run=False,
    )

    returncode = cli.cmd_run_rna_sample(args)

    assert returncode == 0
    assert calls == [f"execute:GSM5102900:{tmp_path}"]
    status = json.loads((output_dir / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] == "success"
    assert status["outputs_complete"] is True


def test_cmd_tune_rna_sample_executes_bounded_tuning_when_not_dry_run(
    tmp_path: Path, monkeypatch
) -> None:
    sample = _make_sample(gse="GSE167363", sample_id="GSM5102900")
    calls: list[str] = []

    monkeypatch.setattr(cli, "find_rna_sample", lambda gse, sample_id: sample)

    output_dir = tmp_path / "GSE167363" / "GSM5102900"

    def fake_execute(current_sample, args):
        calls.append(f"tune:{current_sample.sample_id}:{args.output_root}")
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "run_status.json").write_text(
            json.dumps(
                {
                    "status": "success",
                    "finished_at": "2026-04-19T00:00:00Z",
                    "output_root": str(tmp_path),
                }
            )
            + "\n",
            encoding="utf-8",
        )
        return True

    monkeypatch.setattr(cli, "_execute_tune_rna_sample", fake_execute)

    args = SimpleNamespace(
        gse="GSE167363",
        sample_id="GSM5102900",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=True,
        dry_run=False,
    )

    returncode = cli.cmd_tune_rna_sample(args)

    assert returncode == 0
    assert calls == [f"tune:GSM5102900:{tmp_path}"]
    status = json.loads((output_dir / "run_status.json").read_text(encoding="utf-8"))
    assert status["status"] == "success"


def test_cmd_tune_rna_sample_dry_run_records_tuning_mode_and_flags(
    tmp_path: Path, monkeypatch
) -> None:
    sample = _make_sample(gse="GSE167363", sample_id="GSM5102900")
    monkeypatch.setattr(cli, "find_rna_sample", lambda gse, sample_id: sample)

    args = SimpleNamespace(
        gse="GSE167363",
        sample_id="GSM5102900",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=True,
        dry_run=True,
    )

    returncode = cli.cmd_tune_rna_sample(args)

    assert returncode == 0
    status = json.loads(
        (tmp_path / "GSE167363" / "GSM5102900" / "run_status.json").read_text(
            encoding="utf-8"
        )
    )
    assert status["status"] == "dry_run"
    assert status["mode"] == "tuning"
    assert status["command"][2] == "tune-rna-sample"
    assert "--gse" in status["command"]
    assert "--sample-id" in status["command"]
    assert "--output-root" in status["command"]
    assert "--force" in status["command"]
    assert "--dry-run" in status["command"]


def test_cmd_tune_rna_gse_includes_supported_shared_triplet_when_present(
    tmp_path: Path, monkeypatch
) -> None:
    seen: list[tuple[str, str]] = []
    samples = [
        _make_sample(gse="GSE149689", sample_id="GSM000001", sample_kind="gsm"),
        _make_sample(
            gse="GSE149689",
            sample_id="GSE149689",
            sample_kind="gse_shared",
        ),
        _make_sample(
            gse="GSE149689",
            sample_id="GSE149689-unsupported",
            sample_kind="gse_shared",
            supported=False,
        ),
    ]
    monkeypatch.setattr(cli, "_discover_selected_rna_samples", lambda gse=None: samples)
    monkeypatch.setattr(
        cli,
        "_route_tune_rna_sample",
        lambda sample, args: seen.append((sample.sample_id, sample.sample_kind)) or 0,
    )

    args = SimpleNamespace(
        gse="GSE149689",
        output_root=str(tmp_path),
        python_bin=sys.executable,
        force=False,
        dry_run=True,
    )

    returncode = cli.cmd_tune_rna_gse(args)

    assert returncode == 0
    assert seen == [("GSM000001", "gsm"), ("GSE149689", "gse_shared")]


def _make_sample(
    *,
    gse: str,
    sample_id: str,
    input_type: str = "triplet",
    sample_kind: str = "gsm",
    supported: bool = True,
) -> DiscoveredSample:
    return DiscoveredSample(
        gse=gse,
        sample_id=sample_id,
        input_type=input_type,
        sample_kind=sample_kind,
        supported=supported,
        note="fixture",
        source_name=f"{sample_id}_matrix.mtx.gz",
        matrix_path=Path(f"/tmp/{sample_id}_matrix.mtx.gz"),
        barcodes_path=Path(f"/tmp/{sample_id}_barcodes.tsv.gz"),
        features_path=Path(f"/tmp/{sample_id}_features.tsv.gz"),
    )


def _touch_output_family(
    root: Path, sample: DiscoveredSample, *, include_h5ad: bool
) -> None:
    sample_dir = root / sample.gse / sample.sample_id
    relpaths = [
        "metadata.csv",
        "metadata_qc.csv",
        "qc_summary.csv",
        "validation_result.csv",
        "matrix/matrix.mtx",
        "matrix/barcodes.tsv.gz",
        "matrix/features.tsv.gz",
        "umap_rna_clusters.png",
        "umap_rna_cima_cell_type_l1.png",
        "umap_rna_cima_cell_type_l2.png",
        "umap_rna_cima_cell_type_l1_masked.png",
    ]
    if include_h5ad:
        relpaths.append(f"{sample.sample_id}.h5ad")

    for relpath in relpaths:
        path = sample_dir / relpath
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")

    (sample_dir / "run_status.json").write_text(
        json.dumps(
            {
                "status": "success",
                "finished_at": "2026-04-14T00:00:00Z",
                "output_root": str(root),
            }
        )
        + "\n",
        encoding="utf-8",
    )
