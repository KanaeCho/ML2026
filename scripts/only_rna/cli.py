from __future__ import annotations

import json
import sys
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable

from .annotation import annotate_with_all_versions
from .config import load_default_config
from .doublet import run_doublet_detection
from .discovery import (
    DiscoveredSample,
    discover_rna_samples,
    resolve_data_root,
    selected_rna_gses,
)
from .embedding import run_embedding
from .outputs import write_sample_outputs
from .qc import apply_qc_filters, compute_qc_metrics
from .read_inputs import read_sample_input


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = ROOT / "output" / "rna"
DEFAULT_CONFIG_PATH = ROOT / "scripts" / "only_rna" / "default_config.yaml"


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _resolve_selected_rna_gses(gse: str | None = None) -> tuple[Path, list[str]]:
    data_root = resolve_data_root(ROOT)
    selected = selected_rna_gses(data_root / "reference")
    if gse is not None:
        if gse not in selected:
            raise FileNotFoundError(
                f"GSE not selected in datasets.xlsx for scRNA: {gse}"
            )
        selected = [gse]
    return data_root, selected


def _discover_selected_rna_samples(gse: str | None = None) -> list[DiscoveredSample]:
    data_root, selected = _resolve_selected_rna_gses(gse)
    return discover_rna_samples(data_root / "raw", selected)


def find_rna_sample(gse: str, sample_id: str) -> DiscoveredSample:
    for sample in _discover_selected_rna_samples(gse=gse):
        if sample.sample_id == sample_id:
            return sample
    raise FileNotFoundError(f"RNA sample not found: {gse}/{sample_id}")


def _sample_output_dir(sample: DiscoveredSample, output_root: Path) -> Path:
    return Path(output_root) / sample.gse / sample.sample_id


def _sample_status_file(sample: DiscoveredSample, output_root: Path) -> Path:
    return _sample_output_dir(sample, output_root) / "run_status.json"


def _expected_output_paths(sample: DiscoveredSample, output_root: Path) -> list[Path]:
    output_dir = _sample_output_dir(sample, output_root)
    return [
        output_dir / "metadata.csv",
        output_dir / "metadata_qc.csv",
        output_dir / "qc_summary.csv",
        output_dir / "validation_result.csv",
        output_dir / f"{sample.sample_id}.h5ad",
        output_dir / "matrix" / "matrix.mtx",
        output_dir / "matrix" / "barcodes.tsv.gz",
        output_dir / "matrix" / "features.tsv.gz",
        output_dir / "umap_rna_clusters.png",
        output_dir / "umap_rna_cima_cell_type_l1.png",
        output_dir / "umap_rna_cima_cell_type_l2.png",
        output_dir / "umap_rna_cima_cell_type_l1_masked.png",
    ]


def _outputs_complete(sample: DiscoveredSample, output_root: Path) -> bool:
    return sample.supported and all(
        path.exists() for path in _expected_output_paths(sample, output_root)
    )


def _load_rna_status(sample: DiscoveredSample, output_root: Path) -> dict | None:
    status_file = _sample_status_file(sample, output_root)
    if not status_file.exists():
        return None
    return json.loads(status_file.read_text(encoding="utf-8"))


def _write_rna_status(
    sample: DiscoveredSample, payload: dict, output_root: Path
) -> None:
    output_dir = _sample_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    _sample_status_file(sample, output_root).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _primary_source_path(sample: DiscoveredSample) -> Path | None:
    return (
        sample.matrix_path
        or sample.h5_path
        or sample.archive_path
        or sample.features_path
        or sample.barcodes_path
    )


def _print_rna_discovery(samples: Iterable[DiscoveredSample]) -> None:
    print("gse\tsample_id\tinput_type\tsample_kind\tsupported\tnote\tprimary_path")
    for sample in samples:
        primary = _primary_source_path(sample)
        primary_text = str(primary) if primary is not None else ""
        print(
            f"{sample.gse}\t{sample.sample_id}\t{sample.input_type}\t"
            f"{sample.sample_kind}\t{str(sample.supported).lower()}\t{sample.note}\t{primary_text}"
        )


def _print_rna_status(samples: Iterable[DiscoveredSample], output_root: Path) -> None:
    print(
        "gse\tsample_id\tsample_kind\tstatus\toutputs_complete\tsupported\tlast_finished"
    )
    for sample in samples:
        payload = _load_rna_status(sample, output_root) or {}
        status = payload.get("status", "pending")
        payload_output_root = Path(payload.get("output_root", str(output_root)))
        finished_at = payload.get("finished_at", "")
        print(
            f"{sample.gse}\t{sample.sample_id}\t{sample.sample_kind}\t{status}\t"
            f"{str(_outputs_complete(sample, payload_output_root)).lower()}\t"
            f"{str(sample.supported).lower()}\t{finished_at}"
        )


def _build_placeholder_command(sample: DiscoveredSample, args) -> list[str]:
    python_bin = getattr(args, "python_bin", None) or sys.executable or "python3"
    return [
        python_bin,
        "-m",
        "scripts.only_rna",
        "run-sample",
        "--gse",
        sample.gse,
        "--sample-id",
        sample.sample_id,
        "--input-type",
        sample.input_type,
    ]


def _execute_rna_sample(sample: DiscoveredSample, args) -> bool:
    config = load_default_config(DEFAULT_CONFIG_PATH)
    adata = read_sample_input(sample)
    adata = compute_qc_metrics(adata, config)
    adata = run_doublet_detection(adata, config)
    adata = apply_qc_filters(adata, config)
    adata = run_embedding(adata, config)
    adata = annotate_with_all_versions(
        adata,
        reference_dir=resolve_data_root(ROOT) / "reference",
        azimuth_model_dir=None,
        celltypist_model_path=None,
        singler_model_path=None,
        scanvi_model_path=None,
        methods=["cima"],
    )
    write_sample_outputs(
        adata,
        output_root=Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT)),
        gse=sample.gse,
        sample_id=sample.sample_id,
        config=config,
    )
    return True


def _route_rna_sample(sample: DiscoveredSample, args) -> int:
    output_root = Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT))
    force = bool(getattr(args, "force", False))
    dry_run = bool(getattr(args, "dry_run", False))

    if not sample.supported:
        print(
            f"[skip-rna] {sample.gse}/{sample.sample_id} unsupported: {sample.note}",
            file=sys.stderr,
        )
        return 2

    if _outputs_complete(sample, output_root) and not force:
        print(
            f"[skip-rna] {sample.gse}/{sample.sample_id} already has complete RNA outputs"
        )
        return 0

    command = _build_placeholder_command(sample, args)
    started_at = utc_now()
    status_payload = {
        "sample": {"gse": sample.gse, "sample_id": sample.sample_id},
        "input_type": sample.input_type,
        "sample_kind": sample.sample_kind,
        "command": command,
        "output_root": str(output_root),
        "started_at": started_at,
        "status": "dry_run" if dry_run else "running",
        "outputs_complete": _outputs_complete(sample, output_root),
        "note": "Stable RNA CLI routing and execution.",
    }

    print(f"[run-rna] {sample.gse}/{sample.sample_id}")
    print(" ".join(command))
    if dry_run:
        status_payload["finished_at"] = started_at
        status_payload["status"] = "dry_run"
        _write_rna_status(sample, status_payload, output_root)
        return 0

    _write_rna_status(sample, status_payload, output_root)
    try:
        _execute_rna_sample(sample, args)
    except Exception as exc:
        status_payload["finished_at"] = utc_now()
        status_payload["status"] = "failed"
        status_payload["outputs_complete"] = _outputs_complete(sample, output_root)
        status_payload["note"] = str(exc)
        _write_rna_status(sample, status_payload, output_root)
        raise

    status_payload["finished_at"] = utc_now()
    status_payload["status"] = "success"
    status_payload["outputs_complete"] = _outputs_complete(sample, output_root)
    _write_rna_status(sample, status_payload, output_root)
    return 0


def cmd_discover_rna(args) -> int:
    _print_rna_discovery(_discover_selected_rna_samples(gse=getattr(args, "gse", None)))
    return 0


def cmd_run_rna_sample(args) -> int:
    sample = find_rna_sample(args.gse, args.sample_id)
    if sample.sample_kind == "gse_shared" or sample.sample_id == sample.gse:
        print(
            f"Refusing explicit run-rna-sample target for gse-level shared sample: {sample.gse}/{sample.sample_id}",
            file=sys.stderr,
        )
        return 2
    return _route_rna_sample(sample, args)


def cmd_run_rna_gse(args) -> int:
    samples = [
        sample
        for sample in _discover_selected_rna_samples(gse=args.gse)
        if sample.supported
    ]
    returncode = 0
    for sample in samples:
        current = _route_rna_sample(sample, args)
        if current != 0:
            returncode = current
    return returncode


def cmd_rna_status(args) -> int:
    _print_rna_status(
        _discover_selected_rna_samples(gse=getattr(args, "gse", None)),
        output_root=DEFAULT_OUTPUT_ROOT,
    )
    return 0


__all__ = [
    "cmd_discover_rna",
    "cmd_run_rna_sample",
    "cmd_run_rna_gse",
    "cmd_rna_status",
]
