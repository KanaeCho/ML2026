#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = ROOT / "data" / "raw"
OUTPUT_DIR = ROOT / "output"
R_SCRIPT = ROOT / "scripts" / "process" / "process_single_sample.R"
DOWNLOAD_SCRIPT = ROOT / "scripts" / "process" / "download_from_datasets.py"
FRAGMENT_RE = re.compile(r"^(GSM\d+)_.*fragments.*\.tsv\.gz$")


@dataclass(frozen=True)
class Sample:
    gse: str
    gsm: str
    fragment_file: Path

    @property
    def output_dir(self) -> Path:
        return OUTPUT_DIR / self.gse / self.gsm

    @property
    def log_file(self) -> Path:
        return self.output_dir / "logs" / "sample_qc.log"

    @property
    def status_file(self) -> Path:
        return self.output_dir / "run_status.json"


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def discover_samples(gse: str | None = None) -> list[Sample]:
    if not RAW_DIR.exists():
        raise FileNotFoundError(f"Raw data directory not found: {RAW_DIR}")

    gse_dirs = (
        [RAW_DIR / gse]
        if gse
        else sorted(path for path in RAW_DIR.iterdir() if path.is_dir())
    )
    samples: list[Sample] = []

    for gse_dir in gse_dirs:
        if not gse_dir.exists():
            raise FileNotFoundError(f"GSE directory not found: {gse_dir}")

        by_gsm: dict[str, list[Path]] = {}
        for path in sorted(gse_dir.glob("*.tsv.gz")):
            match = FRAGMENT_RE.match(path.name)
            if not match:
                continue
            gsm = match.group(1)
            by_gsm.setdefault(gsm, []).append(path)

        for gsm, matches in sorted(by_gsm.items()):
            if len(matches) != 1:
                names = ", ".join(path.name for path in matches)
                raise RuntimeError(
                    f"Multiple fragment files found for {gse_dir.name}/{gsm}: {names}"
                )
            samples.append(Sample(gse=gse_dir.name, gsm=gsm, fragment_file=matches[0]))

    return samples


def find_sample(gse: str, gsm: str) -> Sample:
    for sample in discover_samples(gse=gse):
        if sample.gsm == gsm:
            return sample
    raise FileNotFoundError(f"Sample not found: {gse}/{gsm}")


def sample_output_dir(sample: Sample, output_root: Path) -> Path:
    return output_root / sample.gse / sample.gsm


def sample_log_file(sample: Sample, output_root: Path) -> Path:
    return sample_output_dir(sample, output_root) / "logs" / "sample_qc.log"


def sample_status_file(sample: Sample, output_root: Path) -> Path:
    return sample_output_dir(sample, output_root) / "run_status.json"


def expected_outputs(
    sample: Sample, output_profile: str = "full", output_root: Path = OUTPUT_DIR
) -> list[Path]:
    output_dir = sample_output_dir(sample, output_root)
    if output_profile in {"matrix-lite", "validation-lite"}:
        return [
            output_dir / "umap_cima_cell_type_l1.png",
            output_dir / "qc_summary.csv",
            output_dir / "validation_result.csv",
            output_dir / "matrix" / "matrix.mtx",
            output_dir / "matrix" / "barcodes.tsv",
            output_dir / "matrix" / "features.tsv",
        ]

    return [
        output_dir / "qc_overview.png",
        output_dir / "umap_cima_cell_type_l1.png",
        output_dir / "umap_cima_cell_type_l2.png",
        output_dir / "umap_cima_cell_type_l3.png",
        output_dir / "umap_cima_cell_type_l4.png",
        output_dir / "metadata.csv",
        output_dir / "metadata_qc.csv",
        output_dir / "qc_summary.csv",
        output_dir / "matrix" / "matrix.mtx",
        output_dir / "matrix" / "barcodes.tsv.gz",
        output_dir / "matrix" / "features.tsv.gz",
        output_dir / f"{sample.gsm}_seurat_qc.rds",
    ]


def outputs_complete(
    sample: Sample, output_profile: str = "full", output_root: Path = OUTPUT_DIR
) -> bool:
    return all(
        path.exists()
        for path in expected_outputs(
            sample, output_profile=output_profile, output_root=output_root
        )
    )


def load_status(sample: Sample, output_root: Path = OUTPUT_DIR) -> dict | None:
    status_file = sample_status_file(sample, output_root)
    if not status_file.exists():
        return None
    return json.loads(status_file.read_text())


def write_status(sample: Sample, payload: dict, output_root: Path = OUTPUT_DIR) -> None:
    output_dir = sample_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_status_file(sample, output_root).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    )


def run_command(command: list[str], log_file: Path) -> int:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert process.stdout is not None
    with log_file.open("w", encoding="utf-8") as handle:
        for line in process.stdout:
            sys.stdout.write(line)
            handle.write(line)
            handle.flush()

    return process.wait()


def run_download(
    python_bin: str,
    assay: str,
    data_format: str,
    gse: str,
    gsm: str,
    file_kinds: str,
    aria2c: str,
    dry_run: bool,
    skip_network_resolve: bool,
    max_concurrent_downloads: int,
    split: int,
    network_timeout: int,
    links_out: str,
    print_links: bool,
) -> int:
    log_dir = RAW_DIR / "_download_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"download_{assay}_{data_format}.log"

    command = [
        python_bin,
        str(DOWNLOAD_SCRIPT),
        "--assay",
        assay,
        "--data-format",
        data_format,
        "--file-kinds",
        file_kinds,
        "--aria2c",
        aria2c,
        "--max-concurrent-downloads",
        str(max_concurrent_downloads),
        "--split",
        str(split),
        "--network-timeout",
        str(network_timeout),
    ]

    if gse:
        command.extend(["--gse", gse])
    if gsm:
        command.extend(["--gsm", gsm])

    if links_out:
        command.extend(["--links-out", links_out])
    if print_links:
        command.append("--print-links")

    if dry_run:
        command.append("--dry-run")
    if skip_network_resolve:
        command.append("--skip-network-resolve")

    print("[download]")
    print(" ".join(command))
    return run_command(command, log_file)


def run_sample(
    sample: Sample,
    rscript: str,
    nmads: float,
    output_profile: str,
    output_root: Path,
    force: bool,
    dry_run: bool,
) -> int:
    if (
        outputs_complete(sample, output_profile=output_profile, output_root=output_root)
        and not force
    ):
        print(f"[skip] {sample.gse}/{sample.gsm} already has complete outputs")
        return 0

    command = [
        rscript,
        str(R_SCRIPT),
        "--gse",
        sample.gse,
        "--gsm",
        sample.gsm,
        "--nmads",
        str(nmads),
        "--output-profile",
        output_profile,
        "--output-root",
        str(output_root),
    ]

    status_payload = {
        "sample": {"gse": sample.gse, "gsm": sample.gsm},
        "command": command,
        "output_profile": output_profile,
        "output_root": str(output_root),
        "started_at": utc_now(),
        "status": "running",
    }
    write_status(sample, status_payload, output_root=output_root)

    print(f"[run] {sample.gse}/{sample.gsm}")
    print(" ".join(command))

    if dry_run:
        status_payload["status"] = "dry_run"
        status_payload["finished_at"] = utc_now()
        write_status(sample, status_payload, output_root=output_root)
        return 0

    returncode = run_command(command, sample_log_file(sample, output_root))
    status_payload["finished_at"] = utc_now()
    status_payload["returncode"] = returncode
    status_payload["outputs_complete"] = outputs_complete(
        sample, output_profile=output_profile, output_root=output_root
    )
    status_payload["status"] = (
        "success"
        if returncode == 0
        and outputs_complete(
            sample, output_profile=output_profile, output_root=output_root
        )
        else "failed"
    )
    write_status(sample, status_payload, output_root=output_root)

    return returncode


def print_discovery(samples: Iterable[Sample]) -> None:
    print("gse\tgsm\tfragment_file")
    for sample in samples:
        print(f"{sample.gse}\t{sample.gsm}\t{sample.fragment_file.relative_to(ROOT)}")


def print_status(samples: Iterable[Sample]) -> None:
    print("gse\tgsm\tstatus\toutputs_complete\tlast_finished")
    for sample in samples:
        payload = load_status(sample) or {}
        status = payload.get("status", "pending")
        output_profile = payload.get("output_profile", "full")
        output_root = Path(payload.get("output_root", str(OUTPUT_DIR)))
        finished_at = payload.get("finished_at", "")
        print(
            f"{sample.gse}\t{sample.gsm}\t{status}\t"
            f"{str(outputs_complete(sample, output_profile=output_profile, output_root=output_root)).lower()}\t{finished_at}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage scATAC single-sample QC runs")
    subparsers = parser.add_subparsers(dest="command", required=True)

    discover_parser = subparsers.add_parser("discover", help="List discovered samples")
    discover_parser.add_argument("--gse", help="Restrict discovery to one GSE")

    run_sample_parser = subparsers.add_parser("run-sample", help="Run one sample")
    run_sample_parser.add_argument("--gse", required=True)
    run_sample_parser.add_argument("--gsm", required=True)
    run_sample_parser.add_argument("--nmads", type=float, default=4)
    run_sample_parser.add_argument("--rscript", default="Rscript")
    run_sample_parser.add_argument("--output-profile", default="full")
    run_sample_parser.add_argument("--output-root", default=str(OUTPUT_DIR))
    run_sample_parser.add_argument("--force", action="store_true")
    run_sample_parser.add_argument("--dry-run", action="store_true")

    run_gse_parser = subparsers.add_parser(
        "run-gse", help="Run all samples under one GSE"
    )
    run_gse_parser.add_argument("--gse", required=True)
    run_gse_parser.add_argument("--nmads", type=float, default=4)
    run_gse_parser.add_argument("--rscript", default="Rscript")
    run_gse_parser.add_argument("--output-profile", default="full")
    run_gse_parser.add_argument("--output-root", default=str(OUTPUT_DIR))
    run_gse_parser.add_argument("--force", action="store_true")
    run_gse_parser.add_argument("--dry-run", action="store_true")

    download_parser = subparsers.add_parser(
        "download", help="Download GEO files based on datasets.xlsx filtering"
    )
    download_parser.add_argument("--assay", default="scATAC")
    download_parser.add_argument("--data-format", default="fragment")
    download_parser.add_argument("--gse", default="")
    download_parser.add_argument("--gsm", default="")
    download_parser.add_argument("--file-kinds", default="fragment,barcode,singlecell")
    download_parser.add_argument("--aria2c", default="aria2c")
    download_parser.add_argument("--python-bin", default=sys.executable or "python3")
    download_parser.add_argument("--dry-run", action="store_true")
    download_parser.add_argument("--skip-network-resolve", action="store_true")
    download_parser.add_argument("--max-concurrent-downloads", type=int, default=4)
    download_parser.add_argument("--split", type=int, default=8)
    download_parser.add_argument("--network-timeout", type=int, default=10)
    download_parser.add_argument("--links-out", default="")
    download_parser.add_argument("--print-links", action="store_true")

    status_parser = subparsers.add_parser("status", help="Show run status")
    status_parser.add_argument("--gse", help="Restrict status to one GSE")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "discover":
        print_discovery(discover_samples(gse=args.gse))
        return 0

    if args.command == "status":
        print_status(discover_samples(gse=args.gse))
        return 0

    if args.command == "run-sample":
        sample = find_sample(args.gse, args.gsm)
        return run_sample(
            sample=sample,
            rscript=args.rscript,
            nmads=args.nmads,
            output_profile=args.output_profile,
            output_root=Path(args.output_root),
            force=args.force,
            dry_run=args.dry_run,
        )

    if args.command == "run-gse":
        samples = discover_samples(gse=args.gse)
        returncode = 0
        for sample in samples:
            current = run_sample(
                sample=sample,
                rscript=args.rscript,
                nmads=args.nmads,
                output_profile=args.output_profile,
                output_root=Path(args.output_root),
                force=args.force,
                dry_run=args.dry_run,
            )
            if current != 0:
                returncode = current
        return returncode

    if args.command == "download":
        return run_download(
            python_bin=args.python_bin,
            assay=args.assay,
            data_format=args.data_format,
            gse=args.gse,
            gsm=args.gsm,
            file_kinds=args.file_kinds,
            aria2c=args.aria2c,
            dry_run=args.dry_run,
            skip_network_resolve=args.skip_network_resolve,
            max_concurrent_downloads=args.max_concurrent_downloads,
            split=args.split,
            network_timeout=args.network_timeout,
            links_out=args.links_out,
            print_links=args.print_links,
        )

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
