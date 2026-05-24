#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from scripts.co import cli as co_cli
from scripts.longevity import cli as longevity_cli
from scripts.only_rna import backfill_final_celltype
from scripts.only_rna import cli as only_rna_cli
from scripts.process import organize_integrated_products

DEFAULT_DATA_ROOT = Path("/mnt/g/ML2026_data")


def resolve_data_root(
    project_root: Path = ROOT, fallback_candidates: Iterable[Path] | None = None
) -> Path:
    workspace_data = project_root / "data"
    if workspace_data.exists():
        return workspace_data

    candidates: list[Path] = []
    env_root = os.environ.get("ML2026_DATA_ROOT")
    if env_root:
        candidates.append(Path(env_root))
    if fallback_candidates is not None:
        candidates.extend(Path(candidate) for candidate in fallback_candidates)
    else:
        candidates.append(DEFAULT_DATA_ROOT)

    for candidate in candidates:
        if candidate.exists():
            return candidate

    searched = ", ".join(str(path) for path in [workspace_data, *candidates])
    raise FileNotFoundError(f"Unable to resolve data root. Searched: {searched}")


def default_data_root_path(project_root: Path = ROOT) -> Path:
    workspace_data = project_root / "data"
    if workspace_data.exists():
        return workspace_data
    env_root = os.environ.get("ML2026_DATA_ROOT")
    return Path(env_root) if env_root else DEFAULT_DATA_ROOT


DATA_ROOT = default_data_root_path()
RAW_DIR = DATA_ROOT / "raw"
OUTPUT_DIR = ROOT / "output"
ATAC_OUTPUT_DIR = OUTPUT_DIR / "atac"
R_SCRIPT = ROOT / "scripts" / "process" / "process_single_sample.R"
RNA_R_SCRIPT = ROOT / "scripts" / "process" / "process_single_rna_sample.R"
DOWNLOAD_SCRIPT = ROOT / "scripts" / "process" / "download_from_datasets.py"
TEA_SEQ_AUDIT_SCRIPT = ROOT / "scripts" / "process" / "organize_tea_seq_outputs.py"
EXPORT_ATAC_H5AD_SCRIPT = ROOT / "scripts" / "process" / "export_co_atac_h5ad.py"
FRAGMENT_RE = re.compile(r"^(GSM\d+)_.*fragments.*\.tsv\.gz$")
RNA_OUTPUT_DIR = OUTPUT_DIR / "rna"
ATAC_WORKBOOK = DATA_ROOT / "reference" / "atac.xlsx"


@dataclass(frozen=True)
class Sample:
    gse: str
    gsm: str
    fragment_file: Path
    individual_id: str = ""

    @property
    def output_dir(self) -> Path:
        return ATAC_OUTPUT_DIR / self.gse / self.gsm

    @property
    def log_file(self) -> Path:
        return self.output_dir / "logs" / "sample_qc.log"

    @property
    def status_file(self) -> Path:
        return self.output_dir / "run_status.json"


@dataclass(frozen=True)
class RNASample:
    gse: str
    sample_id: str
    input_type: str
    supported: bool
    note: str
    matrix_path: Path | None = None
    barcodes_path: Path | None = None
    features_path: Path | None = None
    h5_path: Path | None = None
    archive_path: Path | None = None

    @property
    def output_dir(self) -> Path:
        return RNA_OUTPUT_DIR / self.gse / self.sample_id

    @property
    def log_file(self) -> Path:
        return self.output_dir / "logs" / "sample_qc.log"

    @property
    def status_file(self) -> Path:
        return self.output_dir / "run_status.json"

    def command_args(self) -> list[str]:
        args = ["--input-type", self.input_type]
        if self.matrix_path is not None:
            args.extend(["--matrix-path", str(self.matrix_path)])
        if self.barcodes_path is not None:
            args.extend(["--barcodes-path", str(self.barcodes_path)])
        if self.features_path is not None:
            args.extend(["--features-path", str(self.features_path)])
        if self.h5_path is not None:
            args.extend(["--h5-path", str(self.h5_path)])
        if self.archive_path is not None:
            args.extend(["--archive-path", str(self.archive_path)])
        return args


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def discover_samples(gse: str | None = None) -> list[Sample]:
    if not RAW_DIR.exists():
        raise FileNotFoundError(f"Raw data directory not found: {RAW_DIR}")
    if not ATAC_WORKBOOK.exists():
        raise FileNotFoundError(f"ATAC workbook not found: {ATAC_WORKBOOK}")

    import pandas as pd

    rules = pd.read_excel(ATAC_WORKBOOK, dtype=str).fillna("")
    required_columns = {"sample", "dataset"}
    missing_columns = required_columns - set(rules.columns)
    if missing_columns:
        raise ValueError(
            f"ATAC workbook missing columns: {', '.join(sorted(missing_columns))}"
        )
    rules["sample"] = rules["sample"].astype(str).str.strip()
    rules["dataset"] = rules["dataset"].astype(str).str.strip()
    rules = rules[
        rules["sample"].str.fullmatch(r"GSM\d+", na=False)
        & rules["dataset"].str.fullmatch(r"GSE\d+", na=False)
    ].copy()
    if gse:
        rules = rules.loc[rules["dataset"] == gse].copy()
    selected = set(zip(rules["dataset"], rules["sample"], strict=False))
    individual_map: dict[tuple[str, str], str] = {}
    if "donor" in rules.columns:
        for _, row in rules.iterrows():
            donor = str(row.get("donor", "")).strip()
            if donor:
                individual_map[(str(row["dataset"]), str(row["sample"]))] = donor

    gse_dirs = (
        [RAW_DIR / gse]
        if gse
        else sorted(RAW_DIR / dataset for dataset in {key[0] for key in selected})
    )
    samples: list[Sample] = []

    for gse_dir in gse_dirs:
        if not gse_dir.exists():
            continue

        by_gsm: dict[str, list[Path]] = {}
        for path in sorted(gse_dir.glob("*.tsv.gz")):
            match = FRAGMENT_RE.match(path.name)
            if not match:
                continue
            gsm = match.group(1)
            if (gse_dir.name, gsm) not in selected:
                continue
            by_gsm.setdefault(gsm, []).append(path)

        for gsm, matches in sorted(by_gsm.items()):
            if len(matches) != 1:
                names = ", ".join(path.name for path in matches)
                raise RuntimeError(
                    f"Multiple fragment files found for {gse_dir.name}/{gsm}: {names}"
                )
            samples.append(
                Sample(
                    gse=gse_dir.name,
                    gsm=gsm,
                    fragment_file=matches[0],
                    individual_id=individual_map.get((gse_dir.name, gsm), ""),
                )
            )

    return samples


def selected_rna_gses(data_root: Path = DATA_ROOT, gse: str | None = None) -> list[str]:
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from scripts.process.download_from_datasets import load_filtered_rows

    rows = load_filtered_rows(
        data_root / "reference" / "datasets.xlsx",
        assay="scRNA",
        gse=gse or "",
        visible_rows_only=True,
    )
    gses = sorted({row.gse for row in rows})
    if gse and gse not in gses:
        raise FileNotFoundError(f"GSE not selected in datasets.xlsx for scRNA: {gse}")
    return gses


def discover_rna_samples_from_local_layout(
    raw_dir: Path, selected_gses: Iterable[str]
) -> list[RNASample]:
    samples: list[RNASample] = []

    for gse in sorted(selected_gses):
        gse_dir = raw_dir / gse
        if not gse_dir.exists():
            continue

        entries = sorted(path for path in gse_dir.iterdir() if path.is_file())
        supported_by_sample: dict[str, RNASample] = {}
        unsupported_by_sample: dict[str, RNASample] = {}

        triplet_parts: dict[str, dict[str, Path]] = {}
        for path in entries:
            name_upper = path.name.upper()
            if gse == "GSE226039" and "PBMC" not in name_upper:
                continue

            gsm_match = re.match(r"^(GSM\d+)", path.name)
            if gsm_match:
                sample_id = gsm_match.group(1)
                if path.suffix == ".h5":
                    supported_by_sample.setdefault(
                        sample_id,
                        RNASample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="h5",
                            supported=True,
                            note="10x h5",
                            h5_path=path,
                        ),
                    )
                    continue
                if path.name.endswith(".tar.gz") and "matrix" in path.name.lower():
                    supported_by_sample.setdefault(
                        sample_id,
                        RNASample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="archive",
                            supported=True,
                            note="matrix tar archive",
                            archive_path=path,
                        ),
                    )
                    continue
                if "matrix" in path.name.lower():
                    triplet_parts.setdefault(sample_id, {})["matrix"] = path
                elif "barcode" in path.name.lower():
                    triplet_parts.setdefault(sample_id, {})["barcodes"] = path
                elif "feature" in path.name.lower() or "gene" in path.name.lower():
                    triplet_parts.setdefault(sample_id, {})["features"] = path
                continue

            if path.name.startswith(f"{gse}_"):
                sample_id = gse
                if "matrix" in path.name.lower() and (
                    ".mtx" in path.name.lower() or path.name.lower().endswith(".tsv.gz")
                ):
                    triplet_parts.setdefault(sample_id, {})["matrix"] = path
                elif "barcode" in path.name.lower():
                    triplet_parts.setdefault(sample_id, {})["barcodes"] = path
                elif "feature" in path.name.lower() or "genes" in path.name.lower():
                    triplet_parts.setdefault(sample_id, {})["features"] = path
                elif (
                    "gene_count" in path.name.lower()
                    or "gene_counts" in path.name.lower()
                ):
                    unsupported_by_sample.setdefault(
                        sample_id,
                        RNASample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="shared-gene-count",
                            supported=False,
                            note="Unsupported shared gene-count matrix for single-cell RNA workflow",
                            matrix_path=path,
                        ),
                    )

        for sample_id, parts in sorted(triplet_parts.items()):
            if {"matrix", "barcodes", "features"} <= parts.keys():
                supported_by_sample.setdefault(
                    sample_id,
                    RNASample(
                        gse=gse,
                        sample_id=sample_id,
                        input_type="triplet",
                        supported=True,
                        note="matrix triplet",
                        matrix_path=parts["matrix"],
                        barcodes_path=parts["barcodes"],
                        features_path=parts["features"],
                    ),
                )

        merged = sorted(
            [*supported_by_sample.values(), *unsupported_by_sample.values()],
            key=lambda item: (item.gse, item.sample_id),
        )
        samples.extend(merged)

    return samples


def discover_rna_samples(
    gse: str | None = None, data_root: Path = DATA_ROOT
) -> list[RNASample]:
    selected_gses = selected_rna_gses(data_root=data_root, gse=gse)
    return discover_rna_samples_from_local_layout(
        raw_dir=data_root / "raw", selected_gses=selected_gses
    )


def find_rna_sample(gse: str, sample_id: str, data_root: Path = DATA_ROOT) -> RNASample:
    for sample in discover_rna_samples(gse=gse, data_root=data_root):
        if sample.sample_id == sample_id:
            return sample
    raise FileNotFoundError(f"RNA sample not found: {gse}/{sample_id}")


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


def rna_sample_output_dir(sample: RNASample, output_root: Path) -> Path:
    return output_root / sample.gse / sample.sample_id


def rna_sample_log_file(sample: RNASample, output_root: Path) -> Path:
    return rna_sample_output_dir(sample, output_root) / "logs" / "sample_qc.log"


def rna_sample_status_file(sample: RNASample, output_root: Path) -> Path:
    return rna_sample_output_dir(sample, output_root) / "run_status.json"


def expected_outputs(
    sample: Sample, output_profile: str = "full", output_root: Path = ATAC_OUTPUT_DIR
) -> list[Path]:
    output_dir = sample_output_dir(sample, output_root)
    if output_profile in {"matrix-lite", "validation-lite"}:
        return [
            output_dir / "umap_cima_cell_type_l1.png",
            output_dir / "umap_cima_cell_type_l2.png",
            output_dir / "qc_summary.csv",
            output_dir / "validation_result.csv",
            output_dir / f"{sample.gsm}.h5ad",
        ]

    return [
        output_dir / "qc_overview.png",
        output_dir / "umap_cima_cell_type_l1.png",
        output_dir / "umap_cima_cell_type_l2.png",
        output_dir / "metadata.csv",
        output_dir / "metadata_qc.csv",
        output_dir / "qc_summary.csv",
        output_dir / f"{sample.gsm}.h5ad",
    ]


def expected_rna_outputs(
    sample: RNASample, output_root: Path = RNA_OUTPUT_DIR
) -> list[Path]:
    discovered = only_rna_cli.DiscoveredSample(
        gse=sample.gse,
        sample_id=sample.sample_id,
        input_type=sample.input_type,
        sample_kind="gse_shared" if sample.sample_id == sample.gse else "gsm",
        supported=sample.supported,
        note=sample.note,
        source_name="",
        matrix_path=sample.matrix_path,
        barcodes_path=sample.barcodes_path,
        features_path=sample.features_path,
        h5_path=sample.h5_path,
        archive_path=sample.archive_path,
    )
    return only_rna_cli._expected_output_paths(discovered, output_root)


def outputs_complete(
    sample: Sample, output_profile: str = "full", output_root: Path = ATAC_OUTPUT_DIR
) -> bool:
    return all(
        path.exists()
        for path in expected_outputs(
            sample, output_profile=output_profile, output_root=output_root
        )
    )


def rna_outputs_complete(sample: RNASample, output_root: Path = RNA_OUTPUT_DIR) -> bool:
    return sample.supported and all(
        path.exists() for path in expected_rna_outputs(sample, output_root=output_root)
    )


def load_status(
    sample: Sample, output_root: Path = ATAC_OUTPUT_DIR
) -> dict[str, Any] | None:
    status_file = sample_status_file(sample, output_root)
    if not status_file.exists():
        return None
    return json.loads(status_file.read_text())


def write_status(
    sample: Sample, payload: dict[str, Any], output_root: Path = ATAC_OUTPUT_DIR
) -> None:
    output_dir = sample_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    sample_status_file(sample, output_root).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    )


def load_rna_status(
    sample: RNASample, output_root: Path = RNA_OUTPUT_DIR
) -> dict[str, Any] | None:
    status_file = rna_sample_status_file(sample, output_root)
    if not status_file.exists():
        return None
    return json.loads(status_file.read_text())


def write_rna_status(
    sample: RNASample, payload: dict[str, Any], output_root: Path = RNA_OUTPUT_DIR
) -> None:
    output_dir = rna_sample_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    rna_sample_status_file(sample, output_root).write_text(
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
    manifest_csv: str,
    sheet_name: str,
    file_kinds: str,
    aria2c: str,
    dry_run: bool,
    include_hidden_rows: bool,
    skip_network_resolve: bool,
    max_concurrent_downloads: int,
    split: int,
    network_timeout: int,
    links_out: str,
    print_links: bool,
) -> int:
    manifest_kind = "atac" if not manifest_csv and assay.lower() == "scatac" else "datasets"
    xlsx_path = ATAC_WORKBOOK if manifest_kind == "atac" else DATA_ROOT / "reference" / "datasets.xlsx"
    log_dir = RAW_DIR / "_download_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    assay_tag = re.sub(r"[^A-Za-z0-9._-]+", "-", assay).strip("-") or "all-assays"
    format_tag = (
        re.sub(r"[^A-Za-z0-9._-]+", "-", data_format).strip("-") or "all-formats"
    )
    log_file = log_dir / f"download_{assay_tag}_{format_tag}.log"

    command = [
        python_bin,
        str(DOWNLOAD_SCRIPT),
        "--xlsx-path",
        str(xlsx_path),
        "--manifest-kind",
        manifest_kind,
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

    if assay:
        command.extend(["--assay", assay])
    if data_format:
        command.extend(["--data-format", data_format])

    if gse:
        command.extend(["--gse", gse])
    if gsm:
        command.extend(["--gsm", gsm])
    if manifest_csv:
        command.extend(["--manifest-csv", manifest_csv])
    if sheet_name:
        command.extend(["--sheet-name", sheet_name])

    if links_out:
        command.extend(["--links-out", links_out])
    if print_links:
        command.append("--print-links")

    if dry_run:
        command.append("--dry-run")
    if include_hidden_rows:
        command.append("--include-hidden-rows")
    if skip_network_resolve:
        command.append("--skip-network-resolve")

    print("[download]")
    print(" ".join(command))
    return run_command(command, log_file)


def run_tea_seq_audit(
    python_bin: str,
    gse: str,
    dry_run: bool,
    delete_legacy: bool,
    score_threshold: float,
    margin_threshold: float,
    low_purity_threshold: float,
) -> int:
    audit_dir = OUTPUT_DIR / gse / "qc_audit"
    log_file = audit_dir / "tea_seq_audit.log"
    command = [
        python_bin,
        str(TEA_SEQ_AUDIT_SCRIPT),
        "--gse",
        gse,
        "--project-root",
        str(ROOT),
        "--score-threshold",
        str(score_threshold),
        "--margin-threshold",
        str(margin_threshold),
        "--low-purity-threshold",
        str(low_purity_threshold),
    ]
    if dry_run:
        command.append("--dry-run")
    if delete_legacy:
        command.append("--delete-legacy")

    print("[tea-seq-audit]")
    print(" ".join(command))
    return run_command(command, log_file)


def run_rna_sample(
    sample: RNASample,
    rscript: str,
    output_root: Path,
    force: bool,
    dry_run: bool,
) -> int:
    if not sample.supported:
        print(f"[skip] {sample.gse}/{sample.sample_id} unsupported: {sample.note}")
        return 2

    if rna_outputs_complete(sample, output_root=output_root) and not force:
        print(
            f"[skip] {sample.gse}/{sample.sample_id} already has complete RNA outputs"
        )
        return 0

    command = [
        rscript,
        str(RNA_R_SCRIPT),
        "--gse",
        sample.gse,
        "--sample-id",
        sample.sample_id,
        "--output-root",
        str(output_root),
        "--project-root",
        str(ROOT),
        *sample.command_args(),
    ]

    status_payload: dict[str, Any] = {
        "sample": {"gse": sample.gse, "sample_id": sample.sample_id},
        "command": command,
        "output_root": str(output_root),
        "input_type": sample.input_type,
        "started_at": utc_now(),
        "status": "running",
    }
    write_rna_status(sample, status_payload, output_root=output_root)

    print(f"[run-rna] {sample.gse}/{sample.sample_id}")
    print(" ".join(command))

    if dry_run:
        status_payload["status"] = "dry_run"
        status_payload["finished_at"] = utc_now()
        write_rna_status(sample, status_payload, output_root=output_root)
        return 0

    returncode = run_command(command, rna_sample_log_file(sample, output_root))
    status_payload["finished_at"] = utc_now()
    status_payload["returncode"] = returncode
    status_payload["outputs_complete"] = rna_outputs_complete(
        sample, output_root=output_root
    )
    status_payload["status"] = (
        "success"
        if returncode == 0 and rna_outputs_complete(sample, output_root=output_root)
        else "failed"
    )
    write_rna_status(sample, status_payload, output_root=output_root)
    return returncode


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
        "--individual-id",
        sample.individual_id,
        "--nmads",
        str(nmads),
        "--output-profile",
        output_profile,
        "--output-root",
        str(output_root),
    ]

    status_payload: dict[str, Any] = {
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

    log_file = sample_log_file(sample, output_root)
    returncode = run_command(command, log_file)
    if returncode == 0:
        export_command = [
            sys.executable,
            str(EXPORT_ATAC_H5AD_SCRIPT),
            "--sample-dir",
            str(sample_output_dir(sample, output_root)),
            "--data-root",
            str(DATA_ROOT),
            "--overwrite",
            "--cleanup",
        ]
        with log_file.open("a", encoding="utf-8") as handle:
            handle.write("\n[export-atac-h5ad]\n")
            handle.write(" ".join(export_command) + "\n")
        export_returncode = run_command(export_command, log_file)
        if export_returncode != 0:
            returncode = export_returncode
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


def print_rna_discovery(samples: Iterable[RNASample]) -> None:
    print("gse\tsample_id\tinput_type\tsupported\tnote\tprimary_path")
    for sample in samples:
        primary = (
            sample.matrix_path
            or sample.h5_path
            or sample.archive_path
            or sample.features_path
            or sample.barcodes_path
        )
        primary_text = str(primary) if primary is not None else ""
        print(
            f"{sample.gse}\t{sample.sample_id}\t{sample.input_type}\t"
            f"{str(sample.supported).lower()}\t{sample.note}\t{primary_text}"
        )


def print_status(samples: Iterable[Sample]) -> None:
    print("gse\tgsm\tstatus\toutputs_complete\tlast_finished")
    for sample in samples:
        payload = load_status(sample) or {}
        status = payload.get("status", "pending")
        output_profile = payload.get("output_profile", "full")
        output_root = Path(payload.get("output_root", str(ATAC_OUTPUT_DIR)))
        finished_at = payload.get("finished_at", "")
        print(
            f"{sample.gse}\t{sample.gsm}\t{status}\t"
            f"{str(outputs_complete(sample, output_profile=output_profile, output_root=output_root)).lower()}\t{finished_at}"
        )


def print_rna_status(samples: Iterable[RNASample]) -> None:
    print("gse\tsample_id\tstatus\toutputs_complete\tsupported\tlast_finished")
    for sample in samples:
        payload = load_rna_status(sample) or {}
        status = payload.get("status", "pending")
        output_root = Path(payload.get("output_root", str(RNA_OUTPUT_DIR)))
        finished_at = payload.get("finished_at", "")
        print(
            f"{sample.gse}\t{sample.sample_id}\t{status}\t"
            f"{str(rna_outputs_complete(sample, output_root=output_root)).lower()}\t"
            f"{str(sample.supported).lower()}\t{finished_at}"
        )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Manage scATAC single-sample QC runs")
    subparsers = parser.add_subparsers(dest="command", required=True)

    discover_parser = subparsers.add_parser("discover", help="List discovered samples")
    discover_parser.add_argument("--gse", help="Restrict discovery to one GSE")

    discover_rna_parser = subparsers.add_parser(
        "discover-rna", help="List discovered RNA samples"
    )
    discover_rna_parser.add_argument("--gse", help="Restrict RNA discovery to one GSE")

    run_sample_parser = subparsers.add_parser("run-sample", help="Run one sample")
    run_sample_parser.add_argument("--gse", required=True)
    run_sample_parser.add_argument("--gsm", required=True)
    run_sample_parser.add_argument("--nmads", type=float, default=4)
    run_sample_parser.add_argument("--rscript", default="Rscript")
    run_sample_parser.add_argument("--output-profile", default="full")
    run_sample_parser.add_argument("--output-root", default=str(ATAC_OUTPUT_DIR))
    run_sample_parser.add_argument("--force", action="store_true")
    run_sample_parser.add_argument("--dry-run", action="store_true")

    run_gse_parser = subparsers.add_parser(
        "run-gse", help="Run all samples under one GSE"
    )
    run_gse_parser.add_argument("--gse", required=True)
    run_gse_parser.add_argument("--nmads", type=float, default=4)
    run_gse_parser.add_argument("--rscript", default="Rscript")
    run_gse_parser.add_argument("--output-profile", default="full")
    run_gse_parser.add_argument("--output-root", default=str(ATAC_OUTPUT_DIR))
    run_gse_parser.add_argument("--force", action="store_true")
    run_gse_parser.add_argument("--dry-run", action="store_true")

    run_rna_sample_parser = subparsers.add_parser(
        "run-rna-sample", help="Run one RNA sample"
    )
    run_rna_sample_parser.add_argument("--gse", required=True)
    run_rna_sample_parser.add_argument("--sample-id", required=True)
    run_rna_sample_parser.add_argument(
        "--python-bin", default=sys.executable or "python3"
    )
    run_rna_sample_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "rna")
    )
    run_rna_sample_parser.add_argument("--force", action="store_true")
    run_rna_sample_parser.add_argument("--dry-run", action="store_true")

    run_rna_gse_parser = subparsers.add_parser(
        "run-rna-gse", help="Run all supported RNA samples under one GSE"
    )
    run_rna_gse_parser.add_argument("--gse", required=True)
    run_rna_gse_parser.add_argument("--python-bin", default=sys.executable or "python3")
    run_rna_gse_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "rna")
    )
    run_rna_gse_parser.add_argument("--force", action="store_true")
    run_rna_gse_parser.add_argument("--dry-run", action="store_true")

    tune_rna_sample_parser = subparsers.add_parser(
        "tune-rna-sample", help="Run baseline-only tuning for one RNA sample"
    )
    tune_rna_sample_parser.add_argument("--gse", required=True)
    tune_rna_sample_parser.add_argument("--sample-id", required=True)
    tune_rna_sample_parser.add_argument(
        "--python-bin", default=sys.executable or "python3"
    )
    tune_rna_sample_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "rna")
    )
    tune_rna_sample_parser.add_argument("--force", action="store_true")
    tune_rna_sample_parser.add_argument("--dry-run", action="store_true")

    tune_rna_gse_parser = subparsers.add_parser(
        "tune-rna-gse",
        help="Run baseline-only tuning for all supported RNA samples under one GSE",
    )
    tune_rna_gse_parser.add_argument("--gse", required=True)
    tune_rna_gse_parser.add_argument(
        "--python-bin", default=sys.executable or "python3"
    )
    tune_rna_gse_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "rna")
    )
    tune_rna_gse_parser.add_argument("--force", action="store_true")
    tune_rna_gse_parser.add_argument("--dry-run", action="store_true")

    download_parser = subparsers.add_parser(
        "download", help="Download GEO files based on datasets.xlsx filtering"
    )
    download_parser.add_argument("--manifest-csv", default="")
    download_parser.add_argument("--sheet-name", default="")
    download_parser.add_argument("--assay", default="")
    download_parser.add_argument("--data-format", default="")
    download_parser.add_argument("--gse", default="")
    download_parser.add_argument("--gsm", default="")
    download_parser.add_argument(
        "--file-kinds",
        default="all",
        help=(
            "Comma-separated kinds: all, count-matrix, fragment, barcode, singlecell, "
            "summary, h5, metadata"
        ),
    )
    download_parser.add_argument("--aria2c", default="aria2c")
    download_parser.add_argument("--python-bin", default=sys.executable or "python3")
    download_parser.add_argument("--dry-run", action="store_true")
    download_parser.add_argument("--include-hidden-rows", action="store_true")
    download_parser.add_argument("--skip-network-resolve", action="store_true")
    download_parser.add_argument("--max-concurrent-downloads", type=int, default=4)
    download_parser.add_argument("--split", type=int, default=8)
    download_parser.add_argument("--network-timeout", type=int, default=10)
    download_parser.add_argument("--links-out", default="")
    download_parser.add_argument("--print-links", action="store_true")

    status_parser = subparsers.add_parser("status", help="Show run status")
    status_parser.add_argument("--gse", help="Restrict status to one GSE")

    rna_status_parser = subparsers.add_parser("rna-status", help="Show RNA run status")
    rna_status_parser.add_argument("--gse", help="Restrict RNA status to one GSE")

    co_discover_parser = subparsers.add_parser(
        "co-discover", help="List co-assay RNA and ATAC samples"
    )
    co_discover_parser.add_argument(
        "--gse", help="Restrict co-assay discovery to one GSE"
    )

    co_status_parser = subparsers.add_parser(
        "co-status", help="Show co-ATAC only_atac run status"
    )
    co_status_parser.add_argument("--gse", help="Restrict co-ATAC status to one GSE")
    co_status_parser.add_argument("--output-root", default=str(ROOT / "output" / "co"))

    co_run_rna_sample_parser = subparsers.add_parser(
        "co-run-rna-sample", help="Run one co-assay RNA sample via only_rna"
    )
    co_run_rna_sample_parser.add_argument("--gse", required=True)
    co_run_rna_sample_parser.add_argument("--sample-id", required=True)
    co_run_rna_sample_parser.add_argument("--python-bin", default=sys.executable or "python3")
    co_run_rna_sample_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "co")
    )
    co_run_rna_sample_parser.add_argument("--force", action="store_true")
    co_run_rna_sample_parser.add_argument("--dry-run", action="store_true")

    co_run_rna_gse_parser = subparsers.add_parser(
        "co-run-rna-gse", help="Run all co-assay RNA samples via only_rna"
    )
    co_run_rna_gse_parser.add_argument("--gse", required=True)
    co_run_rna_gse_parser.add_argument("--python-bin", default=sys.executable or "python3")
    co_run_rna_gse_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "co")
    )
    co_run_rna_gse_parser.add_argument("--force", action="store_true")
    co_run_rna_gse_parser.add_argument("--dry-run", action="store_true")

    co_run_atac_sample_parser = subparsers.add_parser(
        "co-run-atac-sample", help="Run one co-ATAC sample via only_atac"
    )
    co_run_atac_sample_parser.add_argument("--gse", required=True)
    co_run_atac_sample_parser.add_argument("--gsm", required=True)
    co_run_atac_sample_parser.add_argument("--rscript", default="Rscript")
    co_run_atac_sample_parser.add_argument("--output-profile", default="full")
    co_run_atac_sample_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "co")
    )
    co_run_atac_sample_parser.add_argument("--nmads", type=float, default=4)
    co_run_atac_sample_parser.add_argument("--force", action="store_true")
    co_run_atac_sample_parser.add_argument("--dry-run", action="store_true")

    co_run_atac_gse_parser = subparsers.add_parser(
        "co-run-atac-gse", help="Run all co-ATAC samples under one GSE via only_atac"
    )
    co_run_atac_gse_parser.add_argument("--gse", required=True)
    co_run_atac_gse_parser.add_argument("--rscript", default="Rscript")
    co_run_atac_gse_parser.add_argument("--output-profile", default="full")
    co_run_atac_gse_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "co")
    )
    co_run_atac_gse_parser.add_argument("--nmads", type=float, default=4)
    co_run_atac_gse_parser.add_argument("--jobs", type=int, default=1)
    co_run_atac_gse_parser.add_argument("--force", action="store_true")
    co_run_atac_gse_parser.add_argument("--dry-run", action="store_true")

    longevity_discover_parser = subparsers.add_parser(
        "longevity-discover", help="List independent longevity RNA atlas and ATAC samples"
    )
    longevity_discover_parser.add_argument("--sample-id", help="Restrict ATAC discovery to one sample")

    longevity_status_parser = subparsers.add_parser(
        "longevity-status", help="Show independent longevity ATAC run status"
    )
    longevity_status_parser.add_argument("--sample-id", help="Restrict ATAC status to one sample")
    longevity_status_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity")
    )

    longevity_run_atac_sample_parser = subparsers.add_parser(
        "longevity-run-atac-sample", help="Run one independent longevity ATAC sample"
    )
    longevity_run_atac_sample_parser.add_argument("--sample-id", required=True)
    longevity_run_atac_sample_parser.add_argument("--rscript", default="Rscript")
    longevity_run_atac_sample_parser.add_argument("--output-profile", default="full")
    longevity_run_atac_sample_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity")
    )
    longevity_run_atac_sample_parser.add_argument("--nmads", type=float, default=4)
    longevity_run_atac_sample_parser.add_argument("--min-inferred-fragments", type=float, default=None)
    longevity_run_atac_sample_parser.add_argument("--max-inferred-barcodes", type=int, default=None)
    longevity_run_atac_sample_parser.add_argument("--umap-min-dist", type=float, default=None)
    longevity_run_atac_sample_parser.add_argument(
        "--barcode-output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes"),
    )
    longevity_run_atac_sample_parser.add_argument("--barcode-min-fragments", type=int, default=200)
    longevity_run_atac_sample_parser.add_argument("--barcode-max-barcodes", type=int, default=20000)
    longevity_run_atac_sample_parser.add_argument("--barcode-min-tss", type=float, default=2.5)
    longevity_run_atac_sample_parser.add_argument(
        "--barcode-rank-by", choices=["fragments", "tss_then_fragments"], default="fragments"
    )
    longevity_run_atac_sample_parser.add_argument("--barcode-overrides", type=Path, default=None)
    longevity_run_atac_sample_parser.add_argument("--force", action="store_true")
    longevity_run_atac_sample_parser.add_argument("--dry-run", action="store_true")

    longevity_run_atac_all_parser = subparsers.add_parser(
        "longevity-run-atac-all", help="Run all independent longevity ATAC samples"
    )
    longevity_run_atac_all_parser.add_argument("--rscript", default="Rscript")
    longevity_run_atac_all_parser.add_argument("--output-profile", default="full")
    longevity_run_atac_all_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity")
    )
    longevity_run_atac_all_parser.add_argument("--nmads", type=float, default=4)
    longevity_run_atac_all_parser.add_argument("--min-inferred-fragments", type=float, default=None)
    longevity_run_atac_all_parser.add_argument("--max-inferred-barcodes", type=int, default=None)
    longevity_run_atac_all_parser.add_argument("--umap-min-dist", type=float, default=None)
    longevity_run_atac_all_parser.add_argument(
        "--barcode-output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes"),
    )
    longevity_run_atac_all_parser.add_argument("--barcode-min-fragments", type=int, default=200)
    longevity_run_atac_all_parser.add_argument("--barcode-max-barcodes", type=int, default=20000)
    longevity_run_atac_all_parser.add_argument("--barcode-min-tss", type=float, default=2.5)
    longevity_run_atac_all_parser.add_argument(
        "--barcode-rank-by", choices=["fragments", "tss_then_fragments"], default="fragments"
    )
    longevity_run_atac_all_parser.add_argument("--barcode-overrides", type=Path, default=None)
    longevity_run_atac_all_parser.add_argument("--force", action="store_true")
    longevity_run_atac_all_parser.add_argument("--dry-run", action="store_true")

    longevity_ingest_rna_parser = subparsers.add_parser(
        "longevity-ingest-rna",
        help="Write longevity RNA h5ad with 5-class final_celltype and comparison table",
    )
    longevity_ingest_rna_parser.add_argument("--sample-id", default=None)
    longevity_ingest_rna_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity")
    )
    longevity_ingest_rna_parser.add_argument("--dry-run", action="store_true")

    longevity_preprocess_barcodes_parser = subparsers.add_parser(
        "longevity-preprocess-atac-barcodes",
        help="Generate ArchR-based filtered barcode files for longevity ATAC samples",
    )
    longevity_preprocess_barcodes_parser.add_argument("--sample-id", default=None)
    longevity_preprocess_barcodes_parser.add_argument("--rscript", default="Rscript")
    longevity_preprocess_barcodes_parser.add_argument(
        "--output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes"),
    )
    longevity_preprocess_barcodes_parser.add_argument("--min-fragments", type=int, default=200)
    longevity_preprocess_barcodes_parser.add_argument("--max-barcodes", type=int, default=20000)
    longevity_preprocess_barcodes_parser.add_argument("--min-tss", type=float, default=2.5)
    longevity_preprocess_barcodes_parser.add_argument(
        "--rank-by", choices=["fragments", "tss_then_fragments"], default="fragments"
    )
    longevity_preprocess_barcodes_parser.add_argument("--overrides", type=Path, default=None)
    longevity_preprocess_barcodes_parser.add_argument("--force", action="store_true")
    longevity_preprocess_barcodes_parser.add_argument("--dry-run", action="store_true")
    longevity_preprocess_barcodes_parser.add_argument("--keep-going", action="store_true")

    longevity_barcode_status_parser = subparsers.add_parser(
        "longevity-atac-barcode-status",
        help="Show longevity ATAC preprocessed barcode availability",
    )
    longevity_barcode_status_parser.add_argument("--sample-id", default=None)
    longevity_barcode_status_parser.add_argument(
        "--output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes"),
    )
    longevity_barcode_status_parser.add_argument(
        "--format", choices=["table", "csv", "json"], default="table"
    )

    longevity_param_summary_parser = subparsers.add_parser(
        "longevity-summarize-atac-param-contrast",
        help="Summarize longevity ATAC parameter contrast outputs and select best candidate per sample",
    )
    longevity_param_summary_parser.add_argument(
        "--contrast-root", default=str(ROOT / "output" / "longevity_param_contrast")
    )
    longevity_param_summary_parser.add_argument("--min-pass-qc", type=int, default=5000)
    longevity_param_summary_parser.add_argument("--min-frip", type=float, default=0.6)
    longevity_param_summary_parser.add_argument(
        "--format", choices=["table", "csv", "json"], default="table"
    )
    longevity_param_summary_parser.add_argument("--write", action="store_true")

    longevity_param_run_parser = subparsers.add_parser(
        "longevity-run-atac-param-contrast",
        help="Run longevity ATAC samples across predefined barcode parameter candidates",
    )
    longevity_param_run_parser.add_argument("--sample-id", action="append", default=[])
    longevity_param_run_parser.add_argument("--candidate", action="append", default=[])
    longevity_param_run_parser.add_argument("--rscript", default="Rscript")
    longevity_param_run_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity_param_contrast")
    )
    longevity_param_run_parser.add_argument(
        "--barcode-output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes_param_contrast"),
    )
    longevity_param_run_parser.add_argument("--nmads", type=float, default=4)
    longevity_param_run_parser.add_argument("--umap-min-dist", type=float, default=None)
    longevity_param_run_parser.add_argument("--skip-complete", action="store_true")
    longevity_param_run_parser.add_argument("--keep-going", action="store_true")
    longevity_param_run_parser.add_argument("--dry-run", action="store_true")

    longevity_custom_param_run_parser = subparsers.add_parser(
        "longevity-run-atac-custom-param-contrast",
        help="Run one longevity ATAC sample with one custom barcode parameter candidate",
    )
    longevity_custom_param_run_parser.add_argument("--sample-id", required=True)
    longevity_custom_param_run_parser.add_argument("--candidate-name", default=None)
    longevity_custom_param_run_parser.add_argument("--min-fragments", type=int, required=True)
    longevity_custom_param_run_parser.add_argument("--max-barcodes", type=int, required=True)
    longevity_custom_param_run_parser.add_argument("--min-tss", type=float, required=True)
    longevity_custom_param_run_parser.add_argument(
        "--rank-by", choices=["fragments", "tss_then_fragments"], default="fragments"
    )
    longevity_custom_param_run_parser.add_argument("--rscript", default="Rscript")
    longevity_custom_param_run_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity_param_contrast")
    )
    longevity_custom_param_run_parser.add_argument(
        "--barcode-output-root",
        default=str(DATA_ROOT / "reference" / "longevity" / "atac_barcodes_param_contrast"),
    )
    longevity_custom_param_run_parser.add_argument("--nmads", type=float, default=4)
    longevity_custom_param_run_parser.add_argument("--umap-min-dist", type=float, default=None)
    longevity_custom_param_run_parser.add_argument("--skip-complete", action="store_true")
    longevity_custom_param_run_parser.add_argument("--dry-run", action="store_true")

    longevity_param_publish_parser = subparsers.add_parser(
        "longevity-publish-atac-param-contrast",
        help="Copy selected longevity ATAC parameter contrast outputs into the formal longevity output tree",
    )
    longevity_param_publish_parser.add_argument(
        "--selection-csv",
        default=str(ROOT / "output" / "longevity_param_contrast" / "param_contrast_selected_samples.csv"),
    )
    longevity_param_publish_parser.add_argument(
        "--output-root", default=str(ROOT / "output" / "longevity")
    )
    longevity_param_publish_parser.add_argument("--force", action="store_true")
    longevity_param_publish_parser.add_argument("--dry-run", action="store_true")

    backfill_rna_final_parser = subparsers.add_parser(
        "backfill-rna-final-celltype",
        help="Backfill RNA outputs to the 5-class final_celltype vocabulary",
    )
    backfill_rna_final_parser.add_argument(
        "--root", default=str(ROOT / "output" / "co" / "rna" / "GSE224198")
    )
    backfill_rna_final_parser.add_argument("--dry-run", action="store_true")
    backfill_rna_final_parser.add_argument("--fail-on-skip", action="store_true")

    tea_seq_parser = subparsers.add_parser(
        "tea-seq-audit",
        help="Organize TEA-seq accepted outputs and write a dataset-level QC audit",
    )
    tea_seq_parser.add_argument("--gse", required=True)
    tea_seq_parser.add_argument("--python-bin", default=sys.executable or "python3")
    tea_seq_parser.add_argument("--dry-run", action="store_true")
    tea_seq_parser.add_argument("--delete-legacy", action="store_true")
    tea_seq_parser.add_argument("--score-threshold", type=float, default=0.4)
    tea_seq_parser.add_argument("--margin-threshold", type=float, default=0.1)
    tea_seq_parser.add_argument("--low-purity-threshold", type=float, default=0.8)

    organize_products_parser = subparsers.add_parser(
        "organize-products",
        help="Organize annotated branch outputs into integrated product directories",
    )
    organize_products_parser.add_argument(
        "--products",
        default="all",
        help="Comma-separated products: only_atac,only_rna,co_atac,co_rna,all",
    )
    organize_products_parser.add_argument("--output-root", default=str(OUTPUT_DIR))
    organize_products_parser.add_argument(
        "--copy-mode", choices=["symlink", "copy"], default="symlink"
    )
    organize_products_parser.add_argument("--force", action="store_true")
    organize_products_parser.add_argument("--skip-figures", action="store_true")
    organize_products_parser.add_argument("--skip-integration", action="store_true")
    organize_products_parser.add_argument("--integration-n-components", type=int, default=30)
    organize_products_parser.add_argument(
        "--integration-max-umap-fit-cells", type=int, default=100_000
    )
    organize_products_parser.add_argument("--integration-clusters", type=int, default=30)
    organize_products_parser.add_argument("--integration-batch-key", default=None)
    organize_products_parser.add_argument(
        "--integration-method",
        choices=["bbknn", "scanpy_neighbors", "harmony"],
        default=None,
    )
    organize_products_parser.add_argument("--bbknn-neighbors-within-batch", type=int, default=None)
    organize_products_parser.add_argument("--bbknn-trim", type=int, default=None)
    organize_products_parser.add_argument("--leiden-resolution", type=float, default=1.0)
    organize_products_parser.add_argument("--rna-min-cima-l1-score", type=float, default=None)
    organize_products_parser.add_argument("--include-incomplete", action="store_true")

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "discover":
        print_discovery(discover_samples(gse=args.gse))
        return 0

    if args.command == "discover-rna":
        return only_rna_cli.cmd_discover_rna(args)

    if args.command == "status":
        print_status(discover_samples(gse=args.gse))
        return 0

    if args.command == "rna-status":
        return only_rna_cli.cmd_rna_status(args)

    if args.command == "co-discover":
        return co_cli.cmd_discover(args)

    if args.command == "co-status":
        return co_cli.cmd_status(args)

    if args.command == "co-run-rna-sample":
        return co_cli.cmd_run_rna_sample(args)

    if args.command == "co-run-rna-gse":
        return co_cli.cmd_run_rna_gse(args)

    if args.command == "co-run-atac-sample":
        return co_cli.cmd_run_atac_sample(args)

    if args.command == "co-run-atac-gse":
        return co_cli.cmd_run_atac_gse(args)

    if args.command == "longevity-discover":
        return longevity_cli.cmd_discover(args)

    if args.command == "longevity-status":
        return longevity_cli.cmd_status(args)

    if args.command == "longevity-run-atac-sample":
        return longevity_cli.cmd_run_atac_sample(args)

    if args.command == "longevity-run-atac-all":
        return longevity_cli.cmd_run_atac_all(args)

    if args.command == "longevity-ingest-rna":
        return longevity_cli.cmd_ingest_rna(args)

    if args.command == "longevity-preprocess-atac-barcodes":
        return longevity_cli.cmd_preprocess_atac_barcodes(args)

    if args.command == "longevity-atac-barcode-status":
        return longevity_cli.cmd_atac_barcode_status(args)

    if args.command == "longevity-summarize-atac-param-contrast":
        return longevity_cli.cmd_summarize_atac_param_contrast(args)

    if args.command == "longevity-run-atac-param-contrast":
        return longevity_cli.cmd_run_atac_param_contrast(args)

    if args.command == "longevity-run-atac-custom-param-contrast":
        return longevity_cli.cmd_run_atac_custom_param_contrast(args)

    if args.command == "longevity-publish-atac-param-contrast":
        return longevity_cli.cmd_publish_atac_param_contrast(args)

    if args.command == "backfill-rna-final-celltype":
        return backfill_final_celltype.cmd_backfill_rna_final_celltype(args)

    if args.command == "tea-seq-audit":
        return run_tea_seq_audit(
            python_bin=args.python_bin,
            gse=args.gse,
            dry_run=args.dry_run,
            delete_legacy=args.delete_legacy,
            score_threshold=args.score_threshold,
            margin_threshold=args.margin_threshold,
            low_purity_threshold=args.low_purity_threshold,
        )

    if args.command == "organize-products":
        return organize_integrated_products.cmd_organize_products(args)

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

    if args.command == "run-rna-sample":
        return only_rna_cli.cmd_run_rna_sample(args)

    if args.command == "run-rna-gse":
        return only_rna_cli.cmd_run_rna_gse(args)

    if args.command == "tune-rna-sample":
        return only_rna_cli.cmd_tune_rna_sample(args)

    if args.command == "tune-rna-gse":
        return only_rna_cli.cmd_tune_rna_gse(args)

    if args.command == "download":
        return run_download(
            python_bin=args.python_bin,
            assay=args.assay,
            data_format=args.data_format,
            gse=args.gse,
            gsm=args.gsm,
            manifest_csv=args.manifest_csv,
            sheet_name=args.sheet_name,
            file_kinds=args.file_kinds,
            aria2c=args.aria2c,
            dry_run=args.dry_run,
            include_hidden_rows=args.include_hidden_rows,
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
