from __future__ import annotations

import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from scripts.only_rna import cli as only_rna_cli
from scripts.only_rna.discovery import resolve_data_root


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = ROOT / "output" / "co"
ONLY_ATAC_R_SCRIPT = ROOT / "scripts" / "process" / "process_single_sample.R"
EXPORT_CO_ATAC_H5AD = ROOT / "scripts" / "process" / "export_co_atac_h5ad.py"


@dataclass(frozen=True)
class CoSample:
    gse: str
    gsm: str
    assay: str
    individual_id: str
    age: str = ""
    health_status: str = ""
    is_pbmc: str = ""
    max_barcodes_hint: int | None = None
    primary_path: Path | None = None
    barcode_path: Path | None = None
    filtered_metadata_path: Path | None = None
    singlecell_path: Path | None = None
    supported: bool = True
    note: str = ""


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _co_workbook_path() -> Path:
    return resolve_data_root(ROOT) / "reference" / "co.xlsx"


def _raw_root() -> Path:
    return resolve_data_root(ROOT) / "raw"


def _load_co_layout(gse: str | None = None) -> pd.DataFrame:
    workbook = _co_workbook_path()
    if not workbook.exists():
        return pd.DataFrame()
    df = pd.read_excel(workbook, dtype=str).fillna("")
    required_columns = {"sample", "dataset", "age", "health", "donor"}
    missing_columns = required_columns - set(df.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"co.xlsx missing required columns: {missing}")
    for column in required_columns:
        df[column] = df[column].astype(str).str.strip()
    df = df.loc[(df["sample"] != "") & (df["dataset"] != "")].copy()
    if df.empty:
        return pd.DataFrame()
    if gse is not None:
        df = df.loc[df["dataset"] == gse].copy()
    return df


def _discover_co_atac(gse: str | None = None) -> list[CoSample]:
    raw_root = _raw_root()
    samples: list[CoSample] = []
    co_layout = _load_co_layout(gse)
    for _, row in co_layout.iterrows():
        dataset = str(row.get("dataset", "")).strip()
        sample_id = str(row.get("sample", "")).strip()
        atac_dir = raw_root / dataset / sample_id / "ATAC"
        fragment_file = atac_dir / "fragments.tsv.gz"
        barcode_file = atac_dir / "filtered_barcodes.tsv.gz"
        filtered_metadata_file = atac_dir / "filtered_metadata.csv.gz"
        singlecell_file = atac_dir / "singlecell.csv.gz"
        supported = bool(sample_id and fragment_file.exists())
        note = "fragment file" if supported else f"fragment file not found: {fragment_file}"
        samples.append(
            CoSample(
                gse=dataset,
                gsm=sample_id,
                assay="ATAC（RNA+ATAC）",
                individual_id=str(row.get("donor", "")).strip(),
                age=str(row.get("age", "")).strip(),
                health_status=str(row.get("health", "")).strip(),
                is_pbmc="是",
                primary_path=fragment_file if supported else None,
                barcode_path=barcode_file if barcode_file.exists() else None,
                filtered_metadata_path=filtered_metadata_file if filtered_metadata_file.exists() else None,
                singlecell_path=singlecell_file if singlecell_file.exists() else None,
                supported=supported,
                note=note,
            )
        )
    return samples


def _discover_co_rna(gse: str | None = None) -> list[Any]:
    co_layout = _load_co_layout(gse)
    if not co_layout.empty:
        from scripts.only_rna.discovery import DiscoveredSample

        raw_root = _raw_root()
        samples = []
        for _, row in co_layout.iterrows():
            dataset = str(row.get("dataset", "")).strip()
            sample_id = str(row.get("sample", "")).strip()
            rna_dir = raw_root / dataset / sample_id / "RNA"
            matrix = next(
                (path for path in [rna_dir / "matrix.mtx", rna_dir / "matrix.mtx.gz"] if path.exists()),
                rna_dir / "matrix.mtx",
            )
            barcodes = next(
                (path for path in [rna_dir / "barcodes.tsv", rna_dir / "barcodes.tsv.gz"] if path.exists()),
                rna_dir / "barcodes.tsv",
            )
            features = next(
                (path for path in [rna_dir / "features.tsv", rna_dir / "features.tsv.gz"] if path.exists()),
                rna_dir / "features.tsv",
            )
            supported = matrix.exists() and barcodes.exists() and features.exists()
            samples.append(
                DiscoveredSample(
                    gse=dataset,
                    sample_id=sample_id,
                    input_type="triplet",
                    sample_kind="gsm",
                    supported=supported,
                    note="matrix triplet" if supported else f"RNA triplet not found: {rna_dir}",
                    source_name=matrix.name,
                    individual_id=str(row.get("donor", "")).strip(),
                    age=str(row.get("age", "")).strip(),
                    health=str(row.get("health", "")).strip(),
                    donor=str(row.get("donor", "")).strip(),
                    matrix_path=matrix if supported else None,
                    barcodes_path=barcodes if supported else None,
                    features_path=features if supported else None,
                )
            )
        return samples

    return []


def _sample_output_dir(sample: CoSample, output_root: Path) -> Path:
    return output_root / "atac" / sample.gse / sample.gsm


def _status_file(sample: CoSample, output_root: Path) -> Path:
    return _sample_output_dir(sample, output_root) / "run_status.json"


def _expected_atac_outputs(sample: CoSample, output_root: Path, output_profile: str = "full") -> list[Path]:
    out = _sample_output_dir(sample, output_root)
    if output_profile in {"matrix-lite", "validation-lite"}:
        return [
            out / "umap_cima_cell_type_l1.png",
            out / "umap_cima_cell_type_l2.png",
            out / "qc_summary.csv",
            out / "validation_result.csv",
            out / f"{sample.gsm}.h5ad",
        ]

    return [
        out / "metadata.csv",
        out / "metadata_qc.csv",
        out / "qc_summary.csv",
        out / "qc_overview.png",
        out / "umap_cima_cell_type_l1.png",
        out / "umap_cima_cell_type_l2.png",
        out / f"{sample.gsm}.h5ad",
    ]


def _outputs_complete(sample: CoSample, output_root: Path, output_profile: str = "full") -> bool:
    return sample.supported and all(
        path.exists() for path in _expected_atac_outputs(sample, output_root, output_profile)
    )


def _load_status(sample: CoSample, output_root: Path) -> dict[str, Any] | None:
    path = _status_file(sample, output_root)
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_status(sample: CoSample, output_root: Path, payload: dict[str, Any]) -> None:
    output_dir = _sample_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    _status_file(sample, output_root).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _find_co_atac_sample(gse: str, gsm: str) -> CoSample:
    for sample in _discover_co_atac(gse):
        if sample.gsm == gsm:
            return sample
    raise FileNotFoundError(f"co ATAC sample not found: {gse}/{gsm}")


def _find_co_rna_sample(gse: str, sample_id: str):
    for sample in _discover_co_rna(gse):
        if sample.sample_id == sample_id:
            return sample
    raise FileNotFoundError(f"co RNA sample not found: {gse}/{sample_id}")


def _print_discovery(rna_samples: Iterable[Any], atac_samples: Iterable[CoSample]) -> None:
    print("modality\tgse\tsample_id\tindividual_id\tsupported\tnote\tprimary_path")
    for sample in rna_samples:
        primary = sample.archive_path or sample.h5_path or sample.matrix_path or sample.features_path
        individual_id = sample.individual_id
        print(
            f"RNA\t{sample.gse}\t{sample.sample_id}\t"
            f"{individual_id}\t"
            f"{str(sample.supported).lower()}\t{sample.note}\t{primary or ''}"
        )
    for sample in atac_samples:
        print(
            f"ATAC\t{sample.gse}\t{sample.gsm}\t{sample.individual_id}\t"
            f"{str(sample.supported).lower()}\t{sample.note}\t{sample.primary_path or ''}"
        )


def _print_status(atac_samples: Iterable[CoSample], output_root: Path) -> None:
    print("modality\tgse\tsample_id\tindividual_id\tstatus\toutputs_complete\tlast_finished")
    for sample in atac_samples:
        payload = _load_status(sample, output_root) or {}
        status = payload.get("status", "pending")
        output_profile = payload.get("output_profile", "full")
        finished_at = payload.get("finished_at", "")
        print(
            f"ATAC\t{sample.gse}\t{sample.gsm}\t{sample.individual_id}\t{status}\t"
            f"{str(_outputs_complete(sample, output_root, output_profile)).lower()}\t{finished_at}"
        )


def cmd_discover(args) -> int:
    _print_discovery(_discover_co_rna(getattr(args, "gse", None)), _discover_co_atac(getattr(args, "gse", None)))
    return 0


def cmd_status(args) -> int:
    _print_status(_discover_co_atac(getattr(args, "gse", None)), Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT)))
    return 0


def cmd_run_rna_sample(args) -> int:
    args.output_root = str(Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT)) / "rna")
    if not _load_co_layout(args.gse).empty:
        return only_rna_cli._route_rna_sample(_find_co_rna_sample(args.gse, args.sample_id), args)
    return only_rna_cli.cmd_run_rna_sample(args)


def cmd_run_rna_gse(args) -> int:
    args.output_root = str(Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT)) / "rna")
    if not _load_co_layout(args.gse).empty:
        returncode = 0
        for sample in _discover_co_rna(args.gse):
            current = only_rna_cli._route_rna_sample(sample, args)
            if current != 0:
                returncode = current
        return returncode
    return only_rna_cli.cmd_run_rna_gse(args)


def _build_atac_command(sample: CoSample, args, output_root: Path) -> list[str]:
    command = [
        getattr(args, "rscript", "Rscript"),
        str(ONLY_ATAC_R_SCRIPT),
        "--gse", sample.gse,
        "--gsm", sample.gsm,
        "--individual-id", sample.individual_id,
        "--nmads", str(getattr(args, "nmads", 4)),
        "--output-profile", str(getattr(args, "output_profile", "full")),
        "--output-root", str(output_root / "atac"),
    ]
    if sample.primary_path is not None:
        command.extend(["--fragment-file", str(sample.primary_path)])
    if sample.barcode_path is not None:
        command.extend(["--barcode-file", str(sample.barcode_path)])
    if sample.filtered_metadata_path is not None:
        command.extend(["--filtered-metadata-file", str(sample.filtered_metadata_path)])
    if sample.singlecell_path is not None:
        command.extend(["--singlecell-file", str(sample.singlecell_path)])
    if sample.max_barcodes_hint:
        command.extend(["--min-inferred-fragments", "1000"])
        command.extend(["--max-inferred-barcodes", str(sample.max_barcodes_hint)])
    return command


def _run_atac_sample(sample: CoSample, args) -> int:
    output_root = Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT))
    output_profile = str(getattr(args, "output_profile", "full"))
    force = bool(getattr(args, "force", False))
    dry_run = bool(getattr(args, "dry_run", False))
    if not sample.supported:
        print(f"[skip-co-atac] {sample.gse}/{sample.gsm} unsupported: {sample.note}", file=sys.stderr)
        return 2
    if _outputs_complete(sample, output_root, output_profile) and not force:
        print(f"[skip-co-atac] {sample.gse}/{sample.gsm} already complete")
        return 0

    command = _build_atac_command(sample, args, output_root)
    payload = {
        "sample": {"gse": sample.gse, "gsm": sample.gsm, "individual_id": sample.individual_id},
        "command": command,
        "output_profile": output_profile,
        "output_root": str(output_root),
        "started_at": utc_now(),
        "status": "dry_run" if dry_run else "running",
        "outputs_complete": _outputs_complete(sample, output_root, output_profile),
        "note": "Co-ATAC routed through only_atac single-sample workflow.",
    }
    print(f"[run-co-atac] {sample.gse}/{sample.gsm}")
    print(" ".join(command))
    if dry_run:
        payload["finished_at"] = payload["started_at"]
        _write_status(sample, output_root, payload)
        return 0

    _write_status(sample, output_root, payload)
    output_dir = _sample_output_dir(sample, output_root)
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    with (log_dir / "sample_qc.log").open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(command, cwd=ROOT, stdout=log_handle, stderr=subprocess.STDOUT, check=False)
        if result.returncode == 0:
            export_command = [
                sys.executable,
                str(EXPORT_CO_ATAC_H5AD),
                "--sample-dir",
                str(output_dir),
                "--data-root",
                str(resolve_data_root(ROOT)),
                "--overwrite",
                "--cleanup",
            ]
            log_handle.write("\n[export-co-atac-h5ad]\n")
            log_handle.write(" ".join(export_command) + "\n")
            export_result = subprocess.run(
                export_command,
                cwd=ROOT,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                check=False,
            )
            if export_result.returncode != 0:
                result = export_result

    payload["finished_at"] = utc_now()
    payload["returncode"] = str(result.returncode)
    payload["outputs_complete"] = _outputs_complete(sample, output_root, output_profile)
    payload["status"] = "success" if result.returncode == 0 and payload["outputs_complete"] else "failed"
    _write_status(sample, output_root, payload)
    return result.returncode


def cmd_run_atac_sample(args) -> int:
    return _run_atac_sample(_find_co_atac_sample(args.gse, args.gsm), args)


def cmd_run_atac_gse(args) -> int:
    samples = _discover_co_atac(args.gse)
    jobs = max(1, int(getattr(args, "jobs", 1)))
    if jobs == 1:
        returncode = 0
        for sample in samples:
            current = _run_atac_sample(sample, args)
            if current != 0:
                returncode = current
        return returncode

    returncode = 0
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {executor.submit(_run_atac_sample, sample, args): sample for sample in samples}
        for future in as_completed(futures):
            sample = futures[future]
            try:
                current = future.result()
            except Exception as exc:
                print(f"[failed-co-atac] {sample.gse}/{sample.gsm}: {exc}", file=sys.stderr)
                current = 1
            if current != 0:
                returncode = current
    return returncode
