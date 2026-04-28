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
from scripts.only_rna.discovery import discover_rna_samples, resolve_data_root


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT_ROOT = ROOT / "output" / "co"
ONLY_ATAC_R_SCRIPT = ROOT / "scripts" / "process" / "process_single_sample.R"


@dataclass(frozen=True)
class CoSample:
    gse: str
    gsm: str
    assay: str
    individual_id: str
    health_status: str = ""
    is_pbmc: str = ""
    folder_name: str = ""
    max_barcodes_hint: int | None = None
    primary_path: Path | None = None
    supported: bool = True
    note: str = ""


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _manifest_path() -> Path:
    return resolve_data_root(ROOT) / "reference" / "co2_sample_manifest.csv"


def _co1_layout_path() -> Path:
    return resolve_data_root(ROOT) / "raw" / "7555405" / "sample_layout.tsv"


def _raw_root() -> Path:
    return resolve_data_root(ROOT) / "raw"


def _load_manifest(gse: str | None = None) -> pd.DataFrame:
    manifest = _manifest_path()
    if not manifest.exists():
        return pd.DataFrame(columns=["gse", "gsm", "assay", "individual_id", "health_status", "is_pbmc"])
    df = pd.read_csv(manifest).fillna("")
    if gse is not None:
        df = df.loc[df["gse"].astype(str) == gse].copy()
    return df


def _load_co1_layout(gse: str | None = None) -> pd.DataFrame:
    if gse not in {None, "7555405"}:
        return pd.DataFrame()
    layout = _co1_layout_path()
    if not layout.exists():
        return pd.DataFrame()
    return pd.read_csv(layout, sep="\t").fillna("")


def _discover_co_atac(gse: str | None = None) -> list[CoSample]:
    raw_root = _raw_root()
    samples: list[CoSample] = []
    df = _load_manifest(gse) if gse != "7555405" else pd.DataFrame()
    for _, row in df.iterrows():
        assay = str(row.get("assay", "")).strip()
        if not assay.startswith("ATAC"):
            continue
        gse_value = str(row.get("gse", "")).strip()
        gsm = str(row.get("gsm", "")).strip()
        matches = sorted((raw_root / gse_value).glob(f"{gsm}_*fragments*.tsv.gz"))
        supported = len(matches) == 1
        note = "fragment file" if supported else f"expected 1 fragment file, found {len(matches)}"
        samples.append(
            CoSample(
                gse=gse_value,
                gsm=gsm,
                assay=assay,
                individual_id=str(row.get("individual_id", "")).strip(),
                health_status=str(row.get("health_status", "")).strip(),
                is_pbmc=str(row.get("is_pbmc", "")).strip(),
                primary_path=matches[0] if supported else None,
                supported=supported,
                note=note,
            )
        )
    co1 = _load_co1_layout(gse)
    co1_root = raw_root / "7555405"
    for _, row in co1.iterrows():
        folder_name = str(row.get("folder_name", "")).strip()
        sample_id = str(row.get("sample_id", "")).strip()
        fragment_file = co1_root / folder_name / "ATAC" / "fragments.tsv.gz"
        supported = bool(sample_id and fragment_file.exists())
        note = "fragment file" if supported else f"fragment file not found: {fragment_file}"
        try:
            max_barcodes_hint = int(str(row.get("atac_ncells", 0)).strip()) or None
        except (TypeError, ValueError):
            max_barcodes_hint = None
        samples.append(
            CoSample(
                gse="7555405",
                gsm=sample_id,
                assay="ATAC（RNA+ATAC）",
                individual_id=sample_id,
                health_status="健康",
                is_pbmc="是",
                folder_name=folder_name,
                max_barcodes_hint=max_barcodes_hint,
                primary_path=fragment_file if supported else None,
                supported=supported,
                note=note,
            )
        )
    return samples


def _discover_co_rna(gse: str | None = None) -> list[Any]:
    co1 = _load_co1_layout(gse)
    if not co1.empty:
        from scripts.only_rna.discovery import DiscoveredSample

        raw_root = _raw_root() / "7555405"
        samples = []
        for _, row in co1.iterrows():
            folder_name = str(row.get("folder_name", "")).strip()
            sample_id = str(row.get("sample_id", "")).strip()
            rna_dir = raw_root / folder_name / "RNA"
            matrix = rna_dir / "matrix.mtx"
            barcodes = rna_dir / "barcodes.tsv"
            features = rna_dir / "features.tsv"
            supported = matrix.exists() and barcodes.exists() and features.exists()
            samples.append(
                DiscoveredSample(
                    gse="7555405",
                    sample_id=sample_id,
                    input_type="triplet",
                    sample_kind="gsm",
                    supported=supported,
                    note="matrix triplet" if supported else f"RNA triplet not found: {rna_dir}",
                    source_name=matrix.name,
                    individual_id=sample_id,
                    matrix_path=matrix if supported else None,
                    barcodes_path=barcodes if supported else None,
                    features_path=features if supported else None,
                )
            )
        return samples

    data_root = resolve_data_root(ROOT)
    selected = sorted(set(_load_manifest(gse)["gse"].astype(str)))
    samples = discover_rna_samples(data_root / "raw", selected)
    manifest_rna = _load_manifest(gse)
    allowed = set(
        manifest_rna[manifest_rna["assay"].astype(str).str.startswith("RNA")]["gsm"].astype(str)
    )
    return [sample for sample in samples if sample.sample_id in allowed]


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
            out / "matrix" / "matrix.mtx",
            out / "matrix" / "barcodes.tsv",
            out / "matrix" / "features.tsv",
        ]

    return [
        out / "metadata.csv",
        out / "metadata_qc.csv",
        out / "qc_summary.csv",
        out / "qc_overview.png",
        out / "umap_cima_cell_type_l1.png",
        out / "umap_cima_cell_type_l2.png",
        out / "matrix" / "matrix.mtx",
        out / "matrix" / "barcodes.tsv.gz",
        out / "matrix" / "features.tsv.gz",
        out / f"{sample.gsm}_seurat_qc.rds",
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
    if gse == "7555405":
        for sample in _discover_co_rna(gse):
            if sample.sample_id == sample_id:
                return sample
        raise FileNotFoundError(f"co RNA sample not found: {gse}/{sample_id}")
    return only_rna_cli.find_rna_sample(gse, sample_id)


def _print_discovery(rna_samples: Iterable[Any], atac_samples: Iterable[CoSample]) -> None:
    print("modality\tgse\tsample_id\tindividual_id\tsupported\tnote\tprimary_path")
    rna_manifest = _load_manifest(None)
    rna_individual: dict[tuple[str, str], str] = {}
    for _, row in rna_manifest.iterrows():
        if str(row.get("assay", "")).startswith("RNA"):
            rna_individual[(str(row.get("gse", "")), str(row.get("gsm", "")))] = str(row.get("individual_id", ""))
    for sample in rna_samples:
        primary = sample.archive_path or sample.h5_path or sample.matrix_path or sample.features_path
        individual_id = sample.individual_id or rna_individual.get((sample.gse, sample.sample_id), "")
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
    if args.gse == "7555405":
        return only_rna_cli._route_rna_sample(_find_co_rna_sample(args.gse, args.sample_id), args)
    return only_rna_cli.cmd_run_rna_sample(args)


def cmd_run_rna_gse(args) -> int:
    args.output_root = str(Path(getattr(args, "output_root", DEFAULT_OUTPUT_ROOT)) / "rna")
    if args.gse == "7555405":
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
    if sample.folder_name:
        command.extend(["--sample-label", sample.gsm])
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
