from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Iterable, cast

import h5py
import pandas as pd

from scripts.only_rna.discovery import resolve_data_root


ROOT = Path(__file__).resolve().parents[2]
DATASET_ID = "longevity"
DEFAULT_OUTPUT_ROOT = ROOT / "output" / DATASET_ID
ONLY_ATAC_R_SCRIPT = ROOT / "scripts" / "process" / "process_single_sample.R"
EXPORT_ATAC_H5AD = ROOT / "scripts" / "process" / "export_co_atac_h5ad.py"
PREPROCESS_ATAC_BARCODES_R = ROOT / "scripts" / "longevity" / "preprocess_atac_barcodes_archr.R"
LONGEVITY_RNA_MAPPING_VERSION = "longevity_l1_l2_to_5class_v1"
DEFAULT_BARCODE_MIN_FRAGMENTS = 200
DEFAULT_BARCODE_MAX_BARCODES = 20_000
DEFAULT_BARCODE_MIN_TSS = 2.5
HOMEBREW_LIB = Path("/home/linuxbrew/.linuxbrew/lib")
DEFAULT_PARAM_CONTRAST_ROOT = ROOT / "output" / "longevity_param_contrast"
PARAM_CONTRAST_CANDIDATES: dict[str, dict[str, float | int]] = {
    "baseline_50000_tss0": {"min_fragments": 200, "max_barcodes": 50000, "min_tss": 0.0},
    "max30000_tss0": {"min_fragments": 200, "max_barcodes": 30000, "min_tss": 0.0},
    "max30000_tss1p5": {"min_fragments": 200, "max_barcodes": 30000, "min_tss": 1.5},
    "max25000_tss2": {"min_fragments": 200, "max_barcodes": 25000, "min_tss": 2.0},
    "max30000_tss2": {"min_fragments": 200, "max_barcodes": 30000, "min_tss": 2.0},
    "max25000_tss2p5": {"min_fragments": 200, "max_barcodes": 25000, "min_tss": 2.5},
    "max20000_tss2": {"min_fragments": 200, "max_barcodes": 20000, "min_tss": 2.0},
    "max20000_tss2p5": {"min_fragments": 200, "max_barcodes": 20000, "min_tss": 2.5},
}


@dataclass(frozen=True)
class LongevityRnaAtlas:
    dataset: str
    sample_id: str
    h5ad_path: Path
    supported: bool = True
    note: str = "processed h5ad atlas"


@dataclass(frozen=True)
class LongevityAtacSample:
    dataset: str
    sample_id: str
    fragment_file: Path
    supported: bool = True
    note: str = "fragment file"


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def _raw_dir() -> Path:
    return resolve_data_root(ROOT) / "raw" / DATASET_ID


def _barcode_root(output_root: Path | None = None) -> Path:
    if output_root is not None:
        return Path(output_root)
    return resolve_data_root(ROOT) / "reference" / DATASET_ID / "atac_barcodes"


def _default_param_contrast_barcode_root() -> Path:
    return resolve_data_root(ROOT) / "reference" / DATASET_ID / "atac_barcodes_param_contrast"


def _barcode_sample_dir(sample_id: str, output_root: Path | None = None) -> Path:
    return _barcode_root(output_root) / sample_id


def barcode_file_path(sample_id: str, output_root: Path | None = None) -> Path:
    return _barcode_sample_dir(sample_id, output_root) / "filtered_barcodes.tsv.gz"


def barcode_qc_path(sample_id: str, output_root: Path | None = None) -> Path:
    return _barcode_sample_dir(sample_id, output_root) / "barcode_qc.csv.gz"


def barcode_summary_path(sample_id: str, output_root: Path | None = None) -> Path:
    return _barcode_sample_dir(sample_id, output_root) / "summary.json"


def discover_rna() -> list[LongevityRnaAtlas]:
    rna_dir = _raw_dir() / "rna"
    if not rna_dir.exists():
        return []
    return [
        LongevityRnaAtlas(
            dataset=DATASET_ID,
            sample_id=path.stem,
            h5ad_path=path,
        )
        for path in sorted(rna_dir.glob("*.h5ad"))
    ]


def discover_atac(sample_id: str | None = None) -> list[LongevityAtacSample]:
    atac_dir = _raw_dir() / "atac"
    if not atac_dir.exists():
        return []
    samples: list[LongevityAtacSample] = []
    for path in sorted(atac_dir.glob("*_fragments.tsv.gz")):
        current_sample_id = path.name.removesuffix("_fragments.tsv.gz")
        if sample_id is not None and current_sample_id != sample_id:
            continue
        samples.append(
            LongevityAtacSample(
                dataset=DATASET_ID,
                sample_id=current_sample_id,
                fragment_file=path,
            )
        )
    return samples


def _atac_output_dir(sample: LongevityAtacSample, output_root: Path) -> Path:
    return output_root / "atac" / sample.dataset / sample.sample_id


def _rna_output_dir(sample: LongevityRnaAtlas, output_root: Path) -> Path:
    return output_root / "rna" / sample.dataset / sample.sample_id


def _atac_status_file(sample: LongevityAtacSample, output_root: Path) -> Path:
    return _atac_output_dir(sample, output_root) / "run_status.json"


def _expected_atac_outputs(
    sample: LongevityAtacSample, output_root: Path, output_profile: str = "full"
) -> list[Path]:
    out = _atac_output_dir(sample, output_root)
    if output_profile in {"matrix-lite", "validation-lite"}:
        return [
            out / "umap_cima_cell_type_l1.png",
            out / "umap_cima_cell_type_l2.png",
            out / "qc_summary.csv",
            out / "validation_result.csv",
            out / f"{sample.sample_id}.h5ad",
        ]
    return [
        out / "metadata.csv",
        out / "metadata_qc.csv",
        out / "qc_summary.csv",
        out / "qc_overview.png",
        out / "umap_cima_cell_type_l1.png",
        out / "umap_cima_cell_type_l2.png",
        out / f"{sample.sample_id}.h5ad",
    ]


def _atac_outputs_complete(
    sample: LongevityAtacSample, output_root: Path, output_profile: str = "full"
) -> bool:
    return sample.supported and all(
        path.exists() for path in _expected_atac_outputs(sample, output_root, output_profile)
    )


def _nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _count_gzip_lines(path: Path) -> int:
    import gzip

    if not path.exists():
        return 0
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip())


def _override_path() -> Path:
    return resolve_data_root(ROOT) / "reference" / DATASET_ID / "atac_barcode_overrides.csv"


def _load_barcode_overrides(path: Path | None = None) -> dict[str, dict[str, str]]:
    path = path or _override_path()
    if not path.exists():
        return {}
    frame = pd.read_csv(path, dtype=str, keep_default_na=False).fillna("")
    if "sample_id" not in frame.columns:
        raise ValueError(f"Barcode override file missing sample_id column: {path}")
    return {
        str(row["sample_id"]): {str(k): str(v) for k, v in row.items()}
        for _, row in frame.iterrows()
        if str(row["sample_id"]).strip()
    }


def _coerce_int(value: object, default: int) -> int:
    text = "" if value is None else str(value).strip()
    if text == "":
        return default
    return int(float(text))


def _coerce_float(value: object, default: float) -> float:
    text = "" if value is None else str(value).strip()
    if text == "":
        return default
    return float(text)


def _barcode_thresholds_for_sample(sample_id: str, args) -> dict[str, float | int | str]:
    overrides = _load_barcode_overrides(getattr(args, "overrides", None))
    override = overrides.get(sample_id, {})
    min_fragments = _coerce_int(
        override.get("min_fragments"),
        int(getattr(args, "min_fragments", DEFAULT_BARCODE_MIN_FRAGMENTS)),
    )
    max_barcodes = _coerce_int(
        override.get("max_barcodes"),
        int(getattr(args, "max_barcodes", DEFAULT_BARCODE_MAX_BARCODES)),
    )
    min_tss = _coerce_float(
        override.get("min_tss"),
        float(getattr(args, "min_tss", DEFAULT_BARCODE_MIN_TSS)),
    )
    return {
        "min_fragments": min_fragments,
        "max_barcodes": max_barcodes,
        "min_tss": min_tss,
        "rank_by": override.get("rank_by", str(getattr(args, "rank_by", "fragments"))),
        "override_reason": override.get("reason", ""),
    }


def _barcode_status_row(sample: LongevityAtacSample, output_root: Path | None = None) -> dict[str, Any]:
    barcode_file = barcode_file_path(sample.sample_id, output_root)
    qc_file = barcode_qc_path(sample.sample_id, output_root)
    summary_file = barcode_summary_path(sample.sample_id, output_root)
    n_barcodes = _count_gzip_lines(barcode_file) if barcode_file.exists() else 0
    if not barcode_file.exists():
        status = "missing"
    elif n_barcodes == 0:
        status = "empty"
    elif not qc_file.exists() or not summary_file.exists():
        status = "incomplete"
    else:
        status = "ready"
    summary: dict[str, Any] = {}
    if summary_file.exists():
        summary = json.loads(summary_file.read_text(encoding="utf-8"))
    return {
        "sample_id": sample.sample_id,
        "fragment_file": str(sample.fragment_file),
        "barcode_file": str(barcode_file),
        "barcode_file_exists": barcode_file.exists(),
        "n_barcodes": n_barcodes,
        "barcode_qc_exists": qc_file.exists(),
        "summary_exists": summary_file.exists(),
        "min_fragments": summary.get("min_fragments", ""),
        "max_barcodes": summary.get("max_barcodes", ""),
        "min_tss": summary.get("min_tss", ""),
        "rank_by": summary.get("rank_by", ""),
        "status": status,
    }


def _same_numeric_value(left: object, right: object) -> bool:
    try:
        return abs(float(str(left)) - float(str(right))) < 1e-9
    except Exception:
        return False


def _barcode_status_matches_thresholds(
    status: dict[str, Any], thresholds: dict[str, float | int | str]
) -> bool:
    return (
        status.get("status") == "ready"
        and _same_numeric_value(status.get("min_fragments"), thresholds["min_fragments"])
        and _same_numeric_value(status.get("max_barcodes"), thresholds["max_barcodes"])
        and _same_numeric_value(status.get("min_tss"), thresholds["min_tss"])
        and str(status.get("rank_by") or "fragments") == str(thresholds.get("rank_by") or "fragments")
    )


def _run_barcode_thresholds_for_sample(
    sample_id: str, args
) -> dict[str, float | int | str]:
    threshold_args = SimpleNamespace(
        min_fragments=getattr(args, "barcode_min_fragments", DEFAULT_BARCODE_MIN_FRAGMENTS),
        max_barcodes=getattr(args, "barcode_max_barcodes", DEFAULT_BARCODE_MAX_BARCODES),
        min_tss=getattr(args, "barcode_min_tss", DEFAULT_BARCODE_MIN_TSS),
        rank_by=getattr(args, "barcode_rank_by", "fragments"),
        overrides=getattr(args, "barcode_overrides", None),
    )
    return _barcode_thresholds_for_sample(sample_id, threshold_args)


def _print_rows(rows: list[dict[str, Any]], output_format: str) -> None:
    frame = pd.DataFrame(rows)
    if output_format == "json":
        print(json.dumps(rows, indent=2, ensure_ascii=False))
    elif output_format == "csv":
        print(frame.to_csv(index=False), end="")
    else:
        print(frame.to_string(index=False))


def _parse_percent_or_float(value: object) -> float:
    text = "" if value is None else str(value).strip()
    if text == "":
        return 0.0
    if text.endswith("%"):
        return float(text[:-1]) / 100.0
    return float(text)


def _read_metric_csv(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not {"metric", "value"}.issubset(frame.columns):
        raise ValueError(f"Expected metric/value columns in {path}")
    return {str(row["metric"]): str(row["value"]) for _, row in frame.iterrows()}


def _candidate_thresholds_from_barcode_summary(
    sample_id: str, candidate: str, candidate_dir: Path
) -> dict[str, Any]:
    if candidate in PARAM_CONTRAST_CANDIDATES:
        thresholds = PARAM_CONTRAST_CANDIDATES[candidate]
        return {
            "barcode_min_fragments": thresholds["min_fragments"],
            "barcode_max_barcodes": thresholds["max_barcodes"],
            "barcode_min_tss": thresholds["min_tss"],
            "barcode_rank_by": "fragments",
        }
    summary_path = _default_param_contrast_barcode_root() / candidate / sample_id / "summary.json"
    if not summary_path.exists():
        summary_path = candidate_dir / "reference" / DATASET_ID / "atac_barcodes" / sample_id / "summary.json"
    if not summary_path.exists():
        return {
            "barcode_min_fragments": "",
            "barcode_max_barcodes": "",
            "barcode_min_tss": "",
            "barcode_rank_by": "",
        }
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    return {
        "barcode_min_fragments": summary.get("min_fragments", ""),
        "barcode_max_barcodes": summary.get("max_barcodes", ""),
        "barcode_min_tss": summary.get("min_tss", ""),
        "barcode_rank_by": summary.get("rank_by", "fragments"),
    }


def _selection_status_payload(row: pd.Series, source: Path) -> dict[str, Any]:
    def _row_value(name: str) -> str:
        return str(row.get(name, "")).strip()

    def _int_or_none(name: str) -> int | None:
        value = _row_value(name)
        return int(float(value)) if value else None

    def _float_or_none(name: str) -> float | None:
        value = _row_value(name)
        return float(value) if value else None

    candidate = _row_value("candidate")
    sample_id = _row_value("sample_id")
    barcode_summary_path = _default_param_contrast_barcode_root() / candidate / sample_id / "summary.json"
    barcode_summary: dict[str, Any] = {}
    if barcode_summary_path.exists():
        barcode_summary = json.loads(barcode_summary_path.read_text(encoding="utf-8"))

    barcode_parameters = {
        "min_fragments": barcode_summary.get("min_fragments", _int_or_none("barcode_min_fragments")),
        "max_barcodes": barcode_summary.get("max_barcodes", _int_or_none("barcode_max_barcodes")),
        "min_tss": barcode_summary.get("min_tss", _float_or_none("barcode_min_tss")),
        "rank_by": barcode_summary.get("rank_by", _row_value("barcode_rank_by") or "fragments"),
    }
    return {
        "published_from_param_contrast": str(source),
        "published_candidate": candidate,
        "published_selection": {
            "candidate": candidate,
            "source_output_dir": str(source),
            "selection_score": _float_or_none("selection_score"),
            "passes_selection_targets": _row_value("passes_selection_targets").lower() == "true",
            "barcode_parameters": barcode_parameters,
            "metrics": {
                "input_cells": _int_or_none("input_cells"),
                "pass_qc": _int_or_none("pass_qc"),
                "qc_rate": _float_or_none("qc_rate"),
                "median_TSS_enrichment": _float_or_none("median_TSS_enrichment"),
                "median_FRiP": _float_or_none("median_FRiP"),
                "median_fragments": _float_or_none("median_fragments"),
                "cima_l1_low_confidence_frac": _float_or_none("cima_l1_low_confidence_frac"),
                "median_cima_l4_score": _float_or_none("median_cima_l4_score"),
            },
        },
    }


def _score_param_contrast_row(row: dict[str, Any], min_pass_qc: int, min_frip: float) -> float:
    if str(row.get("run_status", "")).lower() != "success" or not bool(row.get("outputs_complete", False)):
        return -1.0
    pass_qc = int(float(row.get("pass_qc", 0) or 0))
    if pass_qc < min_pass_qc:
        return -1.0
    low_conf = float(row.get("cima_l1_low_confidence_frac", 0.0) or 0.0)
    cima_l4 = float(row.get("median_cima_l4_score", 0.0) or 0.0)
    frip = float(row.get("median_FRiP", 0.0) or 0.0)
    if frip < min_frip:
        return -1.0
    qc_rate = float(row.get("qc_rate", 0.0) or 0.0)
    pass_qc_score = min(pass_qc / 15000.0, 1.0)
    return (
        0.35 * (1.0 - low_conf)
        + 0.25 * cima_l4
        + 0.20 * frip
        + 0.10 * qc_rate
        + 0.10 * pass_qc_score
    )


def _collect_param_contrast_rows(contrast_root: Path, min_pass_qc: int, min_frip: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for qc_summary in sorted(contrast_root.glob("*/atac/longevity/*/qc_summary.csv")):
        sample_dir = qc_summary.parent
        sample_id = sample_dir.name
        candidate = qc_summary.parents[3].name
        candidate_dir = contrast_root / candidate
        metrics = _read_metric_csv(qc_summary)
        run_status_path = sample_dir / "run_status.json"
        status_payload: dict[str, Any] = {}
        if run_status_path.exists():
            status_payload = json.loads(run_status_path.read_text(encoding="utf-8"))
        outputs_complete = all(
            (sample_dir / name).exists()
            for name in [
                "metadata.csv",
                "metadata_qc.csv",
                "qc_summary.csv",
                "qc_overview.png",
                "umap_cima_cell_type_l1.png",
                "umap_cima_cell_type_l2.png",
                f"{sample_id}.h5ad",
            ]
        )
        row: dict[str, Any] = {
            "sample_id": sample_id,
            "candidate": candidate,
            "candidate_output_dir": str(sample_dir),
            "run_status": status_payload.get("status", "success" if outputs_complete else "unknown"),
            "outputs_complete": outputs_complete,
            "input_cells": int(float(metrics.get("input_cells", 0) or 0)),
            "pass_qc": int(float(metrics.get("pass_qc", 0) or 0)),
            "singlets": int(float(metrics.get("singlets", 0) or 0)),
            "doublets": int(float(metrics.get("doublets", 0) or 0)),
            "query_cluster_count": int(float(metrics.get("query_cluster_count", 0) or 0)),
            "qc_rate": _parse_percent_or_float(metrics.get("qc_rate", "0")),
            "doublet_fraction": int(float(metrics.get("doublets", 0) or 0))
            / max(int(float(metrics.get("input_cells", 0) or 0)), 1),
            "median_TSS_enrichment": float(metrics.get("median_TSS_enrichment", 0) or 0),
            "median_FRiP": float(metrics.get("median_FRiP", 0) or 0),
            "median_fragments": float(metrics.get("median_fragments", 0) or 0),
            "cima_l1_low_confidence_frac": _parse_percent_or_float(
                metrics.get("cima_l1_low_confidence_frac", "0")
            ),
            "median_cima_l4_score": float(metrics.get("median_cima_l4_score", 0) or 0),
        }
        row.update(_candidate_thresholds_from_barcode_summary(sample_id, candidate, candidate_dir))
        row["passes_selection_targets"] = row["pass_qc"] >= min_pass_qc and row["median_FRiP"] >= min_frip
        row["selection_score"] = _score_param_contrast_row(row, min_pass_qc, min_frip)
        rows.append(row)
    return rows


def _select_best_param_contrast(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    selected: list[dict[str, Any]] = []
    for sample_id, group in pd.DataFrame(rows).groupby("sample_id", sort=True):
        group = group.sort_values(
            ["selection_score", "cima_l1_low_confidence_frac", "median_cima_l4_score", "pass_qc"],
            ascending=[False, True, False, False],
        )
        best = cast(dict[str, Any], group.iloc[0].to_dict())
        best["selected"] = True
        selected.append(best)
    return selected


def _load_status(sample: LongevityAtacSample, output_root: Path) -> dict[str, Any] | None:
    path = _atac_status_file(sample, output_root)
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def _write_status(sample: LongevityAtacSample, output_root: Path, payload: dict[str, Any]) -> None:
    output_dir = _atac_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    _atac_status_file(sample, output_root).write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def _subprocess_env() -> dict[str, str]:
    env = os.environ.copy()
    if HOMEBREW_LIB.exists():
        existing = env.get("LD_LIBRARY_PATH", "")
        paths = [str(HOMEBREW_LIB)]
        paths.extend(path for path in existing.split(":") if path)
        env["LD_LIBRARY_PATH"] = ":".join(dict.fromkeys(paths))
    return env


def _run_command(command: list[str], log_file: Path) -> int:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.Popen(
        command,
        cwd=ROOT,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        env=_subprocess_env(),
    )
    assert process.stdout is not None
    with log_file.open("w", encoding="utf-8") as handle:
        for line in process.stdout:
            sys.stdout.write(line)
            handle.write(line)
            handle.flush()
    return process.wait()


def _print_discovery(rna_samples: Iterable[LongevityRnaAtlas], atac_samples: Iterable[LongevityAtacSample]) -> None:
    print("modality\tdataset\tsample_id\tsupported\tnote\tprimary_path")
    for sample in rna_samples:
        print(
            f"RNA\t{sample.dataset}\t{sample.sample_id}\t"
            f"{str(sample.supported).lower()}\t{sample.note}\t{sample.h5ad_path}"
        )
    for sample in atac_samples:
        print(
            f"ATAC\t{sample.dataset}\t{sample.sample_id}\t"
            f"{str(sample.supported).lower()}\t{sample.note}\t{sample.fragment_file}"
        )


def _print_status(atac_samples: Iterable[LongevityAtacSample], output_root: Path) -> None:
    print("modality\tdataset\tsample_id\tstatus\toutputs_complete\tlast_finished")
    for sample in atac_samples:
        payload = _load_status(sample, output_root) or {}
        status = payload.get("status", "pending")
        output_profile = payload.get("output_profile", "full")
        finished_at = payload.get("finished_at", "")
        print(
            f"ATAC\t{sample.dataset}\t{sample.sample_id}\t{status}\t"
            f"{str(_atac_outputs_complete(sample, output_root, output_profile)).lower()}\t{finished_at}"
        )


def _find_atac_sample(sample_id: str) -> LongevityAtacSample:
    for sample in discover_atac(sample_id):
        return sample
    raise FileNotFoundError(f"Longevity ATAC sample not found: {sample_id}")


def _find_rna_atlas(sample_id: str | None = None) -> LongevityRnaAtlas:
    samples = discover_rna()
    if sample_id is not None:
        samples = [sample for sample in samples if sample.sample_id == sample_id]
    if len(samples) != 1:
        names = ", ".join(sample.sample_id for sample in samples) or "none"
        raise FileNotFoundError(f"Expected one longevity RNA atlas, found: {names}")
    return samples[0]


def _clean_label(value: object) -> str:
    if value is None:
        return ""
    if not isinstance(value, (str, int, float, bool)):
        return str(value).strip()
    if pd.isna(value):
        return ""
    return str(value).strip()


def _infer_longevity_final_celltype(l1: object, l2: object) -> str:
    l2_clean = _clean_label(l2)
    if l2_clean.startswith("NK_"):
        return "NK"
    if l2_clean.startswith("CD4T_"):
        return "CD4_T"
    if l2_clean.startswith("CD8T_") or l2_clean == "CD8T_gdT":
        return "CD8_T"
    if l2_clean.startswith("B_"):
        return "B"
    if l2_clean.startswith("MYE_"):
        return "Myeloid"

    return {
        "NK": "NK",
        "CD4T": "CD4_T",
        "CD8T": "CD8_T",
        "CD8Tn": "CD8_T",
        "B": "B",
        "Bn": "B",
        "Plasma": "B",
        "CD14mono": "Myeloid",
        "CD16mono": "Myeloid",
        "MK": "Myeloid",
        "cDC": "Myeloid",
        "pDC": "Myeloid",
    }.get(_clean_label(l1), "Unknown")


def _add_longevity_final_celltype(obs: pd.DataFrame) -> pd.DataFrame:
    if "L1_annotation" not in obs.columns or "L2_annotation_new" not in obs.columns:
        raise KeyError("Longevity RNA h5ad obs must contain L1_annotation and L2_annotation_new")
    out = obs.copy()
    out["final_celltype"] = pd.Series(
        [
            _infer_longevity_final_celltype(l1, l2)
            for l1, l2 in zip(out["L1_annotation"], out["L2_annotation_new"], strict=True)
        ],
        index=out.index,
        dtype="object",
    )
    return out


def _write_longevity_comparison(obs: pd.DataFrame, output_path: Path) -> None:
    comparison = obs.groupby(
        ["L1_annotation", "L2_annotation_new", "final_celltype"], dropna=False
    ).size()
    comparison = comparison.reset_index()
    comparison = comparison.rename(columns={0: "n_cells"}).sort_values(
        ["L1_annotation", "L2_annotation_new", "final_celltype"]
    )
    comparison.to_csv(output_path, index=False)


def _read_categorical_from_obs(obs_group: h5py.Group, column: str) -> pd.Series:
    group = cast(h5py.Group, obs_group[column])
    categories = [value.decode() if isinstance(value, bytes) else str(value) for value in group["categories"][...]]
    codes = group["codes"][...]
    values = [categories[int(code)] if int(code) >= 0 else "" for code in codes]
    return pd.Series(values, dtype="object")


def _write_obs_categorical(file_path: Path, column: str, values: pd.Series) -> None:
    categories = sorted(set(str(value) for value in values.tolist()))
    category_index = {category: idx for idx, category in enumerate(categories)}
    codes = values.astype(str).map(category_index).astype("int8").to_numpy()
    string_dtype = h5py.string_dtype(encoding="utf-8")
    with h5py.File(file_path, "r+") as handle:
        obs = cast(h5py.Group, handle["obs"])
        if column in obs:
            del obs[column]
        group = obs.create_group(column)
        group.attrs["encoding-type"] = "categorical"
        group.attrs["encoding-version"] = "0.2.0"
        group.attrs["ordered"] = False
        group.create_dataset("categories", data=[category.encode("utf-8") for category in categories], dtype=string_dtype)
        group["categories"].attrs["encoding-type"] = "string-array"
        group["categories"].attrs["encoding-version"] = "0.2.0"
        group.create_dataset("codes", data=codes)
        group["codes"].attrs["encoding-type"] = "array"
        group["codes"].attrs["encoding-version"] = "0.2.0"
        order = [str(value) for value in list(obs.attrs["column-order"])]
        if column not in order:
            order.append(column)
        obs.attrs["column-order"] = order


def _prepare_longevity_obs_from_h5ad(h5ad_path: Path) -> pd.DataFrame:
    with h5py.File(h5ad_path, "r") as handle:
        obs = cast(h5py.Group, handle["obs"])
        frame = pd.DataFrame(
            {
                "L1_annotation": _read_categorical_from_obs(obs, "L1_annotation"),
                "L2_annotation_new": _read_categorical_from_obs(obs, "L2_annotation_new"),
            }
        )
    return _add_longevity_final_celltype(frame)


def cmd_ingest_rna(args) -> int:
    sample = _find_rna_atlas(getattr(args, "sample_id", None))
    output_root = Path(args.output_root)
    output_dir = _rna_output_dir(sample, output_root)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_h5ad = output_dir / sample.h5ad_path.name
    comparison_path = output_dir / "celltype_original_vs_final.csv"

    obs = _prepare_longevity_obs_from_h5ad(sample.h5ad_path)
    counts = obs["final_celltype"].value_counts(dropna=False).sort_index().to_dict()
    print(f"[longevity-ingest-rna] {sample.sample_id}")
    print(f"input={sample.h5ad_path}")
    print(f"output={output_h5ad}")
    print(f"comparison={comparison_path}")
    print(f"final_celltype_counts={counts}")
    if args.dry_run:
        return 0

    shutil.copy2(sample.h5ad_path, output_h5ad)
    _write_obs_categorical(output_h5ad, "final_celltype", cast(pd.Series, obs["final_celltype"]))
    _write_longevity_comparison(obs, comparison_path)
    return 0


def _preprocess_atac_barcodes_sample(sample: LongevityAtacSample, args) -> int:
    output_root = _barcode_root(Path(args.output_root) if getattr(args, "output_root", None) else None)
    output_dir = _barcode_sample_dir(sample.sample_id, output_root)
    thresholds = _barcode_thresholds_for_sample(sample.sample_id, args)
    command = [
        args.rscript,
        str(PREPROCESS_ATAC_BARCODES_R),
        "--fragment-file",
        str(sample.fragment_file),
        "--sample-id",
        sample.sample_id,
        "--output-dir",
        str(output_dir),
        "--min-fragments",
        str(thresholds["min_fragments"]),
        "--max-barcodes",
        str(thresholds["max_barcodes"]),
        "--min-tss",
        str(thresholds["min_tss"]),
        "--rank-by",
        str(thresholds["rank_by"]),
    ]
    if args.force:
        command.append("--force")
    print(f"[longevity-preprocess-atac-barcodes] {sample.sample_id}")
    if thresholds.get("override_reason"):
        print(f"override_reason={thresholds['override_reason']}")
    print(" ".join(command))
    if args.dry_run:
        return 0
    output_dir.mkdir(parents=True, exist_ok=True)
    return _run_command(command, output_dir / "preprocess_archr.log")


def _ensure_atac_barcodes(sample: LongevityAtacSample, args) -> tuple[int, Path | None]:
    barcode_root = Path(getattr(args, "barcode_output_root", _barcode_root()))
    thresholds = _run_barcode_thresholds_for_sample(sample.sample_id, args)
    status = _barcode_status_row(sample, barcode_root)
    if _barcode_status_matches_thresholds(status, thresholds):
        barcode_path = barcode_file_path(sample.sample_id, barcode_root)
        print(f"[longevity-atac-barcodes] {sample.sample_id} ready: {barcode_path}")
        return 0, barcode_path
    if status["status"] == "ready":
        print(
            f"[longevity-atac-barcodes] {sample.sample_id} reprocessing: "
            f"existing thresholds min_fragments={status.get('min_fragments')} "
            f"max_barcodes={status.get('max_barcodes')} min_tss={status.get('min_tss')} "
            f"rank_by={status.get('rank_by') or 'fragments'} "
            f"do not match requested min_fragments={thresholds['min_fragments']} "
            f"max_barcodes={thresholds['max_barcodes']} min_tss={thresholds['min_tss']} "
            f"rank_by={thresholds['rank_by']}"
        )

    preprocess_args = SimpleNamespace(
        rscript=getattr(args, "rscript", "Rscript"),
        output_root=barcode_root,
        min_fragments=thresholds["min_fragments"],
        max_barcodes=thresholds["max_barcodes"],
        min_tss=thresholds["min_tss"],
        rank_by=thresholds["rank_by"],
        overrides=getattr(args, "barcode_overrides", None),
        force=True,
        dry_run=getattr(args, "dry_run", False),
    )
    current = _preprocess_atac_barcodes_sample(sample, preprocess_args)
    if current != 0:
        return current, None
    barcode_path = barcode_file_path(sample.sample_id, barcode_root)
    if getattr(args, "dry_run", False):
        return 0, barcode_path if _nonempty_file(barcode_path) else None
    status = _barcode_status_row(sample, barcode_root)
    if not _barcode_status_matches_thresholds(status, thresholds):
        print(f"[longevity-atac-barcodes] {sample.sample_id} not ready after preprocessing: {status['status']}")
        return 1, None
    return 0, barcode_path


def cmd_preprocess_atac_barcodes(args) -> int:
    samples = discover_atac(getattr(args, "sample_id", None))
    if getattr(args, "sample_id", None) and not samples:
        raise FileNotFoundError(f"Longevity ATAC sample not found: {args.sample_id}")
    returncode = 0
    for sample in samples:
        current = _preprocess_atac_barcodes_sample(sample, args)
        if current != 0:
            returncode = current
            if not getattr(args, "keep_going", False):
                break
    return returncode


def cmd_atac_barcode_status(args) -> int:
    output_root = Path(args.output_root) if getattr(args, "output_root", None) else None
    rows = [_barcode_status_row(sample, output_root) for sample in discover_atac(getattr(args, "sample_id", None))]
    _print_rows(rows, args.format)
    return 0


def cmd_summarize_atac_param_contrast(args) -> int:
    contrast_root = Path(args.contrast_root)
    rows = _collect_param_contrast_rows(contrast_root, args.min_pass_qc, args.min_frip)
    if not rows:
        raise FileNotFoundError(f"No longevity ATAC param contrast qc_summary.csv files under {contrast_root}")
    all_frame = pd.DataFrame(rows).sort_values(["sample_id", "selection_score"], ascending=[True, False])
    selected_frame = pd.DataFrame(_select_best_param_contrast(rows)).sort_values("sample_id")
    if getattr(args, "write", False):
        contrast_root.mkdir(parents=True, exist_ok=True)
        all_frame.to_csv(contrast_root / "param_contrast_summary_all_samples.csv", index=False)
        selected_frame.to_csv(contrast_root / "param_contrast_selected_samples.csv", index=False)
        print(f"wrote={contrast_root / 'param_contrast_summary_all_samples.csv'}")
        print(f"wrote={contrast_root / 'param_contrast_selected_samples.csv'}")
    if args.format == "json":
        print(selected_frame.to_json(orient="records", indent=2))
    elif args.format == "csv":
        print(selected_frame.to_csv(index=False), end="")
    else:
        print(selected_frame.to_string(index=False))
    return 0


def cmd_publish_atac_param_contrast(args) -> int:
    selection_csv = Path(args.selection_csv)
    if not selection_csv.exists():
        raise FileNotFoundError(f"Selection CSV not found: {selection_csv}")
    frame = pd.read_csv(selection_csv, dtype=str, keep_default_na=False)
    required = {"sample_id", "candidate_output_dir"}
    if not required.issubset(frame.columns):
        missing = ", ".join(sorted(required - set(frame.columns)))
        raise ValueError(f"Selection CSV missing required columns: {missing}")
    output_root = Path(args.output_root)
    for _, row in frame.sort_values("sample_id").iterrows():
        sample_id = str(row["sample_id"])
        source = Path(str(row["candidate_output_dir"]))
        dest = output_root / "atac" / DATASET_ID / sample_id
        if not source.exists():
            raise FileNotFoundError(f"Selected source directory not found for {sample_id}: {source}")
        print(f"[publish-longevity-atac] {sample_id}: {source} -> {dest}")
        if args.dry_run:
            continue
        if dest.exists():
            if not args.force:
                raise FileExistsError(f"Destination exists, use --force to replace: {dest}")
            shutil.rmtree(dest)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(source, dest)
        status_path = dest / "run_status.json"
        status_payload: dict[str, Any] = {}
        if status_path.exists():
            status_payload = json.loads(status_path.read_text(encoding="utf-8"))
        status_payload.update(
            {
                "status": "success",
                "published_at": utc_now(),
                "outputs_complete": _atac_outputs_complete(_find_atac_sample(sample_id), output_root),
            }
        )
        status_payload.update(_selection_status_payload(row, source))
        status_path.write_text(json.dumps(status_payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


def _param_contrast_candidate_names(args) -> list[str]:
    requested = getattr(args, "candidate", None) or []
    if not requested:
        return list(PARAM_CONTRAST_CANDIDATES)
    unknown = [name for name in requested if name not in PARAM_CONTRAST_CANDIDATES]
    if unknown:
        raise ValueError(f"Unknown longevity ATAC param contrast candidate(s): {', '.join(unknown)}")
    return list(requested)


def _custom_candidate_name(min_fragments: int, max_barcodes: int, min_tss: float) -> str:
    tss_text = str(min_tss).replace(".", "p")
    return f"custom_min{min_fragments}_max{max_barcodes}_tss{tss_text}"


def cmd_run_atac_custom_param_contrast(args) -> int:
    sample = _find_atac_sample(args.sample_id)
    candidate = args.candidate_name or _custom_candidate_name(
        args.min_fragments, args.max_barcodes, args.min_tss
    )
    output_root = Path(args.output_root) / candidate
    barcode_root = Path(args.barcode_output_root) / candidate
    if args.skip_complete and _param_contrast_sample_complete(sample, output_root):
        print(f"[skip] {candidate}/{sample.sample_id} already complete")
        return 0
    run_args = SimpleNamespace(
        sample_id=sample.sample_id,
        rscript=args.rscript,
        output_profile="full",
        output_root=str(output_root),
        nmads=args.nmads,
        min_inferred_fragments=None,
        max_inferred_barcodes=None,
        umap_min_dist=getattr(args, "umap_min_dist", None),
        barcode_output_root=str(barcode_root),
        barcode_min_fragments=args.min_fragments,
        barcode_max_barcodes=args.max_barcodes,
        barcode_min_tss=args.min_tss,
        barcode_rank_by=args.rank_by,
        barcode_overrides=None,
        force=True,
        dry_run=args.dry_run,
    )
    print(
        f"[longevity-custom-param-contrast] sample={sample.sample_id} candidate={candidate} "
        f"min_fragments={args.min_fragments} max_barcodes={args.max_barcodes} "
        f"min_tss={args.min_tss} rank_by={args.rank_by}"
    )
    return _run_atac_sample(sample, run_args)


def _param_contrast_sample_complete(sample: LongevityAtacSample, output_root: Path) -> bool:
    return _atac_outputs_complete(sample, output_root, "full")


def cmd_run_atac_param_contrast(args) -> int:
    sample_ids = set(getattr(args, "sample_id", None) or [])
    samples = [sample for sample in discover_atac() if not sample_ids or sample.sample_id in sample_ids]
    missing = sample_ids - {sample.sample_id for sample in samples}
    if missing:
        raise FileNotFoundError(f"Longevity ATAC sample(s) not found: {', '.join(sorted(missing))}")
    candidates = _param_contrast_candidate_names(args)
    returncode = 0
    for sample in samples:
        for candidate in candidates:
            thresholds = PARAM_CONTRAST_CANDIDATES[candidate]
            output_root = Path(args.output_root) / candidate
            barcode_root = Path(args.barcode_output_root) / candidate
            if args.skip_complete and _param_contrast_sample_complete(sample, output_root):
                print(f"[skip] {candidate}/{sample.sample_id} already complete")
                continue
            run_args = SimpleNamespace(
                sample_id=sample.sample_id,
                rscript=args.rscript,
                output_profile="full",
                output_root=str(output_root),
                nmads=args.nmads,
                min_inferred_fragments=None,
                max_inferred_barcodes=None,
                umap_min_dist=getattr(args, "umap_min_dist", None),
                barcode_output_root=str(barcode_root),
                barcode_min_fragments=int(thresholds["min_fragments"]),
                barcode_max_barcodes=int(thresholds["max_barcodes"]),
                barcode_min_tss=float(thresholds["min_tss"]),
                barcode_rank_by="fragments",
                barcode_overrides=None,
                force=True,
                dry_run=args.dry_run,
            )
            print(
                f"[longevity-param-contrast] sample={sample.sample_id} candidate={candidate} "
                f"min_fragments={thresholds['min_fragments']} "
                f"max_barcodes={thresholds['max_barcodes']} min_tss={thresholds['min_tss']}"
            )
            current = _run_atac_sample(sample, run_args)
            if current != 0:
                returncode = current
                if not args.keep_going:
                    return returncode
    return returncode


def _run_atac_sample(sample: LongevityAtacSample, args) -> int:
    output_root = Path(args.output_root)
    output_profile = args.output_profile
    if _atac_outputs_complete(sample, output_root, output_profile) and not args.force:
        print(f"[skip] {sample.dataset}/{sample.sample_id} already has complete outputs")
        return 0

    barcode_returncode, generated_barcode = _ensure_atac_barcodes(sample, args)
    if barcode_returncode != 0:
        status_payload = {
            "sample": {"dataset": sample.dataset, "sample_id": sample.sample_id},
            "output_profile": output_profile,
            "output_root": str(output_root),
            "started_at": utc_now(),
            "finished_at": utc_now(),
            "returncode": barcode_returncode,
            "outputs_complete": False,
            "status": "failed",
            "failure_stage": "barcode_preprocessing",
        }
        _write_status(sample, output_root, status_payload)
        return barcode_returncode

    command = [
        args.rscript,
        str(ONLY_ATAC_R_SCRIPT),
        "--gse",
        sample.dataset,
        "--gsm",
        sample.sample_id,
        "--nmads",
        str(args.nmads),
        "--output-profile",
        output_profile,
        "--output-root",
        str(output_root / "atac"),
        "--fragment-file",
        str(sample.fragment_file),
        "--sample-label",
        sample.sample_id,
    ]
    if generated_barcode is not None and _nonempty_file(generated_barcode):
        command.extend(["--barcode-file", str(generated_barcode)])
    if getattr(args, "min_inferred_fragments", None) is not None:
        command.extend(["--min-inferred-fragments", str(args.min_inferred_fragments)])
    if getattr(args, "max_inferred_barcodes", None) is not None:
        command.extend(["--max-inferred-barcodes", str(args.max_inferred_barcodes)])
    if getattr(args, "umap_min_dist", None) is not None:
        command.extend(["--umap-min-dist", str(args.umap_min_dist)])

    status_payload: dict[str, Any] = {
        "sample": {"dataset": sample.dataset, "sample_id": sample.sample_id},
        "command": command,
        "output_profile": output_profile,
        "output_root": str(output_root),
        "started_at": utc_now(),
        "status": "running",
        "barcode_file": str(generated_barcode) if generated_barcode is not None else "",
    }
    _write_status(sample, output_root, status_payload)

    print(f"[longevity-run-atac] {sample.dataset}/{sample.sample_id}")
    print(" ".join(command))
    if args.dry_run:
        status_payload["status"] = "dry_run"
        status_payload["finished_at"] = utc_now()
        _write_status(sample, output_root, status_payload)
        return 0

    log_file = _atac_output_dir(sample, output_root) / "logs" / "sample_qc.log"
    returncode = _run_command(command, log_file)
    if returncode == 0:
        export_command = [
            sys.executable,
            str(EXPORT_ATAC_H5AD),
            "--sample-dir",
            str(_atac_output_dir(sample, output_root)),
            "--data-root",
            str(resolve_data_root(ROOT)),
            "--overwrite",
            "--cleanup",
        ]
        with log_file.open("a", encoding="utf-8") as handle:
            handle.write("\n[export-longevity-atac-h5ad]\n")
            handle.write(" ".join(export_command) + "\n")
        export_returncode = _run_command(export_command, log_file)
        if export_returncode != 0:
            returncode = export_returncode

    status_payload["finished_at"] = utc_now()
    status_payload["returncode"] = returncode
    status_payload["outputs_complete"] = _atac_outputs_complete(sample, output_root, output_profile)
    status_payload["status"] = (
        "success" if returncode == 0 and status_payload["outputs_complete"] else "failed"
    )
    _write_status(sample, output_root, status_payload)
    return returncode


def cmd_discover(args) -> int:
    _print_discovery(discover_rna(), discover_atac(getattr(args, "sample_id", None)))
    return 0


def cmd_status(args) -> int:
    _print_status(discover_atac(getattr(args, "sample_id", None)), Path(args.output_root))
    return 0


def cmd_run_atac_sample(args) -> int:
    return _run_atac_sample(_find_atac_sample(args.sample_id), args)


def cmd_run_atac_all(args) -> int:
    returncode = 0
    for sample in discover_atac():
        current = _run_atac_sample(sample, args)
        if current != 0:
            returncode = current
    return returncode
