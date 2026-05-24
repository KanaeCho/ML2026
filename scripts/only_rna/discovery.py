from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Literal


DEFAULT_DATA_ROOT = Path("/mnt/g/ML2026_data")
ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class DiscoveredSample:
    gse: str
    sample_id: str
    input_type: str
    sample_kind: Literal["gsm", "gse_shared"]
    supported: bool
    note: str
    source_name: str
    individual_id: str = ""
    age: str = ""
    health: str = ""
    donor: str = ""
    matrix_path: Path | None = None
    barcodes_path: Path | None = None
    features_path: Path | None = None
    h5_path: Path | None = None
    archive_path: Path | None = None


def resolve_data_root(cwd: Path | None = None) -> Path:
    workspace_root = Path.cwd() if cwd is None else Path(cwd)
    workspace_data = workspace_root / "data"
    if workspace_data.exists():
        return workspace_data

    candidates: list[Path] = []
    env_root = os.environ.get("ML2026_DATA_ROOT")
    if env_root:
        candidates.append(Path(env_root))
    candidates.append(DEFAULT_DATA_ROOT)

    for candidate in candidates:
        if candidate.exists():
            return candidate

    searched = ", ".join(str(path) for path in [workspace_data, *candidates])
    raise FileNotFoundError(f"Unable to resolve data root. Searched: {searched}")


def selected_rna_gses(reference_root: Path) -> list[str]:
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from scripts.process.download_from_datasets import load_filtered_rows

    rows = load_filtered_rows(
        reference_root / "datasets.xlsx",
        assay="scRNA",
        visible_rows_only=True,
    )
    selected = {row.gse for row in rows}

    co2_manifest = reference_root / "co2_sample_manifest.csv"
    if co2_manifest.exists():
        import pandas as pd

        manifest = pd.read_csv(co2_manifest)
        assay_series = manifest.get("assay", pd.Series(dtype=object)).astype(str)
        gse_series = manifest.get("gse", pd.Series(dtype=object)).astype(str)
        selected.update(
            gse.strip()
            for gse, assay in zip(gse_series, assay_series, strict=False)
            if gse.strip() and assay.strip().startswith("RNA")
        )

    return sorted(selected)


def discover_rna_samples(
    raw_root: Path, selected_gses: list[str]
) -> list[DiscoveredSample]:
    samples: list[DiscoveredSample] = []
    individual_map: dict[tuple[str, str], str] = {}

    co2_manifest = ROOT / "data" / "reference" / "co2_sample_manifest.csv"
    if co2_manifest.exists():
        import pandas as pd

        manifest = pd.read_csv(co2_manifest)
        for _, row in manifest.iterrows():
            gse = str(row.get("gse", "")).strip()
            gsm = str(row.get("gsm", "")).strip()
            individual_id = str(row.get("individual_id", "")).strip()
            assay = str(row.get("assay", "")).strip()
            if gse and gsm and individual_id and assay.startswith("RNA"):
                individual_map[(gse, gsm)] = individual_id

    for gse in sorted(selected_gses):
        gse_dir = raw_root / gse
        if not gse_dir.exists():
            continue

        entries = sorted(path for path in gse_dir.iterdir() if path.is_file())
        supported_by_sample: dict[str, DiscoveredSample] = {}
        unsupported_by_sample: dict[str, DiscoveredSample] = {}
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
                        DiscoveredSample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="h5",
                            sample_kind="gsm",
                            supported=True,
                            note="10x h5",
                            source_name=path.name,
                            individual_id=individual_map.get((gse, sample_id), ""),
                            h5_path=path,
                        ),
                    )
                    continue
                if path.name.endswith(".tar.gz") and "matrix" in path.name.lower():
                    supported_by_sample.setdefault(
                        sample_id,
                        DiscoveredSample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="archive",
                            sample_kind="gsm",
                            supported=True,
                            note="matrix tar archive",
                            source_name=path.name,
                            individual_id=individual_map.get((gse, sample_id), ""),
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
                        DiscoveredSample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="shared-gene-count",
                            sample_kind="gse_shared",
                            supported=False,
                            note="Unsupported shared gene-count matrix for single-cell RNA workflow",
                            source_name=path.name,
                            matrix_path=path,
                        ),
                    )

        for sample_id, parts in sorted(triplet_parts.items()):
            if {"matrix", "barcodes", "features"} <= parts.keys():
                sample_kind: Literal["gsm", "gse_shared"] = (
                    "gse_shared" if sample_id == gse else "gsm"
                )
                supported_by_sample.setdefault(
                    sample_id,
                        DiscoveredSample(
                            gse=gse,
                            sample_id=sample_id,
                            input_type="triplet",
                            sample_kind=sample_kind,
                            supported=True,
                            note="matrix triplet",
                            source_name=parts["matrix"].name,
                            individual_id=individual_map.get((gse, sample_id), ""),
                            matrix_path=parts["matrix"],
                            barcodes_path=parts["barcodes"],
                            features_path=parts["features"],
                    ),
                )

        samples.extend(
            sorted(
                [*supported_by_sample.values(), *unsupported_by_sample.values()],
                key=lambda item: (item.gse, item.sample_id),
            )
        )

    return samples


__all__ = [
    "DEFAULT_DATA_ROOT",
    "DiscoveredSample",
    "discover_rna_samples",
    "resolve_data_root",
    "selected_rna_gses",
]
