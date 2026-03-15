#!/usr/bin/env python3

from __future__ import annotations

import argparse
import html.parser
import json
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from zipfile import ZipFile
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[2]
DATASETS_XLSX = ROOT / "data" / "reference" / "datasets.xlsx"
RAW_DIR = ROOT / "data" / "raw"
ASSAY_COL = "测序数据(scATAC/scRNA)"
FORMAT_COL = "数据格式"
GSM_COL = "样本名(GSM*)"
GSE_COL = "数据集(GSE*)"
NS_MAIN = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
NS_REL = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
DEFAULT_NETWORK_TIMEOUT = 10
USER_AGENT = "ML2026-downloader/1.0"


@dataclass(frozen=True)
class DatasetRow:
    gsm: str
    gse: str
    assay: str
    data_format: str
    raw: dict[str, str]


@dataclass(frozen=True)
class RemoteFile:
    gsm: str
    gse: str
    url: str
    filename: str
    size_bytes: int | None

    @property
    def target_dir(self) -> Path:
        return RAW_DIR / self.gse

    @property
    def target_path(self) -> Path:
        return self.target_dir / self.filename

    @property
    def control_path(self) -> Path:
        return self.target_dir / f"{self.filename}.aria2"


@dataclass(frozen=True)
class LocalFileStatus:
    state: str
    local_size: int | None


@dataclass(frozen=True)
class LocalSampleInventory:
    gsm: str
    gse: str
    files: tuple[Path, ...]
    has_fragment: bool
    has_barcode: bool
    has_tbi: bool


class HrefParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.hrefs: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        for key, value in attrs:
            if key.lower() == "href" and value:
                self.hrefs.append(value)


def normalize_text(value: str) -> str:
    return re.sub(r"\s+", " ", str(value or "")).strip()


def bucket_accession(accession: str) -> str:
    match = re.fullmatch(r"([A-Z]+)(\d+)", accession)
    if not match:
        raise ValueError(f"Unsupported GEO accession: {accession}")
    prefix, digits = match.groups()
    if len(digits) <= 3:
        raise ValueError(f"GEO accession too short to bucket: {accession}")
    return f"{prefix}{digits[:-3]}nnn"


def shared_strings(zip_file: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in zip_file.namelist():
        return []
    root = ET.fromstring(zip_file.read("xl/sharedStrings.xml"))
    ns = {"a": NS_MAIN}
    strings: list[str] = []
    for si in root.findall("a:si", ns):
        text = "".join(node.text or "" for node in si.iter(f"{{{NS_MAIN}}}t"))
        strings.append(text)
    return strings


def parse_xlsx_rows(path: Path) -> list[dict[str, str]]:
    with ZipFile(path) as zip_file:
        workbook = ET.fromstring(zip_file.read("xl/workbook.xml"))
        ns = {"a": NS_MAIN, "r": NS_REL}
        sheets = workbook.find("a:sheets", ns)
        if sheets is None or len(sheets) == 0:
            raise ValueError(f"No sheets found in {path}")

        first_sheet = sheets[0]
        rel_id = first_sheet.attrib[f"{{{NS_REL}}}id"]

        rels_root = ET.fromstring(zip_file.read("xl/_rels/workbook.xml.rels"))
        relmap = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels_root}
        sheet_target = relmap[rel_id]
        if not sheet_target.startswith("worksheets/"):
            sheet_target = f"worksheets/{Path(sheet_target).name}"
        sheet_xml = f"xl/{sheet_target}"

        hyperlinks: dict[str, str] = {}
        rel_path = f"xl/worksheets/_rels/{Path(sheet_target).name}.rels"
        if rel_path in zip_file.namelist():
            rels = ET.fromstring(zip_file.read(rel_path))
            hyperlink_map = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
        else:
            hyperlink_map = {}

        sheet_root = ET.fromstring(zip_file.read(sheet_xml))
        shared = shared_strings(zip_file)
        sheet_hyperlinks = sheet_root.find("a:hyperlinks", ns)
        if sheet_hyperlinks is not None:
            for hyperlink in sheet_hyperlinks.findall("a:hyperlink", ns):
                ref = hyperlink.attrib["ref"]
                rel = hyperlink.attrib.get(f"{{{NS_REL}}}id")
                if rel and rel in hyperlink_map:
                    hyperlinks[ref] = hyperlink_map[rel]

        rows: list[dict[str, str]] = []
        sheet_data = sheet_root.find("a:sheetData", ns)
        if sheet_data is None:
            return rows

        for row in sheet_data.findall("a:row", ns):
            values: dict[str, str] = {}
            for cell in row.findall("a:c", ns):
                ref = cell.attrib["r"]
                col = re.sub(r"\d+", "", ref)
                cell_type = cell.attrib.get("t")
                raw_value = cell.find("a:v", ns)
                value = "" if raw_value is None else raw_value.text or ""
                if cell_type == "s" and value:
                    value = shared[int(value)]
                values[col] = value
                if ref in hyperlinks:
                    values[f"{col}_hyperlink"] = hyperlinks[ref]
            rows.append(values)

    if not rows:
        return []

    header_map = {col: name for col, name in rows[0].items() if not col.endswith("_hyperlink")}
    normalized: list[dict[str, str]] = []
    for row in rows[1:]:
        item: dict[str, str] = {}
        for col, name in header_map.items():
            item[name] = row.get(col, "")
            link_key = f"{col}_hyperlink"
            if link_key in row:
                item[f"{name}_hyperlink"] = row[link_key]
        normalized.append(item)
    return normalized


def load_filtered_rows(
    xlsx_path: Path,
    assay: str,
    data_format: str,
    gse: str = "",
    gsm: str = "",
) -> list[DatasetRow]:
    rows = parse_xlsx_rows(xlsx_path)
    filtered: list[DatasetRow] = []
    for row in rows:
        assay_value = normalize_text(row.get(ASSAY_COL, ""))
        format_value = normalize_text(row.get(FORMAT_COL, ""))
        gsm_value = normalize_text(row.get(GSM_COL, ""))
        gse_value = normalize_text(row.get(GSE_COL, ""))
        gse_match = re.search(r"GSE\d+", gse_value)

        if assay_value != assay:
            continue
        if data_format.lower() not in format_value.lower():
            continue
        if not re.fullmatch(r"GSM\d+", gsm_value):
            continue
        if gse_match is None:
            continue

        item = DatasetRow(
            gsm=gsm_value,
            gse=gse_match.group(0),
            assay=assay_value,
            data_format=format_value,
            raw=row,
        )

        if gse and item.gse != gse:
            continue
        if gsm and item.gsm != gsm:
            continue

        filtered.append(
            DatasetRow(
                gsm=item.gsm,
                gse=item.gse,
                assay=item.assay,
                data_format=item.data_format,
                raw=item.raw,
            )
        )
    return filtered


def listing_url_for_sample(gsm: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{bucket_accession(gsm)}/{gsm}/suppl/"


def listing_url_for_series(gse: str) -> str:
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket_accession(gse)}/{gse}/suppl/"


def urlopen_with_headers(url: str, timeout: int, method: str = "GET"):
    request = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def fetch_listing(url: str, timeout: int = DEFAULT_NETWORK_TIMEOUT) -> list[str]:
    with urlopen_with_headers(url, timeout=timeout) as response:
        content = response.read().decode("utf-8", errors="replace")
    parser = HrefParser()
    parser.feed(content)
    items: list[str] = []
    for href in parser.hrefs:
        if href in {"../", "./"} or href.startswith("?") or href.startswith("/"):
            continue
        items.append(href)
    return sorted(set(items))


def file_matches_sample(gsm: str, filename: str, file_kinds: set[str]) -> bool:
    lower = filename.lower()
    if not lower.startswith(gsm.lower()):
        return False
    if "fragment" in file_kinds and "fragment" in lower and not lower.endswith(".tbi"):
        return True
    if "fragment" in file_kinds and lower.endswith(".tbi") and "fragment" in lower:
        return True
    if "barcode" in file_kinds and "barcode" in lower:
        return True
    if "singlecell" in file_kinds and "singlecell" in lower and lower.endswith(".csv.gz"):
        return True
    if "summary" in file_kinds and "summary" in lower and lower.endswith(".csv.gz"):
        return True
    return False


def remote_file_size(url: str, timeout: int = DEFAULT_NETWORK_TIMEOUT) -> int | None:
    try:
        with urlopen_with_headers(url, timeout=timeout, method="HEAD") as response:
            length = response.headers.get("Content-Length")
            return None if length is None else int(length)
    except Exception:
        try:
            with urlopen_with_headers(url, timeout=timeout) as response:
                length = response.headers.get("Content-Length")
                return None if length is None else int(length)
        except Exception:
            return None


def scan_local_inventory(rows: list[DatasetRow]) -> list[LocalSampleInventory]:
    inventory: list[LocalSampleInventory] = []
    for row in rows:
        gse_dir = RAW_DIR / row.gse
        matches = tuple(sorted(path for path in gse_dir.glob(f"{row.gsm}_*") if path.is_file()))
        names = [path.name.lower() for path in matches]
        inventory.append(
            LocalSampleInventory(
                gsm=row.gsm,
                gse=row.gse,
                files=matches,
                has_fragment=any("fragment" in name and not name.endswith(".tbi") for name in names),
                has_barcode=any("barcode" in name for name in names),
                has_tbi=any(name.endswith(".tbi") and "fragment" in name for name in names),
            )
        )
    return inventory


def resolve_remote_files(
    rows: list[DatasetRow],
    timeout: int,
    file_kinds: set[str],
) -> tuple[list[RemoteFile], list[str]]:
    resolved: list[RemoteFile] = []
    problems: list[str] = []
    listing_cache: dict[str, tuple[list[str], str | None]] = {}
    size_cache: dict[str, int | None] = {}

    def get_listing(url: str, label: str) -> tuple[list[str], str | None]:
        if url in listing_cache:
            return listing_cache[url]

        try:
            result = (fetch_listing(url, timeout=timeout), None)
        except urllib.error.HTTPError as exc:
            result = ([], f"{label} listing unavailable ({exc.code})")
        except Exception as exc:
            result = ([], f"{label} listing unavailable ({exc})")

        listing_cache[url] = result
        return result

    total_rows = len(rows)
    for index, row in enumerate(rows, start=1):
        print(f"[resolve] {index}/{total_rows} {row.gse}/{row.gsm}")
        candidates: list[tuple[str, list[str]]] = []
        local_errors: list[str] = []
        sample_url = listing_url_for_sample(row.gsm)
        sample_listing, sample_error = get_listing(sample_url, "sample")
        candidates.append((sample_url, sample_listing))
        if sample_error:
            local_errors.append(sample_error)

        series_url = listing_url_for_series(row.gse)
        series_listing, series_error = get_listing(series_url, "series")
        candidates.append((series_url, series_listing))
        if series_error:
            local_errors.append(series_error)

        picked: list[RemoteFile] = []
        for base_url, files in candidates:
            matched = [name for name in files if file_matches_sample(row.gsm, name, file_kinds)]
            if not matched:
                continue
            for name in matched:
                file_url = f"{base_url}{name}"
                if file_url not in size_cache:
                    size_cache[file_url] = remote_file_size(file_url, timeout=timeout)
                picked.append(
                    RemoteFile(
                        gsm=row.gsm,
                        gse=row.gse,
                        url=file_url,
                        filename=name,
                        size_bytes=size_cache[file_url],
                    )
                )
            break

        if not picked:
            detail = "; ".join(local_errors) if local_errors else "no matching files"
            requested = ",".join(sorted(file_kinds))
            problems.append(f"{row.gsm}: no requested files ({requested}) found in GEO supplementary listings ({detail})")
            continue

        resolved.extend(picked)

    deduped: dict[tuple[str, str], RemoteFile] = {}
    for item in resolved:
        deduped[(item.gsm, item.filename)] = item
    return list(deduped.values()), problems


def sizeof_fmt(num: int | None) -> str:
    if num is None:
        return "unknown"
    size = float(num)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024 or unit == "TB":
            return f"{size:.2f} {unit}"
        size /= 1024
    return f"{size:.2f} TB"


def local_file_status(item: RemoteFile) -> LocalFileStatus:
    path = item.target_path
    control_exists = item.control_path.exists()

    if control_exists:
        local_size = path.stat().st_size if path.exists() else None
        return LocalFileStatus(state="partial", local_size=local_size)

    if not path.exists():
        return LocalFileStatus(state="missing", local_size=None)

    local_size = path.stat().st_size
    if item.size_bytes is None:
        return LocalFileStatus(state="complete", local_size=local_size)

    if local_size == item.size_bytes:
        return LocalFileStatus(state="complete", local_size=local_size)

    return LocalFileStatus(state="partial", local_size=local_size)


def summarize(
    rows: list[DatasetRow],
    remote_files: list[RemoteFile],
    local_inventory: list[LocalSampleInventory],
) -> dict[str, object]:
    gse_counts: dict[str, int] = {}
    for row in rows:
        gse_counts[row.gse] = gse_counts.get(row.gse, 0) + 1

    local_sample_count = sum(1 for item in local_inventory if item.files)
    local_file_count = sum(len(item.files) for item in local_inventory)
    local_fragment_count = sum(1 for item in local_inventory if item.has_fragment)
    local_barcode_count = sum(1 for item in local_inventory if item.has_barcode)
    local_tbi_count = sum(1 for item in local_inventory if item.has_tbi)

    total_known_bytes = sum(item.size_bytes for item in remote_files if item.size_bytes is not None)
    missing_sizes = sum(1 for item in remote_files if item.size_bytes is None)

    to_download = []
    existing = 0
    partial = 0
    remaining_bytes = 0
    for item in remote_files:
        status = local_file_status(item)
        if status.state == "complete":
            existing += 1
            continue
        if status.state == "partial":
            partial += 1
            to_download.append(item)
            if item.size_bytes is not None and status.local_size is not None:
                remaining_bytes += max(item.size_bytes - status.local_size, 0)
            continue
        to_download.append(item)
        if item.size_bytes is not None:
            remaining_bytes += item.size_bytes

    free_bytes = shutil.disk_usage(RAW_DIR).free if RAW_DIR.exists() else shutil.disk_usage(ROOT).free

    return {
        "filtered_rows": len(rows),
        "unique_gse": len(gse_counts),
        "gse_counts": gse_counts,
        "local_existing_samples": local_sample_count,
        "local_existing_files": local_file_count,
        "local_fragment_samples": local_fragment_count,
        "local_barcode_samples": local_barcode_count,
        "local_tbi_samples": local_tbi_count,
        "resolved_files": len(remote_files),
        "existing_files": existing,
        "partial_files": partial,
        "files_to_download": len(to_download),
        "known_total_bytes": total_known_bytes,
        "remaining_bytes": remaining_bytes,
        "missing_size_files": missing_sizes,
        "free_bytes": free_bytes,
    }


def build_aria2_manifest(remote_files: list[RemoteFile]) -> Path:
    tmp = tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False, prefix="geo-download-", suffix=".txt")
    manifest_path = Path(tmp.name)
    with tmp:
        for item in remote_files:
            if local_file_status(item).state == "complete":
                continue
            item.target_dir.mkdir(parents=True, exist_ok=True)
            tmp.write(item.url + "\n")
            tmp.write(f" dir={item.target_dir}\n")
            tmp.write(f" out={item.filename}\n\n")
    return manifest_path


def write_links_manifest(remote_files: list[RemoteFile], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        handle.write("gse\tgsm\tfilename\tsize_bytes\turl\n")
        for item in sorted(remote_files, key=lambda x: (x.gse, x.gsm, x.filename)):
            size_text = "" if item.size_bytes is None else str(item.size_bytes)
            handle.write(f"{item.gse}\t{item.gsm}\t{item.filename}\t{size_text}\t{item.url}\n")


def print_links(remote_files: list[RemoteFile]) -> None:
    for item in sorted(remote_files, key=lambda x: (x.gse, x.gsm, x.filename)):
        print(f"{item.gse}\t{item.gsm}\t{item.filename}\t{item.url}")


def run_aria2c(aria2c_bin: str, manifest_path: Path, max_concurrent_downloads: int, split: int) -> int:
    command = [
        aria2c_bin,
        "--continue=true",
        "--auto-file-renaming=false",
        "--allow-overwrite=false",
        f"--max-concurrent-downloads={max_concurrent_downloads}",
        f"--split={split}",
        "--min-split-size=16M",
        f"--input-file={manifest_path}",
    ]
    print("[download] " + " ".join(command))
    return subprocess.call(command, cwd=ROOT)


def print_summary(summary: dict[str, object]) -> None:
    print(f"Filtered rows: {summary['filtered_rows']}")
    print(f"Unique GSE: {summary['unique_gse']}")
    for gse, count in sorted(summary["gse_counts"].items()):
        print(f"  {gse}: {count} GSM rows")
    print(f"Local existing GSM samples: {summary['local_existing_samples']}")
    print(f"Local existing files: {summary['local_existing_files']}")
    print(f"Local fragment samples: {summary['local_fragment_samples']}")
    print(f"Local barcode samples: {summary['local_barcode_samples']}")
    print(f"Local fragment index samples: {summary['local_tbi_samples']}")
    print(f"Resolved remote files: {summary['resolved_files']}")
    print(f"Resolved files already present locally: {summary['existing_files']}")
    print(f"Resolved files present but incomplete: {summary['partial_files']}")
    print(f"Files to download: {summary['files_to_download']}")
    print(f"Known total size: {sizeof_fmt(summary['known_total_bytes'])}")
    print(f"Remaining download size: {sizeof_fmt(summary['remaining_bytes'])}")
    print(f"Files with unknown size: {summary['missing_size_files']}")
    print(f"Free disk space: {sizeof_fmt(summary['free_bytes'])}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Filter datasets.xlsx and download GEO fragment/barcode files with aria2c")
    parser.add_argument("--xlsx-path", default=str(DATASETS_XLSX))
    parser.add_argument("--assay", default="scATAC")
    parser.add_argument("--data-format", default="fragment")
    parser.add_argument("--gse", default="")
    parser.add_argument("--gsm", default="")
    parser.add_argument("--file-kinds", default="fragment,barcode")
    parser.add_argument("--aria2c", default="aria2c")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-network-resolve", action="store_true")
    parser.add_argument("--max-concurrent-downloads", type=int, default=4)
    parser.add_argument("--split", type=int, default=8)
    parser.add_argument("--network-timeout", type=int, default=DEFAULT_NETWORK_TIMEOUT)
    parser.add_argument("--links-out", default="")
    parser.add_argument("--print-links", action="store_true")
    args = parser.parse_args()

    xlsx_path = Path(args.xlsx_path)
    if not xlsx_path.exists():
        raise FileNotFoundError(f"datasets workbook not found: {xlsx_path}")

    file_kinds = {
        token.strip().lower()
        for token in args.file_kinds.split(",")
        if token.strip()
    }
    if not file_kinds:
        raise ValueError("At least one file kind must be requested via --file-kinds")

    rows = load_filtered_rows(
        xlsx_path,
        assay=args.assay,
        data_format=args.data_format,
        gse=args.gse,
        gsm=args.gsm,
    )
    print(f"Workbook: {xlsx_path}")
    print(f"Filter: assay={args.assay!r}, data_format contains {args.data_format!r}")
    if args.gse:
        print(f"GSE filter: {args.gse}")
    if args.gsm:
        print(f"GSM filter: {args.gsm}")
    print(f"Requested GEO file kinds: {', '.join(sorted(file_kinds))}")

    if not rows:
        print("No matching rows found")
        return 0

    local_inventory = scan_local_inventory(rows)

    if args.skip_network_resolve:
        summary = {
            "filtered_rows": len(rows),
            "unique_gse": len({row.gse for row in rows}),
            "gse_counts": {gse: sum(1 for row in rows if row.gse == gse) for gse in sorted({row.gse for row in rows})},
            "local_existing_samples": sum(1 for item in local_inventory if item.files),
            "local_existing_files": sum(len(item.files) for item in local_inventory),
            "local_fragment_samples": sum(1 for item in local_inventory if item.has_fragment),
            "local_barcode_samples": sum(1 for item in local_inventory if item.has_barcode),
            "local_tbi_samples": sum(1 for item in local_inventory if item.has_tbi),
            "resolved_files": 0,
            "existing_files": 0,
            "files_to_download": 0,
            "known_total_bytes": 0,
            "remaining_bytes": 0,
            "missing_size_files": 0,
            "free_bytes": shutil.disk_usage(RAW_DIR).free if RAW_DIR.exists() else shutil.disk_usage(ROOT).free,
        }
        print_summary(summary)
        print("Skipped remote resolution")
        return 0

    print(f"Network timeout per request: {args.network_timeout}s")
    remote_files, problems = resolve_remote_files(rows, timeout=args.network_timeout, file_kinds=file_kinds)
    summary = summarize(rows, remote_files, local_inventory)
    print_summary(summary)

    if remote_files:
        if args.print_links:
            print("Direct GEO file links:")
            print_links(remote_files)
        if args.links_out:
            links_out = Path(args.links_out)
            if not links_out.is_absolute():
                links_out = ROOT / links_out
            write_links_manifest(remote_files, links_out)
            print(f"Saved direct GEO links: {links_out}")

    if problems:
        print("Problems:")
        for problem in problems:
            print(f"  - {problem}")

    if not remote_files:
        return 1

    manifest_preview = [
        {
            "gsm": item.gsm,
            "gse": item.gse,
            "filename": item.filename,
            "url": item.url,
            "size_bytes": item.size_bytes,
        }
        for item in remote_files
    ]
    print(json.dumps(manifest_preview[:10], ensure_ascii=False, indent=2))

    if args.dry_run:
        return 0 if not problems else 1

    aria2c_bin = shutil.which(args.aria2c)
    if aria2c_bin is None:
        raise FileNotFoundError("aria2c not found in PATH")

    manifest_path = build_aria2_manifest(remote_files)
    try:
        return run_aria2c(
            aria2c_bin=aria2c_bin,
            manifest_path=manifest_path,
            max_concurrent_downloads=args.max_concurrent_downloads,
            split=args.split,
        )
    finally:
        try:
            manifest_path.unlink()
        except FileNotFoundError:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
