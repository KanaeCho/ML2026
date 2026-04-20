# External Integrations

**Analysis Date:** Mon Apr 20 2026

## APIs & External Services

**Public data sources:**
- NCBI GEO supplementary file hosting - dataset download source for raw sample files
  - SDK/Client: Python standard library `urllib.request` in `scripts/process/download_from_datasets.py`
  - Auth: Not detected
  - Endpoints assembled from `https://ftp.ncbi.nlm.nih.gov/geo/samples/...` and `https://ftp.ncbi.nlm.nih.gov/geo/series/...` in `scripts/process/download_from_datasets.py`

**Reference annotation tooling:**
- Azimuth `pbmcref` reference mapping - RNA cell-type annotation backend invoked during mainline/sample processing
  - SDK/Client: R packages `Azimuth` and `Seurat` invoked from `scripts/only_rna/azimuth.py`
  - Auth: Not detected
  - Invocation: `RunAzimuth(query = query, reference = 'pbmcref', ...)` in `scripts/only_rna/azimuth.py`

**External command-line tools:**
- `Rscript` - required bridge into R-based Azimuth and legacy workflows from `scripts/only_rna/azimuth.py` and `scripts/process/pipeline.py`
  - SDK/Client: subprocess execution via Python `subprocess.run` / `subprocess.Popen`
  - Auth: Not applicable
- `aria2c` - optional high-throughput downloader used by the `download` command in `scripts/process/pipeline.py`
  - SDK/Client: CLI binary path passed through `--aria2c`
  - Auth: Not applicable

## Data Storage

**Databases:**
- Not detected

**File Storage:**
- Local/external filesystem storage only
  - Data root resolved from `./data`, `ML2026_DATA_ROOT`, or `/mnt/g/ML2026_data` in `scripts/only_rna/discovery.py` and `scripts/process/pipeline.py`
  - RNA outputs written under `output/rna/{GSE}/{sample_id}/` by `scripts/only_rna/outputs.py` and `scripts/only_rna/cli.py`
  - Output symlink/external storage assumption documented in `AGENTS.md`

**Caching:**
- In-process listing/size caches for downloader only in `scripts/process/download_from_datasets.py`
- No standalone Redis/Memcached-style cache detected

## Authentication & Identity

**Auth Provider:**
- None detected
  - Implementation: public GEO download access plus local filesystem execution

## Monitoring & Observability

**Error Tracking:**
- None detected

**Logs:**
- Per-sample and per-download log files are written to filesystem by `scripts/process/pipeline.py`
- Run state is tracked with JSON status files such as `run_status.json` from `scripts/only_rna/cli.py` and `scripts/process/pipeline.py`

## CI/CD & Deployment

**Hosting:**
- Not detected as a deployed web service; observed execution model is local CLI/batch processing on a workstation with mounted external storage per `AGENTS.md`

**CI Pipeline:**
- None detected

## Environment Configuration

**Required env vars:**
- `ML2026_DATA_ROOT` - optional external data root override in `scripts/only_rna/discovery.py`, `scripts/process/pipeline.py`, and `scripts/process/process_single_rna_sample.R`

**Secrets location:**
- Not detected

## Reference Assets and Runtime Assumptions

**Reference files:**
- `data_root/reference/datasets.xlsx` - authoritative sample selection source read by `scripts/process/download_from_datasets.py`, `scripts/only_rna/discovery.py`, and `scripts/process/pipeline.py`
- `data_root/reference/cima/*` - local CIMA reference assets loaded by `scripts/only_rna/annotation.py` and described in `AGENTS.md`

**Runtime assumptions:**
- Raw data is organized under `data_root/raw/{GSE}/` and discovered by filename conventions in `scripts/only_rna/discovery.py`
- `GSE226039` is filtered to PBMC files only in `scripts/only_rna/discovery.py` and `scripts/process/pipeline.py`
- Unsupported shared gene-count CSV inputs are explicitly excluded in `scripts/only_rna/discovery.py` and `tests/test_rna_pipeline.py`

## Webhooks & Callbacks

**Incoming:**
- None detected

**Outgoing:**
- GEO HTTP/FTP listing and file fetch requests from `scripts/process/download_from_datasets.py`
- Local subprocess calls to `Rscript` and optional `aria2c` from `scripts/only_rna/azimuth.py` and `scripts/process/pipeline.py`

---

*Integration audit: Mon Apr 20 2026*
