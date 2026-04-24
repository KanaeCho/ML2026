# Technology Stack

**Analysis Date:** Mon Apr 20 2026

## Languages

**Primary:**
- Python 3.11+ - main single-sample RNA pipeline, discovery, tuning, download tooling, and validation in `pyproject.toml`, `scripts/only_rna/*.py`, and `scripts/process/*.py`
- R - Azimuth/Seurat integration and legacy ATAC / TEA-seq workflows in `scripts/only_rna/azimuth.py`, `scripts/process/process_single_rna_sample.R`, and `scripts/process/process_single_sample.R`

**Secondary:**
- YAML - runtime config in `scripts/only_rna/default_config.yaml`
- JSON - status/config artifacts in `opencode.json`, `pyrightconfig.json`, and runtime `run_status.json` paths referenced by `scripts/only_rna/cli.py` and `scripts/process/pipeline.py`
- Markdown - operational documentation in `AGENTS.md`

## Runtime

**Environment:**
- CPython >=3.11 required by `pyproject.toml` and `uv.lock`
- Local R runtime invoked via `Rscript` from `scripts/only_rna/azimuth.py` and `scripts/process/pipeline.py`

**Package Manager:**
- `uv` - project dependency management and tool execution in `pyproject.toml` and `opencode.json`
- Lockfile: present in `uv.lock`

## Frameworks

**Core:**
- AnnData `>=0.8.0` / locked `0.12.6` - annotated matrix container used across `scripts/only_rna/read_inputs.py`, `scripts/only_rna/annotation.py`, and `scripts/only_rna/outputs.py`
- Scanpy `>=1.9.0` - RNA preprocessing, embedding input/output handling, and 10x H5 reading in `pyproject.toml` and `scripts/only_rna/read_inputs.py`
- Seurat + Azimuth (R) - reference mapping with `pbmcref` in `scripts/only_rna/azimuth.py` and `scripts/process/process_single_rna_sample.R`
- Signac + Seurat (R) - retained scATAC workflow in `scripts/process/process_single_sample.R`

**Testing:**
- `pytest>=9.0.3` - dev test runner declared in `pyproject.toml`
- `unittest` - used directly in `tests/test_rna_pipeline.py`

**Build/Dev:**
- BasedPyright `>=1.39.2` - static analysis dependency in `pyproject.toml`, configured by `pyrightconfig.json` and `opencode.json`
- OpenCode LSP wiring - `basedpyright-langserver` launched through `uv run` in `opencode.json`

## Key Dependencies

**Critical:**
- `scanpy` - core single-cell RNA workflow execution referenced in `scripts/only_rna/read_inputs.py` and described in `AGENTS.md`
- `anndata` - `.h5ad` read/write and in-memory data model in `scripts/only_rna/read_inputs.py` and `scripts/only_rna/outputs.py`
- `pandas`, `numpy`, `scipy`, `h5py` - matrix IO, tabular metadata, sparse operations, and HDF5 support across `scripts/only_rna/*.py`
- `igraph` and `leidenalg` - clustering dependencies declared in `pyproject.toml` for embedding/clustering workflows referenced in `AGENTS.md`
- `pyyaml` - runtime config loader in `scripts/only_rna/config.py`

**Infrastructure:**
- `openpyxl` plus custom XLSX parsing logic - dataset selection driven by `datasets.xlsx` in `scripts/process/download_from_datasets.py`
- `geoparse` / `datasets` - downloader-related Python environment dependencies declared in `pyproject.toml`
- `torch`, `torchvision`, `accelerate`, `sentencepiece`, `protobuf` - installed in the Python environment via `pyproject.toml`; present as platform/runtime dependencies even though they are not part of the observed RNA mainline entry flow
- `matplotlib` and `umap-learn` - plotting and UMAP rendering support in `scripts/only_rna/plotting.py`
- `pyside6` - installed dependency in `pyproject.toml`; no direct usage detected in the inspected mainline scripts

## Configuration

**Environment:**
- Data root resolution checks `./data`, then `ML2026_DATA_ROOT`, then `/mnt/g/ML2026_data` in `scripts/only_rna/discovery.py` and `scripts/process/pipeline.py`
- RNA defaults are configured in `scripts/only_rna/default_config.yaml` and loaded by `scripts/only_rna/config.py`
- Python version is pinned by `.python-version` and `pyrightconfig.json`

**Build:**
- Dependency and resolver config: `pyproject.toml`, `uv.lock`
- Type-checking / editor config: `pyrightconfig.json`, `opencode.json`
- Runtime documentation contract: `AGENTS.md`

## CLI / Runtime Entrypoints

**Primary CLI:**
- `scripts/process/pipeline.py` - top-level orchestration CLI for `discover-rna`, `run-rna-sample`, `run-rna-gse`, `tune-rna-sample`, `tune-rna-gse`, `rna-status`, and downloader/audit commands
- `scripts/only_rna/cli.py` - Python-first RNA execution path used by `scripts/process/pipeline.py`

**Auxiliary runtimes:**
- `scripts/only_rna/azimuth.py` - shells out to `Rscript -e` for Azimuth mapping
- `scripts/process/process_single_rna_sample.R` - R implementation for RNA sample processing retained in the process layer
- `scripts/process/process_single_sample.R` - retained ATAC runtime

## Platform Requirements

**Development:**
- Python 3.11 virtual environment expected at `.venv` per `pyrightconfig.json`
- `uv` available for dependency installation and LSP execution per `opencode.json`
- R with packages including `Seurat`, `Azimuth`, `Matrix`, `ggplot2`, `scDblFinder`, `SingleCellExperiment`, `Signac`, and `RhpcBLASctl` required by `scripts/only_rna/azimuth.py`, `scripts/process/process_single_rna_sample.R`, and `scripts/process/process_single_sample.R`
- Optional `aria2c` binary used by downloader CLI in `scripts/process/pipeline.py`

**Production:**
- Local workstation / external-mounted-data execution model, not a packaged service; outputs are written under `output/rna/...` and data is expected under external data roots described in `AGENTS.md`, `scripts/only_rna/discovery.py`, and `scripts/process/pipeline.py`

---

*Stack analysis: Mon Apr 20 2026*
