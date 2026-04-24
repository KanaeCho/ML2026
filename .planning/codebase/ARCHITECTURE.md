# Architecture

**Analysis Date:** Mon Apr 20 2026

## Pattern Overview

**Overall:** Hybrid orchestration with a Python-first RNA subsystem and legacy process scripts.

**Key Characteristics:**
- `scripts/process/pipeline.py` is the top-level CLI router for both legacy scATAC/TEA-seq flows and the current RNA mainline.
- `scripts/only_rna/` is a self-contained sample-processing subsystem built around `AnnData` and pure-Python stage functions.
- Data discovery, execution, and artifact validation are coupled through fixed output contracts under `output/` and `output/rna/`.

## Layers

**CLI Orchestration Layer:**
- Purpose: Parse commands, resolve data roots, route execution, and manage run-status files.
- Location: `scripts/process/pipeline.py`, `scripts/only_rna/cli.py`
- Contains: `argparse` command registration, discovery/status commands, sample routing, dry-run/force handling.
- Depends on: `scripts.only_rna.*`, legacy R scripts in `scripts/process/`, filesystem state under `data/` and `output/`.
- Used by: Direct CLI execution via `python scripts/process/pipeline.py ...`.

**RNA Discovery and Input Layer:**
- Purpose: Resolve external data roots and turn raw GEO file layouts into runnable sample objects.
- Location: `scripts/only_rna/discovery.py`, `scripts/only_rna/read_inputs.py`
- Contains: `DiscoveredSample`, GSE selection from `reference/datasets.xlsx`, triplet/`.h5`/archive readers.
- Depends on: `scripts/process/download_from_datasets.py`, raw files under `data_root/raw/{GSE}/`.
- Used by: `scripts/only_rna/cli.py`, `scripts/only_rna/tuning_orchestrator.py`.

**RNA Processing Layer:**
- Purpose: Apply the single-sample analysis pipeline to one `AnnData` object.
- Location: `scripts/only_rna/qc.py`, `scripts/only_rna/doublet.py`, `scripts/only_rna/embedding.py`, `scripts/only_rna/annotation.py`, `scripts/only_rna/azimuth.py`
- Contains: QC metrics, QC filters, Scrublet integration, embedding/clustering/UMAP, CIMA projection, Azimuth annotation.
- Depends on: `RunConfig` from `scripts/only_rna/models.py`, reference assets under `reference/cima/`, optional R/Azimuth runtime.
- Used by: `scripts/only_rna/cli.py`, `scripts/only_rna/tuning_orchestrator.py`.

**RNA Output and Audit Layer:**
- Purpose: Materialize CSV/H5AD/matrix/PNG artifacts and encode completion checks.
- Location: `scripts/only_rna/outputs.py`, `scripts/only_rna/plotting.py`
- Contains: metadata export, matrix triplet export, validation files, QC plot, dual-annotation UMAP, tuning summaries.
- Depends on: processed `AnnData`, plotting config, output-root conventions.
- Used by: `scripts/only_rna/cli.py`, `scripts/only_rna/tuning_orchestrator.py`.

**Tuning Layer:**
- Purpose: Enumerate bounded preset combinations, execute each candidate, and select a winner.
- Location: `scripts/only_rna/tuning_presets.py`, `scripts/only_rna/tuning_metrics.py`, `scripts/only_rna/tuning_orchestrator.py`
- Contains: preset families, scoring functions, candidate execution loop, selection artifact writing.
- Depends on: the same RNA processing/output pipeline used by mainline sample runs.
- Used by: `scripts/only_rna/cli.py` through `tune-rna-sample` and `tune-rna-gse`.

**Legacy and Auxiliary Processing Layer:**
- Purpose: Preserve non-mainline scATAC, TEA-seq, download, comparison, and reference-building workflows.
- Location: `scripts/process/process_single_sample.R`, `scripts/process/process_single_rna_sample.R`, `scripts/process/organize_tea_seq_outputs.py`, `scripts/process/compare_gse192391_annotation_methods.py`, `scripts/process/build_cima_reference_model.py`, `scripts/process/build_cima_rna_reference_model.py`
- Contains: legacy R entrypoints, dataset download helpers, TEA-seq audits, comparison experiments, reference rebuild scripts.
- Depends on: raw/reference/output directories and external R/Python toolchains.
- Used by: direct script invocation and selected routes from `scripts/process/pipeline.py`.

## Data Flow

**RNA Single-Sample Mainline:**

1. `scripts/process/pipeline.py` parses `run-rna-sample` or `run-rna-gse` and forwards to `scripts/only_rna/cli.py`.
2. `scripts/only_rna/discovery.py` resolves the data root and discovers supported samples from `data_root/raw/{GSE}/` after filtering GSEs from `reference/datasets.xlsx`.
3. `scripts/only_rna/read_inputs.py` converts a discovered triplet, 10x `.h5`, or `matrix.tar.gz` into `AnnData` and normalizes sample metadata columns.
4. `scripts/only_rna/qc.py`, `scripts/only_rna/doublet.py`, and `scripts/only_rna/embedding.py` compute QC fields, mark `pass_qc`, and write `cluster`/`umap_1`/`umap_2` into `adata.obs`.
5. `scripts/only_rna/annotation.py` orchestrates CIMA/Azimuth-style annotation and records backend status in `adata.uns['annotation_method_status']`.
6. `scripts/only_rna/outputs.py` writes sample artifacts under `output/rna/{GSE}/{sample_id}/` and emits `validation_result.csv` against the expected output contract.
7. `scripts/only_rna/cli.py` writes `run_status.json` to the same sample directory before and after execution.

**RNA Bounded Tuning Flow:**

1. `scripts/process/pipeline.py` forwards `tune-rna-sample` and `tune-rna-gse` to `scripts/only_rna/cli.py`.
2. `scripts/only_rna/tuning_orchestrator.py` enumerates candidate IDs from preset families in `scripts/only_rna/tuning_presets.py`.
3. Each candidate reruns the full RNA mainline with a candidate-specific `RunConfig` built through `scripts/only_rna/config.py`.
4. `scripts/only_rna/tuning_metrics.py` scores QC retention, annotation confidence, and embedding quality.
5. `scripts/only_rna/outputs.py` writes candidate outputs below `output/rna/{GSE}/{sample_id}/tuning/{candidate_id}/{GSE}/{sample_id}/` and summary files under `.../tuning/`.

**Legacy scATAC / TEA-seq Flow:**

1. `scripts/process/pipeline.py` handles `discover`, `run-sample`, `run-gse`, `download`, `status`, and `tea-seq-audit` directly.
2. Legacy sample runs shell out to `scripts/process/process_single_sample.R` or `scripts/process/process_single_rna_sample.R` via `run_command(...)` in `scripts/process/pipeline.py`.
3. Dataset-level TEA-seq cleanup/audit runs through `scripts/process/organize_tea_seq_outputs.py` and writes audit artifacts under `output/{GSE}/qc_audit/`.

**State Management:**
- Runtime state is file-based rather than service-based.
- Sample state lives in `run_status.json` under each output directory.
- Analysis state lives inside `AnnData.obs`, `AnnData.var`, `AnnData.obsm`, and `AnnData.uns` during processing.

## Key Abstractions

**Discovered Sample:**
- Purpose: Represent one runnable RNA input target with enough metadata to choose a reader and output path.
- Examples: `scripts/only_rna/discovery.py`, `scripts/only_rna/cli.py`, `tests/only_rna/test_discovery.py`
- Pattern: Immutable dataclass carrying source paths and a `sample_kind` of `gsm` or `gse_shared`.

**RunConfig:**
- Purpose: Carry all sample-processing configuration across stages.
- Examples: `scripts/only_rna/models.py`, `scripts/only_rna/config.py`, `scripts/only_rna/default_config.yaml`
- Pattern: Nested immutable dataclasses loaded from YAML and selectively overridden for tuning.

**Sample Output Contract:**
- Purpose: Define completion criteria for status commands and validation reports.
- Examples: `scripts/only_rna/cli.py`, `scripts/process/pipeline.py`, `scripts/only_rna/outputs.py`
- Pattern: Fixed expected file lists checked with `exists()` and mirrored in validation CSV rows.

**CandidateSpec / CandidateEvaluation:**
- Purpose: Encode bounded tuning search space and scored results.
- Examples: `scripts/only_rna/tuning_orchestrator.py`
- Pattern: Dataclasses with deterministic candidate IDs of the form `qc__azimuth__embedding`.

## Entry Points

**Primary CLI:**
- Location: `scripts/process/pipeline.py`
- Triggers: Manual CLI invocation.
- Responsibilities: Build subcommands, discover samples, route RNA commands into `scripts.only_rna.cli`, shell out to legacy scripts, and maintain run status.

**RNA Command Handlers:**
- Location: `scripts/only_rna/cli.py`
- Triggers: `discover-rna`, `run-rna-sample`, `run-rna-gse`, `tune-rna-sample`, `tune-rna-gse`, `rna-status`.
- Responsibilities: Resolve selected GSEs, forbid explicit `gse_shared` sample runs, execute the RNA processing chain, and emit RNA-specific status rows.

**Reference Builders:**
- Location: `scripts/process/build_cima_rna_reference_model.py`, `scripts/process/build_cima_reference_model.py`
- Triggers: Direct script execution.
- Responsibilities: Rebuild compact CIMA reference assets for runtime annotation.

**Comparison Workflow:**
- Location: `scripts/process/compare_gse192391_annotation_methods.py`
- Triggers: Direct script execution and tests in `tests/test_gse192391_compare.py`.
- Responsibilities: Evaluate multiple annotation backends outside the mainline contract.

## Error Handling

**Strategy:** Fail fast on missing prerequisites, but degrade gracefully inside optional analysis stages.

**Patterns:**
- Discovery and routing raise `FileNotFoundError` or return nonzero status for missing samples and unsupported explicit targets in `scripts/only_rna/cli.py` and `scripts/process/pipeline.py`.
- Doublet detection and embedding fall back to normalized defaults when optional libraries fail in `scripts/only_rna/doublet.py` and `scripts/only_rna/embedding.py`.
- Azimuth wraps external R execution and records `status/detail` instead of crashing downstream consumers in `scripts/only_rna/azimuth.py`.
- Tuning converts candidate exceptions into `reason_code=evaluation_error:...` rows so the search can continue in `scripts/only_rna/tuning_orchestrator.py`.

## Cross-Cutting Concerns

**Logging:** File-backed command logs are written under sample `logs/` directories by `run_command(...)` in `scripts/process/pipeline.py`; Python RNA execution relies primarily on status JSON plus raised exceptions.

**Validation:** Output completeness is validated through expected-path checks in `scripts/process/pipeline.py`, `scripts/only_rna/cli.py`, and `scripts/only_rna/outputs.py`.

**Authentication:** Not detected in the architecture itself. External annotation via Azimuth depends on local R packages, and data selection reads local `reference/datasets.xlsx` rather than a remote service.

---

*Architecture analysis: Mon Apr 20 2026*
