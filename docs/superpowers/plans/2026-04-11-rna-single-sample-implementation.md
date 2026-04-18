# RNA Single-Sample Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the `only_rna` branch mainline for single-sample RNA QC, `L1/L2` CIMA annotation, and UMAP quality outputs on already-downloaded local RNA datasets.

**Architecture:** Keep analysis logic in a dedicated R single-sample entrypoint and keep discovery, dispatch, and status logic in `pipeline.py`. Use a small Python discovery layer that resolves the external data root, classifies RNA input formats, and only runs samples with a supported local input layout.

**Tech Stack:** Python 3, R, Seurat, scDblFinder, Matrix, ggplot2, gzip/tar utilities.

---

### Task 1: Add RNA Discovery Tests

**Files:**
- Create: `tests/test_rna_pipeline.py`
- Modify: `scripts/process/pipeline.py`

- [ ] **Step 1: Write failing tests for RNA data-root resolution and discovery behavior**

- [ ] **Step 2: Run the RNA pipeline tests to confirm they fail before implementation**

Run: `python3 -m unittest -v tests.test_rna_pipeline`

- [ ] **Step 3: Implement minimal Python helpers for RNA data-root resolution and RNA discovery**

- [ ] **Step 4: Re-run the RNA pipeline tests until they pass**

Run: `python3 -m unittest -v tests.test_rna_pipeline`

### Task 2: Add RNA Pipeline Commands

**Files:**
- Modify: `scripts/process/pipeline.py`

- [ ] **Step 1: Add `discover-rna`, `run-rna-sample`, `run-rna-gse`, and `rna-status` command coverage in tests or command-level assertions**

- [ ] **Step 2: Implement parser entries, status handling, and RNA command execution plumbing**

- [ ] **Step 3: Verify the updated Python module compiles**

Run: `python3 -m py_compile scripts/process/pipeline.py`

### Task 3: Add the Single-Sample RNA Workflow

**Files:**
- Create: `scripts/process/process_single_rna_sample.R`

- [ ] **Step 1: Implement RNA input loading for supported local formats**

- [ ] **Step 2: Implement QC, clustering, UMAP, and `L1/L2` CIMA annotation outputs**

- [ ] **Step 3: Write required output files and validation columns for audit-ready review**

- [ ] **Step 4: Verify the script parses successfully**

Run: `Rscript -e "parse(file='scripts/process/process_single_rna_sample.R')"`

### Task 4: Add the RNA Reference Builder

**Files:**
- Create: `scripts/process/build_cima_rna_reference_model.py`

- [ ] **Step 1: Implement an RNA compact reference builder compatible with the current external CIMA RNA asset layout**

- [ ] **Step 2: Verify the script compiles**

Run: `python3 -m py_compile scripts/process/build_cima_rna_reference_model.py`

### Task 5: Update Branch Documentation

**Files:**
- Modify: `AGENTS.md`

- [ ] **Step 1: Rewrite branch positioning from ATAC mainline to RNA single-sample mainline**

- [ ] **Step 2: Document RNA commands, outputs, data-root expectations, and `L1/L2` audit focus**

### Task 6: Verify on Representative Local RNA Samples

**Files:**
- No code changes required if implementation already passes

- [ ] **Step 1: Run Python-side tests**

Run: `python3 -m unittest -v tests.test_rna_pipeline`

- [ ] **Step 2: Run syntax checks**

Run: `python3 -m py_compile scripts/process/pipeline.py scripts/process/build_cima_rna_reference_model.py`

Run: `Rscript -e "parse(file='scripts/process/process_single_rna_sample.R')"`

- [ ] **Step 3: Run RNA discovery**

Run: `python3 scripts/process/pipeline.py discover-rna | head`

- [ ] **Step 4: Run one representative RNA sample into a local temporary output root**

Run: `python3 scripts/process/pipeline.py run-rna-sample --gse GSE167363 --sample-id GSM5102900 --output-root /tmp/ml2026_rna_smoke --force`

- [ ] **Step 5: Confirm required outputs exist for the smoke sample**

Check:
- `/tmp/ml2026_rna_smoke/GSE167363/GSM5102900/qc_summary.csv`
- `/tmp/ml2026_rna_smoke/GSE167363/GSM5102900/validation_result.csv`
- `/tmp/ml2026_rna_smoke/GSE167363/GSM5102900/umap_rna_cima_cell_type_l1.png`
- `/tmp/ml2026_rna_smoke/GSE167363/GSM5102900/umap_rna_cima_cell_type_l2.png`
