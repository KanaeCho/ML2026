# Codebase Structure

**Analysis Date:** Mon Apr 20 2026

## Directory Layout

```text
ML2026/
├── AGENTS.md                 # Current branch operating contract and workflow notes
├── docs/                     # Specs and implementation plans for RNA work
├── scripts/                  # All runnable processing code and notebooks
│   ├── only_rna/             # Python-first RNA subsystem
│   └── process/              # CLI router, legacy R flows, TEA-seq/scATAC utilities
├── tests/                    # Python tests for pipeline and only_rna modules
├── data/                     # Optional in-workspace data root
├── output/                   # Run artifacts; current workflow treats this as generated output
├── pyproject.toml            # Python package/dependency definition
├── uv.lock                   # Locked dependency graph for uv
└── .planning/codebase/       # Generated codebase maps for downstream agents
```

## Directory Purposes

**`scripts/only_rna/`:**
- Purpose: Hold the active single-sample RNA processing subsystem.
- Contains: discovery, input loading, QC, doublet detection, embedding, annotation, plotting, output writing, tuning, and CLI handlers.
- Key files: `scripts/only_rna/cli.py`, `scripts/only_rna/discovery.py`, `scripts/only_rna/read_inputs.py`, `scripts/only_rna/qc.py`, `scripts/only_rna/embedding.py`, `scripts/only_rna/annotation.py`, `scripts/only_rna/outputs.py`, `scripts/only_rna/tuning_orchestrator.py`.

**`scripts/process/`:**
- Purpose: Hold the repository-wide CLI entrypoint plus non-mainline or legacy processing scripts.
- Contains: `argparse` router, R wrappers, reference builders, download utilities, TEA-seq/scATAC scripts, experimental comparison flows.
- Key files: `scripts/process/pipeline.py`, `scripts/process/process_single_sample.R`, `scripts/process/process_single_rna_sample.R`, `scripts/process/download_from_datasets.py`, `scripts/process/organize_tea_seq_outputs.py`, `scripts/process/compare_gse192391_annotation_methods.py`.

**`tests/only_rna/`:**
- Purpose: Verify the Python RNA subsystem at module and workflow boundaries.
- Contains: tests for config loading, discovery, processing stages, CLI routing, outputs, tuning, Azimuth integration, and multimethod annotation.
- Key files: `tests/only_rna/test_cli.py`, `tests/only_rna/test_processing.py`, `tests/only_rna/test_outputs.py`, `tests/only_rna/test_tuning.py`, `tests/only_rna/test_discovery.py`.

**`tests/`:**
- Purpose: Hold higher-level or cross-module tests outside the `only_rna` package folder.
- Contains: pipeline routing tests and comparison-workflow tests.
- Key files: `tests/test_rna_pipeline.py`, `tests/test_gse192391_compare.py`.

**`docs/superpowers/specs/` and `docs/superpowers/plans/`:**
- Purpose: Preserve design docs and implementation plans for the RNA mainline and tuning workflows.
- Contains: markdown specs and matching implementation plans.
- Key files: `docs/superpowers/specs/2026-04-13-only-rna-python-mainline-design.md`, `docs/superpowers/specs/2026-04-14-only-rna-umap-visual-tuning-design.md`.

**`data/`:**
- Purpose: Optional local data root if present.
- Contains: runtime input data and references when the workspace-local layout is used.
- Key files: Not committed as code; resolution logic is defined in `scripts/process/pipeline.py` and `scripts/only_rna/discovery.py`.

**`output/`:**
- Purpose: Generated run artifacts and status files.
- Contains: legacy outputs under `output/{GSE}/{GSM}/` and RNA outputs under `output/rna/{GSE}/{sample_id}/`.
- Key files: Generated files such as `run_status.json`, `metadata.csv`, `qc_summary.csv`, `{sample_id}.h5ad`.

**`.planning/codebase/`:**
- Purpose: Store generated reference documents for planning/execution agents.
- Contains: architecture, structure, stack, testing, conventions, and concern maps.
- Key files: `ARCHITECTURE.md`, `STRUCTURE.md`.

## Key File Locations

**Entry Points:**
- `scripts/process/pipeline.py`: Main CLI router for discovery, runs, download, status, RNA execution, and TEA-seq audit.
- `scripts/process/build_cima_rna_reference_model.py`: Rebuild RNA CIMA reference assets.
- `scripts/process/build_cima_reference_model.py`: Rebuild non-RNA CIMA reference assets.
- `scripts/process/compare_gse192391_annotation_methods.py`: Standalone annotation-comparison workflow.

**Configuration:**
- `pyproject.toml`: Project metadata and Python dependencies.
- `uv.lock`: Resolved dependency lockfile.
- `pyrightconfig.json`: Type-checker configuration.
- `scripts/only_rna/default_config.yaml`: Default QC, plotting, annotation, and Azimuth settings.
- `AGENTS.md`: Current branch behavior contract and output expectations.

**Core Logic:**
- `scripts/only_rna/discovery.py`: Data-root resolution and supported sample discovery.
- `scripts/only_rna/read_inputs.py`: Input readers for triplet, archive, and `.h5` formats.
- `scripts/only_rna/qc.py`: QC metric computation and QC pass/fail flags.
- `scripts/only_rna/doublet.py`: Scrublet-based or normalized doublet handling.
- `scripts/only_rna/embedding.py`: PCA/neighbors/leiden/UMAP pipeline for `pass_qc` cells.
- `scripts/only_rna/annotation.py`: CIMA and Azimuth-oriented annotation orchestration.
- `scripts/only_rna/outputs.py`: Output export and validation report writing.
- `scripts/only_rna/tuning_orchestrator.py`: Bounded tuning execution and winner selection.

**Testing:**
- `tests/test_rna_pipeline.py`: Pipeline-level discovery and data-root tests.
- `tests/only_rna/test_cli.py`: CLI routing and output-contract tests.
- `tests/only_rna/test_processing.py`: Unit tests for read/QC/doublet/embedding/annotation helpers.
- `tests/only_rna/test_tuning.py`: Tuning preset and artifact tests.
- `tests/test_gse192391_compare.py`: Comparison-workflow tests.

## Naming Conventions

**Files:**
- Snake_case Python modules for subsystem stages: `scripts/only_rna/read_inputs.py`, `scripts/only_rna/tuning_metrics.py`.
- Verb-oriented process scripts for one-off workflows: `scripts/process/export_h5ad_obs.py`, `scripts/process/render_integration_umap_panels.py`.
- R entrypoints use snake_case with domain-specific names: `scripts/process/process_single_sample.R`, `scripts/process/process_single_rna_sample.R`.
- Design and plan docs use dated kebab-case names: `docs/superpowers/specs/2026-04-13-only-rna-python-mainline-design.md`.

**Directories:**
- Subsystems are grouped by responsibility rather than by package depth: `scripts/only_rna/`, `scripts/process/`, `tests/only_rna/`.
- Generated or cache directories use tool defaults: `__pycache__/`, `.pytest_cache/`, `.ruff_cache/`.

## Where to Add New Code

**New RNA Feature:**
- Primary code: add stage logic in `scripts/only_rna/` next to the closest concern, such as `scripts/only_rna/annotation.py` for annotation changes or `scripts/only_rna/outputs.py` for artifact changes.
- Tests: mirror the module under `tests/only_rna/`, for example `tests/only_rna/test_outputs.py` or `tests/only_rna/test_processing.py`.

**New CLI Command Related to Existing Workflow:**
- Implementation routing: `scripts/process/pipeline.py` if the command should be globally visible.
- RNA-specific handler: `scripts/only_rna/cli.py` if the command belongs to the RNA command family.

**New Component/Module:**
- RNA implementation: `scripts/only_rna/{new_module}.py`.
- Process/utility implementation: `scripts/process/{new_script}.py` or `.R` if it remains a standalone workflow.

**Utilities:**
- Shared RNA helpers: keep them in `scripts/only_rna/` rather than creating a generic utils folder.
- One-off analysis helpers: keep them in `scripts/process/` close to the invoking script.

**New Documentation:**
- Behavior/spec documents: `docs/superpowers/specs/`.
- Implementation plans: `docs/superpowers/plans/`.
- Agent-facing codebase maps: `.planning/codebase/`.

## Special Directories

**`output/`:**
- Purpose: Runtime artifacts for sample and dataset processing.
- Generated: Yes.
- Committed: Not treated as source-of-truth code; current branch notes in `AGENTS.md` treat it as generated output.

**`data/`:**
- Purpose: Preferred in-workspace data root when present.
- Generated: External/runtime managed.
- Committed: Not relied on as committed source code.

**`scripts/only_rna/__pycache__/` and `tests/**/__pycache__/`:**
- Purpose: Python bytecode caches.
- Generated: Yes.
- Committed: Present in the working tree but not source code.

**`.agents/skills/`:**
- Purpose: Repository-local skill library used by agent tooling.
- Generated: No.
- Committed: Yes.

**`.planning/codebase/`:**
- Purpose: Generated architecture/reference maps consumed by other GSD commands.
- Generated: Yes.
- Committed: Intended to be committed as planning artifacts.

---

*Structure analysis: Mon Apr 20 2026*
