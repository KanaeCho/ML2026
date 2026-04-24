# Codebase Concerns

**Analysis Date:** Mon Apr 20 2026

## Tech Debt

**Duplicated RNA discovery and data-root resolution logic:**
- Issue: RNA sample discovery and data-root fallback logic are implemented twice, once in `scripts/only_rna/discovery.py` and again in `scripts/process/pipeline.py`, with parallel `resolve_data_root(...)`, selected-GSE discovery, sample models, and local-layout parsing.
- Files: `scripts/only_rna/discovery.py`, `scripts/process/pipeline.py`, `tests/test_rna_pipeline.py`, `tests/only_rna/test_discovery.py`
- Impact: Behavior can drift between the Python-first RNA mainline and the older process entrypoint, especially for supported/unsupported sample classification, data-root fallback, and special-case handling.
- Fix approach: Keep one canonical discovery/data-root module in `scripts/only_rna/discovery.py` and have `scripts/process/pipeline.py` delegate to it instead of re-implementing the same rules.

**Documented contract drift between branch docs and sample-root outputs:**
- Issue: branch docs describe sample-root outputs such as `umap_rna_clusters.png` and `umap_rna_azimuth.png`, while `scripts/only_rna/outputs.py` explicitly deletes several legacy sample-root artifacts and currently writes `umap_rna_pbmcref_vs_cima_l1.png` as the main annotation figure.
- Files: `AGENTS.md`, `scripts/only_rna/outputs.py`, `scripts/only_rna/cli.py`, `scripts/process/pipeline.py`
- Impact: Operators and downstream automation can rely on stale output expectations, making run validation and manual review inconsistent with current code behavior.
- Fix approach: Align `AGENTS.md`, status checks, and any orchestration docs with the artifacts actually emitted by `scripts/only_rna/outputs.py`.

**Tuning documentation and implementation are partially misaligned:**
- Issue: branch docs state candidate execution includes `annotate_with_all_versions(..., methods=['cima', 'azimuth'])`, but `scripts/only_rna/tuning_orchestrator.py` currently evaluates candidates with `methods=["azimuth"]` only.
- Files: `AGENTS.md`, `scripts/only_rna/tuning_orchestrator.py`
- Impact: Candidate scoring and audit outputs do not fully match the documented bounded-tuning contract, which can mislead future changes and review of winner selection.
- Fix approach: Either restore CIMA execution in tuning candidates or update branch documentation to state that tuning is currently Azimuth-only for annotation scoring.

**Legacy/placeholder annotation paths remain in active modules:**
- Issue: `scripts/only_rna/annotation.py` still contains `annotate_with_azimuth(...)`, `annotate_with_cell_typist(...)`, `annotate_with_singler(...)`, and `annotate_with_scanvi(...)` best-effort or placeholder helpers, while the main orchestrator uses `run_azimuth_annotation(...)` plus raw-label alignment instead.
- Files: `scripts/only_rna/annotation.py`, `scripts/only_rna/azimuth.py`, `tests/only_rna/test_multimethod.py`
- Impact: The module surface suggests broader supported annotation backends than the branch contract actually guarantees, increasing maintenance burden and ambiguity about which path is authoritative.
- Fix approach: Remove or isolate deprecated/best-effort helpers, and keep one clearly documented production annotation path per backend.

## Known Bugs

**Output status can report success while external annotations degraded silently:**
- Symptoms: sample runs can complete and write outputs even when Azimuth or other optional annotation helpers fail internally, because several code paths convert exceptions into empty/NA outputs rather than hard failures.
- Files: `scripts/only_rna/azimuth.py`, `scripts/only_rna/annotation.py`, `scripts/only_rna/doublet.py`, `scripts/only_rna/embedding.py`
- Trigger: Missing R/Azimuth dependencies, missing `igraph`, scrublet errors, or optional backend import failures.
- Workaround: Inspect `qc_summary.csv`, `validation_result.csv`, and `adata.uns['annotation_method_status']` rather than relying only on exit status or `run_status.json`.

**Shared GSE-level samples are discoverable but not runnable via explicit sample command:**
- Symptoms: discovery returns `sample_kind="gse_shared"` entries, but `run-rna-sample` and `tune-rna-sample` explicitly refuse them.
- Files: `scripts/only_rna/discovery.py`, `scripts/only_rna/cli.py`, `AGENTS.md`
- Trigger: Running `run-rna-sample --gse <GSE> --sample-id <same GSE>` for a shared triplet dataset.
- Workaround: Use `run-rna-gse` or `tune-rna-gse` so the shared sample is included implicitly.

## Security Considerations

**Archive extraction trusts tar contents:**
- Risk: `tarfile.extractall(...)` is used without member path validation when reading archived matrix inputs.
- Files: `scripts/only_rna/read_inputs.py`
- Current mitigation: Archives are limited to discovered local sample bundles, but extraction is still unconditional.
- Recommendations: Validate tar members before extraction or implement safe extraction that rejects absolute paths and parent-directory traversal.

**Subprocess-generated R code embeds file paths directly:**
- Risk: Azimuth execution constructs inline R code strings with interpolated paths and parameters and shells out through `Rscript -e`.
- Files: `scripts/only_rna/azimuth.py`, `scripts/process/compare_gse192391_annotation_methods.py`
- Current mitigation: Paths come from local temp files and configured references, not user-entered shell fragments.
- Recommendations: Prefer writing a temporary `.R` script or passing structured arguments to reduce quoting fragility and improve auditability.

## Performance Bottlenecks

**Bounded tuning reruns the full pipeline per candidate with no caching:**
- Problem: each tuning candidate rereads the sample, recomputes QC, reruns doublet detection, re-embeds, re-annotates, and rewrites full sample outputs.
- Files: `scripts/only_rna/tuning_orchestrator.py`, `scripts/only_rna/read_inputs.py`, `scripts/only_rna/outputs.py`
- Cause: candidate evaluation is fully sequential and does not reuse intermediate artifacts between preset combinations.
- Improvement path: Cache immutable stages per sample, reuse pre-QC or post-QC intermediates where valid, and separate heavy computation from final artifact rendering.

**Tuning overview generation rereads candidate `.h5ad` files after evaluation:**
- Problem: the winner-selection artifact writer loops over candidate output directories and reloads stored `.h5ad` files to render overview panels.
- Files: `scripts/only_rna/outputs.py`
- Cause: `write_tuning_selection_artifacts(...)` reconstructs plotting inputs from disk instead of reusing in-memory candidate results.
- Improvement path: Persist only minimal plotting tables or pass already loaded candidate summaries into the artifact writer.

## Fragile Areas

**Embedding relies on silent fallback modes for missing dependencies and small samples:**
- Files: `scripts/only_rna/embedding.py`, `tests/only_rna/test_processing.py`, `AGENTS.md`
- Why fragile: `leiden` and `umap` downgrade to single-cluster or coordinate-sequence fallbacks on `ImportError`, and small samples bypass the normal graph/UMAP path entirely.
- Safe modification: Change embedding logic together with tests covering `igraph` absence, small-`n_obs` behavior, and expected output columns.
- Test coverage: There are unit tests for fallback branches, but no end-to-end assertion that fallback outputs remain analytically acceptable for real branch datasets.

**Azimuth integration depends on an external R runtime with limited preflight checks:**
- Files: `scripts/only_rna/azimuth.py`, `scripts/only_rna/annotation.py`, `tests/only_rna/test_azimuth.py`
- Why fragile: The mainline assumes `Rscript`, `Azimuth`, `Seurat`, and the `pbmcref` reference are operational at runtime, but dependency validation happens only when annotation is attempted.
- Safe modification: Add an explicit environment preflight command and keep annotation status semantics stable in `adata.uns['annotation_method_status']`.
- Test coverage: Tests mainly monkeypatch `_run_azimuth_r(...)`; they do not exercise a real R/Azimuth environment.

**Branch-specific dataset rules are hard-coded in discovery:**
- Files: `scripts/only_rna/discovery.py`, `scripts/process/pipeline.py`, `AGENTS.md`
- Why fragile: `GSE226039` PBMC-only inclusion and shared gene-count rejection are embedded as literal filename rules, not as configurable dataset policy.
- Safe modification: Move per-dataset exceptions into declarative configuration or a reference manifest before adding more dataset-specific branches.
- Test coverage: Unit tests cover the current exceptions, but each new exceptional dataset would require more code branching.

## Scaling Limits

**Single-sample orchestration does not scale gracefully to larger candidate sets or broad batch execution:**
- Current capacity: tuning is bounded by `config.tuning.max_candidates` and defaults to a small fixed candidate family in `scripts/only_rna/tuning_orchestrator.py`.
- Limit: runtime and storage scale roughly linearly with candidate count because each candidate writes a full nested sample output tree under `output/rna/{GSE}/{sample_id}/tuning/{candidate_id}/{GSE}/{sample_id}/`.
- Scaling path: Separate candidate metrics from full artifact emission, add reuse of intermediates, and reserve full nested outputs for top candidates or debug mode.

**Mainline remains constrained to single-sample processing:**
- Current capacity: branch docs explicitly scope first-stage acceptance to single-sample scRNA-seq only.
- Limit: there is no abstraction in `scripts/only_rna/cli.py` or `scripts/only_rna/tuning_orchestrator.py` for cross-sample normalization, shared embeddings, or integrated annotation consistency.
- Scaling path: Introduce a separate multi-sample integration layer rather than extending the current sample-local pipeline ad hoc.

## Dependencies at Risk

**External annotation runtime stack is operationally heavy:**
- Risk: the main annotation path depends on Python Scanpy/AnnData plus R `Azimuth`/`Seurat`, and branch docs also mention fallback behavior when `igraph` is absent.
- Files: `pyproject.toml`, `scripts/only_rna/azimuth.py`, `scripts/only_rna/embedding.py`, `AGENTS.md`
- Impact: environment drift can change clustering or annotation behavior without code changes, especially across machines where R and Python stacks are provisioned differently.
- Migration plan: add explicit environment validation scripts and version-pinned setup docs for both Python and R dependencies before broadening production use.

## Missing Critical Features

**No explicit preflight for runtime prerequisites:**
- Problem: there is no dedicated command that verifies Python package availability, `Rscript`, Azimuth/Seurat installation, CIMA reference assets, and output/data-root accessibility before a run starts.
- Blocks: reliable operator onboarding and fast diagnosis of environment-specific failures.

**No declarative dataset exception registry:**
- Problem: branch-specific exceptions such as `GSE226039` PBMC-only handling and unsupported shared gene-count matrices are encoded directly in Python conditionals.
- Blocks: safe growth of dataset support without accumulating more special-case branches in discovery code.

## Test Coverage Gaps

**Real external dependency integration is mostly untested:**
- What's not tested: actual `Rscript`/Azimuth execution, real `pbmcref` availability, real scrublet behavior under production environments, and end-to-end runs against the external data-root layout.
- Files: `tests/only_rna/test_azimuth.py`, `tests/only_rna/test_processing.py`, `tests/only_rna/test_tuning.py`, `scripts/only_rna/azimuth.py`
- Risk: environment regressions can pass unit tests while failing on the real branch runtime.
- Priority: High

**Output-contract drift is not asserted against `AGENTS.md`:**
- What's not tested: that branch documentation and orchestration expectations match the exact files emitted by `scripts/only_rna/outputs.py`.
- Files: `AGENTS.md`, `scripts/only_rna/outputs.py`, `scripts/only_rna/cli.py`
- Risk: future agents can rely on stale docs and add more contract divergence unnoticed.
- Priority: Medium

**Repository hygiene issues are present in the tracked workspace:**
- What's not tested: prevention of committed transient artifacts such as `__pycache__` trees, `.ruff_cache/`, `.pytest_cache/`, and `.DS_Store` files that are visible in the working tree.
- Files: `scripts/only_rna/__pycache__/`, `scripts/process/__pycache__/`, `tests/__pycache__/`, `.ruff_cache/`, `.pytest_cache/README.md`, `.DS_Store`, `scripts/.DS_Store`
- Risk: noisy diffs and branch clutter make maintenance harder and can hide meaningful source changes.
- Priority: Low

---

*Concerns audit: Mon Apr 20 2026*
