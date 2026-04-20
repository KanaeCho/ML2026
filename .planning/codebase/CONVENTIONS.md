# Coding Conventions

**Analysis Date:** Mon Apr 20 2026

## Naming Patterns

**Files:**
- Use `snake_case.py` for Python modules under `scripts/only_rna/`, `scripts/process/`, and `tests/`, for example `scripts/only_rna/read_inputs.py`, `scripts/only_rna/tuning_orchestrator.py`, and `tests/only_rna/test_outputs.py`.
- Use `test_*.py` for pytest modules under `tests/`, for example `tests/only_rna/test_cli.py` and `tests/test_rna_pipeline.py`.
- Use descriptive artifact filenames for outputs and validations, for example `scripts/only_rna/outputs.py` writes `metadata.csv`, `qc_summary.csv`, `validation_result.csv`, and `umap_rna_pbmcref_vs_cima_l1.png`.

**Functions:**
- Use `snake_case` for functions and helpers, including private helpers prefixed with `_`, for example `_resolve_selected_rna_gses()` in `scripts/only_rna/cli.py`, `_read_triplet()` in `scripts/only_rna/read_inputs.py`, and `_embedding_parameters()` in `scripts/only_rna/embedding.py`.
- Use verb-first names for pipeline stages and writers, for example `compute_qc_metrics()`, `apply_qc_filters()`, `run_embedding()`, `annotate_with_all_versions()`, and `write_sample_outputs()`.

**Variables:**
- Use `snake_case` for locals and attributes, including explicit domain names such as `pass_qc_mask`, `annotation_method_status`, `expected_paths`, and `selected_gses`; see `scripts/only_rna/outputs.py`, `scripts/only_rna/annotation.py`, and `scripts/process/pipeline.py`.
- Use ALL_CAPS for module constants and path defaults, for example `ROOT`, `DEFAULT_OUTPUT_ROOT`, `DEFAULT_CONFIG_PATH` in `scripts/only_rna/cli.py`, and `DEFAULT_DATA_ROOT`, `RNA_OUTPUT_DIR`, `FRAGMENT_RE` in `scripts/process/pipeline.py`.

**Types:**
- Use `PascalCase` for dataclasses and typed records, for example `RunConfig`, `QcThresholds`, `PlottingConfig`, and `DiscoveredSample` in `scripts/only_rna/models.py` and `scripts/only_rna/discovery.py`.
- Use explicit literal/value enums where helpful, for example `sample_kind: Literal["gsm", "gse_shared"]` in `scripts/only_rna/discovery.py`.

## Code Style

**Formatting:**
- No formatter configuration file was detected at repo root: no `.prettierrc`, `ruff.toml`, `pytest.ini`, `setup.cfg`, or `tox.ini` were found from the inspected manifests.
- Source formatting is consistent with Black-style Python: 4-space indentation, double quotes, trailing commas in multiline calls, and line wrapping around 88-ish chars; see `scripts/only_rna/config.py`, `scripts/only_rna/outputs.py`, and `tests/only_rna/test_tuning.py`.
- Future Python code should match the style already present in `scripts/only_rna/*.py` and `tests/**/*.py`.

**Linting:**
- No active lint runner configuration was detected in `pyproject.toml`, and no dedicated lint config file was found.
- Static typing is partially configured with `pyrightconfig.json`, but it only includes `scripts/process`, not `scripts/only_rna` or `tests`.

## Import Organization

**Order:**
1. `from __future__ import annotations`
2. Standard library imports
3. Third-party imports
4. Local package imports

**Observed examples:**
- `scripts/only_rna/config.py` imports `Path`, `dataclasses`, then `yaml`, then `.models`.
- `scripts/only_rna/outputs.py` imports stdlib, then `anndata`, `pandas`, `scipy`, then `.models` and `.plotting`.
- `tests/only_rna/test_processing.py` imports stdlib, then `anndata`/`numpy`/`pandas`/`scipy`, then `scripts.only_rna.*` modules.

**Path Aliases:**
- No import alias system was detected. Imports use package-relative imports inside `scripts/only_rna/` and absolute package imports like `from scripts.only_rna import cli` or `from scripts.process import pipeline` in tests.

## Type Hints and Dataclasses

**Patterns:**
- Use type annotations broadly on functions and dataclass fields, for example `load_run_config(path: Path) -> RunConfig` in `scripts/only_rna/config.py` and `read_sample_input(sample: DiscoveredSample) -> ad.AnnData` in `scripts/only_rna/read_inputs.py`.
- Use frozen dataclasses for immutable runtime config and sample descriptors, for example `@dataclass(frozen=True)` on all classes in `scripts/only_rna/models.py` and `scripts/only_rna/discovery.py`.
- Use modern unions and generic builtins such as `Path | None`, `list[str]`, and `dict[str, Path]`; see `scripts/only_rna/read_inputs.py` and `scripts/only_rna/outputs.py`.

## Configuration Conventions

**Primary runtime config:**
- Put default runtime config in YAML at `scripts/only_rna/default_config.yaml`.
- Load YAML into dataclasses through `scripts/only_rna/config.py`, specifically `load_run_config()` / `load_default_config()`.
- Organize YAML by top-level sections: `qc`, `plotting`, `annotation`, `azimuth`, and optionally `embedding` and `tuning`; see `scripts/only_rna/default_config.yaml` and tests in `tests/only_rna/test_config.py`.

**Override conventions:**
- Merge CLI/runtime overrides through `merge_cli_overrides()` in `scripts/only_rna/config.py`.
- Prefer nested override keys such as `embedding__n_neighbors` and `azimuth__k_weight`; tests in `tests/only_rna/test_config.py` verify this pattern.
- Preserve untouched config fields when overriding; this behavior is asserted in `tests/only_rna/test_processing.py` and `tests/only_rna/test_config.py`.

**Environment/path resolution:**
- Resolve data roots in code instead of hardcoding per command. Both `scripts/only_rna/discovery.py` and `scripts/process/pipeline.py` prefer `./data`, then `ML2026_DATA_ROOT`, then `/mnt/g/ML2026_data`.
- Keep output roots explicit and overridable through CLI arguments; see parser setup in `scripts/process/pipeline.py` and status routing in `scripts/only_rna/cli.py`.

## CLI and Output Conventions

**CLI design:**
- Define all command interfaces in `scripts/process/pipeline.py` using `argparse` subcommands.
- Route RNA subcommands from `scripts/process/pipeline.py` into `scripts.only_rna.cli`; tests in `tests/only_rna/test_cli.py` assert this contract.
- Use tab-separated plain-text output for discovery/status commands, for example `_print_rna_discovery()` and `_print_rna_status()` in `scripts/only_rna/cli.py`.

**Status tracking:**
- Write `run_status.json` in sample output directories through `_write_rna_status()` in `scripts/only_rna/cli.py`.
- Store structured fields such as `sample`, `input_type`, `sample_kind`, `started_at`, `finished_at`, `status`, and `outputs_complete`.

**Validation artifacts:**
- Every sample write should emit both content artifacts and audit artifacts. `scripts/only_rna/outputs.py` consistently writes `metadata.csv`, `metadata_qc.csv`, `qc_summary.csv`, `validation_result.csv`, `{sample_id}.h5ad`, matrix triplets, and UMAP PNGs.
- Remove stale legacy artifacts before writing current outputs; `scripts/only_rna/outputs.py` deletes entries from `REMOVED_SAMPLE_ROOT_ARTIFACTS`, and `tests/only_rna/test_outputs.py` asserts that behavior.

## DataFrame and AnnData Conventions

**obs column handling:**
- Initialize expected columns explicitly with target dtypes before filling them. `scripts/only_rna/annotation.py` creates string, float, and boolean annotation columns with `pd.Series(..., dtype=...)`.
- Keep failed/unprocessed cells in the full `adata`, then fill pass-QC-derived fields only for `pass_qc` rows; see `scripts/only_rna/embedding.py`, `scripts/only_rna/annotation.py`, and `scripts/only_rna/outputs.py`.
- Use `string`, `boolean`, and float dtypes intentionally rather than generic `object`; tests in `tests/only_rna/test_multimethod.py` assert dtype contracts.

**Copy semantics:**
- Functions generally avoid mutating caller-owned AnnData in place. `compute_qc_metrics()`, `apply_qc_filters()`, `run_doublet_detection()`, `run_embedding()`, `annotate_with_cima()`, and `write_sample_outputs()` all start from copied data or derived frames.

**H5AD sanitization:**
- Sanitize string/object columns before writing `.h5ad` using `_sanitize_dataframe_for_h5ad()` and `_prepare_h5ad_adata()` in `scripts/only_rna/outputs.py`.

## Error Handling

**Patterns:**
- Raise `FileNotFoundError` for missing external resources or sample lookup failures, for example in `scripts/only_rna/discovery.py`, `scripts/only_rna/cli.py`, and `scripts/only_rna/read_inputs.py`.
- Raise `ValueError` for invalid internal state or unsupported input shapes, for example missing triplet paths in `scripts/only_rna/read_inputs.py` and failed tuning candidate generation in `scripts/only_rna/tuning_orchestrator.py`.
- Convert optional tool failures into graceful fallbacks where possible. `scripts/only_rna/doublet.py` catches any scrublet error and normalizes placeholder columns; `scripts/only_rna/embedding.py` catches `ImportError` around Leiden/UMAP and writes deterministic fallback outputs.
- Record runtime failure details into `run_status.json` in `scripts/only_rna/cli.py` before re-raising.

## Logging

**Framework:** `print` + JSON status files

**Patterns:**
- No `logging` module usage was detected in the inspected RNA pipeline modules.
- Operational messaging uses `print()` to stdout/stderr in `scripts/only_rna/cli.py`.
- Durable execution state is stored in JSON (`run_status.json`, `selection_summary.json`, `selected_params.json`) rather than logs alone; see `scripts/only_rna/cli.py` and `scripts/only_rna/outputs.py`.

## Comments

**When to Comment:**
- Inline comments are sparse. Most code relies on descriptive names instead of heavy comments.
- Use section-divider comments for large conceptual areas, as in `scripts/only_rna/annotation.py` and `scripts/only_rna/plotting.py`.

**JSDoc/TSDoc:**
- Not applicable.

**Docstrings:**
- Use short functional docstrings selectively for public helpers, for example `load_run_config()` in `scripts/only_rna/config.py` and `annotate_with_azimuth()` in `scripts/only_rna/annotation.py`.
- Many tests and internal helpers omit docstrings.

## Function Design

**Size:**
- Small stage functions are preferred for core processing, for example `compute_qc_metrics()` in `scripts/only_rna/qc.py` and `run_doublet_detection()` in `scripts/only_rna/doublet.py`.
- Larger orchestration/writer functions exist where contract assembly is centralized, notably `write_sample_outputs()` in `scripts/only_rna/outputs.py` and CLI routing in `scripts/only_rna/cli.py`.

**Parameters:**
- Pass config objects explicitly rather than reading globals inside processing stages, for example `run_embedding(adata, config)` and `write_sample_outputs(..., config=config)`.
- Pass filesystem paths as `Path` objects throughout the codebase.

**Return Values:**
- Return transformed AnnData from stage functions instead of `None`, for example `compute_qc_metrics()`, `apply_qc_filters()`, `run_embedding()`, and `annotate_with_all_versions()`.
- Return output directories from writer/orchestrator functions, for example `write_sample_outputs()` and `write_tuning_selection_artifacts()` in `scripts/only_rna/outputs.py`.

## Module Design

**Exports:**
- Public surface areas are declared with `__all__` in small modules such as `scripts/only_rna/qc.py`, `scripts/only_rna/doublet.py`, `scripts/only_rna/embedding.py`, `scripts/only_rna/read_inputs.py`, `scripts/only_rna/discovery.py`, and `scripts/only_rna/outputs.py`.

**Barrel Files:**
- Not used. `scripts/only_rna/__init__.py` exists, but modules are imported directly rather than through a barrel API.

## Quality Risks Visible From Conventions

- Type checking scope is incomplete because `pyrightconfig.json` only includes `scripts/process`, while the active RNA implementation lives under `scripts/only_rna/`.
- No enforced formatter/linter configuration was detected, so current style consistency depends on contributor discipline and tests.
- Logging/observability is mostly `print()` plus JSON files, which makes structured runtime diagnostics limited outside written artifacts.
- Some behavior is contract-driven by tests rather than centralized config, especially artifact names and validation checks in `tests/only_rna/test_outputs.py`; future changes should update both code and tests together.

---

*Convention analysis: Mon Apr 20 2026*
