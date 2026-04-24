# Testing Patterns

**Analysis Date:** Mon Apr 20 2026

## Test Framework

**Runner:**
- `pytest` is the primary runner, declared in `pyproject.toml` under `[dependency-groups].dev`.
- Config file: Not detected. No `pytest.ini`, `tox.ini`, or `[tool.pytest.ini_options]` section was found in `pyproject.toml`.

**Assertion Library:**
- Native `assert` statements are the main assertion style across `tests/only_rna/*.py`.
- `pytest.approx` is used for float comparisons in files such as `tests/only_rna/test_outputs.py` and `tests/only_rna/test_tuning.py`.
- `unittest.TestCase` is still used in `tests/test_rna_pipeline.py` for older pipeline tests.

**Run Commands:**
```bash
pytest                      # Run all tests
pytest tests/only_rna       # Run RNA-focused suite
pytest tests/only_rna/test_outputs.py -q   # Run a focused file
```

## Test File Organization

**Location:**
- Tests live under `tests/`.
- Active RNA unit/integration-style tests are grouped in `tests/only_rna/`.
- Older parser/discovery coverage also exists at top level in `tests/test_rna_pipeline.py` and `tests/test_gse192391_compare.py`.

**Naming:**
- Use `test_*.py` filenames, for example `tests/only_rna/test_config.py`, `tests/only_rna/test_processing.py`, and `tests/only_rna/test_tuning.py`.
- Use `test_*` function names that describe behavior, often with full contract language, for example `test_write_sample_outputs_validation_contract_excludes_removed_cima_and_comparison_artifacts` in `tests/only_rna/test_outputs.py`.

**Structure:**
```text
tests/
├── only_rna/
│   ├── test_azimuth.py
│   ├── test_cli.py
│   ├── test_config.py
│   ├── test_discovery.py
│   ├── test_embedding.py
│   ├── test_multimethod.py
│   ├── test_outputs.py
│   ├── test_processing.py
│   └── test_tuning.py
├── test_gse192391_compare.py
└── test_rna_pipeline.py
```

## Test Structure

**Suite Organization:**
```python
def test_run_embedding_uses_explicit_embedding_config(monkeypatch):
    adata = ad.AnnData(...)
    config = _make_run_config(...)

    monkeypatch.setattr("scripts.only_rna.embedding.sc.tl.umap", fake_umap)

    out = run_embedding(adata, config)

    assert captured["umap"]["min_dist"] == 0.12
    assert out.obs["cluster"].notna().all()
```

**Patterns:**
- Each file typically defines local factory helpers such as `_make_run_config()`, `_make_output_adata()`, `_make_sample_adata()`, or `_touch_triplet()`; see `tests/only_rna/test_outputs.py`, `tests/only_rna/test_multimethod.py`, and `tests/only_rna/test_discovery.py`.
- Tests are behavior-focused and map closely to production contracts: config loading, CLI routing, artifact presence, dtype guarantees, and fallback semantics.
- Temporary filesystem isolation uses `tmp_path` or `tempfile.TemporaryDirectory()` rather than checked-in fixtures; see `tests/only_rna/test_cli.py` and `tests/test_rna_pipeline.py`.

## Mocking

**Framework:**
- `pytest` `monkeypatch` is the standard mocking mechanism.

**Patterns:**
```python
monkeypatch.setattr(
    "scripts.only_rna.tuning_orchestrator.run_embedding",
    fake_run_embedding,
)
monkeypatch.setattr(
    matplotlib.axes.Axes,
    "legend",
    capture_legend,
)
```

**Observed usage:**
- Patch heavy external dependencies to keep tests local and deterministic, such as `scanpy` functions in `tests/only_rna/test_embedding.py`, `scripts.only_rna.azimuth._run_azimuth_r` in `tests/only_rna/test_azimuth.py`, and full pipeline stage functions in `tests/only_rna/test_tuning.py`.
- Patch rendering primitives to assert plot ergonomics without pixel-perfect snapshots, such as `matplotlib.axes.Axes.legend`, `matplotlib.axes.Axes.scatter`, `matplotlib.axes.Axes.text`, and `plt.subplots` in `tests/only_rna/test_outputs.py`.
- Patch CLI routing targets to verify command dispatch in `tests/only_rna/test_cli.py`.

**What to Mock:**
- External tools and expensive libraries: `scanpy`, Azimuth/R bridges, filesystem-heavy writers, and orchestration stage boundaries.
- Plotting internals when verifying layout/legend/text behavior.

**What NOT to Mock:**
- Core data contracts when cheap to run locally. Many tests instantiate real `AnnData` objects and assert real obs/var transformations in `tests/only_rna/test_processing.py`, `tests/only_rna/test_multimethod.py`, and `tests/only_rna/test_outputs.py`.

## Fixtures and Factories

**Test Data:**
```python
def _make_output_adata() -> ad.AnnData:
    return ad.AnnData(
        X=sparse.csr_matrix(np.array([[1.0, 0.0, 3.0], ...])),
        obs=pd.DataFrame({...}, index=["cell-1", "cell-2", "cell-3"]),
        var=pd.DataFrame({...}, index=["GeneA", "GeneB", "GeneC"]),
    )
```

**Location:**
- Factories are defined inline inside each test module instead of a shared `conftest.py`. Examples include:
  - `_make_run_config()` in `tests/only_rna/test_outputs.py`, `tests/only_rna/test_embedding.py`, and `tests/test_gse192391_compare.py`
  - `_write_triplet()` and `_make_sample()` in `tests/only_rna/test_processing.py`
  - `_touch_triplet()` in `tests/only_rna/test_discovery.py`

## Coverage

**Requirements:**
- No explicit coverage threshold or coverage tool configuration was detected.
- No `coverage`, `pytest-cov`, or CI coverage gate was found in the inspected config files.

**View Coverage:**
```bash
Not configured in repository
```

## Test Types

**Unit Tests:**
- Dominant test type.
- Cover isolated transforms and helpers such as config parsing in `tests/only_rna/test_config.py`, QC/doublet/input reading in `tests/only_rna/test_processing.py`, embedding behavior in `tests/only_rna/test_embedding.py`, and discovery logic in `tests/only_rna/test_discovery.py`.

**Integration Tests:**
- Present as lightweight orchestration tests using monkeypatched stage boundaries.
- `tests/only_rna/test_tuning.py` verifies candidate evaluation and tuning artifact generation across multiple pipeline stages.
- `tests/only_rna/test_cli.py` verifies command dispatch, output-status behavior, and sample execution routing.
- `tests/test_gse192391_compare.py` exercises comparison-batch output generation and method-status contracts.

**E2E Tests:**
- Not detected as a dedicated framework.
- No browser/E2E runner, no subprocess-based full-suite integration harness, and no CI workflow file were observed in inspected files.

## Common Patterns

**Async Testing:**
```python
Not used; the inspected Python test suite is synchronous.
```

**Error Testing:**
```python
def test_run_azimuth_annotation_reports_errors_without_fallback_labels(monkeypatch):
    def _raise(*args, **kwargs):
        raise RuntimeError("Azimuth crashed")

    monkeypatch.setattr(azimuth_module, "_run_azimuth_r", _raise)
    result = azimuth_module.run_azimuth_annotation(...)

    assert result.status == "error"
    assert result.labels is None
```

**Filesystem contract testing:**
```python
sample_dir = write_sample_outputs(...)
assert (sample_dir / "metadata.csv").exists()
validation = pd.read_csv(sample_dir / "validation_result.csv")
assert "output_presence:qc_overview.png" in set(validation["check_name"])
```

**Fallback behavior testing:**
- Dependency-missing fallbacks are explicitly tested, for example:
  - Leiden/UMAP fallback behavior in `tests/only_rna/test_processing.py`
  - disabled/error Azimuth semantics in `tests/only_rna/test_azimuth.py`
  - scrublet fallback/normalization in `tests/only_rna/test_processing.py`

**Visualization contract testing:**
- Plot tests focus on semantic/layout outcomes rather than image snapshots, for example legend columns, marker scales, text labels, aspect ratio, and family-color alignment in `tests/only_rna/test_outputs.py`.

## Validation Mechanisms in Code

**Runtime validation:**
- `scripts/only_rna/outputs.py` writes `validation_result.csv` with checks for output presence and annotation method status.
- `scripts/only_rna/outputs.py` writes `qc_summary.csv` with counts and Azimuth status/detail.
- `scripts/only_rna/cli.py` computes `outputs_complete` by checking a fixed file contract via `_expected_output_paths()`.
- `scripts/only_rna/qc.py`, `scripts/only_rna/embedding.py`, and `scripts/only_rna/annotation.py` raise `KeyError` when required `adata.obs` columns are missing.

**Tested validation contracts:**
- `tests/only_rna/test_outputs.py` verifies required files, exclusion of removed legacy artifacts, status propagation into validation CSVs, and stale artifact cleanup.
- `tests/only_rna/test_cli.py` verifies `rna-status` output and the expected sample-root contract.

## Gaps and Quality Risks Visible From Tests

**Configuration/Test runner gaps:**
- No pytest config, no coverage config, and no CI workflow were detected, so test invocation and enforcement are convention-based.

**Coverage concentration:**
- The RNA path under `tests/only_rna/` is well covered for pure Python logic, but there is little visible automated coverage for non-RNA legacy scripts under `scripts/process/*.py` and R scripts such as `scripts/process/process_single_sample.R` and `scripts/process/process_single_rna_sample.R`.

**Type-check gap:**
- `pyrightconfig.json` excludes the main RNA implementation in `scripts/only_rna/`, so typed code there is validated mainly through tests, not static analysis.

**Contract brittleness:**
- Many tests assert exact filenames and output contracts, especially in `tests/only_rna/test_outputs.py` and `tests/only_rna/test_cli.py`. This is valuable, but it means output renames require coordinated updates across multiple tests and runtime validators.

**Fixture duplication:**
- Test helpers are repeated across files instead of centralized in `tests/conftest.py`, which increases maintenance cost when shared runtime contracts evolve.

**External integration realism:**
- Azimuth, plotting internals, and tuning stage composition are heavily monkeypatched. This keeps tests fast, but leaves a gap for fully integrated real-environment verification.

---

*Testing analysis: Mon Apr 20 2026*
