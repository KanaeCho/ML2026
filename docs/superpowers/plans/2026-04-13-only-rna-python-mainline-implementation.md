# only_rna Python Mainline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a new Python-first single-sample scRNA-seq pipeline under `scripts/only_rna/` that preserves current RNA discovery/CLI contracts where required, supports GSE-level shared triplets through `run-rna-gse`, and produces the full RNA output family with `.h5ad` as the primary object.

**Architecture:** Keep the public RNA command family close to the current interface, but move all new processing into a dedicated `scripts/only_rna/` package with focused submodules for discovery, input reading, QC, doublet handling, embedding, annotation, plotting, outputs, and config. Lock behavior down with layered tests: discovery/CLI contract tests, processing tests, and output contract tests, then wire the package into the existing command surface.

**Tech Stack:** Python 3.11+, pytest, scanpy, anndata, scipy, pandas, numpy, matplotlib, pyyaml.

---

## Planned File Structure

### Create

- `scripts/only_rna/__init__.py` — package exports
- `scripts/only_rna/models.py` — lightweight dataclasses / typed objects for discovery and run config
- `scripts/only_rna/config.py` — YAML loading and CLI override merging
- `scripts/only_rna/discovery.py` — RNA dataset filtering, local-layout discovery, unsupported reason handling
- `scripts/only_rna/read_inputs.py` — triplet / h5 / archive / shared-triplet readers returning `AnnData`
- `scripts/only_rna/qc.py` — QC metric calculation and `pass_qc` assembly
- `scripts/only_rna/doublet.py` — Scrublet integration and normalized doublet fields
- `scripts/only_rna/embedding.py` — normalize / HVG / PCA / neighbors / clustering / UMAP
- `scripts/only_rna/annotation.py` — CIMA L1/L2 annotation on `pass_qc`
- `scripts/only_rna/plotting.py` — standardized UMAP export with readable legend rules
- `scripts/only_rna/outputs.py` — `.h5ad`, CSVs, matrix exports, validation outputs
- `scripts/only_rna/cli.py` — command handlers for `discover-rna`, `run-rna-sample`, `run-rna-gse`, `rna-status`
- `scripts/only_rna/default_config.yaml` — versioned defaults for QC, doublet, embedding, plotting
- `tests/only_rna/test_discovery.py`
- `tests/only_rna/test_processing.py`
- `tests/only_rna/test_outputs.py`
- `tests/only_rna/test_cli.py`

### Modify

- `pyproject.toml` — add `scanpy`, `anndata`, `pyyaml`, and any required scanpy/scrublet runtime deps
- `scripts/process/pipeline.py` — route RNA commands into `scripts.only_rna.cli`
- `AGENTS.md` — update only after behavior is implemented and verified

---

### Task 1: Add Python RNA Dependencies

**Files:**
- Modify: `pyproject.toml`

- [ ] **Step 1: Write the failing dependency smoke test list into the plan execution notes**

Create this checklist in your working notes before editing `pyproject.toml`:

```text
Expected imports after dependency update:
- import anndata
- import scanpy
- import yaml
```

- [ ] **Step 2: Add the minimal runtime dependencies to `pyproject.toml`**

Use entries equivalent to:

```toml
dependencies = [
  "h5py",
  "matplotlib",
  "numpy",
  "openpyxl",
  "pandas",
  "scikit-learn",
  "scipy",
  "anndata",
  "scanpy",
  "pyyaml",
]
```

Keep existing dependencies; only append the missing RNA pipeline requirements.

- [ ] **Step 3: Run a Python import smoke check**

Run:

```bash
python - <<'PY'
import anndata  # noqa: F401
import scanpy  # noqa: F401
import yaml  # noqa: F401
print('imports-ok')
PY
```

Expected: prints `imports-ok`.

- [ ] **Step 4: Commit the dependency update**

Run:

```bash
git add pyproject.toml
git commit -m "build: add only_rna Python dependencies"
```

### Task 2: Create Core Data Models and Config Loader

**Files:**
- Create: `scripts/only_rna/__init__.py`
- Create: `scripts/only_rna/models.py`
- Create: `scripts/only_rna/config.py`
- Create: `scripts/only_rna/default_config.yaml`
- Test: `tests/only_rna/test_processing.py`

- [ ] **Step 1: Write the failing config and model tests**

Add tests like:

```python
from pathlib import Path

from scripts.only_rna.config import load_default_config, merge_cli_overrides
from scripts.only_rna.models import QcThresholds, RunConfig


def test_load_default_config_returns_run_config():
    config = load_default_config(Path('scripts/only_rna/default_config.yaml'))
    assert isinstance(config, RunConfig)
    assert isinstance(config.qc, QcThresholds)
    assert config.plotting.umap_width > 0


def test_merge_cli_overrides_updates_selected_fields_only():
    base = load_default_config(Path('scripts/only_rna/default_config.yaml'))
    merged = merge_cli_overrides(base, min_genes=350, mt_max=18.0)
    assert merged.qc.min_genes == 350
    assert merged.qc.max_pct_mt == 18.0
    assert merged.qc.min_counts == base.qc.min_counts
```

- [ ] **Step 2: Run the new tests to confirm they fail first**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "config or override" -v
```

Expected: FAIL because module or symbols do not exist yet.

- [ ] **Step 3: Implement minimal typed models and config loading**

Create `scripts/only_rna/models.py` with dataclasses along these lines:

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class QcThresholds:
    min_counts: int
    min_genes: int
    max_pct_mt: float
    max_pct_ribo: float


@dataclass(frozen=True)
class PlottingConfig:
    umap_width: float
    umap_height: float
    dpi: int
    point_size: float
    legend_fontsize: float
    legend_title_fontsize: float


@dataclass(frozen=True)
class RunConfig:
    qc: QcThresholds
    plotting: PlottingConfig
```

Create `scripts/only_rna/config.py` with `load_default_config()` and `merge_cli_overrides()` that return a `RunConfig`, not raw dicts.

- [ ] **Step 4: Add a versioned YAML config file**

Create `scripts/only_rna/default_config.yaml` with concrete defaults, for example:

```yaml
qc:
  min_counts: 500
  min_genes: 300
  max_pct_mt: 20.0
  max_pct_ribo: 60.0
plotting:
  umap_width: 10.0
  umap_height: 8.0
  dpi: 180
  point_size: 3.0
  legend_fontsize: 10.0
  legend_title_fontsize: 11.0
```

- [ ] **Step 5: Re-run the config tests until they pass**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "config or override" -v
```

Expected: PASS.

- [ ] **Step 6: Commit the config foundation**

Run:

```bash
git add scripts/only_rna/__init__.py scripts/only_rna/models.py scripts/only_rna/config.py scripts/only_rna/default_config.yaml tests/only_rna/test_processing.py
git commit -m "feat: add only_rna config foundation"
```

### Task 3: Implement Discovery Contract

**Files:**
- Create: `scripts/only_rna/discovery.py`
- Test: `tests/only_rna/test_discovery.py`

- [ ] **Step 1: Write the failing discovery contract tests**

Add tests covering the agreed behavior:

```python
def test_discovers_shared_gse_triplet_as_supported_gse_level_sample(...):
    samples = discover_rna_samples(...)
    assert len(samples) == 1
    assert samples[0].sample_id == 'GSE149689'
    assert samples[0].input_type == 'triplet'
    assert samples[0].sample_kind == 'gse_shared'


def test_marks_shared_gene_count_csv_as_unsupported(...):
    samples = discover_rna_samples(...)
    assert samples[0].supported is False
    assert 'gene-count' in samples[0].note.lower()


def test_prefers_pbmc_files_for_gse226039(...):
    samples = discover_rna_samples(...)
    assert all('PBMC' in sample.source_name.upper() for sample in samples)
```

- [ ] **Step 2: Run the discovery tests to verify RED state**

Run:

```bash
pytest tests/only_rna/test_discovery.py -v
```

Expected: FAIL because `scripts.only_rna.discovery` is not implemented yet.

- [ ] **Step 3: Implement discovery data structures and local-layout parsing**

Create `scripts/only_rna/discovery.py` with functions shaped like:

```python
def resolve_data_root(cwd: Path | None = None) -> Path: ...
def selected_rna_gses(reference_root: Path) -> list[str]: ...
def discover_rna_samples(raw_root: Path, selected_gses: list[str]) -> list[DiscoveredSample]: ...
```

Implement these rules exactly:

- visible `scRNA` rows only
- per-GSM triplet / h5 / archive supported
- GSE-shared triplet supported and tagged as `sample_kind='gse_shared'`
- shared `gene_count` / `gene_counts` CSV unsupported
- `GSE226039` filters to PBMC-only filenames

- [ ] **Step 4: Re-run discovery tests until they pass**

Run:

```bash
pytest tests/only_rna/test_discovery.py -v
```

Expected: PASS.

- [ ] **Step 5: Commit discovery support**

Run:

```bash
git add scripts/only_rna/discovery.py tests/only_rna/test_discovery.py
git commit -m "feat: add only_rna discovery contract"
```

### Task 4: Implement Input Readers

**Files:**
- Create: `scripts/only_rna/read_inputs.py`
- Test: `tests/only_rna/test_processing.py`

- [ ] **Step 1: Write failing reader tests for supported input kinds**

Add tests like:

```python
def test_read_triplet_returns_adata_with_sample_metadata(...):
    adata = read_sample_input(sample)
    assert adata.n_obs == 2
    assert adata.n_vars == 3
    assert adata.obs['sample_id'].nunique() == 1


def test_read_shared_gse_triplet_sets_gse_sample_id(...):
    adata = read_sample_input(sample)
    assert adata.obs['sample_id'].iloc[0] == 'GSE149689'
```

- [ ] **Step 2: Run reader tests to confirm they fail**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "read_triplet or shared_gse_triplet" -v
```

Expected: FAIL because the reader does not exist yet.

- [ ] **Step 3: Implement `read_sample_input()` for triplet, h5, archive, and shared triplet**

Create a function shaped like:

```python
def read_sample_input(sample: DiscoveredSample) -> ad.AnnData:
    ...
```

Ensure it:

- returns `AnnData`
- normalizes `obs_names` and `var_names`
- writes `gse`, `sample_id`, and `input_type` into `adata.obs`
- supports `features.tsv` and `genes.tsv`

- [ ] **Step 4: Re-run the reader tests until they pass**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "read_triplet or shared_gse_triplet" -v
```

Expected: PASS.

- [ ] **Step 5: Commit input reading support**

Run:

```bash
git add scripts/only_rna/read_inputs.py tests/only_rna/test_processing.py
git commit -m "feat: add only_rna input readers"
```

### Task 5: Implement QC and Doublet Filtering

**Files:**
- Create: `scripts/only_rna/qc.py`
- Create: `scripts/only_rna/doublet.py`
- Test: `tests/only_rna/test_processing.py`

- [ ] **Step 1: Write failing tests for QC fields and pass_qc logic**

Add tests like:

```python
def test_compute_qc_metrics_adds_expected_columns(adata, run_config):
    out = compute_qc_metrics(adata, run_config)
    for key in ['n_counts', 'n_genes', 'pct_mt', 'pct_ribo']:
        assert key in out.obs


def test_apply_qc_filters_marks_doublets_as_not_pass_qc(adata, run_config):
    adata.obs['is_doublet'] = [False, True, False]
    out = apply_qc_filters(adata, run_config)
    assert out.obs['pass_qc'].tolist() == [True, False, True]
    assert out.obs['fails_doublet'].tolist() == [False, True, False]
```

- [ ] **Step 2: Run QC tests to confirm failure**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "qc_metrics or pass_qc or doublet" -v
```

Expected: FAIL because QC modules are missing.

- [ ] **Step 3: Implement minimal QC metric calculation**

Create `scripts/only_rna/qc.py` with:

```python
def compute_qc_metrics(adata: ad.AnnData, config: RunConfig) -> ad.AnnData: ...
def apply_qc_filters(adata: ad.AnnData, config: RunConfig) -> ad.AnnData: ...
```

`apply_qc_filters()` must populate at least:

- `pass_qc`
- `fails_count_floor`
- `fails_gene_floor`
- `fails_mt_ceiling`
- `fails_ribo_ceiling`
- `fails_doublet`

- [ ] **Step 4: Implement minimal doublet integration**

Create `scripts/only_rna/doublet.py` with:

```python
def run_doublet_detection(adata: ad.AnnData, config: RunConfig) -> ad.AnnData: ...
```

First version rules:

- if scanpy scrublet runs successfully, write `doublet_score` and `is_doublet`
- if the sample is too small for meaningful scrublet execution in a unit test, allow the test to patch the function and assert output shape/columns only

- [ ] **Step 5: Re-run QC and doublet tests until they pass**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "qc_metrics or pass_qc or doublet" -v
```

Expected: PASS.

- [ ] **Step 6: Commit QC and doublet logic**

Run:

```bash
git add scripts/only_rna/qc.py scripts/only_rna/doublet.py tests/only_rna/test_processing.py
git commit -m "feat: add only_rna qc and doublet filtering"
```

### Task 6: Implement Embedding and CIMA Annotation

**Files:**
- Create: `scripts/only_rna/embedding.py`
- Create: `scripts/only_rna/annotation.py`
- Test: `tests/only_rna/test_processing.py`

- [ ] **Step 1: Write failing tests for pass_qc-only processing**

Add tests like:

```python
def test_embedding_writes_cluster_and_umap_for_pass_qc_only(adata, run_config):
    adata.obs['pass_qc'] = [True, False, True]
    out = run_embedding(adata, run_config)
    assert 'cluster' in out.obs
    assert 'umap_1' in out.obs
    assert out.obs.loc[out.obs['pass_qc'] == False, 'cluster'].isna().all()


def test_annotation_runs_only_on_pass_qc_cells(adata, run_config, tmp_path):
    adata.obs['pass_qc'] = [True, False, True]
    out = annotate_with_cima(adata, reference_dir=tmp_path)
    assert out.obs.loc[out.obs['pass_qc'], 'cima_l1'].notna().all()
    assert out.obs.loc[~out.obs['pass_qc'], 'cima_l1'].isna().all()
```

- [ ] **Step 2: Run the processing tests to verify failure**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "embedding or annotation or pass_qc" -v
```

Expected: FAIL because modules do not exist.

- [ ] **Step 3: Implement the embedding pipeline**

Create `scripts/only_rna/embedding.py` with:

```python
def run_embedding(adata: ad.AnnData, config: RunConfig) -> ad.AnnData:
    """Run normalization, HVG selection, PCA, neighbors, clustering, and UMAP on pass_qc cells only."""
```

Write results back into the full-cell `adata.obs` as:

- `cluster`
- `umap_1`
- `umap_2`

For non-`pass_qc` cells, those values should remain missing.

- [ ] **Step 4: Implement minimal CIMA reference loading and annotation**

Create `scripts/only_rna/annotation.py` with functions shaped like:

```python
def load_cima_reference(reference_dir: Path) -> CimaReference: ...
def annotate_with_cima(adata: ad.AnnData, reference_dir: Path) -> ad.AnnData: ...
```

The first implementation should:

- align genes against `cima_rna_reference_pca_features.tsv.gz`
- produce `cima_l1`, `cima_l2`
- produce `cima_l1_score`, `cima_l1_score_margin`, `cima_l2_score`, `cima_l2_score_margin`
- produce `cima_l1_low_confidence` and `cima_l1_masked`

- [ ] **Step 5: Re-run the processing tests until they pass**

Run:

```bash
pytest tests/only_rna/test_processing.py -k "embedding or annotation or pass_qc" -v
```

Expected: PASS.

- [ ] **Step 6: Commit embedding and annotation support**

Run:

```bash
git add scripts/only_rna/embedding.py scripts/only_rna/annotation.py tests/only_rna/test_processing.py
git commit -m "feat: add only_rna embedding and cima annotation"
```

### Task 7: Implement Plotting and Output Writers

**Files:**
- Create: `scripts/only_rna/plotting.py`
- Create: `scripts/only_rna/outputs.py`
- Test: `tests/only_rna/test_outputs.py`

- [ ] **Step 1: Write failing output contract tests**

Add tests like:

```python
def test_write_outputs_creates_required_rna_files(processed_adata, tmp_path):
    out_dir = write_sample_outputs(processed_adata, tmp_path, sample_id='GSM1')
    expected = [
        'metadata.csv',
        'metadata_qc.csv',
        'qc_summary.csv',
        'validation_result.csv',
        'umap_rna_clusters.png',
        'umap_rna_cima_cell_type_l1.png',
        'umap_rna_cima_cell_type_l2.png',
        'umap_rna_cima_cell_type_l1_masked.png',
        'sample.h5ad',
    ]
    for rel in expected:
        assert (out_dir / rel).exists()


def test_metadata_csv_contains_all_cells(processed_adata, tmp_path):
    out_dir = write_sample_outputs(processed_adata, tmp_path, sample_id='GSM1')
    metadata = pd.read_csv(out_dir / 'metadata.csv')
    assert len(metadata) == processed_adata.n_obs
```

- [ ] **Step 2: Run output tests to confirm failure**

Run:

```bash
pytest tests/only_rna/test_outputs.py -v
```

Expected: FAIL because plotting/output modules are missing.

- [ ] **Step 3: Implement standardized UMAP plotting**

Create `scripts/only_rna/plotting.py` with a helper shaped like:

```python
def save_categorical_umap(
    adata: ad.AnnData,
    color_key: str,
    output_path: Path,
    title: str,
    config: RunConfig,
) -> None:
    ...
```

The first implementation must explicitly set:

- figure size
- point size
- legend font size
- legend title font size
- legend location/layout
- DPI

- [ ] **Step 4: Implement output writers**

Create `scripts/only_rna/outputs.py` with a top-level function shaped like:

```python
def write_sample_outputs(
    adata: ad.AnnData,
    output_root: Path,
    gse: str,
    sample_id: str,
    config: RunConfig,
) -> Path:
    ...
```

It must write:

- `metadata.csv` as all cells
- `metadata_qc.csv` as `pass_qc` subset
- `qc_summary.csv`
- `validation_result.csv`
- `.h5ad`
- matrix export triplet under `matrix/`
- all required UMAP images

- [ ] **Step 5: Re-run output tests until they pass**

Run:

```bash
pytest tests/only_rna/test_outputs.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit plotting and outputs**

Run:

```bash
git add scripts/only_rna/plotting.py scripts/only_rna/outputs.py tests/only_rna/test_outputs.py
git commit -m "feat: add only_rna plotting and output writers"
```

### Task 8: Implement CLI Wiring and Existing Pipeline Integration

**Files:**
- Create: `scripts/only_rna/cli.py`
- Modify: `scripts/process/pipeline.py`
- Test: `tests/only_rna/test_cli.py`

- [ ] **Step 1: Write failing CLI contract tests**

Add tests like:

```python
def test_run_rna_sample_rejects_explicit_gse_shared_target(...):
    result = run_cli(['run-rna-sample', '--gse', 'GSE149689', '--sample-id', 'GSE149689'])
    assert result.exit_code != 0
    assert 'gse-level shared sample' in result.stderr.lower()


def test_run_rna_gse_includes_shared_triplet_when_present(...):
    result = run_cli(['run-rna-gse', '--gse', 'GSE149689'])
    assert result.exit_code == 0
    assert 'GSE149689' in result.stdout
```

- [ ] **Step 2: Run the CLI tests to verify they fail first**

Run:

```bash
pytest tests/only_rna/test_cli.py -v
```

Expected: FAIL because the CLI glue does not exist yet.

- [ ] **Step 3: Implement `scripts/only_rna/cli.py` command handlers**

Add functions shaped like:

```python
def cmd_discover_rna(args) -> int: ...
def cmd_run_rna_sample(args) -> int: ...
def cmd_run_rna_gse(args) -> int: ...
def cmd_rna_status(args) -> int: ...
```

Enforce these rules:

- `run-rna-sample` only accepts sample-level targets
- `run-rna-gse` can implicitly run a GSE-shared triplet
- status checks required files in the new output family

- [ ] **Step 4: Wire `scripts/process/pipeline.py` to delegate RNA commands into `scripts.only_rna.cli`**

Keep the command names stable. Do not rename the public command family.

- [ ] **Step 5: Re-run CLI tests until they pass**

Run:

```bash
pytest tests/only_rna/test_cli.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit CLI integration**

Run:

```bash
git add scripts/only_rna/cli.py scripts/process/pipeline.py tests/only_rna/test_cli.py
git commit -m "feat: wire only_rna commands into pipeline"
```

### Task 9: Run End-to-End Verification and Update Documentation

**Files:**
- Modify: `AGENTS.md`

- [ ] **Step 1: Run the full new test suite**

Run:

```bash
pytest tests/only_rna -v
```

Expected: PASS.

- [ ] **Step 2: Run syntax/import validation on the new package and pipeline entrypoint**

Run:

```bash
python -m py_compile scripts/only_rna/__init__.py scripts/only_rna/models.py scripts/only_rna/config.py scripts/only_rna/discovery.py scripts/only_rna/read_inputs.py scripts/only_rna/qc.py scripts/only_rna/doublet.py scripts/only_rna/embedding.py scripts/only_rna/annotation.py scripts/only_rna/plotting.py scripts/only_rna/outputs.py scripts/only_rna/cli.py scripts/process/pipeline.py
```

Expected: no output.

- [ ] **Step 3: Run discovery against the real workspace layout**

Run:

```bash
python scripts/process/pipeline.py discover-rna
```

Expected:

- lists selected RNA datasets
- retains `GSE226039` PBMC-only behavior
- shows unsupported shared gene-count CSV entries as unsupported

- [ ] **Step 4: Run a representative sample smoke test into a temporary output root**

Run:

```bash
python scripts/process/pipeline.py run-rna-sample --gse GSE167363 --sample-id GSM5102900 --output-root /tmp/ml2026_only_rna_smoke --force
```

Expected:

- exit code `0`
- creates `/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/`

- [ ] **Step 5: Verify required smoke outputs exist**

Check for:

```text
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/metadata.csv
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/metadata_qc.csv
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/qc_summary.csv
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/validation_result.csv
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/umap_rna_cima_cell_type_l1.png
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/umap_rna_cima_cell_type_l2.png
/tmp/ml2026_only_rna_smoke/GSE167363/GSM5102900/matrix/matrix.mtx
```

- [ ] **Step 6: Update `AGENTS.md` to match the implemented reality**

Document only what was actually implemented and verified:

- new `scripts/only_rna/` mainline
- current command behavior
- output object now using `.h5ad`
- any finalized field/output names if they differ from the draft spec

- [ ] **Step 7: Commit verification-backed docs update**

Run:

```bash
git add AGENTS.md
git commit -m "docs: update branch guide for only_rna Python mainline"
```

## Self-Review Checklist

- Spec coverage:
  - CLI stability: Task 8
  - GSE shared triplet support: Tasks 3, 4, 8
  - `run-rna-sample` restriction: Task 8
  - YAML + CLI config: Task 2
  - Python doublet hard filter: Task 5
  - pass_qc-only embedding and annotation: Task 6
  - output family + `.h5ad`: Task 7
  - legend/readability improvements: Task 7
  - layered tests: Tasks 2–8
  - AGENTS sync after implementation: Task 9
- Placeholder scan: no `TODO`, `TBD`, “similar to above”, or undefined task references remain.
- Type consistency: `RunConfig`, `DiscoveredSample`, `read_sample_input()`, `run_doublet_detection()`, `run_embedding()`, `annotate_with_cima()`, and `write_sample_outputs()` are named consistently across tasks.
