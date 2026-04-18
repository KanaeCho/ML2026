# only_rna Bounded Tuning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a bounded, auditable per-dataset tuning workflow for the only_rna pipeline that compares small preset families for QC, Azimuth `pbmcref`, and embedding/UMAP, selects a winner, and preserves the existing `run-rna-*` commands unchanged.

**Architecture:** Keep the current single-candidate RNA runner intact and make it more parameterizable. Add a separate tuning/orchestration layer that enumerates a small family of named presets, runs candidate combinations deterministically, collects audit metrics, selects the best candidate, and writes both the final outputs and selection artifacts. Reuse the robust Azimuth logic from the existing comparison script instead of the lightweight placeholder path in the main annotation module.

**Tech Stack:** Python 3.14, scanpy, anndata, pandas, scipy, matplotlib, pytest, existing R/Azimuth wrapper code in the repo

---

## File Structure

- Modify: `scripts/only_rna/models.py`
  - Add structured config models for embedding, Azimuth, tuning presets, and audit thresholds.
- Modify: `scripts/only_rna/config.py`
  - Load new config sections from YAML and merge CLI overrides into the active run config.
- Modify: `scripts/only_rna/embedding.py`
  - Replace hard-coded embedding heuristics with config-driven parameters while keeping safe defaults.
- Modify: `scripts/only_rna/annotation.py`
  - Keep CIMA logic, but route robust Azimuth execution through a reusable implementation instead of the lightweight placeholder.
- Create: `scripts/only_rna/azimuth.py`
  - Host robust Azimuth `pbmcref` runner and output normalization shared by compare/tuning paths.
- Create: `scripts/only_rna/tuning_presets.py`
  - Define bounded QC, Azimuth, and embedding preset families and named combinations.
- Create: `scripts/only_rna/tuning_metrics.py`
  - Compute candidate-level QC, annotation, and embedding audit scores and reason codes.
- Create: `scripts/only_rna/tuning_orchestrator.py`
  - Run candidate sweeps, select winners, and write tuning summary artifacts.
- Modify: `scripts/only_rna/outputs.py`
  - Add tuning/audit writers without breaking the current sample output contract.
- Modify: `scripts/only_rna/cli.py`
  - Add `tune-rna-sample` and `tune-rna-gse` commands while keeping current `run-*` commands stable.
- Modify: `scripts/process/pipeline.py`
  - Dispatch the new tuning commands to the Python-first only_rna CLI.
- Modify: `scripts/process/compare_gse192391_annotation_methods.py`
  - Reduce duplication by importing shared Azimuth helper(s) from `scripts/only_rna/azimuth.py` if possible.
- Modify/Create tests:
  - `tests/only_rna/test_config.py`
  - `tests/only_rna/test_embedding.py`
  - `tests/only_rna/test_multimethod.py`
  - `tests/only_rna/test_cli.py`
  - `tests/only_rna/test_outputs.py`
  - `tests/test_rna_pipeline.py`
  - `tests/only_rna/test_tuning.py` (new)
- Modify: `AGENTS.md`
  - Update branch-truth documentation only after implementation behavior is real.

The implementation order should be:
1. Config surface and tests
2. Parameterized embedding + robust Azimuth reuse
3. Presets + metrics + tuning orchestrator
4. CLI + pipeline wiring + audit outputs
5. Documentation update

### Task 1: Extend run config for bounded tuning surfaces

**Files:**
- Modify: `scripts/only_rna/models.py`
- Modify: `scripts/only_rna/config.py`
- Test: `tests/only_rna/test_config.py`

- [ ] **Step 1: Write a failing config test for new sections**

Add a test that loads a temporary YAML with `embedding`, `azimuth`, and `tuning` sections and asserts the returned config exposes those values with existing sections unchanged.

```python
def test_load_run_config_reads_embedding_azimuth_and_tuning_sections(tmp_path: Path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        """
qc:
  min_counts: 600
  min_genes: 350
  max_pct_mt: 18.0
  max_pct_ribo: 55.0
embedding:
  n_top_genes: 1500
  n_pcs: 20
  n_neighbors: 12
  resolution: 0.8
  min_dist: 0.35
azimuth:
  enabled: true
  reference: pbmcref
  annotation_levels:
    - l1
    - l2
  k_weight: 50
  mapping_score_k: 100
tuning:
  qc_preset_family: default
  azimuth_preset_family: default
  embedding_preset_family: default
  max_candidates: 9
""".strip()
    )

    config = load_run_config(config_path)

    assert config.qc.min_counts == 600
    assert config.embedding.n_top_genes == 1500
    assert config.embedding.min_dist == 0.35
    assert config.azimuth.enabled is True
    assert config.azimuth.reference == "pbmcref"
    assert config.azimuth.annotation_levels == ["l1", "l2"]
    assert config.tuning.max_candidates == 9
```

- [ ] **Step 2: Run the test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_config.py -k reads_embedding_azimuth_and_tuning_sections -v
```

Expected:
- FAIL because the current models/config loader do not expose these sections.

- [ ] **Step 3: Write a failing override-merging test for new nested fields**

```python
def test_merge_cli_overrides_updates_embedding_and_azimuth_fields():
    config = load_run_config(DEFAULT_CONFIG_PATH)

    updated = merge_cli_overrides(
        config,
        {
            "embedding__n_neighbors": 18,
            "embedding__min_dist": 0.2,
            "azimuth__k_weight": 80,
        },
    )

    assert updated.embedding.n_neighbors == 18
    assert updated.embedding.min_dist == 0.2
    assert updated.azimuth.k_weight == 80
```

- [ ] **Step 4: Run the override test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_config.py -k updates_embedding_and_azimuth_fields -v
```

Expected:
- FAIL because nested override merging currently only supports existing sections.

- [ ] **Step 5: Implement minimal config models**

Add dataclasses in `scripts/only_rna/models.py` for:

```python
@dataclass(frozen=True)
class EmbeddingConfig:
    n_top_genes: int = 2000
    n_pcs: int = 20
    n_neighbors: int = 15
    resolution: float = 1.0
    min_dist: float = 0.5
    spread: float = 1.0
    random_state: int = 0


@dataclass(frozen=True)
class AzimuthConfig:
    enabled: bool = False
    reference: str = "pbmcref"
    annotation_levels: tuple[str, ...] = ("l1", "l2")
    k_weight: int = 50
    n_trees: int = 20
    mapping_score_k: int = 100


@dataclass(frozen=True)
class TuningConfig:
    qc_preset_family: str = "default"
    azimuth_preset_family: str = "default"
    embedding_preset_family: str = "default"
    max_candidates: int = 9
```

and extend `RunConfig` accordingly.

- [ ] **Step 6: Implement minimal YAML loading and override merging**

Update `scripts/only_rna/config.py` so `load_run_config()` reads those sections and `merge_cli_overrides()` can set nested values for them using the existing `section__field` pattern.

- [ ] **Step 7: Run the targeted config tests to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_config.py -k "embedding_azimuth_and_tuning_sections or updates_embedding_and_azimuth_fields" -v
```

Expected:
- PASS for both targeted tests.

### Task 2: Parameterize embedding without breaking default behavior

**Files:**
- Modify: `scripts/only_rna/embedding.py`
- Test: `tests/only_rna/test_embedding.py`

- [ ] **Step 1: Write a failing test that run_embedding consumes explicit embedding config**

```python
def test_run_embedding_uses_explicit_configured_parameters(monkeypatch):
    adata = make_embedding_test_adata()
    config = load_run_config(DEFAULT_CONFIG_PATH)
    config = replace(config, embedding=replace(config.embedding, n_top_genes=123, n_neighbors=7, n_pcs=9, resolution=0.4, min_dist=0.15))

    captured = {}

    monkeypatch.setattr(scanpy.pp, "highly_variable_genes", lambda *args, **kwargs: captured.setdefault("hvg", kwargs))
    monkeypatch.setattr(scanpy.tl, "pca", lambda *args, **kwargs: captured.setdefault("pca", kwargs))
    monkeypatch.setattr(scanpy.pp, "neighbors", lambda *args, **kwargs: captured.setdefault("neighbors", kwargs))
    monkeypatch.setattr(scanpy.tl, "leiden", lambda *args, **kwargs: captured.setdefault("leiden", kwargs))
    monkeypatch.setattr(scanpy.tl, "umap", lambda *args, **kwargs: captured.setdefault("umap", kwargs))

    run_embedding(adata, config=config)

    assert captured["hvg"]["n_top_genes"] == 123
    assert captured["neighbors"]["n_neighbors"] == 7
    assert captured["neighbors"]["n_pcs"] == 9
    assert captured["leiden"]["resolution"] == 0.4
```

- [ ] **Step 2: Run the embedding parameter test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_embedding.py -k uses_explicit_configured_parameters -v
```

Expected:
- FAIL because `embedding.py` still derives most parameters internally.

- [ ] **Step 3: Implement config-driven embedding parameter selection**

Refactor `scripts/only_rna/embedding.py` so the effective parameters come from `config.embedding`, with existing dataset-size heuristics used only as backward-compatible fallbacks when explicit values are absent.

- [ ] **Step 4: Run the targeted embedding test to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_embedding.py -k uses_explicit_configured_parameters -v
```

Expected:
- PASS.

### Task 3: Extract robust Azimuth into a reusable module

**Files:**
- Create: `scripts/only_rna/azimuth.py`
- Modify: `scripts/only_rna/annotation.py`
- Modify: `scripts/process/compare_gse192391_annotation_methods.py`
- Test: `tests/only_rna/test_multimethod.py`

- [ ] **Step 1: Write a failing test that annotation can call the robust Azimuth path**

```python
def test_annotate_with_all_versions_routes_azimuth_through_shared_runner(monkeypatch):
    adata = make_annotation_test_adata()
    called = {}

    def fake_run_azimuth(adata, *, config, sample_id, output_dir):
        called["sample_id"] = sample_id
        adata.obs["azimuth_cell_type"] = "B"
        adata.obs["azimuth_l1"] = "B"
        adata.obs["azimuth_mapping_score"] = 0.9
        return adata

    monkeypatch.setattr("scripts.only_rna.annotation.run_azimuth_annotation", fake_run_azimuth)

    result = annotate_with_all_versions(
        adata,
        config=load_run_config(DEFAULT_CONFIG_PATH),
        sample_id="GSM1",
        output_dir=Path("/tmp/out"),
        methods=["cima", "azimuth"],
    )

    assert called["sample_id"] == "GSM1"
    assert "azimuth_cell_type" in result.obs.columns
```

- [ ] **Step 2: Run the annotation/Azimuth test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_multimethod.py -k routes_azimuth_through_shared_runner -v
```

Expected:
- FAIL because the current annotation module still uses the lightweight internal path.

- [ ] **Step 3: Create `scripts/only_rna/azimuth.py` with the robust runner**

Move the reusable parts of `_run_azimuth_r()` / `run_azimuth_annotation()` from `scripts/process/compare_gse192391_annotation_methods.py` into `scripts/only_rna/azimuth.py`, keeping support for:
- `reference='pbmcref'`
- `annotation_levels`
- `k.weight`
- `n.trees`
- `mapping.score.k`
- output columns for predicted labels and confidence/mapping scores

- [ ] **Step 4: Update `annotation.py` to import and use the shared runner**

Keep the CIMA logic unchanged, but route `methods=["azimuth"]` through the shared robust implementation.

- [ ] **Step 5: Run the targeted Azimuth test to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_multimethod.py -k routes_azimuth_through_shared_runner -v
```

Expected:
- PASS.

### Task 4: Define bounded preset families and scoring metrics

**Files:**
- Create: `scripts/only_rna/tuning_presets.py`
- Create: `scripts/only_rna/tuning_metrics.py`
- Create: `tests/only_rna/test_tuning.py`

- [ ] **Step 1: Write failing tests for preset-family construction**

```python
def test_default_tuning_presets_are_bounded_and_named():
    presets = default_tuning_presets()

    assert set(presets.qc.keys()) == {"baseline", "strict", "lenient"}
    assert set(presets.azimuth.keys()) == {"baseline", "conservative", "smooth"}
    assert set(presets.embedding.keys()) == {"baseline", "separated", "stable"}
```

```python
def test_candidate_score_includes_qc_annotation_embedding_and_reason_code():
    score = summarize_candidate_score(
        qc_score=0.8,
        annotation_score=0.7,
        embedding_score=0.6,
        reason_code="balanced_default",
    )

    assert score.total_score == pytest.approx(0.7)
    assert score.reason_code == "balanced_default"
```

- [ ] **Step 2: Run the tuning unit tests to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_tuning.py -v
```

Expected:
- FAIL because the preset and metrics modules do not exist yet.

- [ ] **Step 3: Implement minimal preset families and score objects**

Create small named preset families only:
- QC: `baseline`, `strict`, `lenient`
- Azimuth: `baseline`, `conservative`, `smooth`
- Embedding: `baseline`, `separated`, `stable`

and a minimal candidate score model that preserves component scores, weighted total, and a reason code.

- [ ] **Step 4: Run the tuning tests to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_tuning.py -v
```

Expected:
- PASS.

### Task 5: Add tuning orchestrator and audit output writers

**Files:**
- Create: `scripts/only_rna/tuning_orchestrator.py`
- Modify: `scripts/only_rna/outputs.py`
- Test: `tests/only_rna/test_tuning.py`

- [ ] **Step 1: Write a failing orchestrator test for candidate selection and artifacts**

```python
def test_run_tuning_selects_best_candidate_and_writes_selection_artifacts(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        "scripts.only_rna.tuning_orchestrator.evaluate_candidate",
        lambda *args, **kwargs: CandidateEvaluation(
            candidate_id=kwargs["candidate_id"],
            total_score={"baseline__baseline__baseline": 0.6, "strict__baseline__stable": 0.9}[kwargs["candidate_id"]],
            reason_code="test",
        ),
    )

    result = run_bounded_tuning(
        sample_id="GSM1",
        gse="GSE1",
        input_sample=object(),
        output_dir=tmp_path,
        config=load_run_config(DEFAULT_CONFIG_PATH),
    )

    assert result.best_candidate_id == "strict__baseline__stable"
    assert (tmp_path / "tuning" / "selection_summary.json").exists()
    assert (tmp_path / "tuning" / "candidates.csv").exists()
```

- [ ] **Step 2: Run the orchestrator test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_tuning.py -k selects_best_candidate_and_writes_selection_artifacts -v
```

Expected:
- FAIL because orchestrator and artifact writers do not yet exist.

- [ ] **Step 3: Implement minimal tuning orchestration**

Create `run_bounded_tuning(...)` that:
- enumerates a bounded cartesian product capped by `config.tuning.max_candidates`
- evaluates each candidate deterministically
- writes `tuning/candidates.csv`
- writes `tuning/selection_summary.json`
- writes `tuning/selected_params.json`
- returns the selected candidate metadata

- [ ] **Step 4: Run the targeted orchestrator test to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_tuning.py -k selects_best_candidate_and_writes_selection_artifacts -v
```

Expected:
- PASS.

### Task 6: Add tuning commands without breaking current run commands

**Files:**
- Modify: `scripts/only_rna/cli.py`
- Modify: `scripts/process/pipeline.py`
- Test: `tests/only_rna/test_cli.py`
- Test: `tests/test_rna_pipeline.py`

- [ ] **Step 1: Write a failing CLI test for `tune-rna-sample` dispatch**

```python
def test_main_dispatches_tune_rna_sample(monkeypatch):
    called = {}

    def fake_run(args):
        called["sample_id"] = args.sample_id
        return 0

    monkeypatch.setattr("scripts.only_rna.cli.handle_tune_rna_sample", fake_run)

    exit_code = main(["tune-rna-sample", "--gse", "GSE167363", "--sample-id", "GSM5102900"])

    assert exit_code == 0
    assert called["sample_id"] == "GSM5102900"
```

- [ ] **Step 2: Run the CLI test to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_cli.py -k dispatches_tune_rna_sample -v
```

Expected:
- FAIL because the command does not exist yet.

- [ ] **Step 3: Implement new CLI and pipeline commands**

Add:
- `tune-rna-sample`
- `tune-rna-gse`

and wire them through `scripts/process/pipeline.py` while leaving existing `discover-rna`, `run-rna-sample`, `run-rna-gse`, and `rna-status` semantics unchanged.

- [ ] **Step 4: Run targeted CLI and pipeline tests to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_cli.py -k tune_rna -v
uv run --with pytest python -m pytest tests/test_rna_pipeline.py -k tune_rna -v
```

Expected:
- PASS for the new command-path coverage.

### Task 7: Verify the integrated slice and update branch docs

**Files:**
- Modify: `AGENTS.md`
- Test: integration-focused existing test files listed above

- [ ] **Step 1: Run the focused only_rna test slice**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_config.py tests/only_rna/test_embedding.py tests/only_rna/test_multimethod.py tests/only_rna/test_tuning.py tests/only_rna/test_cli.py tests/test_rna_pipeline.py -v
```

Expected:
- PASS.

- [ ] **Step 2: Run syntax verification for changed Python files**

Run:

```bash
python3 -m py_compile scripts/only_rna/models.py scripts/only_rna/config.py scripts/only_rna/embedding.py scripts/only_rna/annotation.py scripts/only_rna/azimuth.py scripts/only_rna/tuning_presets.py scripts/only_rna/tuning_metrics.py scripts/only_rna/tuning_orchestrator.py scripts/only_rna/outputs.py scripts/only_rna/cli.py scripts/process/pipeline.py tests/only_rna/test_config.py tests/only_rna/test_embedding.py tests/only_rna/test_multimethod.py tests/only_rna/test_tuning.py tests/only_rna/test_cli.py tests/test_rna_pipeline.py
```

Expected:
- no output.

- [ ] **Step 3: Update `AGENTS.md` with the real implemented behavior only**

Document only what is true after code and tests pass:
- new tuning commands
- bounded preset approach
- Azimuth integration status
- new tuning outputs/audit artifacts

- [ ] **Step 4: Re-run the focused only_rna test slice after the doc update**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_config.py tests/only_rna/test_embedding.py tests/only_rna/test_multimethod.py tests/only_rna/test_tuning.py tests/only_rna/test_cli.py tests/test_rna_pipeline.py -v
```

Expected:
- PASS.

## Self-Review

- Spec coverage: This plan covers the approved bounded workflow: config surfaces, parameterized embedding, robust Azimuth reuse, bounded preset families, candidate scoring, orchestration, audit outputs, stable current run commands, new tuning commands, and AGENTS sync.
- Placeholder scan: No `TODO`/`TBD` placeholders remain; each task contains concrete files, tests, and commands.
- Type consistency: The plan uses `EmbeddingConfig`, `AzimuthConfig`, `TuningConfig`, `CandidateEvaluation`, and `run_bounded_tuning(...)` consistently across later tasks.
