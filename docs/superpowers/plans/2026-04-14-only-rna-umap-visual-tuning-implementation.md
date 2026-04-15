# only_rna UMAP Visual Tuning Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Improve only_rna UMAP readability by making L2 legends readable, increasing point visibility, and preventing the main UMAP panel from being visually squeezed into a narrow rectangle.

**Architecture:** Keep the existing single-PNG output contract and tune only the display layer. Concentrate the implementation in `scripts/only_rna/plotting.py`, with any config extension kept backward-compatible and tightly scoped. Validate behavior through focused regression tests in `tests/only_rna/test_outputs.py`.

**Tech Stack:** Python 3.14, matplotlib, pandas, anndata, pytest

---

## File Structure

- Modify: `scripts/only_rna/plotting.py`
  - Add layout logic that preserves a readable main UMAP panel under large legends.
  - Add moderate point-visibility improvements without changing biology/output contract.
  - Keep legend behavior adaptive for high-cardinality categorical plots.
- Modify: `tests/only_rna/test_outputs.py`
  - Add regression tests for L2-like many-category legend layout and panel-preservation behavior.
- Optional modify only if strictly needed: `scripts/only_rna/models.py`, `scripts/only_rna/default_config.yaml`
  - Only extend config if the plotting changes require a new explicit field that cannot be kept safely behind `getattr(...)` defaults.

The plan intentionally avoids touching clustering, annotation, output filenames, CLI wiring, or smoke-run biology.

### Task 1: Add regression coverage for many-category readability

**Files:**
- Modify: `tests/only_rna/test_outputs.py`
- Test: `tests/only_rna/test_outputs.py`

- [ ] **Step 1: Write the failing test for L2-style panel squeeze prevention**

Add this test near the existing legend-readability test:

```python
def test_save_categorical_umap_keeps_readable_panel_for_many_categories(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = ad.AnnData(
        X=sparse.csr_matrix(np.eye(27)),
        obs=pd.DataFrame(
            {
                "umap_1": np.linspace(-5.0, 5.0, 27),
                "umap_2": np.linspace(-3.0, 3.0, 27),
                "cima_l2": [f"L2_{i}" for i in range(27)],
            },
            index=[f"cell-{i}" for i in range(27)],
        ),
    )

    captured: dict[str, object] = {}
    original_subplots = plt.subplots

    def capture_subplots(*args, **kwargs):
        captured["figsize"] = kwargs.get("figsize")
        captured["dpi"] = kwargs.get("dpi")
        return original_subplots(*args, **kwargs)

    original_legend = matplotlib.axes.Axes.legend

    def capture_legend(self, *args, **kwargs):
        captured.update(
            {
                "ncol": kwargs.get("ncol"),
                "markerscale": kwargs.get("markerscale"),
                "bbox_to_anchor": kwargs.get("bbox_to_anchor"),
                "loc": kwargs.get("loc"),
            }
        )
        return original_legend(self, *args, **kwargs)

    monkeypatch.setattr(plt, "subplots", capture_subplots)
    monkeypatch.setattr(matplotlib.axes.Axes, "legend", capture_legend)

    output_path = tmp_path / "l2_many_categories.png"
    save_categorical_umap(
        adata,
        color_key="cima_l2",
        output_path=output_path,
        title="L2 readability",
        config=config,
    )

    assert output_path.exists()
    assert captured["markerscale"] >= 4.0
    assert captured["ncol"] >= 3
    assert captured["loc"] == "center left"
    assert captured["bbox_to_anchor"] is not None
```

- [ ] **Step 2: Run the new regression to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_outputs.py -k keeps_readable_panel_for_many_categories -v
```

Expected:
- FAIL because current plotting logic still uses the same base figure size for heavy L2 legends and has no stronger panel-preservation policy.

- [ ] **Step 3: Add a second failing test for point visibility defaults**

Append this test below the first one:

```python
def test_save_categorical_umap_uses_visibility_tuned_scatter_defaults(
    tmp_path: Path, monkeypatch
):
    config = _make_run_config()
    adata = _make_output_adata()

    captured: dict[str, object] = {}
    original_scatter = matplotlib.axes.Axes.scatter

    def capture_scatter(self, *args, **kwargs):
        captured.update(
            {
                "s": kwargs.get("s"),
                "alpha": kwargs.get("alpha"),
                "linewidths": kwargs.get("linewidths"),
            }
        )
        return original_scatter(self, *args, **kwargs)

    monkeypatch.setattr(matplotlib.axes.Axes, "scatter", capture_scatter)

    output_path = tmp_path / "visibility_tuned.png"
    save_categorical_umap(
        adata,
        color_key="cluster",
        output_path=output_path,
        title="Visibility tuning",
        config=config,
    )

    assert output_path.exists()
    assert captured["s"] > config.plotting.point_size
    assert captured["alpha"] < 0.95
```

- [ ] **Step 4: Run the second regression to verify RED**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_outputs.py -k visibility_tuned_scatter_defaults -v
```

Expected:
- FAIL because the current implementation still uses the raw `config.plotting.point_size` directly and does not apply a dedicated visibility enhancement policy.

- [ ] **Step 5: Commit the failing tests once they are both verified RED**

```bash
GIT_MASTER=1 git add tests/only_rna/test_outputs.py
GIT_MASTER=1 git commit -m "Add UMAP readability regression tests" -m "Ultraworked with [Sisyphus](https://github.com/code-yeongyu/oh-my-openagent)" -m "Co-authored-by: Sisyphus <clio-agent@sisyphuslabs.ai>"
```

### Task 2: Implement single-figure readability tuning in plotting.py

**Files:**
- Modify: `scripts/only_rna/plotting.py`
- Test: `tests/only_rna/test_outputs.py`

- [ ] **Step 1: Add a helper that computes adaptive legend layout for many-category plots**

Add this helper above `save_categorical_umap`:

```python
def _legend_layout(category_count: int, plotting: RunConfig | object) -> tuple[int, tuple[float, float]]:
    configured_ncols = int(getattr(plotting, "legend_ncols", 0))
    configured_anchor = getattr(plotting, "legend_bbox_to_anchor", None)

    if configured_ncols > 0:
        ncols = configured_ncols
    elif category_count >= 24:
        ncols = 3
    elif category_count >= 12:
        ncols = 2
    else:
        ncols = 1

    if configured_anchor is not None:
        anchor = configured_anchor
    elif ncols >= 3:
        anchor = (1.02, 0.5)
    elif ncols == 2:
        anchor = (1.01, 0.5)
    else:
        anchor = (1.0, 0.5)

    return ncols, anchor
```

- [ ] **Step 2: Replace the current inline legend-column logic with the helper**

Change the current block from:

```python
if legend_ncols <= 0:
    if len(categories) >= 24:
        legend_ncols = 3
    elif len(categories) >= 12:
        legend_ncols = 2
    else:
        legend_ncols = 1
```

to:

```python
legend_ncols, legend_bbox_to_anchor = _legend_layout(len(categories), plotting)
```

- [ ] **Step 3: Add moderate point-visibility tuning without changing biological output**

Inside the category loop, replace the current scatter call:

```python
ax.scatter(
    subset["umap_1"],
    subset["umap_2"],
    s=plotting.point_size,
    c=[cmap(idx)],
    label=category,
    linewidths=0,
    alpha=0.9,
)
```

with:

```python
scatter_size = float(getattr(plotting, "display_point_size", plotting.point_size * 2.0))
scatter_alpha = float(getattr(plotting, "point_alpha", 0.82))
scatter_linewidth = float(getattr(plotting, "point_edge_width", 0.1))

ax.scatter(
    subset["umap_1"],
    subset["umap_2"],
    s=scatter_size,
    c=[cmap(idx)],
    label=category,
    linewidths=scatter_linewidth,
    alpha=scatter_alpha,
)
```

This keeps the change moderate: visibly easier to read, but not a full style reset.

- [ ] **Step 4: Strengthen the legend block for L2-like category pressure**

Keep the existing legend call but ensure it uses the tuned values:

```python
legend = ax.legend(
    title=color_key,
    loc=legend_location,
    bbox_to_anchor=legend_bbox_to_anchor,
    ncol=legend_ncols,
    frameon=False,
    fontsize=plotting.legend_fontsize,
    title_fontsize=plotting.legend_title_fontsize,
    markerscale=legend_markerscale,
)
```

and immediately after `ax.set_ylabel("umap_2")`, add:

```python
ax.set_aspect("auto")
```

plus reserve more deterministic space before save:

```python
fig.tight_layout(rect=(0.0, 0.0, 0.82, 1.0))
```

Replace the existing plain `fig.tight_layout()` call with that `rect=...` form so the right-side legend area stops crushing the scatter panel.

- [ ] **Step 5: Run the targeted plotting regressions to verify GREEN**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_outputs.py -k "readable_legend or keeps_readable_panel_for_many_categories or visibility_tuned_scatter_defaults" -v
```

Expected:
- PASS for all targeted readability tests.

- [ ] **Step 6: Run the full plotting/output test file**

Run:

```bash
uv run --with pytest python -m pytest tests/only_rna/test_outputs.py -v
```

Expected:
- PASS for the entire file.

- [ ] **Step 7: Run syntax verification for the changed plotting files**

Run:

```bash
python3 -m py_compile scripts/only_rna/plotting.py tests/only_rna/test_outputs.py
```

Expected:
- no output

- [ ] **Step 8: Commit the readability implementation with its paired tests**

```bash
GIT_MASTER=1 git add scripts/only_rna/plotting.py tests/only_rna/test_outputs.py
GIT_MASTER=1 git commit -m "Improve only_rna UMAP readability" -m "Ultraworked with [Sisyphus](https://github.com/code-yeongyu/oh-my-openagent)" -m "Co-authored-by: Sisyphus <clio-agent@sisyphuslabs.ai>"
```

### Task 3: Regenerate preview outputs and verify the visible result

**Files:**
- Modify: `tmp_preview/GSE167363/GSM5102900/*.png` (generated artifacts, not source)
- Test: smoke-run outputs under `/tmp`

- [ ] **Step 1: Regenerate the smoke output with the readability-tuned plotting layer**

Run:

```bash
rm -rf "/tmp/ml2026_only_rna_smoke_readable_v2"
uv run python scripts/process/pipeline.py run-rna-sample --gse GSE167363 --sample-id GSM5102900 --output-root /tmp/ml2026_only_rna_smoke_readable_v2 --force
```

Expected:
- command exits successfully
- output directory exists under `/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900/`

- [ ] **Step 2: Verify the four required UMAP PNGs exist in the new smoke output**

Run:

```bash
python3 - <<'PY'
from pathlib import Path
base = Path('/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900')
for name in [
    'umap_rna_clusters.png',
    'umap_rna_cima_cell_type_l1.png',
    'umap_rna_cima_cell_type_l2.png',
    'umap_rna_cima_cell_type_l1_masked.png',
]:
    path = base / name
    print(name, path.exists())
PY
```

Expected:
- all four print `True`

- [ ] **Step 3: Refresh the project-local preview directory**

Run:

```bash
mkdir -p "/mnt/f/ydd/ML2026/tmp_preview/GSE167363/GSM5102900"
cp \
  "/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900/umap_rna_clusters.png" \
  "/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900/umap_rna_cima_cell_type_l1.png" \
  "/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900/umap_rna_cima_cell_type_l2.png" \
  "/tmp/ml2026_only_rna_smoke_readable_v2/GSE167363/GSM5102900/umap_rna_cima_cell_type_l1_masked.png" \
  "/mnt/f/ydd/ML2026/tmp_preview/GSE167363/GSM5102900/"
```

Expected:
- preview directory contains the refreshed four PNGs

- [ ] **Step 4: Confirm the preview directory contents**

Run:

```bash
ls -1 "/mnt/f/ydd/ML2026/tmp_preview/GSE167363/GSM5102900"
```

Expected output includes exactly:
- `umap_rna_clusters.png`
- `umap_rna_cima_cell_type_l1.png`
- `umap_rna_cima_cell_type_l2.png`
- `umap_rna_cima_cell_type_l1_masked.png`

- [ ] **Step 5: Commit the refreshed preview artifacts if the user wants them tracked**

Only if the user explicitly wants the preview PNGs committed:

```bash
GIT_MASTER=1 git add tmp_preview/GSE167363/GSM5102900/*.png
GIT_MASTER=1 git commit -m "Refresh only_rna UMAP preview images" -m "Ultraworked with [Sisyphus](https://github.com/code-yeongyu/oh-my-openagent)" -m "Co-authored-by: Sisyphus <clio-agent@sisyphuslabs.ai>"
```

If the user does **not** want preview PNGs tracked, skip this commit and leave the refreshed files uncommitted.

## Self-review checklist

- Spec coverage: This plan covers all approved design goals: L2 legend readability, larger legend markers, moderate point visibility enhancement, and prevention of the L2 panel being visually crushed by the legend.
- Placeholder scan: No `TODO`, `TBD`, `implement later`, or vague filler steps remain.
- Type consistency: The plan keeps `save_categorical_umap(...)` as the public API and does not invent a second plotting entrypoint.

## Execution note

This is a small display-layer tuning plan. It should be executed before any further biological or clustering retuning so that visual regressions remain easy to isolate.
