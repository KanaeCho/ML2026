# only_rna UMAP Visual Tuning Design

## Goal

Improve the readability of the generated only_rna UMAP figures without changing clustering results, annotation logic, output filenames, or the biological meaning of the plots.

This tuning pass specifically targets two user-visible issues already observed in the smoke outputs:

1. L2 legends become hard to read when category count is high.
2. The effective UMAP plotting region can be visually compressed into a narrow rectangle because the legend consumes too much horizontal space.

## Scope

This design applies only to the visualization layer used by `scripts/only_rna/plotting.py` and the tests that validate that layer.

In scope:

- improve L2 legend layout for many categories
- improve legend marker readability
- improve main scatter-point visibility
- reduce the chance that the main UMAP panel is squeezed into an unreadable rectangle
- add regression coverage for these readability guarantees

Out of scope:

- changing clustering results
- changing CIMA annotation logic
- changing output file names or output family
- changing smoke-run biological conclusions
- redesigning the full plotting stack or introducing a new plotting library

## User-approved direction

The approved direction is a **single-figure, moderate enhancement** approach.

That means:

- keep the current output contract and one-PNG-per-plot behavior
- improve readability substantially, but avoid a dramatic style change
- tune both the legend and the scatter layer together
- treat L2 as the highest-pressure plot and allow it to use more aggressive legend layout than simpler plots

The user additionally clarified that the L2 UMAP is currently being visually compressed into a rectangle, so the design must address both legend readability and layout pressure on the main plotting panel.

## Design

### 1. Main plotting area takes priority

The plotting function should stop treating the legend as an afterthought that can freely consume space from the scatter area. The main UMAP panel is the primary information surface, so the layout must preserve a readable plotting region even when the legend is large.

Concretely, the plotting function should:

- keep the scatter panel visually dominant
- reserve a more predictable legend area rather than letting a long single-column legend squeeze the axes
- avoid layouts where `tight_layout()` plus a dense legend turns the visible embedding into a narrow strip

### 2. Automatic multi-column legend policy

The legend should use an adaptive column strategy based on the number of categories.

Baseline policy:

- small category count: 1 column
- medium category count: 2 columns
- high category count: 3 columns

The current threshold logic already moves in this direction. This tuning pass should preserve that behavior and refine it as needed for L2-heavy plots.

The important requirement is not the exact threshold numbers but the resulting effect: the L2 legend must stop presenting as a long unreadable vertical wall that crushes the plotting area.

### 3. Legend marker readability must be decoupled from plot point size

The main root cause already diagnosed was that the legend marker size inherited the very small scatter point size, making the legend color swatches almost useless.

The design therefore requires:

- explicit legend marker scaling independent of main point size
- legend markers large enough to be visually matched to categories
- preservation of category text labels, but with color chips/markers that are actually usable

This remains a hard requirement for all categorical UMAPs and is especially important for L2.

### 4. Moderate point-visibility enhancement

The scatter points should become easier to see, but not so large that dense regions become muddy.

This pass should apply a moderate enhancement strategy:

- modestly increase point size from the current tiny baseline where needed
- keep alpha in a readable but not fully opaque range
- optionally allow a light edge/contrast treatment if needed to preserve point identity against dense local neighborhoods

The purpose is readability, not stylistic dramatization.

### 5. L2-specific tolerance for stronger layout tuning

Because L2 has the largest category counts, the plotting system should allow the L2 plot to be rendered with stricter legend/layout behavior than cluster or L1 plots if necessary.

This does not require different filenames or a different plotting API. It only means the plotting implementation may branch on category pressure and use a more generous legend strategy when the plotted label cardinality is large.

### 6. Configuration strategy

The existing `PlottingConfig` is currently minimal. This design allows small, focused extension of the display-layer config only if needed to support the improved layout.

Candidate optional plotting fields:

- `legend_markerscale`
- `legend_ncols`
- `legend_location`
- `legend_bbox_to_anchor`
- `point_alpha`
- `point_edge_width`
- `point_edge_color`

These fields are optional tuning controls. The implementation should remain backward-compatible with existing config instances and defaults.

### 7. Tests

Regression tests must explicitly validate readability-oriented behavior rather than just file existence.

At minimum the tests should cover:

1. legend marker scale is larger than raw scatter marker usage for many-category plots
2. multi-column legend behavior activates for high-cardinality category sets
3. the generated image is still produced successfully under many-category conditions
4. the L2-like many-category case does not regress into a single-column layout that predictably compresses the panel

The tests should remain implementation-aware but not overfit to incidental matplotlib internals beyond the legend kwargs the plotting code intentionally controls.

## Files expected to change

- `scripts/only_rna/plotting.py`
- `tests/only_rna/test_outputs.py`

Possible small supporting changes if absolutely necessary:

- `scripts/only_rna/models.py`
- `scripts/only_rna/default_config.yaml`

No other functional modules should change for this display-only tuning pass.

## Acceptance criteria

This design is satisfied when all of the following are true:

1. The L2 legend is visibly easier to read than the current version.
2. The legend markers are no longer effectively invisible relative to the legend text.
3. The UMAP panel is no longer visually squeezed into an obvious narrow rectangle under many-category L2 rendering.
4. Main points are easier to see without causing a major style shift.
5. Existing output filenames and output-writing code paths remain unchanged.
6. Biological results and annotation outputs remain unchanged; this pass only affects display.
7. Regression tests explicitly cover the intended readability behavior.

## Non-goals and guardrails

- Do not change cluster assignments.
- Do not change CIMA labels or confidence thresholds.
- Do not introduce separate “preview-only” rendering logic.
- Do not split one plot into multiple files.
- Do not redesign the project around a new visualization backend.

## Summary

This is a focused display-layer refinement. The approved strategy is to keep the current one-image output contract, but make the plots actually usable by improving three things together: legend marker readability, adaptive legend layout for high-cardinality L2 labels, and protection of the main UMAP panel from being visually crushed by the legend.
