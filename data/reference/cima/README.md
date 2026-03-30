# CIMA reference assets

Kept reference input:

- `CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`
  - Source: `https://ftp.cngb.org/pub/SciRAID/trueblood/cima/CIMA_Resource/Cell_Atlas/CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`
  - Role: authoritative CIMA scATAC reference atlas with per-cell `cell_type_l1` through `cell_type_l4` in AnnData `obs`

Derived runtime assets used by the single-sample pipeline:

- `cima_atac_celltype_hierarchy.csv`
  - Derived from the reference `.h5ad`
  - Role: compact L1/L2/L3/L4 hierarchy table

- `cima_atac_reference_lsi_features.tsv.gz`
  - Derived from a balanced reference subset sampled from the `.h5ad`
  - Role: compact ATAC reference model containing selected feature indices, feature IDs, IDF weights, and LSI loadings for query projection and peak-identity validation

- `cima_atac_reference_l1_centroids.tsv`
- `cima_atac_reference_l2_centroids.tsv`
- `cima_atac_reference_l3_centroids.tsv`
- `cima_atac_reference_l4_centroids.tsv`
  - Derived from the same sampled reference subset
  - Role: per-level centroids in the shared reference LSI space used for hierarchical nearest-centroid label transfer

- `cima_atac_reference_model.json`
  - Role: metadata for the compact reference model build parameters

Generation script:

- `scripts/process/build_cima_reference_model.py`
  - Builds the compact reference model directly from the `.h5ad`

Removed exploratory / superseded assets:

- sample-level metadata tables
- scATAC proportion summary tables
- marker summary tables
- the old peakset-index transfer asset after replacing it with the centroid-based reference model

Rationale:

- The pipeline now uses a compact reference LSI model derived from real CIMA scATAC cells instead of direct peakset-overlap scoring.
- Keeping the reduced runtime files avoids carrying redundant downloads while preserving reproducible annotation inputs for this branch.
