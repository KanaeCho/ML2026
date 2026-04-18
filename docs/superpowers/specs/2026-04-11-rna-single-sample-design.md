# RNA Single-Sample Design

## Goal

On branch `only_rna`, the first-phase mainline is single-sample scRNA-seq processing for the RNA datasets already filtered in `datasets.xlsx` and already present under `/mnt/g/ML2026_data/raw/`.

Each supported RNA sample should complete:

- QC
- query-native RNA clustering and UMAP
- CIMA RNA `L1` and `L2` annotation
- sample-level validation outputs for UMAP quality review

ATAC and TEA-seq code stays in the repository but is no longer the branch acceptance target.

## Acceptance Scope

First-phase acceptance target:

- cover the RNA samples already selected in `datasets.xlsx` and downloaded locally
- focus on single-sample execution only
- prioritize `L1` and `L2` annotation quality

The branch does not need cross-sample RNA integration in this phase.

## Quality Standard

“Celltype UMAP quality is good enough” means both the plots and simple audit metrics are acceptable.

Primary evaluation:

- `L1` labels should be visually readable on query-native RNA UMAP
- `L2` labels should be mostly stable without broad mixed clouds
- low-confidence labels should be explicit and reviewable
- cluster-to-label consistency should be summarized numerically

Audit categories:

- sample-level QC retention
- low-confidence fraction
- `L1` cluster purity
- `L2` cluster purity
- label mixing within clusters

## Supported Inputs

The RNA workflow should discover and process these local input styles:

- 10x Matrix Market triplets: `matrix.mtx(.gz)` + `barcodes.tsv(.gz)` + `features.tsv(.gz)` or `genes.tsv(.gz)`
- 10x `.h5`
- tar archives containing `matrix.mtx`, `barcodes.tsv`, and `genes.tsv`/`features.tsv`

Discovery must prefer PBMC files when one GSM has multiple tissues, such as `GSE226039`.

Shared count matrices that are not cell-level scRNA-seq inputs should be discoverable as unsupported rather than silently treated as valid single-cell inputs.

## Runtime Assets

The workflow should support the current workspace layout where `data/` may be absent locally and the real data root is `/mnt/g/ML2026_data`.

RNA reference assets are already available in `/mnt/g/ML2026_data/reference/cima/`:

- `CIMA_RNA_6484974cells_36326genes_compressed.h5ad`
- `cima_rna_reference_pca_features.tsv.gz`
- `cima_rna_celltype_hierarchy.csv`
- `cima_rna_reference_l1_centroids.tsv`
- `cima_rna_reference_l2_centroids.tsv`
- `cima_rna_reference_model.json`

The first phase should use the existing RNA compact reference assets directly.

## Public Interfaces

`pipeline.py` should expose RNA-specific commands:

- `discover-rna`
- `run-rna-sample`
- `run-rna-gse`
- `rna-status`

The single-sample RNA entrypoint should be `scripts/process/process_single_rna_sample.R`.

## Output Layout

RNA outputs should live under:

- `output/rna/{GSE}/{sample_id}/`

Required outputs:

- `qc_overview.png`
- `metadata.csv`
- `metadata_qc.csv`
- `qc_summary.csv`
- `validation_result.csv`
- `umap_rna_clusters.png`
- `umap_rna_cima_cell_type_l1.png`
- `umap_rna_cima_cell_type_l2.png`
- `umap_rna_cima_cell_type_l1_masked.png`
- `matrix/matrix.mtx`
- `matrix/barcodes.tsv.gz`
- `matrix/features.tsv.gz`
- `{sample_id}_seurat_qc.rds`

## Validation Columns

`validation_result.csv` should include at least:

- `cell_barcode`
- `seurat_clusters`
- `umap_rna_1`
- `umap_rna_2`
- `cima_cell_type_l1`
- `cima_cell_type_l2`
- `cima_cell_type_l1_masked`
- `cima_l1_low_confidence`
- `cima_l1_cluster_purity`
- `cima_l2_cluster_purity`
- `cima_l2_score`
- `cima_l2_score_margin`

## Implementation Order

1. Update branch documentation so `AGENTS.md` matches the RNA mainline.
2. Add Python-side RNA sample discovery and pipeline commands.
3. Add the single-sample RNA R workflow using current local reference assets.
4. Verify the workflow on a small set of representative local RNA datasets.
