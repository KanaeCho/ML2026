# Output Contracts

本文件列出 AI 判断 RNA/ATAC/workflow 是否完成时应检查的输出文件和关键字段。

## RNA 样本输出

目录：

```text
output/rna/{GSE}/{sample_id}/
```

必需文件：

```text
metadata.csv
metadata_qc.csv
qc_summary.csv
qc_thresholds.json
validation_result.csv
qc_overview.png
umap_rna_pbmcref_vs_cima_l1.png
umap_rna_pbmcref_highlight.png
umap_rna_cima_l1.png
{sample_id}.h5ad
run_status.json
logs/sample_qc.log
```

关键字段：

- `pass_qc`
- `cluster`
- `umap_1`, `umap_2`
- `azimuth_cell_type`
- `azimuth_score`
- `azimuth_low_confidence`
- `azimuth_cima_l1`
- `final_celltype`
- `final_celltype_mapping`

成功判断：

- `run_status.json` 存在且 return code 为 0 或 status 为 success。
- expected outputs 完整。
- `metadata_qc.csv` 与 `.h5ad.obs` 对齐到 pass-QC 且 final celltype 已知细胞。
- UMAP 图可打开且 label 使用当前字段。

## ATAC 样本输出

目录之一：

```text
output/atac/{GSE}/{GSM}/
output/co/atac/{dataset}/{sample}/
output/longevity/atac/longevity/{sample_id}/
```

必需文件：

```text
metadata.csv
metadata_qc.csv
qc_summary.csv
validation_result.csv
qc_overview.png
umap_cima_cell_type_l1.png
umap_cima_cell_type_l2.png
{sample_id}.h5ad
run_status.json
logs/sample_qc.log
```

关键字段：

- `pass_qc`
- `nCount_ATAC`
- `nFeature_ATAC`
- `TSS.enrichment`
- `FRiP`
- `scDblFinder.score`
- `scDblFinder.class`
- `seurat_clusters`
- `umap_atac_1`, `umap_atac_2`
- `cima_ref_umap_1`, `cima_ref_umap_2`
- `cima_cell_type_l1`
- `cima_cell_type_l2`
- `cima_l1_low_confidence`
- `final_celltype`

成功判断：

- `run_status.json` 记录 success 或 outputs_complete true。
- `.h5ad` 存在。
- `metadata_qc.csv` 只含 pass-QC 细胞。
- L1/L2 UMAP 图存在。

## Product-level 输出

目录：

```text
output/1.only_atac/
output/2.only_rna/
output/3.co_atac/
output/4.co_rna/
```

必需文件：

```text
product_status.json
manifests/samples.csv
manifests/cells_metadata.csv
manifests/output_files.csv
qc/sample_qc_summary.csv
qc/validation_summary.csv
{product}.h5ad
integration/integration_summary.json
integration/integration_metrics.csv
integration/sample_mixing_summary.csv
figures/*_panels.png
```

关键字段：

- `integrated_umap_1`, `integrated_umap_2`
- `integrated_cluster`
- `integration_method`
- `integration_feature_space`
- `integration_included`
- `integration_exclusion_reason`
- RNA product: `integrated_cima_l1`, `integrated_cima_l1_score`, `integrated_cima_l2`, `integrated_cima_l2_score`

成功判断：

- `integrated_equals_original_umap` 应为 false。
- panel 图使用 `integrated_umap_1/2`，不是 per-sample UMAP。
- product `.h5ad` 是 metadata-level 对象，不应包含全量 count/peak matrix。
