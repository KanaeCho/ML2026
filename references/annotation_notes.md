# 注释说明

本文档整理当前项目中的 RNA/ATAC 注释逻辑和注意事项。

## RNA 注释

当前 RNA 主线注释方法：

```text
Azimuth pbmcref
```

相关代码：

- `scripts/only_rna/azimuth.py`
- `scripts/only_rna/annotation.py`
- `scripts/only_rna/final_celltype.py`
- `scripts/only_rna/outputs.py`

主要输出字段：

- `azimuth_cell_type`
- `azimuth_cell_type_l1_raw`
- `azimuth_cell_type_l2_raw`
- `azimuth_cima_l1_raw`
- `azimuth_cima_l1`
- `azimuth_score`
- `azimuth_score_margin`
- `azimuth_low_confidence`
- `final_celltype`
- `final_celltype_mapping`

当前 RNA final celltype 统一为五类：

- `CD4_T`
- `CD8_T`
- `B`
- `Myeloid`
- `NK`

注意：当前主线默认使用 Azimuth/pbmcref。`annotation.py` 中保留其他注释 helper，不代表这些都是当前主线必需产物。

## ATAC 注释

当前 ATAC 注释方法：

```text
CIMA ATAC compact feature model + centroid assignment
```

相关代码：

- `scripts/process/process_single_sample.R`
- `scripts/process/integrate_product_embeddings.py`

主要输出字段：

- `cima_cell_type_l1`
- `cima_cell_type_l2`
- `cima_cell_type_l3`
- `cima_cell_type_l4`
- `cima_cell_type_l1_masked`
- `cima_l1_low_confidence`
- `cima_l4_score`
- `cima_l4_score_margin`
- `cima_ref_umap_1`
- `cima_ref_umap_2`

## product-level RNA CIMA 标签

product-level RNA integration 会在 CIMA RNA PCA compact feature space 中计算诊断标签：

- `integrated_cima_l1`
- `integrated_cima_l1_score`
- `integrated_cima_l2`
- `integrated_cima_l2_score`

这些标签用于诊断 product-level projection，不替代单样本主线 Azimuth/pbmcref 注释。

## 注释结果排查

RNA 优先检查：

- `azimuth_status`
- `azimuth_score_mean`
- `azimuth_score_margin_mean`
- `azimuth_low_confidence_fraction`
- `final_celltype`
- `umap_rna_pbmcref_vs_cima_l1.png`

ATAC 优先检查：

- `cima_l1_low_confidence_frac`
- `median_cima_l4_score`
- `cima_cell_type_l1_masked`
- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`

## 需要人工判断的情况

- label 和预期生物学不一致。
- RNA 和 ATAC 注释不一致。
- low confidence 比例高。
- cluster 内多个 label 混杂。
- marker 证据和自动注释冲突。
