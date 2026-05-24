# 结果文件索引

本文档整理当前项目中已确认的结果文件类型和路径模式。由于真实 `output/` 目录是外部挂载或软链接，本文档优先记录输出契约和路径模式；具体每个样本的真实文件清单后续需要只读索引外部结果目录后补充。

## RNA QC 结果

| 路径模式 | 分析步骤 | 生成脚本 | 是否最终结果 | 是否适合 PPT | 备注 |
| --- | --- | --- | --- | --- | --- |
| `output/rna/{GSE}/{sample_id}/qc_summary.csv` | RNA QC summary | `scripts/only_rna/outputs.py` | 是 | 适合做表格 | 样本级 QC 和注释摘要 |
| `output/rna/{GSE}/{sample_id}/qc_thresholds.json` | RNA QC threshold audit | `scripts/only_rna/outputs.py` | 是 | 可用于方法补充 | 记录动态阈值 |
| `output/rna/{GSE}/{sample_id}/qc_overview.png` | RNA QC 可视化 | `scripts/only_rna/plotting.py` | 是 | 适合 | QC overview |

## ATAC QC 结果

| 路径模式 | 分析步骤 | 生成脚本 | 是否最终结果 | 是否适合 PPT | 备注 |
| --- | --- | --- | --- | --- | --- |
| `output/atac/{GSE}/{GSM}/qc_summary.csv` | ATAC QC summary | `process_single_sample.R` | 是 | 适合做表格 | only ATAC |
| `output/co/atac/{dataset}/{sample}/qc_summary.csv` | ATAC QC summary | `process_single_sample.R` | 是 | 适合 | co ATAC |
| `output/longevity/atac/longevity/{sample_id}/qc_summary.csv` | ATAC QC summary | `process_single_sample.R` | 是 | 适合 | longevity ATAC |
| `*/qc_overview.png` | ATAC QC 可视化 | `process_single_sample.R` | 是 | 适合 | QC overview |

## RNA 聚类和投影结果

| 路径模式 | 内容 | 是否适合 PPT | 备注 |
| --- | --- | --- | --- |
| `metadata.csv` / `metadata_qc.csv` | 包含 `cluster`, `umap_1`, `umap_2` | 表格不直接放 PPT | 用于审计和后续绘图 |
| `umap_rna_pbmcref_vs_cima_l1.png` | pbmcref vs CIMA L1 UMAP | 适合 | RNA 主结果图 |
| `umap_rna_pbmcref_highlight.png` | pbmcref highlight UMAP | 适合 | 注释可读性检查 |
| `umap_rna_cima_l1.png` | RNA CIMA L1 UMAP | 适合 | 5 类 final celltype 展示 |

## ATAC 聚类和投影结果

| 路径模式 | 内容 | 是否适合 PPT | 备注 |
| --- | --- | --- | --- |
| `metadata.csv` / `metadata_qc.csv` | 包含 `seurat_clusters`, `umap_atac_1/2`, `cima_ref_umap_1/2` | 表格不直接放 PPT | 用于审计 |
| `umap_cima_cell_type_l1.png` | ATAC CIMA L1 UMAP | 适合 | ATAC 主结果图 |
| `umap_cima_cell_type_l2.png` | ATAC CIMA L2 UMAP | 适合 | ATAC 细分注释图 |

## 注释结果

| 路径模式 | 内容 | 备注 |
| --- | --- | --- |
| RNA `metadata.csv` | Azimuth/pbmcref 和 final celltype 字段 | 全细胞审计 |
| RNA `metadata_qc.csv` | QC 后且 final celltype 已知细胞 | 下游结果主要使用 |
| ATAC `metadata.csv` | CIMA ATAC L1-L4 字段 | 全细胞审计 |
| ATAC `validation_result.csv` | ATAC pass-QC 关键字段 | legacy/product 可能使用 |

## marker 结果

当前项目中未找到 RNA marker gene 或 ATAC marker 主线输出。需要人工确认 notebook 或外部结果中是否存在。

## metadata 文件

| 文件 | 内容 |
| --- | --- |
| `metadata.csv` | 全细胞 metadata |
| `metadata_qc.csv` | QC 后细胞 metadata，RNA 当前还会对 final celltype 已知细胞做过滤 |
| `validation_result.csv` | 输出验证或关键字段导出 |
| product `manifests/cells_metadata.csv` | product-level cell metadata |
| product `manifests/samples.csv` | product-level sample manifest |

## 批处理输出

| 文件 | 内容 |
| --- | --- |
| `run_status.json` | 单样本运行状态 |
| `logs/sample_qc.log` | 单样本运行日志 |
| product `product_status.json` | product organization/integration 状态 |
| `tuning/candidates.csv` | RNA tuning candidate 分数 |
| `tuning/selection_summary.json` | RNA tuning 选择结果 |
| `tuning/selected_params.json` | RNA tuning 选中参数 |

## h5ad/RDS/loom 等对象

| 类型 | 当前状态 | Git 建议 |
| --- | --- | --- |
| `.h5ad` | RNA/ATAC 输出对象 | 不上传 Git，只索引路径 |
| `.rds` | ATAC/ArchR/Seurat 相关对象可能存在 | 不上传 Git |
| `.loom` | 当前项目中未确认 | 如果存在，不上传 Git |

## fragment 文件

fragment 是 ATAC 原始或近原始输入，不应上传 Git。

路径模式：

```text
data_root/raw/{GSE}/{GSM...fragments...}.tsv.gz
data_root/raw/{dataset}/{sample}/ATAC/fragments.tsv.gz
data_root/raw/longevity/atac/*_fragments.tsv.gz
```

## peak matrix

ATAC peak-by-cell matrix 由 `process_single_sample.R` 生成，可能位于：

```text
{sample_output}/matrix/matrix.mtx
{sample_output}/matrix/barcodes.tsv(.gz)
{sample_output}/matrix/features.tsv(.gz)
```

这些是大型中间/结果文件，不建议上传 Git。

## gene activity matrix

当前项目中未找到 gene activity matrix 主线输出。

## motif 或 peak annotation 结果

当前项目中未找到 motif enrichment 主线输出。peak annotation 在 ATAC 脚本中用于构建 ChromatinAssay，但未确认有单独结果文件。
