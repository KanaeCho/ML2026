# 输入输出说明

本文档整理当前 RNA/ATAC 工作流已经确认的输入和输出结构。目的是在可复用文档中替代本地绝对路径，让以后新项目可以直接参考这个结构准备数据、配置样本表和查找结果。

## 数据根目录解析规则

当前代码按以下优先级解析数据根目录：

1. 仓库根目录下的 `./data`，如果存在。
2. 环境变量 `ML2026_DATA_ROOT`，如果设置。
3. 当前工作区 fallback 路径 `/mnt/g/ML2026_data`。

第三项是当前机器上的本地路径，不应作为可移植 Git 项目的默认假设。后续复用时建议设置 `ML2026_DATA_ROOT`，或复制 `configs/config_template.yaml` 为本地配置文件。

## 输出根目录

默认输出根目录通常为：

```text
output/
```

当前工作区中的 `output/` 被记录为外部输出挂载或软链接。完整输出内容不应提交到 Git。

## reference 输入

预期 reference 目录结构：

```text
data_root/reference/
├── datasets.xlsx
├── atac.xlsx
├── co.xlsx
├── peak.bed
├── longevity/
│   ├── atac_barcodes/
│   ├── atac_barcodes_param_contrast/
│   └── atac_barcode_overrides.csv
└── cima/
    ├── CIMA_RNA_6484974cells_36326genes_compressed.h5ad
    ├── CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad
    ├── cima_rna_celltype_hierarchy.csv
    ├── cima_rna_reference_pca_features.tsv.gz
    ├── cima_rna_reference_l1_centroids.tsv
    ├── cima_rna_reference_l2_centroids.tsv
    ├── cima_atac_celltype_hierarchy.csv
    ├── cima_atac_reference_lsi_features.tsv.gz
    ├── cima_atac_reference_l1_centroids.tsv
    └── cima_atac_reference_l2_centroids.tsv
```

部分 RNA/ATAC CIMA L3/L4 centroid 文件也会被相关流程使用。公开发布 reference 设置说明前，需要再确认完整必需文件清单。

## RNA 输入

only RNA 分支从以下文件发现已筛选 scRNA 样本：

```text
data_root/reference/datasets.xlsx
```

支持的 raw 输入模式包括：

```text
data_root/raw/{GSE}/{GSM...}_matrix.mtx(.gz)
data_root/raw/{GSE}/{GSM...}_barcodes.tsv(.gz)
data_root/raw/{GSE}/{GSM...}_features.tsv(.gz)
data_root/raw/{GSE}/{GSM...}_genes.tsv(.gz)
data_root/raw/{GSE}/{GSM...}.h5
data_root/raw/{GSE}/{GSM...}_matrix.tar.gz
data_root/raw/{GSE}/{GSE}_matrix.mtx(.gz)
```

已确认的特殊规则：

- `GSE226039`：当前 discovery 逻辑只保留 PBMC 文件。
- `GSE198533`：共享 gene-count CSV 当前标记为不支持，不作为单细胞矩阵输入。
- GSE-level shared triplet 可以被发现，但显式 `run-rna-sample` 可能拒绝 `gse_shared` 目标；这类情况应使用 GSE-level 命令。

## RNA 输出

默认输出结构：

```text
output/rna/{GSE}/{sample_id}/
├── metadata.csv
├── metadata_qc.csv
├── qc_summary.csv
├── qc_thresholds.json
├── validation_result.csv
├── qc_overview.png
├── umap_rna_pbmcref_vs_cima_l1.png
├── umap_rna_pbmcref_highlight.png
├── umap_rna_cima_l1.png
├── {sample_id}.h5ad
├── run_status.json
└── logs/sample_qc.log
```

RNA tuning 输出结构：

```text
output/rna/{GSE}/{sample_id}/tuning/
├── candidates.csv
├── selection_summary.json
├── selected_params.json
├── umap_rna_candidates_overview_cima_l1.png
├── umap_rna_candidates_overview_pbmcref.png
└── umap_rna_candidates_overview_pbmcref_highlight.png
```

candidate 嵌套输出路径：

```text
output/rna/{GSE}/{sample_id}/tuning/{candidate_id}/{GSE}/{sample_id}/
```

## ATAC 输入

only ATAC 样本发现读取：

```text
data_root/reference/atac.xlsx
```

预期 raw 输入模式：

```text
data_root/raw/{GSE}/{GSM...fragments...}.tsv.gz
```

ATAC 脚本也支持显式传入路径：

- `--fragment-file`
- `--barcode-file`
- `--filtered-metadata-file`
- `--singlecell-file`

barcode 来源优先级：

1. 显式或自动发现的 `filtered_barcodes`。
2. 显式或自动发现的 `filtered_metadata`。
3. 显式或自动发现的 `singlecell` metadata。
4. 从 fragment counts 推断 barcode。

## ATAC 输出

ATAC 输出可能位于以下分支之一：

```text
output/atac/{GSE}/{GSM}/
output/co/atac/{dataset}/{sample}/
output/longevity/atac/longevity/{sample_id}/
```

典型 full output 包括：

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

根据 output profile 和清理逻辑，部分路径还可能存在：

```text
matrix/matrix.mtx
matrix/barcodes.tsv(.gz)
matrix/features.tsv(.gz)
{sample_id}_seurat_qc.rds
```

大型 matrix 和 RDS 文件不要提交到 Git。

## 共测输入和输出

共测样本发现读取：

```text
data_root/reference/co.xlsx
```

预期列：

```text
sample,dataset,age,health,donor
```

预期 raw layout：

```text
data_root/raw/{dataset}/{sample}/RNA/matrix.mtx(.gz)
data_root/raw/{dataset}/{sample}/RNA/barcodes.tsv(.gz)
data_root/raw/{dataset}/{sample}/RNA/features.tsv(.gz)
data_root/raw/{dataset}/{sample}/ATAC/fragments.tsv.gz
```

可选 ATAC 辅助文件：

```text
filtered_barcodes.tsv.gz
filtered_metadata.csv.gz
singlecell.csv.gz
```

共测输出：

```text
output/co/rna/{dataset}/{sample}/
output/co/atac/{dataset}/{sample}/
```

当前项目说明：共测分支目前不做 paired-cell matching，也不做联合 RNA+ATAC embedding。

## longevity 输入和输出

longevity 是独立额外数据通道，不属于 only RNA、only ATAC 或 co 自动发现。

RNA 输入：

```text
data_root/raw/longevity/rna/*.h5ad
```

ATAC 输入：

```text
data_root/raw/longevity/atac/*_fragments.tsv.gz
```

barcode 预处理输出：

```text
data_root/reference/longevity/atac_barcodes/{sample_id}/
├── filtered_barcodes.tsv.gz
├── barcode_qc.csv.gz
└── summary.json
```

longevity 正式输出：

```text
output/longevity/rna/longevity/{sample_id}/
output/longevity/atac/longevity/{sample_id}/
```

参数对照输出：

```text
output/longevity_param_contrast/{candidate}/atac/longevity/{sample_id}/
data_root/reference/longevity/atac_barcodes_param_contrast/{candidate}/{sample_id}/
```

## product-level 输出

product organization 输出：

```text
output/1.only_atac/
output/2.only_rna/
output/3.co_atac/
output/4.co_rna/
```

每个 product 目录预期包含：

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

product-level 整合只使用低维 embedding，不合并全量 RNA count matrix 或 ATAC peak matrix。

## 不进入 Git 的文件

不要提交：

- `data/` 真实内容。
- `output/` 真实内容。
- 完整 `QualityControl/` 和 `ArchRLogs/` 内容，除非人工筛选小型示例。
- `.h5ad`、`.rds`、`.RData`、`.loom`、`.bam`、`.fastq.gz`、fragment 文件、大型 matrix 和大型压缩表。
- 包含本地绝对路径的 local config。

应在文档中描述这些文件，而不是把它们直接提交到 Git。
