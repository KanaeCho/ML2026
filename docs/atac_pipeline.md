# ATAC Pipeline 文档

本文档整理当前项目中已经实现的 ATAC 分析流程。内容基于当前仓库中的 R/Python 脚本和项目说明，不把未确认步骤写成已完成。

## 总览

当前 ATAC 主线脚本为：

```text
scripts/process/process_single_sample.R
```

该脚本被以下分支复用：

- only ATAC
- co ATAC
- longevity ATAC

主要 Python 路由入口为：

```text
scripts/process/pipeline.py
scripts/co/cli.py
scripts/longevity/cli.py
```

## 1. 输入数据

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/pipeline.py`, `scripts/process/process_single_sample.R`, `scripts/co/cli.py`, `scripts/longevity/cli.py` |
| 输入类型 | fragments.tsv.gz、barcode 文件、filtered metadata、singlecell metadata、peak reference |
| only ATAC 样本表 | `data_root/reference/atac.xlsx` |
| co ATAC 样本表 | `data_root/reference/co.xlsx` |
| longevity ATAC | 扫描 `data_root/raw/longevity/atac/*_fragments.tsv.gz` |
| 当前状态 | 已完成 |
| 是否支持批处理 | 支持 |
| PPT 适用性 | 适合放数据结构页 |

## 2. 样本信息和 metadata

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/pipeline.py`, `scripts/co/cli.py`, `scripts/longevity/cli.py` |
| 样本字段 | `gse`, `gsm` 或 `sample_id`, `individual_id`, `age`, `health`, `donor` |
| 输出 metadata | `metadata.csv`, `metadata_qc.csv`, `validation_result.csv`, `.h5ad.obs` |
| 当前状态 | 已完成 |
| 不确定项 | co 分支 RNA 和 ATAC barcode 是否一一对应未确认；当前项目不做 paired-cell matching |

## 3. fragments 或 peak matrix 输入

| 项目 | 内容 |
| --- | --- |
| fragments | 默认从 raw GSE/GSM 文件名发现，也可用 `--fragment-file` 显式传入 |
| barcode 来源 | `filtered_barcodes`, `filtered_metadata`, `singlecell`, fragment-count inference |
| peak 文件 | `data_root/reference/peak.bed` |
| peak-by-cell matrix | 由 Signac `FeatureMatrix` 构建 |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合解释 fragments 和 peak matrix 关系 |

## 4. ATAC QC 指标

| 指标 | 含义 |
| --- | --- |
| `nCount_ATAC` | 每个细胞 peak counts |
| `nFeature_ATAC` | 每个细胞检测到的 peak 数 |
| `TSS.enrichment` | TSS enrichment，衡量开放染色质信号质量 |
| `FRiP` | fragments in peaks fraction |
| `nucleosome_signal` | nucleosome signal |
| `blacklist_fraction` | 落在 blacklist 区域的比例 |
| `fragments` | fragment 总量 |
| `unique_ratio` | 如果可用，表示 unique fragment/read 比例 |
| `scDblFinder.score` | doublet score |
| `scDblFinder.class` | singlet/doublet 分类 |

对应脚本：`scripts/process/process_single_sample.R`。

## 5. 过滤标准

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/process_single_sample.R` |
| 方法 | MAD outlier + doublet exclusion |
| 参数 | `--nmads`，当前默认 4 |
| 过滤字段 | `nCount_ATAC`, `TSS.enrichment`, `FRiP`, `scDblFinder.class` |
| 输出字段 | `pass_qc`, `outlier` |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合 QC 方法页 |

## 6. LSI 或其他降维

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/process_single_sample.R` |
| 方法 | ATAC query-native LSI/UMAP 流程，以及 CIMA reference-space LSI projection |
| 输出字段 | `umap_atac_1`, `umap_atac_2`, `cima_ref_umap_1`, `cima_ref_umap_2` |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合解释 ATAC LSI 和 UMAP |

## 7. 聚类

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/process_single_sample.R` |
| 输出字段 | `seurat_clusters` |
| 当前状态 | 已完成 |
| PPT 适用性 | 可和 UMAP/注释图一起展示 |

## 8. UMAP/t-SNE 等投影

| 项目 | 内容 |
| --- | --- |
| 当前投影 | UMAP |
| 输出图 | `umap_cima_cell_type_l1.png`, `umap_cima_cell_type_l2.png` |
| 坐标字段 | `umap_atac_1/2`, `cima_ref_umap_1/2` |
| t-SNE | 当前项目中未找到主线输出 |
| PPT 适用性 | 非常适合 |

## 9. gene activity 分析

| 项目 | 内容 |
| --- | --- |
| 当前状态 | 当前项目中未找到 gene activity matrix 主线输出 |
| 需要人工确认 | 是否存在于 notebook 或外部结果目录 |
| PPT 适用性 | 如果后续确认存在，可用于 RNA/ATAC 对照解释 |

## 10. motif 或 peak annotation

| 项目 | 内容 |
| --- | --- |
| motif | 当前项目中未找到 motif enrichment 主线输出 |
| peak annotation | 脚本中使用 EnsDb 和 peak 注释构建 ChromatinAssay，但未确认单独导出 peak annotation 文件 |
| 需要人工确认 | 是否有外部 motif/peak annotation 结果 |

## 11. cluster/cell type 注释

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/process/process_single_sample.R` |
| 方法 | CIMA ATAC compact feature model + centroid assignment |
| 注释层级 | L1, L2, L3, L4 |
| 输出字段 | `cima_cell_type_l1`, `cima_cell_type_l2`, `cima_cell_type_l3`, `cima_cell_type_l4` |
| 低置信字段 | `cima_cell_type_l1_masked`, `cima_l1_low_confidence`, `cima_l4_score`, `cima_l4_score_margin` |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合放 ATAC 注释结果页 |

## 12. 多样本合并或批处理

| 分支 | 命令 |
| --- | --- |
| only ATAC | `run-sample`, `run-gse`, `status` |
| co ATAC | `co-run-atac-sample`, `co-run-atac-gse`, `co-status` |
| longevity ATAC | `longevity-run-atac-sample`, `longevity-run-atac-all`, parameter contrast commands |
| product integration | `organize-products` |

当前状态：已完成。

## 13. 结果导出

典型输出：

| 输出 | 说明 |
| --- | --- |
| `metadata.csv` | 全细胞 metadata |
| `metadata_qc.csv` | pass-QC 细胞 metadata |
| `qc_summary.csv` | 样本级 QC 和注释摘要 |
| `validation_result.csv` | pass-QC 细胞关键字段和验证结果 |
| `qc_overview.png` | QC 总览图 |
| `umap_cima_cell_type_l1.png` | CIMA L1 UMAP |
| `umap_cima_cell_type_l2.png` | CIMA L2 UMAP |
| `{sample_id}.h5ad` | ATAC h5ad 输出 |
| `run_status.json` | 运行状态 |
| `logs/sample_qc.log` | 样本运行日志 |

## 14. 可视化

主要可视化输出：

- `qc_overview.png`
- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`
- product-level `figures/*_panels.png`

## 15. 最终结果和中间对象

建议后续结果审阅优先查看：

- `qc_summary.csv`
- `validation_result.csv`
- `metadata_qc.csv`
- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`
- `{sample_id}.h5ad`
- product-level manifests 和 integration metrics

大型 `.h5ad`、`.rds`、matrix、fragment 文件不建议上传 Git。
