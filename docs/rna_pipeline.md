# RNA Pipeline 文档

本文档整理当前项目中已经实现的 RNA 分析流程。内容只基于当前仓库中已找到的脚本、配置和项目说明，不把未确认步骤写成已完成。

## 总览

当前 RNA 主线位于：

```text
scripts/only_rna/
```

主要入口通过：

```text
scripts/process/pipeline.py
```

当前主线是 Python-first 单样本 RNA 流程，支持按 GSE 批处理，并支持 baseline-only tuning。

## 1. 输入数据

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/discovery.py`, `scripts/only_rna/read_inputs.py`, `scripts/process/pipeline.py` |
| 输入来源 | `data_root/reference/datasets.xlsx` 中可见且 `assay=scRNA` 的样本，以及 `data_root/raw/{GSE}/` 下本地已下载矩阵 |
| 支持格式 | matrix triplet、10x `.h5`、`matrix.tar.gz`、部分 GSE shared triplet |
| 当前状态 | 已完成 |
| 是否支持批处理 | 支持 `run-rna-gse` 和 `tune-rna-gse` |
| PPT 适用性 | 可用于“数据来源和样本结构”页 |
| 不确定项 | `GSE198533` shared gene-count CSV 当前不支持；marker 相关输入未确认 |

支持的输入路径模式见 `docs/io_summary.md`。

## 2. 样本信息和 metadata

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/discovery.py`, `scripts/only_rna/outputs.py`, `scripts/co/cli.py` |
| 输入 metadata | `datasets.xlsx`；共测分支为 `co.xlsx` |
| 关键字段 | `gse`, `sample_id`, `input_type`, `sample_kind`, `age`, `health`, `donor`, `individual_id` |
| 输出 metadata | `metadata.csv`, `metadata_qc.csv`, `.h5ad.obs` |
| 当前状态 | 已完成 |
| 是否支持批处理 | 支持 |
| PPT 适用性 | 可用于样本表和流程图 |
| 不确定项 | 不同分支的 sample ID/barcode 是否完全一致需要具体样本确认 |

## 3. QC 指标

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/qc.py` |
| 指标 | `n_counts`, `n_genes`, `pct_mt`, `pct_ribo` |
| doublet 指标 | `doublet_score`, `is_doublet`，由 `scripts/only_rna/doublet.py` 处理 |
| 输出文件 | `metadata.csv`, `metadata_qc.csv`, `qc_summary.csv`, `qc_thresholds.json`, `qc_overview.png` |
| 当前状态 | 已完成 |
| 是否支持批处理 | 支持，每个样本单独计算 |
| PPT 适用性 | 适合放 QC 指标解释和样本 QC 图 |

## 4. 过滤标准

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/qc.py`, `scripts/only_rna/default_config.yaml` |
| 方法 | dynamic hybrid MAD |
| lower-tail 指标 | `n_counts`, `n_genes`，在 `log10(x + 1)` 空间计算动态阈值 |
| upper-tail 指标 | `pct_mt`, `pct_ribo`，在原值空间计算动态阈值 |
| doublet 过滤 | `is_doublet == True` 不进入 `pass_qc` |
| 审计输出 | `qc_thresholds.json` 记录阈值、MAD、guardrail 和最终阈值 |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合放 QC 过滤逻辑页 |

当前默认 QC 参数见 `docs/parameters.md`。

## 5. 标准化和归一化

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/embedding.py` |
| 上游依赖 | `pass_qc == True` 的细胞 |
| 方法 | `scanpy.pp.normalize_total(target_sum=1e4)` 和 `scanpy.pp.log1p` |
| 当前状态 | 已完成 |
| 是否支持批处理 | 支持，每个样本单独执行 |
| PPT 适用性 | 适合放 RNA 流程图，不一定单独成页 |

## 6. 高变基因选择

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/embedding.py` |
| 方法 | `scanpy.pp.highly_variable_genes(flavor='seurat')` |
| 默认数量 | 根据样本规模自动选择；大样本默认最多 1000，小样本默认最多 2000，受基因数限制 |
| 当前状态 | 已完成 |
| PPT 适用性 | 可放在方法流程页 |

## 7. PCA 和降维

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/embedding.py` |
| 方法 | `scanpy.pp.scale(max_value=10)`, `scanpy.tl.pca` |
| 默认 PC | 配置默认为 30，并受细胞数和基因数限制 |
| 当前状态 | 已完成 |
| PPT 适用性 | 可用于说明 PCA、neighbors、UMAP 的关系 |

## 8. 聚类

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/embedding.py` |
| 方法 | `scanpy.pp.neighbors`, `scanpy.tl.leiden` |
| 输出字段 | `cluster` |
| fallback | 如果环境缺少 `igraph`，会回退到最小 deterministic fallback |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合和 UMAP/注释图一起展示 |

## 9. UMAP/t-SNE 等投影

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/embedding.py`, `scripts/only_rna/plotting.py`, `scripts/only_rna/outputs.py` |
| 当前投影 | UMAP |
| 输出字段 | `umap_1`, `umap_2` |
| 输出图 | `umap_rna_pbmcref_vs_cima_l1.png`, `umap_rna_pbmcref_highlight.png`, `umap_rna_cima_l1.png` |
| t-SNE | 当前项目中未找到主线输出 |
| PPT 适用性 | 非常适合 |

## 10. marker gene 分析

| 项目 | 内容 |
| --- | --- |
| 当前状态 | 当前项目中未找到主线 marker gene 表或 marker 绘图输出 |
| 是否支持批处理 | 未确认 |
| PPT 适用性 | 如果后续找到 marker 结果，适合放入注释验证页 |
| 需要人工确认 | 是否存在于 notebook 或外部结果目录 |

## 11. cluster/cell type 注释

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/annotation.py`, `scripts/only_rna/azimuth.py`, `scripts/only_rna/final_celltype.py`, `scripts/only_rna/outputs.py` |
| 当前主线方法 | Azimuth `pbmcref` |
| 输出字段 | `azimuth_cell_type`, `azimuth_cell_type_l2_raw`, `azimuth_cima_l1_raw`, `azimuth_cima_l1`, `final_celltype` |
| final celltype | `CD4_T`, `CD8_T`, `B`, `Myeloid`, `NK` |
| 审计字段 | `azimuth_score`, `azimuth_score_margin`, `azimuth_low_confidence`, `annotation_method_status` |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合放注释策略和结果图 |

注意：`annotation.py` 中保留一些 alternative annotation helper，但当前主线默认是 Azimuth `pbmcref`。

## 12. 多样本合并或批处理

| 项目 | 内容 |
| --- | --- |
| 对应脚本 | `scripts/only_rna/cli.py`, `scripts/process/pipeline.py`, `scripts/process/organize_integrated_products.py`, `scripts/process/integrate_product_embeddings.py` |
| 单样本命令 | `run-rna-sample` |
| GSE 批处理命令 | `run-rna-gse` |
| tuning 命令 | `tune-rna-sample`, `tune-rna-gse` |
| 状态命令 | `rna-status` |
| product integration | `organize-products` |
| 当前状态 | 已完成 |
| PPT 适用性 | 适合放多样本流程图 |

## 13. 结果导出

| 输出 | 说明 |
| --- | --- |
| `metadata.csv` | 全细胞 metadata 审计表 |
| `metadata_qc.csv` | QC 后且 final celltype 已知的细胞表 |
| `qc_summary.csv` | 样本级 QC 和注释摘要 |
| `qc_thresholds.json` | QC 动态阈值审计 |
| `validation_result.csv` | 输出完整性和注释状态检查 |
| `{sample_id}.h5ad` | pass-QC 且 final celltype 已知的细胞矩阵和 obs |
| `run_status.json` | 运行状态、命令、输出完整性 |
| UMAP PNG | RNA 注释和高亮图 |

## 14. 可视化

主要可视化输出：

- `qc_overview.png`
- `umap_rna_pbmcref_vs_cima_l1.png`
- `umap_rna_pbmcref_highlight.png`
- `umap_rna_cima_l1.png`
- tuning overview UMAPs
- product-level panel figures

## 15. 最终结果和中间对象

最终结果建议以以下文件为主：

- `metadata_qc.csv`
- `qc_summary.csv`
- `validation_result.csv`
- `{sample_id}.h5ad`
- UMAP PNG
- product-level `manifests/cells_metadata.csv`
- product-level `integration/integration_metrics.csv`
- product-level `figures/*_panels.png`

大型 `.h5ad` 不建议上传 Git，只在文档中索引路径。
