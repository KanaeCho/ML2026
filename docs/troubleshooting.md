# Troubleshooting 问题排查手册

本文档给出常见 RNA/ATAC 项目问题的排查入口。它不是让你直接重跑，而是告诉你先看哪些文件、按什么顺序判断问题来源。

## 快速定位

| 现象 | 优先检查 |
| --- | --- |
| 样本运行失败 | `run_status.json`, `logs/sample_qc.log` |
| 输出文件缺失 | expected output 列表、`validation_result.csv`, `run_status.json` |
| QC 后细胞太少 | `qc_summary.csv`, `qc_thresholds.json`, `metadata.csv` |
| RNA 注释不好 | `qc_summary.csv`, `metadata_qc.csv`, `umap_rna_*.png`, `annotation_method_status` |
| ATAC 注释不好 | `qc_summary.csv`, `validation_result.csv`, `umap_cima_cell_type_*.png` |
| UMAP 混乱 | UMAP 图、metadata 中 sample/GSE/donor/cluster/label 字段 |
| product 整合异常 | `integration/integration_summary.json`, `integration/integration_metrics.csv`, `sample_mixing_summary.csv` |
| 版本混乱 | `run_status.json`, 输出目录时间、tuning 目录、product manifest |

## RNA 常见问题

### RNA QC 异常

优先检查：

- `metadata.csv`
- `metadata_qc.csv`
- `qc_summary.csv`
- `qc_thresholds.json`
- `qc_overview.png`

检查顺序：

1. 查看 `n_cells_total` 和 `n_cells_pass_qc`。
2. 查看 `pass_qc_fraction`。
3. 查看 `final_min_counts`, `final_min_genes`, `final_max_pct_mt`, `final_max_pct_ribo`。
4. 查看 `qc_thresholds.json` 中 guardrail 是否触发。
5. 查看 doublet 比例。
6. 判断是样本质量问题还是阈值问题。

### RNA 注释图不好

优先检查：

- `umap_rna_pbmcref_vs_cima_l1.png`
- `umap_rna_pbmcref_highlight.png`
- `umap_rna_cima_l1.png`
- `metadata_qc.csv`
- `qc_summary.csv`
- `validation_result.csv`

检查顺序：

1. `azimuth_status` 是否为 ok。
2. `azimuth_score_mean` 是否偏低。
3. `azimuth_low_confidence_fraction` 是否偏高。
4. `final_celltype` 是否大量 Unknown 或被过滤。
5. 绘图字段是否是最新的 `azimuth_cima_l1` / `final_celltype`。
6. 判断是注释质量问题、label 映射问题还是绘图字段问题。

## ATAC 常见问题

### ATAC QC 异常

优先检查：

- `qc_summary.csv`
- `metadata.csv`
- `metadata_qc.csv`
- `qc_overview.png`
- `logs/sample_qc.log`

检查顺序：

1. 查看 input cells 和 pass QC cells。
2. 查看 `median_TSS_enrichment`。
3. 查看 `median_FRiP`。
4. 查看 `median_fragments`。
5. 查看 doublet 数量。
6. 查看 barcode_source 是 filtered barcode 还是 fragment inference。
7. 如果 barcode 来自 inference，再检查 barcode threshold 和 rank plot。

### ATAC fragment 或 peak matrix 异常

优先检查：

- fragment 文件路径。
- tabix index 是否存在或能生成。
- `peak.bed` 是否存在。
- `logs/sample_qc.log` 中 FeatureMatrix 相关错误。
- `matrix/matrix.mtx` 是否生成。

不要直接重跑全体样本，先对单样本确认 fragment 和 peak reference 是否正确。

### ATAC 注释不好

优先检查：

- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`
- `validation_result.csv`
- `qc_summary.csv`
- `metadata_qc.csv`

检查顺序：

1. 查看 `cima_l1_low_confidence_frac`。
2. 查看 `median_cima_l4_score`。
3. 查看 L1/L2 label 是否大量 Unknown。
4. 查看 UMAP basis 是 reference projected 还是 query native。
5. 判断是 QC 问题、reference projection 问题还是 label 解释问题。

## product-level 整合常见问题

### batch effect 明显或样本分开

优先检查：

- `manifests/cells_metadata.csv`
- `integration/integration_summary.json`
- `integration/integration_metrics.csv`
- `integration/sample_mixing_summary.csv`
- `figures/*_gse_panels.png`
- `figures/*_sample_panels.png`

检查顺序：

1. 确认使用的是 `integrated_umap_1/2`，不是 per-sample UMAP。
2. 查看 integration method 是 Harmony、BBKNN 还是 scanpy neighbors。
3. 查看 batch key。
4. 查看每个样本细胞数是否极端不平衡。
5. 查看低置信 RNA 细胞是否被排除出整合。
6. 判断是生物差异、batch effect 还是输入样本组成差异。

### 多样本合并后某个样本主导结果

优先检查：

- `manifests/samples.csv`
- `manifests/cells_metadata.csv`
- `qc/sample_qc_summary.csv`
- `sample_mixing_summary.csv`

检查顺序：

1. 统计每个样本细胞数。
2. 查看大样本是否在 UMAP 上占据主要区域。
3. 查看样本是否来自不同组织或条件。
4. 检查是否应分层展示而不是强行合并解释。

## 批处理失败

优先检查：

- 失败样本的 `run_status.json`。
- 失败样本的 `logs/sample_qc.log`。
- 输入文件路径。
- 样本是否在 manifest 中。
- expected output 是否缺失。

常见来源：

- 输入文件缺失。
- reference 文件缺失。
- R/Python 环境缺依赖。
- barcode 文件格式不对。
- fragment index 缺失。
- Azimuth/Rscript 调用失败。
- 输出目录权限或挂载问题。

## 结果版本混乱

优先检查：

- sample-root 输出目录。
- tuning candidate 嵌套目录。
- product manifest 中的 `source_output_dir`。
- `run_status.json` 中的命令和时间。
- 是否存在旧 legacy 输出。

处理建议：

1. 不要直接删除旧结果。
2. 先记录哪个目录是当前主线输出。
3. 在 `docs/results_index.md` 中标注版本关系。
4. 如需清理，先征得确认。
