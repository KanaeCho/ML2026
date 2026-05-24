# RNA Workflow SOP

本文件是 AI 执行 RNA 单样本和多样本处理时的操作参考。它基于当前项目已经实现的 `scripts/only_rna/` 主线，不引入未确认步骤。

## 适用输入

RNA 样本发现来自：

```text
data_root/reference/datasets.xlsx
```

支持的 raw 输入包括：

- 10x triplet：`matrix.mtx(.gz)`、`barcodes.tsv(.gz)`、`features.tsv(.gz)` 或 `genes.tsv(.gz)`。
- 10x `.h5`。
- `matrix.tar.gz`。
- 支持的 GSE-level shared triplet。

当前不把 shared gene-count CSV 当作单细胞矩阵主线输入。

## 发现和状态

```bash
uv run python scripts/process/pipeline.py discover-rna
uv run python scripts/process/pipeline.py discover-rna --gse <GSE>
uv run python scripts/process/pipeline.py rna-status
uv run python scripts/process/pipeline.py rna-status --gse <GSE>
```

发现结果中应确认：

- `supported` 是否为 true。
- `input_type` 是否符合实际文件。
- `matrix_path` / `h5_path` / `archive_path` 是否存在。
- GSE shared sample 是否需要通过 GSE-level 命令运行。

## 单样本运行

```bash
uv run python scripts/process/pipeline.py run-rna-sample --gse <GSE> --sample-id <SAMPLE>
```

常用控制参数：

- `--force`：覆盖已完成样本。
- `--dry-run`：只打印计划，不实际运行。
- `--output-root <PATH>`：覆盖默认 RNA 输出根。

## GSE 批处理

```bash
uv run python scripts/process/pipeline.py run-rna-gse --gse <GSE>
```

批处理策略：

1. 先 `discover-rna --gse <GSE>`。
2. 再 `rna-status --gse <GSE>`。
3. 先跑一个代表性 smoke sample。
4. smoke sample 输出完整后再运行 `run-rna-gse`。
5. 运行后再次 `rna-status --gse <GSE>`。

## RNA 主线步骤

1. 读取 count matrix 到 AnnData。
2. 写入样本 metadata：`sample`、`dataset`、`age`、`health`、`donor` 等。
3. 计算 QC 指标：`n_counts`、`n_genes`、`pct_mt`、`pct_ribo`。
4. 执行 doublet detection，生成 `doublet_score`、`is_doublet`。
5. 应用 dynamic hybrid MAD QC：
   - `n_counts`、`n_genes` 在 `log10(x + 1)` 空间做 lower-tail threshold。
   - `pct_mt`、`pct_ribo` 在原值空间做 upper-tail threshold。
   - `is_doublet == True` 不进入 `pass_qc`。
6. 只对 `pass_qc` 细胞执行 embedding：normalize_total、log1p、HVG、scale、PCA、neighbors、Leiden、UMAP。
7. 运行 Azimuth `pbmcref` 注释。
8. 将 Azimuth/CIMA L1 映射到 5 类 final celltype：`CD4_T`、`CD8_T`、`B`、`Myeloid`、`NK`。
9. 写出审计输出。

## RNA 输出契约

默认输出目录：

```text
output/rna/{GSE}/{sample_id}/
```

每个完成样本至少应包含：

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

`{sample_id}.h5ad` 只保留 pass-QC 且 final celltype 已知的细胞和表达矩阵。

## 关键 metadata 字段

RNA metadata 至少应关注：

- `pass_qc`
- `fails_count_floor`
- `fails_gene_floor`
- `fails_mt_ceiling`
- `fails_ribo_ceiling`
- `fails_doublet`
- `doublet_score`
- `is_doublet`
- `cluster`
- `umap_1`, `umap_2`
- `azimuth_cell_type`
- `azimuth_score`
- `azimuth_score_margin`
- `azimuth_low_confidence`
- `azimuth_cima_l1_raw`
- `azimuth_cima_l1`
- `final_celltype`
- `final_celltype_mapping`

## Tuning

当前 tuning 是 baseline-only，用于复核单 candidate 审计产物，不是大规模参数搜索。

```bash
uv run python scripts/process/pipeline.py tune-rna-sample --gse <GSE> --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py tune-rna-gse --gse <GSE>
```

输出位于：

```text
output/rna/{GSE}/{sample_id}/tuning/
```

## RNA 审计重点

运行后必须检查：

- `run_status.json`：命令、return code、`outputs_complete`。
- `qc_summary.csv`：总细胞数、pass-QC 数、final output 数、annotation status。
- `qc_thresholds.json`：动态阈值是否合理。
- UMAP 图：是否有严重塌缩、桥状结构或大面积 Unknown。
- `validation_result.csv`：输出完整性检查是否通过。

## 常见 RNA 问题

- 找不到样本：先看 `discover-rna` 是否支持该 input layout。
- GSE shared sample 不能用 `run-rna-sample`：改用 `run-rna-gse`。
- Azimuth 失败：检查 R、Seurat、Azimuth 环境和 `annotation_method_status`。
- pass-QC 过低：检查 `qc_thresholds.json` 和 `qc_overview.png`，不要立即放宽阈值。
- UMAP 不可读：先确认 pass-QC 细胞数、HVG/PC 数和日志，再考虑参数。
