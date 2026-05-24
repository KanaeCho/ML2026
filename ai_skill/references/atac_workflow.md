# ATAC Workflow SOP

本文件是 AI 执行 ATAC 单样本、co-ATAC 和 longevity ATAC 处理时的操作参考。ATAC 主线脚本是 `scripts/process/process_single_sample.R`。

## 适用输入

ATAC 主输入是：

```text
fragments.tsv.gz
```

可选辅助输入：

- `filtered_barcodes.tsv.gz`
- `filtered_metadata.csv.gz`
- `singlecell.csv.gz`
- `data_root/reference/peak.bed`
- CIMA ATAC reference assets

barcode 优先级：

1. 显式 `--barcode-file`。
2. 自动发现或显式 `filtered_metadata`。
3. 自动发现或显式 `singlecell` metadata。
4. 从 fragment counts 推断 barcode。

## only ATAC 发现、状态和运行

```bash
uv run python scripts/process/pipeline.py discover
uv run python scripts/process/pipeline.py discover --gse <GSE>
uv run python scripts/process/pipeline.py status
uv run python scripts/process/pipeline.py status --gse <GSE>
uv run python scripts/process/pipeline.py run-sample --gse <GSE> --gsm <GSM>
uv run python scripts/process/pipeline.py run-gse --gse <GSE>
```

only ATAC 样本表来自：

```text
data_root/reference/atac.xlsx
```

## co-ATAC 发现、状态和运行

```bash
uv run python scripts/process/pipeline.py co-discover
uv run python scripts/process/pipeline.py co-status
uv run python scripts/process/pipeline.py co-run-atac-sample --gse <DATASET> --gsm <SAMPLE>
uv run python scripts/process/pipeline.py co-run-atac-gse --gse <DATASET> --jobs <N>
```

co 样本表来自：

```text
data_root/reference/co.xlsx
```

raw layout：

```text
data_root/raw/{dataset}/{sample}/ATAC/fragments.tsv.gz
```

如果同目录存在 barcode 辅助文件，co 包装层会显式传给 ATAC 主线，避免误用过宽的 fragment-count inference。

## longevity ATAC 发现、barcode 和运行

```bash
uv run python scripts/process/pipeline.py longevity-discover
uv run python scripts/process/pipeline.py longevity-atac-barcode-status
uv run python scripts/process/pipeline.py longevity-preprocess-atac-barcodes --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py longevity-run-atac-sample --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py longevity-run-atac-all
```

longevity raw layout：

```text
data_root/raw/longevity/atac/*_fragments.tsv.gz
```

barcode 预处理输出：

```text
data_root/reference/longevity/atac_barcodes/{sample_id}/filtered_barcodes.tsv.gz
data_root/reference/longevity/atac_barcodes/{sample_id}/barcode_qc.csv.gz
data_root/reference/longevity/atac_barcodes/{sample_id}/summary.json
```

默认 barcode 参数：

- `barcode_min_fragments=200`
- `barcode_max_barcodes=20000`
- `barcode_min_tss=2.5`
- `barcode_rank_by=fragments`

## ATAC 主线步骤

1. 定位 fragment 文件。
2. 确定初始 barcode 集合。
3. 用 `peak.bed` 构建 peak-by-cell matrix。
4. 构建 ChromatinAssay / Seurat object。
5. 计算 QC 指标：`nCount_ATAC`、`nFeature_ATAC`、`TSS.enrichment`、`FRiP`、`nucleosome_signal`、`blacklist_fraction`、`fragments`。
6. 执行 scDblFinder doublet detection，生成 `scDblFinder.score` 和 `scDblFinder.class`。
7. 应用 MAD outlier + doublet exclusion，生成 `pass_qc`。
8. 对 pass-QC 细胞执行 TF-IDF、SVD/LSI、neighbors、Leiden、UMAP。
9. 投影到 CIMA ATAC compact feature model。
10. 写入 CIMA L1-L4 注释和低置信审计字段。
11. 写出 metadata、QC summary、validation、UMAP、`.h5ad` 和 status。

## ATAC 输出契约

ATAC 输出可能位于：

```text
output/atac/{GSE}/{GSM}/
output/co/atac/{dataset}/{sample}/
output/longevity/atac/longevity/{sample_id}/
```

每个完成样本至少应包含：

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

正式样本输出当前不依赖保留 `matrix/` 目录或 RDS 文件；这些大型中间对象不要提交 Git。

## 关键 metadata 字段

ATAC metadata 至少应关注：

- `pass_qc`
- `outlier`
- `nCount_ATAC`
- `nFeature_ATAC`
- `TSS.enrichment`
- `FRiP`
- `nucleosome_signal`
- `blacklist_fraction`
- `fragments`
- `scDblFinder.score`
- `scDblFinder.class`
- `seurat_clusters`
- `umap_atac_1`, `umap_atac_2`
- `cima_ref_umap_1`, `cima_ref_umap_2`
- `cima_cell_type_l1`
- `cima_cell_type_l1_masked`
- `cima_cell_type_l2`, `cima_cell_type_l3`, `cima_cell_type_l4`
- `cima_l1_low_confidence`
- `cima_l4_score`
- `cima_l4_score_margin`
- `final_celltype`

## ATAC 审计重点

运行后必须检查：

- `run_status.json`：return code 和 `outputs_complete`。
- `qc_summary.csv`：pass-QC 数、FRiP、TSS、低置信比例。
- `validation_result.csv`：关键字段和输出完整性。
- `qc_overview.png`：QC 分布和阈值。
- `umap_cima_cell_type_l1.png`、`umap_cima_cell_type_l2.png`：注释是否可读。
- `.h5ad` 是否存在且能被 product integration 读取。

## 常见 ATAC 问题

- barcode 数过多：检查是否缺少 `filtered_barcodes.tsv.gz`，是否回退到了 fragment-count inference。
- barcode 数过少：检查 fragment 文件、min fragments、TSS 阈值和 barcode 预处理 summary。
- FRiP/TSS 低：先检查样本质量和 barcode 选择，不要直接放宽 QC。
- UMAP 失败：检查 pass-QC 细胞数、peak 数、LSI 日志和 R 包环境。
- CIMA 低置信高：检查 reference projection、样本组织来源和 QC 质量。
