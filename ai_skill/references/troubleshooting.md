# Troubleshooting SOP

本文件用于指导 AI 排查 RNA/ATAC 批处理失败或结果异常。

## 基本原则

- 不要一上来建议全量重跑。
- 先定位通道、样本、阶段、输入、输出和日志。
- 区分“已出现问题”“潜在风险”“当前项目中未找到”。
- 只在有证据时下结论。

## 定位问题阶段

先判断错误属于哪一类：

- 样本发现失败。
- 输入文件缺失或格式不匹配。
- reference 资产缺失。
- QC 指标计算失败。
- doublet 检测失败。
- embedding / clustering / UMAP 失败。
- Azimuth 或 CIMA 注释失败。
- 输出导出失败。
- product integration 失败。

## 必查文件

单样本：

```text
run_status.json
logs/sample_qc.log
qc_summary.csv
validation_result.csv
metadata.csv
metadata_qc.csv
```

product：

```text
product_status.json
manifests/samples.csv
qc/sample_qc_summary.csv
qc/validation_summary.csv
integration/integration_summary.json
integration/integration_metrics.csv
```

## RNA 常见问题

### 找不到样本

检查：

- `discover-rna` 输出。
- `datasets.xlsx` 行是否可见且 assay 为 scRNA。
- raw 文件是否符合支持 layout。
- GSE shared sample 是否应通过 `run-rna-gse`。

### QC 后细胞过少

检查：

- `qc_thresholds.json`。
- `qc_overview.png`。
- `fails_count_floor`、`fails_gene_floor`、`fails_mt_ceiling`、`fails_ribo_ceiling`、`fails_doublet` 比例。
- doublet score 是否异常。

不要直接放宽阈值；先确认样本质量和输入是否正确。

### Azimuth 失败

检查：

- R 是否可调用。
- Seurat / Azimuth 是否安装。
- `annotation_method_status`。
- `logs/sample_qc.log` 中 Azimuth 阶段报错。

### UMAP 不可读

检查：

- pass-QC 细胞数。
- HVG 数和 PC 数是否足够。
- 是否 fallback 到 minimal clustering。
- annotation low-confidence fraction。

## ATAC 常见问题

### barcode 过多

检查：

- 是否缺失 `filtered_barcodes.tsv.gz`。
- 是否误用 fragment-count inference。
- co/longevity wrapper 是否传入了 `--barcode-file`。
- `summary.json` 中 selected barcode 数。

### barcode 过少

检查：

- fragment 文件是否正确。
- `min_fragments`、`min_tss`、`max_barcodes`。
- `barcode_qc.csv.gz`。
- ArchR preprocessing 是否成功。

### FRiP/TSS 低

检查：

- sample raw quality。
- barcode selection。
- peak reference 是否匹配。
- `qc_overview.png`。

不要直接把低质量样本纳入正式结果。

### CIMA 低置信高

检查：

- QC 后细胞质量。
- CIMA reference 资产。
- 样本是否是 reference 覆盖不佳的组织或状态。
- cluster purity 和 score margin。

## Product integration 常见问题

### 样本未进入 product

检查：

- 样本 `.h5ad` 是否存在。
- `validation_result.csv` 是否存在。
- output root 是否正确。
- product organizer 是否排除了该分支或 dataset。

### 低置信 RNA 细胞没有 integrated UMAP

这是当前 only_rna product 的预期行为：低于默认 `integrated_cima_l1_score` 阈值的细胞保留 metadata，但不参与 product-level UMAP/Leiden。

### panel 图缺失

检查：

- 是否使用了 `--skip-figures`。
- `integrated_umap_1/2` 是否存在。
- `product_status.json` 中是否记录跳过原因。

## 何时需要询问用户

以下情况先问用户：

- 需要删除或覆盖旧结果。
- 需要全量重跑。
- 需要改变 QC 阈值或 barcode 阈值。
- 需要把某些低质量样本纳入正式结果。
- 需要修改主线注释策略。
