# 多样本批处理文档

本文档整理当前项目中的多样本批处理逻辑，包括样本来源、输入路径、输出结构、状态记录、跳过已完成样本和失败定位方式。

## 1. 样本列表从哪里来

| 分支 | 样本来源 | 读取脚本 |
| --- | --- | --- |
| only RNA | `data_root/reference/datasets.xlsx`，可见且 `assay=scRNA` 的行 | `scripts/only_rna/discovery.py`, `scripts/process/pipeline.py` |
| only ATAC | `data_root/reference/atac.xlsx` | `scripts/process/pipeline.py` |
| co RNA/ATAC | `data_root/reference/co.xlsx` | `scripts/co/cli.py` |
| longevity RNA | `data_root/raw/longevity/rna/*.h5ad` | `scripts/longevity/cli.py` |
| longevity ATAC | `data_root/raw/longevity/atac/*_fragments.tsv.gz` | `scripts/longevity/cli.py` |
| product integration | 已完成的样本输出目录 | `scripts/process/organize_integrated_products.py` |

## 2. 每个样本的输入路径如何定义

RNA only 分支：

```text
data_root/raw/{GSE}/{GSM...}_matrix.mtx(.gz)
data_root/raw/{GSE}/{GSM...}_barcodes.tsv(.gz)
data_root/raw/{GSE}/{GSM...}_features.tsv(.gz)
```

ATAC only 分支：

```text
data_root/raw/{GSE}/{GSM...fragments...}.tsv.gz
```

co 分支：

```text
data_root/raw/{dataset}/{sample}/RNA/
data_root/raw/{dataset}/{sample}/ATAC/
```

longevity 分支：

```text
data_root/raw/longevity/rna/*.h5ad
data_root/raw/longevity/atac/*_fragments.tsv.gz
```

## 3. 每个样本生成哪些结果

RNA 样本典型输出：

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

ATAC 样本典型输出：

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

product 输出：

```text
product_status.json
manifests/samples.csv
manifests/cells_metadata.csv
manifests/output_files.csv
qc/sample_qc_summary.csv
qc/validation_summary.csv
integration/integration_summary.json
integration/integration_metrics.csv
integration/sample_mixing_summary.csv
figures/*_panels.png
```

## 4. RNA 和 ATAC 是否使用相同 ID 或 barcode

当前结论：不完全相同，需要按分支区分。

| 分支 | ID 规则 |
| --- | --- |
| only RNA | 通常为 `GSE/sample_id`，多数 sample_id 为 GSM，也可能是 GSE shared sample |
| only ATAC | 通常为 `GSE/GSM` |
| co | 使用 `co.xlsx` 中的英文 `sample`，例如 `donorA_Day0` |
| longevity | 使用独立 sample_id，例如 W 开头样本或 processed atlas 名称 |

当前项目不执行 co RNA+ATAC paired-cell matching，因此不要默认假设 RNA barcode 和 ATAC barcode 一一对应。

## 5. 是否存在硬编码路径

已确认存在部分本机路径或 dataset 特例：

- `/mnt/g/ML2026_data`：当前 data root fallback。
- `/mnt/g/ML2026_output`：当前 output 挂载说明。
- `/home/linuxbrew/.linuxbrew/lib`：longevity 代码中的系统库路径。
- `GSE226039` PBMC-only 规则。
- `GSE198533` shared gene-count CSV unsupported 规则。
- product integration 中排除 `GSE206284`、`GSE282769` 的规则。

这些暂时不建议直接改动，应该先记录在文档和参数表中，后续再考虑配置化。

## 6. 是否有日志记录

有。

每个样本通常有：

```text
run_status.json
logs/sample_qc.log
```

`run_status.json` 通常记录：

- 样本 ID。
- 执行命令。
- 输出根目录。
- 开始和结束时间。
- return code。
- `outputs_complete`。
- `status`。

## 7. 是否能单独重跑某个样本

支持。

示例：

```bash
uv run python scripts/process/pipeline.py run-rna-sample --gse GSE_ID --sample-id SAMPLE_ID
uv run python scripts/process/pipeline.py co-run-atac-sample --gse DATASET_ID --gsm SAMPLE_ID
uv run python scripts/process/pipeline.py longevity-run-atac-sample --sample-id SAMPLE_ID
```

## 8. 是否能跳过已完成样本

支持。

当前多个入口都会根据 expected outputs 判断样本是否完成。如果已完成且没有使用 `--force`，则跳过。

常见参数：

- `--force`：强制重跑。
- `--dry-run`：只显示命令，不实际运行。
- `--skip-complete`：部分 longevity parameter contrast 命令支持。

## 9. 批处理失败时如何定位

排查顺序：

1. 查看命令行输出。
2. 查看样本目录下的 `run_status.json`。
3. 查看样本目录下的 `logs/sample_qc.log`。
4. 检查 expected output 是否缺失。
5. 根据失败阶段查看对应脚本：RNA 看 `scripts/only_rna/`，ATAC 看 `process_single_sample.R`。
6. 检查输入路径和 metadata/sample ID/barcode 是否一致。

## 10. 是否适合整理成通用 workflow

适合。

当前批处理已经具备：

- 样本发现。
- 单样本运行。
- GSE/批量运行。
- 输出完整性检查。
- 状态记录。
- 日志记录。
- product-level 汇总。

后续需要进一步配置化的部分：

- dataset 特例。
- data_root/output_root。
- ATAC barcode 参数。
- integration 参数。
- R 环境依赖。
