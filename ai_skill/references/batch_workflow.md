# Batch Processing SOP

本文件描述 AI 如何用当前项目流程做多样本批处理和审计。

## 批处理前检查

1. 明确通道：only RNA、only ATAC、co、longevity、product integration。
2. 确认 data root：优先 `./data`，其次 `ML2026_DATA_ROOT`，最后才是项目本机 fallback。
3. 确认样本表或 raw layout 存在。
4. 运行 discover 命令。
5. 运行 status 命令。
6. 选择一个 smoke sample。
7. 确认用户是否允许批量运行和覆盖旧结果。

## Discover 命令

```bash
uv run python scripts/process/pipeline.py discover-rna
uv run python scripts/process/pipeline.py discover
uv run python scripts/process/pipeline.py co-discover
uv run python scripts/process/pipeline.py longevity-discover
```

## Status 命令

```bash
uv run python scripts/process/pipeline.py rna-status
uv run python scripts/process/pipeline.py status
uv run python scripts/process/pipeline.py co-status
uv run python scripts/process/pipeline.py longevity-status
uv run python scripts/process/pipeline.py longevity-atac-barcode-status
```

## 单样本 smoke test

RNA：

```bash
uv run python scripts/process/pipeline.py run-rna-sample --gse <GSE> --sample-id <SAMPLE>
```

only ATAC：

```bash
uv run python scripts/process/pipeline.py run-sample --gse <GSE> --gsm <GSM>
```

co：

```bash
uv run python scripts/process/pipeline.py co-run-rna-sample --gse <DATASET> --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py co-run-atac-sample --gse <DATASET> --gsm <SAMPLE>
```

longevity：

```bash
uv run python scripts/process/pipeline.py longevity-ingest-rna --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py longevity-run-atac-sample --sample-id <SAMPLE>
```

## 批量运行

RNA GSE：

```bash
uv run python scripts/process/pipeline.py run-rna-gse --gse <GSE>
```

only ATAC GSE：

```bash
uv run python scripts/process/pipeline.py run-gse --gse <GSE>
```

co RNA/ATAC：

```bash
uv run python scripts/process/pipeline.py co-run-rna-gse --gse <DATASET>
uv run python scripts/process/pipeline.py co-run-atac-gse --gse <DATASET> --jobs <N>
```

longevity ATAC：

```bash
uv run python scripts/process/pipeline.py longevity-run-atac-all
```

## 参数对照

当前 longevity ATAC 支持 barcode 参数对照：

```bash
uv run python scripts/process/pipeline.py longevity-run-atac-param-contrast --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py longevity-run-atac-custom-param-contrast --sample-id <SAMPLE> --min-fragments <N> --max-barcodes <N> --min-tss <X>
uv run python scripts/process/pipeline.py longevity-summarize-atac-param-contrast --write
uv run python scripts/process/pipeline.py longevity-publish-atac-param-contrast --dry-run
```

发布参数对照结果到正式目录前必须确认选择表和目标输出。

## Product-level 整合

```bash
uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink
```

可选产品：

- `only_atac`
- `only_rna`
- `co_atac`
- `co_rna`
- `all`

默认执行低维整合。只有在生成 metadata-only 检查时才使用 `--skip-integration`。

## 批处理后审计

每批完成后汇总：

- 样本总数。
- `success`、`failed`、`dry_run`、skipped 数量。
- 缺失的 expected output。
- RNA `n_cells_total`、`n_cells_pass_qc`、`n_cells_final_output`。
- ATAC pass-QC 细胞数、median FRiP、median TSS、CIMA low-confidence fraction。
- 每个失败样本的日志路径和首个错误阶段。
- product-level `integration_summary.json` 和 `integration_metrics.csv`。

## 失败处理策略

- 单样本失败：只重跑该样本。
- 同一 GSE 多样本失败：先找共同输入或 reference 问题。
- 所有样本失败：优先检查环境、data root、reference 文件和 CLI 参数。
- UMAP 或注释质量差：不要直接全量重跑，先检查 QC summary、低置信比例和样本来源。
