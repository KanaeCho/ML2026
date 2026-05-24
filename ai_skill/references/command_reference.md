# Command Reference

统一入口：

```bash
uv run python scripts/process/pipeline.py <command> [options]
```

## RNA

```bash
uv run python scripts/process/pipeline.py discover-rna [--gse <GSE>]
uv run python scripts/process/pipeline.py rna-status [--gse <GSE>]
uv run python scripts/process/pipeline.py run-rna-sample --gse <GSE> --sample-id <SAMPLE> [--force] [--dry-run]
uv run python scripts/process/pipeline.py run-rna-gse --gse <GSE> [--force] [--dry-run]
uv run python scripts/process/pipeline.py tune-rna-sample --gse <GSE> --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py tune-rna-gse --gse <GSE>
```

## only ATAC

```bash
uv run python scripts/process/pipeline.py discover [--gse <GSE>]
uv run python scripts/process/pipeline.py status [--gse <GSE>]
uv run python scripts/process/pipeline.py run-sample --gse <GSE> --gsm <GSM> [--nmads 4] [--force] [--dry-run]
uv run python scripts/process/pipeline.py run-gse --gse <GSE> [--nmads 4] [--force] [--dry-run]
```

## co RNA/ATAC

```bash
uv run python scripts/process/pipeline.py co-discover [--gse <DATASET>]
uv run python scripts/process/pipeline.py co-status [--gse <DATASET>]
uv run python scripts/process/pipeline.py co-run-rna-sample --gse <DATASET> --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py co-run-rna-gse --gse <DATASET>
uv run python scripts/process/pipeline.py co-run-atac-sample --gse <DATASET> --gsm <SAMPLE>
uv run python scripts/process/pipeline.py co-run-atac-gse --gse <DATASET> --jobs <N>
```

## longevity

```bash
uv run python scripts/process/pipeline.py longevity-discover [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-status [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-ingest-rna [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-atac-barcode-status [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-preprocess-atac-barcodes [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-run-atac-sample --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py longevity-run-atac-all
uv run python scripts/process/pipeline.py longevity-run-atac-param-contrast [--sample-id <SAMPLE>]
uv run python scripts/process/pipeline.py longevity-run-atac-custom-param-contrast --sample-id <SAMPLE> --min-fragments <N> --max-barcodes <N> --min-tss <X>
uv run python scripts/process/pipeline.py longevity-summarize-atac-param-contrast --write
uv run python scripts/process/pipeline.py longevity-publish-atac-param-contrast --dry-run
```

## product integration

```bash
uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink
uv run python scripts/process/pipeline.py organize-products --products only_rna --copy-mode symlink
uv run python scripts/process/pipeline.py organize-products --products only_atac --copy-mode symlink
uv run python scripts/process/pipeline.py organize-products --products co_rna --copy-mode symlink
uv run python scripts/process/pipeline.py organize-products --products co_atac --copy-mode symlink
```

默认不要使用 `--skip-integration`，除非用户只想生成 metadata-level product。
