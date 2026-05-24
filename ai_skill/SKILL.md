# Single-Cell RNA/ATAC Batch Workflow Skill

## 什么时候使用

当用户需要复用当前项目风格的单细胞 RNA / ATAC / RNA+ATAC 共测多样本处理流程时，使用本 skill。

本 skill 的目标是让 AI 可以按既有流程执行或指导：样本发现、输入检查、RNA QC、RNA 聚类和 UMAP、RNA 注释、ATAC barcode 选择、ATAC QC、ATAC LSI/UMAP、ATAC CIMA 投影注释、多样本批处理、输出审计和失败排查。

本 skill 不用于重新设计新的生信方法。除非用户明确要求，不要替换当前 QC、聚类、投影或注释策略。

## 核心入口

优先使用统一 CLI：

```bash
uv run python scripts/process/pipeline.py <command> [options]
```

核心代码入口：

- RNA 主线：`scripts/only_rna/`
- ATAC 主线：`scripts/process/process_single_sample.R`
- 共测分支：`scripts/co/cli.py`
- longevity 分支：`scripts/longevity/cli.py`
- product-level 整合：`scripts/process/organize_integrated_products.py`, `scripts/process/integrate_product_embeddings.py`

## AI 执行总原则

- 先发现样本，再运行样本。
- 先跑或检查一个 smoke sample，再批处理全量样本。
- 每个样本单独 QC、聚类、投影和注释；不要默认合并全量 count matrix 或 peak matrix。
- 多样本层面只在 product integration 阶段合并低维 embedding 和 metadata。
- 每次运行后必须检查 `run_status.json`、summary CSV、UMAP 图和 `.h5ad` 是否完整。
- 不擅自删除原始数据、旧结果或大型输出。
- 不擅自全量重跑；如果会重跑大量样本，先说明预计影响并询问用户。
- 不把 marker gene、gene activity、motif enrichment 写成主线步骤；当前主线中这些不是必需产物或未确认。

## 标准工作顺序

### 1. 判断数据通道

先确认用户要处理的是哪一类数据：

- only RNA：`datasets.xlsx` 中 scRNA 样本，对应 `run-rna-*`。
- only ATAC：`atac.xlsx` 中 ATAC fragment 样本，对应 `run-sample` / `run-gse`。
- co RNA/ATAC：`co.xlsx` 中 paired RNA/ATAC 样本，对应 `co-run-*`。
- longevity：独立 `data_root/raw/longevity/` 通道，对应 `longevity-*`。
- product integration：基于已完成样本输出组织 `output/1.only_atac/` 等 product。

### 2. 发现样本

按通道运行 discover，不要直接猜样本路径。

```bash
uv run python scripts/process/pipeline.py discover-rna
uv run python scripts/process/pipeline.py discover
uv run python scripts/process/pipeline.py co-discover
uv run python scripts/process/pipeline.py longevity-discover
```

如果只处理一个 GSE 或 dataset，优先加 `--gse <ID>` 或对应 sample 过滤参数。

### 3. 检查已有状态

```bash
uv run python scripts/process/pipeline.py rna-status
uv run python scripts/process/pipeline.py status
uv run python scripts/process/pipeline.py co-status
uv run python scripts/process/pipeline.py longevity-status
uv run python scripts/process/pipeline.py longevity-atac-barcode-status
```

状态检查重点：

- `status` 是否为 `success`。
- `outputs_complete` 是否为 true。
- 是否缺少 `.h5ad`、UMAP、summary 或 validation 文件。
- 日志是否存在错误。

### 4. Smoke sample

对新数据或新参数，先运行一个代表性样本。

RNA 示例：

```bash
uv run python scripts/process/pipeline.py run-rna-sample --gse <GSE> --sample-id <SAMPLE>
```

ATAC 示例：

```bash
uv run python scripts/process/pipeline.py run-sample --gse <GSE> --gsm <GSM>
```

co 示例：

```bash
uv run python scripts/process/pipeline.py co-run-rna-sample --gse <DATASET> --sample-id <SAMPLE>
uv run python scripts/process/pipeline.py co-run-atac-sample --gse <DATASET> --gsm <SAMPLE>
```

longevity ATAC 示例：

```bash
uv run python scripts/process/pipeline.py longevity-run-atac-sample --sample-id <SAMPLE>
```

### 5. 批处理

smoke sample 成功后再批处理。

```bash
uv run python scripts/process/pipeline.py run-rna-gse --gse <GSE>
uv run python scripts/process/pipeline.py run-gse --gse <GSE>
uv run python scripts/process/pipeline.py co-run-rna-gse --gse <DATASET>
uv run python scripts/process/pipeline.py co-run-atac-gse --gse <DATASET> --jobs <N>
uv run python scripts/process/pipeline.py longevity-run-atac-all
```

批处理前后都要重新运行 status。失败样本只针对失败样本排查，不要默认全量重跑。

### 6. Product-level 整合

样本级输出完成后再组织 product。

```bash
uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink
```

当前 product-level 整合只合并低维 embedding 和 metadata，不合并全量 RNA count matrix 或 ATAC peak matrix。RNA 默认 Harmony，ATAC 默认 BBKNN。只有用户明确要求时才使用 `--skip-integration`。

## RNA 主线 SOP

参考：`ai_skill/references/rna_workflow.md`。

RNA 单样本流程：

1. 从 `datasets.xlsx` 发现已筛选 scRNA 样本。
2. 读取 triplet、10x h5、`matrix.tar.gz` 或 supported shared triplet。
3. 计算 `n_counts`、`n_genes`、`pct_mt`、`pct_ribo`。
4. 执行 doublet detection。
5. 用 dynamic hybrid MAD 生成 `pass_qc`。
6. 对 pass-QC 细胞执行 normalize、log1p、HVG、PCA、neighbors、Leiden、UMAP。
7. 执行 Azimuth `pbmcref` 注释。
8. 映射到 5 类 final celltype：`CD4_T`、`CD8_T`、`B`、`Myeloid`、`NK`。
9. 写出 metadata、QC summary、validation、UMAP、`.h5ad` 和 status。

## ATAC 主线 SOP

参考：`ai_skill/references/atac_workflow.md`。

ATAC 单样本流程：

1. 定位 `fragments.tsv.gz`。
2. 按优先级确定 barcode：显式 `--barcode-file`、filtered metadata、singlecell metadata、fragment-count inference。
3. 用 `peak.bed` 构建 peak-by-cell matrix。
4. 计算 fragments、TSS enrichment、FRiP、nucleosome、blacklist、peak counts 等 QC。
5. 执行 scDblFinder doublet detection。
6. 用 MAD outlier + doublet exclusion 生成 `pass_qc`。
7. 执行 TF-IDF、SVD/LSI、neighbors、Leiden、UMAP。
8. 投影到 CIMA ATAC reference compact feature space。
9. 写出 CIMA L1-L4 注释、metadata、QC summary、validation、UMAP、`.h5ad` 和 status。

## 批处理和审计 SOP

参考：`ai_skill/references/batch_workflow.md` 和 `ai_skill/references/output_contracts.md`。

每批样本完成后必须汇总：

- 样本总数、成功数、失败数、跳过数。
- 每个样本 `n_cells_total`、`n_cells_pass_qc` 或 ATAC pass-QC 数。
- annotation status 和低置信比例。
- UMAP 图是否存在。
- `.h5ad` 是否存在。
- `validation_result.csv` 是否通过输出完整性检查。
- 失败样本对应 `logs/sample_qc.log` 的错误阶段。

## 失败排查 SOP

参考：`ai_skill/references/troubleshooting.md`。

遇到失败或结果异常时按顺序定位：

1. 通道：RNA、ATAC、co、longevity、product integration。
2. 阶段：发现、读取、QC、doublet、embedding、注释、导出、绘图、整合。
3. 输入：样本表、raw 文件、reference 文件、barcode 文件是否匹配。
4. 输出：status、summary、validation、UMAP、`.h5ad` 是否缺失。
5. 日志：`logs/sample_qc.log` 或命令输出中第一个实际错误。
6. 影响范围：单样本、单 GSE、单通道还是 product-level。
7. 修复策略：优先小范围重跑，避免直接全量重跑。

## 输出不能进 Git

不要把以下文件加入 Git：

- `data/` 真实内容。
- `output/` 真实内容。
- `.h5ad`、`.rds`、`.RData`、`.loom`。
- fragment、FASTQ、BAM、matrix、压缩 metadata 大文件。
- `QualityControl/`、`ArchRLogs/` 全量结果。
- 本地绝对路径配置和隐私 metadata。

## 如何复用到新项目

如果要把本 skill 用到另一个 RNA/ATAC 项目：

1. 保留本 skill 的流程结构。
2. 替换项目专用路径、样本表、reference 资产和输出目录。
3. 先用新项目代码验证 discover/status/smoke sample 命令。
4. 删除或标注不适用的 ML2026 特例，例如特定 GSE 规则、CIMA 资产路径或 longevity 分支。
5. 不要直接照搬已完成样本数、真实输出路径或本地 fallback 路径作为新项目事实。
