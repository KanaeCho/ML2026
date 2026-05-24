# 项目现状总结

本文档基于当前工作区已经审计到的仓库文件，整理 ML2026 RNA/ATAC 多样本分析项目的现状。这里只记录当前项目中有证据支持的内容。未找到或未核实的内容会标注为“当前项目中未找到”或“未确认”。

## 当前项目目的

本项目正在被整理成一个可复用的 Git 项目，用于后续处理类似 RNA/ATAC 多样本数据。当前目标不是重新设计生信流程，而是把已经实现的分析逻辑、输入输出关系、参数、批处理命令、结果检查和问题排查经验保存下来。

主要目标：

- 复用当前已有 RNA/ATAC 处理代码。
- 将代码、配置模板、文档和 troubleshooting 说明纳入 Git。
- 将原始数据、大型中间对象和完整结果目录排除在 Git 之外。
- 同时提供给人看的文档和给 AI 使用的 workflow/skill。
- 沉淀 QC、聚类、UMAP、注释、batch effect、批处理失败等异常结果的排查经验。

## 当前顶层目录

```text
ML2026/
├── AGENTS.md
├── README.md
├── .gitignore
├── pyproject.toml
├── uv.lock
├── configs/
├── scripts/
│   ├── only_rna/
│   ├── process/
│   ├── co/
│   └── longevity/
├── tests/
├── docs/
├── references/
├── templates/
├── ai_skill/
├── results_examples/
├── notebooks/
├── QualityControl/
├── ArchRLogs/
├── data
└── output
```

## 当前主要代码模块

| 模块 | 路径 | 当前作用 |
| --- | --- | --- |
| RNA 主线 | `scripts/only_rna/` | Python-first RNA 单样本流程，包括发现、读取、QC、doublet、embedding、注释、输出和 tuning。 |
| 总 CLI | `scripts/process/pipeline.py` | 项目主命令入口，负责 RNA、ATAC、co、longevity、download、status 和 product organization 路由。 |
| ATAC 主线 | `scripts/process/process_single_sample.R` | R/Signac ATAC 单样本处理脚本，被 only ATAC、co ATAC 和 longevity ATAC 复用。 |
| product-level 整合 | `scripts/process/organize_integrated_products.py` 和 `scripts/process/integrate_product_embeddings.py` | 将已完成样本输出整理成四类 product，并执行低维整合。 |
| 共测分支 | `scripts/co/cli.py` | 从 `co.xlsx` 发现并运行共测 RNA/ATAC 样本。 |
| longevity 分支 | `scripts/longevity/cli.py` | 处理 longevity RNA atlas ingest、ATAC barcode 预处理和参数对照。 |
| 测试 | `tests/` | 测试 RNA discovery、输出契约、product organization 和 longevity barcode preprocessing。 |

## RNA 流程当前状态

RNA 工作流已经实现，并且可通过 CLI 运行。

已确认步骤：

1. 从 `datasets.xlsx` 的可见 scRNA 行和本地 raw layout 发现样本。
2. 读取 matrix triplet、10x `.h5`、matrix archive 和部分 GSE-level shared triplet。
3. 计算 QC 指标：`n_counts`、`n_genes`、`pct_mt`、`pct_ribo`。
4. 执行 doublet 处理。
5. 使用 dynamic hybrid MAD 和 guardrail 进行 QC 过滤。
6. 对 pass-QC 细胞执行 Scanpy normalization、log transform、HVG、PCA、neighbors、Leiden 和 UMAP。
7. 使用 Azimuth `pbmcref` 作为当前主线注释方法。
8. 将 RNA 结果映射到 `CD4_T`、`CD8_T`、`B`、`Myeloid`、`NK` 五类 final celltype。
9. 导出 sample-level metadata、QC summary、validation、h5ad 和 UMAP 图。
10. 支持 baseline-only tuning，candidate 为 `baseline__baseline__baseline`。
11. 支持 product-level RNA 低维整合，使用 CIMA compact feature projection 和 Harmony。

当前项目中未找到或未确认：

- RNA marker gene 主线输出。
- t-SNE 主线输出。

## ATAC 流程当前状态

ATAC 工作流已经通过 R/Signac 脚本实现，并由 Python CLI 包装成可批处理流程。

已确认步骤：

1. 根据分支从 `atac.xlsx`、`co.xlsx` 或 longevity raw fragment 目录发现 ATAC 样本。
2. 支持默认 GSE/GSM fragment 文件，也支持显式 `--fragment-file`。
3. barcode 来源优先级为 filtered barcode、filtered metadata、singlecell metadata、fragment-count inference。
4. 使用 `peak.bed` 和 Signac `FeatureMatrix` 构建 peak-by-cell matrix。
5. 计算 ATAC QC 指标，包括 count、feature、TSS enrichment、FRiP、nucleosome signal、blacklist fraction、fragments 和可选 unique ratio。
6. 使用 `scDblFinder` 进行 doublet 检测。
7. 使用 MAD outlier 和 doublet exclusion 进行 QC 过滤。
8. 在 ATAC 脚本中执行 LSI、UMAP 和 clustering。
9. 使用 CIMA ATAC compact feature model 和 centroids 进行 L1-L4 注释。
10. 导出 metadata、QC summary、validation、UMAP 图、matrix、部分 profile 下的 RDS、h5ad 和 run status。
11. 支持 product-level ATAC 低维整合，使用 CIMA ATAC LSI compact feature projection 和 BBKNN。

当前项目中未找到或未确认：

- ATAC gene activity matrix 主线输出。
- motif enrichment 主线输出。
- 独立 peak annotation 结果文件。当前脚本中存在用于构建 ChromatinAssay 的基因注释，但未确认有单独导出结果。

## 多样本批处理当前状态

多样本批处理已经实现。

| 分支 | 样本来源 | 主要命令 | 状态 |
| --- | --- | --- | --- |
| only RNA | `data_root/reference/datasets.xlsx` | `discover-rna`、`run-rna-sample`、`run-rna-gse`、`tune-rna-*`、`rna-status` | 已实现 |
| only ATAC | `data_root/reference/atac.xlsx` | `status`、`run-sample`、`run-gse` | 已实现 |
| co RNA | `data_root/reference/co.xlsx` | `co-discover`、`co-run-rna-sample`、`co-run-rna-gse` | 已实现 |
| co ATAC | `data_root/reference/co.xlsx` | `co-run-atac-sample`、`co-run-atac-gse` | 已实现 |
| longevity RNA | `data_root/raw/longevity/rna/*.h5ad` | `longevity-ingest-rna` | 已实现，但属于特殊分支 |
| longevity ATAC | `data_root/raw/longevity/atac/*_fragments.tsv.gz` | `longevity-run-atac-*`、`longevity-preprocess-atac-barcodes`、参数对照命令 | 已实现，但属于特殊分支 |
| integrated products | 已完成的样本输出目录 | `organize-products` | 已实现 |

代码中已确认的批处理能力：

- 支持单样本重跑。
- 支持 GSE 或全样本运行。
- 已完成输出会被跳过，除非使用 `--force`。
- 多个命令支持 `--dry-run`。
- 每个样本写 `run_status.json`。
- 每个样本写 `logs/sample_qc.log`。
- 根据 expected output contract 检查输出完整性。

## 可以直接复用的代码

当前可直接复用的主要代码：

- `scripts/process/pipeline.py`
- `scripts/only_rna/`
- `scripts/process/process_single_sample.R`
- `scripts/co/cli.py`
- `scripts/longevity/cli.py`
- `scripts/process/organize_integrated_products.py`
- `scripts/process/integrate_product_embeddings.py`
- `scripts/process/render_product_umap_panels.py`
- `tests/test_rna_pipeline.py`
- `tests/test_integrated_products.py`
- `tests/test_longevity_atac_barcode_preprocessing.py`

## 后续需要整理或配置化的地方

这些不是当前文档整理阶段必须立即修改的内容，但应该在后续 Git 项目化过程中持续跟踪。

| 问题 | 证据 | 建议 |
| --- | --- | --- |
| 数据根 fallback 是本机路径 | Python/R 代码中存在 `/mnt/g/ML2026_data` | 当前可保留作为本机 fallback，但文档和模板中应使用占位符。 |
| dataset 特例写在代码里 | 例如 `GSE226039`、`GSE198533`、`GSE206284`、`GSE282769` | 后续可整理成 dataset policy 表或 YAML。 |
| RNA discovery 逻辑不止一处 | `scripts/only_rna/discovery.py` 和 `scripts/process/pipeline.py` | 后续让 `scripts/only_rna/discovery.py` 成为唯一 canonical 实现。 |
| ATAC R 脚本较长 | `scripts/process/process_single_sample.R` | 当前先不拆，后续在测试充分后再按 input/QC/annotation/export 拆分。 |
| R 环境未统一记录 | R 包直接在脚本里 `library(...)` | 后续补充 R 环境说明或 `environment.yml`。 |
| 测试中存在本地绝对路径 | `tests/only_rna/test_tuning.py` | Git 发布前应改成相对路径或临时目录。 |
| notebook 状态不清楚 | 多个 notebook 仍在 `scripts/` 下 | 需要人工判断哪些是正式、探索或旧版本。 |

## Git 上传适配性总结

适合上传 Git：

- 源代码。
- 测试。
- 配置模板。
- 文档。
- 人工筛选后的小型示例图或结果说明。
- AI workflow/skill 说明。

不适合普通 Git 上传：

- 原始数据。
- 完整输出目录。
- 大型 `.h5ad`、`.rds`、`.RData`、`.loom`、matrix、fragment、FASTQ、BAM 文件。
- 完整 ArchR 日志。
- 本地路径配置。
- 隐私 metadata 或未确认可公开的样本表。

## 当前文档整理阶段

第一批项目化文件已包括：

- `README.md`
- `configs/config_template.yaml`
- `configs/sample_sheet_template.csv`
- `docs/project_summary.md`
- `docs/io_summary.md`
- `docs/todo_and_uncertainties.md`
- `git_preparation_checklist.md`

后续计划继续补充 RNA、ATAC、批处理、参数、结果索引、问题排查、PPT 提纲和 AI skill 文档。
