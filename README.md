# ML2026 RNA/ATAC 多样本分析工作流

本仓库用于整理当前 ML2026 项目中已经完成或正在进行的 RNA/ATAC 多样本分析流程。目标不是重新设计新的生信流程，而是把已有代码、参数、输入输出规范、批处理逻辑、结果检查方法和问题排查经验整理成一个以后可以长期复用、维护和上传 Git 的项目。

本项目后续主要服务于三个对象：

- 给自己和其他人看的项目说明、流程文档和复用指南。
- 给 Git 仓库使用的代码、配置模板、文档和测试。
- 给 AI 使用的 workflow/skill，让 AI 后续能理解项目结构、检查结果和辅助排查问题。

## 项目范围

当前已确认的工作流包括：

- 单样本 scRNA-seq 处理：样本发现、计数矩阵读取、QC、doublet 处理、过滤、Scanpy 降维聚类、Azimuth `pbmcref` 注释、RNA 5 类 final celltype 映射、UMAP 可视化和结果导出。
- RNA baseline-only tuning：固定 `baseline__baseline__baseline` candidate 的审计型调参流程。
- 单样本 scATAC-seq 处理：fragment/barcode 输入、peak-by-cell matrix 构建、QC、doublet 处理、过滤、LSI/UMAP、CIMA ATAC 注释和结果导出。
- 共测 RNA/ATAC 样本处理：基于共测样本表组织 RNA 和 ATAC 分支运行。
- longevity 独立分支：processed RNA atlas ingest、ATAC barcode 预处理、ATAC 参数对照和结果发布。
- product-level 低维整合：整理已完成样本输出，生成 only RNA、only ATAC、co RNA、co ATAC 四类整合产品；该步骤只合并低维 embedding，不合并全量 count matrix。

当前项目中尚未确认属于主线产物的内容：

- RNA marker gene 表。
- ATAC gene activity matrix。
- ATAC motif enrichment 结果。
- t-SNE 结果。

如果这些内容存在于 notebook 或外部结果目录中，后续需要人工确认后再写入正式文档。

## 当前目录结构

```text
ML2026/
├── README.md
├── .gitignore
├── pyproject.toml
├── uv.lock
├── configs/
│   ├── config_template.yaml
│   └── sample_sheet_template.csv
├── scripts/
│   ├── only_rna/
│   ├── process/
│   ├── co/
│   └── longevity/
├── docs/
│   ├── project_summary.md
│   ├── io_summary.md
│   └── todo_and_uncertainties.md
├── references/
├── templates/
├── ai_skill/
├── results_examples/
├── notebooks/
└── tests/
```

当前阶段保留已有代码结构，不大规模移动 `scripts/` 下的核心脚本。第一步只新增项目说明、配置模板、Git 上传检查和文档骨架，避免破坏现有流程。

## 数据和结果管理原则

真实原始数据和大型结果不要提交到 Git。

以下内容应保留在本地或外部数据盘：

- 原始 count matrix、fragment 文件、FASTQ/BAM 文件和大型压缩表。
- `.h5ad`、`.rds`、`.RData`、`.loom`、完整 `matrix/` 目录和完整输出目录。
- 完整 `data/`、`output/`、`QualityControl/`、`ArchRLogs/` 内容，除非人工挑选极少量小型示例并确认可公开。
- 本地绝对路径配置、隐私 metadata、未脱敏样本信息。

Git 中应保存：

- 源代码。
- 测试。
- 配置模板。
- 文档。
- 小型示例说明。
- 结果索引。
- 问题排查流程。
- AI workflow/skill 草案。

## 环境依赖

Python 依赖目前写在 `pyproject.toml`，锁定文件为 `uv.lock`。项目使用 Python 3.11+。

主要 Python 包包括：

- `anndata`
- `scanpy`
- `pandas`
- `numpy`
- `scipy`
- `scikit-learn`
- `matplotlib`
- `bbknn`
- `harmonypy`
- `igraph`
- `leidenalg`
- `openpyxl`

ATAC 和 Azimuth 相关流程依赖 R 环境。当前 R 依赖尚未整理成独立 `environment.yml` 或安装说明，已从脚本中确认使用过：

- `Signac`
- `Seurat`
- `GenomeInfoDb`
- `EnsDb.Hsapiens.v86`
- `scDblFinder`
- `SingleCellExperiment`
- `Matrix`
- `rtracklayer`
- `ggplot2`
- `patchwork`
- `optparse`

R 环境仍需要后续人工确认和文档化。

## 快速开始

以下命令需要在仓库根目录运行，并要求数据根和参考文件已经准备好。

发现 RNA 样本：

```bash
uv run python scripts/process/pipeline.py discover-rna
```

运行单个 RNA 样本：

```bash
uv run python scripts/process/pipeline.py run-rna-sample --gse GSE_ID --sample-id SAMPLE_ID
```

运行一个 RNA GSE：

```bash
uv run python scripts/process/pipeline.py run-rna-gse --gse GSE_ID
```

查看 RNA 状态：

```bash
uv run python scripts/process/pipeline.py rna-status
```

运行一个共测 ATAC 样本：

```bash
uv run python scripts/process/pipeline.py co-run-atac-sample --gse DATASET_ID --gsm SAMPLE_ID
```

整理 product-level 整合输出：

```bash
uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink
```

更多输入输出路径说明见 `docs/io_summary.md`。

## 主要文档

- `docs/project_summary.md`：当前项目结构、流程状态和可复用代码总结。
- `docs/io_summary.md`：输入、输出、reference、raw data 和 result 目录说明。
- `docs/todo_and_uncertainties.md`：需要人工确认的内容和后续整理任务。
- `git_preparation_checklist.md`：上传 Git 前的检查清单。

后续会继续补充 RNA pipeline、ATAC pipeline、批处理、参数、结果索引、问题排查、PPT 提纲和 AI skill 文档。

## Git 上传前必须检查

上传 Git 前请至少检查：

- `git status --short`
- `git diff`
- 是否有原始数据被 stage。
- 是否有大型结果文件被 stage。
- 是否有本地绝对路径、隐私信息或未脱敏 metadata。
- notebook 是否包含大输出或本地路径。
- `.gitignore` 是否覆盖 `data/`、`output/`、`.h5ad`、`.rds`、fragment、matrix、log 等文件。

不要在未确认远程仓库和分支前执行 `git push`。
