# ML2026 Agent Guide

## 当前分支定位

当前分支：`fix`。

本分支不是 `main`，修改本文件只影响 `fix` 分支。不要在未得到用户确认前切换、合并或推送到 `main`。

本分支有两个目标：

- 维护 ML2026 当前已实现的单细胞 RNA/ATAC 多样本分析流程。
- 将当前 RNA/ATAC QC、聚类、投影、注释、批处理、审计流程整理成可复用的 AI workflow skill。

## AI Workflow Skill 入口

可复用 workflow skill 位于：

```text
ai_skill/SKILL.md
```

当用户目标是“以后遇到类似 scRNA-seq / scATAC-seq / RNA+ATAC 共测数据时，让 AI 直接照着流程处理”时，AI 应优先读取：

```text
ai_skill/SKILL.md
ai_skill/references/rna_workflow.md
ai_skill/references/atac_workflow.md
ai_skill/references/batch_workflow.md
ai_skill/references/output_contracts.md
ai_skill/references/command_reference.md
ai_skill/references/troubleshooting.md
ai_skill/references/reuse_guide.md
```

`AGENTS.md` 记录 ML2026 当前项目事实和操作边界；`ai_skill/` 记录可复用 SOP。不要把完整 SOP 重复写进 `AGENTS.md`。

复用到新项目时，只复用 workflow 顺序和审计逻辑，不要直接照搬 ML2026 的本地路径、已完成样本数、GSE 特例或真实输出状态。

## 项目事实

当前项目聚焦单细胞 RNA/ATAC 多样本处理，主要通道包括：

- only RNA
- only ATAC
- RNA+ATAC 共测 co
- longevity 独立通道
- product-level low-dimensional integration

数据根解析优先级：

1. 仓库根目录下的 `./data`
2. 环境变量 `ML2026_DATA_ROOT`
3. 本机 fallback `/mnt/g/ML2026_data`

输出默认位于：

```text
output/
```

真实数据、完整输出和大型中间对象不属于 Git 长期追踪对象。

## 主要入口

统一 CLI：

```text
scripts/process/pipeline.py
```

常用调用方式：

```bash
uv run python scripts/process/pipeline.py <command> [options]
```

核心代码入口：

- RNA 主线：`scripts/only_rna/`
- ATAC 主线：`scripts/process/process_single_sample.R`
- co 分支：`scripts/co/cli.py`
- longevity 分支：`scripts/longevity/cli.py`
- product organization：`scripts/process/organize_integrated_products.py`
- product embedding integration：`scripts/process/integrate_product_embeddings.py`
- product panel rendering：`scripts/process/render_product_umap_panels.py`

脚本目录说明见：

```text
scripts/README.md
```

## 主线流程摘要

### RNA

RNA 主线是 Python-first 单样本流程，支持按 GSE 批处理。

流程摘要：

1. 从 `data_root/reference/datasets.xlsx` 发现可见且 `assay=scRNA` 的样本。
2. 读取 triplet、10x h5、`matrix.tar.gz` 或 supported shared triplet。
3. 计算 `n_counts`、`n_genes`、`pct_mt`、`pct_ribo`。
4. 执行 doublet detection。
5. 使用 dynamic hybrid MAD 生成 `pass_qc`。
6. 对 pass-QC 细胞执行 Scanpy normalization、HVG、PCA、neighbors、Leiden、UMAP。
7. 执行 Azimuth `pbmcref` 注释。
8. 映射到 5 类 `final_celltype`：`CD4_T`、`CD8_T`、`B`、`Myeloid`、`NK`。
9. 输出 metadata、QC、validation、UMAP、`.h5ad`、`run_status.json`。

RNA 详细 SOP 见：

```text
ai_skill/references/rna_workflow.md
docs/rna_pipeline.md
```

### ATAC

ATAC 主线复用 `scripts/process/process_single_sample.R`，only ATAC、co ATAC 和 longevity ATAC 都围绕该主线或包装逻辑执行。

流程摘要：

1. 定位 `fragments.tsv.gz`。
2. 按优先级确定 barcode：显式 barcode、filtered metadata、singlecell metadata、fragment-count inference。
3. 使用 `peak.bed` 构建 peak-by-cell matrix。
4. 计算 TSS、FRiP、fragment、nucleosome、blacklist 等 QC 指标。
5. 执行 scDblFinder doublet detection。
6. 使用 MAD outlier + doublet exclusion 生成 `pass_qc`。
7. 执行 TF-IDF、SVD/LSI、neighbors、Leiden、UMAP。
8. 投影到 CIMA ATAC compact feature model 并生成 L1-L4 注释。
9. 输出 metadata、QC、validation、UMAP、`.h5ad`、`run_status.json`。

ATAC 详细 SOP 见：

```text
ai_skill/references/atac_workflow.md
docs/atac_pipeline.md
```

### 批处理

批处理标准顺序：

1. `discover`
2. `status`
3. smoke sample
4. batch run
5. audit
6. product integration

批处理 SOP 和命令速查见：

```text
ai_skill/references/batch_workflow.md
ai_skill/references/command_reference.md
docs/batch_processing.md
```

### Product-Level Integration

Product integration 只合并低维 embedding 和 metadata，不合并全量 RNA count matrix 或 ATAC peak matrix。

默认输出：

```text
output/1.only_atac/
output/2.only_rna/
output/3.co_atac/
output/4.co_rna/
```

默认命令：

```bash
uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink
```

默认不要使用 `--skip-integration`，除非用户只想生成 metadata-level product。

## 当前输出契约

RNA 样本完成后至少应有：

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

ATAC 样本完成后至少应有：

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

完整 output contract 见：

```text
ai_skill/references/output_contracts.md
docs/io_summary.md
```

## 禁止或需确认操作

未经用户确认，不要：

- 删除原始数据、本地结果或旧输出。
- 全量重跑大批量样本。
- 修改核心 QC、barcode、annotation 或 integration 方法。
- 把 marker gene、gene activity、motif enrichment 写成当前主线必需步骤。
- 上传真实数据、大型结果、中间对象或隐私 metadata。
- push 到 `main`。
- merge `fix` 到 `main`。
- 修改 GitHub 默认分支。
- 把 ML2026 的本地路径、GSE 特例、完成样本数当作新项目事实照搬。

如果需要覆盖旧结果、全量重跑、改变阈值、发布参数对照结果或推送远程，必须先向用户确认。

## Git 上传边界

不要上传或追踪：

- `data/` 真实内容
- `output/` 真实内容
- `QualityControl/` 全量内容
- `ArchRLogs/` 全量内容
- `.h5ad`
- `.rds`
- `.RData`
- `.rda`
- `.loom`
- fragment 文件
- FASTQ/BAM/CRAM 文件
- matrix 大文件
- `.csv.gz` / `.tsv.gz` 大型压缩表
- `.env` 和本地配置
- `.opencode/`
- `.planning/`
- notebook
- 本地学习笔记

提交前至少检查：

```bash
git status --short
git diff --check
git ls-files | rg -i '\.(h5ad|rds|RData|rda|loom|bam|bai|cram|fastq|fq|mtx|ipynb)$|fragments\.tsv\.gz$|\.(csv|tsv)\.gz$|(^|/)data/|(^|/)output/'
```

## 文档维护规则

如果修改以下内容，需要同步更新对应文档：

- RNA/ATAC 输入发现规则
- QC 阈值或过滤逻辑
- barcode 选择逻辑
- annotation 层级或 final celltype 映射
- 输出文件契约
- CLI 命令接口
- product integration 行为
- `ai_skill/` workflow SOP

项目事实更新写入 `AGENTS.md` 或 `docs/`；可复用执行流程更新写入 `ai_skill/`。

## 当前分支工作建议

在 `fix` 分支继续开发时，优先保持：

- `AGENTS.md` 简洁记录项目事实和 AI 入口。
- `ai_skill/` 记录可复用 RNA/ATAC workflow SOP。
- `docs/` 记录给人看的项目说明。
- 代码变更按功能分组 commit。
- push 前由用户确认。
