# Git 上传准备检查清单

这个文件用于在提交或上传 Git 仓库前检查项目是否安全、清晰、可复用。不要在完成高优先级检查前执行 `git push`。

## 1. 确认仓库范围

- [ ] 确认 Git 仓库只包含代码、测试、文档、配置模板、AI workflow 和少量示例说明。
- [ ] 确认原始数据和完整输出目录不进入 Git。
- [ ] 确认目标仓库是公开仓库还是私有仓库。
- [ ] 确认 `AGENTS.md` 是保留为内部说明，还是拆分成公开文档。
- [x] 默认按公开仓库标准处理：`.planning/`、`.opencode/`、`opencode.json` 和未清理 notebook 暂不提交。

## 2. 检查工作区状态

提交前先运行：

```bash
git status --short
git diff
```

需要确认：

- [ ] 没有原始数据被 stage。
- [ ] 没有大型结果文件被 stage。
- [ ] 没有本地私有配置被 stage。
- [ ] 没有缓存文件被误加入。
- [ ] 所有准备提交的文件都是有意加入的。

## 3. 不应该上传到普通 Git 的文件

以下内容不要作为普通 Git 文件提交：

- [ ] `data/` 真实内容。
- [ ] `output/` 真实内容。
- [ ] 完整 `QualityControl/` 内容。
- [ ] 完整 `ArchRLogs/` 内容。
- [ ] `.h5ad` 文件。
- [ ] `.rds`、`.RData`、`.rda` 文件。
- [ ] `.loom` 文件。
- [ ] `*.bam`、`*.bai`、`*.cram` 文件。
- [ ] `*.fastq`、`*.fastq.gz`、`*.fq.gz` 文件。
- [ ] `*_fragments.tsv.gz` 文件。
- [ ] 完整 `matrix/` 输出目录。
- [ ] 大型 `*.mtx`、`*.mtx.gz`、`*.csv.gz`、`*.tsv.gz` 文件。
- [ ] 完整 per-sample 日志，除非是人工挑选的排查案例。
- [ ] 包含本地路径的 local config。
- [ ] notebook checkpoint 目录。
- [ ] 未清理的 `*.ipynb` notebook。
- [ ] `.planning/` 本地规划/debug 文件，除非人工筛选后确认可公开。
- [ ] `.opencode/` 和 `opencode.json` 本地工具配置，除非团队明确需要共享。

## 4. `.gitignore` 检查

当前 `.gitignore` 已经排除 `data` 和 `output`，并补充了大型数据和结果文件类型。上传前仍需确认这些规则是否满足当前仓库需求：

```gitignore
QualityControl/
ArchRLogs/
*.h5ad
*.rds
*.RData
*.rda
*.loom
*.bam
*.bai
*.cram
*.fastq
*.fastq.gz
*.fq
*.fq.gz
*_fragments.tsv.gz
*.mtx
*.mtx.gz
*.csv.gz
*.tsv.gz
*.log
logs/
local_config*.yaml
.env
.env.*
.claude_resources.json
.planning/
.opencode/
opencode.json
*.ipynb
.ipynb_checkpoints/
```

如果确实需要提交小型示例图或小型示例表，请放入 `results_examples/`，并先人工检查内容是否可公开。

## 5. 本地路径检查

公开上传前建议搜索本地绝对路径：

```bash
rg "/mnt/|/home/|C:|D:|F:|G:" .
```

当前已知路径风险：

- [ ] `/mnt/g/ML2026_data`：当前数据根 fallback。
- [ ] `/mnt/g/ML2026_output`：当前 output 挂载说明。
- [x] `/mnt/f/...`：`tests/only_rna/test_tuning.py` 中的本地测试路径已移除；仍需人工确认 `.planning/` 和本地工具文件是否上传。
- [ ] `/home/linuxbrew/.linuxbrew/lib`：longevity 代码中的系统库路径。

如果这些路径是代码中的可覆盖 fallback，可以保留并在文档中说明；如果出现在模板或公开文档中，应替换成占位符。

## 6. Notebook 检查

提交 notebook 前需要：

- [ ] 清理大型输出。
- [ ] 删除本地绝对路径。
- [ ] 删除不应公开的样本信息或 metadata。
- [ ] 确认 notebook 是正式流程、探索记录还是历史版本。
- [ ] 根据用途移动到 `notebooks/rna/`、`notebooks/atac/` 或 `notebooks/exploration/`。

当前已发现但默认不上传的 notebook：

- [ ] `scripts/scATAC.ipynb`
- [ ] `scripts/step01_signac_single_sample_covid-redo.ipynb`
- [ ] `scripts/process/GSM8671454.ipynb`

这些 notebook 需要清空输出、删除本地路径和敏感 metadata 后，再移动到 `notebooks/rna/`、`notebooks/atac/` 或 `notebooks/exploration/` 并显式取消 ignore。

## 7. 环境依赖检查

Python：

- [ ] 确认 `pyproject.toml` 里的依赖足够支撑当前文档化流程。
- [ ] 确认 `uv.lock` 是否随仓库提交。
- [ ] 确认 CUDA 相关 Torch index 是否适合目标使用者。

R：

- [ ] 记录 R 版本。
- [ ] 记录 R 包版本。
- [ ] 补充 Bioconductor 包安装说明。
- [ ] 后续考虑添加 `environment.yml` 或 `docs/environment.md`。

脚本中已确认使用的 R 包包括：

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

## 8. 文档完整性检查

第一批已完成：

- [x] `README.md`
- [x] `configs/config_template.yaml`
- [x] `configs/sample_sheet_template.csv`
- [x] `docs/project_summary.md`
- [x] `docs/io_summary.md`
- [x] `docs/todo_and_uncertainties.md`

后续仍需人工复核或继续补充：

- [x] `docs/rna_pipeline.md`
- [x] `docs/atac_pipeline.md`
- [x] `docs/batch_processing.md`
- [x] `docs/results_index.md`
- [x] `docs/parameters.md`
- [x] `docs/problem_index.md`
- [x] `docs/troubleshooting.md`
- [x] `docs/ppt_knowledge_outline.md`
- [x] `references/knowledge_map.md`
- [x] `references/qc_metrics.md`
- [x] `references/annotation_notes.md`
- [x] `references/diagnostic_workflow.md`
- [x] `templates/case_note_template.md`
- [x] `templates/ppt_outline_template.md`
- [x] `templates/report_template.md`
- [x] `ai_skill/SKILL.md`

## 9. 测试和验证

代码改动提交前建议运行轻量测试：

```bash
uv run pytest tests/test_rna_pipeline.py tests/test_integrated_products.py tests/test_longevity_atac_barcode_preprocessing.py
```

Git 准备阶段不要为了验证文档而重跑大规模数据处理。

运行测试前需要确认：

- [ ] 测试不依赖私有外部数据。
- [ ] 测试不会写入大型输出。
- [ ] 测试里的本地绝对路径已修复或已记录为待处理。

## 10. 推荐 commit 分组

建议按小提交组织，不要一次性提交所有内容：

1. `init reusable project documentation structure`
2. `add RNA and ATAC pipeline documentation`
3. `add batch processing config templates`
4. `add troubleshooting and diagnostic workflow docs`
5. `add AI workflow skill draft`
6. `add PPT and report templates`
7. `update gitignore for data and result artifacts`
8. `document runtime environment requirements`

## 11. push 前最终检查

- [ ] `git status --short` 只包含预期文件。
- [ ] `git diff --cached` 已完整审阅。
- [ ] `.gitignore` 已排除大型数据和结果文件。
- [ ] README 说明数据和输出都在 Git 外部。
- [ ] 模板中没有私有本地路径。
- [ ] 没有 secrets 或隐私 metadata 被 stage。
- [ ] 没有原始数据或完整输出文件被 stage。
- [ ] 没有大型二进制文件被 stage，除非明确使用 Git LFS 或 DVC 管理。
- [ ] commit message 按用途分组。
- [ ] 远程仓库和目标分支确认后，再执行 `git push`。
