# RNA/ATAC 多样本项目整理与复用 Skill

## 用途

当用户需要整理、复用、维护或排查当前 ML2026 风格 RNA/ATAC 多样本分析项目时，使用本 skill。

本 skill 的目标不是重新设计生信流程，而是帮助 AI 理解已有项目结构、已有 RNA/ATAC 工作流、输入输出关系、批处理方式、Git 上传边界、结果排查逻辑和文档整理方式。

## 适用场景

当用户提出以下需求时，应使用本 skill：

- 整理 RNA/ATAC 多样本分析项目。
- 复用已有 RNA/ATAC 批处理流程。
- 检查 RNA/ATAC 输出是否完整。
- 回顾 QC、过滤、降维、聚类、UMAP、注释逻辑。
- 整理 Git 上传前的代码和文档。
- 整理 README、docs、references、templates。
- 生成 PPT 提纲或论文方法整理。
- 排查 QC 异常、聚类异常、UMAP 混乱、注释不清楚、batch effect 或批处理失败。
- 将一次问题处理过程记录成 case note。

## 项目结构读取顺序

AI 处理该项目时应优先读取：

1. `README.md`
2. `docs/project_summary.md`
3. `docs/io_summary.md`
4. `docs/rna_pipeline.md`
5. `docs/atac_pipeline.md`
6. `docs/batch_processing.md`
7. `docs/parameters.md`
8. `docs/results_index.md`
9. `docs/troubleshooting.md`
10. `references/diagnostic_workflow.md`
11. `references/knowledge_map.md`
12. `git_preparation_checklist.md`

必要时再读取代码：

- RNA：`scripts/only_rna/`
- ATAC：`scripts/process/process_single_sample.R`
- 总入口：`scripts/process/pipeline.py`
- co：`scripts/co/cli.py`
- longevity：`scripts/longevity/cli.py`
- product integration：`scripts/process/organize_integrated_products.py`, `scripts/process/integrate_product_embeddings.py`

## 使用原则

- 优先整理已有内容，不擅自新增不存在的分析步骤。
- 不擅自替换 RNA/ATAC 方法。
- 不重跑大规模计算，除非用户明确要求。
- 不删除已有文件。
- 不移动或重命名核心脚本，除非用户确认。
- 不把潜在问题写成当前项目已经发生的问题。
- 对不确定内容标注“未确认”或“当前项目中未找到”。
- 结论尽量追溯到脚本、配置、metadata、日志或结果文件。
- Git 上传前必须保护原始数据、大型结果和本地路径。

## 输入

本 skill 可使用的输入包括：

- 原始数据路径说明。
- 样本信息表。
- metadata。
- RNA/ATAC 脚本。
- notebook。
- config 文件。
- QC 阈值。
- 聚类和降维参数。
- 注释规则。
- marker gene 表，如果项目中存在。
- 已有结果目录。
- 日志文件。
- `.h5ad`、`.rds` 等中间对象路径。
- UMAP/QC/annotation 图。

## 输出

本 skill 可帮助生成或维护：

- 项目结构说明。
- RNA pipeline 文档。
- ATAC pipeline 文档。
- 批处理说明。
- 参数表。
- 输入输出索引。
- 结果文件索引。
- PPT 知识整理提纲。
- 异常问题索引。
- troubleshooting 文档。
- case note 模板。
- todo 和不确定项列表。
- Git 上传检查清单。

## RNA 工作流概览

当前 RNA 主线：

1. 从 `datasets.xlsx` 发现样本。
2. 读取 triplet、h5、archive 或 shared triplet。
3. 计算 `n_counts`、`n_genes`、`pct_mt`、`pct_ribo`。
4. doublet 检测。
5. dynamic hybrid MAD 过滤。
6. Scanpy normalization、log1p、HVG、PCA、neighbors、Leiden、UMAP。
7. Azimuth `pbmcref` 注释。
8. 映射到 5 类 final celltype。
9. 输出 metadata、QC、validation、h5ad、UMAP 和 status。
10. 可执行 baseline-only tuning。
11. 可进入 product-level Harmony 低维整合。

参考：`docs/rna_pipeline.md`。

## ATAC 工作流概览

当前 ATAC 主线：

1. 从 `atac.xlsx`、`co.xlsx` 或 longevity raw 目录发现样本。
2. 读取 fragments 和 barcode 辅助文件。
3. 使用 `peak.bed` 构建 peak-by-cell matrix。
4. 计算 TSS、FRiP、nucleosome、blacklist、fragments 等 QC。
5. scDblFinder doublet 检测。
6. MAD outlier 和 doublet 过滤。
7. LSI、聚类和 UMAP。
8. CIMA ATAC L1-L4 注释。
9. 输出 metadata、QC、validation、UMAP、h5ad、status 和日志。
10. 可进入 product-level BBKNN 低维整合。

参考：`docs/atac_pipeline.md`。

## 批处理检查规则

检查批处理时优先查看：

- 样本表来源是否正确。
- 输入路径是否存在。
- `run_status.json` 是否存在。
- `outputs_complete` 是否为 true。
- expected output 是否完整。
- `logs/sample_qc.log` 是否有错误。
- 是否使用了 `--force`、`--dry-run` 或 `--skip-complete`。

参考：`docs/batch_processing.md`。

## 异常诊断规则

遇到异常结果时，不要直接建议全量重跑。按顺序检查：

1. 是 RNA 还是 ATAC。
2. 是单样本还是多样本整合。
3. 是 QC、过滤、降维、聚类、投影、注释、导出还是批处理问题。
4. 输入文件是否正确。
5. metadata、barcode、sample ID 是否匹配。
6. QC 后细胞数是否合理。
7. 是否存在低质量样本或低质量 cluster。
8. 聚类参数是否记录。
9. 投影是否被 sample/batch 主导。
10. marker 是否存在；若项目中未找到，明确标注。
11. 注释规则是否可追溯。
12. 绘图使用的 label 是否是最新版本。
13. 日志是否报错。
14. 是否需要小范围重跑。
15. 是否需要人工确认。

参考：`references/diagnostic_workflow.md`, `docs/troubleshooting.md`, `docs/problem_index.md`。

## Git 上传安全规则

AI 不应建议上传以下文件：

- `data/` 真实内容。
- `output/` 真实内容。
- `QualityControl/` 全量内容。
- `ArchRLogs/` 全量内容。
- `.h5ad`, `.rds`, `.RData`, `.loom`。
- FASTQ/BAM/fragment/matrix 大文件。
- 隐私 metadata。
- 本地绝对路径配置。

AI 可以帮助上传或整理：

- 代码。
- 测试。
- 文档。
- 配置模板。
- 小型示例说明。
- AI skill 自身。

参考：`git_preparation_checklist.md`。

## case note 记录

每次异常问题排查后，建议使用：

```text
templates/case_note_template.md
```

记录内容包括：问题类型、现象、涉及数据、涉及文件、初步判断、检查顺序、检查结果、处理尝试、是否重跑、当前结论和后续待确认。

## 禁止或需确认的操作

未得到用户确认前，不要：

- 删除旧结果。
- 移动核心脚本。
- 重命名输出目录。
- 修改核心 QC 或注释逻辑。
- 全量重跑数据。
- 把大文件加入 Git。
- 把潜在问题写成已发生问题。
- 把当前项目中未找到的步骤写成已完成。
