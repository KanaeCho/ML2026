# PPT 提纲模板

## 1. 项目背景

- 研究问题：
- 数据类型：
- 样本来源：
- 当前分析目标：
- 对应文件：

## 2. 数据类型和样本结构

- RNA 样本：
- ATAC 样本：
- 共测样本：
- longevity 或特殊分支：
- 对应文件：`docs/io_summary.md`

## 3. RNA 分析流程

- 输入：
- QC：
- 降维/聚类/UMAP：
- 注释：
- 输出：
- 对应文件：`docs/rna_pipeline.md`

## 4. ATAC 分析流程

- 输入：
- QC：
- peak matrix：
- LSI/聚类/UMAP：
- CIMA 注释：
- 输出：
- 对应文件：`docs/atac_pipeline.md`

## 5. 多样本批处理设计

- 样本表来源：
- 单样本重跑：
- 跳过已完成：
- 日志和状态：
- 对应文件：`docs/batch_processing.md`

## 6. RNA QC 指标解释

- `n_counts`
- `n_genes`
- `pct_mt`
- `pct_ribo`
- doublet
- 对应图：`qc_overview.png`

## 7. ATAC QC 指标解释

- `nCount_ATAC`
- `TSS.enrichment`
- `FRiP`
- `nucleosome_signal`
- `blacklist_fraction`
- doublet
- 对应图：`qc_overview.png`

## 8. 降维、聚类和投影的关系

- PCA/LSI：
- neighbors graph：
- Leiden/Seurat clusters：
- UMAP：

## 9. 注释策略

- RNA：Azimuth `pbmcref`
- ATAC：CIMA ATAC
- product-level 诊断标签：
- 对应文件：`references/annotation_notes.md`

## 10. RNA 和 ATAC 结果对应关系

- 当前是否共测：
- 是否做 paired-cell matching：当前项目不做
- 需要注意的解释限制：

## 11. 当前主要结果

- RNA UMAP：
- ATAC UMAP：
- product-level panels：
- QC summary：

## 12. 适合放入 PPT 的图和表

- `qc_overview.png`
- `umap_rna_pbmcref_vs_cima_l1.png`
- `umap_rna_cima_l1.png`
- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`
- `figures/*_panels.png`

## 13. 当前遇到或潜在的异常问题

- 已出现：
- 潜在问题：
- 对应文件：`docs/problem_index.md`

## 14. 异常问题排查逻辑

- 先定位数据类型。
- 再定位样本和阶段。
- 查输入、metadata、QC、UMAP、注释、日志。
- 小范围检查优先于全量重跑。

## 15. 当前项目不足和待确认内容

- marker gene：
- gene activity：
- motif：
- notebook：
- R 环境：

## 16. 后续分析和整理方向

- 继续补文档：
- 清理 Git 上传风险：
- 整理 AI skill：
- 准备论文方法或组会材料：
