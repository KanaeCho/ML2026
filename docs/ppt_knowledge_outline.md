# PPT 知识提纲

本文档为后续组会汇报、论文方法整理和个人复习准备 PPT 结构。每一页尽量标注可用结果文件路径；如果当前项目中未找到，则明确标注。

## 1. 项目背景

- 说明 RNA/ATAC 多样本分析项目目的。
- 强调当前项目已形成 RNA、ATAC、co、longevity 和 product-level 整合分支。
- 对应文件：`README.md`, `docs/project_summary.md`

## 2. 数据类型和样本结构

- RNA only。
- ATAC only。
- co RNA/ATAC。
- longevity 独立通道。
- 对应文件：`docs/io_summary.md`, `docs/batch_processing.md`

## 3. RNA 分析流程

- 输入 count matrix。
- QC 和过滤。
- normalize/log/HVG/PCA。
- clustering 和 UMAP。
- Azimuth `pbmcref` 注释。
- final 5 类 celltype。
- 对应文件：`docs/rna_pipeline.md`

## 4. ATAC 分析流程

- 输入 fragments。
- barcode 选择。
- peak-by-cell matrix。
- QC 和 doublet。
- LSI/UMAP。
- CIMA ATAC L1-L4 注释。
- 对应文件：`docs/atac_pipeline.md`

## 5. 多样本批处理设计

- 样本表来源。
- 单样本和 GSE 批处理。
- 跳过已完成样本。
- 日志和 run status。
- product organization。
- 对应文件：`docs/batch_processing.md`

## 6. RNA QC 指标解释

- `n_counts`。
- `n_genes`。
- `pct_mt`。
- `pct_ribo`。
- doublet。
- 对应图：`output/rna/{GSE}/{sample_id}/qc_overview.png`

## 7. ATAC QC 指标解释

- `nCount_ATAC`。
- `nFeature_ATAC`。
- `TSS.enrichment`。
- `FRiP`。
- `nucleosome_signal`。
- `blacklist_fraction`。
- 对应图：`output/atac/{GSE}/{GSM}/qc_overview.png` 或 `output/co/atac/{dataset}/{sample}/qc_overview.png`

## 8. 降维、聚类和投影的关系

- RNA：HVG -> PCA -> neighbors -> Leiden -> UMAP。
- ATAC：peak matrix -> LSI -> clustering -> UMAP。
- product：低维 embedding -> integration -> UMAP/cluster。
- 对应文件：`references/knowledge_map.md`

## 9. 注释策略

- RNA：Azimuth `pbmcref`。
- ATAC：CIMA ATAC centroid assignment。
- product RNA：CIMA projection-space 诊断标签。
- 对应文件：`references/annotation_notes.md`

## 10. RNA 和 ATAC 结果对应关系

- co 分支存在 paired sample 结构。
- 当前项目不做 paired-cell matching。
- 不应默认解释为 barcode 级一一对应。
- 对应文件：`docs/io_summary.md`

## 11. 当前主要结果

可用图：

- `output/rna/{GSE}/{sample_id}/umap_rna_pbmcref_vs_cima_l1.png`
- `output/rna/{GSE}/{sample_id}/umap_rna_cima_l1.png`
- `output/atac/{GSE}/{GSM}/umap_cima_cell_type_l1.png`
- `output/co/atac/{dataset}/{sample}/umap_cima_cell_type_l1.png`
- `output/1.only_atac/figures/*_panels.png`
- `output/2.only_rna/figures/*_panels.png`
- `output/3.co_atac/figures/*_panels.png`
- `output/4.co_rna/figures/*_panels.png`

## 12. 适合放入 PPT 的图和表

- QC overview 图。
- RNA CIMA L1 UMAP。
- RNA pbmcref highlight 图。
- ATAC CIMA L1/L2 UMAP。
- product-level CIMA L1/L2 panels。
- product-level integrated cluster panels。
- `qc_summary.csv` 汇总表。
- `integration_metrics.csv` 指标表。

## 13. 当前遇到或潜在的异常问题

已出现的整理风险：

- 本地路径和外部挂载依赖。
- 大文件和结果文件不适合上传 Git。
- dataset 特例硬编码。
- marker/gene activity/motif 主线未确认。

潜在问题：

- QC 异常。
- UMAP 混乱。
- 注释不清楚。
- batch effect。
- 批处理失败。
- 结果版本混乱。

## 14. 异常问题排查逻辑

- 先看数据类型。
- 再看样本和阶段。
- 查输入文件、metadata、QC summary、UMAP、annotation、log。
- 判断问题来源。
- 小范围检查优先。
- 需要修改流程时先人工确认。
- 对应文件：`references/diagnostic_workflow.md`, `docs/troubleshooting.md`

## 15. 当前项目不足和待确认内容

- R 环境未完整文档化。
- notebook 需要分类。
- marker gene 输出未确认。
- ATAC gene activity 未确认。
- motif 输出未确认。
- dataset exceptions 需要配置化。

## 16. 后续分析和整理方向

- 完善 RNA/ATAC pipeline 文档。
- 完善 results index。
- 整理 AI skill。
- 清理 Git 上传风险。
- 补充小型公开示例。
- 准备组会 PPT 和论文方法描述。
