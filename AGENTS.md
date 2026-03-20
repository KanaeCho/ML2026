# ML2026 项目说明

## 项目定位
本项目当前只做两件事：

1. 对单样本 scATAC-seq 数据进行质控，输出 QC 图，并生成可用于后续分析的矩阵结果。
2. 将所有样本矩阵整合起来，输出整合后的可视化结果，用于判断整合质量和批次效应情况。

补充说明：
- 当前 R/Signac 流程内部主要使用 peak×cell 计数矩阵。
- 如果后续分析必须统一成 cell×peak 方向，应在输出规范中明确写清并统一处理，不要在不同脚本里混用。

## 当前状态
- 当前数据集为 2 个 GSE，共 28 个 GSM：`GSE190992` 和 `GSE283744`。
- 28 个样本的单样本 QC 已完成。
- 28 个样本的整合前 hard-QC 已完成。
- `integration_merged/` 与 `integration_sketch/` 已完成构建。
- `integration_sketch_analysis/` 与 `integration_merged_analysis/` 已完成 LSI / Harmony / UMAP / 聚类和整合质量评估输出。
- 当前已提供整合后细胞类型注释辅助脚本 `annotate_integration_celltypes.R`，用于在整合结果上追加 broad cell type / subtype 标注。

## 目录结构

### `data/`
原始数据和共享参考文件。

#### `data/reference/`
- `datasets.xlsx`: 数据集统计表。
- `peak.bed`: 统一 peak 定义文件。
- `peaks.csv`: 带 peak_id 的参考索引文件。

#### `data/raw/`
- 按 GSE 组织原始输入数据。
- 每个样本至少应包含 fragment 文件。
- 如果有 `*_filtered_barcodes.tsv.gz`，优先作为初始细胞集合使用。
- 如果有 `*singlecell*.csv.gz`，优先作为无官方 filtered barcodes 样本的 cell-calling 来源。
- 文件发现依赖命名约定：
  - fragment 文件名以 `GSM` 开头，且包含 `fragments`
  - barcode 文件名以 `GSM` 开头，且包含 `barcodes`

### `output/`
- 运行时生成的结果目录，不是项目的长期真值来源。
- 当前按 `output/{GSE}/{GSM}` 组织单样本结果。
- 当前还包含以下整合相关目录：
  - `output/integration_merged/`
  - `output/integration_sketch/`
  - `output/integration_sketch_analysis/`
  - `output/integration_merged_analysis/`
  - `output/hard_qc_review/`

### `scripts/process/`
当前主要工作目录。

关键文件：
- `GSM8671454.ipynb`: 单样本基准 notebook，用于确认流程逻辑。
- `process_single_sample.R`: 单样本 QC 与矩阵生成脚本。
- `regenerate_qc_overview.R`: 基于已有 `*_seurat_qc.rds` 重生 `qc_overview.png`。
- `apply_integration_hard_qc.py`: 在单样本 `pass_qc` 结果上追加统一 hard-QC。
- `merge_integration_matrices.py`: 严格检查 feature 一致性后合并多样本矩阵。
- `build_integration_sketch.py`: 从各样本 `integration_qc` 结果中均衡抽样，构建 sketch 输入。
- `filter_integration_cells.py`: 基于整合后注释 metadata 按指定列过滤细胞，生成新的 merged 输入目录。
- `run_batch_integration.R`: 对 sketch 或全量 merged 输入运行 TF-IDF、LSI、Harmony、UMAP、聚类和整合质量评估。
- `run_single_sample_umap.R`: 基于已有单样本 `*_seurat_qc.rds` 对 hard-QC 前的 `pass_qc` 细胞重建单样本 LSI / UMAP / 聚类，并可叠加 GEO 标签映射或 scRNA->scATAC label transfer 结果。
- `annotate_integration_celltypes.R`: 基于整合结果中的 cluster，为细胞追加 broad cell type / subtype 标注。
- `de_novo_annotate_integration_celltypes.py`: 仅基于当前 cluster marker peaks、UMAP 邻近关系和样本偏置，重算一版 de novo 细胞类型注释表。
- `review_celltype_annotation_validation.py`: 基于当前整合结果、cluster marker 和样本组成，输出注释验证汇总与人工复核报告。
- `review_hard_qc_and_integration_readiness.py`: 汇总 hard-QC 前后表现，生成整合准备报告。
- `pipeline.py`: Python 流程管理入口，负责样本发现、调度和日志。
- `download_from_datasets.py`: 从 `datasets.xlsx` 过滤样本并组织 GEO supplementary 下载任务。
- `export_h5ad_obs.py`: 将 GEO 提供的 `.h5ad` 文件中的 `obs` 元数据导出为 CSV/CSV.GZ，供单样本标签映射使用。

## 当前文件分工
- `process_single_sample.R`
  - 输入：`GSE`、`GSM`
  - 行为：自动发现 fragment 和 barcode 文件，运行单样本 QC，保存最终产物
  - 输出目录：`output/{GSE}/{GSM}`

- `run_single_sample_umap.R`
  - 用途：读取已有 `output/{GSE}/{GSM}/{GSM}_seurat_qc.rds`，仅对 `pass_qc` 细胞运行单样本 TF-IDF、LSI、聚类和 UMAP
  - 当前支持：
    - 直接读取外部注释 CSV，按 barcode/sample 映射 GEO 细胞类型并绘制单样本 UMAP
    - 使用外部 scRNA Seurat reference RDS 做 `GeneActivity + FindTransferAnchors/TransferData` 的 RNA->ATAC 标签转移
  - 当前输出：
    - `single_sample_umap_by_cluster.png`
    - `single_sample_umap_by_qc.png`
    - `single_sample_umap_metadata.csv.gz`
    - 如提供参考标签，还会额外输出按 GEO label 或 transferred label 着色的 UMAP 图和 `single_sample_umap_report.md`

- `regenerate_qc_overview.R`
  - 用途：只重绘总 QC 图，不重跑整套单样本 QC

- `apply_integration_hard_qc.py`
  - 用途：基于现有 `metadata_qc.csv` 与 `matrix/` 再做一轮统一硬阈值筛选
  - 输出：每个样本单独的 `integration_qc/`

- `merge_integration_matrices.py`
  - 用途：将各样本 `integration_qc` 矩阵严格对齐后横向合并
  - 约束：只有所有 `features.tsv.gz` 内容与顺序完全一致时才允许直接合并
  - 输出：`output/integration_merged/`

- `build_integration_sketch.py`
  - 用途：从 hard-QC 后细胞中按样本均衡抽样，构建 sketch 输入
  - 当前默认：每样本抽取 `1000` 个细胞
  - 输出：`output/integration_sketch/`

- `filter_integration_cells.py`
  - 用途：基于已有整合后注释结果，对 `integration_merged/` 的 cell 子集做过滤，并重建新的 merged 输入目录
  - 当前可用于删除 `celltype == Unknown / non-PBMC-like` 的细胞，再重跑整合
  - 当前也支持按 `seurat_clusters` 整群排除，或按 `scDblFinder.score` 做全局/cluster 定向细胞剔除，用于敏感性重跑
  - 输入：
    - `output/integration_merged/`
    - `output/integration_merged_analysis/integrated_metadata_celltyped.csv.gz`
  - 输出：例如 `output/integration_merged_without_unknown/`

- `run_batch_integration.R`
  - 用途：对 sketch 或全量 merged 输入运行正式整合与质量评估
  - 当前方法：`Signac + Seurat + harmony`
  - 当前默认：按 peak 可及细胞数筛选 top `30000` 个 peaks
  - 当前支持通过 `--drop-lsi-dims` 或自动 QC 相关性规则，从 downstream Harmony / 邻居图 / UMAP 中排除受 QC 污染的 LSI 维度
  - 当前支持：`--skip-mixing-metrics`，用于在只关心重整合 / 聚类 / UMAP 时跳过 batch mixing 计算
  - 输出：`output/integration_sketch_analysis/` 或 `output/integration_merged_analysis/`

- `diagnose_bridge_clusters.py`
  - 用途：基于整合后的 `integrated_metadata.csv.gz`、残余 `scDblFinder.score`、QC 指标、batch 构成和注释验证结果，生成 bridge / dirty cluster 诊断报告
  - 当前输出：
    - `bridge_cluster_diagnostic_summary.csv`
    - `bridge_candidate_doublet_cells.csv.gz`
    - `bridge_diagnostic_report.md`
    - `bridge_diagnostic_summary.json`

- `annotate_integration_celltypes.R`
  - 用途：在整合结果的 `integrated_metadata.csv.gz` 基础上追加细胞类型标签
  - 当前支持：
    - 直接使用内置的 `cluster -> celltype/subtype` 映射
    - 通过 `--reference-metadata` 基于共享 `global_cell_id` 对重整合后的新 cluster 做多数票标签转移
    - 通过 `--annotation-map` 读取外部注释表，并通过 `--output-suffix` 输出平行版本结果而不覆盖原文件
    - broad cell type UMAP 使用固定颜色映射，保证不同注释版本之间同名 celltype 颜色一致，便于横向对比
  - 输出：
    - `cluster_celltype_annotation_map.csv`
    - `integrated_metadata_celltyped.csv.gz`
    - `post_harmony_umap_by_celltype.png`
    - `post_harmony_umap_by_celltype_subtype.png`
    - `celltype_annotation_notes.md`

- `de_novo_annotate_integration_celltypes.py`
  - 用途：不使用旧人工标签，仅依据当前 `cluster_top_accessible_peaks.csv` 的 nearest genes、当前 UMAP cluster 邻近关系和样本偏置，生成一版 de novo 注释建议
  - 输入：
    - `cluster_top_accessible_peaks.csv`
    - `integrated_metadata.csv.gz`
    - `cluster_by_sample_counts.csv`
  - 输出：
    - `cluster_celltype_annotation_map_de_novo.csv`
    - `cluster_marker_annotation_scores_de_novo.csv`
    - `celltype_annotation_notes_de_novo.md`
    - 可配合 `annotate_integration_celltypes.R --annotation-map ... --output-suffix _de_novo` 生成 `integrated_metadata_celltyped_de_novo.csv.gz` 和对应 UMAP 图

- `review_celltype_annotation_validation.py`
  - 用途：复核当前 celltype/subtype 注释是否受到 marker 不支持、样本偏置或旧标签混合的影响
  - 输入：
    - `integrated_metadata.csv.gz`
    - `cluster_celltype_annotation_map.csv`
    - `cluster_top_accessible_peaks.csv`
    - `cluster_annotation_reference_overlap.csv`（如果存在）
  - 输出：
    - `cluster_annotation_validation_summary.csv`
    - `annotation_validation_report.md`

- `review_hard_qc_and_integration_readiness.py`
  - 用途：复核 hard-QC 固定阈值在全队列上的实际效果，并输出整合准备报告
  - 输出：`output/hard_qc_review/`

- `pipeline.py`
  - 当前支持：`discover`、`run-sample`、`run-gse`、`download`、`status`

- `download_from_datasets.py`
  - 当前默认过滤：`scATAC` + `fragment`
  - 当前支持通过 `--file-kinds` 指定 GEO 文件类别，例如 `fragment`、`barcode`、`singlecell`、`summary`

- `export_h5ad_obs.py`
  - 用途：从 GEO 的 `.h5ad` processed object 中提取 `obs` 级别注释，导出为后续 `run_single_sample_umap.R --annotation-csv` 可用的表格

## 单样本处理流程
1. 读取样本参数。
   - 输入至少包括 `GSE`、`GSM`。
   - fragment 和 barcode 文件按命名规则自动发现。

2. 确定初始 barcode 集合。
   - 优先使用 `*_filtered_barcodes.tsv.gz`。
   - 若缺失，则优先使用 `*singlecell*.csv.gz` 中的官方 cell-calling 结果。
   - 若两者都缺失，则在脚本内部完成 fallback barcode 预筛。
   - `singlecell.csv.gz` 中的附加逐 barcode 指标会带 `singlecell_` 前缀写入 metadata。

3. 读取统一参考。
   - 加载 `data/reference/peak.bed`。
   - 准备 hg38 注释并统一染色体命名。

4. 构建 peak×cell 矩阵。
   - 检查 tabix 索引，缺失时自动创建。
   - 使用 `FeatureMatrix()` 在统一 peak 集合上生成计数矩阵。

5. 构建 Seurat / Signac 对象并计算 QC 指标。
   - `nCount_ATAC`
   - `nFeature_ATAC`
   - `TSS.enrichment`
   - `nucleosome_signal`
   - `total_fragments`
   - `FRiP`
   - `unique_ratio`
   - `blacklist_fraction`

6. 进行 doublet 检测。
   - 当前使用 `scDblFinder`。

7. 进行 QC 过滤。
   - 当前主流程使用 MAD 方式识别异常值。
   - 当前主要基于 `nCount_ATAC`、`TSS.enrichment`、`FRiP`。
   - doublet 检测结果参与最终过滤。

8. 输出单样本结果。
   - `qc_overview.png`
   - `matrix/matrix.mtx`
   - `matrix/barcodes.tsv.gz`
   - `matrix/features.tsv.gz`
   - `metadata.csv`
   - `metadata_qc.csv`
   - `qc_summary.csv`
   - `GSM*_seurat_qc.rds`

9. 整合前统一 hard-QC。
   - 这一步不覆盖原始单样本 QC 结果。
   - 当前默认阈值：
     - `nCount_ATAC >= 1000`
     - `nCount_ATAC <= 100000`
     - `TSS.enrichment >= 4`
     - `FRiP >= 0.35`
     - `blacklist_fraction <= 0.05`
     - `nucleosome_signal <= 4`
   - 当前不把 `unique_ratio` 作为统一 hard-QC 阈值。
   - 输出目录为每个样本下的 `integration_qc/`。

## 多样本整合流程
1. 对 28 个 `integration_qc` peak×cell 矩阵做严格 feature 对齐检查，并生成总合并矩阵。
2. 基于 sketch 或全量 merged 输入运行 LSI / Harmony / UMAP / 聚类。
3. 输出整合质量评估结果，并生成整合 metadata。
4. 在需要时追加 cluster 到 cell type 的标注结果。

当前整合输出至少包括：
- Harmony 前后按 `GSE`、`GSM` 着色的低维图
- 按 `nCount_ATAC`、`FRiP`、`TSS.enrichment`、`blacklist_fraction` 着色的图
- `cluster × sample` 热图
- LSI 与 QC 指标相关性
- 至少一种 batch mixing 指标
- `cluster_top_accessible_peaks.csv`
- `integrated_metadata.csv.gz`

## QC 原则
- 优先使用原始数据附带的 filtered barcodes 作为初筛。
- 缺失官方 filtered barcodes 时，必须提供可重复的 fallback 流程。
- QC 判断以样本分布和可解释性为主，不在主流程中提前写死全局硬阈值。
- `blacklist_fraction` 作为常规 QC 指标保留。
- 所有过滤结果都应记录在 metadata 中，而不是只保留最终细胞列表。
- 如果统一阈值发生变化，必须同步更新脚本和本文件。

## 文档同步约定
- 只要修改了项目结构、处理流程、脚本职责或输出规范，就必须同步更新 `AGENTS.md`。
- `AGENTS.md` 记录的是当前实际状态，不保留已经删除、停用或仅具历史意义的说明。
