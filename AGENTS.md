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
- `run_batch_integration.R`: 对 sketch 或全量 merged 输入运行 TF-IDF、LSI、Harmony、UMAP、聚类和整合质量评估。
- `annotate_integration_celltypes.R`: 基于整合结果中的 cluster，为细胞追加 broad cell type / subtype 标注。
- `review_hard_qc_and_integration_readiness.py`: 汇总 hard-QC 前后表现，生成整合准备报告。
- `pipeline.py`: Python 流程管理入口，负责样本发现、调度和日志。
- `download_from_datasets.py`: 从 `datasets.xlsx` 过滤样本并组织 GEO supplementary 下载任务。

## 当前文件分工
- `process_single_sample.R`
  - 输入：`GSE`、`GSM`
  - 行为：自动发现 fragment 和 barcode 文件，运行单样本 QC，保存最终产物
  - 输出目录：`output/{GSE}/{GSM}`

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

- `run_batch_integration.R`
  - 用途：对 sketch 或全量 merged 输入运行正式整合与质量评估
  - 当前方法：`Signac + Seurat + harmony`
  - 当前默认：按 peak 可及细胞数筛选 top `30000` 个 peaks
  - 输出：`output/integration_sketch_analysis/` 或 `output/integration_merged_analysis/`

- `annotate_integration_celltypes.R`
  - 用途：在整合结果的 `integrated_metadata.csv.gz` 基础上追加细胞类型标签
  - 输出：
    - `cluster_celltype_annotation_map.csv`
    - `integrated_metadata_celltyped.csv.gz`
    - `post_harmony_umap_by_celltype.png`
    - `post_harmony_umap_by_celltype_subtype.png`
    - `celltype_annotation_notes.md`

- `review_hard_qc_and_integration_readiness.py`
  - 用途：复核 hard-QC 固定阈值在全队列上的实际效果，并输出整合准备报告
  - 输出：`output/hard_qc_review/`

- `pipeline.py`
  - 当前支持：`discover`、`run-sample`、`run-gse`、`download`、`status`

- `download_from_datasets.py`
  - 当前默认过滤：`scATAC` + `fragment`
  - 当前支持通过 `--file-kinds` 指定 GEO 文件类别，例如 `fragment`、`barcode`、`singlecell`、`summary`

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
