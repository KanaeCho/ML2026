# ML2026 项目说明

## 项目定位
本分支当前聚焦单样本 scATAC-seq 流程，不再以多样本整合作为主目标。

本分支当前验收目标只有两项：

1. 对 `GSE190992` 和 `GSE283744` 的全部样本统一对齐到 `data/reference/peak.bed`，逐样本独立完成质控。
2. 为每个样本输出 QC 图、矩阵结果和 4 张 UMAP 图；UMAP 分别按 CIMA transfer 得到的 `cell_type_l1`、`cell_type_l2`、`cell_type_l3`、`cell_type_l4` 上色，并使用统一颜色体系。

补充说明：
- 当前 R / Signac 流程内部继续以 peak x cell 计数矩阵为主。
- `data/reference/peak.bed` 是两个数据集全部样本共用的 reference peak 集，不在样本之间混用不同 peak 定义。
- 多样本整合相关旧脚本已从本分支移除，不作为当前维护对象。

## 当前状态
- 当前数据集仍为 2 个 GSE，共 28 个 GSM：`GSE190992` 和 `GSE283744`。
- `process_single_sample.R` 已确认使用 `data/reference/peak.bed` 构建单样本 peak x cell 矩阵，并在 QC 后生成 query-native ATAC UMAP 坐标、reference-projected CIMA UMAP 图以及 CIMA L1-L4 注释结果。
- 当前单样本主流程支持两种输出模式：默认 `full` 输出完整矩阵 / metadata / RDS；`matrix-lite` 输出面向外部验证的精简矩阵与一个 CIMA L1 UMAP。
- CIMA scATAC reference 已确认来自 `data/reference/cima/CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`。
- 该参考的 AnnData `obs` 已确认包含逐细胞 `cell_type_l1`、`cell_type_l2`、`cell_type_l3`、`cell_type_l4`。
- 当前分支运行时保留的 CIMA 资产只有：上述 `.h5ad`、`data/reference/cima/cima_atac_celltype_hierarchy.csv`、`data/reference/cima/cima_atac_reference_lsi_features.tsv.gz`、`data/reference/cima/cima_atac_reference_l1_centroids.tsv`、`data/reference/cima/cima_atac_reference_l2_centroids.tsv`、`data/reference/cima/cima_atac_reference_l3_centroids.tsv`、`data/reference/cima/cima_atac_reference_l4_centroids.tsv`、`data/reference/cima/cima_atac_reference_model.json`、`data/reference/cima/README.md`。

## 目录结构

### `data/`
原始数据和共享参考文件。

#### `data/reference/`
- `datasets.xlsx`: 数据集统计表。
- `peak.bed`: 两个 GSE 全部样本统一使用的 peak 定义文件。
- `peaks.csv`: 带 peak_id 的参考索引文件。

#### `data/reference/cima/`
- `CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`: 当前单样本注释使用的 CIMA scATAC 参考图谱。
- `cima_atac_celltype_hierarchy.csv`: 从参考 `.h5ad` 提取的 L1/L2/L3/L4 层级映射表。
- `cima_atac_reference_lsi_features.tsv.gz`: 从参考 `.h5ad` 抽样构建的紧凑参考 LSI 模型特征表，包含 feature index、feature ID、IDF 和 LSI loadings。
- `cima_atac_reference_l1_centroids.tsv` / `cima_atac_reference_l2_centroids.tsv` / `cima_atac_reference_l3_centroids.tsv` / `cima_atac_reference_l4_centroids.tsv`: 同一参考模型在各层级上的 centroid 表。
- `cima_atac_reference_model.json`: 参考模型构建参数记录。
- `README.md`: 当前保留资产与删除策略说明。
- 其中 `.h5ad` 是注释来源与字段定义的权威参考，单样本运行时实际直接消费的是层级表和参考 LSI centroid 模型。

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
- 当前工作区中的 `output/` 映射到外部结果根，当前实际指向 `/mnt/g/ML2026_output`。
- 当前按 `output/{GSE}/{GSM}` 组织单样本结果。
- 当前跨样本的 ATAC-only 汇总 / 整合产物集中在 `output/1.only_atac/`。
- 当前 `output/reference/` 保留供外部输出驱动脚本直接读取的 `datasets.xlsx` 与 `cima/` 参考资产副本。
- 当前已确认的单样本基础输出包括：
  - `qc_overview.png`
  - `matrix/matrix.mtx`
  - `matrix/barcodes.tsv.gz`
  - `matrix/features.tsv.gz`
  - `metadata.csv`
  - `metadata_qc.csv`
  - `qc_summary.csv`
  - `GSM*_seurat_qc.rds`
- 当前已确认的单样本 UMAP / annotation 输出包括：
  - `umap_cima_cell_type_l1.png`
  - `umap_cima_cell_type_l2.png`
  - `umap_cima_cell_type_l3.png`
  - `umap_cima_cell_type_l4.png`
- 当前支持轻量输出模式 `matrix-lite`，仅保留：
  - `matrix/matrix.mtx`
  - `matrix/barcodes.tsv`
  - `matrix/features.tsv`
  - `validation_result.csv`
  - `umap_cima_cell_type_l1.png`
  - `qc_summary.csv`
- 当前已确认写回单样本 metadata 的新增列包括：
  - `cima_cell_type_l1`
  - `cima_cell_type_l2`
  - `cima_cell_type_l3`
  - `cima_cell_type_l4`
  - `cima_l4_score`
  - `cima_l4_score_margin`
  - `umap_atac_1`
  - `umap_atac_2`
  - `cima_ref_umap_1`
  - `cima_ref_umap_2`
- 当前已确认的跨样本 ATAC-only 目录包括：
  - `output/1.only_atac/qc_reports/`
  - `output/1.only_atac/accepted_integration_bbknn_gsm/`

### `scripts/process/`
当前主要工作目录。

当前主线相关文件：
- `GSM8671454.ipynb`: 单样本基准 notebook，用于确认流程逻辑。
- `process_single_sample.R`: 当前单样本主入口，负责 reference peak 对齐、QC 指标计算、CIMA L1-L4 注释、per-sample UMAP 绘制以及矩阵 / metadata 输出。
- `regenerate_qc_overview.R`: 基于已有 `*_seurat_qc.rds` 重绘 `qc_overview.png`。
- `integrate_accepted_matrix_lite.py`: 读取按数据集筛选通过的 matrix-lite 样本，复用 CIMA compact feature model 将各样本投影到统一 reference-derived LSI 空间；当前支持基线 `pooled-umap` 和低内存 `bbknn` 两种整合模式，可显式指定 `--input-root` 与 `--reference-dir`，输出一个整合 UMAP 及按 CIMA cell type、GSE/GSM 来源、健康状态和 QC 指标着色的可视化图、整合 metadata 与内存监控日志。
- `render_integration_umap_panels.py`: 基于已有整合 metadata 直接补画 panel UMAP 大图，不重跑整合；当前用于把 `accepted_integration_metadata.csv` 按 `CIMA L1`、`CIMA L2`、`GSE`、健康状态拆成“每类一个子图”的面板图。
- `summarize_matrix_lite_qc.R`: 汇总 `matrix-lite` 样本的 `qc_summary.csv`，生成跨样本 QC 指标表与基线 QC 图。
- `export_baseline_qc_excel.py`: 把基线 QC 汇总导出为 Excel，并按数据集分别判断“是否可以接受”。
- `run_single_sample_umap.R`: 基于已有单样本 `*_seurat_qc.rds` 对 `pass_qc` 细胞重建单样本 LSI / UMAP / 聚类，并可叠加外部标签映射或 RNA->ATAC label transfer 结果。
- `filter_integration_cells.py`: 基于已有整合后注释 metadata 对旧整合输入做 cell 子集过滤。
- `de_novo_annotate_integration_celltypes.py`: 基于 cluster marker、UMAP 邻近关系和样本偏置生成 de novo 细胞类型注释建议。
- `diagnose_bridge_clusters.py`: 基于整合 metadata、QC 指标和 batch 构成生成 bridge / dirty cluster 诊断报告。
- `review_celltype_annotation_validation.py`: 复核当前整合注释是否受到 marker 不支持、样本偏置或旧标签混合的影响。
- `pipeline.py`: Python 流程管理入口，负责样本发现、调度和日志。
- `download_from_datasets.py`: 从 `datasets.xlsx` 过滤样本并组织 GEO supplementary 下载任务。
- `export_h5ad_obs.py`: 将 GEO 提供的 `.h5ad` 文件中的 `obs` 元数据导出为 CSV/CSV.GZ，供单样本标签映射使用。
- `build_cima_reference_model.py`: 从 CIMA `.h5ad` 构建当前运行时使用的紧凑参考 LSI centroid 模型。

## 当前文件分工
- `process_single_sample.R`
  - 输入：`GSE`、`GSM`
  - 当前行为：自动发现 fragment 和 barcode 文件，读取统一 `peak.bed`，运行单样本 QC，并使用 CIMA scATAC 参考的紧凑 LSI centroid 模型对 QC 后细胞进行分层注释
  - 当前输出目录：`output/{GSE}/{GSM}`
  - 当前会输出 `cima_cell_type_l1` ~ `cima_cell_type_l4`、`cima_l4_score`、`cima_l4_score_margin`、query-native `umap_atac_1` / `umap_atac_2`、reference-space `cima_ref_umap_1` / `cima_ref_umap_2`，以及 4 张按分层标签上色的 reference-projected CIMA UMAP 图

- `build_cima_reference_model.py`
  - 输入：`data/reference/cima/CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`
  - 当前行为：按 L4 平衡抽样 CIMA 参考细胞，构建紧凑 LSI 参考模型并输出各层级 centroid 文件

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

- `integrate_accepted_matrix_lite.py`
  - 用途：对 matrix-lite 模式下按数据集分别判定“可以接受”的样本做跨样本联合可视化；`--method pooled-umap` 保留当前 pooled baseline，`--method bbknn` 走 GSM/GSE batch-aware 的低内存图整合；当前默认输入样本根目录为工作区 `output/`，默认接受样本清单位于 `output/1.only_atac/qc_reports/`，默认参考 LSI 资产目录位于 `output/reference/cima/`，输出目录可由 `--output-dir` 指向例如 `output/1.only_atac/accepted_integration_bbknn_gsm/`

- `render_integration_umap_panels.py`
  - 用途：对已生成的整合 metadata 复用现有 UMAP 坐标补画 panel 大图；当前默认输入为 `output/1.only_atac/accepted_integration_bbknn_gsm/accepted_integration_metadata.csv` 与 `output/reference/cima/cima_atac_celltype_hierarchy.csv`，输出 `accepted_integration_cima_l1_panels.png`、`accepted_integration_cima_l2_panels.png`、`accepted_integration_gse_panels.png`、`accepted_integration_health_status_panels.png`

- `summarize_matrix_lite_qc.R`
  - 用途：对全部 matrix-lite 样本的 `qc_summary.csv` 做跨样本汇总；当前默认从工作区 `output/` 读取单样本结果，并把 baseline QC 指标总表、异常样本表和基线 QC 图写到 `output/1.only_atac/qc_reports/`

- `export_baseline_qc_excel.py`
  - 用途：把 `output/1.only_atac/qc_reports/` 下的 baseline QC 汇总导出为 Excel，并按数据集内部分布生成“是否可以接受”列

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

- `filter_integration_cells.py`
  - 用途：基于已有整合后注释结果，对旧整合输入目录的 cell 子集做过滤并重建新的 merged 输入目录

- `de_novo_annotate_integration_celltypes.py`
  - 用途：不使用旧人工标签，仅依据当前 cluster marker peaks、UMAP 邻近关系和样本偏置生成一版 de novo 注释建议

- `diagnose_bridge_clusters.py`
  - 用途：基于整合 metadata、残余 `scDblFinder.score`、QC 指标、batch 构成和注释验证结果生成 bridge / dirty cluster 诊断报告

- `review_celltype_annotation_validation.py`
  - 用途：复核当前 celltype/subtype 注释是否受到 marker 不支持、样本偏置或旧标签混合的影响

- `export_baseline_qc_excel.py`
  - 用途：把 `output/1.only_atac/qc_reports/` 下的 baseline QC 汇总导出为 Excel，并按数据集内部分布生成“是否可以接受”列

- `pipeline.py`
  - 当前支持：`discover`、`run-sample`、`run-gse`、`download`、`status`
  - 当前 `run-sample` / `run-gse` 支持 `--output-profile` 与 `--output-root`
  - 其中 `--output-profile matrix-lite` 用于磁盘受限场景；当前默认输出根是工作区 `output/`（已映射到外部结果根），如有需要仍可用 `--output-root` 改写

- `download_from_datasets.py`
  - 当前默认过滤：`scATAC` + `fragment`
  - 当前默认下载文件类别包括：`fragment`、`barcode`、`singlecell`
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

4. 构建 peak x cell 矩阵。
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

8. 输出基础单样本结果。
   - `qc_overview.png`
   - `matrix/matrix.mtx`
   - `matrix/barcodes.tsv.gz`
   - `matrix/features.tsv.gz`
   - `metadata.csv`
   - `metadata_qc.csv`
   - `qc_summary.csv`
   - `GSM*_seurat_qc.rds`

9. 对 QC 后细胞做单样本降维与 UMAP。
   - 当前路线会同时保留 query-native ATAC UMAP 与 reference-projected CIMA UMAP。
   - query-native 路线为 `RunTFIDF()` → `FindTopFeatures(min.cutoff = "q0")` → `RunSVD()` → `RunUMAP(reduction = "lsi")`。
   - CIMA 绘图路线会复用 query 细胞投影到参考 LSI 空间后的 embedding 再计算 reference-space UMAP。
   - UMAP 仅针对 `pass_qc == TRUE` 的细胞计算。
   - 当前写回 metadata 的坐标列为 query-native `umap_atac_1` / `umap_atac_2` 与 reference-space `cima_ref_umap_1` / `cima_ref_umap_2`。

10. 基于 CIMA scATAC 参考做 L1-L4 注释。
    - 当前参考文件为 `data/reference/cima/CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`。
    - 当前参考 `obs` 中使用的关键列为 `cell_type_l1`、`cell_type_l2`、`cell_type_l3`、`cell_type_l4`。
    - 当前实现不是跨模态 RNA-to-ATAC，而是同模态 ATAC-to-ATAC 的 reference-model transfer：先从 CIMA `.h5ad` 按 L4 平衡抽样构建紧凑 TF-IDF/LSI 参考模型，再把 query 细胞投影到同一参考 LSI 空间中，按层级 centroid 做最近邻匹配。
    - 当前注释按层级逐层进行：先直接预测 L1，再在该 L1 的允许子集中预测 L2，再在对应 L2 子集中预测 L3，最后在对应 L3 子集中预测 L4。
    - 当前输出的 4 张 UMAP 分别按 `cima_cell_type_l1`、`cima_cell_type_l2`、`cima_cell_type_l3`、`cima_cell_type_l4` 上色，且主线绘图基于 reference-projected CIMA UMAP。
    - 当前 4 张 UMAP 使用统一层级配色：L1 先定义主色，L2/L3/L4 在各自 L1 主色下生成同色系渐变。

## CIMA 参考约束
- 当前唯一明确的公开入口是 `https://db.cngb.org/trueblood/cima/resource`。
- 当前已确认的公开下载路径为：`https://ftp.cngb.org/pub/SciRAID/trueblood/cima/CIMA_Resource/Cell_Atlas/CIMA_ATAC_3762242cells_338036peaks_compressed.h5ad`。
- 当前已确认其模态为 scATAC，并且与本项目的 `data/reference/peak.bed` 共享 338036 个 peak。
- 当前已确认的注释 metadata 列名为：`cell_type_l1`、`cell_type_l2`、`cell_type_l3`、`cell_type_l4`。
- 当前项目运行时不再保留与主线无关的 CIMA exploratory 文件；如需额外保留，只能是被当前主线直接消费的文件或从参考 `.h5ad` 派生出的紧凑参考模型资产。

## QC 原则
- 优先使用原始数据附带的 filtered barcodes 作为初筛。
- 缺失官方 filtered barcodes 时，必须提供可重复的 fallback 流程。
- QC 判断以样本分布和可解释性为主，不在主流程中提前写死不必要的全局硬阈值。
- `blacklist_fraction` 作为常规 QC 指标保留。
- 所有过滤结果都应记录在 metadata 中，而不是只保留最终细胞列表。
- 单样本 UMAP 与 transfer 注释都必须建立在明确记录的 QC 后细胞集合之上。

## 文档同步约定
- 只要修改了项目结构、处理流程、脚本职责、输出规范或参考注释来源，就必须同步更新 `AGENTS.md`。
- `AGENTS.md` 记录的是当前分支的实际目标、已确认事实和明确待定项；不能把未验证的假设写成既成事实。
