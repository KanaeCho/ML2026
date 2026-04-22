# ML2026 项目说明

## 项目定位
当前分支聚焦 **单样本 scRNA-seq** 主线。

本分支第一阶段验收目标：

1. 对 `datasets.xlsx` 中已筛选、且原始计数矩阵已下载到本地的数据，逐样本完成 RNA `QC + 聚类 + CIMA RNA L1/L2 注释 + UMAP 可视化`。
2. 在当前主线上，保留一个 **baseline-only tuning 工作流**，用于以固定 `baseline__baseline__baseline` 参数执行单样本 RNA 主线，并输出可复核的单 candidate 审计产物。
3. 重点保证 `L1 + L2` 的 celltype UMAP 质量可读，并同时输出可复核的简单审计指标。

补充说明：
- 当前以 **单样本** 为主，不把跨样本 RNA 整合作为第一阶段验收目标。
- CIMA RNA 注释是当前主线；marker-based 注释目前不作为主流程必需产物。
- 仓库中的 ATAC / TEA-seq 流程继续保留，但在本分支降级为非主线维护对象。

## 当前状态
- 当前工作区默认通过外部数据根读取数据：优先 `./data`，否则读取环境变量 `ML2026_DATA_ROOT`，再否则回退到 `/mnt/g/ML2026_data`。
- 当前工作区中的 `output/` 是软链接，实际指向 `/mnt/g/ML2026_output`。
- `pipeline.py discover-rna` 当前可识别 10 个已筛选 RNA GSE，共 55 个样本条目，其中 54 个可直接运行，1 个当前仅记为不支持：
  - 支持：`GSE149689`、`GSE157007`、`GSE167363`、`GSE192391`、`GSE198891`、`GSE213516`、`GSE226039`、`GSE231794`、`GSE268936`
  - 不支持：`GSE198533` 的共享 `gene_counts_matrix.csv.gz`，当前不按单细胞矩阵处理
- 已验证的 RNA smoke sample：
  - `GSE167363/GSM5102900`
  - 成功输出 `metadata.csv`、`metadata_qc.csv`、`qc_summary.csv`、`validation_result.csv`、`qc_overview.png`、`umap_rna_pbmcref_vs_cima_l1.png`、矩阵导出、`run_status.json` 和 `GSM5102900.h5ad`
- 当前主线默认 QC 已切换到更严格的一档：`counts_lower_nmads=2.5`、`genes_lower_nmads=2.5`、`pct_mt_upper_nmads=2.5`、`pct_ribo_upper_nmads=3.5`
- 当前代码中的 baseline-only tuning 路径已经落地：
  - `pipeline.py` / `scripts.only_rna.cli` 已支持 `tune-rna-sample` 与 `tune-rna-gse`
  - tuning 当前只执行固定 candidate：`baseline__baseline__baseline`
  - 该 candidate 会真实执行单样本主线（读取矩阵、QC、doublet、embedding、Azimuth、样本输出），不是占位 stub
- `GSE226039` 当前按文件名只保留 `PBMC` 样本参与 RNA 发现与运行，不把 Ileum / Rectum 组织一起纳入本分支主线。

## 目录结构

### `data/`
- 当前仓库内通常不存在真实 `data/` 目录；运行时数据默认来自外部数据根。
- 共享参考位于 `data_root/reference/`，原始 RNA 输入位于 `data_root/raw/`。

### `data_root/reference/`
- `datasets.xlsx`
  - RNA 样本发现基于 Excel 第一张表中当前可见、且 `assay=scRNA` 的行。
- `cima/`
  - `CIMA_RNA_6484974cells_36326genes_compressed.h5ad`
  - `cima_rna_celltype_hierarchy.csv`
  - `cima_rna_reference_pca_features.tsv.gz`
  - `cima_rna_reference_l1_centroids.tsv`
  - `cima_rna_reference_l2_centroids.tsv`
  - `cima_rna_reference_l3_centroids.tsv`
  - `cima_rna_reference_l4_centroids.tsv`
  - `cima_rna_reference_model.json`

### `data_root/raw/`
- 按 `GSE` 组织原始 RNA 输入。
- 当前 RNA 单样本发现支持：
  - 每个 `GSM` 一套 `matrix.mtx(.gz) + barcodes.tsv(.gz) + features/genes.tsv(.gz)`
  - 每个 `GSM` 一个 10x `.h5`
  - 每个 `GSM` 一个 `matrix.tar.gz`
  - 一个 `GSE` 共享的 Matrix Market triplet，此时 `sample_id = GSE`
- 当前不支持把共享 gene-count CSV 直接当作单细胞矩阵主线输入。
- 对 `GSE226039`，只接受文件名包含 `PBMC` 的样本。

### `output/`
- 运行结果不是项目长期真值来源。
- RNA 单样本默认输出目录：`output/rna/{GSE}/{sample_id}/`
- 当前每个 RNA 样本至少输出：
  - `metadata.csv`
  - `metadata_qc.csv`
  - `qc_summary.csv`
  - `qc_thresholds.json`
  - `validation_result.csv`
  - `qc_overview.png`
  - `{sample_id}.h5ad`
  - `matrix/matrix.mtx`
  - `matrix/barcodes.tsv.gz`
  - `matrix/features.tsv.gz`
  - `umap_rna_pbmcref_vs_cima_l1.png`
  - `umap_rna_pbmcref_highlight.png`
  - `umap_rna_cima_l1.png`
  - `run_status.json`
- 当前 baseline-only tuning 额外输出：
  - 根目录：`output/rna/{GSE}/{sample_id}/tuning/`
  - 至少包含：
    - `candidates.csv`
    - `selection_summary.json`
    - `selected_params.json`
    - `umap_rna_candidates_overview_cima_l1.png`
    - `umap_rna_candidates_overview_pbmcref.png`
    - `umap_rna_candidates_overview_pbmcref_highlight.png`
  - 单 candidate 的样本级输出当前写到：
    - `output/rna/{GSE}/{sample_id}/tuning/baseline__baseline__baseline/{GSE}/{sample_id}/`
    - 其内部仍沿用常规单样本输出族（`metadata.csv`、`qc_summary.csv`、`.h5ad`、UMAP 图等）
- 仓库中旧有的 ATAC / TEA-seq 输出仍位于：
  - `output/{GSE}/{GSM}/`
  - `output/1.only_atac/`
  - `output/GSE214546/qc_audit/`
  但这些目录不属于本分支 RNA 第一阶段验收主线。

## 当前主线脚本

### `scripts/only_rna/`
- 当前 RNA 主线已切换为 Python-first 子系统。
- 当前已实现模块：
  - `config.py`：YAML 默认配置与 CLI override 合并；当前支持 `qc`、`plotting`、可选 `annotation.methods`，以及可选 `embedding` / `azimuth` / `tuning` 配置面；默认 tuning 为 baseline-only
  - `discovery.py`：数据根解析、`datasets.xlsx` 可见 `scRNA` 行发现、样本本地布局发现
  - `read_inputs.py`：triplet / `.h5` / `matrix.tar.gz` / GSE-shared triplet 读取为 `AnnData`
  - `qc.py`：`n_counts`、`n_genes`、`pct_mt`、`pct_ribo` 计算与 QC fail flags / `pass_qc`
  - `doublet.py`：Python doublet 结果归一化与 `scanpy.pp.scrublet` 路径
  - `embedding.py`：仅对 `pass_qc` 细胞执行 embedding / clustering / UMAP，并写回 `cluster`、`umap_1`、`umap_2`；当前会优先使用显式 `config.embedding`，若保持 dataclass 默认值则回退到原有 PBMC 启发式默认参数
  - `annotation.py`：保留最小 CIMA 资产加载与多方法编排；当前 mainline 默认执行 shared Azimuth `pbmcref` 注释，并把 `annotation_method_status` 写入 `adata.uns`
  - `azimuth.py`：shared Azimuth `pbmcref` 执行模块；提供 `run_azimuth_annotation(...)`、显式 `status/detail` 语义，以及真实 R/Seurat/Azimuth 路径
  - `qc_calibration.py`：跨样本执行 baseline 与 `stricter_v1` QC 对照的小型校准脚本；当前 `stricter_v1` 使用 `counts_lower_nmads=2.5`、`genes_lower_nmads=2.5`、`pct_mt_upper_nmads=2.5`
  - `tuning_presets.py`：baseline-only preset 定义（QC / Azimuth / embedding 均只保留 `baseline`）
  - `tuning_metrics.py`：单 candidate 的 `qc_score` / `annotation_score` / `embedding_score` 以及汇总逻辑
  - `tuning_orchestrator.py`：固定 `baseline__baseline__baseline` candidate 的执行与审计写出
  - `plotting.py`：按配置输出可读的 categorical UMAP
  - `outputs.py`：写出 `.h5ad`、metadata/QC/validation CSV、矩阵导出、UMAP 图；当前也负责 tuning 选择产物写出
  - `cli.py`：`discover-rna` / `run-rna-sample` / `run-rna-gse` / `tune-rna-sample` / `tune-rna-gse` / `rna-status` 的 Python 主线实现

### `scripts/process/pipeline.py`
- 当前支持的 RNA 命令：
  - `discover-rna`
  - `run-rna-sample`
  - `run-rna-gse`
  - `tune-rna-sample`
  - `tune-rna-gse`
  - `rna-status`
- 当前行为：
  - 保持上述命令名稳定
  - RNA 命令分发到 `scripts.only_rna.cli`
  - `run-rna-sample` 显式拒绝 `gse_shared` 目标
  - `run-rna-gse` 会隐式包含受支持的 GSE-shared triplet
  - `tune-rna-sample` 当前同样显式拒绝 `gse_shared` 目标，并只执行 `baseline__baseline__baseline`
  - `tune-rna-gse` 当前会像 `run-rna-gse` 一样隐式包含受支持的 GSE-shared triplet，并只执行 `baseline__baseline__baseline`
  - `rna-status` 当前按新的 RNA 输出族（包含 `{sample_id}.h5ad`）检查完成度

### `scripts/process/build_cima_rna_reference_model.py`
- 输入：
  - `CIMA_RNA_6484974cells_36326genes_compressed.h5ad`
- 当前行为：
  - 按 `L4` 平衡抽样参考细胞
  - 依据 `variances_norm` 选择高变基因
  - 对参考表达矩阵做标准化 PCA
  - 输出当前运行时使用的 RNA compact PCA feature model、L1-L4 centroid、层级表和 model json
- 当前主要用途：
  - 重建全量 RNA reference 资产
  - 做 PBMC-focused / 子集 reference 试验

## RNA 单样本处理流程
1. 发现样本。
   - 只读取 `datasets.xlsx` 中当前可见的 scRNA 行。
   - 只处理本地已存在矩阵文件的样本。

2. 读取计数矩阵。
   - 当前支持 triplet、`.h5`、`matrix.tar.gz`。
   - 若命中共享 GSE 级 Matrix Market triplet，则用 `GSE` 作为 `sample_id`。

3. 计算基础 QC。
   - `n_counts`
   - `n_genes`
   - `pct_mt`
   - `pct_ribo`
   - `doublet_score`
   - `is_doublet`

4. 执行 QC 过滤。
   - 当前主线使用 dynamic hybrid MAD 阈值：
     - `n_counts` / `n_genes` 在 `log10(x + 1)` 空间做 lower-tail dynamic threshold
     - `pct_mt` / `pct_ribo` 在原值空间做 upper-tail dynamic threshold
     - 最终阈值会受 guardrail 约束并写入 `qc_thresholds.json`
   - `is_doublet == True` 的细胞不会进入 `pass_qc`

5. 对 `pass_qc` 细胞做 RNA 降维与聚类。
   - `scanpy.pp.normalize_total`
   - `scanpy.pp.log1p`
   - `scanpy.pp.highly_variable_genes`
   - `scanpy.pp.scale`
   - `scanpy.tl.pca`
   - `scanpy.pp.neighbors`
   - `scanpy.tl.leiden`（若环境缺 `igraph`，当前会回退到最小 deterministic fallback）
   - `scanpy.tl.umap`（或小样本 / fallback 写回）

6. 基于 shared Azimuth `pbmcref` 做主线 RNA 注释。
   - 当前 mainline 默认通过 `run_azimuth_annotation(...)` 调用 R/Seurat/Azimuth `RunAzimuth(..., reference='pbmcref')`
   - 默认输出 `azimuth_cell_type`，并在 `adata.uns['annotation_method_status']` 中记录 `status/detail`
   - CIMA / 其他方法仍保留在多方法与 tuning 路径中，但不再是单样本 sample-root 完成度的主线契约

7. 写回 metadata 和 UMAP 图。
   - `metadata.csv` 保留全细胞信息
   - `metadata_qc.csv` / `validation_result.csv` 聚焦 `pass_qc` 细胞
- 当前 sample-root 主线输出 `cluster` 与 `azimuth_cell_type` 对应的 query-native RNA UMAP
  - 另会输出 `qc_overview.png`，按每个 metric 画出样本实际使用的 dynamic MAD 阈值；其中 `n_counts` / `n_genes` 面板在 `log10(x + 1)` 空间展示 median、raw MAD cutoff 与 final applied cutoff，`pct_mt` / `pct_ribo` 面板在原值空间展示同类信息
  - 另会输出 `umap_rna_pbmcref_highlight.png`
  - 另会输出 `umap_rna_cima_l1.png`，且该图严格使用 `azimuth_cima_l1`，不回退到 `cima_l1`，并省略 `Unknown` 标签

## 当前 baseline-only tuning 流程
1. 固定 candidate。
   - 当前只保留一个 candidate：`baseline__baseline__baseline`。

2. 为该 candidate 构造运行配置。
   - 当前通过 `default_tuning_presets()` 取 `baseline` preset，并用 `merge_cli_overrides(...)` 覆盖到 base config。

3. 执行真实单样本主线。
   - 当前顺序为：
     - `read_sample_input(...)`
     - `compute_qc_metrics(...)`
     - `run_doublet_detection(...)`
     - `apply_qc_filters(...)`
     - `run_embedding(...)`
     - `annotate_with_all_versions(..., methods=['azimuth'])`
     - `write_sample_outputs(...)`
   - 当前 reference 目录与主线保持一致，来自 `resolve_data_root(ROOT) / 'reference'`。

4. 计算单 candidate 分数。
   - `qc_score`：当前按 `n_cells_pass_qc / n_cells_total` 计算并 clamp 到 `[0,1]`
   - `annotation_score`：当前按 `confidence_mean * (1 - low_confidence_fraction)` 计算；若方法状态不是 `ok`，则为 `0.0`
   - `embedding_score`：当前按 `separation_score - fragmentation_penalty` 计算并 clamp 到 `[0,1]`
   - `total_score`：当前为三项算术平均

5. 写出 tuning 审计。
   - 当前会写出 `candidates.csv`
   - 若该 candidate 成功，则写出 `selection_summary.json` 和 `selected_params.json`
   - 若该 candidate 失败，当前会保留 `candidates.csv`，但不会伪造 selection summary；此时 orchestration 会报错并停止

## 当前 RNA 输出字段
- `metadata.csv` / `metadata_qc.csv` 当前会写回：
  - `pass_qc`
  - `fails_count_floor`
  - `fails_gene_floor`
  - `fails_mt_ceiling`
  - `fails_ribo_ceiling`
  - `fails_doublet`
  - `doublet_score`
  - `is_doublet`
  - `cluster`
  - `umap_1`
  - `umap_2`
  - `azimuth_cell_type`
  - `azimuth_score`
  - `azimuth_score_margin`
  - `azimuth_low_confidence`
  - `cima_l1`
  - `cima_l2`
  - `cima_l1_masked`
  - `cima_l1_low_confidence`
  - `cima_l1_score`
  - `cima_l1_score_margin`
  - `cima_l2_score`
  - `cima_l2_score_margin`

## 当前 QC / validation 审计输出
- `qc_summary.csv` 当前至少包含：
  - `sample_id`
  - `gse`
  - `n_cells_total`
  - `n_cells_pass_qc`
  - `n_cells_fail_qc`
  - `pass_qc_fraction`
  - `azimuth_status`
  - `azimuth_detail`
  - `azimuth_score_mean`
  - `azimuth_score_margin_mean`
  - `azimuth_low_confidence_fraction`
  - `annotation_score`
  - `qc_threshold_method`
  - `final_min_counts`
  - `final_min_genes`
  - `final_max_pct_mt`
  - `final_max_pct_ribo`
- `qc_thresholds.json` 当前记录样本级动态阈值审计，包括 method、每个 metric 的 MAD 统计、raw threshold、final threshold 和 guardrail 应用情况
- `validation_result.csv` 当前至少包含：
  - `completion`
  - `metadata_all_cells`
  - `metadata_qc_pass_qc_only`
  - `output_presence:*` 系列检查
  - `annotation_status:{method}` 系列检查（若 `annotation_method_status` 存在）

## 非主线但保留的流程
- `process_single_sample.R`
  - 单样本 scATAC 主入口
- `organize_tea_seq_outputs.py`
  - TEA-seq accepted 输出整理与 QC 审计
- `render_tea_seq_cima_cluster_umap.py`
  - TEA-seq 的 cluster-level CIMA L1 辅助 UMAP
- `render_tea_seq_adt_broad_umap.py`
  - TEA-seq 的 ADT broad-label 辅助 UMAP

这些脚本当前继续保留，但不属于 `only_rna` 分支第一阶段验收主线。

## 文档同步约定
- 只要修改了 RNA 主线目标、输入发现规则、QC 逻辑、注释层级、输出规范或命令接口，就必须同步更新 `AGENTS.md`。
- `AGENTS.md` 只记录当前分支已确认的真实行为；不要把未来计划写成既成事实。
