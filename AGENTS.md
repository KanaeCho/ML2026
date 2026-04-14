# ML2026 项目说明

## 项目定位
当前分支聚焦 **单样本 scRNA-seq** 主线。

本分支第一阶段验收目标：

1. 对 `datasets.xlsx` 中已筛选、且原始计数矩阵已下载到本地的数据，逐样本完成 RNA `QC + 聚类 + CIMA RNA L1/L2 注释 + UMAP 可视化`。
2. 重点保证 `L1 + L2` 的 celltype UMAP 质量可读，并同时输出可复核的简单审计指标。

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
  - 成功输出 `metadata.csv`、`metadata_qc.csv`、`qc_summary.csv`、`validation_result.csv`、四张 RNA UMAP 图、矩阵导出、`run_status.json` 和 `GSM5102900.h5ad`
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
  - `validation_result.csv`
  - `{sample_id}.h5ad`
  - `matrix/matrix.mtx`
  - `matrix/barcodes.tsv.gz`
  - `matrix/features.tsv.gz`
  - `umap_rna_clusters.png`
  - `umap_rna_cima_cell_type_l1.png`
  - `umap_rna_cima_cell_type_l2.png`
  - `umap_rna_cima_cell_type_l1_masked.png`
  - `run_status.json`
- 仓库中旧有的 ATAC / TEA-seq 输出仍位于：
  - `output/{GSE}/{GSM}/`
  - `output/1.only_atac/`
  - `output/GSE214546/qc_audit/`
  但这些目录不属于本分支 RNA 第一阶段验收主线。

## 当前主线脚本

### `scripts/only_rna/`
- 当前 RNA 主线已切换为 Python-first 子系统。
- 当前已实现模块：
  - `config.py`：YAML 默认配置与 CLI override 合并
  - `discovery.py`：数据根解析、`datasets.xlsx` 可见 `scRNA` 行发现、样本本地布局发现
  - `read_inputs.py`：triplet / `.h5` / `matrix.tar.gz` / GSE-shared triplet 读取为 `AnnData`
  - `qc.py`：`n_counts`、`n_genes`、`pct_mt`、`pct_ribo` 计算与 QC fail flags / `pass_qc`
  - `doublet.py`：Python doublet 结果归一化与 `scanpy.pp.scrublet` 路径
  - `embedding.py`：仅对 `pass_qc` 细胞执行 embedding / clustering / UMAP，并写回 `cluster`、`umap_1`、`umap_2`
  - `annotation.py`：最小 CIMA 资产加载、`pass_qc` 细胞 `L1/L2` 注释、低置信 masking
  - `plotting.py`：按配置输出可读的 categorical UMAP
  - `outputs.py`：写出 `.h5ad`、metadata/QC/validation CSV、矩阵导出、UMAP 图
  - `cli.py`：`discover-rna` / `run-rna-sample` / `run-rna-gse` / `rna-status` 的 Python 主线实现

### `scripts/process/pipeline.py`
- 当前支持的 RNA 命令：
  - `discover-rna`
  - `run-rna-sample`
  - `run-rna-gse`
  - `rna-status`
- 当前行为：
  - 保持上述命令名稳定
  - RNA 命令分发到 `scripts.only_rna.cli`
  - `run-rna-sample` 显式拒绝 `gse_shared` 目标
  - `run-rna-gse` 会隐式包含受支持的 GSE-shared triplet
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
   - 当前已验证的默认 YAML 阈值：
     - `min_counts = 500`
     - `min_genes = 300`
     - `max_pct_mt = 20.0`
     - `max_pct_ribo = 60.0`
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

6. 基于 CIMA RNA 参考做 L1/L2 注释。
   - 使用 `cima_rna_reference_pca_features.tsv.gz` 中的 `gene_mean/gene_std/pc_dim_*`
   - 先预测 L1，再在允许的 L1 子集中预测 L2
   - 输出 `cima_l1_score`、`cima_l1_score_margin`、`cima_l2_score`、`cima_l2_score_margin`
   - 当前低置信定义已验证为：`cima_l1_score < 0.35` 或 `cima_l1_score_margin < 0.05`
   - 低置信细胞在 `cima_l1_masked` 中标成 `Unknown`

7. 写回 metadata 和 UMAP 图。
   - `metadata.csv` 保留全细胞信息
   - `metadata_qc.csv` / `validation_result.csv` 聚焦 `pass_qc` 细胞
   - 输出 `cluster`、`cima_l1`、`cima_l2`、`cima_l1_masked` 对应的 query-native RNA UMAP

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
  - `cima_l1`
  - `cima_l2`
  - `cima_l1_masked`
  - `cima_l1_low_confidence`
  - `cima_l1_score`
  - `cima_l1_score_margin`
  - `cima_l2_score`
  - `cima_l2_score_margin`

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
