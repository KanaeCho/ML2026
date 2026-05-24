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
  - 成功输出 `metadata.csv`、`metadata_qc.csv`、`qc_summary.csv`、`validation_result.csv`、`qc_overview.png`、`umap_rna_pbmcref_vs_cima_l1.png`、`run_status.json` 和 `GSM5102900.h5ad`
- 当前主线默认 QC 已切换到更严格的一档：`counts_lower_nmads=2.5`、`genes_lower_nmads=2.5`、`pct_mt_upper_nmads=2.5`、`pct_ribo_upper_nmads=3.5`
- 当前代码中的 baseline-only tuning 路径已经落地：
  - `pipeline.py` / `scripts.only_rna.cli` 已支持 `tune-rna-sample` 与 `tune-rna-gse`
  - tuning 当前只执行固定 candidate：`baseline__baseline__baseline`
  - 该 candidate 会真实执行单样本主线（读取矩阵、QC、doublet、embedding、Azimuth、样本输出），不是占位 stub
- `GSE226039` 当前按文件名只保留 `PBMC` 样本参与 RNA 发现与运行，不把 Ileum / Rectum 组织一起纳入本分支主线。
- 当前新增 `scripts/co/` 共测子系统；`co-run-atac-*` 会把共测 ATAC 样本包装后复用 `scripts/process/process_single_sample.R` 的 only_atac 主线，并输出到 `output/co/atac/`。
- `co2` / `GSE206284` 已从当前工作区清理：`data/reference/co2_sample_manifest.csv`、`data/reference/co2_rna_atac_pairing.csv`、`data/raw/GSE206284/` 以及对应 RNA/ATAC 输出目录均不再保留。
- 当前 `co1` 共测数据集 `7555405` 已接入 `scripts/co/`：样本发现读取 `data_root/reference/co.xlsx`，原始样本目录与 pipeline 发现 / 输出均使用英文 `sample`（如 `donorA_Day0`）。
- 当前 `co1` / `7555405` 共测分支任务已完成：24 个样本的 RNA 输出位于 `output/co/rna/7555405/`，24 个样本的 ATAC 输出位于 `output/co/atac/7555405/`，co-ATAC 状态均为 `success` 且 `outputs_complete=true`。
- 当前 `GSE224198` 已按共测目录结构接入 `scripts/co/`：8 个 paired RNA/ATAC 样本写入 `data/reference/co.xlsx`，样本目录为 `data/raw/GSE224198/{sample}/RNA` 与 `data/raw/GSE224198/{sample}/ATAC`；GEO 公开元数据提供 `SLE` 与 `healthy` donor group，但未提供年龄，因此 `age` 留空。
- 当前新增独立 `longevity` 额外数据通道，不属于 only_rna、only_atac GSE/GSM 自动发现，也不属于 `co` 共测分支；原始数据位于 `data/raw/longevity/rna/` 与 `data/raw/longevity/atac/`，发现和运行通过 `pipeline.py longevity-*` 命令完成，不依赖 `datasets.xlsx`、`atac.xlsx` 或 `co.xlsx`。
- 当前整合分支新增内存受控的 product-level low-dimensional integration 工作流：通过 `uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink` 生成 `output/1.only_atac/`、`output/2.only_rna/`、`output/3.co_atac/`、`output/4.co_rna/`。该工作流严格不合并全量 count matrix；RNA 逐样本按全基因 library size 做 `normalize_total(1e4)`，对 CIMA feature 做 `log1p` 后投影到 CIMA RNA PCA compact feature 空间，ATAC 逐样本投影到 CIMA ATAC LSI compact feature 空间，只在 product 层合并低维 `float32` embedding。RNA product 默认对低维 embedding 执行 Harmony 校正后，再用 Scanpy neighbors + UMAP + Leiden 写回真实 `integrated_umap_1/2` 与 `integrated_cluster`；ATAC product 默认仍使用 BBKNN + Scanpy UMAP + Scanpy Leiden。only_rna 默认只让 `integrated_cima_l1_score >= 0.5` 的高置信 RNA 细胞参与 product-level UMAP / Leiden，以减少低置信 CIMA 投影细胞造成的桥状结构；低置信细胞仍保留在 product metadata 中，并通过 `integration_included=false` 与 `integration_exclusion_reason=low_integrated_cima_l1_score` 审计。若输出挂载点不允许 symlink，则写入 `SOURCE_OUTPUT_DIR.txt` pointer 指向原始样本目录，避免复制大矩阵导致内存/磁盘压力。

## 目录结构

### `data/`
- 当前仓库内通常不存在真实 `data/` 目录；运行时数据默认来自外部数据根。
- 共享参考位于 `data_root/reference/`，原始 RNA 输入位于 `data_root/raw/`。

### `data_root/reference/`
- `datasets.xlsx`
  - RNA 样本发现基于 Excel 第一张表中当前可见、且 `assay=scRNA` 的行。
- `longevity/atac_barcodes/`
  - 独立 longevity ATAC 的 ArchR barcode 预处理输出目录；每个样本子目录至少包含 `filtered_barcodes.tsv.gz`、`barcode_qc.csv.gz`、`summary.json`。`longevity-run-atac-*` 会先确保该目录中存在 ready 状态的 `filtered_barcodes.tsv.gz`，若缺失则先调用 ArchR 预处理生成 barcode，再把该文件传给 only_atac 主线的 `--barcode-file`。
- `longevity/atac_barcode_overrides.csv`
  - 可选的 longevity ATAC barcode 预处理阈值覆盖表，列为 `sample_id,min_fragments,max_barcodes,min_tss,reason`。
- `co.xlsx`
  - 共测数据集运行清单；当前只保留 `sample`、`dataset`、`age`、`health`、`donor` 五列，`scripts/co/` 运行时用 `sample` 映射到原始目录与输出样本名。
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
- 独立 longevity 额外数据不按 GSE/GSM 发现：`data_root/raw/longevity/rna/` 当前包含处理后的 `scRNA_v3_L2.h5ad` atlas；`data_root/raw/longevity/atac/` 当前包含 10 个 `W*_fragments.tsv.gz` 单样本 ATAC fragment 文件。

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
  - `umap_rna_pbmcref_vs_cima_l1.png`
  - `umap_rna_pbmcref_highlight.png`
  - `umap_rna_cima_l1.png`
  - `run_status.json`
  - `{sample_id}.h5ad` 只保留 pass-QC 且 `final_celltype` 为已知 5 类 RNA 大类的细胞和表达矩阵；`.obs` 至少包含 `sample`、`dataset`、`age`、`health`、`donor`、`final_celltype`、`final_celltype_mapping`、`pbmcref_celltype`、`azimuth_cima_l1_raw`、`azimuth_cima_l1`、`azimuth_cell_type_l2_raw`、`umap_1/2`。RNA `final_celltype` 统一使用 5 类：`CD4_T`、`CD8_T`、`B`、`Myeloid`、`NK`；`azimuth_cima_l1` 与 `final_celltype` 保持一致，原始 CIMA L1 审计标签保留在 `azimuth_cima_l1_raw`。映射时 `ILC -> NK`、`unconvensional_T -> CD8_T`，且若 pbmcref L2 包含 `NK` 则优先映射为 `NK`。`metadata.csv` 仍保留全细胞审计，`metadata_qc.csv` / `validation_result.csv` 与 h5ad 的已知 final celltype 细胞集合对齐。
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
- 仓库中 ATAC / TEA-seq 输出当前位于：
  - ATAC 单样本默认输出目录：`output/atac/{GSE}/{GSM}/`
  - `output/1.only_atac/`（旧 accepted integration 输出已由当前 only_atac product 接管）
  - `output/GSE214546/qc_audit/`
  但这些目录不属于本分支 RNA 第一阶段验收主线。
- 共测 ATAC 输出目录：`output/co/atac/{GSE}/{GSM}/`
- 对 `7555405`，共测输出目录使用英文 sample id，例如 RNA `output/co/rna/7555405/donorA_Day0/` 与 ATAC `output/co/atac/7555405/donorA_Day0/`。
- 当前 co-ATAC 样本至少输出：
  - `metadata.csv`
  - `metadata_qc.csv`
  - `qc_summary.csv`
  - `validation_result.csv`
  - `qc_overview.png`
  - `umap_cima_cell_type_l1.png`
  - `umap_cima_cell_type_l2.png`
  - `{GSM}.h5ad`
  - `run_status.json`
  - `{GSM}.h5ad` 只保留 pass-QC 细胞和 peak-by-cell 矩阵；`.obs` 至少包含 `sample`、`dataset`、`age`、`health`、`donor`、`final_celltype`、`cima_cell_type_l1`、`cima_cell_type_l1_masked`、`cima_cell_type_l1_cluster_consensus`、`cima_cell_type_l2/3/4`、`cima_l1_low_confidence`、`cima_l1_cluster_purity`、`cima_l4_score`、`cima_l4_score_margin`、`umap_atac_1/2`、`cima_ref_umap_1/2`。ATAC `final_celltype` 使用 `cima_cell_type_l1`；单样本最终不保留 `matrix/` 目录和 `{GSM}_seurat_qc.rds`。
- 当前 integrated product 输出目录：
  - `output/1.only_atac/`：聚合 only_atac 样本输出 `output/atac/{GSE}/{GSM}/`，排除 co2 `GSE206284`；legacy ATAC cell metadata 来源为 `validation_result.csv`。
  - `output/2.only_rna/`：聚合 `output/rna/{GSE}/{sample_id}/` 下已完成的 only_rna RNA 样本，跳过 tuning candidate 嵌套输出，并在 product 层排除 co2 `GSE206284`。
  - `output/3.co_atac/`：聚合 `output/co/atac/7555405/{sample_id}/` 的 co1 ATAC 样本。
  - `output/4.co_rna/`：聚合 `output/co/rna/7555405/{sample_id}/` 的 co1 RNA 样本。
- 每个 integrated product 目录至少包含：
  - `product_status.json`
  - `manifests/samples.csv`
  - `manifests/cells_metadata.csv`
  - `manifests/output_files.csv`
  - `qc/sample_qc_summary.csv`
  - `qc/validation_summary.csv`
  - `{product}.h5ad`（metadata-only 对象，不含表达/peak count matrix）
  - `integration/integration_summary.json`
  - `integration/integration_metrics.csv`
  - `integration/sample_mixing_summary.csv`
  - `figures/*_cima_l1_panels.png`、`figures/*_cima_l2_panels.png`、`figures/*_integrated_cluster_panels.png`、`figures/*_gse_panels.png`、`figures/*_sample_panels.png`
- integrated product 的 cell metadata 会追加并保留追踪字段：`product`、`branch`、`modality`、`gse`、`sample_id`、`individual_id`、`source_output_dir`、`source_cell_id`、`global_cell_id`、`is_co_sample`、`co_dataset`、`co_dataset_id`。
- integrated product 的 cell metadata 还会写入真实 product-level 整合字段：`integrated_umap_1`、`integrated_umap_2`、`integrated_cluster`、`integration_method`、`integration_feature_space`、`integration_included`、`integration_exclusion_reason`。RNA product 还会基于同一个 CIMA RNA PCA embedding 和 reference centroids 写入 `integrated_cima_l1`、`integrated_cima_l1_score`、`integrated_cima_l2`、`integrated_cima_l2_score`，用于检查与整合空间一致的 CIMA label。only_rna 中低于默认 CIMA L1 置信阈值的细胞保留 `integrated_cima_l1/l2` 诊断标签，但不写入 product-level `integrated_umap_1/2` 和 `integrated_cluster`。panel 图优先使用 `integrated_umap_1/2`，不再使用 per-sample UMAP 或 reference UMAP 作为主验收坐标，且会自然跳过缺少 integrated 坐标的低置信细胞。
- integrated product 当前不执行 co1 RNA+ATAC 跨模态整合，不做联合 embedding，不做 paired-cell matching，也不做 RNA-to-ATAC label transfer。

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
  - `outputs.py`：写出 `.h5ad`、metadata/QC/validation CSV、UMAP 图；当前也负责 tuning 选择产物写出。样本级 `.h5ad` 保留 pass-QC 且 `final_celltype != Unknown` 的细胞矩阵，不再额外导出 `matrix/` 目录。
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
- 当前支持的共测命令：
  - `co-discover`
  - `co-status`
  - `co-run-rna-sample`
  - `co-run-rna-gse`
  - `co-run-atac-sample`
  - `co-run-atac-gse`
- 当前支持的整合产品命令：
  - `organize-products`
- 当前支持的独立 longevity 命令：
  - `longevity-discover`
  - `longevity-status`
  - `longevity-ingest-rna`
  - `longevity-preprocess-atac-barcodes`
  - `longevity-atac-barcode-status`
  - `longevity-run-atac-sample`
  - `longevity-run-atac-all`
  - `longevity-run-atac-param-contrast`
  - `longevity-run-atac-custom-param-contrast`
  - `longevity-summarize-atac-param-contrast`
  - `longevity-publish-atac-param-contrast`
- 共测命令当前行为：
- `co-discover` 当前基于 `data/reference/co.xlsx` 列出 co1 RNA 与 ATAC 样本，并显示 `individual_id`。
  - `co-run-rna-*` 复用 `scripts.only_rna.cli` 的 RNA 单样本流程，默认输出根为 `output/co/rna/`
  - `co-run-atac-*` 复用 `scripts/process/process_single_sample.R` 的 only_atac 主线，默认输出根为 `output/co/atac/`；`co-run-atac-gse` 支持 `--jobs` 并行处理样本
  - `co-status` 当前只检查 co-ATAC only_atac 输出完整性
- longevity 命令当前行为：
  - `longevity-discover` 直接扫描 `data_root/raw/longevity/rna/*.h5ad` 与 `data_root/raw/longevity/atac/*_fragments.tsv.gz`，不读取 Excel 清单。
  - `longevity-preprocess-atac-barcodes` 调用 `scripts/longevity/preprocess_atac_barcodes_archr.R` 生成 ArchR-based filtered barcode 文件，默认写到 `data_root/reference/longevity/atac_barcodes/{sample_id}/`；默认阈值为 `min_fragments=200`、`max_barcodes=20000`、`min_tss=2.5`，可通过 CLI 或 `data_root/reference/longevity/atac_barcode_overrides.csv` 覆盖。
    - barcode 预处理默认按 fragment 数排序选择 top barcodes；实验性 `--rank-by tss_then_fragments` 可按 ArchR TSS enrichment 优先排序，用于高 FRiP 目标下的单样本参数探索。
  - `longevity-atac-barcode-status` 汇总每个 longevity ATAC 样本的预处理 barcode 文件状态、barcode 数量以及 QC/summary 文件存在性。
  - `longevity-run-atac-*` 将独立 longevity ATAC fragment 文件通过显式 `--fragment-file` 传给 `scripts/process/process_single_sample.R`，默认输出根为 `output/longevity/atac/longevity/{sample_id}/`；运行 ATAC 主流程前会先检查 `data_root/reference/longevity/atac_barcodes/{sample_id}/filtered_barcodes.tsv.gz` 是否 ready，若未 ready 则先调用 ArchR barcode 预处理，随后把生成的 `filtered_barcodes.tsv.gz` 显式传给 only_atac 主线的 `--barcode-file`；longevity ATAC 支持可选传递 `--umap-min-dist` 覆盖 query-native `RunUMAP(min.dist=...)`，但默认不传该参数，保持 Seurat/uwot 默认 UMAP 行为。
  - `longevity-run-atac-param-contrast` 在隔离目录中按预定义候选参数运行 longevity ATAC barcode 预处理与 only_atac 主线；barcode 输出写到 `data_root/reference/longevity/atac_barcodes_param_contrast/{candidate}/{sample_id}/`，样本输出写到 `output/longevity_param_contrast/{candidate}/atac/longevity/{sample_id}/`，不会写正式 `output/longevity/atac/longevity/`。
  - `longevity-run-atac-custom-param-contrast` 在前台对单个 longevity ATAC 样本运行一个自定义 barcode 参数 candidate，用于根据上一次结果逐步调整参数；输出目录规则同 `longevity-run-atac-param-contrast`。
  - `longevity-summarize-atac-param-contrast` 扫描 `output/longevity_param_contrast/{candidate}/atac/longevity/{sample_id}/qc_summary.csv`，汇总各样本 candidate 指标，并按 `median_FRiP >= 0.6` 与 `pass_qc >= --min-pass-qc` 的硬目标先过滤，再按 `cima_l1_low_confidence_frac`、`median_cima_l4_score`、`median_FRiP`、`qc_rate`、`pass_qc` 的加权分数选择每个样本的 best candidate；`--write` 会写出 `param_contrast_summary_all_samples.csv` 与 `param_contrast_selected_samples.csv`。
  - `longevity-publish-atac-param-contrast` 读取 `param_contrast_selected_samples.csv`，把每个样本选中的 candidate 输出目录复制到正式 `output/longevity/atac/longevity/{sample_id}/`；默认不覆盖已有目录，需显式 `--force`，并支持 `--dry-run`。
  - `longevity-ingest-rna` 读取 processed RNA atlas，按 `L1_annotation` / `L2_annotation_new` 写入 5 类 `final_celltype`，默认输出到 `output/longevity/rna/longevity/{sample_id}/`，且只写修改后的 `.h5ad` 与 `celltype_original_vs_final.csv` 前后对比表；longevity RNA 不执行 only_rna 单样本 QC/Azimuth 主线，也不按 Unknown/空过滤 barcode。
- `organize-products` 当前行为：
  - 默认命令为 `uv run python scripts/process/pipeline.py organize-products --products all --copy-mode symlink`
  - `--products` 支持 `only_atac`、`only_rna`、`co_atac`、`co_rna`、`all`
  - `--copy-mode` 支持 `symlink` 与 `copy`，默认 `symlink`；若挂载点禁止 symlink，会自动降级为 `SOURCE_OUTPUT_DIR.txt` pointer，不复制大文件
  - 默认执行真实 product-level 低维整合，并写入 `integrated_umap_1/2` 与 `integrated_cluster`
  - 默认整合方法按 product 区分：RNA product 默认 `--integration-method harmony`，ATAC product 默认 `--integration-method bbknn`；`--integration-method scanpy_neighbors` 可用于普通 Scanpy neighbors fallback
  - `--skip-integration` 才会跳过整合；此时只生成 metadata-level product，不作为当前整合验收通过状态
  - `--skip-figures` 只跳过 panel 渲染，不跳过 integration
  - `--integration-n-components`、`--integration-max-umap-fit-cells`、`--integration-clusters`、`--integration-batch-key`、`--bbknn-neighbors-within-batch`、`--bbknn-trim`、`--leiden-resolution`、`--rna-min-cima-l1-score` 控制低维整合参数
  - 为控制内存，该命令按 product 串行处理，不合并表达矩阵或 peak matrix；RNA/ATAC 均按样本逐个投影为低维 `float32` embedding，product 层只合并低维 embedding

### `scripts/process/organize_integrated_products.py`
- 当前整合产品组织器，负责发现四类已注释样本输出、校验必要文件、写出 product manifest、调用低维整合引擎、生成 metadata-only `.h5ad` 和调用 panel renderer。
- 当前真实生成状态：`output/1.only_atac/` 为 70/70 样本、409532 个 QC 后细胞；`output/2.only_rna/` 为 35/35 样本、386313 个 QC 后细胞；`output/3.co_atac/` 为 24/24 样本、69813 个 QC 后细胞；`output/4.co_rna/` 为 24/24 样本、98670 个 QC 后细胞。

### `scripts/process/integrate_product_embeddings.py`
- 当前 product-level 低维整合引擎。
- RNA product 使用样本 `.h5ad` 与 `cima_rna_reference_pca_features.tsv.gz`，逐样本对 pass-QC 细胞按全基因 count 总量执行 `normalize_total(1e4)`，对 CIMA feature 子集执行 `log1p`，再用 CIMA reference feature model 的 `gene_mean/gene_std` 标准化后投影到 CIMA RNA PCA compact feature 空间。only_rna 默认 `--integration-method harmony`、`--integration-batch-key gse`、`--rna-min-cima-l1-score 0.5`；co_rna 默认 `--integration-method harmony`、`--integration-batch-key sample_id`，不启用 only_rna 的 CIMA L1 置信过滤；ATAC product 默认仍使用 `--integration-method bbknn`、`sample_id` 批次和原 BBKNN 参数。
- ATAC product 使用样本 `.h5ad` 内的 peak-by-cell 矩阵与 `cima_atac_reference_lsi_features.tsv.gz`，逐样本投影到 CIMA ATAC LSI compact feature 空间。
- product 层对低维 embedding 执行所选整合后写出 Scanpy UMAP 和 Scanpy Leiden cluster；RNA product 当前默认是 Harmony 校正后的 `X_harmony` 邻接图，ATAC 当前默认是 BBKNN graph。整合输出包括 `integration/integration_summary.json`、`integration/integration_metrics.csv`、`integration/sample_mixing_summary.csv`。
- RNA product 会额外用 CIMA RNA L1/L2 centroid 对低维 embedding 做余弦最近质心赋值，写入 `integrated_cima_l1/l2`。RNA 单样本主注释仍来自 Azimuth `pbmcref`；`integrated_cima_l1/l2` 是 CIMA projection-space nearest-centroid 诊断标签，不替代 pbmcref 主注释。当前诊断显示 only_rna 中低 `integrated_cima_l1_score` 细胞会造成明显桥状结构，因此 only_rna 默认仅用高置信投影细胞拟合 product-level UMAP / Leiden，低置信细胞不进入整合图但保留审计字段。
- 当前真实整合指标：`only_atac` 混合邻居比例 0.9988、15 clusters；`only_rna` 共有 386313 个 QC 后细胞，其中 288324 个高置信细胞参与 Harmony product-level 整合、97989 个低置信细胞保留 metadata 但不参与 UMAP/Leiden，混合邻居比例 0.9750、14 clusters、`integrated_cima_l1` cluster purity 0.8271、`integrated_cima_l2` cluster purity 0.5911；`co_atac` 混合邻居比例 1.0000、12 clusters；`co_rna` 98670 个 QC 后细胞全部参与 Harmony product-level 整合，混合邻居比例 1.0000、13 clusters、`integrated_cima_l1` cluster purity 0.8013、`integrated_cima_l2` cluster purity 0.5803。四个 product 的 `integrated_equals_original_umap` 均为 `false`。

### `scripts/process/render_product_umap_panels.py`
- 当前 integrated product panel 图渲染器沿用旧 accepted integration 的视觉规范：灰色背景点、按类别高亮、panel 标题包含 `n=`、rasterized scatter、CIMA L1 基础配色和 L2 阶梯色。
- 坐标列优先级为 `integrated_umap_1/2`、`cima_ref_umap_1/2`、`umap_atac_1/2`、`umap_1/2`；当前主验收要求使用 `integrated_umap_1/2`。若缺坐标列则在 `product_status.json` 中记录跳过原因，不伪造 UMAP。
- 当前每个 product 额外输出 `*_integrated_cluster_panels.png`，用于检查真实整合后的 cluster 结构。

### `scripts/co/`
- 当前共测子系统已实现：
  - `cli.py`：共测 manifest 发现、RNA 路由复用、ATAC only_atac 包装调度与状态检查
  - `process_co_atac_sample.R`：旧独立 co-ATAC 质控与 query-native 可视化脚本，当前不再作为 `co-run-atac-*` 主路由
- 当前共测数据接入方式：
  - 从 `data/reference/co.xlsx` 读取 `sample`、`dataset`、`age`、`health`、`donor`
  - RNA 输入来自 `data/raw/{dataset}/{sample}/RNA/matrix.mtx(.gz) + barcodes.tsv(.gz) + features.tsv(.gz)`
  - ATAC 输入来自 `data/raw/{dataset}/{sample}/ATAC/fragments.tsv.gz`，若同目录存在 `filtered_barcodes.tsv.gz`、`filtered_metadata.csv.gz` 或 `singlecell.csv.gz`，`co-run-atac-*` 会显式传给 only_atac 主线用于初始 barcode 选择
- 当前 co-ATAC 处理流程：
  - `co-run-atac-*` 通过显式 `--fragment-file` 把共测 fragment 文件传给 only_atac 主线，并优先显式传递共测样本目录中的 barcode 辅助文件，避免嵌套目录样本误回退到过宽的 fragment-count barcode inference
  - only_atac 从 `fragments.tsv.gz` 读取 ATAC 片段；初始 barcode 优先级为显式/发现到的 `filtered_barcodes`、`filtered_metadata`、`singlecell`，若均不存在才按主线 barcode inference 从 fragment counts 推断初始 barcode 集合
  - 使用 `peak.bed` 构建 peak-by-cell matrix，计算 QC 指标，使用 `scDblFinder` 做 doublet 检测，并执行 only_atac 既有 QC / CIMA ATAC L1/L2 注释 / reference-space UMAP 输出流程

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
  - `azimuth_cima_l1_raw`
  - `azimuth_cima_l1`
  - `final_celltype`
  - `final_celltype_mapping`
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
  - `n_cells_final_output`
  - `n_cells_unknown_final_celltype_removed`
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
