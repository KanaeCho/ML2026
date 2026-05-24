# scripts 目录说明

本目录保存当前项目的主要分析脚本和辅助工具。当前仓库保留既有代码结构，不在 Git 整理阶段大规模移动脚本；建议通过本文档识别主线入口、分支入口和历史/辅助脚本。

## 总入口

- `scripts/process/pipeline.py` 是当前推荐的统一 CLI 入口。
- RNA、共测、longevity 和 product-level 整合命令都通过该文件分发。
- 常用运行方式为 `uv run python scripts/process/pipeline.py <command> ...`。

## RNA 主线

- `scripts/only_rna/` 是当前 Python-first RNA 单样本主线。
- 主线功能包括样本发现、矩阵读取、QC、doublet 检测、embedding、Azimuth `pbmcref` 注释、UMAP 渲染和输出审计。
- `scripts/only_rna/cli.py` 实现 `discover-rna`、`run-rna-sample`、`run-rna-gse`、`tune-rna-sample`、`tune-rna-gse` 和 `rna-status` 等命令。
- `scripts/only_rna/default_config.yaml` 保存 RNA 主线默认参数。
- `scripts/only_rna/tuning_orchestrator.py` 当前只保留 baseline-only tuning 工作流。

## ATAC 主线

- `scripts/process/process_single_sample.R` 是当前 only_atac 单样本主线脚本。
- 主线功能包括 fragment/barcode 输入、peak-by-cell matrix 构建、QC、doublet 检测、CIMA ATAC 注释、UMAP 输出和审计文件写出。
- ATAC 主线通常由 `scripts/process/pipeline.py`、`scripts/co/cli.py` 或 `scripts/longevity/cli.py` 间接调用。
- `scripts/process/run_single_sample_umap.R`、`scripts/process/regenerate_qc_overview.R` 和 `scripts/process/summarize_matrix_lite_qc.R` 是 ATAC/R 输出相关辅助脚本。

## 共测分支

- `scripts/co/` 保存 paired RNA/ATAC 共测分支逻辑。
- `scripts/co/cli.py` 读取 `data_root/reference/co.xlsx`，并将 RNA 样本转发到 `scripts.only_rna.cli`，将 ATAC 样本包装后复用 only_atac 主线。
- `scripts/co/process_co_atac_sample.R` 是旧的独立 co-ATAC 脚本，当前不作为 `co-run-atac-*` 的主路由。

## Longevity 分支

- `scripts/longevity/` 保存独立 longevity 数据通道逻辑。
- `scripts/longevity/cli.py` 实现 longevity RNA ingest、ATAC barcode 预处理、ATAC 单样本运行、参数对照汇总和发布等命令。
- `scripts/longevity/preprocess_atac_barcodes_archr.R` 使用 ArchR 生成 longevity ATAC 的 filtered barcode 文件。
- longevity 分支不读取 `datasets.xlsx`、`atac.xlsx` 或 `co.xlsx`，而是扫描 `data_root/raw/longevity/`。

## Product-Level 整合

- `scripts/process/organize_integrated_products.py` 负责发现已完成样本输出、组织四类 product 目录、写 manifest、调用整合与 panel 渲染。
- `scripts/process/integrate_product_embeddings.py` 负责 RNA/ATAC product-level 低维 embedding 整合，不合并全量 count matrix 或 peak matrix。
- `scripts/process/render_product_umap_panels.py` 负责基于 product-level 坐标生成 CIMA L1/L2、cluster、GSE 和 sample panel 图。
- `scripts/process/append_atac_product_covariates.py` 是 ATAC product metadata 补充工具。

## 数据下载与候选辅助脚本

- `scripts/process/download_from_datasets.py` 用于按参考表下载或整理公开数据，实际可用性依赖外部网络和数据源状态。
- `scripts/process/fragment_top_barcodes.py`、`scripts/process/gse214546_barcode_candidates.py` 和 `scripts/process/gse282769_barcode_candidates.py` 是 barcode 候选或数据集特异辅助脚本。
- 这类脚本应在运行前先检查参数、输入路径和输出目录，避免覆盖正式结果。

## 历史或探索脚本

- `scripts/Function_scRNA.R` 和 `scripts/process/process_single_rna_sample.R` 属于保留的 R 版 RNA/历史辅助逻辑，当前不是 RNA 第一阶段主线入口。
- `scripts/process/render_covid19_minidata_pooled_umap.py` 属于特定数据或探索性可视化脚本，不作为通用主线入口。
- 对不确定用途的脚本，建议先查阅 `docs/` 和代码调用关系，再决定是否纳入正式工作流。

## 维护约定

- 修改 RNA 主线目标、输入发现规则、QC 逻辑、注释层级、输出规范或命令接口时，需要同步更新 `AGENTS.md` 和 `docs/` 中对应说明。
- 新增脚本应优先通过 `scripts/process/pipeline.py` 或分支 CLI 暴露稳定入口。
- 不建议在脚本中写入个人机器绝对路径；应使用 CLI 参数、配置文件或 `ML2026_DATA_ROOT` 等环境变量。
