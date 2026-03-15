# ML2026 项目说明

## 项目定位
本项目当前只做两件事：

1. 对单样本 scATAC-seq 数据进行质控，输出 QC 图，并生成可用于后续分析的矩阵结果。
2. 将所有样本的矩阵整合起来，输出整合后的可视化结果，用于判断整合质量和批次效应情况。

补充说明：
- 当前 R/Signac 流程内部主要使用 peak×cell 计数矩阵。
- 如果后续分析必须统一成 cell×peak 方向，应在输出规范中明确写清并统一处理，不要在不同脚本里混用。

## 当前任务
1. 固定单样本 QC 流程。
   - 输入：按 `GSE`、`GSM` 自动发现 fragment 文件和可选的 filtered barcodes 文件。
   - 输出：最终矩阵、最终 metadata、QC 图。
   - 目的：判断单样本数据质量是否可接受。

2. 固定多样本整合流程。
   - 将所有样本矩阵合并。
   - 输出整合后的降维图、聚类图或批次分布图。
   - 用于判断样本之间是否仍存在明显批次效应。

## 目录结构

### `data/`
原始数据和共享参考文件。

#### `data/reference/`
- `datasets.xlsx`: 只读数据集统计表。
- `peak.bed`: 统一 peak 定义文件。
- `peaks.csv`: 带 peak_id 的参考索引文件。

#### `data/raw/`
- 按 GSE 组织原始输入数据。
- 每个样本至少应包含 fragment 文件。
- 如果有 `*_filtered_barcodes.tsv.gz`，优先作为初始细胞集合使用。
- 文件发现依赖命名约定：
  - fragment 文件名以 `GSM` 开头，且包含 `fragments`
  - barcode 文件名以 `GSM` 开头，且包含 `barcodes`

### `output/`
- 运行时生成的结果目录，不是项目的长期真值来源。
- 结果可以删除后重新生成。
- 当前按 `output/{GSE}/{GSM}` 组织。
- 单样本脚本只保留最终产物，不额外堆积中间文件。
- 当前额外包含：
  - `output/integration_merged/`: 28 个样本 hard-QC 后的总合并矩阵与合并 metadata。
  - `output/integration_sketch/`: 用于整合参数探索的 balanced sketch 矩阵与 metadata。
  - `output/integration_sketch_analysis/`: sketch 矩阵上的 LSI / Harmony / UMAP / 聚类与整合质量评估结果。
  - `output/integration_merged_analysis/`: 全量 merged 矩阵上的正式 LSI / Harmony / UMAP / 聚类与整合质量评估结果。
  - `output/hard_qc_review/`: hard-QC 前后评估图、汇总表与整合准备报告。

### `scripts/process/`
当前主要工作目录。

关键文件：
- `GSM8671454.ipynb`: 当前唯一保留的单样本处理 notebook，用于确认流程是否正确。
- `process_single_sample.R`: 从 notebook 抽出的单样本可调用脚本。
- `regenerate_qc_overview.R`: 基于现有 `*_seurat_qc.rds` 重生单样本 `qc_overview.png`，用于只调整总图版式而不重跑整套 QC。
- `apply_integration_hard_qc.py`: 基于现有 `metadata_qc.csv` 和 `matrix/` 对样本做整合前统一硬阈值二次筛选，生成 integration-ready 输出。
- `merge_integration_matrices.py`: 严格检查 28 个样本 `integration_qc/matrix/features.tsv.gz` 是否完全一致且顺序一致，随后将 28 个 peak×cell 矩阵按列合并为一份总矩阵。
- `build_integration_sketch.py`: 从各样本 `integration_qc` 中按样本均衡抽取细胞，构建用于 LSI / Harmony 参数探索的 sketch 矩阵。
- `run_batch_integration.R`: 正式整合脚本，读取 `integration_sketch/` 或 `integration_merged/`，运行 TF-IDF、LSI、Harmony、UMAP、聚类、batch mixing 评估、cluster × sample 热图、LSI-QC 相关性与 marker peak 摘要。
- `review_hard_qc_and_integration_readiness.py`: 汇总 28 个样本 `metadata_qc.csv` 与 `integration_qc/metadata_integration_qc.csv`，验证 hard-QC 结果、输出前后对比图，并生成整合准备报告。
- `pipeline.py`: Python 流程管理入口，负责样本发现、日志和调度。
- `download_from_datasets.py`: 从 `datasets.xlsx` 过滤样本并组织 GEO 下载任务。

## 当前文件分工
- `GSM8671454.ipynb`
  - 用途：先在单一样本上调通流程。
  - 输出要求：图在 notebook 中直接显示即可。
  - 当前作为项目中的基准流程。
- `process_single_sample.R`
  - 用途：将单样本 notebook 逻辑脚本化。
  - 输入：`GSE`、`GSM`。
  - 行为：自动发现 fragment 和 barcode 文件，运行单样本 QC，保存最终产物。
  - 对缺失 `filtered_barcodes` 的样本，当前在脚本内部完成 barcode 候选集合预筛，不再往 `data/raw` 回写新的 `filtered_barcodes` 文件。
  - 输出目录：`output/{GSE}/{GSM}`。
- `regenerate_qc_overview.R`
  - 用途：读取已有的 `GSM*_seurat_qc.rds`，重新生成 `qc_overview.png`。
  - 适用场景：只修改总 QC 图内容或布局，不重新跑 `FeatureMatrix`、QC 指标和 doublet。
- `apply_integration_hard_qc.py`
  - 用途：在现有单样本 `pass_qc` 结果基础上，再做一轮统一硬阈值筛选。
  - 输出：每个样本单独的 `integration_qc/` 目录，用于后续批次整合。
- `merge_integration_matrices.py`
  - 用途：将 28 个样本的 `integration_qc` 矩阵严格对齐后横向合并。
  - 约束：只有所有 `features.tsv.gz` 内容与顺序完全一致时才允许直接合并；否则应先重排对齐，不能盲目 `cbind`。
  - 当前实现采用流式写入，并先写临时文件后再替换正式输出，避免中途退出时留下伪完成的合并矩阵。
  - 输出：`output/integration_merged/` 下的总 `matrix.mtx`、总 `barcodes.tsv.gz`、总 `features.tsv.gz`、合并 metadata 与样本列偏移清单。
- `build_integration_sketch.py`
  - 用途：在不直接加载全量整合对象的前提下，先构建一个按样本均衡抽样的 sketch 输入，用于去批次参数探索。
  - 行为：当前默认每个样本抽取 `1000` 个 hard-QC 后细胞，输出 sketch 版 `matrix.mtx`、`barcodes.tsv.gz`、`features.tsv.gz`、`merged_metadata.csv`、`sample_manifest.csv` 与 `sketch_summary.json`。
  - 输出：`output/integration_sketch/`。
- `run_batch_integration.R`
  - 用途：对 sketch 或全量 merged 输入运行正式批次整合与质量评估。
  - 行为：当前使用 `Signac + Seurat + harmony`，先按 peak 在细胞中的可及数预筛 top peaks（当前默认 `30000`），再运行 TF-IDF、LSI、Harmony、邻居图、聚类、UMAP，并输出 batch mixing 指标、cluster × sample 热图、LSI 与 QC 指标相关性、以及 cluster top accessible peaks 与最近基因摘要。
  - 环境：脚本会按系统自动追加项目本地库路径（Windows 使用 `.r-win-library/`，Linux/WSL 使用 `.r-linux-library/`），用于优先加载项目侧 R 包依赖。
  - 输入：`output/integration_sketch/` 或 `output/integration_merged/`。
  - 输出：`output/integration_sketch_analysis/` 或 `output/integration_merged_analysis/`。
- `review_hard_qc_and_integration_readiness.py`
  - 用途：复核 hard-QC 固定阈值在全队列上的实际效果，并评估是否适合进入去批次阶段。
  - 行为：从现有 `metadata_qc.csv` 重算 hard-QC 掩码，与 `integration_qc/metadata_integration_qc.csv` 做逐样本一致性校验，统计失败原因，并生成 SVG 图与 markdown 报告。
  - 输出：`output/hard_qc_review/`。
- `pipeline.py`
  - 用途：管理单样本脚本执行。
  - 当前支持：样本发现、单样本运行、按 GSE 批量运行、状态查看、下载入口。
  - 调用方式：由 Python 调用 `process_single_sample.R`，并记录日志与状态。
- `download_from_datasets.py`
  - 用途：读取 `data/reference/datasets.xlsx`，按条件过滤样本，解析 GEO supplementary 下载任务。
  - 当前默认过滤：`scATAC` + `fragment`。
  - 当前支持额外按 `GSE`、`GSM` 过滤。
  - 当前支持通过 `--file-kinds` 指定下载的 GEO 文件类别，例如 `fragment`、`barcode`、`singlecell`、`summary`。
  - 下载前先统计过滤结果、本地已存在样本、目标文件和预计空间，再调用 `aria2c` 下载到 `data/raw/{GSE}`。
  - 当前支持网络超时参数，并对 GEO listing 请求做缓存，避免重复抓取同一 GSE 的 series 目录。
  - 当前支持先输出待下载的 GEO 直链清单，再执行实际下载。
  - 本地文件状态当前按 `.aria2` 控制文件和远端大小区分 `complete/partial/missing`；只要 `.aria2` 仍在，就一律视为未完成并在重跑时继续续传。

## 当前实现方向
- 项目环境已安装 `PySide6`，后续将以 Python + Qt 的桌面应用形式管理流程执行。
- Qt 应用负责流程管理、参数组织、日志展示和状态记录。
- 具体计算仍应由稳定的脚本入口执行，不应把全部分析逻辑直接塞进界面事件里。
- 预期职责分离：
  - Python 负责流程管理、样本发现、状态记录、日志管理
  - R 负责单样本 QC、矩阵生成和多样本整合计算
- 当前已落地的脚本化入口就是 `scripts/process/pipeline.py` + `scripts/process/process_single_sample.R`。
- 当前下载入口是 `scripts/process/pipeline.py download`，其内部调用 `scripts/process/download_from_datasets.py`。

## 单样本处理流程
1. 读取样本参数。
   - 输入至少包括 `GSE`、`GSM`。
   - `peak.bed` 和输出根目录由项目约定固定，不作为用户日常输入参数。
   - fragment 文件根据 `GSE`、`GSM` 和命名规则自动发现。
   - barcode 文件根据 `GSE`、`GSM` 和命名规则自动发现。

2. 确定初始 barcode 集合。
   - 如果存在 `*_filtered_barcodes.tsv.gz`，优先直接读取。
   - 如果不存在，但存在 `*singlecell*.csv.gz`，优先使用其中的官方 cell barcode 调用结果作为初始细胞集合。
   - `singlecell.csv.gz` 当前优先尝试使用以下字段判定 cell barcode：
     - `is_cell_barcode` / `is__cell_barcode`
     - `cell_id`
     - `passed_filters`（仅作为更弱的回退）
   - 如果不存在，不预先在 `data/raw` 下生成新的 `filtered_barcodes` 文件。
   - 缺失官方 `filtered_barcodes` 且缺失 `singlecell.csv.gz` 的样本，barcode 过滤在单样本 QC 流程内部完成。
   - 当前预筛方法：先计算 barcode-rank 的 `T_knee`，再从前段 rank 曲线中计算更严格的 `T_inflection`，同时从 `log10(fragment_count + 1)` 密度分布中自动找两峰之间的谷底 `T_floor_auto`。
   - 候选 barcode 当前同时受两类约束：
     - count 阈值：`max(T_knee, T_inflection, T_floor_auto)`
     - rank 约束：只保留前 `N_candidate_auto` 个 barcode，其中 `N_candidate_auto` 当前由 `T_inflection` 对应的 rank 给出
   - 如果密度谷底无法稳定识别，则退回到与 `T_knee` 挂钩的保底阈值。
   - 初始 barcode 集合用于进入矩阵构建，不等同于最终 QC 保留细胞。
   - 当前 `singlecell.csv.gz` 中的附加逐 barcode 指标会带 `singlecell_` 前缀写入 metadata，供后续对照和复查。

3. 读取统一参考。
   - 加载 `data/reference/peak.bed`。
   - 准备 hg38 注释并统一染色体命名。

4. 构建 fragment 对象与 peak×cell 矩阵。
   - 使用初始 barcode 集合作为输入细胞集合。
   - 检查 tabix 索引，缺失时自动创建。
   - 使用 `FeatureMatrix()` 生成统一 peak 坐标下的计数矩阵。

5. 构建 Seurat/Signac 对象。
   - 创建 `ChromatinAssay`。
   - 写入注释、样本编号和数据集编号。

6. 计算 QC 指标。
   - `nCount_ATAC`
   - `nFeature_ATAC`
   - `TSS.enrichment`
   - 当前环境下如果需要生成 `TSSPlot`，应确保 `TSSEnrichment` 保留位置富集矩阵。
   - `nucleosome_signal`
   - `total_fragments`
   - `FRiP`
   - `unique_ratio`
   - `blacklist_fraction`

7. 进行 doublet 检测。
   - 当前使用 `scDblFinder`。
   - 输出 doublet 分类和分数，写入 metadata。

8. 进行 QC 过滤。
   - 当前主流程使用 MAD 方式识别异常值。
   - 当前 MAD 主要基于 `nCount_ATAC`、`TSS.enrichment`、`FRiP`。
   - doublet 检测结果参与最终过滤。
   - 当前不在 notebook 中启用固定阈值二次过滤。
   - 如果后续需要统一阈值，应先比较多样本 QC 分布，再决定是否在整合前增加一轮标准化过滤。

9. 输出单样本结果。
   - 每个样本保存 1 张总 QC 图 `qc_overview.png`，内部包含多个子图
- 当前 `qc_overview` 不再展示前置 barcode 诊断图；即使 barcode 集合来自 fragment 推断或 `singlecell.csv.gz`，总图也只保留正式 QC 图。
   - `matrix/matrix.mtx`
   - `matrix/barcodes.tsv.gz`
   - `matrix/features.tsv.gz`
   - `metadata.csv`
   - `metadata_qc.csv`
   - `qc_summary.csv`
   - `GSM*_seurat_qc.rds`
   - QC 可视化步骤应尽量具备容错性，避免单个子图失败导致整批图输出中断

10. 整合前统一硬阈值筛选。
   - 这一步不覆盖原始单样本 QC 结果。
   - 当前基于已有 `metadata_qc.csv` 与 `matrix/` 再做一次二次筛选。
   - 当前默认保守阈值为：
     - `nCount_ATAC >= 1000`
     - `nCount_ATAC <= 100000`
     - `TSS.enrichment >= 4`
     - `FRiP >= 0.35`
     - `blacklist_fraction <= 0.05`
     - `nucleosome_signal <= 4`
   - 当前不把 `unique_ratio` 作为统一硬阈值，因为不同数据集定义和缺失情况不一致。
   - 输出目录为每个样本下的 `integration_qc/`，其中包含：
     - `matrix/`
     - `metadata_integration_qc.csv`
     - `integration_hard_qc_summary.csv`

单样本阶段的核心目标不是立即做下游分析，而是先回答两个问题：
- 这个样本质量是否足够好？
- 这个样本能否稳定转换成统一格式的矩阵结果？

## 当前脚本化状态
- 单样本脚本 `process_single_sample.R` 已经落地，并通过 `GSE283744/GSM8671454` 真实运行验证。
- 流程管理入口 `pipeline.py` 已经落地。
- 下载脚本 `download_from_datasets.py` 已经落地。
- 当前 `pipeline.py` 支持：
  - `discover`
  - `run-sample`
  - `run-gse`
  - `download`
  - `status`
- 脚本接口保持收窄，日常运行优先只传 `GSE`、`GSM`。
- `peak.bed`、输出根目录、文件自动发现规则在脚本内部按项目约定统一处理。
- 下载阶段当前依据 `datasets.xlsx` 过滤 `scATAC + fragment`。
- 按当前 workbook 内容，这个过滤结果是 28 条 GSM 记录，对应 2 个 GSE：`GSE190992` 和 `GSE283744`。
- 当前机器已安装 `aria2c`，下载入口可直接调用。
- 真实联网解析 GEO supplementary 时可能受站点响应速度影响，当前脚本通过每请求超时和 listing 缓存控制等待时间。
- `datasets.xlsx` 当前只提供 GSE 页面链接，不直接提供 fragment/barcode 文件直链；精确文件 URL 仍由脚本从 GEO supplementary 目录解析。
- 下载验证当前以 `gzip -t` 为准，不以文件大小是否匹配作为最终完整性判断。
- 目前 28 个样本的 fragment 文件都已经完成下载，并补齐了 `.tbi`。
- `GSE190992` 的 18 个样本原始 GEO 中没有 `filtered_barcodes` 文件；此前自动生成的 fallback barcode 文件已删除，后续改为在 `process_single_sample.R` 内部完成 barcode 过滤。
- `GSE283744` 的多个 fragment 原始为普通 `gzip`，当前已自动重压成 `BGZF` 后建好 `.tbi`。
- `GSM8671454` 与 `GSM8671455` 曾发现 fragment 压缩损坏，因此旧 QC 结果不可再视为可信；当前已重新下载 fragment，并重新跑完单样本 QC。
- 当前 `process_single_sample.R` 已支持缺失官方 `filtered_barcodes` 时的内部 barcode 预筛，但 `qc_overview` 只保留正式 QC 图，不再展示 barcode 诊断图。
- `GSM5737281` 的首次无白名单验证虽然跑通，但候选 barcode 仍过多；当前已继续收紧为 `T_knee + T_inflection + T_floor_auto + N_candidate_auto` 联合约束，并已重新验证通过。
- 当前发现 `GSE190992` 的 `FRiP / TSS` 整体仍明显弱于 `GSE283744`，因此后续优先尝试下载该数据集在 GEO supplementary 中的 `singlecell.csv.gz`，并改用其官方 cell barcode 调用结果重新进行单样本 QC。
- `GSE190992` 的 `singlecell.csv.gz` 已下载完成，并确认包含官方 barcode 级 cell-calling 字段（如 `is__cell_barcode`、`passed_filters`、`TSS_fragments`、`on_target_fragments`）。
- `GSM5737281` 已用 `singlecell.csv.gz` 成功验证新流程：
  - `barcode_source = singlecell_csv`
  - `input_cells = 8834`（旧版为 `10002`）
  - `pass_qc = 7499`（旧版为 `6912`）
  - `qc_rate = 84.89%`（旧版为 `69.11%`）
  - `median_TSS_enrichment = 5.39`
  - `median_FRiP = 0.5897`
- 基于该验证结果，后续如需重跑 `GSE190992`，应统一改用 `singlecell.csv.gz` 入口，而不是旧的 fragment 预筛入口。
- 当前 `GSM5737281` 的收紧后结果为：
  - `barcode_candidate_threshold = 4887`
  - `input_cells = 10002`
  - `pass_qc = 6912`
  - `median_TSS_enrichment = 5.41`
  - `median_FRiP = 0.5923`
- 当前已确认这版无官方 `filtered_barcodes` 的预筛逻辑可以继续用于 `GSE190992` 的批量运行。
- 当前如需仅调整 `qc_overview.png` 的内容或布局，优先使用 `regenerate_qc_overview.R` 基于现有 `*_seurat_qc.rds` 重生，不重新运行整套单样本 QC。
- 当前整合前硬阈值筛选已独立成 `apply_integration_hard_qc.py`，默认以保守统一阈值对子样本现有 `pass_qc` 结果再筛一轮，不覆盖原始单样本输出。
- 当前 28 个样本的单样本 QC 已全部完成，所有样本目录下都已生成 `qc_summary.csv`、`metadata_qc.csv`、`matrix/`、`qc_overview.png` 与 `GSM*_seurat_qc.rds`。
- 当前 28 个样本的整合前 hard-QC 也已全部完成，所有样本都已生成 `integration_qc/metadata_integration_qc.csv`、`integration_qc/matrix/` 与 `integration_hard_qc_summary.csv`。
- 当前 hard-QC 总体保留 `173,527 / 188,216` 个细胞（`92.20%`）；按数据集分层为：
  - `GSE190992`: `119,342 / 132,871`（`89.82%`）
  - `GSE283744`: `54,185 / 55,345`（`97.90%`）
- 当前 hard-QC 失败主因以 `FRiP_low`、`TSS.enrichment_low` 与 `nCount_ATAC_low` 为主；`blacklist_fraction` 在当前队列中没有成为主要过滤瓶颈。
- 当前最需要关注的 hard-QC 低保留样本包括：`GSM5737291`、`GSM5737286`、`GSM5737292`、`GSM5737298`、`GSM5737296`，且主要由低 `FRiP` 驱动。
- 当前已完成 `integration_merged/` 总矩阵构建，合并结果为 `338,036` peaks、`173,527` cells、`1,038,538,612` non-zero entries。
- 当前已完成 `integration_sketch/` 参数探索输入构建，结果为 `338,036` peaks、`28,000` cells（每样本 `1000`）、`178,242,460` non-zero entries。
- 当前已完成 `integration_sketch_analysis/` sketch 整合结果输出：使用 `30,000` 个 peaks、`28,000` 个 cells、`16` 个 clusters，`mixing_source_gse_mean = 0.5414`，`mixing_source_gsm_mean = 0.9485`。
- 当前已完成 `integration_merged_analysis/` 全量正式整合结果输出：使用 `30,000` 个 peaks、`173,527` 个 cells、`22` 个 clusters，`mixing_source_gse_mean = 0.4202`，`mixing_source_gsm_mean = 0.9266`。
- 当前 sketch / full 两套整合分析目录都已输出：
  - `pre_harmony_lsi_by_gse.png`
  - `pre_harmony_lsi_by_gsm.png`
  - `post_harmony_umap_by_gse.png`
  - `post_harmony_umap_by_gsm.png`
  - `post_harmony_umap_by_qc.png`
  - `post_harmony_umap_by_cluster.png`
  - `cluster_by_sample_heatmap.png`
  - `lsi_qc_correlation_heatmap.png`
  - `batch_mixing_metrics.csv`
  - `cluster_top_accessible_peaks.csv`
  - `integrated_metadata.csv.gz`
  - `integration_summary.json`
  - `integration_report.md`
- 当前已新增 `output/hard_qc_review/`，其中包含：
  - `hard_qc_before_after_metrics.svg`
  - `hard_qc_metric_summary.csv`
  - `hard_qc_failure_reasons.csv`
  - `hard_qc_failure_combinations_top20.csv`
  - `hard_qc_sample_summary.csv`
  - `hard_qc_dataset_summary.csv`
  - `hard_qc_validation.csv`
  - `integration_readiness_report.md`
- 当前 `hard_qc_validation.csv` 已确认 28 个样本重算得到的 hard-QC 结果与现有 `integration_qc` 输出完全一致。
- 当前机器资源概况：
  - `16` 个逻辑 CPU
  - `47.05 GB` 内存（当前可用约 `42.21 GB`）
  - 工作目录所在磁盘当前可用空间约 `110.57 GB`
- 当前批量单样本运行策略仍保持 `pipeline.py` 的串行 `run-gse`；当前瓶颈已不再是单样本 QC，而是如何在不额外制造全矩阵副本的前提下进入整合。
- 当前单样本脚本对 `unique_ratio` 做兼容处理：如果原始 fragment 统计里没有 `reads_count`，则 `unique_ratio` 留空，不再因为该字段缺失导致汇总失败。
- 当前已手动清理 2 个已完整文件遗留的 `.aria2`，并删除 6 个损坏的 fragment 文件及其 `.aria2`，待后续重下。
- 随后又定位到 4 个待补文件，其中 3 个损坏文件及其 `.aria2` 已删除，并生成了仅针对这 4 个文件的重下清单 `data/raw/_download_logs/retry_missing_4.txt`。
- 如后续新增整合脚本或桌面界面，必须继续同步本文件。

## 多样本整合任务
单样本流程稳定后，进入整合阶段。

目标：
- 合并多个样本矩阵
- 保留样本来源和批次信息
- 输出整合后的图，判断整合是否成功
- 评估并选择批次效应校正方案

候选方法：
- Harmony
- Seurat
- scVI 或其他 scverse 方案

前提：
- 单样本 QC 逻辑稳定
- 所有样本输出格式一致
- 样本级 metadata 足够完整

当前整合阶段拆成两步：
1. 先对 28 个 `integration_qc` peak×cell 矩阵做严格 feature 对齐检查，并生成一份总合并矩阵。
2. 再基于合并矩阵做 LSI / Harmony / UMAP 与整合质量验证。

当前整合阶段的实际状态：
- 第 1 步已经完成，`output/integration_merged/` 可直接作为对齐完成后的统一输入与 bookkeeping 产物。
- 第 2 步已经落地为正式脚本 `scripts/process/run_batch_integration.R`，并已对 `integration_sketch/` 与 `integration_merged/` 各跑通一套 LSI / Harmony / UMAP / 聚类与整合质量评估输出。
- 当前流程先用 `output/integration_sketch/` 做参数探索，再将同一整合思路扩展到 `output/integration_merged/` 全量对象。
- 当前需要避免把 `338,036 × 173,527` 的总矩阵及其多个 TF-IDF / SVD 中间副本一次性全部常驻内存。

当前整合验证必须包含：
- Harmony 前后按 `GSE`、`GSM` 着色的低维图
- 按 `nCount_ATAC`、`FRiP`、`TSS.enrichment`、`blacklist_fraction` 着色的图
- `cluster × sample` 组成热图
- LSI 维度与 QC 变量相关性
- 至少一类 batch mixing 定量指标
- marker / gene activity 的可解释性检查

整合阶段至少需要输出以下图之一或多项：
- PCA / LSI 图
- UMAP 图
- 按样本着色的整合图
- 按批次着色的整合图
- 聚类结果图

## QC 原则
- 先使用原始数据附带的 filtered barcodes 作为初筛。
- 如果样本没有 filtered barcodes 文件，必须提供可重复的 fallback 流程生成初始 barcode 集合。
- QC 判断以样本分布和可解释性为主，不在项目说明中提前写死全局硬阈值。
- `blacklist_fraction` 作为常规 QC 指标保留，用于观察异常开放区域信号。
- 所有过滤结果都应记录在 metadata 中，而不是只保留最终细胞列表。
- 如果后续确定统一阈值，应在脚本和 AGENTS 中同步更新，不要只改其中一处。

## 文档同步约定
- 只要修改了项目结构、处理流程、脚本职责或输出规范，就必须同步更新 `AGENTS.md`。
- `AGENTS.md` 记录的是当前实际状态，不应保留已经删除或停用的脚本说明。

## 当前建议的工作顺序
1. 先用 `GSM8671454.ipynb` 确认单样本流程逻辑没问题。
2. 用 `process_single_sample.R` 固化 notebook 逻辑。
3. 用 `pipeline.py` 批量跑多个样本并统一输出格式。
4. 最后进行多样本整合和整合质量评估。
