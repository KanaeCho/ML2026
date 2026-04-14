# only_rna Python 主线设计

## 1. 目标

本设计定义 `only_rna` 分支的新主线：在 `scripts/only_rna/` 下建立一个**全新、Python-first、单样本 scRNA-seq 处理子系统**，用于处理 `datasets.xlsx` 中已筛选、且原始输入已完整下载到本地的数据。

第一阶段目标是让每个受支持样本完成以下流程：

- 单样本输入发现与读取
- QC 指标计算
- Python doublet 检测
- `pass_qc` 判定
- 基于 `pass_qc` 细胞的降维、聚类、UMAP
- CIMA RNA `L1/L2` 注释
- 可复核的审计输出与更易读的 UMAP 图

新主线的目标是**功能等价**，而不是复现当前 `Seurat/scDblFinder` 产物的逐值一致性。

## 2. 非目标

本设计明确不覆盖以下内容：

- 跨样本 RNA integration
- ATAC 或 TEA-seq 主线改造
- 对旧 `scripts/process/` 代码做增量重构后复用其核心处理逻辑
- 第一阶段内追求与旧 R 主线完全一致的聚类编号、doublet 结果或 UMAP 坐标
- 把共享 `gene_count` / `gene_counts` CSV 视为受支持单细胞输入

旧流程仍可暂时保留在仓库中作为历史实现，但新主线的设计与实现边界应独立清晰。

## 3. 总体架构决策

### 3.1 架构方向

采用 **“兼容外壳 + 全新 Python 主线”**：

- 尽量保留现有 CLI 命令族与数据发现契约，降低分支内使用方式的迁移成本。
- 单样本处理、注释、绘图、输出写出全部在 `scripts/only_rna/` 中用新 Python 代码实现。
- 新主线不复用旧 `scripts/process/process_single_rna_sample.R` 的分析主体。

### 3.2 技术路线

处理栈采用：

- `AnnData` 作为样本级对象
- `Scanpy` 作为基础单细胞处理主线
- Python doublet 工具首版即接入，当前首选为 `scanpy.pp.scrublet()`
- `matplotlib` / `scanpy.pl` 为 UMAP 输出基础

### 3.3 结果等价定义

“功能等价”定义为：

- 同类输入能被发现并运行
- 能输出 QC、聚类、UMAP、CIMA L1/L2、审计文件
- 输出的主文件集与当前主线保持同一工作目的
- 图可读性与审计信息达到当前分支验收目标

不要求：

- 聚类标签数值与旧流程相同
- UMAP 坐标相同
- doublet 判定逐细胞一致
- CIMA score 数值完全一致

## 4. 运行环境与依赖约束

### 4.1 数据根解析

新主线必须继续遵守当前仓库约定的数据根解析顺序：

1. `./data`
2. 环境变量 `ML2026_DATA_ROOT`
3. `/mnt/g/ML2026_data`

这项规则属于外部运行契约，不应在新主线中改变。

### 4.2 输出根

默认输出根继续对齐当前主线语义，即使用 `output/rna/` 作为默认输出入口。当前工作区中的 `output/` 为软链接这一事实不属于新主线逻辑的一部分，但新主线必须兼容这种工作区布局。

### 4.3 Python 依赖

现有 `pyproject.toml` 中尚未声明 `scanpy`、`anndata` 和 Python doublet 相关依赖，因此新主线实现时需要显式补充这些依赖。第一阶段设计假定以下依赖会被加入项目依赖集中：

- `scanpy`
- `anndata`
- `umap-learn`（如未被上游依赖间接提供，也应显式声明）
- `pyyaml`

若 `scanpy.pp.scrublet()` 的运行还需要额外依赖，则应一并声明，而不是让运行时隐式失败。

## 5. CLI 设计

### 5.1 保持接近现有格式

CLI 命令组保持接近当前格式，继续提供：

- `discover-rna`
- `run-rna-sample`
- `run-rna-gse`
- `rna-status`

这组命令在新主线中可由新的 Python 入口实现，但外部接口风格应尽量维持稳定。

### 5.2 命令职责

#### `discover-rna`

负责：

- 读取 `datasets.xlsx` 中当前可见、且 `assay=scRNA` 的行
- 结合本地原始目录布局发现可运行样本
- 标记受支持与不受支持输入
- 输出供用户审查的发现结果

#### `run-rna-sample`

负责：

- 运行一个显式的样本级对象
- 仅支持真正的 sample-level 目标

明确**禁止**把 GSE 级共享矩阵样本作为 `run-rna-sample` 的显式目标。

#### `run-rna-gse`

负责：

- 运行一个 GSE 下的新主线样本集合
- 在内部隐式处理 GSE 级共享 triplet 的情况

也就是说，**GSE 级共享矩阵是必须支持的，但主要通过 `run-rna-gse` 暴露，不作为单独的显式 sample 命令入口。**

#### `rna-status`

负责：

- 汇总发现状态
- 汇总执行状态
- 检查每个样本输出是否完整

它应继续承担“是否跑完、产物是否齐全”的轻量审计职责。

## 6. 输入发现契约

### 6.1 数据集选择来源

发现逻辑必须继续以 `reference/datasets.xlsx` 第一张表中：

- 当前可见行
- `assay = scRNA`

作为候选数据集来源。

### 6.2 支持的本地原始输入格式

新主线必须沿用当前支持范围，识别以下输入：

1. **每个 GSM 一套 triplet**
   - `matrix.mtx(.gz)`
   - `barcodes.tsv(.gz)`
   - `features.tsv(.gz)` 或 `genes.tsv(.gz)`

2. **每个 GSM 一个 10x `.h5`**

3. **每个 GSM 一个 `matrix.tar.gz`**

4. **GSE 级共享 Matrix Market triplet**
   - 以 GSE 为命名主体
   - 发现后作为一个共享矩阵运行对象

### 6.3 GSE 级共享矩阵支持策略

对 GSE 级共享 triplet：

- 它是新主线中的**一级支持输入类型**
- 发现阶段必须显式识别
- 运行阶段必须可执行
- 但外部 CLI 上主要通过 `run-rna-gse` 隐式纳入

### 6.4 不支持输入

对于共享 `gene_count` / `gene_counts` CSV：

- 发现阶段必须显式标记为 unsupported
- 不能静默当作单细胞输入运行
- 状态输出中应保留“不支持”的原因

### 6.5 `GSE226039` 特例

对 `GSE226039`：

- 仅分析文件名中可判定为 `PBMC` 的样本
- Ileum / Rectum 等其他组织文件必须在新主线 discovery 中继续被排除

这是已确认的业务规则，不应在第一阶段改变。

## 7. 样本处理流程设计

### 7.1 总体顺序

单样本处理流程固定为：

1. 读取输入矩阵为 `AnnData`
2. 计算基础 QC 指标
3. 执行 Python doublet 检测
4. 计算 `pass_qc`
5. 仅对 `pass_qc` 细胞执行：
   - 归一化
   - 特征选择
   - PCA
   - 邻居图
   - 聚类
   - UMAP
   - CIMA `L1/L2` 注释
6. 将结果回写到样本对象与导出文件

### 7.2 输入读入语义

不同输入格式最终都要被统一转换成具有以下最小语义的 `AnnData`：

- `X` 为细胞 × 基因计数矩阵
- `obs_names` 为唯一细胞条码
- `var_names` 为唯一基因名或规范化后基因字段
- `obs` 中写入样本级别元信息，如 `gse`、`sample_id`、`input_type`

### 7.3 基础 QC 指标

至少应计算并保留：

- `n_counts`
- `n_genes`
- `pct_mt`
- `pct_ribo`

若后续实现中需要兼容旧命名，可在导出时提供映射，但内部设计应以清晰的 Python 风格字段为主。

### 7.4 Doublet 设计

第一阶段默认方案：

- 使用 Python 侧 doublet 工具
- 当前首选：`scanpy.pp.scrublet()`

输出至少应写回：

- `doublet_score`
- `is_doublet`

在设计语义上，doublet 结果属于 `pass_qc` 的**硬过滤条件之一**。

### 7.5 `pass_qc` 语义

`pass_qc` 为一个总布尔字段，由下列条件组合得到：

- 计数阈值通过
- feature 阈值通过
- `pct_mt` 阈值通过
- `pct_ribo` 阈值通过
- doublet 过滤通过

同时应保留各分项失败原因字段，便于审计，例如：

- `fails_count_floor`
- `fails_gene_floor`
- `fails_mt_ceiling`
- `fails_ribo_ceiling`
- `fails_doublet`

具体字段名可在实现阶段进一步标准化，但必须满足“全细胞可审计”的设计意图。

### 7.6 配置化阈值

QC 阈值不应硬编码在处理逻辑中，而应采用：

- **YAML 默认配置**
- **CLI 局部覆盖**

配置内容至少需要覆盖：

- count 下限
- gene 下限
- `pct_mt` 上限
- `pct_ribo` 上限
- doublet 相关参数
- 降维/聚类关键参数
- 绘图导出参数

CLI 覆盖仅用于少量运行级调整，不应替代 YAML 作为主配置源。

### 7.7 降维与聚类

对 `pass_qc` 细胞执行新主线嵌入流程。第一阶段不要求完全复刻 Seurat 行为，但需要形成稳定、可复核的 query-native RNA 结构。设计上应至少包括：

- 标准化
- 高变特征或等价特征选择
- PCA
- 邻居图
- 聚类
- UMAP

聚类结果与 UMAP 结果必须回写到 `pass_qc` 子集对应的全细胞元数据中，其余被过滤细胞相关列可为空。

### 7.8 CIMA 注释

第一阶段只计算：

- `L1`
- `L2`

且只对 `pass_qc` 细胞执行注释。新主线必须继续消费当前外部 CIMA reference 资产，包括：

- `cima_rna_reference_pca_features.tsv.gz`
- `cima_rna_reference_l1_centroids.tsv`
- `cima_rna_reference_l2_centroids.tsv`
- `cima_rna_celltype_hierarchy.csv`
- `cima_rna_reference_model.json`

注释输出至少应包含：

- `cima_l1`
- `cima_l2`
- `cima_l1_score`
- `cima_l1_score_margin`
- `cima_l2_score`
- `cima_l2_score_margin`

并保留低置信逻辑及 masked L1 视图所需字段，例如：

- `cima_l1_low_confidence`
- `cima_l1_masked`

## 8. 输出设计

### 8.1 输出目录

单样本默认输出目录继续为：

- `output/rna/{GSE}/{sample_id}/`

### 8.2 产物全集

新主线第一阶段保留当前主线的**产物全集**，包括：

- `qc_overview.png`
- `metadata.csv`
- `metadata_qc.csv`
- `qc_summary.csv`
- `validation_result.csv`
- `umap_rna_clusters.png`
- `umap_rna_cima_cell_type_l1.png`
- `umap_rna_cima_cell_type_l2.png`
- `umap_rna_cima_cell_type_l1_masked.png`
- `matrix/matrix.mtx`
- `matrix/barcodes.tsv.gz`
- `matrix/features.tsv.gz`
- 样本级主对象文件 `.h5ad`

其中，旧 `*_seurat_qc.rds` 的角色由 `.h5ad` 取代。

### 8.3 字段命名风格

导出字段允许比旧流程更规范，采用**清晰、Python 风格、`snake_case` 优先**的命名方式。设计原则是：

- 先保证语义清晰
- 再考虑与旧字段的可映射性

不要求为了兼容历史而保留全部旧字段名。

### 8.4 `metadata.csv`

`metadata.csv` 定义为：

- **全细胞总表**

所有输入细胞都必须保留在该表中。对被过滤细胞：

- QC、doublet、`pass_qc` 等列必须有值
- 聚类、UMAP、CIMA 相关列可以为空

### 8.5 `metadata_qc.csv`

`metadata_qc.csv` 定义为：

- **`pass_qc` 子集视图**

该文件应聚焦后续嵌入与注释真正参与的细胞。

### 8.6 `qc_summary.csv`

`qc_summary.csv` 定义为：

- **样本级 QC / 统计摘要**

应汇总至少以下信息：

- 输入总细胞数
- `pass_qc` 细胞数与比例
- doublet 数与比例
- 各类 QC 失败计数
- 低置信 L1 比例
- cluster purity 类统计

### 8.7 `validation_result.csv`

`validation_result.csv` 定义为：

- **流程完成性与关键检查结果**

它不再承担“全细胞主表”角色，而应聚焦验证和审计，例如：

- 关键输出是否存在
- 关键列是否完整
- 注释与 UMAP 是否成功生成
- 样本级统计检查是否通过

如实现上需要混合少量样本级指标与布尔检查也可以，但职责必须与 `qc_summary.csv` 区分清楚。

## 9. 绘图与可读性设计

### 9.1 问题定义

当前旧流程中的 RNA UMAP 图 legend 过小，主要因为图例仅做了轻量文本字号设置，没有建立统一绘图规范。新主线必须从设计上解决这一问题。

### 9.2 绘图规范

UMAP 绘图模块必须提供统一的可读性配置，至少覆盖：

- 图尺寸
- DPI
- 点大小与透明度
- 标题字号
- legend 文字字号
- legend title 字号
- legend key 尺寸
- legend 间距
- 多列 legend 或合适的布局策略

设计目标不是“简单放大一点”，而是形成一套可持续复用的 RNA UMAP 导出规范。

### 9.3 必须产出的图

第一阶段必须导出：

- cluster UMAP
- CIMA L1 UMAP
- CIMA L2 UMAP
- CIMA L1 masked UMAP

### 9.4 低置信显示

低置信 L1 细胞必须能在 masked 图中被统一显示为 `Unknown` 或等价低置信类别，且图例配色与标注方式应可直接支持人工审阅。

## 10. 模块边界与包结构

`scripts/only_rna/` 采用**包结构 + 子模块**，并按流程细拆，而不是少量大文件。

建议结构如下：

- `scripts/only_rna/__init__.py`
- `scripts/only_rna/cli.py`
- `scripts/only_rna/discovery.py`
- `scripts/only_rna/read_inputs.py`
- `scripts/only_rna/qc.py`
- `scripts/only_rna/doublet.py`
- `scripts/only_rna/embedding.py`
- `scripts/only_rna/annotation.py`
- `scripts/only_rna/plotting.py`
- `scripts/only_rna/outputs.py`
- `scripts/only_rna/config.py`
- `scripts/only_rna/models.py`

### 10.1 `cli.py`

负责命令定义、参数解析与调度，不承载分析实现细节。

### 10.2 `discovery.py`

负责数据根解析、`datasets.xlsx` 筛选、输入布局识别、PBMC 特例、unsupported 原因写出。

### 10.3 `read_inputs.py`

负责把不同原始格式读成统一 `AnnData`，并处理必要的格式归一化。

### 10.4 `qc.py`

负责基础 QC 指标与 `pass_qc` 逻辑。

### 10.5 `doublet.py`

负责 doublet 检测逻辑及其结果标准化。

### 10.6 `embedding.py`

负责 `pass_qc` 子集上的嵌入流程。

### 10.7 `annotation.py`

负责 CIMA L1/L2 注释与相关置信字段生成。

### 10.8 `plotting.py`

负责统一的 UMAP 导出规范。

### 10.9 `outputs.py`

负责 CSV、矩阵导出、对象写出和验证文件生成。

### 10.10 `config.py`

负责 YAML 默认配置读取、CLI 覆盖合并、模块级配置分发。

### 10.11 `models.py`

负责轻量 typed data structures，例如：

- discovery 结果
- sample run 描述
- 配置对象
- validation 结果对象

## 11. 测试结构设计

测试纳入本次设计，采用**实现模块细拆、测试按层次合并**的方式。

建议测试目录：

- `tests/only_rna/test_discovery.py`
- `tests/only_rna/test_processing.py`
- `tests/only_rna/test_outputs.py`
- `tests/only_rna/test_cli.py`

### 11.1 `test_discovery.py`

覆盖：

- 数据根解析
- `datasets.xlsx` 可见 `scRNA` 行筛选
- per-GSM triplet / h5 / archive 识别
- shared GSE triplet 识别
- shared gene-count CSV unsupported 判定
- `GSE226039` PBMC-only 过滤

### 11.2 `test_processing.py`

覆盖：

- 读入后的 `AnnData` 基础结构
- QC 指标与 `pass_qc` 逻辑
- doublet 结果接入
- `pass_qc` 子集嵌入
- CIMA 注释只作用于 `pass_qc`

### 11.3 `test_outputs.py`

覆盖：

- `.h5ad` 写出
- `metadata.csv` 全细胞语义
- `metadata_qc.csv` 子集语义
- `qc_summary.csv` 与 `validation_result.csv` 的职责边界
- UMAP 文件写出
- 矩阵导出存在性

### 11.4 `test_cli.py`

覆盖：

- `discover-rna`
- `run-rna-sample`
- `run-rna-gse`
- `rna-status`
- `run-rna-sample` 不允许显式运行 GSE 级共享样本
- `run-rna-gse` 可隐式纳入共享矩阵

## 12. 状态与验证设计

新主线应保留“状态文件 + 输出完整性检查”的思路，使 `rna-status` 仍能回答两个问题：

1. 这个样本是否已经执行完成？
2. 这个样本的关键产物是否齐全？

输出完整性检查至少应覆盖产物全集中的关键文件。

## 13. 与旧主线的关系

第一阶段实现过程中，旧 `scripts/process/` 路径可以继续保留，但新主线不应继续扩大对其内部分析逻辑的依赖。迁移目标应是让 `scripts/only_rna/` 逐步成为该分支默认的可维护 Python 子系统。

这意味着：

- 旧流程可作为行为参考
- 旧 discovery / CLI 契约可作为兼容外壳参考
- 但新主线分析实现应在新目录中自洽闭环

## 14. 风险与阶段边界

### 14.1 第一阶段主要风险

- Python doublet 结果与旧流程差异较大
- Scanpy 主线下的聚类/UMAP 形态与当前 R 实现不同
- 不同原始输入格式归一化时出现字段不一致
- CIMA 参考资产接入过程中出现基因对齐问题

### 14.2 应对原则

- 以“功能等价 + 可审计”优先，而不是追求旧流程数值复刻
- 在测试层把 discovery 契约、输出职责和 CLI 行为锁死
- 把配置与处理逻辑解耦，降低后续调参成本

## 15. 实施后文档同步要求

只要实现阶段对以下任一项产生实际变更，就必须同步更新 `AGENTS.md`：

- RNA 主线目标
- discovery 规则
- QC 逻辑
- 注释层级
- 输出规范
- 命令接口

`AGENTS.md` 只能记录实现后已确认的真实行为，不能把本设计中的未来意图提前写成既成事实。

## 16. 设计结论

本设计最终确定：

- 新主线是 `scripts/only_rna/` 下的全新 Python-first 子系统
- CLI 保持接近现有命令族
- discovery 规则完全沿用当前契约，并明确支持 GSE 级共享 triplet
- `run-rna-gse` 隐式处理 GSE 级共享矩阵，`run-rna-sample` 禁止显式运行该类目标
- 处理对象为 `AnnData`
- 处理主线为 `Scanpy`
- 首版接入 Python doublet，并作为 `pass_qc` 硬过滤条件
- CIMA 仅对 `pass_qc` 细胞做 `L1/L2` 注释
- 输出保留当前产物全集，但主对象改为 `.h5ad`
- 字段命名允许更规范的 Python 风格
- UMAP 图例与版式采用统一可读性规范
- 代码模块细拆，测试按层次合并

该设计已经足以支撑下一步实现计划编写。
