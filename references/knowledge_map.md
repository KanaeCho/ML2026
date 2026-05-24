# RNA/ATAC 知识框架

本文档用于整理当前项目涉及的 RNA/ATAC 分析知识，方便后续复习、写方法、做 PPT 和排查问题。

## RNA 知识整理

### RNA 分析主要步骤

1. 准备原始 count matrix。
2. 读取样本 metadata。
3. 计算 QC 指标。
4. 根据 QC 指标和 doublet 结果过滤细胞。
5. 归一化和 log transform。
6. 选择高变基因。
7. PCA 降维。
8. 构建 neighbors graph。
9. Leiden 聚类。
10. UMAP 投影。
11. cell type 注释。
12. 导出 metadata、QC summary、h5ad 和图。
13. 多样本 product-level 低维整合。

### 每一步目的

| 步骤 | 目的 |
| --- | --- |
| 输入读取 | 将不同格式的 count matrix 统一成 AnnData |
| QC | 识别低质量细胞、异常细胞和技术噪音 |
| doublet | 识别可能由多个细胞混合形成的 barcode |
| 过滤 | 保留可靠细胞进入下游分析 |
| 归一化 | 减少 library size 差异影响 |
| HVG | 选择最能反映细胞差异的基因 |
| PCA | 降低维度、减少噪音 |
| 聚类 | 将相似细胞分组 |
| UMAP | 可视化细胞关系 |
| 注释 | 将 cluster/cell 映射到生物学细胞类型 |
| 导出 | 保存可复核的结果和审计文件 |

### 当前项目对应代码

| 知识点 | 代码位置 |
| --- | --- |
| 样本发现 | `scripts/only_rna/discovery.py` |
| 输入读取 | `scripts/only_rna/read_inputs.py` |
| QC | `scripts/only_rna/qc.py` |
| doublet | `scripts/only_rna/doublet.py` |
| embedding | `scripts/only_rna/embedding.py` |
| annotation | `scripts/only_rna/annotation.py`, `scripts/only_rna/azimuth.py` |
| final celltype | `scripts/only_rna/final_celltype.py` |
| output | `scripts/only_rna/outputs.py` |
| product integration | `scripts/process/integrate_product_embeddings.py` |

### 当前项目关键参数

详见 `docs/parameters.md`。核心参数包括：

- RNA QC nMAD 参数。
- count/gene/mt/ribo guardrail。
- normalization target sum 10000。
- HVG flavor `seurat`。
- PCA 默认 30 components。
- Azimuth reference `pbmcref`。
- final celltype 5 类。

### 当前项目已有结果图

- `qc_overview.png`
- `umap_rna_pbmcref_vs_cima_l1.png`
- `umap_rna_pbmcref_highlight.png`
- `umap_rna_cima_l1.png`
- product-level `*_cima_l1_panels.png`
- product-level `*_cima_l2_panels.png`
- product-level `*_integrated_cluster_panels.png`

## ATAC 知识整理

### ATAC 分析主要步骤

1. 准备 fragment 文件和 barcode 文件。
2. 读取 peak reference。
3. 构建 fragment object。
4. 构建 peak-by-cell matrix。
5. 创建 ChromatinAssay 和 Seurat object。
6. 计算 TSS enrichment、FRiP、nucleosome signal 等 QC 指标。
7. doublet 检测。
8. QC 过滤。
9. LSI 降维。
10. 聚类和 UMAP。
11. CIMA ATAC 注释。
12. 导出 QC、metadata、UMAP、matrix、h5ad。
13. product-level 低维整合。

### ATAC QC 指标含义

| 指标 | 含义 |
| --- | --- |
| `nCount_ATAC` | peak 区域总 counts |
| `nFeature_ATAC` | 检测到的 peak 数 |
| `TSS.enrichment` | TSS 附近信号富集程度，反映文库质量 |
| `FRiP` | fragments in peaks fraction，反映 fragments 落入 peak 的比例 |
| `nucleosome_signal` | nucleosome pattern 相关质量指标 |
| `blacklist_fraction` | 落入 blacklist 区域的比例 |
| `fragments` | fragment 总数 |
| `scDblFinder.class` | singlet/doublet 判断 |

### fragments、peak matrix、gene activity 的关系

- fragments 是 ATAC 原始或近原始输入，记录每条 fragment 的基因组位置和 barcode。
- peak matrix 是将 fragments 计数到 peak 区域后得到的 peak-by-cell 矩阵。
- gene activity 通常是将 ATAC 信号汇总到基因附近得到的 gene-level proxy。
- 当前项目中已确认 fragments 和 peak matrix；gene activity 主线输出当前项目中未找到。

### LSI、聚类、投影关系

- ATAC 数据稀疏，通常先做 TF-IDF/LSI 类降维。
- 聚类基于低维表示或 graph。
- UMAP 用于展示低维空间中的细胞结构。
- 当前项目还会将样本投影到 CIMA ATAC reference-space，用于 CIMA 注释和可视化。

### 当前项目对应代码

| 知识点 | 代码位置 |
| --- | --- |
| ATAC 主线 | `scripts/process/process_single_sample.R` |
| only ATAC 路由 | `scripts/process/pipeline.py` |
| co ATAC 路由 | `scripts/co/cli.py` |
| longevity ATAC 路由 | `scripts/longevity/cli.py` |
| h5ad export | `scripts/process/export_co_atac_h5ad.py` |
| product integration | `scripts/process/integrate_product_embeddings.py` |

### 当前项目已有结果图

- `qc_overview.png`
- `umap_cima_cell_type_l1.png`
- `umap_cima_cell_type_l2.png`
- product-level `*_cima_l1_panels.png`
- product-level `*_cima_l2_panels.png`
- product-level `*_integrated_cluster_panels.png`

## 多样本批处理知识整理

### 为什么要批处理

多样本项目需要保证每个样本使用一致的流程、参数、输出格式和状态记录。批处理可以减少手工操作错误，并让失败样本可定位、可重跑。

### sample sheet 或 metadata 如何组织

当前项目使用不同分支的样本表：

- only RNA：`datasets.xlsx`
- only ATAC：`atac.xlsx`
- co：`co.xlsx`
- longevity：扫描 raw 目录

通用模板见 `configs/sample_sheet_template.csv`。

### 批处理结果如何汇总

样本级结果通过 `run_status.json`、`qc_summary.csv` 和 `validation_result.csv` 检查。product-level 通过 manifests、QC summary、integration metrics 和 panel figures 汇总。

## 异常问题知识整理

### 常见异常问题类型

- QC 异常。
- 过滤后细胞数过少。
- batch effect 明显。
- UMAP 混乱。
- 聚类不清楚。
- marker 或注释不一致。
- RNA/ATAC 注释不一致。
- fragment/peak matrix 异常。
- 批处理失败。
- 输出缺失。
- 参数和版本不可追溯。

### 哪些问题只需要改文档、label 或图

- 图标题不清楚。
- label 字段说明不清楚。
- 当前结果正确但路径或版本没有记录。
- 需要解释 low confidence label 的含义。

### 哪些问题可能需要小范围重跑

- 单样本 barcode 参数不合适。
- 单样本 QC 参数需要调整。
- UMAP 绘图字段错误，需要重新导出图。
- annotation backend 失败后需要补跑单样本。

### 哪些问题必须人工判断

- marker 和预期不一致。
- RNA/ATAC 注释不一致。
- batch effect 与真实生物差异的区分。
- 是否删除旧结果。
- 是否调整核心阈值或 dataset rule。
