# 师姐的 scATAC-seq 单样本质控流程学习笔记

## 一、整体流程概览

师姐的质控流程分为以下步骤：
1. 数据准备（fragment文件索引）
2. 创建Seurat对象
3. 计算QC指标
4. Doublet检测
5. 保存结果

---

## 二、详细代码解析

### 步骤1: 准备工作

```r
# 加载必需的包
library(Signac)        # scATAC-seq分析核心包
library(Seurat)        # 单细胞数据分析
library(EnsDb.Hsapiens.v86)  # 人类基因组注释
library(scDblFinder)   # Doublet检测

# 获取基因组注释
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"  # 设置染色体命名风格为chr1, chr2...
genome(annotations) <- "hg38"

# 读取blacklist区域（需要过滤的基因组区域）
blacklist_hg38 <- rtracklayer::import("hg38-blacklist.v2.bed")
```

**关键点：**
- `annotations` 包含基因位置信息，用于计算TSS enrichment
- `blacklist` 是ENCODE项目定义的需要排除的基因组区域（重复序列、高信号区域）

---

### 步骤2: 为fragment文件创建索引

```r
# 给每个fragment文件生成.tbi索引（加速读取）
fragment_files <- list.files(fragment_dir, pattern = "\\.fragments.tsv.gz$", full.names = TRUE)
for (f in fragment_files) {
  system(paste("tabix -p bed", shQuote(f)))
}
```

**为什么需要索引？**
- fragment文件通常很大（几GB）
- tabix索引可以快速定位特定区域的fragments
- Signac的很多函数需要indexed fragment文件

---

### 步骤3: 创建Seurat对象（核心函数）

```r
Create_scATAC_object <- function(barcodes_file, fragment_file, peaks) {
    # 1. 读取barcode列表（真实细胞）
    cell <- readLines(barcodes_file)

    # 2. 创建Fragment对象（关联fragment文件）
    frags <- CreateFragmentObject(path = fragment_file, cells = cell)

    # 3. 生成peak-by-cell矩阵
    counts <- FeatureMatrix(
        fragments = frags,
        features = peaks,    # peaks区域
        cells = cell         # 只统计这些细胞
    )

    # 4. 创建ChromatinAssay（包含counts和fragment信息）
    atac.assay <- CreateChromatinAssay(counts, fragments = frags)

    # 5. 创建Seurat对象
    atac.obj <- CreateSeuratObject(atac.assay, assay = "ATAC")

    return(atac.obj)
}
```

**关键理解：**
- `CreateFragmentObject`: 不读取全部数据，只是建立连接
- `FeatureMatrix`: 统计每个peak在每个细胞中的fragment数
- `CreateChromatinAssay`: 同时保存counts矩阵和fragment文件路径

---

### 步骤4: 计算QC指标（最重要！）

```r
# 添加基因组注释（必须在计算TSS之前）
Annotation(atac_obj) <- annotations

# 1. 核小体信号 (Nucleosome Signal)
atac_obj <- NucleosomeSignal(atac_obj)
```
**核小体信号解释：**
- 测量fragment长度分布
- 好的细胞：有明显的核小体周期性（147bp, 294bp...）
- 低质量细胞：fragment长度分布混乱
- **标准：< 2**

```r
# 2. TSS富集分数 (TSS Enrichment)
atac_obj <- TSSEnrichment(atac_obj)
```
**TSS富集分数解释：**
- TSS (Transcription Start Site) = 转录起始位点
- 好的ATAC数据：TSS附近有很强的信号（开放染色质）
- 计算方法：TSS±2kb区域的信号 / 背景信号
- **标准：> 5**

```r
# 3. 统计总fragment数
total_fragments <- CountFragments(fragment_file)
total_fragments <- subset(total_fragments, total_fragments$CB %in% colnames(atac_obj))
row.names(total_fragments) <- total_fragments$CB
atac_obj$fragments <- total_fragments[colnames(atac_obj), "frequency_count"]
```
**为什么需要总fragment数？**
- `nCount_ATAC` = peaks中的fragments
- `fragments` = 所有fragments（包括peaks外的）
- 用于计算FRiP

```r
# 4. FRiP (Fraction of Reads in Peaks)
atac_obj <- FRiP(object = atac_obj, assay = 'ATAC', total.fragments = 'fragments')
```
**FRiP解释：**
- FRiP = peaks中的reads / 总reads
- 好的细胞：大部分reads落在peaks中
- **标准：> 0.4 (40%)**

```r
# 5. Blacklist比例
atac_obj$blacklist_fraction <- FractionCountsInRegion(
    object = atac_obj,
    assay = "ATAC",
    regions = blacklist_hg38
)
```
**Blacklist比例解释：**
- 落在blacklist区域的reads比例
- 这些区域通常是重复序列、高信号假阳性区域
- **标准：< 0.05 (5%)**

```r
# 6. Doublet检测
atac_obj <- runscDblFinder(atac_obj, assay = 'ATAC')
```
**Doublet检测解释：**
- 识别两个细胞混在一起的"假细胞"
- scDblFinder基于模拟doublets来识别
- 输出：`scDblFinder.class` (singlet/doublet) 和 `scDblFinder.score` (0-1)

---

## 三、质控标准总结

师姐使用的过滤标准（从metadata推断）：

| QC指标 | 标准 | 含义 |
|--------|------|------|
| nCount_ATAC | 2,000 - 100,000 | peaks中的fragment数 |
| TSS.enrichment | > 5 | TSS区域富集程度 |
| nucleosome_signal | < 2 | 核小体信号强度 |
| FRiP | > 0.4 | peaks中reads比例 |
| blacklist_fraction | < 0.05 | blacklist区域比例 |
| scDblFinder.class | singlet | 非doublet |

---

## 四、与你之前流程的对比

### 你的流程：
```r
# 1. 统计fragments
fragment_counts <- CountFragments(fragment_file)

# 2. 过滤细胞（只用fragment数）
valid_cells <- fragment_counts$CB[fragment_counts$frequency_count > 1000]

# 3. 生成矩阵
counts <- FeatureMatrix(fragments = frags, features = peaks_gr, cells = valid_cells)
```
**问题：**
- ❌ 只用fragment数过滤，标准太宽松
- ❌ 没有计算TSS、nucleosome、FRiP
- ❌ 没有doublet检测

### 师姐的流程：
```r
# 1. 先创建对象（包含所有细胞）
atac_obj <- Create_scATAC_object(...)

# 2. 计算完整QC指标
atac_obj <- NucleosomeSignal(atac_obj)
atac_obj <- TSSEnrichment(atac_obj)
atac_obj <- FRiP(atac_obj)
atac_obj <- runscDblFinder(atac_obj)

# 3. 用多个标准联合过滤
atac_obj <- subset(atac_obj,
  nCount_ATAC > 2000 & nCount_ATAC < 100000 &
  TSS.enrichment > 5 &
  nucleosome_signal < 2 &
  FRiP > 0.4 &
  blacklist_fraction < 0.05 &
  scDblFinder.class == "singlet"
)
```
**优势：**
- ✓ 多维度评估细胞质量
- ✓ 更严格、更准确
- ✓ 符合scATAC-seq领域标准

---

## 五、关键函数说明

### 1. CreateFragmentObject()
```r
frags <- CreateFragmentObject(path = fragment_file, cells = cell)
```
- 不加载全部数据到内存
- 只建立到fragment文件的连接
- 需要.tbi索引文件

### 2. FeatureMatrix()
```r
counts <- FeatureMatrix(fragments = frags, features = peaks, cells = cell)
```
- 统计每个peak在每个细胞中的fragment数
- 返回稀疏矩阵（peaks × cells）

### 3. NucleosomeSignal()
```r
atac_obj <- NucleosomeSignal(atac_obj)
```
- 需要fragment文件
- 计算每个细胞的核小体周期性
- 添加`nucleosome_signal`列到metadata

### 4. TSSEnrichment()
```r
atac_obj <- TSSEnrichment(atac_obj)
```
- 需要基因组注释（Annotation）
- 计算TSS±2kb区域的富集
- 添加`TSS.enrichment`列到metadata

### 5. FRiP()
```r
atac_obj <- FRiP(object = atac_obj, assay = 'ATAC', total.fragments = 'fragments')
```
- 需要总fragment数（metadata中的`fragments`列）
- 计算peaks中reads占比
- 添加`FRiP`列到metadata

---

## 六、学习建议

1. **先理解概念**：
   - 什么是TSS？为什么TSS区域应该是开放的？
   - 什么是核小体？为什么fragment长度有周期性？
   - 什么是FRiP？为什么高FRiP代表高质量？

2. **实践步骤**：
   - 用你的GSM8671454样本，按师姐的流程重新处理
   - 对比过滤前后的细胞数和质量
   - 可视化QC指标分布

3. **参考资料**：
   - Signac官方教程：https://stuartlab.org/signac/
   - ATAC-seq原理：Buenrostro et al., 2013, Nature Methods
   - scATAC-seq质控标准：10x Genomics官方文档

---

## 七、下一步：应用到你的数据

准备好后，我们可以：
1. 为GSM8671454创建完整的质控脚本
2. 计算所有QC指标
3. 用师姐的标准过滤
4. 对比与之前结果的差异
