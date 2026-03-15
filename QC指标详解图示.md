# scATAC-seq QC指标详解

## 1. TSS Enrichment（TSS富集分数）

### 原理图示
```
基因结构：
    ↓ TSS (转录起始位点)
====|===================>  基因
    |-2kb-|  |+2kb|

好的ATAC数据：
    ████████              ← TSS附近信号很强（开放染色质）
    ▁▁▁▁▁▁▁▁              ← 远离TSS信号弱（背景）

差的ATAC数据：
    ▃▃▃▃▃▃▃▃              ← TSS附近信号不强
    ▃▃▃▃▃▃▃▃              ← 背景信号也高（噪音）
```

### 计算方法
```
TSS enrichment = TSS±100bp区域平均信号 / TSS±2kb区域背景信号
```

### 判断标准
- **优秀**: > 10
- **良好**: 5-10
- **可接受**: 3-5
- **差**: < 3

---

## 2. Nucleosome Signal（核小体信号）

### 原理图示
```
DNA缠绕在核小体上：

好的细胞 - fragment长度分布有周期性：
数量
 |     ╱╲              ╱╲
 |    ╱  ╲            ╱  ╲
 |___╱____╲__________╱____╲___
     147bp  294bp   441bp      fragment长度
     (单核小体)(双核小体)

差的细胞 - fragment长度分布混乱：
数量
 |   ▃▅▃▅▃▅▃▅▃▅▃▅▃▅▃▅
 |___________________________
                fragment长度
```

### 计算方法
```
nucleosome_signal = 147-294bp fragments / <147bp fragments
```

### 判断标准
- **优秀**: < 1
- **良好**: 1-2
- **可接受**: 2-4
- **差**: > 4

---

## 3. FRiP (Fraction of Reads in Peaks)

### 原理图示
```
基因组：
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Peaks（开放染色质区域）：
    ████      ████    ████

好的细胞 - 大部分reads落在peaks中：
    ████      ████    ████     ← 80%的reads
    ▃▃        ▃▃      ▃▃       ← 20%的reads

差的细胞 - reads分散在全基因组：
    ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
```

### 计算方法
```
FRiP = peaks中的reads数 / 总reads数
```

### 判断标准
- **优秀**: > 0.6 (60%)
- **良好**: 0.4-0.6
- **可接受**: 0.2-0.4
- **差**: < 0.2

---

## 4. nCount_ATAC vs fragments

### 区别说明
```
Fragment文件中的所有fragments：
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
    ↑                                  ↑
    fragments = 50,000 (总数)

只统计peaks中的fragments：
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    ████      ████    ████
    ↑                      ↑
    nCount_ATAC = 20,000 (peaks中)

FRiP = 20,000 / 50,000 = 0.4
```

---

## 5. Doublet（双细胞）

### 原理图示
```
正常细胞（Singlet）：
    ┌─────────┐
    │  Cell A │  → 正常的peak模式
    └─────────┘

Doublet（两个细胞混在一起）：
    ┌─────────┐
    │  Cell A │
    │    +    │  → 异常高的信号
    │  Cell B │     peaks数量多
    └─────────┘
```

### scDblFinder检测方法
1. 模拟人工doublets（随机组合两个细胞）
2. 训练分类器区分真实细胞和人工doublets
3. 给每个细胞打分（0-1）

### 判断标准
- **scDblFinder.class = "singlet"**: 保留
- **scDblFinder.class = "doublet"**: 过滤
- **scDblFinder.score < 0.3**: 可信的singlet
- **scDblFinder.score > 0.5**: 可能是doublet

---

## 6. Blacklist Fraction

### 原理图示
```
基因组：
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Blacklist区域（需要排除）：
    ████              ████
    ↑                 ↑
    重复序列          高信号假阳性区域

好的细胞：
    ▃▃                ▃▃     ← 很少reads落在blacklist
    ████████████████████     ← 大部分reads在正常区域

差的细胞：
    ████              ████   ← 很多reads落在blacklist
    ▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃▃
```

### 判断标准
- **优秀**: < 0.01 (1%)
- **良好**: 0.01-0.05
- **可接受**: 0.05-0.1
- **差**: > 0.1

---

## 完整质控流程图

```
原始数据
   ↓
[1] 创建Fragment对象
   ├─ fragment文件（.tsv.gz）
   ├─ peaks文件（.bed）
   └─ barcodes文件（.tsv）
   ↓
[2] 生成peak-by-cell矩阵
   ↓
[3] 创建Seurat对象
   ├─ 添加基因组注释
   └─ 添加blacklist区域
   ↓
[4] 计算QC指标
   ├─ nCount_ATAC (自动)
   ├─ nFeature_ATAC (自动)
   ├─ NucleosomeSignal()
   ├─ TSSEnrichment()
   ├─ FRiP()
   ├─ blacklist_fraction
   └─ scDblFinder()
   ↓
[5] 质控过滤
   ├─ nCount_ATAC: 2,000-100,000
   ├─ TSS.enrichment > 5
   ├─ nucleosome_signal < 2
   ├─ FRiP > 0.4
   ├─ blacklist_fraction < 0.05
   └─ scDblFinder.class == "singlet"
   ↓
高质量细胞
```

---

## 师姐 vs 你的流程对比

### 你的流程（简化版）
```
Fragment文件
   ↓
统计每个barcode的fragment数
   ↓
过滤：fragments > 1000
   ↓
生成矩阵（只用过滤后的细胞）
   ↓
创建Seurat对象
   ↓
整合、聚类
```
**问题：**
- 只用一个指标（fragment数）
- 标准太宽松（1000太低）
- 没有评估数据质量（TSS、nucleosome）
- 没有doublet检测

### 师姐的流程（标准版）
```
Fragment文件
   ↓
生成矩阵（包含所有细胞）
   ↓
创建Seurat对象
   ↓
计算6个QC指标：
  - nCount_ATAC
  - TSS enrichment
  - nucleosome signal
  - FRiP
  - blacklist fraction
  - doublet score
   ↓
多维度过滤（6个标准）
   ↓
高质量细胞
   ↓
整合、聚类
```
**优势：**
- 多维度评估
- 标准严格
- 符合领域规范

---

## 实际案例对比

假设你的GSM8671454样本：

### 你的流程结果：
- 原始barcodes: 426,131
- 过滤后细胞: 6,988
- 过滤标准: fragments > 1000
- **问题**: 可能包含很多低质量细胞

### 师姐的流程预期结果：
- 原始barcodes: 426,131
- fragments > 2000: ~10,000
- + TSS > 5: ~8,000
- + nucleosome < 2: ~7,000
- + FRiP > 0.4: ~6,000
- + 非doublet: ~5,500
- **最终**: ~5,500 高质量细胞

**差异分析：**
- 你的6,988细胞中，可能有1,500个低质量细胞
- 这些低质量细胞会：
  - 形成低质量clusters
  - 影响整合效果
  - 降低下游分析准确性
