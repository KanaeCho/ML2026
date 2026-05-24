# QC 指标说明

本文档解释当前项目中常见 RNA/ATAC QC 指标的含义和排查用途。

## RNA QC 指标

| 指标 | 含义 | 常见异常含义 |
| --- | --- | --- |
| `n_counts` | 每个细胞总 UMI/count 数 | 太低可能是低质量细胞；太高可能是 doublet 或高 RNA content |
| `n_genes` | 每个细胞检测到的基因数 | 太低可能是低质量；太高可能是 doublet |
| `pct_mt` | 线粒体基因比例 | 高值可能提示细胞受损或死亡 |
| `pct_ribo` | ribosomal gene 比例 | 高值可能提示技术或细胞状态差异 |
| `doublet_score` | doublet 评分 | 高值提示混合 barcode |
| `is_doublet` | doublet 判断 | True 的细胞不进入 pass-QC |

当前项目 RNA QC 使用 dynamic hybrid MAD，并将阈值写入 `qc_thresholds.json`。

## ATAC QC 指标

| 指标 | 含义 | 常见异常含义 |
| --- | --- | --- |
| `nCount_ATAC` | peak counts 总数 | 太低可能细胞质量差；太高可能 doublet |
| `nFeature_ATAC` | 检测到的 peak 数 | 太低可能信号不足；太高可能 doublet |
| `TSS.enrichment` | TSS 附近信号富集 | 低值常提示 ATAC 文库质量差 |
| `FRiP` | fragments in peaks fraction | 低值常提示有效开放区域信号低 |
| `nucleosome_signal` | nucleosome pattern 指标 | 异常可能提示文库质量问题 |
| `blacklist_fraction` | blacklist 区域比例 | 高值提示噪音或异常 mapping |
| `fragments` | fragment 总数 | 太低说明信息不足；太高可能 doublet |
| `scDblFinder.score` | doublet 分数 | 高值提示 doublet |
| `scDblFinder.class` | singlet/doublet 分类 | doublet 会被过滤 |

## QC 图怎么看

优先看：

- 过滤前后细胞数。
- 每个关键指标的分布。
- 是否有极端离群值。
- 过滤阈值是否明显过严或过松。
- doublet 比例是否异常。
- ATAC 中 TSS 和 FRiP 是否同时偏低。

## QC 问题排查顺序

1. 看 `qc_summary.csv`。
2. 看 `qc_overview.png`。
3. RNA 看 `qc_thresholds.json`。
4. ATAC 看 barcode source 和 fragment 统计。
5. 看 doublet 比例。
6. 看日志中是否有 fallback 或错误。
7. 再决定是否需要单样本重跑。
