# 参数整理

本文档整理当前项目中已经确认的关键参数。参数来源包括 YAML 配置、Python argparse 默认值和 R 脚本参数。未确认是否应统一到 config 的参数会单独标注。

## RNA QC 参数

来源：`scripts/only_rna/default_config.yaml`

| 参数 | 当前值 | 影响步骤 | 是否建议统一到 config | 备注 |
| --- | --- | --- | --- | --- |
| `qc.method` | `dynamic_hybrid_mad` | RNA QC 过滤 | 已在 config | 当前主线方法 |
| `counts_lower_nmads` | 2.5 | `n_counts` lower-tail 阈值 | 已在 config | log10p1 空间 |
| `genes_lower_nmads` | 2.5 | `n_genes` lower-tail 阈值 | 已在 config | log10p1 空间 |
| `pct_mt_upper_nmads` | 2.5 | `pct_mt` upper-tail 阈值 | 已在 config | 原值空间 |
| `pct_ribo_upper_nmads` | 3.5 | `pct_ribo` upper-tail 阈值 | 已在 config | 原值空间 |
| `min_cells_for_dynamic` | 50 | 小样本规则 | 已在 config | 小样本会放宽 nMAD |
| `count_floor_min` | 100 | count guardrail | 已在 config | 最小 counts 下界 |
| `count_floor_max` | 1500 | count guardrail | 已在 config | 最大 counts 下界 |
| `gene_floor_min` | 100 | gene guardrail | 已在 config | 最小 genes 下界 |
| `gene_floor_max` | 1200 | gene guardrail | 已在 config | 最大 genes 下界 |
| `pct_mt_ceiling_min` | 5.0 | mt guardrail | 已在 config | 最小 mt ceiling |
| `pct_mt_ceiling_max` | 40.0 | mt guardrail | 已在 config | 最大 mt ceiling |
| `pct_ribo_ceiling_min` | 20.0 | ribo guardrail | 已在 config | 最小 ribo ceiling |
| `pct_ribo_ceiling_max` | 80.0 | ribo guardrail | 已在 config | 最大 ribo ceiling |

## RNA embedding 参数

来源：`scripts/only_rna/embedding.py`

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| `normalize_total target_sum` | 10000 | RNA 标准化 | 固定写在代码中 |
| `log1p` | enabled | RNA log transform | 固定 |
| `hvg flavor` | `seurat` | 高变基因选择 | 固定 |
| `n_top_genes` | auto | HVG | 大样本最多 1000，小样本最多 2000 |
| `scale max_value` | 10 | PCA 前 scale | 固定 |
| `n_pcs` | 默认 30，受细胞数/基因数限制 | PCA/neighbors | 来自 config model 默认值 |
| `n_neighbors` | auto | neighbors | 大样本最多 10，小样本最多 15 |
| `leiden resolution` | auto | clustering | 大样本 0.5，小样本 1.0 |
| `umap min_dist` | 0.5 | UMAP | 来自 embedding config 默认值 |
| `umap spread` | 1.0 | UMAP | 来自 embedding config 默认值 |
| `random_state` | 42 | UMAP | 来自 embedding config 默认值 |

## RNA annotation 参数

来源：`scripts/only_rna/default_config.yaml`, `scripts/only_rna/azimuth.py`, `scripts/only_rna/outputs.py`

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| `annotation.methods` | `azimuth` | RNA 主注释 | 当前主线默认 |
| `azimuth.reference` | `pbmcref` | Azimuth reference | 当前主线 reference |
| `annotation_levels` | `l1`, `l2` | Azimuth 输出层级 | 配置中记录 |
| final celltype vocabulary | `CD4_T`, `CD8_T`, `B`, `Myeloid`, `NK` | 输出 h5ad 和 metadata | 当前 RNA 输出统一 5 类 |

## RNA tuning 参数

来源：`scripts/only_rna/default_config.yaml`, `scripts/only_rna/tuning_presets.py`

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| `qc_preset_family` | `baseline_only` | tuning candidate | 当前只保留 baseline |
| `azimuth_preset_family` | `baseline_only` | tuning candidate | 当前只保留 baseline |
| `embedding_preset_family` | `baseline_only` | tuning candidate | 当前只保留 baseline |
| `max_candidates` | 1 | tuning candidate 数 | 固定单 candidate |
| candidate ID | `baseline__baseline__baseline` | tuning 输出目录 | 当前唯一 candidate |

## ATAC QC 参数

来源：`scripts/process/process_single_sample.R`, `scripts/process/pipeline.py`

| 参数 | 当前值 | 影响步骤 | 是否建议统一到 config | 备注 |
| --- | --- | --- | --- | --- |
| `--nmads` | 4 | ATAC outlier filtering | 建议 | MAD multiplier |
| `doublet_method` | `scDblFinder` | doublet 检测 | 可记录 | R 脚本中固定使用 |
| `scDblFinder nfeatures` | 25 | doublet 检测 | 可记录 | R 脚本中固定 |
| `scDblFinder processing` | `normFeatures` | doublet 检测 | 可记录 | R 脚本中固定 |

## ATAC barcode 参数

来源：`scripts/process/process_single_sample.R`, `scripts/longevity/cli.py`

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| barcode priority | filtered_barcodes, filtered_metadata, singlecell, fragment inference | barcode 选择 | 当前主线逻辑 |
| `--min-inferred-fragments` | null | fragment-count inference | 可选覆盖 |
| `--max-inferred-barcodes` | null | fragment-count inference | 可选覆盖 |
| longevity `min_fragments` | 200 | ArchR barcode preprocessing | 默认值 |
| longevity `max_barcodes` | 20000 | ArchR barcode preprocessing | 默认值 |
| longevity `min_tss` | 2.5 | ArchR barcode preprocessing | 默认值 |
| longevity `rank_by` | `fragments` | barcode ranking | 可选 `tss_then_fragments` |

## ATAC annotation 和 UMAP 参数

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| `peak_file` | `data_root/reference/peak.bed` | FeatureMatrix | 需要 reference 文件 |
| CIMA levels | L1-L4 | ATAC 注释 | CIMA centroid assignment |
| `--umap-min-dist` | null | query-native UMAP | null 表示使用 Seurat/uwot 默认 |

## product-level integration 参数

来源：`scripts/process/organize_integrated_products.py`, `scripts/process/integrate_product_embeddings.py`

| 参数 | 当前值 | 影响步骤 | 备注 |
| --- | --- | --- | --- |
| products | `only_atac`, `only_rna`, `co_atac`, `co_rna` | product organization | `--products all` 处理全部 |
| copy mode | `symlink` | 样本输出引用 | 失败时可降级 pointer |
| n_components | 30 | CIMA compact embedding | 默认 |
| max_umap_fit_cells | 100000 | UMAP fit | 控制内存 |
| RNA method | `harmony` | RNA product integration | only_rna/co_rna 默认 |
| ATAC method | `bbknn` | ATAC product integration | only_atac/co_atac 默认 |
| only_rna batch_key | `gse` | Harmony batch | 默认 |
| other batch_key | `sample_id` | batch correction | 默认 |
| bbknn neighbors_within_batch | 1 | BBKNN | 默认 |
| bbknn trim | 60 | BBKNN | 默认 |
| leiden_resolution | 1.0 | integrated cluster | 默认 |
| only_rna min CIMA L1 score | 0.5 | inclusion filter | 低置信细胞保留 metadata，但不参与 UMAP/Leiden |

## 需要后续统一到 config 的参数

优先级较高：

- data root 和 output root。
- dataset exception rules。
- ATAC `--nmads`。
- ATAC barcode inference 参数。
- longevity barcode preprocessing 参数。
- product integration 参数。
- R 环境和 reference 文件路径。
