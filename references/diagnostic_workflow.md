# 异常结果诊断流程

本文档用于以后遇到“结果不好”“图不好”“注释不对”“聚类奇怪”“样本分开”“marker 不明显”“某一步失败”等问题时，按顺序排查。原则是先定位问题来源，再决定是否需要小范围重跑；不要一开始就直接重跑全流程。

## 总原则

1. 先确认是 RNA、ATAC、co、longevity 还是 product-level 整合问题。
2. 先确认是单样本问题还是多样本合并问题。
3. 先看输入、metadata、barcode/sample ID、QC 后细胞数和日志。
4. 再看聚类、UMAP、注释、marker 或 batch effect。
5. 修改前先记录证据和判断。
6. 能小范围检查就不要大范围重跑。
7. 需要人工判断的地方必须标注。

## 第一步：识别问题类型

常见类型：

- RNA QC 异常。
- ATAC QC 异常。
- 过滤后细胞数过少。
- 样本之间差异过大。
- batch effect 明显。
- 聚类过粗或过细。
- UMAP/t-SNE 投影混乱。
- cluster 分得不清楚。
- marker gene 不明显。
- marker gene 和预期细胞类型不一致。
- 注释图不好或 label 不清楚。
- annotation table 和 metadata 对不上。
- RNA 和 ATAC 注释不一致。
- ATAC gene activity 异常。
- peak matrix 或 fragment 文件异常。
- 多样本合并后某个样本主导结果。
- 批处理部分样本失败。
- 输出文件缺失。
- 路径硬编码。
- 参数没有记录。
- 结果版本混乱。

## 第二步：定位相关文件

每次诊断至少定位以下文件：

| 文件类型 | RNA | ATAC | product-level |
| --- | --- | --- | --- |
| 运行状态 | `run_status.json` | `run_status.json` | `product_status.json` |
| 日志 | `logs/sample_qc.log` | `logs/sample_qc.log` | 命令输出或 product status |
| QC summary | `qc_summary.csv` | `qc_summary.csv` | `qc/sample_qc_summary.csv` |
| metadata | `metadata.csv`, `metadata_qc.csv` | `metadata.csv`, `metadata_qc.csv` | `manifests/cells_metadata.csv` |
| validation | `validation_result.csv` | `validation_result.csv` | `qc/validation_summary.csv` |
| UMAP 图 | `umap_rna_*.png` | `umap_cima_cell_type_*.png` | `figures/*_panels.png` |
| 中间对象 | `{sample_id}.h5ad` | `{sample_id}.h5ad`, 可能有 RDS | `{product}.h5ad` |
| 参数 | `qc_thresholds.json`, config | `qc_summary.csv`, command, log | `integration/integration_summary.json` |

## 第三步：按顺序排查

1. 数据类型：RNA 还是 ATAC。
2. 分析范围：单样本还是多样本整合。
3. 发生阶段：QC、过滤、降维、聚类、投影、注释、导出还是批处理。
4. 输入文件是否存在且路径正确。
5. metadata、barcode、sample ID 是否匹配。
6. QC 后细胞数是否合理。
7. 是否存在低质量样本或低质量 cluster。
8. 聚类参数是否记录。
9. 投影是否被 sample、GSE、donor 或 batch 主导。
10. marker 是否清楚，若当前项目中没有 marker 输出则标注“当前项目中未找到”。
11. 注释规则是否可追溯。
12. 绘图使用的 label 是否是最新版本。
13. 日志是否有报错或 fallback。
14. 是否需要小范围重跑。
15. 是否需要人工确认。

## 第四步：区分问题来源

| 来源 | 判断依据 |
| --- | --- |
| 数据质量问题 | QC 指标异常、过滤后细胞数少、FRiP/TSS/mt/ribo 异常 |
| 参数设置问题 | 阈值过严/过松、cluster resolution 不合适、barcode 参数不合适 |
| 批处理逻辑问题 | 某些样本失败、路径错误、跳过逻辑异常 |
| metadata 对应关系问题 | sample ID、barcode、donor、health 对不上 |
| 注释规则问题 | label 映射不清楚、low confidence 高、annotation status 非 ok |
| 绘图或导出问题 | metadata 有值但图中 label 不对、输出缺失 |
| 版本混乱问题 | 同一样本有多个输出目录、tuning candidate 和 sample-root 混淆 |
| 证据不足 | 缺少日志、缺少参数、缺少旧版本记录，需要人工确认 |

## 第五步：处理建议分级

| 等级 | 含义 | 示例 |
| --- | --- | --- |
| 可以直接文档化 | 不改结果，只补说明 | 记录某个字段含义、记录路径模式 |
| 可以小范围检查 | 读取 metadata/QC/log，不重跑 | 检查某样本 `qc_summary.csv` 和 UMAP |
| 需要用户确认后修改 | 可能改变流程或解释 | 修改 label 映射、调整 dataset rule |
| 需要小范围重跑 | 只重跑单样本或少数样本 | 修改某样本 barcode 参数后重跑 |
| 不建议当前执行 | 成本高或证据不足 | 全量重跑、重建 reference、重做整合 |

## 诊断记录要求

每次排查都建议用 `templates/case_note_template.md` 记录：

- 问题类型。
- 问题现象。
- 涉及样本和文件。
- 初步判断。
- 检查顺序。
- 检查结果。
- 是否重跑。
- 当前结论。
- 后续待确认。
