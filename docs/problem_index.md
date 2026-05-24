# 问题索引

本文档是 RNA/ATAC 项目的问题类型索引。这里列出的很多问题是未来可能遇到的潜在问题，不代表当前项目都已经发生。

| 问题类型 | 常见表现 | 可能原因 | 优先检查文件 | 检查顺序 | 可尝试处理 | 是否可能需要重跑 | 是否需要人工确认 | 当前项目是否出现 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| RNA QC 异常 | pass-QC 比例过低或 QC 图异常 | 样本质量差、阈值过严、doublet 高、mt/ribo 异常 | `qc_summary.csv`, `qc_thresholds.json`, `qc_overview.png` | 先看细胞数，再看阈值和 fail flags | 记录问题；必要时单样本检查参数 | 可能 | 是 | 潜在问题 |
| ATAC QC 异常 | TSS/FRiP 低、pass-QC 少 | fragment 质量差、barcode 选择不合适、peak reference 问题 | `qc_summary.csv`, `metadata.csv`, `qc_overview.png`, log | 先看 barcode_source，再看 TSS/FRiP/fragments | 小范围检查 barcode 参数 | 可能 | 是 | 潜在问题 |
| 细胞数异常减少 | 过滤后细胞数过少 | QC 阈值、doublet、输入矩阵不对 | QC summary、metadata、run log | 比较 total/pass/final cells | 单样本检查，不先全量重跑 | 可能 | 是 | 潜在问题 |
| 样本间差异过大 | UMAP 按样本分开 | 生物差异、batch effect、样本质量差 | product figures, sample QC summary | 先看 QC，再看 batch/sample panels | 分层解释或检查 integration 参数 | 可能 | 是 | 潜在问题 |
| batch effect 明显 | product UMAP 被 GSE/sample 主导 | 批次效应、batch key 不合适 | integration metrics, sample mixing summary | 看 method、batch key、mixing | 可能小范围重整合 | 可能 | 是 | 潜在问题 |
| 聚类过粗或过细 | cluster 太少或太碎 | resolution、neighbors、样本质量 | metadata, UMAP, 参数 | 查 cluster 数和 resolution | 调参数需先确认 | 可能 | 是 | 潜在问题 |
| UMAP/投影混乱 | 细胞云不清晰或桥状结构 | 低质量细胞、低置信投影、batch effect | UMAP 图、metadata、integration summary | 确认坐标来源和 inclusion filter | 可能小范围重算 UMAP | 可能 | 是 | 潜在问题 |
| cluster 无明显 marker | marker 不清楚 | 聚类不佳、marker 流程缺失、样本组成复杂 | marker 表、metadata、UMAP | 先确认 marker 文件是否存在 | 当前先记录缺失 | 可能 | 是 | 当前项目中未找到 marker 主线 |
| marker gene 和预期不一致 | marker 指向其他细胞类型 | 注释错误、cluster 混合、marker 表版本问题 | marker 表、annotation table | 查 marker 生成脚本和 label | 需人工判断 | 可能 | 是 | 潜在问题 |
| 注释图不好或 label 不清楚 | UMAP label 混乱、Unknown 多 | annotation 低置信、label 映射问题 | UMAP、metadata_qc、qc_summary | 查 annotation status 和 score | 先修文档/label 解释 | 可能 | 是 | 潜在问题 |
| annotation table 和 metadata 不一致 | 表格和图 label 对不上 | 字段版本混乱、旧输出 | metadata, validation, plotting script | 查绘图使用字段 | 可能只需重画图 | 可能 | 是 | 潜在问题 |
| RNA/ATAC 注释不一致 | 同样本 RNA 和 ATAC label 不一致 | 模态差异、barcode 不匹配、注释体系不同 | co RNA/ATAC metadata | 先确认是否 paired matching | 多数需解释而非重跑 | 可能 | 是 | 潜在问题 |
| ATAC gene activity 异常 | gene activity 不符合预期 | 当前未确认有 gene activity 输出 | gene activity matrix | 先确认文件是否存在 | 未确认 | 是 | 当前项目中未找到 |
| peak matrix 或 fragment 文件异常 | FeatureMatrix 失败或 cells 很少 | fragment/index/peak 文件错误 | log, fragment path, peak.bed | 查文件存在和 tabix | 可能单样本重跑 | 可能 | 是 | 潜在问题 |
| 多样本合并后某个样本主导 | UMAP 主要由一个样本组成 | 样本细胞数不均衡、质量差异 | cells_metadata, sample_mixing_summary | 查每样本细胞数 | 可能分层展示 | 可能 | 是 | 潜在问题 |
| 批处理部分样本失败 | 部分 run_status failed | 输入缺失、环境错误、输出权限 | run_status, log | 定位样本和阶段 | 单样本重跑 | 可能 | 是 | 潜在问题 |
| 输出文件缺失 | expected output 不完整 | 脚本失败、写权限、旧 profile | validation_result, run_status | 查 missing 文件 | 补跑或修导出 | 可能 | 是 | 潜在问题 |
| 路径硬编码 | 新机器无法运行 | 本地路径写死 | grep `/mnt/`, config | 查代码和模板 | 配置化 | 不一定 | 是 | 已出现 |
| 参数没有记录 | 无法复现实验 | 参数分散在脚本/CLI/log | config, run_status, logs | 找参数来源 | 补文档 | 不一定 | 是 | 已出现部分风险 |
| 结果版本混乱 | 多个目录不知道用哪个 | tuning/legacy/product 混合 | run_status, manifests | 查 source_output_dir | 记录版本索引 | 不一定 | 是 | 潜在问题 |
