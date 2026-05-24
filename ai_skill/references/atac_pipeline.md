# AI 参考：ATAC Pipeline

详见人类文档 `docs/atac_pipeline.md`。

核心事实：

- 当前 ATAC 主线是 `scripts/process/process_single_sample.R`。
- only ATAC、co ATAC 和 longevity ATAC 都复用该主线或其包装逻辑。
- 当前 ATAC 注释方法是 CIMA ATAC compact feature model + centroid assignment。
- gene activity 和 motif 主线输出当前项目中未找到。
