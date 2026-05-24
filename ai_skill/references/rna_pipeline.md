# AI 参考：RNA Pipeline

详见人类文档 `docs/rna_pipeline.md`。AI 回答或排查 RNA 问题时，应优先引用该文档中的当前项目事实。

核心事实：

- 当前 RNA 主线在 `scripts/only_rna/`。
- 总入口在 `scripts/process/pipeline.py`。
- 当前主注释方法是 Azimuth `pbmcref`。
- 输出 final celltype 为 5 类：`CD4_T`, `CD8_T`, `B`, `Myeloid`, `NK`。
- marker gene 主线输出当前项目中未找到。
