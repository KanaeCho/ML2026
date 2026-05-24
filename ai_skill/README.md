# AI Skill 说明

本目录存放给 AI 使用的 single-cell RNA/ATAC workflow skill。它的用途不是普通项目介绍，而是把当前项目已经实现的 RNA、ATAC、co、longevity 和 product integration 流程整理成可复用的执行手册。

如果以后处理类似数据，可以把本目录作为上下文给 AI，让 AI 按 `SKILL.md` 中的 discover -> status -> smoke sample -> batch run -> audit 顺序执行。

当前结构：

```text
ai_skill/
├── SKILL.md
├── references/
│   ├── rna_workflow.md
│   ├── atac_workflow.md
│   ├── batch_workflow.md
│   ├── output_contracts.md
│   ├── command_reference.md
│   ├── troubleshooting.md
│   └── reuse_guide.md
└── templates/
    ├── batch_run_plan.template.md
    └── batch_audit_report.template.md
```

AI skill 应告诉 AI：

- 如何判断数据通道。
- 如何发现样本和检查状态。
- RNA 如何执行 QC、聚类、UMAP、Azimuth/CIMA 映射和 final celltype 输出。
- ATAC 如何执行 barcode 选择、QC、LSI/UMAP、CIMA 投影和注释输出。
- 多样本如何批处理。
- product-level integration 何时运行。
- 如何审计输出完整性。
- 如何排查 QC、聚类、UMAP、注释和批处理失败。

注意：`AGENTS.md` 是当前 ML2026 项目专用上下文；`ai_skill/` 是可复用 workflow skill。不要把已完成样本数、本地绝对路径或 ML2026 特例直接当作新项目事实照搬。
