# AI Skill 说明

本目录用于存放给 AI 使用的 workflow/skill 说明。

AI skill 只是这个 Git 项目的一部分，不是项目主体。给人看的主要说明应放在 `README.md` 和 `docs/` 中。

计划结构：

```text
ai_skill/
├── SKILL.md
└── references/
    ├── rna_pipeline.md
    ├── atac_pipeline.md
    ├── batch_processing.md
    ├── diagnostic_workflow.md
    ├── troubleshooting.md
    ├── knowledge_map.md
    └── case_notes_template.md
```

AI skill 应告诉 AI：

- 如何读取本项目结构。
- 哪些文件是主线代码。
- 哪些内容不能在未确认时修改。
- 如何检查 RNA/ATAC 输出是否完整。
- 如何排查 QC、聚类、UMAP、注释、batch effect 和批处理失败。
- 如何记录 case note。
- 如何辅助生成 PPT 提纲和项目文档。
