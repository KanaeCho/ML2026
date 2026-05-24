# 待办事项和不确定项

本文档记录当前项目中需要人工确认、后续整理或后续重构的内容。这里的内容不代表流程错误，只表示在项目上传 Git 或长期复用前需要进一步确认。

## 需要人工确认的内容

| 内容 | 当前状态 | 为什么需要确认 |
| --- | --- | --- |
| Git 仓库是公开还是私有 | 未确认 | 决定本地路径、样本状态和内部说明能保留多少。 |
| `QualityControl/` 的正式用途 | 未确认 | 目录中包含 `.rds` 和 PDF，可能适合做示例，但不能整体上传。 |
| `ArchRLogs/` 的正式用途 | 未确认 | 日志可用于排查问题，但可能包含本地路径，不应整体上传。 |
| notebook 的状态 | 未确认 | 需要判断哪些是正式流程、探索记录、旧版本或废弃文件。 |
| RNA marker gene 输出 | 当前项目中未找到 | 需要确认是否存在于 notebook 或外部结果目录。 |
| ATAC gene activity 输出 | 当前项目中未找到 | 需要确认是否在主线之外生成过。 |
| ATAC motif 输出 | 当前项目中未找到 | 需要确认是否在主线之外生成过。 |
| t-SNE 输出 | 当前项目中未找到 | 当前主线以 UMAP 为主。 |
| R 环境配置 | 未确认 | R 包已在脚本中使用，但还没有统一安装说明或环境文件。 |
| 完整 reference 文件清单 | 部分确认 | CIMA 和 peak 文件是必需的，但公开 setup 文档前需要确认完整清单。 |

## Git 准备待办

| 待办 | 优先级 | 说明 |
| --- | --- | --- |
| 扩展 `.gitignore` | 已完成第一版 | 已排除大数据、大结果、日志和本地配置；后续可继续补充。 |
| 检查 `.claude_resources.json` | 高 | 该文件包含本地路径，已加入 ignore。 |
| 检查 `.opencode/` 和 `opencode.json` | 暂不上传 | 已按公开仓库安全标准加入 ignore；若团队需要共享，再人工审查后取消。 |
| 清理或移动 notebook | 暂不上传 | 已默认 ignore `*.ipynb`；清理输出、本地路径和隐私信息后再选择性提交。 |
| 修复测试里的本地绝对路径 | 已完成 | `tests/only_rna/test_tuning.py` 已改为基于仓库位置解析 `default_config.yaml`。 |
| 补充 R 环境说明 | 已完成第一版 | 已新增 `docs/environment.md`，但 R 包版本尚未锁定。 |
| 编写 RNA pipeline 文档 | 高 | 计划文件为 `docs/rna_pipeline.md`。 |
| 编写 ATAC pipeline 文档 | 高 | 计划文件为 `docs/atac_pipeline.md`。 |
| 编写批处理文档 | 高 | 计划文件为 `docs/batch_processing.md`。 |
| 编写 troubleshooting 和 diagnostic workflow | 高 | 用于长期复用和 AI 辅助排查。 |
| 编写 AI skill 草案 | 中 | 计划放在 `ai_skill/`，这是辅助内容，不是项目主体。 |

## 后续代码整理候选项

以下内容不建议在当前文档整理阶段立即修改，应在有测试保护后逐步处理。

| 候选项 | 当前证据 | 后续建议 |
| --- | --- | --- |
| 统一 data root 解析 | 多个模块都有 data root 解析逻辑 | 当前先保留，后续集中为一个 canonical helper。 |
| dataset 特例配置化 | GSE 特例写在代码里 | 后续整理成 dataset policy 表或 YAML。 |
| 拆分 ATAC R 脚本 | `process_single_sample.R` 功能较多 | 先文档化，后续按 input/QC/annotation/export 拆分。 |
| 区分 legacy 和 mainline 脚本 | 已完成第一版 | 已新增 `scripts/README.md` 分类说明；后续可继续细化每个辅助脚本的状态。 |
| 统一参数记录 | 参数分散在 YAML、argparse 和 R options 中 | 先写 `docs/parameters.md`，再考虑集中到 config。 |

## 当前已知风险

| 风险 | 状态 | 建议处理 |
| --- | --- | --- |
| 原始数据误提交 | 可预防 | 保持 `data/` ignore，并用扩展名 ignore 大文件。 |
| 大型结果误提交 | 可预防 | 保持 `output/` ignore，并 ignore `.h5ad`、`.rds`、matrix、fragment 等。 |
| 公开文档中出现本地路径 | 当前已有内部文档/代码包含本地路径 | 模板使用占位符，公开前审阅 README 和 docs。 |
| 清理时破坏已有流程 | 可预防 | 在文档和测试完善前，不移动或重构核心脚本。 |
| 结果版本混乱 | 潜在风险 | 后续补充 results index 和 run-status 约定。 |

## 后续建议新增文档

第一批文档之后，建议继续补充：

1. `docs/rna_pipeline.md`
2. `docs/atac_pipeline.md`
3. `docs/batch_processing.md`
4. `docs/parameters.md`
5. `docs/results_index.md`
6. `docs/problem_index.md`
7. `docs/troubleshooting.md`
8. `references/knowledge_map.md`
9. `references/qc_metrics.md`
10. `references/annotation_notes.md`
11. `references/diagnostic_workflow.md`
12. `templates/case_note_template.md`
13. `templates/ppt_outline_template.md`
14. `ai_skill/SKILL.md`
