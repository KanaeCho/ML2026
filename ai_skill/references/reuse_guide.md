# Reuse Guide

本文件说明如何把本 skill 迁移到新 RNA/ATAC 项目。

## 可以复用的部分

- discover -> status -> smoke sample -> batch run -> audit 的执行顺序。
- RNA QC、embedding、annotation 的审计思路。
- ATAC barcode、QC、LSI、CIMA projection 的审计思路。
- 输出契约检查方式。
- troubleshooting 顺序。
- Git 中不上传大型数据和结果的规则。

## 不能直接照搬的部分

- 具体 GSE/GSM 样本列表。
- 本机 data root fallback。
- 已完成样本数和真实输出状态。
- CIMA reference 资产路径，除非新项目也使用同一 reference。
- `GSE226039`、`GSE198533`、`GSE206284` 等数据集特例。
- longevity 分支，除非新项目也有同类独立通道。

## 迁移步骤

1. 复制 `ai_skill/` 到新项目或作为上下文提供给 AI。
2. 在新项目写一个项目专用 `AGENTS.md`，只记录新项目事实。
3. 修改样本发现规则和 data root。
4. 确认新项目是否有 RNA、ATAC、co、longevity、product integration 中的哪些通道。
5. 先实现 discover/status。
6. 用一个 smoke sample 验证每条通道。
7. 再运行批处理。
8. 用 output contracts 审计结果。

## 给 AI 的复用提示词

```text
请读取 ai_skill/SKILL.md，并按其中的 single-cell RNA/ATAC batch workflow 帮我处理当前项目。
先只做样本发现、输入检查和 smoke sample 计划，不要直接全量运行。
不要上传 data/output/h5ad/rds/fragment/matrix 文件。
```
