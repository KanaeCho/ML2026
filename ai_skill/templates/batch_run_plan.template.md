# Batch Run Plan Template

## 目标

- 数据通道：`RNA / ATAC / co / longevity / product integration`
- 数据集：`<GSE or dataset>`
- 样本范围：`<all or selected samples>`

## 输入检查

- data root：`<path or env>`
- 样本表：`<datasets.xlsx / atac.xlsx / co.xlsx / raw scan>`
- reference：`<peak.bed / CIMA / Azimuth>`
- raw 文件：`<present / missing>`

## Discover 命令

```bash
<command>
```

## Status 命令

```bash
<command>
```

## Smoke Sample

- 样本：`<sample>`
- 命令：

```bash
<command>
```

## 批处理命令

```bash
<command>
```

## 审计清单

- `run_status.json`
- `qc_summary.csv`
- `validation_result.csv`
- UMAP PNG
- `.h5ad`
- `logs/sample_qc.log`

## 风险和确认点

- 是否覆盖旧输出：`yes/no`
- 是否运行大批量任务：`yes/no`
- 是否需要用户确认：`yes/no`
