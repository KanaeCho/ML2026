# Batch Audit Report Template

## 批次信息

- 数据通道：`<RNA / ATAC / co / longevity / product>`
- 数据集：`<GSE or dataset>`
- 输出根：`<output root>`
- 审计时间：`<timestamp>`

## 样本状态

| 状态 | 数量 |
| --- | ---: |
| total |  |
| success |  |
| failed |  |
| skipped |  |
| dry_run |  |

## 关键指标

| sample | n_total | n_pass_qc | final_output | annotation_status | outputs_complete | note |
| --- | ---: | ---: | ---: | --- | --- | --- |

## 缺失输出

| sample | missing files |
| --- | --- |

## 失败样本

| sample | failed stage | log path | first error |
| --- | --- | --- | --- |

## 可视化检查

- QC overview：`pass / flag / fail`
- UMAP readability：`pass / flag / fail`
- annotation confidence：`pass / flag / fail`

## 结论

- 是否可进入下一步：`yes/no`
- 是否需要小范围重跑：`yes/no`
- 是否需要用户确认：`yes/no`
