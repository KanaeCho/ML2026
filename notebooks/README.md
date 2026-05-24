# Notebooks 说明

本目录用于存放经过整理和确认的 notebook。

当前仓库中已有 notebook 仍需要人工分类后再移动：

- `scripts/scATAC.ipynb`
- `scripts/step01_signac_single_sample_covid-redo.ipynb`
- `scripts/process/GSM8671454.ipynb`

建议分类：

- `notebooks/rna/`：RNA 相关的正式或说明性 notebook。
- `notebooks/atac/`：ATAC 相关的正式或说明性 notebook。
- `notebooks/exploration/`：探索性、临时分析或历史版本 notebook。

提交 notebook 前必须检查：

- 清理大型输出。
- 删除本地绝对路径。
- 删除不应公开的样本信息。
- 确认 notebook 是正式流程说明、探索记录还是历史文件。
