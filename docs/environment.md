# 运行环境说明

本文档记录当前项目代码运行所需的主要 Python 与 R 环境。这里列出的是当前仓库中可见代码的依赖线索，不等同于完全锁定的生产环境快照。

## Python 环境

- 当前 `pyproject.toml` 要求 Python `>=3.11`。
- 推荐使用 `uv` 管理环境和运行命令。
- 常用命令格式：`uv run python scripts/process/pipeline.py <command> ...`。
- 开发测试依赖包含 `pytest`。

当前 `pyproject.toml` 中记录的主要 Python 包包括：

- `anndata`
- `scanpy`
- `bbknn`
- `harmonypy`
- `numpy`
- `pandas`
- `scipy`
- `scikit-learn`
- `matplotlib`
- `umap-learn`
- `h5py`
- `openpyxl`
- `pyyaml`
- `igraph`
- `leidenalg`
- `pyranges`
- `torch`
- `torchvision`
- `tqdm`

## R 环境

当前仓库中没有发现 `renv.lock` 或等价 R lockfile，因此 R 包版本尚未在仓库中锁定。以下依赖来自 `scripts/**/*.R` 中的 `library(...)` 调用。

主要 ATAC/RNA 脚本用到的 R 包包括：

- `optparse`
- `Signac`
- `Seurat`
- `GenomeInfoDb`
- `EnsDb.Hsapiens.v86`
- `scDblFinder`
- `SingleCellExperiment`
- `Matrix`
- `rtracklayer`
- `ggplot2`
- `patchwork`
- `data.table`
- `jsonlite`
- `ArchR`
- `RhpcBLASctl`
- `harmony`
- `SoupX`
- `tidyverse`

## 外部参考和数据根

- 当前代码优先读取工作区下的 `./data`。
- 若 `./data` 不存在，则读取环境变量 `ML2026_DATA_ROOT`。
- 若环境变量也不存在，则回退到项目约定的外部数据根。
- Git 仓库不应上传真实 raw data、`.h5ad` 大对象、fragment 文件、完整 output 目录或本地绝对路径配置。

## Azimuth 与 Seurat 依赖

- RNA 主线的主要注释路径调用 R/Seurat/Azimuth 的 `pbmcref` reference。
- 运行 RNA 注释前需确认 R 中可加载 Seurat/Azimuth 相关依赖，并且本机能访问所需 reference。
- 若 Azimuth 失败，先检查 R 包、reference 下载/缓存、输入矩阵 gene symbol 和运行日志，不建议直接全量重跑。

## ArchR 依赖

- longevity ATAC barcode 预处理使用 `scripts/longevity/preprocess_atac_barcodes_archr.R`。
- 运行前需确认 R 中可加载 `ArchR`、`data.table`、`jsonlite` 和 `optparse`。
- ArchR 相关日志目录和大中间文件不应上传 Git。

## 测试

常用测试命令：

```bash
uv run pytest
```

若只验证 RNA tuning 相关改动，可运行：

```bash
uv run pytest tests/only_rna/test_tuning.py
```

## 当前未完全确认事项

- R 包版本尚未锁定。
- notebooks 的实际依赖尚未逐个审计。
- 不同机器上的 Azimuth reference 缓存位置尚未标准化。
- GPU 不是当前 RNA/ATAC 主线验收的必要条件，但 `torch` 已在 Python 依赖中出现。
