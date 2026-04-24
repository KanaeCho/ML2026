# Dynamic scRNA QC 升级计划

## 目标

将当前单样本 scRNA-seq 主线中的静态 QC 阈值，升级为**数据驱动的 dynamic hybrid MAD-based QC**，并直接替代现有 static QC 主线行为。

升级后必须满足：

1. `n_counts` / `n_genes` 在 `log10(x + 1)` 空间做 lower-tail 动态阈值。
2. `pct_mt` / `pct_ribo` 在原值空间做 upper-tail 动态阈值。
3. `doublet` 继续作为独立 hard gate。
4. 保留现有 fail flags / `pass_qc` 契约不变。
5. 增加样本级阈值审计产物，保证可复核。
6. 保持下游 embedding / annotation / outputs 的主语义稳定。

## 已锁定决策

### 主线行为

- 不再保留 static QC 作为主线路径。
- 主线直接切换为 `dynamic_hybrid_mad`。
- 不做“static 与 dynamic 并行模式”作为长期主设计。

### 阈值语义

- `n_counts`: `log10(n_counts + 1)` lower-tail MAD 阈值。
- `n_genes`: `log10(n_genes + 1)` lower-tail MAD 阈值。
- `pct_mt`: 原值空间 upper-tail MAD 阈值。
- `pct_ribo`: 原值空间 upper-tail MAD 阈值。

### QC 契约

以下列名必须保持不变：

- `fails_count_floor`
- `fails_gene_floor`
- `fails_mt_ceiling`
- `fails_ribo_ceiling`
- `fails_doublet`
- `pass_qc`

### 输出契约

- `metadata.csv` 继续保留全细胞。
- `metadata_qc.csv` 继续仅保留 `pass_qc == True` 细胞。
- 新增样本级阈值审计文件：`qc_thresholds.json`。
- `qc_summary.csv` 增加最终生效阈值字段。

## 当前基线

当前 QC 由以下文件驱动：

- `scripts/only_rna/qc.py`
- `scripts/only_rna/models.py`
- `scripts/only_rna/config.py`
- `scripts/only_rna/default_config.yaml`

当前静态默认阈值：

- `min_counts = 500`
- `min_genes = 300`
- `max_pct_mt = 20.0`
- `max_pct_ribo = 60.0`

当前主线顺序：

1. `read_sample_input(...)`
2. `compute_qc_metrics(...)`
3. `run_doublet_detection(...)`
4. `apply_qc_filters(...)`
5. `run_embedding(...)`
6. `annotate_with_all_versions(...)`
7. `write_sample_outputs(...)`

## 非目标

本阶段不做：

- 跨样本联合 QC
- doublet 方法学重做
- mixture model / 概率模型 QC
- 跨 assay 统一 QC 框架
- 输出目录结构的大规模重构

## 目标方案

### 1. 动态阈值核心

#### `n_counts` / `n_genes`

对每个样本分别计算：

- `log10(n_counts + 1)`
- `log10(n_genes + 1)`

使用：

- `median`
- `MAD = median(|x - median(x)|)`

得到 lower-tail 阈值：

`lower_bound = median - k * MAD`

然后映射回原始空间得到：

- `final_min_counts`
- `final_min_genes`

#### `pct_mt` / `pct_ribo`

在原值空间使用：

`upper_bound = median + k * MAD`

得到：

- `final_max_pct_mt`
- `final_max_pct_ribo`

### 2. guardrail 机制

动态阈值不能直接裸用，需要做 bounded guardrails：

- `min_counts` 不得低于全局下限，也不得高于合理上限。
- `min_genes` 不得低于全局下限，也不得高于合理上限。
- `max_pct_mt` 必须被 clamp 到预设区间。
- `max_pct_ribo` 必须被 clamp 到预设区间。

此外必须处理：

- 小样本
- `MAD == 0`
- 指标全常数/近常数
- 没有 mt/ribo 基因匹配
- 极端筛选导致 `pass_qc_fraction` 过低

注意：这些 fallback 仍然属于 dynamic family 内部稳定化，不是回退到旧 static 主线。

### 3. doublet 规则

保留现有独立 hard gate：

- `fails_doublet = is_doublet`
- `pass_qc` 仍然是所有 fail flag 的否定合取

## 配置设计

建议将 `qc:` 配置切换为动态参数风格：

```yaml
qc:
  method: dynamic_hybrid_mad
  counts_lower_nmads: 3.0
  genes_lower_nmads: 3.0
  pct_mt_upper_nmads: 3.0
  pct_ribo_upper_nmads: 3.5
  min_cells_for_dynamic: 50

  count_floor_min: 100
  count_floor_max: 1500
  gene_floor_min: 100
  gene_floor_max: 1200
  pct_mt_ceiling_min: 5.0
  pct_mt_ceiling_max: 40.0
  pct_ribo_ceiling_min: 20.0
  pct_ribo_ceiling_max: 80.0
```

说明：

- 不再用 `min_counts=500` 这类静态主线字段表达主设计。
- 如果需要兼容旧字段，最多作为 legacy 兼容入口，不再作为主线语义。

## 审计输出设计

### 新增文件

- `output/rna/{GSE}/{sample_id}/qc_thresholds.json`

### 建议结构

```json
{
  "sample_id": "GSMxxxx",
  "gse": "GSExxxx",
  "method": "dynamic_hybrid_mad",
  "n_cells_total": 12345,
  "metrics": {
    "n_counts": {
      "transform": "log10p1",
      "direction": "lower",
      "center": 3.12,
      "mad": 0.18,
      "nmads": 3.0,
      "raw_threshold": 2.58,
      "final_threshold_original_scale": 379,
      "guardrails_applied": []
    },
    "n_genes": {},
    "pct_mt": {},
    "pct_ribo": {}
  },
  "fallbacks": {
    "small_sample_rule_used": false,
    "zero_mad_metrics": [],
    "retention_guard_triggered": false
  }
}
```

### `qc_summary.csv` 增加字段

- `qc_threshold_method`
- `final_min_counts`
- `final_min_genes`
- `final_max_pct_mt`
- `final_max_pct_ribo`

## 文件改动清单

### 1. `scripts/only_rna/models.py`

目标：新增动态 QC 配置与样本级阈值审计数据结构。

改动：

- 重构/扩展 `QcThresholds`，使其表达动态 MAD + guardrails 参数。
- 新增如 `ComputedQcThresholds` 之类的 dataclass，承载样本实际计算结果。

### 2. `scripts/only_rna/config.py`

目标：让 YAML 和 override 正确加载新的 dynamic QC 配置。

改动：

- `load_run_config()` 改为读取新的 `qc:` 结构。
- `merge_cli_overrides()` 支持新的 nested keys。

### 3. `scripts/only_rna/default_config.yaml`

目标：默认配置直接切换到 dynamic 主线。

改动：

- 用 dynamic MAD + guardrail 参数替代旧静态阈值字段。

### 4. `scripts/only_rna/qc.py`

目标：实现动态阈值计算和 fail flags 契约保持。

改动：

- 保留 `compute_qc_metrics()` 只负责算基础指标。
- 重写 `apply_qc_filters()`：
  - 先计算样本级动态阈值
  - 写入 `adata.uns["qc_thresholds"]`
  - 再生成原契约 fail flags 和 `pass_qc`

新增 helper：

- `_median_and_mad(...)`
- `_compute_lower_tail_mad_threshold_log10(...)`
- `_compute_upper_tail_mad_threshold(...)`
- `_compute_sample_qc_thresholds(...)`

### 5. `scripts/only_rna/outputs.py`

目标：把动态阈值输出为正式审计文件。

改动：

- 新增 `qc_thresholds.json` writer。
- `qc_summary.csv` 写入最终生效阈值字段。
- `validation_result.csv` 可增加阈值审计存在性检查。

### 6. `scripts/only_rna/cli.py`

目标：样本完成度检查纳入新文件。

改动：

- `_expected_output_paths()` 新增 `qc_thresholds.json`。

### 7. `scripts/only_rna/plotting.py`

目标：QC 图要使用样本实际阈值，而不是旧静态配置值。

改动：

- `save_qc_overview()` 优先从 `adata.uns["qc_thresholds"]` 读取阈值。
- 仅在缺失时才 fallback 到 config。

### 8. `scripts/only_rna/tuning_presets.py`

目标：baseline-only tuning 路径和新主线保持一致。

改动：

- baseline preset 切换为 dynamic QC 配置，而不是旧静态阈值 dataclass。

### 9. `scripts/only_rna/tuning_orchestrator.py`

目标：确保 tuning 执行的就是新的 dynamic 主线。

改动：

- 校验 baseline candidate 用的是新 QC 配置。
- 确保新 audit 文件出现在 tuning/mainline 路径中（按当前输出契约而定）。

### 10. `AGENTS.md`

目标：文档反映真实行为。

改动：

- 删除静态阈值主线描述。
- 更新为 dynamic hybrid MAD 真实行为。
- 记录新增 `qc_thresholds.json` 和 summary 字段。

## 执行 Wave

### Wave 1：配置与数据结构

文件：

- `models.py`
- `config.py`
- `default_config.yaml`
- `tests/only_rna/test_config.py`

目标：

- 建好 dynamic QC 配置层
- 保证 YAML / override 可以驱动新方案

验证：

```bash
pytest tests/only_rna/test_config.py -q
```

### Wave 2：QC 核心数学逻辑

文件：

- `qc.py`
- `tests/only_rna/test_processing.py`

目标：

- 动态阈值正确计算
- guardrails 生效
- fail flags / `pass_qc` 契约保持

验证：

```bash
pytest tests/only_rna/test_processing.py -q
```

### Wave 3：审计输出与契约同步

文件：

- `outputs.py`
- `cli.py`
- `plotting.py`
- `tests/only_rna/test_outputs.py`
- `tests/only_rna/test_cli.py`

目标：

- 输出 `qc_thresholds.json`
- `qc_summary.csv` 扩字段
- 输出完整性契约纳入新文件
- 图上阈值与实际应用阈值一致

验证：

```bash
pytest tests/only_rna/test_outputs.py tests/only_rna/test_cli.py -q
```

### Wave 4：baseline-only tuning 与文档收尾

文件：

- `tuning_presets.py`
- `tuning_orchestrator.py`
- `tests/only_rna/test_tuning.py`
- `AGENTS.md`

目标：

- tuning 与主线保持一致
- 文档更新为真实行为

验证：

```bash
pytest tests/only_rna/test_tuning.py -q
pytest tests/only_rna -q
```

## 必须新增/修改的测试

### `tests/only_rna/test_processing.py`

新增：

- `test_apply_qc_filters_computes_counts_threshold_in_log_space`
- `test_apply_qc_filters_computes_genes_threshold_in_log_space`
- `test_apply_qc_filters_computes_mt_threshold_in_original_space`
- `test_apply_qc_filters_computes_ribo_threshold_in_original_space`
- `test_apply_qc_filters_applies_guardrail_bounds`
- `test_apply_qc_filters_handles_zero_mad_without_crash`
- `test_apply_qc_filters_handles_small_sample_without_crash`
- `test_apply_qc_filters_preserves_doublet_hard_gate`

### `tests/only_rna/test_outputs.py`

新增/修改：

- `test_write_sample_outputs_writes_qc_thresholds_json`
- `test_qc_summary_includes_final_dynamic_thresholds`
- `test_qc_overview_uses_sample_specific_thresholds`

### `tests/only_rna/test_tuning.py`

新增/修改：

- baseline-only tuning 仍使用新 dynamic QC config
- 审计输出存在

## 风险与缓解

### 风险 1：小样本阈值不稳定

缓解：

- `min_cells_for_dynamic`
- `MAD == 0` fallback
- threshold clamp

### 风险 2：动态阈值过严导致样本几乎被清空

缓解：

- guardrail 上下限
- retention 过低时触发内部 fallback 规则
- 写入 `fallbacks.retention_guard_triggered`

### 风险 3：图和真实阈值不一致

缓解：

- plotting 只读 `adata.uns["qc_thresholds"]`

### 风险 4：旧测试大量依赖静态默认值

缓解：

- 优先更新 `test_config.py` 与 `test_processing.py`
- 其余测试逐步去除对旧固定数值的依赖

## 验收标准

1. 主线不再以静态阈值为默认行为。
2. `n_counts / n_genes` 在 `log10(x+1)` 空间算 lower-tail 动态阈值。
3. `pct_mt / pct_ribo` 在原值空间算 upper-tail 动态阈值。
4. `doublet` 仍然是独立 hard gate。
5. fail flags / `pass_qc` 列名与布尔语义保持不变。
6. 每个样本输出 `qc_thresholds.json`。
7. `qc_summary.csv` 包含最终动态阈值字段。
8. `metadata.csv` / `metadata_qc.csv` 的 all-cells vs pass-QC-only 语义保持不变。
9. embedding / annotation 仍只依赖 `pass_qc`。
10. `pytest tests/only_rna -q` 通过。
11. `AGENTS.md` 已更新为真实行为。

## 执行建议

按 Wave 顺序执行，不建议跳步。

推荐顺序：

1. Wave 1
2. Wave 2
3. Wave 3
4. Wave 4

只有在前一 Wave 的短测试通过后，才进入下一 Wave。
