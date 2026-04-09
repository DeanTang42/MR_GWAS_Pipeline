# MR_GWAS_Pipeline

这是一个可独立运行、可用 `git` 管理的 MR 分析项目骨架，包含三部分核心能力：

1. GWAS 摘要统计标准化与 canonical 位点 ID 生成。
2. 基于 `TwoSampleMR` 的单暴露单结局 MR 分析。
3. 基于 `TwoSampleMR + MVMR` 的 `X-M-Y` 中介 MVMR 分析。
4. 可配置 exposure 数量的通用 MVMR 分析。

项目已经按当前机器环境整理好默认配置，代码、配置、文档可以直接纳入 `git`；原始数据、临时文件和结果文件默认被 `.gitignore` 排除。

## 目录结构

```text
MR_GWAS_Pipeline/
├── bin/
├── config/
├── data/
│   ├── Org/
│   ├── reference/
│   ├── standardized/
│   ├── exp/
│   └── out/
├── docs/
├── logs/
├── results/
├── scripts/
├── .gitignore
├── environment.yml
└── requirements.txt
```

## 当前默认环境

- Python: `3.8.20`
- R: `4.3.1`
- Python 包:
  - `polars==1.8.2`
  - `questionary==2.1.0`
  - `rich==14.3.3`

R 包不再由项目自动安装。需要的包名清单见 [config/r_packages.txt](config/r_packages.txt)。

## 配置文件

所有常用参数都集中放在 [config/defaults.env](config/defaults.env)：

- `ORG_DATA_DIR`
- `STANDARDIZED_OUTPUT_DIR`
- `EXP_OUTPUT_DIR`
- `OUT_OUTPUT_DIR`
- `RESULTS_DIR`
- `REFERENCE_DIR`
- `PLINK_BIN`
- `CLUMP_BFILE`
- `FRQ_REFERENCE`
- `CLUMP_R2`
- `CLUMP_KB`
- `CLUMP_P1`
- `CLUMP_POP`

## 快速开始

### 1. 安装 Python 依赖

```bash
bash bin/bootstrap.sh
```

R 包请按 [config/r_packages.txt](config/r_packages.txt) 自行安装。

### 2. 原始 GWAS 的放置位置

建议把原始下载文件放在：

```text
data/Org/
```

交互模式下，标准化脚本会优先把这个目录作为输入路径提示。

### 3. 标准化原始 GWAS

第一步只做标准化，统一生成 `data/standardized/*.tsv.gz`。同一个标准化文件后续可以重复转换成 exposure 或 outcome。

```bash
bash bin/standardize_gwas.sh \
  --non-interactive \
  --input /path/to/gwas.tsv.gz \
  --output /home/ding/MR_GWAS_Pipeline/data/standardized/EXPO.tsv.gz \
  --mode B \
  --snp-col variant_id \
  --chr-col chromosome \
  --pos-col base_pair_location \
  --allele1-col effect_allele \
  --allele2-col other_allele \
  --stat-type BETA \
  --stat-col beta \
  --se-col standard_error \
  --p-col p_value \
  --pval-format raw \
  --freq-type EAF \
  --freq-col effect_allele_frequency
```

### 4. 如有需要，为缺失 EAF 做 MAF 定向注释

如果某个 GWAS 只提供 `MAF`，标准化阶段会保留 `MAF` 并让 `EAF=NA`。这时可以先从 1KG bfile 生成频率参考：

```bash
bash bin/build_1kg_frq_reference.sh
```

默认会使用 `CLUMP_BFILE`，输出：

```text
data/reference/g1000_eur_colon.frq
data/reference/g1000_eur_colon.frq.tsv.gz
```

然后对标准化文件补 EAF：

```bash
bash bin/annotate_eaf.sh \
  --non-interactive \
  --input data/standardized/EXPO.tsv.gz \
  --frq data/reference/g1000_eur_colon.frq.tsv.gz
```

该步骤默认直接替换原标准化文件，并额外生成：

- `*.eaf_matched.tsv.gz`: 仅包含 `BIM_ID` 匹配到 1KG `.frq` 的行
- `*.eaf_frequency_diff.tsv.gz`: `EAF_1KG` 与 `EAF_FROM_MAF` 差异超过阈值的位点
- `*.eaf_annotation.log`: 匹配数、注释数和频率差异统计

注释后的标准化文件会新增：

- `EAF_1KG`: 使用 1KG `.frq` 的 A1/A2/MAF 定向到 `EFFECT_ALLELE` 的频率
- `EAF_FROM_MAF`: 使用原 GWAS 文件的 `MAF`，再按 1KG A1/A2 定向到 `EFFECT_ALLELE` 的频率
- `EAF_ABS_DIFF`: `EAF_1KG` 和 `EAF_FROM_MAF` 的绝对差

### 5. 转换成 exposure 或 outcome

第二步从 `data/standardized` 中读取文件，生成 TwoSampleMR 输入。`role=exposure` 会先按 `clump_p1` 预筛候选位点，再调用 PLINK clump；`role=outcome` 只做格式转换。

生成 outcome：

```bash
bash bin/format_tsmr.sh \
  --non-interactive \
  --input data/standardized/EXPO.tsv.gz \
  --output data/out/EXPO_outcome.csv \
  --role outcome \
  --phenotype EXPO
```

生成 exposure：

```bash
bash bin/format_tsmr.sh \
  --non-interactive \
  --input data/standardized/EXPO.tsv.gz \
  --output data/exp/EXPO_exposure.csv \
  --role exposure \
  --phenotype EXPO \
  --clump-r2 0.1 \
  --clump-kb 500 \
  --clump-p1 5e-8
```

如果 exposure 只希望在指定基因/染色体区间内 clump，可以加：

```bash
--region 1:55000000-56000000
```

非 `--non-interactive` 模式下，`format_tsmr.sh` 会交互式询问：

- `clump_r2`
- `clump_kb`
- `clump_p1`
- `pop`
- `plink` 路径
- `clump` 参考 `bfile` 前缀

如果不显式写 `--output`，脚本会优先使用 `config/defaults.env` 中的输出目录：

- `STANDARDIZED_OUTPUT_DIR`
- `OUT_OUTPUT_DIR`
- `EXP_OUTPUT_DIR`

### 6. 运行 MR 分析

将 `EXP_OUTPUT_DIR` 和 `OUT_OUTPUT_DIR` 对应目录中的文件准备好之后，运行：

```bash
bash bin/run_mr.sh EXPO HF
```

结果会写入：

```text
results/EXPO_HF/
```

包括：

- `mr_results.txt`
- `iv_table.txt`
- `heterogeneity.txt`
- `pleiotropy.txt`
- `report.txt`
- `scatter.png`
- `funnel.png`
- `leaveoneout.png`

### 7. 运行 X-M-Y 中介 MVMR

新入口会交互式逐项询问：

- `X` 对应的 exposure 文件，固定从 `data/exp/` 选择
- `M` 对应的 exposure 文件，固定从 `data/exp/` 选择
- `X / M / Y` 对应的 standardized `.gz` 文件，固定从 `data/standardized/` 选择
- `X / M / Y` 的展示名称
- 结果输出目录

然后会把选项完整列出，确认后开始分析：

```bash
bash bin/run_mvmr_mediation.sh
```

该流程当前约定：

- `X/M` 的 IV 清单来自各自的 exposure 文件，并取 SNP 并集
- `X/M/Y` 的 `beta/se/alleles/eaf` 从 standardized 文件中提取
- `X-Y` 与 `M-Y` 的方向协调使用 `TwoSampleMR::harmonise_data()`
- 在进入 `format_mvmr()` 前，会对协调后的最终 MVMR SNP list 再做一次全局 clump
- `MVMR` 部分使用 `MVMR::format_mvmr()`、`strength_mvmr()`、`pleiotropy_mvmr()`、`ivw_mvmr()`

输出目录下会生成：

- `iv_union.txt`
- `mvmr_input.txt`
- `harmonised_xy_*.txt`
- `harmonised_my_*.txt`
- `harmonised_total_xm_*.txt`
- `harmonised_total_xy_*.txt`
- `mvmr_strength.txt`
- `mvmr_pleiotropy.txt`
- `mvmr_ivw.txt`
- `mediation_summary.txt`
- `report.txt`

### 8. 运行可配置 exposure 数量的通用 MVMR

如果你希望一次纳入 2 个以上 exposure，而不是固定 `X-M-Y` 中介结构，使用：

```bash
bash bin/run_mvmr.sh
```

交互流程会先询问：

- 需要纳入的 exposure 数量
- 每个 exposure 对应的 `data/exp/` 文件
- 每个 exposure 对应的 `data/standardized/` 文件
- outcome 对应的 `data/standardized/` 文件
- 最终输出目录

该流程的规则是：

- 所有 exposure 的 IV 先取 SNP 并集
- 每个 exposure 都分别和 outcome 做一次 `harmonise_data(action = 2)`
- 以第 1 个 exposure 为锚点统一 outcome 方向
- 如果其他 exposure 的方向与锚点相反，则自动翻转对应 beta 后保留
- 在进入 `format_mvmr()` 前，会对协调后的最终 MVMR SNP list 再做一次全局 clump

输出目录下会生成：

- `iv_union.txt`
- `pairwise_harmonise_summary.txt`
- `orientation_summary.txt`
- `mvmr_input.txt`
- `mvmr_strength.txt`
- `mvmr_pleiotropy.txt`
- `mvmr_ivw.txt`
- `report.txt`

## ID 规则

标准化后的 MR 文件会同时保留三套位点标识：

- `SNP`: `bim_id`，格式为 `CHR:POS:A1:A2`，用于 MR 匹配和 clump
- `variant_id`: `CHR:POS:A1:A2`，其中 `A1/A2` 按字典序排列
- `rsid`: 注释列

标准化 TSV 额外会保留：

- `EFFECT_ALLELE` / `OTHER_ALLELE`: 效应统计量实际对应的等位基因方向
- `EAF` / `MAF`: 能解析时分别保留；缺失值显式写为 `NA`
- `EAF_1KG` / `EAF_FROM_MAF` / `EAF_ABS_DIFF`: 仅在运行 EAF 注释脚本后出现，用于记录 MAF 定向注释结果

## 标准化边界

- `Mode C (A1/A2)` 必须显式提供效应方向是 `A1` 还是 `A2`；`Unknown` 不再支持。
- 如果输入只有 `MAF` 而不是 allele-specific frequency，标准化脚本不会强行把它转换成 `EAF`；可在标准化后运行 `bin/annotate_eaf.sh` 用 1KG A1/A2 为 `MAF` 定向。
- 标准化脚本本身不做 panel-based strand flip、回文位点消歧和方向纠正；这些会留给后续注释步骤或下游 `harmonise_data()`。

## Git 管理建议

项目已经配置好 `.gitignore`，默认不跟踪：

- `data/Org/`
- `data/exp/`
- `data/out/`
- `data/standardized/`
- `results/`
- `logs/`

因此适合把代码、配置和文档放进 `git`，而把原始数据和结果文件留在工作目录本地。

## 更多说明

详细中文操作说明见 [docs/使用说明.md](docs/使用说明.md)。
