# MR_GWAS_Pipeline

这是一个可独立运行、可用 `git` 管理的 MR 分析项目骨架，包含两部分核心能力：

1. GWAS 摘要统计标准化与等位基因对齐。
2. 基于 `TwoSampleMR` 的单暴露单结局 MR 分析。

项目已经按当前机器环境整理好默认配置，代码、配置、文档可以直接纳入 `git`；大体积参考数据、临时文件和结果文件默认被 `.gitignore` 排除。

## 目录结构

```text
MR_GWAS_Pipeline/
├── bin/
├── config/
├── data/
│   ├── standardized/
│   ├── exposure/
│   ├── outcome/
│   └── reference/
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

- `REFERENCE_PANEL`
- `STANDARDIZED_OUTPUT_DIR`
- `EXPOSURE_OUTPUT_DIR`
- `OUTCOME_OUTPUT_DIR`
- `RESULTS_DIR`
- `PLINK_BIN`
- `CLUMP_BFILE`
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

### 2. 准备 reference panel

项目不再提供 reference panel 构建脚本。你只需要把自己构建好的 panel 放到 `REFERENCE_PANEL` 指定的位置即可。格式要求见 [docs/reference_panel格式说明.md](docs/reference_panel格式说明.md)。

### 3. 标准化 outcome GWAS

如果不加 `--non-interactive`，脚本会进入交互模式。现在对于 exposure 的 clump 参数也会在交互过程中询问，并自动带出 `config/defaults.env` 里的默认值。

```bash
bash bin/standardize_gwas.sh \
  --non-interactive \
  --input /path/to/outcome.tsv.gz \
  --output /home/ding/MR_GWAS_Pipeline/data/outcome/HF.csv \
  --output-format mr \
  --mr-role out \
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
  --freq-col effect_allele_frequency \
  --phenotype HF \
  --sample-size 344182
```

### 4. 标准化 exposure GWAS 并自动 clump

`--mr-role exp` 会自动执行 `clump_data()`，并使用 `bim_id` 作为 `SNP` 与参考 `.bim` 对接。非 `--non-interactive` 模式下会交互式询问：

- `clump_r2`
- `clump_kb`
- `clump_p1`
- `pop`
- `plink` 路径
- `clump` 参考 `bfile` 前缀

```bash
bash bin/standardize_gwas.sh \
  --non-interactive \
  --input /path/to/exposure.tsv.gz \
  --output-format mr \
  --mr-role exp \
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
  --freq-col effect_allele_frequency \
  --phenotype EXPO \
  --sample-size 344182
```

如果不显式写 `--output`，脚本会优先使用 `config/defaults.env` 中的输出目录：

- `STANDARDIZED_OUTPUT_DIR`
- `OUTCOME_OUTPUT_DIR`
- `EXPOSURE_OUTPUT_DIR`

### 5. 运行 MR 分析

将 `EXPOSURE_OUTPUT_DIR` 和 `OUTCOME_OUTPUT_DIR` 对应目录中的文件准备好之后，运行：

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

## ID 规则

标准化后的 MR 文件会同时保留三套位点标识：

- `SNP`: `bim_id`，格式为 `CHR_POS_A1_A2`，用于 MR 匹配和 clump
- `variant_id`: `CHR:POS:REF:ALT`，用于可读性和位点追踪
- `rsid`: 注释列

## reference panel 的作用

reference panel 不是拿来直接做 MR 的，而是用来：

1. 统一位点编码
2. 统一等位基因方向
3. 用频率消解回文 SNP
4. 剔除方向不可靠或无法匹配的位点

## Git 管理建议

项目已经配置好 `.gitignore`，默认不跟踪：

- `data/exposure/`
- `data/outcome/`
- `data/standardized/`
- `data/reference/`
- `results/`
- `logs/`

因此适合把代码、配置和文档放进 `git`，而把大 reference 和结果文件留在工作目录本地。

## 更多说明

详细中文操作说明见 [docs/使用说明.md](docs/使用说明.md)，reference panel 格式说明见 [docs/reference_panel格式说明.md](docs/reference_panel格式说明.md)。
