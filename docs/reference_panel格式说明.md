# 1KG `.frq` 参考频率说明

当前项目的标准化流程本身不再依赖 reference panel。reference 现在只用于一个独立步骤：当 GWAS 只有 `MAF`、没有 `EAF` 时，用 1KG `.frq` 的 A1/A2 方向为原文件 `MAF` 定向，从而补充 `EAF`。

## 1. 构建命令

默认使用 [defaults.env](/home/ding/MR_GWAS_Pipeline/config/defaults.env) 里的 `CLUMP_BFILE`：

```bash
bash bin/build_1kg_frq_reference.sh
```

也可以手动指定：

```bash
bash bin/build_1kg_frq_reference.sh \
  --bfile /path/to/g1000_eur_colon \
  --out-prefix /home/ding/MR_GWAS_Pipeline/data/reference/g1000_eur_colon
```

## 2. 输出文件

构建脚本会生成两个文件：

- `data/reference/g1000_eur_colon.frq`
- `data/reference/g1000_eur_colon.frq.tsv.gz`

其中 `.frq` 是 PLINK 原生输出，`.frq.tsv.gz` 是给 Python 注释脚本读取的规范化版本。

## 3. 必需列

规范化版本必须包含以下列：

```text
SNP	A1	A2	MAF
```

语义如下：

- `SNP`: 必须和标准化 GWAS 的 `BIM_ID` 一致，例如 `CHR:POS:A1:A2`
- `A1`: PLINK `.frq` 里的 A1
- `A2`: PLINK `.frq` 里的 A2
- `MAF`: A1 的频率

## 4. EAF 注释逻辑

注释脚本读取标准化 GWAS：

```text
BIM_ID	CHR	POS	VARIANT_ID	RSID	EFFECT_ALLELE	OTHER_ALLELE	BETA	SE	P	EAF	MAF
```

只处理 `EAF=NA` 且 `MAF` 有值的行。

如果 `BIM_ID` 能匹配 `.frq`，且等位基因方向一致：

- `EFFECT_ALLELE == A1`: `EAF_FROM_MAF = MAF`
- `EFFECT_ALLELE == A2`: `EAF_FROM_MAF = 1 - MAF`
- `EAF_1KG` 同理使用 1KG `.frq` 中的 `MAF` 计算
- 主列 `EAF` 用 `EAF_FROM_MAF` 填充

如果匹配不到 reference、`MAF` 缺失或等位基因对不上，`EAF` 继续保持 `NA`。

## 5. 审计列

运行 `annotate_eaf.sh` 后，标准化文件会新增：

- `EAF_1KG`
- `EAF_FROM_MAF`
- `EAF_ABS_DIFF`

## 6. 配套输出

默认命令：

```bash
bash bin/annotate_eaf.sh \
  --non-interactive \
  --input data/standardized/GWAS.tsv.gz
```

会替换原标准化文件，并生成：

- `GWAS.eaf_matched.tsv.gz`
- `GWAS.eaf_frequency_diff.tsv.gz`
- `GWAS.eaf_annotation.log`

如果不希望替换原文件：

```bash
bash bin/annotate_eaf.sh \
  --non-interactive \
  --input data/standardized/GWAS.tsv.gz \
  --no-replace-input \
  --output data/standardized/GWAS.eaf_annotated.tsv.gz
```
