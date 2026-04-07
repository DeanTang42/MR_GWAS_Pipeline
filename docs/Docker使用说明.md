# Docker 使用说明

这个项目可以打包成 Docker 镜像，镜像内固定 Python、R、PLINK 和 R 包依赖；原始 GWAS、处理结果、reference bfile 不打进镜像，运行时通过 volume 挂载。

## 构建镜像

在项目根目录执行：

```bash
cd /path/to/MR_GWAS_Pipeline
docker build -t mr-gwas-pipeline:latest .
```

构建时会从 conda-forge、bioconda、MRCIEU r-universe 和 GitHub 下载依赖，其中 `MRPRESSO` 来自 GitHub。因此构建环境需要能访问这些源。

## 推荐目录挂载

假设服务器上的项目目录是：

```text
/hwmaster/tangyuzhong/project/MR_GWAS_Pipeline
```

假设 PLINK reference 目录是：

```text
/path/to/g1000_eur/
```

并且该目录里包含：

```text
g1000_eur_colon.bed
g1000_eur_colon.bim
g1000_eur_colon.fam
```

进入容器：

```bash
docker run --rm -it \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/data:/opt/MR_GWAS_Pipeline/data \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/results:/opt/MR_GWAS_Pipeline/results \
  -v /path/to/g1000_eur:/ref \
  -e CLUMP_BFILE=/ref/g1000_eur_colon \
  mr-gwas-pipeline:latest bash
```

容器内默认路径：

```text
PROJECT_DIR=/opt/MR_GWAS_Pipeline
PLINK_BIN=/opt/conda/envs/mr/bin/plink
R_LIB_PATH=/opt/conda/envs/mr/lib/R/library
CLUMP_BFILE=/ref/g1000_eur_colon
```

## 标准化 exposure 并 clump

可以直接在宿主机执行：

```bash
docker run --rm -it \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/data:/opt/MR_GWAS_Pipeline/data \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/results:/opt/MR_GWAS_Pipeline/results \
  -v /path/to/g1000_eur:/ref \
  -e CLUMP_BFILE=/ref/g1000_eur_colon \
  mr-gwas-pipeline:latest \
  bash bin/standardize_gwas.sh \
    --non-interactive \
    --input /opt/MR_GWAS_Pipeline/data/Org/input.tsv.gz \
    --output /opt/MR_GWAS_Pipeline/data/exp/EXPOSURE.csv \
    --output-format mr \
    --mr-role exp \
    --mode B \
    --snp-col rsids \
    --chr-col Chrom \
    --pos-col Pos_hg19 \
    --allele1-col effectAllele \
    --allele2-col otherAllele \
    --stat-type BETA \
    --stat-col Beta \
    --se-col SE \
    --p-col Pval \
    --pval-format raw \
    --freq-type EAF \
    --freq-col ImpMAF \
    --phenotype ANGPTL3 \
    --sample-size 35682 \
    --clump-p1 5e-8
```

说明：如果输入只有 `ImpMAF`，但你希望最终 `eaf.exposure` 不为空，可以像上面这样用 `--freq-type EAF --freq-col ImpMAF`，即把 `ImpMAF` 作为 effect allele frequency 使用。

## 运行 MR

准备好：

```text
data/exp/EXPOSURE.csv
data/out/OUTCOME.csv
```

然后执行：

```bash
docker run --rm -it \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/data:/opt/MR_GWAS_Pipeline/data \
  -v /hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/results:/opt/MR_GWAS_Pipeline/results \
  mr-gwas-pipeline:latest \
  bash bin/run_mr.sh EXPOSURE OUTCOME
```

结果会写入宿主机：

```text
/hwmaster/tangyuzhong/project/MR_GWAS_Pipeline/results/EXPOSURE_OUTCOME/
```

## 常见问题

- 如果 clump 结果为 0，先检查 `SNP` 格式是否和 `CLUMP_BFILE.bim` 第二列一致。本项目默认使用 `CHR:POS:A1:A2`，因此推荐挂载 `g1000_eur_colon` 版本。
- 如果提示找不到 R 包，说明你没有使用 Docker 镜像内的 `Rscript`，或者镜像构建时 R 包安装失败。先在容器里运行 `Rscript -e 'library(data.table); library(TwoSampleMR); library(MRPRESSO); library(mr.raps)'` 检查。
- 不要把大型 GWAS、reference 和结果文件复制进镜像；这些文件应通过 `-v` 挂载。
