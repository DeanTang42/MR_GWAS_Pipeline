# syntax=docker/dockerfile:1
FROM mambaorg/micromamba:1.5.10

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends ca-certificates bash git procps \
    && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER
WORKDIR /opt/MR_GWAS_Pipeline

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba create -y -f /tmp/environment.yml \
    && micromamba clean --all --yes

SHELL ["micromamba", "run", "-n", "mr", "/bin/bash", "-lc"]
RUN Rscript -e 'options(repos=c(MRCIEU="https://mrcieu.r-universe.dev", CRAN="https://cloud.r-project.org")); install.packages(c("TwoSampleMR","mr.raps"), Ncpus=max(1, parallel::detectCores()-1))' \
    && Rscript -e 'remotes::install_github("rondolab/MR-PRESSO", upgrade="never")' \
    && Rscript -e 'library(data.table); library(TwoSampleMR); library(MRPRESSO); library(mr.raps); library(dplyr); library(ggplot2); cat("R packages OK\n")'

SHELL ["/bin/bash", "-lc"]
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/MR_GWAS_Pipeline

USER root
RUN install -m 0755 /opt/MR_GWAS_Pipeline/docker/entrypoint.sh /usr/local/bin/mr-pipeline-entrypoint \
    && chmod +x /opt/MR_GWAS_Pipeline/bin/*.sh

USER $MAMBA_USER
ENV PROJECT_DIR=/opt/MR_GWAS_Pipeline \
    CONDA_ENV=mr \
    CONDA_PREFIX=/opt/conda/envs/mr \
    R_LIB_PATH=/opt/conda/envs/mr/lib/R/library \
    PLINK_BIN=/opt/conda/envs/mr/bin/plink \
    CLUMP_BFILE=/ref/g1000_eur_colon \
    ORG_DATA_DIR=/opt/MR_GWAS_Pipeline/data/Org \
    STANDARDIZED_OUTPUT_DIR=/opt/MR_GWAS_Pipeline/data/standardized \
    EXP_OUTPUT_DIR=/opt/MR_GWAS_Pipeline/data/exp \
    OUT_OUTPUT_DIR=/opt/MR_GWAS_Pipeline/data/out \
    RESULTS_DIR=/opt/MR_GWAS_Pipeline/results \
    SKIP_CONDA_ACTIVATE=1

ENTRYPOINT ["mr-pipeline-entrypoint"]
CMD ["bash"]
