#!/usr/bin/env bash
set -euo pipefail

if command -v micromamba >/dev/null 2>&1; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "${CONDA_ENV:-mr}"
fi

export PROJECT_DIR="${PROJECT_DIR:-/opt/MR_GWAS_Pipeline}"
export R_LIB_PATH="${R_LIB_PATH:-${CONDA_PREFIX:-/opt/conda/envs/mr}/lib/R/library}"
export PLINK_BIN="${PLINK_BIN:-${CONDA_PREFIX:-/opt/conda/envs/mr}/bin/plink}"
export CLUMP_BFILE="${CLUMP_BFILE:-/ref/g1000_eur_colon}"
export ORG_DATA_DIR="${ORG_DATA_DIR:-${PROJECT_DIR}/data/Org}"
export STANDARDIZED_OUTPUT_DIR="${STANDARDIZED_OUTPUT_DIR:-${PROJECT_DIR}/data/standardized}"
export EXP_OUTPUT_DIR="${EXP_OUTPUT_DIR:-${PROJECT_DIR}/data/exp}"
export OUT_OUTPUT_DIR="${OUT_OUTPUT_DIR:-${PROJECT_DIR}/data/out}"
export RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
export SKIP_CONDA_ACTIVATE="${SKIP_CONDA_ACTIVATE:-1}"

exec "$@"
