#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config/defaults.env"
CONFIG_LOCAL_FILE="${PROJECT_DIR}/config/defaults.local.env"

if [[ -f "${CONFIG_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
fi
if [[ -f "${CONFIG_LOCAL_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CONFIG_LOCAL_FILE}"
fi

CONDA_SH="${CONDA_SH:-/home/ding/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-GWAS}"

if [[ "${SKIP_CONDA_ACTIVATE:-0}" != "1" ]]; then
    if [[ -f "${CONDA_SH}" ]]; then
        set +u
        source "${CONDA_SH}"
        conda activate "${CONDA_ENV}"
        set -u
    else
        echo "[WARN] CONDA_SH 不存在，跳过 conda activate: ${CONDA_SH}" >&2
    fi
fi

export MR_PIPELINE_R_LIB_PATH="${R_LIB_PATH:-/home/ding/R/4.4.1_MR}"
export MR_PIPELINE_EXP_DIR="${EXP_OUTPUT_DIR:-${EXPOSURE_OUTPUT_DIR:-${PROJECT_DIR}/data/exp}}"
export MR_PIPELINE_EXPOSURE_DIR="${MR_PIPELINE_EXP_DIR}"
export MR_PIPELINE_STANDARDIZED_OUTPUT_DIR="${STANDARDIZED_OUTPUT_DIR:-${PROJECT_DIR}/data/standardized}"
export MR_PIPELINE_RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"
export MR_PIPELINE_MVMR_CLUMP_PLINK="${MVMR_GLOBAL_CLUMP_PLINK:-${PLINK_BIN:-/home/ding/miniconda3/envs/GWAS/bin/plink}}"
export MR_PIPELINE_MVMR_CLUMP_BFILE="${MVMR_GLOBAL_CLUMP_BFILE:-${CLUMP_BFILE:-/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur_colon}}"
export MR_PIPELINE_MVMR_CLUMP_R2="${MVMR_GLOBAL_CLUMP_R2:-${CLUMP_R2:-0.01}}"
export MR_PIPELINE_MVMR_CLUMP_KB="${MVMR_GLOBAL_CLUMP_KB:-${CLUMP_KB:-1000}}"
export MR_PIPELINE_MVMR_CLUMP_P1="${MVMR_GLOBAL_CLUMP_P1:-1}"

mkdir -p "${MR_PIPELINE_RESULTS_DIR}"

python "${PROJECT_DIR}/scripts/run_mvmr_mediation.py" "$@"
