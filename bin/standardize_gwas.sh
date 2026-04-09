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
export MR_PIPELINE_ORG_DIR="${ORG_DATA_DIR:-${PROJECT_DIR}/data/Org}"
export MR_PIPELINE_STANDARDIZED_OUTPUT_DIR="${STANDARDIZED_OUTPUT_DIR:-${PROJECT_DIR}/data/standardized}"
export MR_PIPELINE_EXP_DIR="${EXP_OUTPUT_DIR:-${EXPOSURE_OUTPUT_DIR:-${PROJECT_DIR}/data/exp}}"
export MR_PIPELINE_OUT_DIR="${OUT_OUTPUT_DIR:-${OUTCOME_OUTPUT_DIR:-${PROJECT_DIR}/data/out}}"
export MR_PIPELINE_EXPOSURE_DIR="${MR_PIPELINE_EXP_DIR}"
export MR_PIPELINE_OUTCOME_DIR="${MR_PIPELINE_OUT_DIR}"
export MR_PIPELINE_CLUMP_PLINK="${PLINK_BIN:-/home/ding/miniconda3/envs/GWAS/bin/plink}"
export MR_PIPELINE_CLUMP_BFILE="${CLUMP_BFILE:-/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur_colon}"
export MR_PIPELINE_CLUMP_R2="${CLUMP_R2:-0.1}"
export MR_PIPELINE_CLUMP_KB="${CLUMP_KB:-500}"
export MR_PIPELINE_CLUMP_P1="${CLUMP_P1:-1e-4}"
export MR_PIPELINE_CLUMP_POP="${CLUMP_POP:-EUR}"

mkdir -p "${MR_PIPELINE_ORG_DIR}" "${MR_PIPELINE_STANDARDIZED_OUTPUT_DIR}" "${MR_PIPELINE_EXP_DIR}" "${MR_PIPELINE_OUT_DIR}"

python "${PROJECT_DIR}/scripts/gwas_standardizer.py" "$@"
