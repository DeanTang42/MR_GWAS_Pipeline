#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config/defaults.env"

if [[ -f "${CONFIG_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
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

mkdir -p "${MR_PIPELINE_RESULTS_DIR}"

python "${PROJECT_DIR}/scripts/run_mvmr.py" "$@"
