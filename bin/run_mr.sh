#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_FILE="${PROJECT_DIR}/config/defaults.env"

if [[ -f "${CONFIG_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${CONFIG_FILE}"
fi

R_LIB_PATH="${R_LIB_PATH:-/home/ding/R/4.4.1_MR}"
export MR_PIPELINE_R_LIB_PATH="${R_LIB_PATH}"
export MR_PIPELINE_ORG_DIR="${ORG_DATA_DIR:-${PROJECT_DIR}/data/Org}"
export MR_PIPELINE_EXP_DIR="${EXP_OUTPUT_DIR:-${EXPOSURE_OUTPUT_DIR:-${PROJECT_DIR}/data/exp}}"
export MR_PIPELINE_OUT_DIR="${OUT_OUTPUT_DIR:-${OUTCOME_OUTPUT_DIR:-${PROJECT_DIR}/data/out}}"
export MR_PIPELINE_EXPOSURE_DIR="${MR_PIPELINE_EXP_DIR}"
export MR_PIPELINE_OUTCOME_DIR="${MR_PIPELINE_OUT_DIR}"
export MR_PIPELINE_RESULTS_DIR="${RESULTS_DIR:-${PROJECT_DIR}/results}"

mkdir -p "${MR_PIPELINE_RESULTS_DIR}"

Rscript "${PROJECT_DIR}/scripts/MR_single_pipeline.R" \
    --project-dir "${PROJECT_DIR}" \
    --r-lib-path "${R_LIB_PATH}" \
    "$@"
