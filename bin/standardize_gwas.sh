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

set +u
source "${CONDA_SH}"
conda activate "${CONDA_ENV}"
set -u

export MR_PIPELINE_R_LIB_PATH="${R_LIB_PATH:-/home/ding/R/4.4.1_MR}"
export MR_PIPELINE_REFERENCE_PANEL="${REFERENCE_PANEL:-}"
export MR_PIPELINE_STANDARDIZED_OUTPUT_DIR="${STANDARDIZED_OUTPUT_DIR:-${PROJECT_DIR}/data/standardized}"
export MR_PIPELINE_EXPOSURE_DIR="${EXPOSURE_OUTPUT_DIR:-${PROJECT_DIR}/data/exposure}"
export MR_PIPELINE_OUTCOME_DIR="${OUTCOME_OUTPUT_DIR:-${PROJECT_DIR}/data/outcome}"
export MR_PIPELINE_CLUMP_PLINK="${PLINK_BIN:-/home/ding/miniconda3/envs/GWAS/bin/plink}"
export MR_PIPELINE_CLUMP_BFILE="${CLUMP_BFILE:-/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur}"
export MR_PIPELINE_CLUMP_R2="${CLUMP_R2:-0.1}"
export MR_PIPELINE_CLUMP_KB="${CLUMP_KB:-500}"
export MR_PIPELINE_CLUMP_P1="${CLUMP_P1:-1e-4}"
export MR_PIPELINE_CLUMP_POP="${CLUMP_POP:-EUR}"

mkdir -p "${MR_PIPELINE_STANDARDIZED_OUTPUT_DIR}" "${MR_PIPELINE_EXPOSURE_DIR}" "${MR_PIPELINE_OUTCOME_DIR}"

args=("$@")
has_reference=false
for ((i = 0; i < ${#args[@]}; i++)); do
    if [[ "${args[i]}" == "--reference" ]]; then
        has_reference=true
        break
    fi
done

if [[ "${has_reference}" == false && -n "${REFERENCE_PANEL:-}" && -f "${REFERENCE_PANEL}" ]]; then
    args=(--reference "${REFERENCE_PANEL}" "${args[@]}")
fi

python "${PROJECT_DIR}/scripts/gwas_standardizer.py" "${args[@]}"
