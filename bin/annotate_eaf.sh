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

export MR_PIPELINE_STANDARDIZED_OUTPUT_DIR="${STANDARDIZED_OUTPUT_DIR:-${PROJECT_DIR}/data/standardized}"
export MR_PIPELINE_REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/data/reference}"
export MR_PIPELINE_FRQ_REFERENCE="${FRQ_REFERENCE:-${MR_PIPELINE_REFERENCE_DIR}/$(basename "${CLUMP_BFILE:-g1000_eur_colon}").frq.tsv.gz}"

mkdir -p "${MR_PIPELINE_REFERENCE_DIR}"

python "${PROJECT_DIR}/scripts/annotate_eaf_from_frq.py" "$@"
