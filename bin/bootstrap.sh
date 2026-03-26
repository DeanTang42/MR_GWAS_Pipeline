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

python -m pip install -r "${PROJECT_DIR}/requirements.txt"
echo "Python dependencies installed."
echo "R package list: ${PROJECT_DIR}/config/r_packages.txt"
