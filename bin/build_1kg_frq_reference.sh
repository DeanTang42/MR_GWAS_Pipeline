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

BFILE="${CLUMP_BFILE:-/home/ding/MR_LPA/Ref/g1000_eur/g1000_eur_colon}"
PLINK="${PLINK_BIN:-/home/ding/miniconda3/envs/GWAS/bin/plink}"
REFERENCE_DIR="${REFERENCE_DIR:-${PROJECT_DIR}/data/reference}"
OUT_PREFIX="${REFERENCE_DIR}/$(basename "${BFILE}")"

usage() {
    cat <<'EOF'
用法:
  bash bin/build_1kg_frq_reference.sh [--bfile PREFIX] [--out-prefix PREFIX] [--plink PATH]

输出:
  OUT_PREFIX.frq         PLINK 原生频率文件，包含 SNP/A1/A2/MAF 等列
  OUT_PREFIX.frq.tsv.gz  规范化 TSV，仅保留 SNP/A1/A2/MAF，供 EAF 注释脚本读取
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --bfile)
            BFILE="$2"
            shift 2
            ;;
        --out-prefix)
            OUT_PREFIX="$2"
            shift 2
            ;;
        --plink)
            PLINK="$2"
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "[ERROR] 未知参数: $1" >&2
            usage >&2
            exit 1
            ;;
    esac
done

mkdir -p "$(dirname "${OUT_PREFIX}")"

echo "[INFO] PLINK: ${PLINK}"
echo "[INFO] bfile: ${BFILE}"
echo "[INFO] output prefix: ${OUT_PREFIX}"

"${PLINK}" --bfile "${BFILE}" --freq --out "${OUT_PREFIX}"

FRQ_FILE="${OUT_PREFIX}.frq"
FRQ_TSV_GZ="${OUT_PREFIX}.frq.tsv.gz"

awk '
    BEGIN { OFS = "\t" }
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            idx[$i] = i
        }
        if (!("SNP" in idx) || !("A1" in idx) || !("A2" in idx) || !("MAF" in idx)) {
            print "[ERROR] .frq 缺少 SNP/A1/A2/MAF 列" > "/dev/stderr"
            exit 1
        }
        print "SNP", "A1", "A2", "MAF"
        next
    }
    NF > 0 {
        print $idx["SNP"], $idx["A1"], $idx["A2"], $idx["MAF"]
    }
' "${FRQ_FILE}" | gzip -c > "${FRQ_TSV_GZ}"

gzip -t "${FRQ_TSV_GZ}"

echo "[INFO] PLINK .frq: ${FRQ_FILE}"
echo "[INFO] normalized reference: ${FRQ_TSV_GZ}"
