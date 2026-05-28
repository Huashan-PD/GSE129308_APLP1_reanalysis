#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

mkdir -p data/raw/SRA data/raw/FASTQ outputs/tables logs
LOG="logs/download_sra_optional.log"
ERROR_LOG="logs/error.log"
exec >> "${LOG}" 2>&1

trap 'echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] ERROR 01b_download_sra_optional.sh failed at line ${LINENO}: ${BASH_COMMAND}" >> "${ERROR_LOG}"' ERR

if [[ "${DOWNLOAD_FASTQ:-0}" != "1" ]]; then
  echo "DOWNLOAD_FASTQ is not 1; skipping optional SRA/FASTQ download."
  exit 0
fi

for cmd in prefetch fasterq-dump curl; do
  if ! command -v "${cmd}" >/dev/null 2>&1; then
    echo "Required command not found for optional FASTQ download: ${cmd}" >&2
    exit 1
  fi
done

RUNINFO="outputs/tables/GSE129308_SRA_RunInfo.csv"
RUNS="outputs/tables/GSE129308_SRR_accessions.txt"

echo "Fetching SRA RunInfo for GSE129308"
curl -L --fail "https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=GSE129308" -o "${RUNINFO}"

awk -F',' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="Run") run_col=i; next} NR>1 && run_col>0 && $run_col ~ /^SRR/ {print $run_col}' "${RUNINFO}" \
  | sort -u > "${RUNS}"

if [[ ! -s "${RUNS}" ]]; then
  echo "No SRR accessions parsed from ${RUNINFO}" >&2
  exit 1
fi

while read -r srr; do
  [[ -z "${srr}" ]] && continue
  echo "Downloading ${srr}"
  prefetch "${srr}" --output-directory data/raw/SRA
  fasterq-dump "data/raw/SRA/${srr}" --outdir data/raw/FASTQ --split-files --threads "${FASTQ_THREADS:-4}"
done < "${RUNS}"

echo "Optional FASTQ download complete. Cell Ranger reprocessing is intentionally not run by default."
