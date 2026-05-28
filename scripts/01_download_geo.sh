#!/usr/bin/env bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${PROJECT_ROOT}"

mkdir -p data/raw/GSE129308_RAW outputs/tables logs
LOG="logs/download.log"
ERROR_LOG="logs/error.log"
exec >> "${LOG}" 2>&1

trap 'echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] ERROR 01_download_geo.sh failed at line ${LINENO}: ${BASH_COMMAND}" >> "${ERROR_LOG}"' ERR

BASE_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129308/suppl"
RAW_TAR="data/raw/GSE129308_RAW.tar"
METRICS="data/raw/GSE129308_Sequencing_metrics.csv.gz"

file_size() {
  if stat -f%z "$1" >/dev/null 2>&1; then
    stat -f%z "$1"
  else
    stat -c%s "$1"
  fi
}

remote_size() {
  curl -L --fail --silent --head "$1" \
    | awk 'tolower($1)=="content-length:" {print $2}' \
    | tr -d '\r' \
    | tail -n 1
}

download_with_resume() {
  local url="$1"
  local dest="$2"
  local expected=""
  expected="$(remote_size "${url}" || true)"

  if [[ -s "${dest}" && -n "${expected}" ]]; then
    local actual
    actual="$(file_size "${dest}")"
    if [[ "${actual}" == "${expected}" ]]; then
      echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Found complete file, skipping: ${dest}"
      return 0
    fi
    echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Resuming incomplete file: ${dest} (${actual}/${expected} bytes)"
  elif [[ -s "${dest}" ]]; then
    echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Found existing file without remote size check, skipping: ${dest}"
    return 0
  else
    echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Downloading: ${url}"
  fi

  curl -L --fail -C - -o "${dest}" "${url}"

  if [[ -n "${expected}" ]]; then
    local actual_after
    actual_after="$(file_size "${dest}")"
    if [[ "${actual_after}" != "${expected}" ]]; then
      echo "Downloaded size mismatch for ${dest}: ${actual_after}/${expected}" >&2
      exit 1
    fi
  fi
}

download_with_resume "${BASE_URL}/GSE129308_RAW.tar" "${RAW_TAR}"
download_with_resume "${BASE_URL}/GSE129308_Sequencing_metrics.csv.gz" "${METRICS}"

echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Writing tar file list"
tar -tf "${RAW_TAR}" > outputs/tables/GSE129308_RAW_filelist.txt

if ! find data/raw/GSE129308_RAW -type f -name "*.h5" | grep -q .; then
  echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Extracting RAW tar"
  tar -xf "${RAW_TAR}" -C data/raw/GSE129308_RAW/
else
  echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Existing extracted H5 files found, skipping extraction"
fi

find data/raw/GSE129308_RAW -type f | sort > outputs/tables/downloaded_raw_files.txt
echo "[$(date "+%Y-%m-%d %H:%M:%S %Z")] Download step complete"
