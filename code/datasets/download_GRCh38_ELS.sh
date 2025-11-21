#!/bin/bash
# Download ENCODE SCREEN candidate enhancers (GRCh38)

set -euo pipefail

OUTPUT_DIR="data/genome"
URL="https://downloads.wenglab.org/cCREs/GRCh38-ELS.bed"
OUTPUT_FILE="${OUTPUT_DIR}/GRCh38-ELS.bed"

mkdir -p "${OUTPUT_DIR}"

if [ -f "${OUTPUT_FILE}" ]; then
  echo "GRCh38-ELS.bed already exists, skipping download"
  exit 0
fi

echo "Downloading GRCh38-ELS.bed from ${URL}..."
curl -L "${URL}" -o "${OUTPUT_FILE}"

if [ $? -eq 0 ]; then
  echo "Successfully downloaded GRCh38-ELS.bed"
  wc -l "${OUTPUT_FILE}" || true
else
  echo "Failed to download GRCh38-ELS.bed" >&2
  exit 1
fi


