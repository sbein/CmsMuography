#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
EXPACS_DIR="${EXPACS_DIR:-${REPO_DIR}/../External/EXPACS}"
CXX="${CXX:-g++}"
OUTPUT="${OUTPUT:-${REPO_DIR}/expacs/expacs_muon_grid}"

if [[ ! -f "${EXPACS_DIR}/subroutines.cpp" ]]; then
  echo "Missing ${EXPACS_DIR}/subroutines.cpp; set EXPACS_DIR to unpacked parma_cpp." >&2
  exit 2
fi

"${CXX}" -O2 -std=c++17 \
  "${REPO_DIR}/expacs/expacs_muon_grid.cpp" \
  "${EXPACS_DIR}/subroutines.cpp" \
  -o "${OUTPUT}"

echo "built ${OUTPUT}"
