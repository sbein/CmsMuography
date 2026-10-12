#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
WORKSPACE_DIR="$(cd "${REPO_DIR}/.." && pwd)"
CMSSW_SRC="${CMSSW_SRC:-${WORKSPACE_DIR}/CMSSW_16_1_0_pre4/src}"
INPUT_LIST="${INPUT_LIST:-${REPO_DIR}/data/OORun2025_AOD_100.txt}"
RUN_ID="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
PARTS_DIR="${PARTS_DIR:-${REPO_DIR}/test_output/voxel_parts_${RUN_ID}}"
PART_LOGS_DIR="${PART_LOGS_DIR:-${REPO_DIR}/logs/voxel_parts_${RUN_ID}}"
MERGED_ROOT="${MERGED_ROOT:-${REPO_DIR}/test_output/OORun2025_AOD_voxel_map_zup_robust.root}"
PEDESTAL_ROOT="${PEDESTAL_ROOT:-${REPO_DIR}/test_output/OORun2025_AOD_voxel_map_zup_robust_pedestals.root}"
PER_FILE_TIMEOUT="${PER_FILE_TIMEOUT:-1200}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scram runtime -sh)"
cd "${REPO_DIR}"
mkdir -p "${PARTS_DIR}" "${PART_LOGS_DIR}" "$(dirname "${MERGED_ROOT}")"

success_list="${PARTS_DIR}/successful_parts.txt"
failure_list="${PARTS_DIR}/failed_inputs.tsv"
: > "${success_list}"
printf "index\tstatus\tinput\tpart_log\n" > "${failure_list}"

index=0
while IFS= read -r input; do
  [[ -z "${input}" || "${input}" =~ ^[[:space:]]*# ]] && continue
  index=$((index + 1))
  part_root="$(printf "%s/part_%03d.root" "${PARTS_DIR}" "${index}")"
  part_log="$(printf "%s/part_%03d.log" "${PART_LOGS_DIR}" "${index}")"

  echo "file ${index}: ${input}"
  set +e
  timeout --kill-after=60s "${PER_FILE_TIMEOUT}" \
    env PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 \
    python3 scripts/cosmic_voxel_map.py \
      --input "${input}" \
      --output "${part_root}" \
      > "${part_log}" 2>&1
  status=$?
  set -e

  if [[ "${status}" -eq 0 ]] && grep -q '^Wrote ' "${part_log}" && [[ -s "${part_root}" ]]; then
    echo "${part_root}" >> "${success_list}"
  else
    printf "%s\t%s\t%s\t%s\n" "${index}" "${status}" "${input}" "${part_log}" >> "${failure_list}"
    tail -20 "${part_log}" || true
  fi
done < "${INPUT_LIST}"

if [[ ! -s "${success_list}" ]]; then
  echo "No successful per-file outputs; not merging." >&2
  exit 1
fi

python3 scripts/merge_voxel_maps.py --input-list "${success_list}" --output "${MERGED_ROOT}" --skip-bad-files
python3 scripts/suppress_voxel_pedestals.py --input "${MERGED_ROOT}" --output "${PEDESTAL_ROOT}" --write-source

echo "wrote ${MERGED_ROOT}"
echo "wrote ${PEDESTAL_ROOT}"
