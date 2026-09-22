#!/usr/bin/env bash
set -euo pipefail

cmssw_src="${CMSSW_SRC:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src}"
project_root="${PROJECT_ROOT:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern}"
input_list="${INPUT_LIST:-CosmicMuography/data/OORun2025_AOD_disk_root_one_pfn_per_file_xrdfs_readable.txt}"
run_id="${RUN_ID:-$(date +%Y%m%d_%H%M%S)}"
parts_dir="${PARTS_DIR:-CosmicMuography/test_output/OORun2025_AOD_disk58_parts_${run_id}}"
part_logs_dir="${PART_LOGS_DIR:-CosmicMuography/logs/OORun2025_AOD_disk58_parts_${run_id}}"
merged_root="${MERGED_ROOT:-CosmicMuography/test_output/OORun2025_AOD_disk58_voxel_map_zup_robust.root}"
pedestal_root="${PEDESTAL_ROOT:-CosmicMuography/test_output/OORun2025_AOD_disk58_voxel_map_zup_robust_pedestals.root}"
per_file_timeout="${PER_FILE_TIMEOUT:-1200}"
progress_every="${PROGRESS_EVERY:-10000}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${cmssw_src}"
eval "$(scram runtime -sh)"

cd "${project_root}"
mkdir -p "${parts_dir}" "${part_logs_dir}" "$(dirname "${merged_root}")"

success_list="${parts_dir}/successful_parts.txt"
failure_list="${parts_dir}/failed_pfns.tsv"
: > "${success_list}"
printf "index\tstatus\tpfn\tpart_log\n" > "${failure_list}"

echo "start: $(date -Is)"
echo "host: $(hostname -f 2>/dev/null || hostname)"
echo "input_list: ${input_list}"
echo "parts_dir: ${parts_dir}"
echo "part_logs_dir: ${part_logs_dir}"
echo "merged_root: ${merged_root}"
echo "pedestal_root: ${pedestal_root}"
echo "per_file_timeout: ${per_file_timeout}"

index=0
while IFS= read -r pfn; do
  [[ -z "${pfn}" || "${pfn}" =~ ^[[:space:]]*# ]] && continue
  index=$((index + 1))
  part_root="$(printf "%s/part_%03d.root" "${parts_dir}" "${index}")"
  part_log="$(printf "%s/part_%03d.log" "${part_logs_dir}" "${index}")"

  echo "file ${index}: start $(date -Is)"
  echo "file ${index}: ${pfn}"
  set +e
  timeout --kill-after=60s "${per_file_timeout}" \
    env PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 \
    python3 CosmicMuography/scripts/cosmic_voxel_map.py \
      --input "${pfn}" \
      --output "${part_root}" \
      --progress-every "${progress_every}" \
      > "${part_log}" 2>&1
  status=$?
  set -e

  if [[ "${status}" -eq 0 ]] && grep -q '^Wrote ' "${part_log}" && [[ -s "${part_root}" ]]; then
    echo "${part_root}" >> "${success_list}"
    echo "file ${index}: success $(date -Is)"
  else
    printf "%s\t%s\t%s\t%s\n" "${index}" "${status}" "${pfn}" "${part_log}" >> "${failure_list}"
    echo "file ${index}: failed status=${status} $(date -Is)"
    tail -20 "${part_log}" || true
  fi
done < "${input_list}"

n_success=$(wc -l < "${success_list}")
n_failure=$(( $(wc -l < "${failure_list}") - 1 ))
echo "per-file loop done: success=${n_success} failure=${n_failure}"

if [[ "${n_success}" -le 0 ]]; then
  echo "No successful per-file outputs; not merging." >&2
  exit 1
fi

PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 python3 CosmicMuography/scripts/merge_voxel_maps.py \
  --input-list "${success_list}" \
  --output "${merged_root}" \
  --skip-bad-files

PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 python3 CosmicMuography/scripts/suppress_voxel_pedestals.py \
  --input "${merged_root}" \
  --output "${pedestal_root}" \
  --write-source

echo "finish: $(date -Is)"
