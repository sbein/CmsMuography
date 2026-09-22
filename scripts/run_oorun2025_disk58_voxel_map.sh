#!/usr/bin/env bash
set -euo pipefail

cmssw_src="${CMSSW_SRC:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src}"
project_root="${PROJECT_ROOT:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern}"
input_list="${INPUT_LIST:-CosmicMuography/data/OORun2025_AOD_disk_root_one_pfn_per_file_xrdfs_readable.txt}"
output_root="${OUTPUT_ROOT:-CosmicMuography/test_output/OORun2025_AOD_disk58_voxel_map_zup.root}"
pedestal_root="${PEDESTAL_ROOT:-CosmicMuography/test_output/OORun2025_AOD_disk58_voxel_map_zup_pedestals.root}"
progress_every="${PROGRESS_EVERY:-10000}"

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${cmssw_src}"
eval "$(scram runtime -sh)"

cd "${project_root}"
mkdir -p "$(dirname "${output_root}")"

echo "start: $(date -Is)"
echo "host: $(hostname -f 2>/dev/null || hostname)"
echo "input_list: ${input_list}"
echo "output_root: ${output_root}"
echo "pedestal_root: ${pedestal_root}"
echo "n_input_pfns: $(grep -vc '^[[:space:]]*#' "${input_list}")"

PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 python3 CosmicMuography/scripts/cosmic_voxel_map.py \
  --input "${input_list}" \
  --output "${output_root}" \
  --skip-bad-files \
  --progress-every "${progress_every}"

PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 python3 CosmicMuography/scripts/suppress_voxel_pedestals.py \
  --input "${output_root}" \
  --output "${pedestal_root}" \
  --write-source

echo "finish: $(date -Is)"
