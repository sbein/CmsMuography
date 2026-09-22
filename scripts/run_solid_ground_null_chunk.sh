#!/usr/bin/env bash
set -euo pipefail

CHUNK_INDEX="${1:?chunk index}"
NCHUNKS="${2:?number of chunks}"

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern}"
CMSSW_SRC="${CMSSW_SRC:-${PROJECT_DIR}/CMSSW_16_1_0_pre4/src}"
DATA_ROOT="${DATA_ROOT:-${PROJECT_DIR}/CosmicMuography/test_output/OORun2025_AOD_disk55_no_kipt_voxel_map_zup_120m.root}"
TRACKS_ROOT="${TRACKS_ROOT:-${PROJECT_DIR}/CosmicMuography/test_output/OORun2025_AOD_disk55_no_kipt_voxel_map_zup_robust.root}"
FLUX_TABLE="${FLUX_TABLE:-${PROJECT_DIR}/CosmicMuography/test_output/expacs_muon_flux_grid_cms_p5.csv}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/CosmicMuography/test_output/solid_ground_null_chunks_120m}"
WORK_DIR="${_CONDOR_SCRATCH_DIR:-${PWD}}"
LOCAL_OUT="${WORK_DIR}/solid_ground_expected_chunk_${CHUNK_INDEX}.root"
FINAL_OUT="${OUT_DIR}/solid_ground_expected_chunk_${CHUNK_INDEX}.root"
TMP_FINAL="${FINAL_OUT}.tmp.${HOSTNAME:-host}.$$"

cleanup() {
  rm -f "${LOCAL_OUT}" "${TMP_FINAL}"
}
trap cleanup EXIT

source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_SRC}"
eval "$(scram runtime -sh)"
cd "${PROJECT_DIR}"

mkdir -p "${OUT_DIR}"
python3 CosmicMuography/scripts/make_solid_ground_expectation.py \
  --data-root "${DATA_ROOT}" \
  --tracks-root "${TRACKS_ROOT}" \
  --flux-table "${FLUX_TABLE}" \
  --chunk-index "${CHUNK_INDEX}" \
  --nchunks "${NCHUNKS}" \
  --output "${LOCAL_OUT}"

python3 - <<PY
import ROOT
path = "${LOCAL_OUT}"
root_file = ROOT.TFile.Open(path)
hist = root_file.Get("solid_ground_expected_voxels") if root_file and not root_file.IsZombie() else None
if not hist:
    raise SystemExit(f"local chunk output is not readable: {path}")
print("validated local chunk", path, "integral", hist.Integral())
root_file.Close()
PY

cp -f "${LOCAL_OUT}" "${TMP_FINAL}"
python3 - <<PY
import ROOT
path = "${TMP_FINAL}"
root_file = ROOT.TFile.Open(path)
hist = root_file.Get("solid_ground_expected_voxels") if root_file and not root_file.IsZombie() else None
if not hist:
    raise SystemExit(f"copied chunk output is not readable: {path}")
print("validated copied chunk", path, "integral", hist.Integral())
root_file.Close()
PY
mv -f "${TMP_FINAL}" "${FINAL_OUT}"
echo "published final chunk ${FINAL_OUT}"
