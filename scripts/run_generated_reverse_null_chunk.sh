#!/usr/bin/env bash
set -euo pipefail

CHUNK_INDEX="${1:?chunk index}"
NCHUNKS="${2:?number of chunks}"

PROJECT_DIR="${PROJECT_DIR:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern}"
CMSSW_SRC="${CMSSW_SRC:-${PROJECT_DIR}/CMSSW_16_1_0_pre4/src}"
TEMPLATE_ROOT="${TEMPLATE_ROOT:-${PROJECT_DIR}/CosmicMuography/test_output/OORun2025_AOD_disk55_no_kipt_zonly_smoothcalib_lima_alpha0p05_120m_rhcavern.root}"
TEMPLATE_HIST="${TEMPLATE_HIST:-observed_voxels}"
FLUX_TABLE="${FLUX_TABLE:-${PROJECT_DIR}/CosmicMuography/test_output/expacs_muon_flux_grid_cms_p5.csv}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/CosmicMuography/test_output/generated_reverse_null_chunks_120m}"
NTHROWS="${NTHROWS:-20000000}"
MIN_SURFACE_ENERGY_MEV="${MIN_SURFACE_ENERGY_MEV:-30000}"
MIN_COSZ="${MIN_COSZ:-0.1}"
SEED="${SEED:-20260822}"
WORK_DIR="${_CONDOR_SCRATCH_DIR:-${PWD}}"
LOCAL_OUT="${WORK_DIR}/generated_reverse_expected_chunk_${CHUNK_INDEX}.root"
FINAL_OUT="${OUT_DIR}/generated_reverse_expected_chunk_${CHUNK_INDEX}.root"
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
python3 CosmicMuography/scripts/make_generated_solid_ground_expectation.py \
  --sampling-mode reverse-cms-volume \
  --template-root "${TEMPLATE_ROOT}" \
  --template-hist "${TEMPLATE_HIST}" \
  --flux-table "${FLUX_TABLE}" \
  --chunk-index "${CHUNK_INDEX}" \
  --nchunks "${NCHUNKS}" \
  --nthrows "${NTHROWS}" \
  --seed "${SEED}" \
  --min-surface-energy-mev "${MIN_SURFACE_ENERGY_MEV}" \
  --min-cosz "${MIN_COSZ}" \
  --output "${LOCAL_OUT}"

python3 - <<PY
import ROOT
path = "${LOCAL_OUT}"
root_file = ROOT.TFile.Open(path)
hist = root_file.Get("generator_expected_voxels") if root_file and not root_file.IsZombie() else None
if not hist:
    raise SystemExit(f"local chunk output is not readable: {path}")
print("validated local generated chunk", path, "integral", hist.Integral())
root_file.Close()
PY

cp -f "${LOCAL_OUT}" "${TMP_FINAL}"
python3 - <<PY
import ROOT
path = "${TMP_FINAL}"
root_file = ROOT.TFile.Open(path)
hist = root_file.Get("generator_expected_voxels") if root_file and not root_file.IsZombie() else None
if not hist:
    raise SystemExit(f"copied chunk output is not readable: {path}")
print("validated copied generated chunk", path, "integral", hist.Integral())
root_file.Close()
PY
mv -f "${TMP_FINAL}" "${FINAL_OUT}"
echo "published final generated chunk ${FINAL_OUT}"
