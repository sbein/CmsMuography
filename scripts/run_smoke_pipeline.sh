#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
PYTHON="${PYTHON:-python3}"
OUT_DIR="${OUT_DIR:-${PROJECT_DIR}/test_output/smoke_pipeline}"
RUNTIME="${RUNTIME:-current}"
CMSSW_SRC="${CMSSW_SRC:-/afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src}"
RUN_AOD="${RUN_AOD:-0}"
AOD_INPUT="${AOD_INPUT:-${PROJECT_DIR}/data/test_OORun2025_AOD_nonempty.txt}"
REDIRECTOR="${REDIRECTOR:-root://cmsxrootd.fnal.gov/}"

if [[ "${RUNTIME}" == "cmssw" ]]; then
  set +u
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  set -u
  cd "${CMSSW_SRC}"
  eval "$(scram runtime -sh)"
elif [[ "${RUNTIME}" != "current" ]]; then
  echo "Unsupported RUNTIME=${RUNTIME}; use current or cmssw." >&2
  exit 2
fi

cd "${PROJECT_DIR}"
mkdir -p "${OUT_DIR}"

echo "[0/6] environment"
echo "project=${PROJECT_DIR}"
echo "out=${OUT_DIR}"
echo "python=$(${PYTHON} -c 'import sys; print(sys.executable)')"
"${PYTHON}" -V
"${PYTHON}" - <<'PY'
import importlib.util
print("ROOT_available", importlib.util.find_spec("ROOT") is not None)
PY

echo "[1/6] EXPACS-format flux table and synthetic track fixtures"
"${PYTHON}" scripts/make_smoke_fixtures.py -o "${OUT_DIR}/fixtures"
FLUX_TABLE="${OUT_DIR}/fixtures/smoke_expacs_flux.csv"
TRACKS_ROOT="${OUT_DIR}/fixtures/smoke_tracks.root"

echo "[2/6] reconstruction/voxel-painting smoke via tracks-tree rebuild"
"${PYTHON}" scripts/rebuild_voxel_map_from_tracks.py \
  --input "${TRACKS_ROOT}" \
  --output "${OUT_DIR}/observed_rebuilt_voxels.root" \
  --xmin -2000 --xmax 2000 \
  --ymin -3000 --ymax 3000 \
  --zmin 0 --zmax 4000 \
  --voxel-size 500 \
  --progress-every 10 \
  --copy-tracks

if [[ "${RUN_AOD}" == "1" ]]; then
  echo "[2b/6] optional direct AOD reconstruction smoke"
  "${PYTHON}" scripts/cosmic_voxel_map.py \
    --input "${AOD_INPUT}" \
    --redirector "${REDIRECTOR}" \
    --max-events 200 \
    --output "${OUT_DIR}/observed_aod_voxels.root" \
    --skip-bad-files \
    --progress-every 50
fi

echo "[3/6] propagation/null smoke with reverse CMS-volume solid-ground model"
"${PYTHON}" scripts/make_generated_solid_ground_expectation.py \
  --sampling-mode reverse-cms-volume \
  --template-root "${OUT_DIR}/observed_rebuilt_voxels.root" \
  --template-hist cavern_muon_voxels \
  --flux-table "${FLUX_TABLE}" \
  --nthrows 2500 \
  --seed 12345 \
  --min-surface-energy-mev 30000 \
  --min-cosz 0.2 \
  --range-grid-points 400 \
  --progress-every 1000 \
  --output "${OUT_DIR}/generated_solid_ground_smoke.root"

echo "[4/6] Li-Ma/profile-likelihood volume smoke"
"${PYTHON}" scripts/make_lima_volume.py \
  --observed-root "${OUT_DIR}/observed_rebuilt_voxels.root" \
  --observed-hist cavern_muon_voxels \
  --expected-root "${OUT_DIR}/generated_solid_ground_smoke.root" \
  --expected-hist generator_expected_voxels \
  --method lima \
  --alpha 0.05 \
  --normalize global \
  --output "${OUT_DIR}/lima_raw_smoke.root"

echo "[5/6] broad smooth calibration smoke"
"${PYTHON}" scripts/calibrate_lima_expectation.py \
  -i "${OUT_DIR}/lima_raw_smoke.root" \
  -o "${OUT_DIR}/lima_smoothcalib_smoke.root" \
  --method lima \
  --alpha 0.05 \
  --min-expected 1e-12 \
  --calibration-model product1d \
  --x-smooth-cm 1000 \
  --y-smooth-cm 1000 \
  --z-smooth-cm 1000

echo "[6/6] visualization smoke"
"${PYTHON}" scripts/make_lima_slice_stacks.py \
  -i "${OUT_DIR}/lima_smoothcalib_smoke.root" \
  -o "${OUT_DIR}/lima_smoothcalib_smoke_slice_stacks.root" \
  --hist-names cavity_excess_significance_calibrated,signed_significance_calibrated,log_transmission_excess \
  --copy-hists observed_voxels,expected_voxels_global_scaled,calibrated_expected_voxels,smooth_calibration_factor \
  --plot-units m \
  --display-rebin 1 \
  --log-3d-color \
  --z-slices-cm 250,750,1250,1750,2250,2750,3250,3750 \
  --y-slices-cm=-2500,-1500,-750,0,750,1500,2500 \
  --png-dir "${OUT_DIR}/png"

echo "smoke pipeline complete"
echo "outputs: ${OUT_DIR}"
