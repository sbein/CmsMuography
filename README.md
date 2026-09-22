# CMS Cavern Muography

This repository contains a reproducible starter pipeline for CMS cavern
muography with 2025 CMS cosmic AOD data. It builds a z-up voxel occupancy map
from reconstructed cosmic muons, constructs an EXPACS/PARMA-inspired solid
overburden null, computes Li-Ma/profile-likelihood cavity evidence volumes,
and makes decorated ROOT canvases with CMS/cavern/ground/compass overlays.

The code is intentionally lightweight: CMSSW/FWLite/PyROOT for AOD and ROOT
histograms, plus small Python scripts and Condor wrappers. Large generated
ROOT files and PNG displays are not committed.

## Coordinate Convention

The analysis display frame is the cavern frame:

- `x = CMS global x`
- `y = CMS global z`, the beam-pipe axis, horizontal through the barrel
- `z = CMS global y`, physically upward out of the cavern

Histograms are stored internally in cm unless a display script is asked for
`--plot-units m`. The final display convention uses meters. Muon ray painting
is clipped to the upward side by default, so the map represents trajectories
propagated toward positive cavern `z`.

## Repository Layout

- `scripts/`: AOD reading, voxel painting, null generation, Li-Ma volumes,
  display canvases, and Condor wrappers.
- `data/`: collaborator-local DAS/Rucio/XRootD file lists.
- `expacs/`: local wrapper notes for EXPACS/PARMA surface flux tables.
- `docs/`: smoke-test route and collaborator notes.

## Quick Smoke Test

For a new collaborator, start with the synthetic end-to-end smoke pipeline:

```bash
RUNTIME=current ./scripts/run_smoke_pipeline.sh
```

It exercises the EXPACS-format flux interface, propagation/null generation,
track-based voxel reconstruction, Li-Ma calculation, broad calibration,
and decorated visualization. See `docs/smoke_tests.md` for the lxplus setup
and optional direct-AOD smoke.

## Environment On Lxplus

```bash
cd /afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src && source /cvmfs/cms.cern.ch/cmsset_default.sh && eval "$(scram runtime -sh)" && cd /afs/cern.ch/work/s/sbein/Muography/CmsCavern
```

For Rucio:

```bash
source /cvmfs/cms.cern.ch/rucio/setup.sh && rucio rule list --account bbilin && rucio rule list --account iyenioca
```

For remote AOD reads, the FNAL redirector was often robust:

```text
root://cmsxrootd.fnal.gov/
```

## Input Cosmic Dataset

The main first-pass sample was:

```text
/Cosmics/OORun2025-PromptReco-v1/AOD
```

Per-user DAS/Rucio/XRootD file lists are generated into `data/` and are not
tracked by Git. In Samuel's lxplus working area one useful local file list has
been:

```text
data/OORun2025_AOD_disk_root_one_pfn_per_file_xrdfs_readable_no_kipt.txt
```

The code also has helpers to regenerate or validate file lists:

```bash
./scripts/make_das_filelist.sh /Cosmics/OORun2025-PromptReco-v1/AOD 100 > data/OORun2025_AOD_100.txt && python3 scripts/check_xrootd_filelist.py data/OORun2025_AOD_100.txt --timeout 20
```

## Step 1: Build A Raw Voxel Map From AOD

Small direct test:

```bash
python3 scripts/cosmic_voxel_map.py --input data/OORun2025_AOD_100.txt --redirector root://cmsxrootd.fnal.gov/ --max-events 1000 --output test_output/OORun2025_small_voxel_map_zup.root --skip-bad-files --progress-every 100
```

Robust per-file run, useful when individual disk replicas stall:

```bash
PROJECT_ROOT=$PWD INPUT_LIST=data/OORun2025_AOD_disk_root_one_pfn_per_file_xrdfs_readable_no_kipt.txt MERGED_ROOT=test_output/OORun2025_AOD_voxel_map_zup_robust.root PEDESTAL_ROOT=test_output/OORun2025_AOD_voxel_map_zup_robust_pedestals.root ./scripts/run_oorun2025_disk58_per_file.sh
```

The raw map file contains:

- `cavern_muon_voxels`: TH3 of voxelized muon-line crossings.
- `tracks`: accepted-track diagnostics, including cavern-frame reference
  point and direction.
- `cutflow`: event/track/voxel counts.
- `config_json`: exact geometry and cuts.

If a run produced a good `tracks` tree but the TH3 needs to be rebuilt:

```bash
python3 scripts/rebuild_voxel_map_from_tracks.py --input test_output/OORun2025_AOD_voxel_map_zup_robust.root --output test_output/OORun2025_AOD_voxel_map_zup_rebuilt.root
```

## Step 2: Make The Solid-Ground Null

The preferred null is independent of the real angular distribution except for
broad later calibration. It samples EXPACS/PARMA surface muons, samples random
impact points through the finite CMS cylinder, propagates straight lines
through a uniform solid overburden to the ground surface, and computes survival
with a Bethe/radiative range model.

The local EXPACS wrapper is in `expacs/`. The project expects a CSV flux table
like:

```text
test_output/expacs_muon_flux_grid_cms_p5.csv
```

Generate one null file directly:

```bash
python3 scripts/make_generated_solid_ground_expectation.py --sampling-mode reverse-cms-volume --template-root test_output/OORun2025_AOD_voxel_map_zup_robust.root --template-hist cavern_muon_voxels --flux-table test_output/expacs_muon_flux_grid_cms_p5.csv --nthrows 2000000 --min-surface-energy-mev 30000 --min-cosz 0.1 --output test_output/generated_solid_ground_bethe_reverse_2M.root
```

For large statistics, use the Condor wrapper:

```bash
condor_submit scripts/submit_generated_reverse_null_120m.condor
```

If chunked outputs are produced, merge their `generator_expected_voxels`
histograms with:

```bash
python3 scripts/merge_voxel_maps.py --hist-name generator_expected_voxels --input-list test_output/generated_reverse_null_chunks_120m/successful_parts.txt --output test_output/generated_solid_ground_bethe_reverse_merged.root --skip-bad-files
```

## Step 3: Compute Li-Ma / Cavity Evidence Volumes

Raw Li-Ma or fixed-Poisson significance:

```bash
python3 scripts/make_lima_volume.py --observed-root test_output/OORun2025_AOD_voxel_map_zup_robust.root --observed-hist cavern_muon_voxels --expected-root test_output/generated_solid_ground_bethe_reverse_2M.root --expected-hist generator_expected_voxels --method lima --alpha 0.05 --normalize global --output test_output/OORun2025_classic_lima_alpha0p05.root
```

Broad smooth calibration, deliberately too smooth to erase cavity-scale
features:

```bash
python3 scripts/calibrate_lima_expectation.py -i test_output/OORun2025_classic_lima_alpha0p05.root -o test_output/OORun2025_smoothcalib_lima_alpha0p05.root --method lima --alpha 0.05 --calibration-model zxy-poly --x-smooth-cm 2500 --y-smooth-cm 2500 --z-smooth-cm 1500
```

Important output histograms include:

- `observed_voxels`
- `expected_voxels_scaled`
- `expected_voxels_calibrated`
- `signed_significance` or `signed_significance_calibrated`
- `cavity_excess_significance` or `cavity_excess_significance_calibrated`
- `log_transmission_excess`
- `observed_over_expected_calibrated`
- `smooth_calibration_factor`

## Step 4: Make Decorated Displays

```bash
python3 scripts/make_lima_slice_stacks.py -i test_output/OORun2025_smoothcalib_lima_alpha0p05.root -o test_output/OORun2025_smoothcalib_lima_alpha0p05_slice_stacks.root --hist-names cavity_excess_significance_calibrated,signed_significance_calibrated,log_transmission_excess,observed_over_expected_calibrated --copy-hists observed_voxels,expected_voxels_scaled,expected_voxels_calibrated,smooth_calibration_factor --plot-units m --display-rebin 2 --log-3d-color --png-dir test_output/OORun2025_smoothcalib_lima_alpha0p05_slice_stacks_png
```

The display script writes ROOT canvases and optional PNGs with:

- CMS detector wireframe in black.
- Approximate CMS cavern hall and PX56 shaft wireframe in black.
- Approximate ground surface overlay.
- 3D compass rose that rotates with the 3D scene.
- XY slices in z layers and XZ slices in y layers.
- 3D threshold grids using loose-to-tighter positive-value quantiles.

## Create A GitHub Remote

This local repo can be pushed from a machine with GitHub authentication:

```bash
cd cms-cavern-muography && gh repo create cms-cavern-muography --private --source=. --remote=origin --push
```

Without `gh`, create an empty GitHub repository in the browser, then:

```bash
git remote add origin git@github.com:USER_OR_ORG/cms-cavern-muography.git && git push -u origin main
```

## Reproducibility Notes

- Do not commit large outputs; keep ROOT and PNG products in `test_output/`,
  EOS, or another analysis storage area.
- The display geometry is approximate, documented in the ROOT metadata, and
  meant for interpretation rather than survey-grade CAD.
- Internal ROOT axes are usually cm. Use `--plot-units m` for final displays.
- Avoid fitting calibration fields so flexibly that they absorb 50 m-scale
  cavities or vertically segmented variations.
