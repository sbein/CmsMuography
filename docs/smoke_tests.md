# Smoke Tests

The fastest verified route for a new collaborator is the synthetic smoke
pipeline:

```bash
cd /afs/cern.ch/work/s/sbein/Muography/CmsMuography && source /cvmfs/cms.cern.ch/cmsset_default.sh && cd /afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src && eval "$(scram runtime -sh)" && cd /afs/cern.ch/work/s/sbein/Muography/CmsMuography && RUNTIME=current ./scripts/run_smoke_pipeline.sh
```

This tests, in order:

1. EXPACS/PARMA table interface using a tiny EXPACS-format CSV.
2. Reconstruction-side voxel painting using a synthetic accepted-track ROOT
   tree and `rebuild_voxel_map_from_tracks.py`.
3. Straight-line reverse propagation through the finite CMS cylinder and
   solid-ground survival model.
4. Li-Ma/profile-likelihood volume construction.
5. Broad smooth null calibration.
6. Decorated ROOT visualization with CMS/cavern/ground/compass overlays.

The smoke output defaults to:

```text
test_output/smoke_pipeline
```

The synthetic route is deliberately independent of grid proxies and AOD
replica availability. Once that passes, run the optional direct-AOD smoke:

```bash
RUNTIME=current RUN_AOD=1 AOD_INPUT=data/OORun2025_AOD_100.txt ./scripts/run_smoke_pipeline.sh
```

That reads the file list named by `AOD_INPUT` through
`root://cmsxrootd.fnal.gov/`, so it needs a valid CMS runtime and whatever
proxy/access setup is required by the site. See `data/README.md` for creating
a small local list; PFN/LFN replica reports are intentionally not tracked in
the public repository.

For an lxplus batch-style environment that has not already sourced CMSSW:

```bash
RUNTIME=cmssw CMSSW_SRC=/afs/cern.ch/work/s/sbein/Muography/CmsCavern/CMSSW_16_1_0_pre4/src ./scripts/run_smoke_pipeline.sh
```

Successful completion should produce:

- `observed_rebuilt_voxels.root`
- `generated_solid_ground_smoke.root`
- `lima_raw_smoke.root`
- `lima_smoothcalib_smoke.root`
- `lima_smoothcalib_smoke_slice_stacks.root`
- `png/*.png`
