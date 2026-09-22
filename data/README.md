# Local Data Lists

This directory is reserved for collaborator-local DAS, Rucio, and XRootD file
lists. The public repository intentionally does not track PFN/LFN staging
reports or site-specific replica metadata.

To create a small AOD input list on lxplus:

```bash
./scripts/make_das_filelist.sh /Cosmics/OORun2025-PromptReco-v1/AOD 100 > data/OORun2025_AOD_100.txt && python3 scripts/check_xrootd_filelist.py data/OORun2025_AOD_100.txt --timeout 20
```

Then pass the readable file list explicitly, for example:

```bash
RUN_AOD=1 AOD_INPUT=data/OORun2025_AOD_100.txt ./scripts/run_smoke_pipeline.sh
```
