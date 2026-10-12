# Input lists

Generate a small XRootD list from DAS after loading CMSSW:

```bash
./scripts/make_das_filelist.sh /Cosmics/OORun2025-PromptReco-v1/AOD 100 > data/OORun2025_AOD_100.txt
```

The analyzer accepts ROOT URLs, `/store/...` LFNs, or text files containing
one input per line. Site-specific PFN inventories, Rucio reports, and large
production lists are intentionally ignored by Git because they can expire or
depend on an individual account. Keep those local under `data/`.

The two checked-in examples are only lightweight direct-AOD smoke inputs; they
may eventually require regeneration as replicas move.
