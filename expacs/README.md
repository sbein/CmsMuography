# EXPACS/PARMA wrapper

This directory contains the local CMS-cavern wrapper used to tabulate the
surface muon flux from the official EXPACS/PARMA C++ package.

Upstream package used on lxplus:

- EXPACS: https://phits.jaea.go.jp/expacs/
- C++ package: https://phits.jaea.go.jp/expacs/data/parma_cpp.zip

Please cite the upstream model if these outputs are used beyond exploratory
studies:

- T. Sato, PLOS ONE 10(12): e0144679 (2015), doi:10.1371/journal.pone.0144679
- T. Sato, PLOS ONE 11(8): e0160390 (2016), doi:10.1371/journal.pone.0160390

The first-pass null model in this analysis uses EXPACS/PARMA only for the
surface muon flux shape. The observed accepted-track sample provides the
detector acceptance kernel, and the expected volume is globally normalized to
the observed voxel integral before significance maps are computed.

When running `expacs_muon_grid`, pass `--expacs-dir` pointing at the unpacked
EXPACS/PARMA directory so the upstream `input/` tables are found reliably.
