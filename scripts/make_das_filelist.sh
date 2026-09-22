#!/usr/bin/env bash
set -euo pipefail

dataset="${1:-/Cosmics/pORun2025-PromptReco-v1/AOD}"
limit="${2:-100}"
redirector="${3:-root://cms-xrd-global.cern.ch/}"

if ! command -v dasgoclient >/dev/null 2>&1; then
  source /cvmfs/cms.cern.ch/cmsset_default.sh
  if command -v scram >/dev/null 2>&1; then
    eval "$(scram runtime -sh 2>/dev/null || true)"
  fi
fi

dasgoclient --query="file dataset=${dataset} instance=prod/global" --limit="${limit}" |
  awk -v redirector="${redirector}" '
    /^\/store\// { print redirector $0; next }
    NF { print }
  '
