#!/usr/bin/env bash
set -euo pipefail

did="${1:?usage: list_file_replicas.sh /store/...root|cms:/store/...root}"
protocol="${2:-root}"

source /cvmfs/cms.cern.ch/rucio/setup.sh

case "${did}" in
  cms:*) ;;
  /store/*) did="cms:${did}" ;;
esac

rucio list-file-replicas "${did}" --pfns --protocol "${protocol}"
