#!/bin/bash
# Backward-compatible alias; submission always goes through preparation/validation.
set -euo pipefail
cd "$(dirname "$0")/.."
config=${1:-config/pipeline.json}
array=${2:-1-22%4}
chromosomes=${array%%%*}
concurrency=4
if [[ "$array" == *%* ]]; then concurrency=${array#*%}; fi
if [[ $# -ge 2 ]]; then shift 2; else shift "$#"; fi
exec bash scripts/run_pipeline.sh --config "$config" --chromosomes "$chromosomes" \
  --max-concurrent "$concurrency" "$@"
