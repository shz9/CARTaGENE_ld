#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
config=${1:-config/pipeline.json}
array=${2:-1-22%4}
[[ -f "$config" ]] || { echo "Missing $config" >&2; exit 1; }
mkdir -p logs
sbatch --array="$array" scripts/chromosome.sbatch "$config"
