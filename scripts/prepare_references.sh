#!/bin/bash
# Run directly on the login node. Never submit this script with sbatch.
set -euo pipefail
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "Run reference preparation directly on the login node, outside a SLURM job." >&2
  exit 1
fi
cd "$(dirname "$0")/.."
config=config/pipeline.json
if [[ "${1:-}" == "--config" ]]; then
  config=${2:?Supply a configuration file after --config}
  shift 2
fi
vcf=${1:?Usage: bash scripts/prepare_references.sh [--config pipeline.json] GRCh38.vcf.gz [panel1 ...]}
shift
args=(--config "$config" --prepare-only --dbsnp-vcf "$vcf")
if [[ $# -gt 0 ]]; then args+=(--rsids "$@"); fi
exec bash scripts/run_pipeline.sh "${args[@]}"
