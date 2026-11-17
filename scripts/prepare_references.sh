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
vcf=${1:?Usage: bash scripts/prepare_references.sh [--config pipeline.json] GRCh38.vcf.gz panel1 [panel2 ...]}
shift
[[ $# -gt 0 ]] || { echo "Supply at least one rsID panel" >&2; exit 1; }
[[ -f "$vcf" ]] || { echo "Download the GRCh38 dbSNP VCF on the login node first: $vcf" >&2; exit 1; }
# Use the same Python module for venv creation and execution.
module load "${PYTHON_MODULE:-python/3.13}"
# Small downloads first; the local dbSNP scan can take considerably longer.
env/ld/bin/python scripts/prepare_genetic_map.py
env/ld/bin/python scripts/prepare_blocks.py --config "$config"
env/ld/bin/python scripts/prepare_panel.py --assembly GRCh38 --dbsnp-vcf "$vcf" --rsids "$@"
