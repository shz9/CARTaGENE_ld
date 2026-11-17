#!/bin/bash
# Unified login-node entry point: bootstrap, prepare, validate, then submit.
set -euo pipefail
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "Run this script directly on the login node, outside a SLURM job." >&2
  exit 1
fi
cd "$(dirname "$0")/.."
if [[ "${1:-}" == --help || "${1:-}" == -h ]]; then
  cat <<'HELP'
Usage: bash scripts/run_pipeline.sh [options]
  --config FILE          Pipeline configuration (default: config/pipeline.json)
  --chromosomes LIST     Autosomes, e.g. 22 or 1-3,22 (default: 1-22)
  --max-concurrent N     Maximum simultaneous chromosome jobs (default: 4)
  --rsids FILE [FILE...] Override config rsids (default: metadata/rsids.txt)
  --dbsnp-vcf FILE       Local GRCh38 VCF (default: data/reference/GCF_000001405.40.gz)
  --dbsnp-url URL        Download URL if the VCF is missing; fetches URL.md5 as well
  --check-only          Validate existing inputs; no installation/download/submission
  --prepare-only        Prepare and validate inputs without submitting
  --help                Show this help
All preparation runs on the login node. Successful normal runs submit a SLURM array.
Set PYTHON_MODULE to override python/3.13 consistently for setup and jobs.
HELP
  exit 0
fi
check_only=0
for arg in "$@"; do [[ "$arg" != --check-only ]] || check_only=1; done
module load "${PYTHON_MODULE:-python/3.13}"
module load plink
if ! [[ -x env/ld/bin/python ]] || ! env/ld/bin/python -c \
  'import sys, magenpy, numpy, pandas, zarr; assert magenpy.__version__ == "0.2.5"; assert zarr.__version__.split(".")[0] == "2"' 2>/dev/null; then
  if [[ "$check_only" == 1 ]]; then
    echo "Missing or invalid env/ld. Run without --check-only to set it up." >&2
    exit 1
  fi
  bash scripts/setup_env.sh
fi
# A module change must not silently reuse a venv based on another Python version.
[[ "$(python3 -c 'import sys; print(sys.version_info[:2])')" == \
   "$(env/ld/bin/python -c 'import sys; print(sys.version_info[:2])')" ]] || {
  echo "env/ld uses a different Python version; recreate it with the selected PYTHON_MODULE." >&2
  exit 1
}
env/ld/bin/python -m pip check
env/ld/bin/python -u scripts/run_pipeline.py "$@"
