#!/bin/bash
set -euo pipefail
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
  echo "Create the environment directly on the login node, outside a SLURM job." >&2
  exit 1
fi
cd "$(dirname "$0")/.."
# Use the same Python module for venv creation and execution.
module load "${PYTHON_MODULE:-python/3.13}"
python3 -m venv env/ld
env/ld/bin/python -m pip install --upgrade pip
env/ld/bin/python -m pip install 'magenpy==0.2.5' 'zarr<3'
env/ld/bin/python -m pip freeze > env/ld/requirements.lock.txt
