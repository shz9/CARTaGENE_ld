#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
python3 -m venv env/ld
env/ld/bin/python -m pip install --upgrade pip
env/ld/bin/python -m pip install 'magenpy==0.2.5' 'zarr<3'
env/ld/bin/python -m pip freeze > env/ld/requirements.lock.txt
