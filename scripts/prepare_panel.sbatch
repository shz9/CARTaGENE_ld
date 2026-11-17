#!/bin/bash
#SBATCH --account=def-sgravel
#SBATCH --job-name=cag-panel
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=05:00:00
#SBATCH --output=logs/panel-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
module load "${PYTHON_MODULE:-python/3.13}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
vcf=${1:?Missing local GRCh38 dbSNP VCF}
panel=${2:?Missing panel output path}
shift 2
[[ $# -gt 0 ]] || { echo "Supply at least one rsID list" >&2; exit 1; }
# All inputs must be local: compute nodes have no internet access.
env/ld/bin/python -u scripts/prepare_panel.py --assembly GRCh38 \
  --dbsnp-vcf "$vcf" --output "$panel" --rsids "$@"
