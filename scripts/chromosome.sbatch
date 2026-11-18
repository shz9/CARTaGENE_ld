#!/bin/bash
#SBATCH --account=def-sgravel
#SBATCH --job-name=cag-ld
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=05:00:00
#SBATCH --output=logs/ld-%A_%a.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
# Use the same Python module for venv creation and execution.
module load "${PYTHON_MODULE:-python/3.13}"
module load plink
plink_bin=$(command -v plink2 || command -v plink)
"$plink_bin" --version
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
env/ld/bin/python -u scripts/run_chromosome.py --config "${1:-config/pipeline.json}" \
  --chrom "${SLURM_ARRAY_TASK_ID:?Missing array chromosome}" --plink "$plink_bin" \
  --threads "$OMP_NUM_THREADS"
