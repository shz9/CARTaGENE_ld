#!/bin/bash
#SBATCH --account=def-sgravel
#SBATCH --job-name=cag-reference
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/reference-%j.out
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?Submit from the repository root}"
vcf=${1:?Usage: sbatch scripts/prepare_references.sbatch GRCh38.vcf.gz panel1 [panel2 ...]}
shift
[[ $# -gt 0 ]] || { echo "Supply at least one rsID panel" >&2; exit 1; }
env/ld/bin/python scripts/prepare_panel.py --assembly GRCh38 --dbsnp-vcf "$vcf" --rsids "$@"
env/ld/bin/python scripts/prepare_genetic_map.py
