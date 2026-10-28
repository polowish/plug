#!/bin/bash
# the compute-node half of the mdcath stream: extract every domain mdcath_fetch.py stages,
# then delete it (scripts/mdcath_extract.py --watch). CPU only; ~36 core-hours for all
# 5,398 domains at the measured rate (largest domain 8 GB RAM). Resumable: resubmit and it
# puts back anything a killed run had claimed.
#
#   sbatch scripts/mdcath_extract.sbatch [staging] [out]
#SBATCH --job-name=mdcath-extract
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/%x_%j.out
set -euo pipefail
STAGING="${1:-$HOME/scratch/mdcath_staging}"
OUT="${2:-$HOME/fac/data/mdcath}"
PY="${MDCATH_PY:-$HOME/scratch/venvs/mdcath/bin/python}"
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")/..}"
echo "[$(date -u +%FT%TZ)] host=$(hostname)  plug $(git rev-parse --short HEAD)  staging=$STAGING  out=$OUT"
# 8 extractions at once: the largest domains peak at ~8 GB each
"$PY" -u scripts/mdcath_extract.py --watch "$STAGING" --out "$OUT" --workers "${WORKERS:-8}"
