#!/bin/bash
# fetch the two large reservoirs, one per array task, so uniref90 and the pdb mirror run
# concurrently and each gets its own log + its own re-run handle (sbatch --array=1 ...).
#
#   sbatch scripts/download_reservoirs.sbatch          # both
#   sbatch --array=0 scripts/download_reservoirs.sbatch  # uniref90 only
#
# both are many-hour transfers; --time below is sized for that.

#SBATCH --job-name=plug-reservoirs
#SBATCH --array=0-1
#SBATCH --output=logs/reservoirs_%A_%a.out
#SBATCH --error=logs/reservoirs_%A_%a.out
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=48:00:00
### ---- UCSF HPC Core: fill these in for your cluster ----
### #SBATCH --partition=CHANGEME
### #SBATCH --account=CHANGEME
### -------------------------------------------------------

set -euo pipefail

# resolve the repo from this script's location, so sbatch can be run from anywhere
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
mkdir -p "$repo/logs"

targets=(uniref pdb)
target="${targets[${SLURM_ARRAY_TASK_ID:-0}]}"

echo "[$(date -u +%FT%TZ)] host=$(hostname) task=${SLURM_ARRAY_TASK_ID:-0} target=$target"
sh "$here/download_reservoirs.sh" "$target"
echo "[$(date -u +%FT%TZ)] $target finished"
