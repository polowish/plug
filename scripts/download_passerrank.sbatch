#!/bin/bash
# passerrank on its own, calling the repo's data/download_passerrank.sh unmodified.
#
#   sbatch scripts/download_passerrank.sbatch
#
# it is the long pole of the benchmark set: one serial curl to the uniprot rest api per
# accession, thousands of them, with no resume — so a timeout restarts it from zero.
# isolating it means that timeout doesn't take the other 11 benchmarks down with it.

#SBATCH --job-name=plug-passerrank
#SBATCH --output=logs/passerrank_%j.out
#SBATCH --error=logs/passerrank_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
### ---- UCSF HPC Core: fill these in for your cluster ----
### #SBATCH --partition=CHANGEME
### #SBATCH --account=CHANGEME
### -------------------------------------------------------

set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
mkdir -p "$repo/logs"

command -v curl >/dev/null 2>&1 || { echo "missing required binary: curl" >&2; exit 1; }

echo "[$(date -u +%FT%TZ)] ==> passerrank (serial uniprot fetches, expect hours)"
sh "$repo/data/download_passerrank.sh"

fa="$repo/data/passerrank/passerrank.fasta"
echo "[$(date -u +%FT%TZ)] passerrank finished: $(grep -c '^>' "$fa") sequences in $fa"
