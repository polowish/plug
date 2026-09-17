#!/bin/bash
# fetch the 11 benchmarks + the sifts table into plug/data/, where the package hardcodes
# them (benchmarks/*.py all read c.data / "<name>").
#
#   sbatch scripts/download_benchmarks.sbatch
#
# deliberately ONE sequential job rather than an array: the whole set is ~1.4 GB and most
# scripts finish in seconds, so there is nothing to gain from parallelism — and running
# them in sequence means the three scripts that call download_sifts.sh (standalone,
# allobench, passerrank) cannot race on the same output file.
#
# skipped here: download_pdb.sh (the mirror goes to $DEST, see download_reservoirs.sh) and
# download_passerrank.sh (thousands of serial uniprot calls — its own job, so a timeout
# doesn't take the other benchmarks with it).

#SBATCH --job-name=plug-benchmarks
#SBATCH --output=logs/benchmarks_%j.out
#SBATCH --error=logs/benchmarks_%j.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=20:00
#SBATCH --partition=cpu

set -euo pipefail
module load CBI
module load miniforge3
conda activate /home/remote/smadejski/scratch/venvs/plug


here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(dirname "$here")"
mkdir -p "$repo/logs"

# fail now rather than partway through the set
for b in curl unzip python3; do
  command -v "$b" >/dev/null 2>&1 || { echo "missing required binary: $b" >&2; exit 1; }
done
# cafa5 needs the huggingface cli; download_cafa5.sh prefers .venv/bin/hf, else PATH
if [ ! -x "$repo/.venv/bin/hf" ] && ! command -v hf >/dev/null 2>&1; then
  echo "missing 'hf' (huggingface-hub cli) — activate the plug env first, or cafa5 will fail" >&2
  exit 1
fi

for s in "$repo/data"/download_*.sh; do
  case "$s" in
    *download_pdb.sh|*download_passerrank.sh) continue ;;
  esac
  name="$(basename "$s" .sh)"; name="${name#download_}"
  echo "[$(date -u +%FT%TZ)] ==> $name"
  sh "$s"
done

echo "[$(date -u +%FT%TZ)] benchmarks finished; data in $repo/data"
du -sh "$repo/data"/*/ 2>/dev/null || true
