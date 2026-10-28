#!/bin/sh
# download cryptobench (cryptic binding sites: residues that line a pocket in a holo
# structure but not in the apo one) into data/cryptobench, then build the per-residue
# table plug reads (scripts/cryptobench_build.py: mmcif -> full chain sequence + labels).
# osf.io/pz4a9, MIT. ~1.2 GB, almost all of it the apo/holo mmcif files.
set -e
d="$(dirname "$0")/cryptobench"
mkdir -p "$d/folds"
get() { curl -fsSL --retry 3 "https://osf.io/download/$1/" -o "$2"; }
get ta2ju "$d/dataset.json"
get 5s93p "$d/folds.json"
get q2n43 "$d/folds/test.json"
get 7bvs4 "$d/folds/train-fold-0.json"
get akfgn "$d/folds/train-fold-1.json"
get nmj39 "$d/folds/train-fold-2.json"
get btdqy "$d/folds/train-fold-3.json"
[ -s "$d/cif-files.zip" ] || get c5vp8 "$d/cif-files.zip"
[ -d "$d/cif-files" ] || (cd "$d" && unzip -q cif-files.zip -d cif-files)
echo "downloaded; now build the table (biopython): python scripts/cryptobench_build.py"
