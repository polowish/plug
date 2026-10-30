#!/bin/sh
# download CheZOD (per-residue protein order/disorder from NMR chemical shifts: a Z-score
# per residue, high = ordered, ~<8 = disordered; Nielsen & Mulder 2016/2020) into
# data/chezod, then build the table plug reads (scripts/chezod_build.py).
#
# the copy used is ADOPT's (Redl et al. 2023, github.com/PeptoneLtd/ADOPT, MIT), pinned to
# the commit below; the raw sets come from github.com/protein-nmr/CheZOD (no licence file).
#   117_dataset_raw.json   CheZOD117, the test set (117 BMRB entries)
#   chezod_1159.json       CheZOD1325 with the 166 entries that duplicate or resemble the
#                          117 removed (85 of the 117 BMRB ids are IN the 1325 set) -> train/val
#   1325_dataset_raw.json  the unfiltered 1325, kept for reference only; never train on it
# 999.0 in a z-score list means no chemical shifts for that residue. ~2.5 MB in all.
set -e
d="$(dirname "$0")/chezod"
mkdir -p "$d"
sha=5053812037e215eb4d8aa54aef597892fcd27ccd
for f in 117_dataset_raw.json chezod_1159.json 1325_dataset_raw.json; do
    curl -fsSL --retry 3 "https://raw.githubusercontent.com/PeptoneLtd/ADOPT/$sha/datasets/$f" -o "$d/$f"
done
echo "downloaded; now build the table: python scripts/chezod_build.py"
