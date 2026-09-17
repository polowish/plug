#!/bin/sh
# download BindingDB (curated protein-ligand binding affinities, articles subset) into data/bindingdb.
# the tsv's target sequence column holds each protein; affinity columns hold Ki/Kd/IC50/EC50.
set -e
d="$(dirname "$0")/bindingdb"
mkdir -p "$d"
# bindingdb rotates old releases off the server, so the YYYYMM in the filename goes stale.
# walk back from this month until one resolves.
base=https://www.bindingdb.org/rwd/bind/downloads/BindingDB_BindingDB_Articles
for back in 0 1 2 3 4 5 6 7 8 9 10 11; do
  m=$(date -u -v-${back}m +%Y%m 2>/dev/null || date -u -d "$back months ago" +%Y%m)
  curl -fsSL "${base}_${m}_tsv.zip" -o "$d/bindingdb.zip" && break
done
[ -s "$d/bindingdb.zip" ] || { echo "no bindingdb release found in the last 12 months" >&2; exit 1; }
cd "$d" && unzip -q -o bindingdb.zip && rm bindingdb.zip
