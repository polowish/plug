#!/bin/sh
# download passerrank allosteric proteins: pull the asd table, then fetch each unique
# pdb_uniprot accession's sequence from uniprot into one fasta.
set -e
d="$(dirname "$0")/passerrank"
mkdir -p "$d"
curl -fsSL https://raw.githubusercontent.com/smu-tao-group/PASSerRank/HEAD/data/source_data/ASD_Release_201909_AS.txt -o "$d/ASD_Release_201909_AS.txt"
awk -F'\t' 'NR>1 && $4 ~ /^[A-Z0-9]+$/ {print $4}' "$d/ASD_Release_201909_AS.txt" | sort -u > "$d/accessions.txt"
# accessions left out on purpose (../known_missing/passerrank.txt): inactive uniprot entries
# return an EMPTY fasta with http 200, so curl -f alone never noticed them missing
km="$(dirname "$0")/../known_missing/passerrank.txt"
: > "$d/passerrank.fasta"
: > "$d/failed_accessions.txt"
while read -r acc; do
  if [ -f "$km" ] && cut -f1 "$km" | grep -qx "$acc"; then continue; fi
  ok=""
  for t in 1 2 3; do
    if seq=$(curl -fsSL "https://rest.uniprot.org/uniprotkb/$acc.fasta") && [ -n "$seq" ]; then
      printf '%s\n' "$seq" >> "$d/passerrank.fasta"; ok=1; break
    fi
    sleep $((5 * t))
  done
  [ -n "$ok" ] || echo "$acc" >> "$d/failed_accessions.txt"
done < "$d/accessions.txt"
if [ -s "$d/failed_accessions.txt" ]; then
  echo "passerrank: $(wc -l < "$d/failed_accessions.txt") accessions gave no sequence after retries" \
       "(see $d/failed_accessions.txt) -- add them to known_missing/passerrank.txt with the reason" >&2
  exit 1
fi
# its site residues are numbered like the pdb file, so we need sifts to renumber them
sh "$(dirname "$0")/download_sifts.sh"
