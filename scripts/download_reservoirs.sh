#!/bin/sh
# download the two large, reusable reservoirs to an external path — uniref90 (the sequence
# reservoir) and a pdb mmcif mirror (the structure reservoir). these live outside the repo
# because they are general-purpose datasets worth sharing across projects; the 11 benchmarks
# stay in plug/data where the package hardcodes them (see download_benchmarks.sbatch).
#
# usage:  sh download_reservoirs.sh [uniref|pdb]      # no arg does both
#         DEST=/other/path sh download_reservoirs.sh  # override the path for one run
#         FORCE=1 sh download_reservoirs.sh uniref    # re-download even if present
set -e

# ---- EDIT ME ----  (a DEST= in the environment overrides this, e.g. for a test run)
DEST="${DEST:-/mnt/fac/CX500059_DS1/smadejski/data}"
# -----------------

# fail before a multi-hour transfer rather than partway through it
need() {
  command -v "$1" >/dev/null 2>&1 || { echo "missing required binary: $1" >&2; exit 1; }
}

# a reservoir is "done" if its sentinel exists and is non-empty; FORCE=1 re-does it anyway
done_already() {
  [ -z "$FORCE" ] && [ -e "$1" ] && [ -s "$1" ]
}

uniref() {
  need curl
  d="$DEST/uniref"
  f="$d/uniref90.fasta.gz"
  if done_already "$f"; then
    echo "uniref90 already at $f (FORCE=1 to re-download)"
    return 0
  fi
  mkdir -p "$d"
  part="$f.part"
  [ -n "$FORCE" ] && rm -f "$part"   # FORCE means start over, not resume
  echo "==> uniref90 -> $f"
  # download into .part and rename only on success. a job that hits its time limit then
  # leaves a truncated .part that the next run resumes — and never a truncated
  # uniref90.fasta.gz that the skip guard above would mistake for a finished download
  # (which would silently feed build_unlabeled_trainset a partial reservoir).
  # -C - resumes: this is tens of gb and will likely span more than one attempt. left
  # gzipped on purpose — config.iter_fasta opens .gz transparently, so decompressing
  # would cost ~3x the disk for nothing.
  # -# is the compact progress bar: the default meter writes a wide status line every second,
  # which makes an 8-hour slurm log painful to tail.
  curl -fL -C - -# https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz \
    -o "$part"
  mv "$part" "$f"
  echo "uniref90 complete: $(du -h "$f" | cut -f1)"
  # keep the release note next to it so the reservoir version is on record
  curl -fsSL https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.release_note \
    -o "$d/uniref90.release_note" || true
}

pdb() {
  need rsync
  d="$DEST/pdb"
  if done_already "$d/.rsync_complete"; then
    echo "pdb mirror already at $d (FORCE=1 to re-sync)"
    return 0
  fi
  mkdir -p "$d"
  echo "==> pdb mmcif mirror -> $d"
  # --delete keeps the mirror in sync with rcsb, which is only safe because the target is a
  # dedicated pdb/ subdir — never point this at $DEST itself.
  # one stream only: rcsb throttles concurrent connections, and rsync resumes on its own.
  rsync -rlpt -z --delete --port=33444 \
    rsync.rcsb.org::ftp_data/structures/divided/mmCIF/ "$d/"
  touch "$d/.rsync_complete"
}

case "${1:-all}" in
  uniref) uniref ;;
  pdb)    pdb ;;
  all)    uniref; pdb ;;
  *)      echo "usage: $0 [uniref|pdb]" >&2; exit 2 ;;
esac

echo "done. point plug at them in .env:"
echo "  UNIREF_FASTA=$DEST/uniref/uniref90.fasta.gz"
echo "  PDB_DIR=$DEST/pdb"
