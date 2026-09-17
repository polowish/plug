# scripts/ — cluster download jobs

SLURM wrappers for fetching PLUG's data on a cluster. Nothing here changes the PLUG
package — the benchmark downloaders in `data/` are called unmodified, and the two large
reservoirs use `UNIREF_FASTA` / `PDB_DIR`, which `config.py` already supports.

## Where each dataset goes, and why

| data | size | destination | why |
|---|---|---|---|
| 10 benchmarks + sifts | ~1.4 GB | `plug/data/` | the package hardcodes these — every `benchmarks/*.py` reads `c.data / "<name>"` |
| passerrank | small | `plug/data/` | same, but slow enough to need its own job |
| uniref90 | ~30 GB gzipped | `$DEST` | large and general-purpose — worth sharing across projects |
| PDB mmCIF mirror | tens of GB | `$DEST` | same |

`$DEST` defaults to `/mnt/fac/CX500059_DS1/smadejski/data`. Edit the marked line at the top
of `download_reservoirs.sh`, or set `DEST=` in the environment for a one-off run.

## Before the first submit

1. **Create the log dir**: `mkdir -p plug/logs`. `#SBATCH --output` is resolved by SLURM
   *before* the job starts, so the scripts cannot create it themselves — without it SLURM
   silently drops the output.
2. **Fill in the SLURM placeholders**: `--partition` and `--account` are commented out in
   each `.sbatch` header, since queue names are cluster-specific.
3. **Activate the environment** (for `download_benchmarks.sbatch` only): the cafa5 step
   needs the `hf` CLI. Either activate the `plug` env from `yamls/plug.yml` or make `hf`
   resolvable on `PATH`. The job preflight-checks this and exits immediately if missing.

## Usage

```bash
mkdir -p logs

sbatch scripts/download_benchmarks.sbatch     # ~1.4 GB, start here — fastest feedback
sbatch scripts/download_passerrank.sbatch     # hours, runs independently
sbatch scripts/download_reservoirs.sbatch     # both reservoirs, one array task each

sbatch --array=0 scripts/download_reservoirs.sbatch   # uniref90 only
sbatch --array=1 scripts/download_reservoirs.sbatch   # pdb only
```

Outside SLURM the reservoir script runs standalone:

```bash
sh scripts/download_reservoirs.sh uniref
DEST=/tmp/test sh scripts/download_reservoirs.sh uniref   # override for a test
FORCE=1 sh scripts/download_reservoirs.sh uniref          # re-download from scratch
```

## After the reservoirs land

Add to `plug/.env` (the script prints these lines when it finishes):

```
UNIREF_FASTA=/mnt/fac/CX500059_DS1/smadejski/data/uniref/uniref90.fasta.gz
PDB_DIR=/mnt/fac/CX500059_DS1/smadejski/data/pdb
```

Then the normal pipeline works:

```bash
python -m plug.fastas
RESERVOIR=seq python -m plug.build_unlabeled_trainset
```

## Design notes

**Why the benchmarks are one sequential job, not an array.** The whole set is ~1.4 GB and
most scripts finish in seconds, so parallelism buys nothing — and `download_sifts.sh` is
invoked by three different scripts (standalone, allobench, passerrank). Run concurrently
they would race on the same output file. Sequential execution removes the race for free.

**Why passerrank is separate.** It issues one serial UniProt REST call per accession,
thousands of them, with no resume — a timeout restarts it from zero. On its own job, that
timeout doesn't take the other 11 benchmarks with it.

**Why uniref90 downloads to `.part`.** The final `uniref90.fasta.gz` name appears only
after a *complete* download, so the skip guard can never mistake a truncated file for a
finished one. This matters: an earlier version checked the real filename, and a job killed
at 33 MB was reported as "already downloaded" on the next run — which would have fed
`build_unlabeled_trainset` a silently truncated reservoir. `curl -C -` resumes the `.part`.

**Why uniref90 stays gzipped.** `config.iter_fasta` opens `.gz` transparently, so
decompressing costs ~3× the disk for no benefit.

**Why the PDB rsync keeps `--delete`.** It is correct for a mirror you re-sync, and safe
only because the target is a dedicated `pdb/` subdirectory. Never point it at `$DEST`
itself. One stream only — RCSB throttles concurrent connections.

**Skip guards.** Every target is skipped if already present, overridable with `FORCE=1`.
Worth knowing: of the repo's own downloaders only `download_atlas.sh` is natively
idempotent — the rest re-download unconditionally.

## Not covered

`STRUCT_TESTS` (`data/structs`) — the test-set structures for the foldseek structural
leakage check. No downloader exists for it upstream; it is an open TODO in the main README,
so `RESERVOIR=struct|both` cannot run end-to-end until that is solved separately.
