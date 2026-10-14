#!/bin/sh
# download atlas MD per-residue RMSF (conformational dynamics) into data/atlas.
# atlas ships each protein chain as one ~70MB analysis zip that's almost all trajectory —
# we http-range-extract just the ~5KB <pdb>_RMSF.tsv out of each remote zip, no trajectories pulled.
# each RMSF.tsv is one row per residue: seq (1-letter aa), RMSF_R1, RMSF_R2, RMSF_R3 (Å, 3 MD replicas).
set -e
d="$(dirname "$0")/atlas"
mkdir -p "$d/rmsf"
# the list of PDB chains in ATLAS v2 (Nov 2024)
curl -fsSL https://www.dsimb.inserm.fr/ATLAS/data/download/distributions/2024_11_18_ATLAS_pdb.txt -o "$d/atlas_pdb.txt"
python3 - "$d" <<'PY'
import io, sys, time, urllib.request, zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

d = Path(sys.argv[1]); out = d / "rmsf"
api = "https://www.dsimb.inserm.fr/ATLAS/api/ATLAS/analysis/{}"
chains = [ln.strip() for ln in open(d / "atlas_pdb.txt") if ln.strip()]
# chains left out on purpose (../known_missing/atlas.txt): not fetched, not failures
km = d.resolve().parent.parent / "known_missing" / "atlas.txt"
known = {ln.split("\t")[0] for ln in open(km) if ln.strip() and not ln.startswith("#")} if km.exists() else set()
chains = [ch for ch in chains if ch not in known]


class RangeFile(io.RawIOBase):  # seekable file backed by http range requests (zipfile only reads the bytes it needs)
    def __init__(self, url):
        self.url, self.pos = url, 0
        r = urllib.request.urlopen(urllib.request.Request(url, headers={"Range": "bytes=0-0"}))
        self.size = int(r.headers["Content-Range"].split("/")[1])

    def seek(self, off, whence=0): self.pos = {0: off, 1: self.pos + off, 2: self.size + off}[whence]; return self.pos
    def tell(self): return self.pos
    def seekable(self): return True

    def read(self, n=-1):
        n = self.size - self.pos if n < 0 else n
        if not n: return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(req).read(); self.pos += len(data); return data


def pull(chain):  # extract only <chain>_RMSF.tsv from the remote analysis zip
    dst = out / f"{chain}.tsv"
    if dst.exists(): return chain, "skip", ""
    err = ""
    for attempt in range(3):  # a failed fetch used to drop the chain silently
        try:
            zf = zipfile.ZipFile(RangeFile(api.format(chain)))
            name = next(n for n in zf.namelist() if n.endswith("_RMSF.tsv"))
            dst.write_bytes(zf.read(name)); return chain, "ok", ""
        except Exception as e:
            err = f"{type(e).__name__}: {e}"
            time.sleep(5 * 2 ** attempt)
    print(f"  {chain}: {err}", flush=True); return chain, "fail", err


results = []
with ThreadPoolExecutor(max_workers=8) as ex:
    for i, r in enumerate(ex.map(pull, chains), 1):
        results.append(r)
        if i % 100 == 0: print(f"  {i}/{len(chains)}", flush=True)
# what was fetched and what failed, so a partial download is visible after the fact
with open(d / "download_manifest.tsv", "w") as f:
    f.write("chain\tstatus\terror\n")
    for ch, st, err in results: f.write(f"{ch}\t{st}\t{err}\n")
failed = [ch for ch, st, _ in results if st == "fail"]
print(f"atlas: {len(list(out.glob('*.tsv')))} RMSF profiles -> {out} ({len(known)} known-missing left out)")
if failed:
    # a partial download changes the benchmark, and a saved split no longer matches it
    sys.exit(f"atlas: {len(failed)} chains failed after retries: {failed[:10]} -- re-run to resume, "
             "or add them to known_missing/atlas.txt with the reason")
PY
