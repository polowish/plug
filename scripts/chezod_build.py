#!/usr/bin/env python
# chezod -> data/chezod/chezod.tsv, the table plug's Chezod benchmark reads. one row per
# BMRB entry:
#
#   brmid      the BMRB entry id
#   split      test = CheZOD117; train/val = ADOPT's chezod_1159 (CheZOD1325 minus the
#              entries that duplicate or resemble the 117), val = the entries whose
#              sha1("chezod:" + brmid) is 0 mod 10 (~10%, fixed by id, not by draw order)
#   sequence   the entry's sequence, as the chemical shifts were assigned on it
#   zscore     per position of `sequence`: the CheZOD Z-score, "nan" where the entry has
#              no shifts for that residue (999.0 in the source)
#
#   python scripts/chezod_build.py [data/chezod]
import csv
import hashlib
import json
import sys
from pathlib import Path

MISSING = 999.0


def records(path: Path) -> list[dict]:
    """[{brmid, sequence, zscore}] from either of ADOPT's layouts: a list of records
    (117_dataset_raw.json) or a pandas column dict (chezod_1159.json), whose score
    field is spelt `z-score`."""
    d = json.load(open(path))
    if isinstance(d, dict):
        cols = list(d)
        d = [{k: d[k][i] for k in cols} for i in d[cols[0]]]
    out = []
    for r in d:
        z = r.get("zscore", r.get("z-score"))
        if len(z) != len(r["sequence"]):
            raise SystemExit(f"{path.name} {r['brmid']}: {len(z)} z-scores for {len(r['sequence'])} residues")
        out.append({"brmid": str(r["brmid"]), "sequence": r["sequence"].upper(), "zscore": z})
    return out


def is_val(brmid: str) -> bool:
    return int(hashlib.sha1(f"chezod:{brmid}".encode()).hexdigest()[:8], 16) % 10 == 0


def main() -> int:
    d = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "data" / "chezod")
    test = records(d / "117_dataset_raw.json")
    pool = records(d / "chezod_1159.json")
    clash = {r["brmid"] for r in test} & {r["brmid"] for r in pool}
    clash |= {r["brmid"] for r in pool if r["sequence"] in {t["sequence"] for t in test}}
    if clash:
        raise SystemExit(f"{len(clash)} train/val entries are also test entries: {sorted(clash)[:10]}")
    rows = [dict(r, split="test") for r in test]
    rows += [dict(r, split="val" if is_val(r["brmid"]) else "train") for r in pool]
    out = d / "chezod.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["brmid", "split", "sequence", "zscore"], delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({**r, "zscore": ",".join("nan" if x == MISSING else f"{x:g}" for x in r["zscore"])})
    for s in ("train", "val", "test"):
        rs = [r for r in rows if r["split"] == s]
        n = sum(len(r["sequence"]) for r in rs)
        k = sum(x != MISSING for r in rs for x in r["zscore"])
        print(f"{s:5s} {len(rs):5d} entries {n:7d} residues, {k} with a z-score ({k / n:.1%})")
    print(f"-> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
