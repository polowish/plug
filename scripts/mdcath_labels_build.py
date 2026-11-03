#!/usr/bin/env python
# mdcath extraction summaries -> data/mdcath/mdcath_labels.tsv, the table plug's Mdcath
# benchmark reads (dyn_models docs/plans/dynamics_benchmarks.md, items 2 and 3). one row
# per domain:
#
#   domain, split, sequence, L        split = ESMDynamic's mdcath train/val/test lists,
#                                     never re-dealt (so ESMDynamic is scored only on
#                                     domains it did not train on)
#   n_folded_320, n_ok_379, n_ok_413  replicas behind the labels below, from the shared
#                                     stability table (disco-lab-5f, mdcath_stability.py):
#                                     "folded" at 320 K; "folded" or "partial" (= not
#                                     unfolded) at the hot temperature
#   rmsf_320                          per-residue rmsf (A, mdcath's own) averaged over the
#                                     replicas folded at 320 K; needs >= 3 of them
#   drmsf_379, drmsf_413              thermal response: mean rmsf over the hot replicas
#                                     that did not unfold minus rmsf_320; needs >= 3 of
#                                     each. empty where a domain is not eligible
#   drmsf_379_resid                   drmsf_379 minus its per-domain least-squares fit on
#                                     rmsf_320: the part of the response that baseline
#                                     flexibility does not predict
#   pairs_contact                     up to 96 residue pairs (|i-j| >= 6) in contact
#                                     (ca-ca < 8 A) in >= 5% of 320 K frames, with that
#                                     frequency: "i,j,v;..." 1-based
#   pairs_dccm                        up to 96 pairs, 32 per separation band
#                                     ([6,12), [12,24), [24,inf)), with the 320 K ca
#                                     dynamic cross-correlation
#   the pair columns only for domains folded in all five 320 K replicas, since the
#   summaries average contacts and dccm over every replica at a temperature.
#
# per replica, for noise ceilings and replica-level uncertainty (the columns above stay the
# probe targets, unchanged):
#   calls_<T>, T in 320/379/413       each replica's stability call, in replica order
#   rmsf_reps_<T>                     each replica's per-residue rmsf, "|" between
#                                     replicas, whatever its call (filter with calls_<T>)
#   pairs_contact_reps                the pairs_contact pairs with each 320 K replica's
#                                     frequency, from the per-frame contact maps
#                                     (frames/<dom>.h5): "i,j,v0/v1/v2/v3/v4;..." 1-based
#   there is no per-replica dccm: the extractor averaged it over replicas and kept no
#   coordinates, so it needs the raw trajectories again.
#
# pairs are drawn with a generator seeded by the domain id, so a rebuild draws the same
# pairs. a sidecar mdcath_labels.meta.json records the sources and rules.
#
#   python scripts/mdcath_labels_build.py --mdcath ~/fac/data/mdcath \
#       [--frames ~/fac/data/mdcath/frames]   (needs h5py: ~/scratch/venvs/mdcath)
#       --splits ~/scratch/github/esmdynamic-investigation/investigation/data/mdcath \
#       [--out data/mdcath]
import argparse
import csv
import json
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np

HOT = (379, 413)
MIN_REPLICAS = 3
N_CONTACT, CONTACT_MIN = 96, 0.05
SEP_BANDS, PER_BAND = ((6, 12), (12, 24), (24, 10**6)), 32
MIN_SEP = 6


def calls(path: Path) -> dict:
    """{domain: {(temperature, replica): call}} from the stability table."""
    out = defaultdict(dict)
    with open(path) as f:
        for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"):
            out[r["domain"]][(int(r["temperature"]), int(r["replica"]))] = r["call"]
    return out


def fmt(v) -> str:
    return ",".join(f"{x:.3f}" for x in v)


def fmt_pairs(pairs) -> str:
    return ";".join(f"{i + 1},{j + 1},{v:.3f}" for i, j, v in pairs)


def replica_contacts(frames: Path, dom: str, L: int, pairs, reps) -> list | None:
    """(len(pairs), n_reps) 320 K contact frequency per replica, read off the bit-packed
    per-frame maps (np.packbits of the i<j upper triangle, row-major); None without them."""
    import h5py

    path = frames / f"{dom}.h5"
    if not path.exists():
        return None
    k = np.array([i * (2 * L - i - 1) // 2 + (j - i - 1) for i, j in pairs], np.int64)
    out = np.full((len(pairs), len(reps)), np.nan)
    with h5py.File(path, "r") as f:
        for ri, r in enumerate(reps):
            key = f"320/{r}/contacts"
            if key not in f:
                continue
            bits = f[key][()]                                   # (F, ceil(L(L-1)/2 / 8))
            out[:, ri] = ((bits[:, k >> 3] >> (7 - (k & 7))) & 1).mean(0)
    return out


def draw(rng, idx, n):
    return idx if len(idx) <= n else idx[np.sort(rng.choice(len(idx), n, replace=False))]


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--mdcath", type=Path, required=True, help="dir with summaries/ and stability.tsv")
    p.add_argument("--splits", type=Path, required=True, help="ESMDynamic's {train,val,test}.csv")
    p.add_argument("--frames", type=Path, default=None, help="per-frame files (default <mdcath>/frames)")
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "data" / "mdcath")
    a = p.parse_args()

    frames = a.frames or a.mdcath / "frames"
    split = {d: s for s in ("train", "val", "test") for d in (a.splits / f"{s}.csv").read_text().split()}
    stab = calls(a.mdcath / "stability.tsv")
    a.out.mkdir(parents=True, exist_ok=True)
    rows, counts = [], defaultdict(int)
    for path in sorted((a.mdcath / "summaries").glob("*.npz")):
        dom = path.stem
        if dom not in split:
            counts["not_in_esmdynamic_split"] += 1; continue
        z = np.load(path)
        temps = [int(t) for t in z["temperatures"]]
        reps = [int(r) for r in z["replicas"]]
        seq = str(z["sequence"])
        L = len(seq)
        c = stab.get(dom, {})

        def ok(T, allowed):
            return [ri for ri, r in enumerate(reps) if c.get((T, r)) in allowed]

        cold = ok(320, {"folded"})
        row = {"domain": dom, "split": split[dom], "sequence": seq, "L": L,
               "n_folded_320": len(cold)}
        rmsf = z["rmsf"].astype(np.float64)                       # (T, R, L)
        for T in (320,) + HOT:
            ti = temps.index(T)
            row[f"calls_{T}"] = ",".join(c.get((T, r), "na") for r in reps)
            row[f"rmsf_reps_{T}"] = "|".join(fmt(rmsf[ti, ri]) for ri in range(len(reps)))
        base = rmsf[temps.index(320), cold].mean(0) if len(cold) >= MIN_REPLICAS else None
        row["rmsf_320"] = fmt(base) if base is not None else ""
        for T in HOT:
            hot = ok(T, {"folded", "partial"})
            row[f"n_ok_{T}"] = len(hot)
            if base is not None and len(hot) >= MIN_REPLICAS:
                row[f"drmsf_{T}"] = fmt(rmsf[temps.index(T), hot].mean(0) - base)
                counts[f"drmsf_{T}"] += 1
            else:
                row[f"drmsf_{T}"] = ""
        if row["drmsf_379"]:
            d = np.array([float(x) for x in row["drmsf_379"].split(",")])
            A = np.c_[np.ones(L), base]
            coef, *_ = np.linalg.lstsq(A, d, rcond=None)
            row["drmsf_379_resid"] = fmt(d - A @ coef)
        else:
            row["drmsf_379_resid"] = ""

        row["pairs_contact"] = row["pairs_dccm"] = row["pairs_contact_reps"] = ""
        if len(cold) == len(reps) == 5:
            t0 = temps.index(320)
            freq = z["contact_freq"][t0].astype(np.float64)
            dccm = z["dccm"][t0].astype(np.float64)
            rng = np.random.default_rng(zlib.crc32(dom.encode()))
            i, j = np.triu_indices(L, MIN_SEP)
            con = np.flatnonzero(freq[i, j] >= CONTACT_MIN)
            pick = draw(rng, con, N_CONTACT)
            row["pairs_contact"] = fmt_pairs((i[k], j[k], freq[i[k], j[k]]) for k in pick)
            per = replica_contacts(frames, dom, L, [(i[k], j[k]) for k in pick], reps)
            if per is not None:
                # the summary pools frames over replicas: the frame-weighted replica mean
                # must give it back (to float16 rounding), or the bit indexing is wrong
                w = z["n_frames"][t0].astype(np.float64)
                pooled = (per * w).sum(1) / w.sum()
                if np.abs(pooled - freq[i[pick], j[pick]]).max(initial=0) > 2e-3:
                    raise SystemExit(f"{dom}: per-replica contacts do not average to the summary's")
                row["pairs_contact_reps"] = ";".join(
                    f"{i[k] + 1},{j[k] + 1}," + "/".join(f"{v:.3f}" for v in per[n])
                    for n, k in enumerate(pick))
                counts["pairs_contact_reps"] += 1
            sep = j - i
            picks = [draw(rng, np.flatnonzero((sep >= lo) & (sep < hi)), PER_BAND) for lo, hi in SEP_BANDS]
            pick = np.sort(np.concatenate(picks))
            row["pairs_dccm"] = fmt_pairs((i[k], j[k], dccm[i[k], j[k]]) for k in pick if np.isfinite(dccm[i[k], j[k]]))
            counts["pairs"] += 1
        counts["rmsf_320"] += base is not None
        rows.append(row)

    out = a.out / "mdcath_labels.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    meta = {"sources": {"summaries": str(a.mdcath / "summaries"), "stability": str(a.mdcath / "stability.tsv"),
                        "frames": str(frames),
                        "splits": str(a.splits)},
            "rules": {"min_replicas": MIN_REPLICAS, "cold": "folded at 320 K", "hot": "folded or partial",
                      "hot_temperatures": HOT, "contact_min_frequency": CONTACT_MIN, "n_contact_pairs": N_CONTACT,
                      "dccm_bands": SEP_BANDS, "dccm_per_band": PER_BAND, "min_separation": MIN_SEP,
                      "pairs_require": "all 5 replicas folded at 320 K",
                      "per_replica": "calls_<T>, rmsf_reps_<T> (T = 320, 379, 413), pairs_contact_reps "
                                     "(from frames); no per-replica dccm"},
            "counts": dict(counts), "n_rows": len(rows)}
    (a.out / "mdcath_labels.meta.json").write_text(json.dumps(meta, indent=1) + "\n")
    by = defaultdict(lambda: defaultdict(int))
    for r in rows:
        for k in ("rmsf_320", "drmsf_379", "drmsf_413", "pairs_dccm"):
            by[k][r["split"]] += bool(r[k])
    print(f"{out}: {len(rows)} domains; " + "; ".join(f"{k} {dict(v)}" for k, v in by.items()))
    print("  " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
