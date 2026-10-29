#!/usr/bin/env python
# mdcath: per-frame descriptors + per-domain summaries out of one mdcath_dataset_<dom>.h5,
# so the ~0.2-2.5 GB trajectory file can be deleted afterwards (the whole dataset is 3.6 TB).
# one pass serves two users:
#   summaries (dyn_models items 2-3): per-residue rmsf per temperature/replica, ca contact
#     frequency and the dynamical cross-correlation matrix (dccm) per temperature
#   frames (the ensemble pilot, dyn_models wt/ensemble-pilot): per frame, a bit-packed ca
#     contact map, backbone phi/psi, and distances to 8 fixed nonlocal partners per residue
#
#   python scripts/mdcath_extract.py data/mdcath/raw/mdcath_dataset_1ss3A00.h5 --out data/mdcath
#
# format (also written into every output file's attrs, so a file explains itself):
#   frames/<dom>.h5    /<T>/<rep>/ frame (int32), time_ns (float32), contacts (uint8,
#                      frames x ceil(L(L-1)/2 / 8): np.packbits of the upper triangle i<j,
#                      row-major, ca-ca < CONTACT_A), phi, psi (float16 radians, nan where
#                      undefined), knn_dist (float16 angstrom, frames x L x K), rmsd, rg
#                      / knn_idx (int16 L x K), sequence, resid, attrs incl. knn_reference
#   summaries/<dom>.npz  sequence, temperatures, rmsf (T x rep x L, nm, mdcath's own),
#                      rmsf_ca (T x rep x L, angstrom, recomputed here as a check),
#                      contact_freq (T x L x L float16), dccm (T x L x L float16), n_frames
#
# knn partners are chosen ONCE per domain from the shared starting structure
# (pdbProteinAtoms: every temperature and replica starts from it), so a difference between
# temperatures compares the same residue pairs. contacts use ESMDynamic's 8 A ca cutoff.
import argparse
import json
import os
import sys
from pathlib import Path

import h5py
import numpy as np

CONTACT_A = 8.0
K = 8
MIN_SEP = 6          # |i - j| > 5
FRAME_NS = 1.0       # mdcath writes one frame per ns
THREE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
         "GLY": "G", "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H", "HIE": "H", "HID": "H",
         "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
         "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def _text(ds):
    v = ds[()]
    return v.decode() if isinstance(v, bytes) else str(v)


def backbone(g):
    """Atom indices of N, CA, C per residue, the sequence, residue ids, and the starting
    structure's ca coordinates, all from the domain's pdbProteinAtoms (atom order = coords)."""
    names, resids, resn, xyz, keys = [], [], [], [], []
    for line in _text(g["pdbProteinAtoms"]).splitlines():
        if line.startswith(("ATOM", "HETATM")):
            names.append(line[12:16].strip())
            resn.append(line[17:20].strip())
            resids.append(int(line[22:26]))
            xyz.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            # a residue is a run of atoms sharing number, insertion code and name: an
            # insertion (52, 52A) repeats the number, so the number alone merges residues
            keys.append((line[22:27], line[17:20]))
    n_atoms = int(g.attrs["numProteinAtoms"])
    if len(names) != n_atoms:
        raise SystemExit(f"pdbProteinAtoms has {len(names)} atoms, numProteinAtoms says {n_atoms}")
    if not np.array_equal(np.asarray(resids), g["resid"][()]):
        raise SystemExit("pdbProteinAtoms residue ids disagree with the resid dataset")
    order, idx, resnum = [], {}, []
    for a, (nm, r, rn, key) in enumerate(zip(names, resids, resn, keys)):
        if a == 0 or key != keys[a - 1]:
            order.append(len(order))
            idx[order[-1]] = {"resn": rn}
            resnum.append(r)
        idx[order[-1]][nm] = a
    for need in ("N", "CA", "C"):
        missing = [resnum[r] for r in order if need not in idx[r]]
        if missing:
            raise SystemExit(f"residues {missing[:5]} have no {need} atom")
    n = np.array([idx[r]["N"] for r in order])
    ca = np.array([idx[r]["CA"] for r in order])
    c = np.array([idx[r]["C"] for r in order])
    seq = "".join(THREE.get(idx[r]["resn"], "X") for r in order)
    start_ca = np.asarray(xyz, np.float32)[ca]
    return n, ca, c, seq, np.asarray(resnum), start_ca


def dihedral(p0, p1, p2, p3):
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1, axis=-1, keepdims=True)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    x = (v * w).sum(-1)
    y = (np.cross(b1, v) * w).sum(-1)
    return np.arctan2(y, x)


def phi_psi(X, n, ca, c):
    """(frames, L) phi and psi in radians; nan at the chain ends."""
    N, CA, C = X[:, n], X[:, ca], X[:, c]
    F, L = CA.shape[:2]
    phi = np.full((F, L), np.nan, np.float32)
    psi = np.full((F, L), np.nan, np.float32)
    if L > 1:
        phi[:, 1:] = dihedral(C[:, :-1], N[:, 1:], CA[:, 1:], C[:, 1:])
        psi[:, :-1] = dihedral(N[:, :-1], CA[:, :-1], C[:, :-1], N[:, 1:])
    return phi, psi


def superpose(CA, ref):
    """Every frame's ca (F, L, 3) rigidly fitted onto ref (L, 3) (Kabsch), so rmsf and dccm
    measure internal motion, whatever the trajectory's own frame alignment is."""
    A = CA - CA.mean(1, keepdims=True)
    B = ref - ref.mean(0)
    U, _, Vt = np.linalg.svd(np.einsum("fli,lj->fij", A, B))
    d = np.sign(np.linalg.det(np.einsum("fij,fjk->fik", U, Vt)))
    U[:, :, -1] *= d[:, None]
    R = np.einsum("fij,fjk->fik", U, Vt)
    return np.einsum("fli,fij->flj", A, R)


def knn_partners(start_ca, k=K, min_sep=MIN_SEP):
    """(L, k) indices of each residue's k nearest residues at |i - j| >= min_sep, by
    distance in the starting structure; -1 where a short chain has fewer candidates."""
    L = len(start_ca)
    d = np.linalg.norm(start_ca[:, None] - start_ca[None], axis=-1)
    sep = np.abs(np.arange(L)[:, None] - np.arange(L)[None])
    d = np.where(sep >= min_sep, d, np.inf)
    idx = np.argsort(d, axis=1, kind="stable")[:, :k]
    idx = np.where(np.isfinite(np.take_along_axis(d, idx, 1)), idx, -1)
    if idx.shape[1] < k:
        idx = np.pad(idx, ((0, 0), (0, k - idx.shape[1])), constant_values=-1)
    return idx.astype(np.int16)


def extract(path: Path, out: Path, keep_frames: bool = True) -> dict:
    f = h5py.File(path, "r")
    dom = next(iter(f.keys()))
    g = f[dom]
    n, ca, c, seq, resid, start_ca = backbone(g)
    L = len(seq)
    iu = np.triu_indices(L, 1)
    knn = knn_partners(start_ca)
    temps = sorted((k for k in g.keys() if k.isdigit()), key=int)
    reps = sorted({r for T in temps for r in g[T].keys()}, key=int)

    rmsf = np.full((len(temps), len(reps), L), np.nan, np.float32)
    rmsf_ca = np.full_like(rmsf, np.nan)
    freq = np.zeros((len(temps), L, L), np.float32)
    dccm = np.zeros((len(temps), L, L), np.float32)
    n_frames = np.zeros((len(temps), len(reps)), np.int32)

    fr = None
    if keep_frames:
        (out / "frames").mkdir(parents=True, exist_ok=True)
        fr = h5py.File(out / "frames" / f"{dom}.h5.part", "w")
        fr.attrs.update({"domain": dom, "format": "mdcath_frames_v1", "contact_cutoff_A": CONTACT_A,
                         "contacts": "np.packbits of the upper triangle (i<j, row-major) of ca-ca < cutoff",
                         "knn_k": K, "knn_min_separation": MIN_SEP,
                         "knn_reference": "starting structure (pdbProteinAtoms), shared by every temperature and replica",
                         "frame_interval_ns": FRAME_NS, "source": path.name})
        fr["sequence"] = seq
        fr["resid"] = resid.astype(np.int32)
        fr["knn_idx"] = knn
    safe = np.where(knn >= 0, knn, 0)

    for ti, T in enumerate(temps):
        cnt = 0
        for ri, r in enumerate(reps):
            if r not in g[T]:
                continue
            rg_ = g[T][r]
            X = rg_["coords"][()]                       # (F, atoms, 3), angstrom
            F = X.shape[0]
            n_frames[ti, ri] = F
            rmsf[ti, ri] = rg_["rmsf"][()]
            CA = X[:, ca].astype(np.float64)
            fit = superpose(CA, start_ca.astype(np.float64))
            dev = fit - fit.mean(0)
            rmsf_ca[ti, ri] = np.sqrt((dev ** 2).sum(-1).mean(0))
            # dccm per replica, fluctuations about that replica's mean after fitting every
            # frame onto the starting structure; averaged over replicas below
            cov = np.einsum("fid,fjd->ij", dev, dev) / F
            sd = np.sqrt(np.diag(cov))
            dccm[ti] += cov / np.outer(sd, sd)
            D = np.linalg.norm(CA[:, :, None] - CA[:, None, :], axis=-1)
            contact = D < CONTACT_A
            freq[ti] += contact.sum(0)
            cnt += F
            if fr is not None:
                grp = fr.create_group(f"{T}/{r}")
                grp["frame"] = np.arange(F, dtype=np.int32)
                grp["time_ns"] = (np.arange(F) * FRAME_NS).astype(np.float32)
                grp.create_dataset("contacts", data=np.packbits(contact[:, iu[0], iu[1]], axis=1),
                                   compression="gzip", compression_opts=4)
                phi, psi = phi_psi(X, n, ca, c)
                grp.create_dataset("phi", data=phi.astype(np.float16), compression="gzip", compression_opts=4)
                grp.create_dataset("psi", data=psi.astype(np.float16), compression="gzip", compression_opts=4)
                kd = np.take_along_axis(D, np.broadcast_to(safe, (F, L, K)), axis=2)
                kd = np.where(knn[None] >= 0, kd, np.nan)
                grp.create_dataset("knn_dist", data=kd.astype(np.float16), compression="gzip", compression_opts=4)
                grp["rmsd"] = rg_["rmsd"][()].astype(np.float32)
                grp["rg"] = rg_["gyrationRadius"][()].astype(np.float32)
        nrep = int((n_frames[ti] > 0).sum())
        freq[ti] /= max(cnt, 1)
        dccm[ti] /= max(nrep, 1)
    f.close()

    (out / "summaries").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "summaries" / f"{dom}.npz", sequence=np.array(seq), resid=resid,
                        temperatures=np.array([int(t) for t in temps]), replicas=np.array([int(r) for r in reps]),
                        rmsf=rmsf, rmsf_ca=rmsf_ca, contact_freq=freq.astype(np.float16),
                        dccm=dccm.astype(np.float16), n_frames=n_frames,
                        contact_cutoff_A=CONTACT_A)
    if fr is not None:
        fr.close()
        (out / "frames" / f"{dom}.h5.part").rename(out / "frames" / f"{dom}.h5")
    return {"domain": dom, "L": L, "frames": int(n_frames.sum()),
            "rmsf_check_max_abs_A": float(np.nanmax(np.abs(rmsf * 10 - rmsf_ca)))}


def _one(job):
    path, out, keep_frames = job
    try:
        return extract(path, out, keep_frames=keep_frames), None
    except BaseException as e:                              # recorded; the worker carries on
        return {"domain": path.stem}, f"{type(e).__name__}: {e}"


def watch(staging: Path, out: Path, workers: int, keep_frames: bool) -> int:
    """Consume <staging>/*.h5 as mdcath_fetch.py delivers them: claim (atomic rename into
    claimed/), extract, append to <out>/done.txt, delete. A failed domain moves to failed/
    with its error in <out>/extract_failed.txt. Stops when FETCH_DONE exists and nothing is
    left. Claimed files from a killed run are put back first, so a rerun resumes."""
    import time
    from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait

    claimed, failed = staging / "claimed", staging / "failed"
    for d in (claimed, failed, out):
        d.mkdir(parents=True, exist_ok=True)
    for f in claimed.glob("*.h5"):
        os.replace(f, staging / f.name)
    running = {}
    n_ok = n_bad = 0
    with ProcessPoolExecutor(max_workers=workers) as pool:
        while True:
            while len(running) < workers:
                todo = sorted(staging.glob("*.h5"))
                if not todo:
                    break
                src = todo[0]
                dst = claimed / src.name
                try:
                    os.replace(src, dst)
                except FileNotFoundError:
                    continue
                running[pool.submit(_one, (dst, out, keep_frames))] = dst
            if not running:
                if (staging / "FETCH_DONE").exists() and not list(staging.glob("*.h5")):
                    break
                time.sleep(15)
                continue
            done, _ = wait(list(running), timeout=60, return_when=FIRST_COMPLETED)
            for fut in done:
                dst = running.pop(fut)
                info, err = fut.result()
                if err is None:
                    with open(out / "done.txt", "a") as fh:
                        fh.write(info["domain"] + "\n")
                    dst.unlink()
                    n_ok += 1
                    print(json.dumps(info), flush=True)
                else:
                    os.replace(dst, failed / dst.name)
                    with open(out / "extract_failed.txt", "a") as fh:
                        fh.write(f"{dst.stem}\t{err}\n")
                    n_bad += 1
                    print(f"FAILED {dst.stem}: {err}", flush=True)
    print(f"done: {n_ok} extracted, {n_bad} failed", flush=True)
    return 1 if n_bad else 0


def main() -> int:
    p = argparse.ArgumentParser(description="mdcath h5 -> per-frame descriptors + summaries")
    p.add_argument("h5", nargs="*", type=Path)
    p.add_argument("--out", type=Path, default=Path(__file__).resolve().parent.parent / "data" / "mdcath")
    p.add_argument("--no-frames", action="store_true", help="summaries only")
    p.add_argument("--watch", type=Path, metavar="STAGING",
                   help="consume what mdcath_fetch.py stages, deleting each file once extracted")
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    if args.watch:
        return watch(args.watch, args.out, args.workers, not args.no_frames)
    for h in args.h5:
        info = extract(h, args.out, keep_frames=not args.no_frames)
        print(json.dumps(info), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
