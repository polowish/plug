#!/usr/bin/env python
# cryptobench -> data/cryptobench/cryptobench.tsv, the table plug's Cryptobench and Apoholo
# benchmarks read. one row per apo CHAIN (an apo entry can contribute several chains):
#
#   apo_pdb, chain, uniprot, fold      the apo structure and its cryptobench fold
#                                      (test, train-0..3)
#   sequence                           the chain's full sequence from the apo mmcif
#                                      (pdbx_poly_seq_scheme), modelled or not
#   cryptic_site                       1-based positions in `sequence` of every residue in
#                                      any of the chain's cryptic pockets (union over its
#                                      holo structures)
#   n_site_unmapped                    pocket residues whose author number is not in the
#                                      chain's sequence scheme (dropped, never guessed)
#   main_holo, main_prmsd              the holo structure cryptobench marks as the main one
#                                      (highest pocket rmsd)
#   displacement                       per position of `sequence`: c-alpha displacement (A)
#                                      between the apo chain and the main holo chain after a
#                                      kabsch fit over every residue modelled in both;
#                                      "nan" where either side has no c-alpha
#
# residues are matched apo<->holo by aligning the two chains' full sequences, not by author
# number, which the two entries need not share.
#
#   python scripts/cryptobench_build.py [data/cryptobench]
import csv
import json
import sys
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1


def chain_scheme(cif: dict, chain: str):
    """([(one-letter, auth key or None)] in sequence order, author chain id) for `chain`,
    from pdbx_poly_seq_scheme; the key is the author residue number + insertion code.
    `chain` is matched as an author chain id first, then as a label (asym) id, since an
    entry can name a chain either way."""
    s = "_pdbx_poly_seq_scheme."
    strands, asyms = cif[s + "pdb_strand_id"], cif[s + "asym_id"]
    use = strands if chain in strands else asyms
    auth_chain = chain if use is strands else next((st for st, a in zip(strands, asyms) if a == chain), chain)
    rows = zip(use, cif[s + "mon_id"], cif[s + "pdb_seq_num"],
               cif[s + "auth_seq_num"], cif.get(s + "pdb_ins_code", ["."] * len(cif[s + "mon_id"])))
    out = []
    for strand, mon, pdb_num, auth_num, ins in rows:
        if strand != chain:
            continue
        aa = seq1(mon, custom_map={"MSE": "M"}) or "X"
        aa = aa if len(aa) == 1 else "X"
        observed = auth_num not in ("?", ".")
        key = (pdb_num + ("" if ins in ("?", ".") else ins)) if observed else None
        out.append((aa, key))
    return out, auth_chain


def ca_coords(cif: dict, chain: str) -> dict:
    """{auth key: c-alpha xyz} for the chain's first model."""
    a = "_atom_site."
    model = cif.get(a + "pdbx_PDB_model_num")
    first = model[0] if model else None
    ins = cif.get(a + "pdbx_PDB_ins_code", ["?"] * len(cif[a + "label_atom_id"]))
    out = {}
    for i, (name, ch, num) in enumerate(zip(cif[a + "label_atom_id"], cif[a + "auth_asym_id"], cif[a + "auth_seq_id"])):
        if name != "CA" or ch != chain or (model and model[i] != first):
            continue
        alt = cif[a + "label_alt_id"][i]
        key = num + ("" if ins[i] in ("?", ".") else ins[i])
        if key in out and alt not in (".", "?", "A"):
            continue                                   # keep the first altloc
        out[key] = (float(cif[a + "Cartn_x"][i]), float(cif[a + "Cartn_y"][i]), float(cif[a + "Cartn_z"][i]))
    return out


def kabsch_rmsd_fit(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """P fitted onto Q (both n x 3); returns the fitted P."""
    p0, q0 = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - p0).T @ (Q - q0))
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    return (P - p0) @ U @ D @ Vt + q0


def main() -> int:
    d = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "data" / "cryptobench")
    data = json.load(open(d / "dataset.json"))
    folds = json.load(open(d / "folds.json"))
    fold_of = {pdb: f for f, ids in folds.items() for pdb in ids}
    cifs = {p.stem.lower(): p for p in (d / "cif-files").rglob("*.cif")}
    cache: dict = {}

    def cif_of(pdb):
        if pdb not in cache:
            cache.clear()                              # one apo + one holo at a time
            cache[pdb] = MMCIF2Dict(str(cifs[pdb.lower()]))
        return cache[pdb]

    aligner = PairwiseAligner(mode="global", match_score=2, mismatch_score=-1,
                              open_gap_score=-5, extend_gap_score=-0.5)
    rows, skipped, no_holo_chain = [], [], []
    for apo, holos in sorted(data.items()):
        for chain in sorted({h["apo_chain"] for h in holos}):
            hs = [h for h in holos if h["apo_chain"] == chain]
            if apo.lower() not in cifs:
                skipped.append((apo, chain, "no apo mmcif")); continue
            apo_cif = cif_of(apo)
            scheme, apo_auth = chain_scheme(apo_cif, chain)
            if not scheme:
                skipped.append((apo, chain, "chain not in pdbx_poly_seq_scheme")); continue
            seq = "".join(a for a, _ in scheme)
            pos_of = {k: i for i, (_, k) in enumerate(scheme) if k is not None}
            site, unmapped = set(), 0
            for h in hs:
                for sel in h["apo_pocket_selection"]:
                    ch, _, key = sel.partition("_")
                    if ch in (chain, apo_auth) and key in pos_of:
                        site.add(pos_of[key] + 1)
                    else:
                        unmapped += 1
            apo_ca = ca_coords(apo_cif, apo_auth)
            main = max(hs, key=lambda h: (h.get("is_main_holo_structure", False), h.get("pRMSD", 0)))
            disp = [float("nan")] * len(seq)
            hp = main["holo_pdb_id"].lower()
            if hp in cifs:
                hscheme, holo_auth = chain_scheme(cif_of(hp), main["holo_chain"])
                holo_ca = ca_coords(cif_of(hp), holo_auth)
                hseq = "".join(a for a, _ in hscheme)
                if not hseq:
                    no_holo_chain.append(f"{apo}_{chain} -> {main['holo_pdb_id']}_{main['holo_chain']}")
                    hp = None
            if hp in cifs:
                aln = aligner.align(seq, hseq)[0]
                pairs = []
                for (a0, a1), (b0, b1) in zip(*aln.aligned):
                    for i, j in zip(range(a0, a1), range(b0, b1)):
                        ka, kb = scheme[i][1], hscheme[j][1]
                        if ka in apo_ca and kb in holo_ca:
                            pairs.append((i, apo_ca[ka], holo_ca[kb]))
                if len(pairs) >= 10:
                    P = np.array([p[1] for p in pairs]); Q = np.array([p[2] for p in pairs])
                    dist = np.linalg.norm(kabsch_rmsd_fit(P, Q) - Q, axis=1)
                    for (i, _, _), x in zip(pairs, dist):
                        disp[i] = float(x)
            rows.append({"apo_pdb": apo, "chain": chain, "uniprot": hs[0].get("uniprot_id", ""),
                         "fold": fold_of.get(apo, ""), "sequence": seq,
                         "cryptic_site": ",".join(map(str, sorted(site))),
                         "n_site_unmapped": unmapped, "main_holo": f"{main['holo_pdb_id']}_{main['holo_chain']}",
                         "main_prmsd": main.get("pRMSD", ""),
                         "displacement": ",".join("nan" if x != x else f"{x:.2f}" for x in disp)})
    out = d / "cryptobench.tsv"
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    n_disp = sum(any(x != "nan" for x in r["displacement"].split(",")) for r in rows)
    print(f"{out}: {len(rows)} apo chains, {sum(bool(r['cryptic_site']) for r in rows)} with a site, "
          f"{n_disp} with displacement, {sum(r['n_site_unmapped'] for r in rows)} pocket residues unmapped; "
          f"{len(skipped)} skipped")
    for s in skipped[:20]:
        print("  skipped", *s)
    if no_holo_chain:
        print(f"  {len(no_holo_chain)} main holo chains missing from their mmcif's sequence scheme "
              f"(no displacement): {', '.join(no_holo_chain[:10])}{' ...' if len(no_holo_chain) > 10 else ''}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
