# cryptobench: cryptic binding sites -- residues that line a ligand pocket in a holo
# structure but not in the apo structure of the same protein (osf.io/pz4a9, MIT). one item
# per apo chain, with the site counted along that chain's full sequence. ships its own
# split: a test set and four train folds; here folds 0-2 train, fold 3 validates.
# the table is built from the mmcif files by scripts/cryptobench_build.py
# (data/download_cryptobench.sh says how).
#
# apoholo is the same chains with a different target: each residue's c-alpha displacement
# between the apo structure and its main holo structure (the one with the largest pocket
# rmsd), after a fit over every residue modelled in both. nan where either structure has
# no c-alpha for the residue.
import csv
import sys

from . import Benchmark
from .. import config as c

csv.field_size_limit(sys.maxsize)
SPLIT = {"train-0": "train", "train-1": "train", "train-2": "train", "train-3": "val", "test": "test"}


def _rows():
    path = c.data / "cryptobench" / "cryptobench.tsv"
    with open(path, newline="") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            split = SPLIT.get(r["fold"])
            if split is None: continue
            yield r, split


class Cryptobench(Benchmark):
    name = "cryptobench"
    # every item is about one apo chain
    key = "pdb_chain"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        for r, split in _rows():
            yield {"sequence": r["sequence"],
                   "pdb_chain": f"{r['apo_pdb']}_{r['chain']}",
                   "uniprot": r["uniprot"],
                   # 1-based positions along the sequence above
                   "cryptic_site": [int(p) for p in r["cryptic_site"].split(",") if p],
                   "split": split}


class Apoholo(Benchmark):
    name = "apoholo"
    key = "pdb_chain"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        for r, split in _rows():
            disp = [float(x) for x in r["displacement"].split(",")]
            if all(x != x for x in disp): continue          # no main holo structure to compare
            yield {"sequence": r["sequence"],
                   "pdb_chain": f"{r['apo_pdb']}_{r['chain']}",
                   "uniprot": r["uniprot"],
                   "main_holo": r["main_holo"],
                   # one value per residue of the sequence; nan = not modelled in both
                   "displacement": disp,
                   "split": split}
