# chezod: per-residue order/disorder measured by NMR -- the CheZOD Z-score, from how far a
# residue's backbone chemical shifts sit from random-coil values (Nielsen & Mulder; high =
# ordered, below ~8 = disordered). one item per BMRB entry, with a score per residue of
# the sequence the shifts were assigned on. ships its own split: CheZOD117 is the test
# set; ADOPT's 1159 (CheZOD1325 minus the entries that duplicate or resemble the 117)
# is train, with ~10% carved off as val by a hash of the entry id.
# the table is built by scripts/chezod_build.py (data/download_chezod.sh says how).
import csv

from . import Benchmark
from .. import config as c


class Chezod(Benchmark):
    name = "chezod"
    # every item is about one BMRB entry
    key = "brmid"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        with open(c.data / "chezod" / "chezod.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                yield {"sequence": r["sequence"],
                       "brmid": r["brmid"],
                       # one value per residue of the sequence; nan = no chemical shifts
                       "zscore": [float(x) for x in r["zscore"].split(",")],
                       "split": r["split"]}
