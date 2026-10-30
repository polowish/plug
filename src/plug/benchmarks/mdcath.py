# mdcath: CATH domains simulated at 320-450 K, five replicas each (Mirarchi et al. 2024,
# huggingface compsciencelab/mdCATH, CC BY 4.0). one item per domain, read from the label
# table scripts/mdcath_labels_build.py makes from the extraction summaries
# (scripts/mdcath_fetch.py + mdcath_extract.py). ships its split as ESMDynamic's
# train/val/test domain lists, so a model trained on mdcath is scored on held-out domains.
#
# per-residue fields are lists as long as the sequence, or None where the domain is not
# eligible for that label (too few folded replicas; the table's header comment has the
# rules). pair fields are lists of (i, j, value) with 1-based positions, or None.
import csv
import sys

from . import Benchmark
from .. import config as c

csv.field_size_limit(sys.maxsize)
RESIDUE = ("rmsf_320", "drmsf_379", "drmsf_413", "drmsf_379_resid")
PAIRS = ("pairs_contact", "pairs_dccm")


def _floats(s):
    return [float(x) for x in s.split(",")] if s else None


def _pairs(s):
    if not s: return None
    out = []
    for t in s.split(";"):
        i, j, v = t.split(",")
        out.append((int(i), int(j), float(v)))
    return out


class Mdcath(Benchmark):
    name = "mdcath"
    # every item is about one cath domain
    key = "domain"
    published = ("train", "val", "test")

    @classmethod
    def rows(cls):
        with open(c.data / "mdcath" / "mdcath_labels.tsv", newline="") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                yield {"sequence": r["sequence"], "domain": r["domain"], "split": r["split"],
                       **{k: _floats(r[k]) for k in RESIDUE},
                       **{k: _pairs(r[k]) for k in PAIRS}}
