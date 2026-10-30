# one benchmark = a list of items. each benchmark file only says how to read its own
# rows; this class does the rest: which split to use, and which one to hand back.
from torch.utils.data import Dataset
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from datetime import date
from pathlib import Path

from .. import config as c
from .. import homology

# items each benchmark deliberately leaves out: known_missing/<name>.txt, "key<TAB>why" per
# line. dropped on read, and not counted as failures by the download scripts, so what a
# benchmark contains is decided there on purpose, not by which downloads happened to succeed.
KNOWN_MISSING = Path(__file__).resolve().parents[3] / "known_missing"


class Benchmark(Dataset):
    # what this benchmark is called, also the name of its saved split file
    name = ""
    # the field that names the protein an item is about
    # (two fields, as a tuple, when an item is about a pair of proteins)
    key = "sequence"
    # the field holding that protein's sequence (also two, for pairs)
    seq = "sequence"
    # the splits this benchmark came with, if it came with any
    published = ()
    # worked out the first time a homology split is asked for, then kept here
    _split = None

    # each benchmark writes this one and nothing else: yield one dict per item.
    # if the benchmark came with splits, give every item the split it belongs to.
    @classmethod
    def rows(cls):
        raise NotImplementedError

    # every sequence in this benchmark, for the searches that keep the training set clean.
    # a benchmark holding sequences its items don't carry adds them on top of these.
    def sequences(self):
        for it in self.items:
            for f in self.seqfields(): yield it[f]

    # split is the split you want: "train", "val", "test", or "all".
    # how is where that split comes from:
    #   published — the split the benchmark came with (if any)
    #   homology  — a split worked out here, so no two splits share similarity
    #   whole     — no split at all, the whole benchmark is one test set (for self/unsurpervised)
    def __init__(self, split="test", how=None):
        # if not specified, use the benchmark's own split when it has one or the whole thing
        how = how or ("published" if self.published else "whole")
        # read every item first — the split has to be decided over all of them
        items, self.n_excluded = self.filtered(self.rows())
        self.n_unsplit = 0

        if how == "published":
            # if published requested but there is no published throw error
            if not self.published:
                raise ValueError(f"{self.name} ships no published split — "
                                 "use how='homology' to create a split or how='whole'"
                                 " to use data as one big test set")
            # if published requested but published doesnt have split throw error
            if split not in self.published + ("all",):
                raise ValueError(f"{self.name} publishes {'/'.join(self.published)}"
                                 f", {split!r} is not allowed")

        # label every item from the homology split, working it out if it isn't cached yet
        elif how == "homology":
            split_map = self.split_map()
            self.check_split(items, split_map)
            # items whose protein the split does not know (only reachable with
            # PLUG_ALLOW_SPLIT_MISMATCH=1): left out and counted, never a KeyError
            self.n_unsplit = sum(any(it[k] not in split_map for k in self.keyfields()) for it in items)
            items = self.stamp(items, split_map)

        # no splits, so everything counts as test
        elif how == "whole":
            for it in items: it["split"] = "test"

        else:
            raise ValueError(f"how: 'published', 'homology' or 'whole', not {how}")

        # keep only the splits that were asked for
        self.items = items if split == "all" else [it for it in items if it["split"] == split]

    # the number of items in this benchmark
    def __len__(self): return len(self.items)

    # get an item by index
    def __getitem__(self, i): return self.items[i]

    # key and seq are one field for most benchmarks and two for pairs. these always
    # give back a tuple, so nothing further down has to check which it is.
    @classmethod
    def keyfields(cls): return (cls.key,) if isinstance(cls.key, str) else tuple(cls.key)

    @classmethod
    def seqfields(cls): return (cls.seq,) if isinstance(cls.seq, str) else tuple(cls.seq)

    # the items this benchmark leaves out on purpose, {key: why}
    @classmethod
    def excluded(cls):
        f = KNOWN_MISSING / f"{cls.name}.txt"
        if not f.exists(): return {}
        out = {}
        for line in open(f):
            line = line.rstrip("\n")
            if line and not line.startswith("#"):
                k, _, why = line.partition("\t")
                out[k] = why
        return out

    # rows minus the excluded items -> (items, how many were left out)
    @classmethod
    def filtered(cls, rows):
        excl = cls.excluded()
        items = list(rows)
        if not excl: return items, 0
        kept = [it for it in items if not any(it[k] in excl for k in cls.keyfields())]
        return kept, len(items) - len(kept)

    # {protein: sha1 of its sequence} and one sha256 over all of it: what a split was built from
    @classmethod
    def fingerprint(cls, items):
        seqs = {}
        for it in items:
            for k, s in zip(cls.keyfields(), cls.seqfields()):
                seqs.setdefault(it[k], it[s])
        per = {k: hashlib.sha1(seqs[k].encode()).hexdigest()[:16] for k in sorted(seqs)}
        h = hashlib.sha256("".join(f"{k}\t{v}\n" for k, v in per.items()).encode()).hexdigest()
        return h, per

    # a saved split is only valid for exactly the data it was dealt from. compare the data
    # being loaded with the split (and with its fingerprint sidecar, when there is one) and
    # stop with the difference rather than quietly use a split built for other data.
    # PLUG_ALLOW_SPLIT_MISMATCH=1 downgrades this to a warning (unknown proteins are then
    # left out and counted in self.n_unsplit).
    @classmethod
    def check_split(cls, items, split_map):
        fp, per = cls.fingerprint(items)
        side = sidecar(cls.split_file())
        added = sorted(set(per) - set(split_map))
        removed = sorted(set(split_map) - set(per))
        changed = []
        if side.exists():
            meta = json.loads(side.read_text())
            if meta.get("fingerprint") == fp: return
            old = meta.get("proteins", {})
            changed = sorted(k for k in per if k in old and old[k] != per[k])
        elif not added and not removed:
            return              # a split saved before fingerprints existed; membership agrees
        if not (added or removed or changed): return
        show = lambda xs: f"{len(xs)} ({', '.join(xs[:8])}{', ...' if len(xs) > 8 else ''})"
        msg = (f"{cls.name}: the data no longer matches its saved split {cls.split_file()}: "
               f"{show(added)} proteins not in the split, {show(removed)} in the split but not "
               f"in the data, {show(changed)} with a changed sequence. A split is only valid for "
               f"the data it was dealt from -- restore that data (see known_missing/), or build a "
               f"new split on purpose (PLUG_BUILD_SPLITS=1, after moving the old one away).")
        if os.environ.get("PLUG_ALLOW_SPLIT_MISMATCH") == "1":
            print("WARNING: " + msg, file=sys.stderr); return
        raise ValueError(msg)

    # location of cached splots. the settings are part of the file name, so a split made
    # at 20% identity is never reused after you change that setting to 30%.
    @classmethod
    def split_file(cls):
        return c.splits / f"{cls.name}__{c.split_tag()}.csv"

    # which split each protein is in, as {protein: "train"/"val"/"test"}. worked out
    # once, cached, and read back from then on — grouping the sequences takes
    # minutes on the big benchmarks, so it should only ever happen once per setting.
    @classmethod
    def split_map(cls):
        if cls._split is None:
            path = cls.split_file()
            # read it if saved. building one is opt-in: a re-dealt split reshuffles which
            # proteins are train/val/test (two extra atlas chains once moved 2/3 of its test
            # set), so a missing file must stop the run, not silently produce a new split
            if path.exists():
                cls._split = _read(path)
            elif os.environ.get("PLUG_BUILD_SPLITS") == "1":
                cls._split = _build(cls, path)
            else:
                raise FileNotFoundError(
                    f"no saved homology split {path}: refusing to deal a new one, which would "
                    f"reshuffle train/val/test. Point SPLITS_DIR at the saved splits, or set "
                    f"PLUG_BUILD_SPLITS=1 to build one on purpose.")
        return cls._split

    # label each item with the split its protein is in. for a pair, both proteins have
    # to be in the same split — if one is in train and the other in test, the pair
    # belongs to neither, so it is left out.
    @classmethod
    def stamp(cls, items, split_map):
        out = []
        for it in items:
            names = {split_map[it[k]] for k in cls.keyfields()}
            if len(names) == 1:
                it["split"] = names.pop()
                out.append(it)
        return out

# work out fresh homology splits: collect every protein and its sequence, count how many items
# each protein has, group the similar ones together, then deal out those groups out 80/10/10
# so that a protein and everything like it always end up in the same split.
def _build(cls, path):
    seqs, weight = {}, Counter()
    items, _ = cls.filtered(cls.rows())
    for it in items:
        for k, s in zip(cls.keyfields(), cls.seqfields()):
            seqs.setdefault(it[k], it[s])
            weight[it[k]] += 1
    # nothing to group, which almost always means the data was never downloaded
    if not seqs:
        raise FileNotFoundError(f"{cls.name} has no items to split — "
                                "download it first, see data/download_*.sh")
    split_map = homology.deal(homology.clusters(seqs), weight)
    _write(path, split_map)
    write_sidecar(cls, path, items)
    return split_map


# the fingerprint record saved next to a split: <split>.json
def sidecar(path):
    return Path(path).with_suffix(".json")


def write_sidecar(cls, path, items):
    fp, per = cls.fingerprint(items)
    sidecar(path).write_text(json.dumps({
        "benchmark": cls.name, "split_file": Path(path).name, "tag": c.split_tag(),
        "seed": c.seed, "min_id": c.min_id, "cov": c.cov, "n_proteins": len(per),
        "excluded": sorted(cls.excluded()), "fingerprint": fp, "written": date.today().isoformat(),
        "proteins": per}, indent=1))

# read a saved split back in
def _read(path):
    with open(path, newline="") as f:
        return {r["protein"]: r["split"] for r in csv.DictReader(f)}

# save a split as a two column csv: the protein, and the split it is in
def _write(path, split_map):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("protein", "split"))
        w.writerows(sorted(split_map.items()))

from .megascale import Megascale
from .fireprot import Fireprot
from .dms import Dms
from .cafa5 import Cafa5
from .allobench import Allobench
from .ppi import Ppi
from .shs27k import Shs27k
from .passerrank import Passerrank
from .lp_pdbbind import LpPdbbind
from .bindingdb import Bindingdb
from .atlas import Atlas
from .cryptobench import Apoholo, Cryptobench
from .mdcath import Mdcath


# every benchmark, looked up by name — this is what fastas.py and build_trainset loop
# over. a benchmark only turns up here if it was imported above, so keep that list current.
REGISTRY = {b.name: b for b in Benchmark.__subclasses__()}
# what `from plug.benchmarks import *` hands you: the benchmark classes, nothing else
__all__ = ["Benchmark", "REGISTRY", *(b.__name__ for b in REGISTRY.values())]