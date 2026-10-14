# write the fingerprint record (<split>.json) for saved homology splits that predate it.
#
#   SPLITS_DIR=/path/to/splits python scripts/fingerprint_splits.py atlas allobench ...
#
# refuses a benchmark whose current data does not match its saved split exactly (same
# proteins), because the record would then vouch for data the split was not dealt from.
# items in known_missing/<name>.txt are left out first, as everywhere else.
import sys

from plug import config as c
from plug.benchmarks import REGISTRY, _read, sidecar, write_sidecar

bad = 0
for name in sys.argv[1:]:
    cls = REGISTRY[name]
    path = cls.split_file()
    if not path.exists():
        print(f"{name}: no saved split at {path}"); bad += 1; continue
    items, n_excl = cls.filtered(cls.rows())
    _, per = cls.fingerprint(items)
    split = _read(path)
    added, removed = sorted(set(per) - set(split)), sorted(set(split) - set(per))
    if added or removed:
        print(f"{name}: data does not match {path.name}: {len(added)} not in split {added[:5]}, "
              f"{len(removed)} in split but not in data {removed[:5]} -- no record written"); bad += 1
        continue
    write_sidecar(cls, path, items)
    print(f"{name}: {len(per)} proteins ({n_excl} items excluded) -> {sidecar(path).name}")
sys.exit(1 if bad else 0)
