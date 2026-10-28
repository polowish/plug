#!/usr/bin/env python
# fetch mdcath domain files into a size-capped staging directory, for mdcath_extract.py
# --watch to consume (extract, then delete). runs on the LOGIN node: compute nodes have
# no internet. the login node kills every process when its ssh session ends (tmux too),
# so this is resumable by construction -- rerun it and it skips whatever is already
# extracted (<out>/done.txt), staged, or being extracted:
#
#   tmux new -s mdcath
#   ~/scratch/venvs/mdcath/bin/python scripts/mdcath_fetch.py --staging ~/scratch/mdcath_staging \
#       --out ~/fac/data/mdcath --cap-gb 300
#
# writes <staging>/FETCH_DONE when every domain is accounted for, which tells the
# extraction worker it may stop once the staging directory is empty.
import argparse
import os
import shutil
import sys
import time
from pathlib import Path

REPO = "compsciencelab/mdCATH"


def domain_list(cache: Path) -> list[str]:
    """Every domain in the dataset, from the hub listing (cached in <out>/domains.txt)."""
    if cache.exists():
        return [d for d in cache.read_text().split() if d]
    from huggingface_hub import list_repo_files

    files = list_repo_files(REPO, repo_type="dataset")
    doms = sorted(f[len("data/mdcath_dataset_"):-3] for f in files
                  if f.startswith("data/mdcath_dataset_") and f.endswith(".h5"))
    cache.parent.mkdir(parents=True, exist_ok=True)
    cache.write_text("\n".join(doms) + "\n")
    return doms


def staged_bytes(staging: Path) -> int:
    return sum(p.stat().st_size for p in staging.rglob("*.h5*") if p.is_file())


def fetch(dom: str, staging: Path, tries: int = 3) -> bool:
    from huggingface_hub import hf_hub_download

    tmp = staging / ".incoming"
    tmp.mkdir(exist_ok=True)
    for k in range(tries):
        try:
            p = hf_hub_download(REPO, f"mdcath_dataset_{dom}.h5", subfolder="data",
                                repo_type="dataset", local_dir=tmp)
            os.replace(p, staging / f"{dom}.h5")          # whole file or nothing
            return True
        except Exception as e:                              # network, 5xx, rate limit
            print(f"  {dom}: attempt {k + 1} failed: {e}", file=sys.stderr, flush=True)
            time.sleep(30 * (k + 1))
    return False


def main() -> int:
    p = argparse.ArgumentParser(description="fetch mdcath h5 files into a capped staging dir")
    p.add_argument("--staging", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True, help="the extractor's output dir (done.txt)")
    p.add_argument("--cap-gb", type=float, default=300.0)
    p.add_argument("--only", type=Path, help="a file of domain ids to fetch instead of all")
    args = p.parse_args()

    staging, out = args.staging, args.out
    staging.mkdir(parents=True, exist_ok=True)
    (staging / "FETCH_DONE").unlink(missing_ok=True)
    doms = (args.only.read_text().split() if args.only else domain_list(out / "domains.txt"))
    failed = []
    t0, got, nbytes = time.time(), 0, 0
    print(f"{len(doms)} domains; staging {staging}; cap {args.cap_gb:g} GB", flush=True)
    for i, dom in enumerate(doms):
        done = set((out / "done.txt").read_text().split()) if (out / "done.txt").exists() else set()
        if dom in done or any((staging / sub / f"{dom}.h5").exists() for sub in ("", "claimed", "failed")):
            continue
        while staged_bytes(staging) > args.cap_gb * 1e9:
            time.sleep(20)                                  # the worker is behind; let it drain
        if fetch(dom, staging):
            got += 1
            nbytes += (staging / f"{dom}.h5").stat().st_size if (staging / f"{dom}.h5").exists() else 0
        else:
            failed.append(dom)
        if got and got % 50 == 0:
            rate = nbytes / max(time.time() - t0, 1) / 1e6
            print(f"  [{time.strftime('%H:%M:%S')}] {i + 1}/{len(doms)} seen, {got} fetched, "
                  f"{rate:.0f} MB/s, staged {staged_bytes(staging) / 1e9:.0f} GB, {len(failed)} failed", flush=True)
    shutil.rmtree(staging / ".incoming", ignore_errors=True)
    if failed:
        (out / "fetch_failed.txt").write_text("\n".join(failed) + "\n")
        print(f"{len(failed)} domains failed: {out / 'fetch_failed.txt'} (rerun to retry)", flush=True)
        return 1
    (staging / "FETCH_DONE").write_text(time.strftime("%Y-%m-%dT%H:%M:%S") + "\n")
    print(f"all domains accounted for; {got} fetched this run", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
