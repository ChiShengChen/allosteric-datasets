#!/usr/bin/env python
"""Mechanical check that a set conforms to FORMAT.md.

A format is only as good as the thing that refuses to load a file breaking it. Every
check here corresponds to an error that has actually occurred in one of the three source
repositories: an anchor index out of range after a node set was filtered, a label array
left as float, residue numbers that ascend everywhere except across a chain break, and a
rule tag that was absent because the writer knew what the set was.

Exit code is 0 only if every file passes.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np

REQUIRED_ARRAYS = ("cb", "resnums", "anchor", "y", "chain_id")
RULE_TAGS = {
    "coord_type": {"CB"},
    "label_rule": {"expert-curated", "proxy-distal-druglike",
                   "4A-heavy-atom-to-modulator", "holo-derived-apo-paired",
                   "none-annotated"},
    "polarity": {"positive", "negative"},
    "conformation": {"holo", "apo"},
    "source": {"curated", "allobench", "quantum-allostery"},
    "licence_tier": {"redistributable", "build-only"},
}


def check(path):
    """Return a list of complaints about one file. Empty means it conforms."""
    bad = []
    z = np.load(path, allow_pickle=True)
    have = set(z.files)

    for f in REQUIRED_ARRAYS:
        if f not in have:
            bad.append(f"missing array {f!r}")
    for tag, allowed in RULE_TAGS.items():
        if tag not in have:
            bad.append(f"missing rule tag {tag!r} — FORMAT.md gives it no default")
        elif str(z[tag]) not in allowed:
            bad.append(f"{tag}={str(z[tag])!r} is not one of {sorted(allowed)}")
    if bad:
        return bad                       # the rest would just cascade

    cb, resnums = z["cb"], z["resnums"]
    anchor, y, chain = z["anchor"], z["y"], z["chain_id"]
    n = len(cb)

    if cb.ndim != 2 or cb.shape[1] != 3:
        bad.append(f"cb has shape {cb.shape}, expected (N, 3)")
    if not np.isfinite(cb).all():
        bad.append("cb contains non-finite coordinates")
    for name, arr in (("resnums", resnums), ("y", y), ("chain_id", chain)):
        if len(arr) != n:
            bad.append(f"{name} has length {len(arr)}, expected {n}")

    if anchor.size == 0:
        bad.append("anchor is empty — the seeded formulation needs a seed")
    elif anchor.min() < 0 or anchor.max() >= n:
        bad.append(f"anchor index out of range: [{anchor.min()}, {anchor.max()}] "
                   f"against N = {n}")
    elif len(set(anchor.tolist())) != len(anchor):
        bad.append("anchor contains duplicate indices")

    if not np.isin(y, (0, 1)).all():
        bad.append(f"y is not binary: values {sorted(set(y.tolist()))[:6]}")
    elif int(y.sum()) == 0 and str(z["polarity"]) == "positive":
        bad.append("no positive residue on a target declared positive")
    elif int(y.sum()) > 0 and str(z["polarity"]) == "negative":
        bad.append(f"{int(y.sum())} positives on a target declared negative")

    # Ascending within a chain, not across one: a multimer restarts numbering.
    for c in np.unique(chain):
        r = resnums[chain == c]
        if len(r) > 1 and not (np.diff(r) > 0).all():
            bad.append(f"resnums not strictly ascending within chain {c!r}")

    if str(z["conformation"]) == "apo" and "holo_pdb" not in have:
        bad.append("apo file without holo_pdb — labels must name the structure they "
                   "came from")
    return bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", help="npz files or directories of them")
    ap.add_argument("--quiet", action="store_true")
    a = ap.parse_args()

    files = []
    for p in a.paths:
        files += (sorted(glob.glob(os.path.join(p, "*.npz"))) if os.path.isdir(p)
                  else [p])
    if not files:
        print("no .npz files found", file=sys.stderr)
        return 2

    failed = 0
    for f in files:
        try:
            bad = check(f)
        except Exception as e:                       # a file that will not open at all
            bad = [f"unreadable: {type(e).__name__}: {e}"]
        if bad:
            failed += 1
            print(f"FAIL {f}")
            for b in bad:
                print(f"       {b}")
        elif not a.quiet:
            print(f"ok   {os.path.basename(f)}")

    print(f"\n{len(files) - failed}/{len(files)} conform to spec/FORMAT.md")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
