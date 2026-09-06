#!/usr/bin/env python
"""Convert the AlloBench route into the unified format — locally, never committed.

**This set is not redistributable and this builder is the only thing in this repository
that touches it.** The coordinates are RCSB and public domain, but the residue
annotations reach us through AlloBench from the Allosteric Site Database, whose terms are
"research use only; cite AlloBench and ASD, and do not redistribute them onward". So the
machinery ships and the labels do not, exactly as upstream does it. Output lands in
`sets/allobench/`, which `.gitignore` excludes, and every file it writes carries
`licence_tier = build-only` so a downstream script can refuse to publish one.

Run it against a checkout of leo07010/allosteric-dataset-pipeline:

    git clone https://github.com/leo07010/allosteric-dataset-pipeline
    python3 builders/from_allobench.py --src allosteric-dataset-pipeline

## The reason this builder exists rather than a pointer to the upstream files

**Upstream now ships Cβ.** `dataset/cb/` carries `cb_coords` (Cα substituted at glycine)
keyed identically to `dataset/samples/`, which carries Cα in `coords`. Every other set
in this repository is Cβ, and the qasc_plus analysis of this same route was run on Cα
because that directory did not exist when it was vendored. That mismatch is not
cosmetic there: its ALPS `RADIUS = 12.0` was tuned on Cβ contact geometry, and its
README lists the coordinate difference as one of three load-bearing incompatibilities
between the two sets.

Converting from `cb/` removes one of the three. It does not remove the other two — the
label rule is still 4 Å-to-modulator against expert curation, and both sets are still
holo — so these files are still never pooled with `sets/curated/`. What it does mean is
that a number produced from this builder and a number produced on the curated set now
differ by two stated things instead of three.

`resnums` is asserted identical between the two directories rather than assumed: they
are separate files and a silent misalignment would put every label on the wrong residue.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def load_folds(src):
    """UniProt-grouped fold assignment, keyed by sample, from the upstream manifest."""
    p = os.path.join(src, "metadata", "manifest.json")
    if not os.path.exists(p):
        return {}
    m = json.load(open(p))
    # The manifest keeps folds as {fold_id: [sample keys]}, not as a field on each
    # sample. Inverting it here rather than guessing: an earlier version of this
    # loader looked for a "fold" key on the manifest entries, found none, and
    # attached folds to nothing while reporting success.
    out = {}
    for fold_id, keys in m.get("folds", {}).items():
        for k in keys:
            out[k] = int(fold_id)
    return out


def convert_one(sample_path, cb_path, folds):
    key = os.path.basename(sample_path)[:-4]
    s = np.load(sample_path, allow_pickle=True)
    c = np.load(cb_path, allow_pickle=True)

    if not np.array_equal(s["resnums"], c["resnums"]):
        return key, None, "resnums differ between samples/ and cb/"

    meta = json.loads(str(s["meta"])) if "meta" in s.files else {}
    y = np.asarray(s["allo_labels"], np.int64)
    anchor = np.where(np.asarray(s["active_site_mask"]).astype(bool))[0].astype(np.int64)
    if anchor.size == 0:
        return key, None, "no active-site residue; the seeded formulation needs a seed"
    if int(y.sum()) == 0:
        return key, None, "no positive residue"

    chain = str(meta.get("chain", "A"))
    n = len(y)
    rec = dict(
        cb=np.asarray(c["cb_coords"], np.float32),
        resnums=np.asarray(s["resnums"], np.int32),
        anchor=anchor,
        y=y,
        chain_id=np.array([chain] * n, dtype="U4"),
        resnames=np.asarray(c["resnames"], dtype="U3"),
        coord_type=np.asarray("CB"),
        label_rule=np.asarray("4A-heavy-atom-to-modulator"),
        polarity=np.asarray("positive"),
        conformation=np.asarray("holo"),
        source=np.asarray("allobench"),
        licence_tier=np.asarray("build-only"),
        anchor_rule=np.asarray("asd-annotated-active-site"),
        chain_id_source=np.asarray("recorded"),
        pdb=np.asarray(str(meta.get("name", meta.get("pdb", "")))),
        uniprot=np.asarray(str(meta.get("uniprot", ""))),
        modulator=np.asarray(str(meta.get("modulator", ""))),
        terms=np.asarray("research use only; cite AlloBench and ASD; "
                         "do not redistribute onward"),
    )
    if key in folds:
        rec["fold"] = np.asarray(int(folds[key]))
    return key, rec, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True,
                    help="checkout of leo07010/allosteric-dataset-pipeline")
    ap.add_argument("--out", default=os.path.join(HERE, "sets", "allobench"))
    a = ap.parse_args()

    sys.path.insert(0, os.path.join(HERE, "spec"))
    from validate import check                                      # noqa: E402

    samples = sorted(glob.glob(os.path.join(a.src, "dataset", "samples", "*.npz")))
    if not samples:
        print(f"no samples under {a.src}/dataset/samples — see this file's docstring",
              file=sys.stderr)
        return 2
    folds = load_folds(os.path.join(a.src))
    os.makedirs(a.out, exist_ok=True)

    written, reasons, no_cb = 0, {}, 0
    for sp in samples:
        cp = os.path.join(a.src, "dataset", "cb", os.path.basename(sp))
        if not os.path.exists(cp):
            no_cb += 1
            continue
        key, rec, why = convert_one(sp, cp, folds)
        if rec is None:
            reasons[why] = reasons.get(why, 0) + 1
            continue
        p = os.path.join(a.out, f"{key}.npz")
        np.savez_compressed(p, **rec)
        bad = check(p)
        if bad:
            os.remove(p)
            reasons[bad[0]] = reasons.get(bad[0], 0) + 1
            continue
        written += 1

    print(f"{written} of {len(samples)} samples converted -> {a.out}")
    print(f"folds attached: {sum(1 for f in glob.glob(os.path.join(a.out, '*.npz')) if 'fold' in np.load(f).files)}")
    if no_cb:
        print(f"  {no_cb} had no matching file in dataset/cb/")
    for why, n in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"  {n:4d} dropped: {why}")
    print("\nNOT REDISTRIBUTABLE. These files carry licence_tier = build-only and the "
          "output directory is gitignored.")


if __name__ == "__main__":
    main()
