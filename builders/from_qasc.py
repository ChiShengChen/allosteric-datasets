#!/usr/bin/env python
"""Convert the redistributable target sets of qasc_plus into the unified format.

Six sets come across. They are *not* one dataset in six pieces: they use two different
label rules and two different anchor rules, and one of them has no positives at all by
construction. The conversion's only real work is attaching those tags, because in the
source repository they live in build-script docstrings rather than in the files.

  qasc_plus/data/          ->  sets/                 label rule              polarity
  targets_curated (97)         curated/              expert-curated          positive
  targets (11)                 proxy_tier_a/         proxy-distal-druglike   positive
  targets_b (90)               proxy_tier_b/         proxy-distal-druglike   positive
  targets_multimer (88)        multimer/             proxy-distal-druglike   positive
  targets_negative (90)        negatives/            none-annotated          negative
  matched_pos + matched_neg    matched/              (per member)            per member

`targets_curated_small` is not converted: it is the N <= 520 subset of `curated` and is
recoverable from it with one filter, so shipping it separately would give the same
target two identities.

**Anchor rule travels with the file too.** Section 10.1 of the source repository found
that its protein-level experiment was not measuring what it claimed, because positives
took the anchor from the curated table and negatives took it from a 4.5 A cofactor
contact, and every method is seeded at the anchor. The `matched/` set exists because of
that finding. A reader who cannot see which anchor rule produced a file cannot see that
class of bug, so `anchor_rule` is written into every file here.

Chain is recovered from the filename where the source did not store it: those builders
write `<pdb>_<chain>.npz` and keep a single chain per file.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "spec"))
from validate import check                                        # noqa: E402

# (source dir, output dir, label_rule, polarity, anchor_rule)
SETS = [
    ("targets_curated", "curated", "expert-curated", "positive", "curated-table"),
    ("targets", "proxy_tier_a", "proxy-distal-druglike", "positive",
     "cofactor-contact"),
    ("targets_b", "proxy_tier_b", "proxy-distal-druglike", "positive",
     "cofactor-contact"),
    ("targets_multimer", "multimer", "proxy-distal-druglike", "positive",
     "cofactor-contact"),
    ("targets_negative", "negatives", "none-annotated", "negative",
     "cofactor-contact"),
    ("matched_pos", "matched", "expert-curated", "positive", "uniform-matched"),
    ("matched_neg", "matched", "none-annotated", "negative", "uniform-matched"),
]

KEY = re.compile(r"^(?P<pdb>[0-9A-Za-z]{4})(?:_(?P<chain>[A-Za-z0-9]+))?")


def convert_one(src, label_rule, polarity, anchor_rule):
    z = np.load(src, allow_pickle=True)
    stem = os.path.basename(src)[:-4]
    m = KEY.match(stem)
    pdb = m.group("pdb").upper() if m else ""
    n = len(z["cb"])

    if "chain_id" in z.files:
        chain = np.asarray(z["chain_id"], dtype="U4")
        chain_source = "recorded"
    else:
        c = (m.group("chain") if m and m.group("chain") else "A")
        chain = np.array([c] * n, dtype="U4")
        chain_source = "from-filename" if m and m.group("chain") else "assumed-A"

    rec = dict(
        cb=np.asarray(z["cb"], np.float32),
        resnums=np.asarray(z["resnums"], np.int32),
        anchor=np.asarray(z["anchor"], np.int64),
        y=np.asarray(z["y"], np.int64),
        chain_id=chain,
        coord_type=np.asarray("CB"),
        label_rule=np.asarray(label_rule),
        polarity=np.asarray(polarity),
        conformation=np.asarray("holo"),
        source=np.asarray("curated"),
        licence_tier=np.asarray("redistributable"),
        anchor_rule=np.asarray(anchor_rule),
        chain_id_source=np.asarray(chain_source),
        pdb=np.asarray(pdb),
    )
    return stem, rec


def main():
    ap = argparse.ArgumentParser()
    HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--qasc", default=os.path.join(os.path.dirname(HERE), "qasc_plus"),
                    help="checkout of the qasc_plus repository")
    ap.add_argument("--out", default=os.path.join(HERE, "sets"))
    a = ap.parse_args()

    total, allheld = 0, []
    for src_dir, out_dir, rule, polarity, anchor_rule in SETS:
        src = os.path.join(a.qasc, "data", src_dir)
        files = sorted(glob.glob(os.path.join(src, "*.npz")))
        if not files:
            print(f"  {src_dir:18s} MISSING at {src}")
            continue
        dst = os.path.join(a.out, out_dir)
        os.makedirs(dst, exist_ok=True)
        quarantine = os.path.join(a.out, os.pardir, "quarantine", out_dir)
        kept, held = 0, []
        for f in files:
            stem, rec = convert_one(f, rule, polarity, anchor_rule)
            # matched/ merges two source dirs; keep the polarity in the filename so a
            # positive and a negative of the same PDB cannot overwrite each other.
            name = (f"{stem}_{polarity[:3]}" if out_dir == "matched" else stem)
            tmp = os.path.join(dst, f"{name}.npz")
            np.savez_compressed(tmp, **rec)
            bad = check(tmp)
            if bad:
                # Non-conforming files are held, never silently repaired and never
                # quietly shipped. Repairing them means changing the node set, which
                # shifts every anchor index, and that is the source repository's call.
                os.makedirs(quarantine, exist_ok=True)
                os.replace(tmp, os.path.join(quarantine, f"{name}.npz"))
                held.append((name, bad))
            else:
                kept += 1
        note = f"  ({len(held)} quarantined)" if held else ""
        print(f"  {src_dir:18s} -> sets/{out_dir:14s} {kept:4d} files "
              f"[{rule}, {polarity}]{note}")
        for name, bad in held:
            print(f"        HELD {name}: {bad[0]}")
        total += kept
        allheld.extend((out_dir, n, b) for n, b in held)
    print(f"\n{total} targets converted, {len(allheld)} quarantined")
    if allheld:
        with open(os.path.join(a.out, os.pardir, "quarantine", "README.md"), "w") as fh:
            fh.write("# Quarantine\n\nFiles that came out of the conversion without "
                     "conforming to `spec/FORMAT.md`.\n\nThey are held rather than "
                     "repaired: every repair here means dropping residues from the node "
                     "set, which shifts every `anchor` index, and that is a decision for "
                     "the source repository rather than for a converter.\n\n")
            for out_dir, name, bad in allheld:
                fh.write(f"- `{out_dir}/{name}` — {'; '.join(bad)}\n")


if __name__ == "__main__":
    main()
