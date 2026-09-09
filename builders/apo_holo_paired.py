#!/usr/bin/env python
"""Build node-for-node matched apo and holo files, so the pair isolates one variable.

`sets/apo_holo/` answers "how does a method do on apo input". It cannot answer "how much
does apo input *cost*", because its 15 targets are a different population from the holo
sets — hand-picked, heavily evidenced systems, three of whose organiser-mandated pairs
their own audit rejected. ALPS scores 0.694 there against 0.592 and 0.612 on the two
holo sets, and reading that as "apo is easier" would be reading a target-selection
difference as a conformational one.

The control that isolates the variable is the same protein scored twice. This builder
writes both members of each arm with:

* the **same node set** — the intersection of the frozen apo node set with the residues
  the holo entry models, by author numbering
* the **same anchor**, the arm's frozen active site restricted to that intersection
* the **same labels**, the frozen `label_residues`, which the source already derived by
  transferring the holo footprint onto the apo

so the only thing that differs between `<arm>__apo.npz` and `<arm>__holo.npz` is which
coordinate file the positions came from. A paired test over those is a measurement of
the conformational cost and of nothing else.

Numbering agreement is the thing that can silently ruin this, so it is reported per arm
and an arm is dropped rather than written when the two entries do not share enough
residues. The ABL1 arms are the known risk: the manifest records that those entries use
different numbering conventions.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from from_quantum_allostery import parse_cif_atoms                # noqa: E402

MIN_SHARE = 0.60      # of the frozen apo node set


def build(arm, frozen, spec, structures):
    apo_pdb, apo_ch = spec["apo"]["pdb"], spec["apo"]["chain"]
    holo_pdb, holo_ch = spec["holo"]["pdb"], spec["holo"]["chain"]
    t = frozen["targets"][arm]
    ids = list(t["residue_ids"])

    pa = os.path.join(structures, "apo", f"{apo_pdb}.cif.gz")
    ph = os.path.join(structures, "holo", f"{holo_pdb}.cif.gz")
    if not (os.path.exists(pa) and os.path.exists(ph)):
        return None, "missing structure"
    apo, _, apo_het = parse_cif_atoms(pa, apo_ch)
    holo, _, holo_het = parse_cif_atoms(ph, holo_ch)

    shared = [r for r in ids if r in apo and r in holo]
    frac = len(shared) / len(ids)
    if frac < MIN_SHARE:
        return None, (f"only {len(shared)}/{len(ids)} ({frac:.0%}) of the frozen node "
                      f"set is modelled in both {apo_pdb}:{apo_ch} and "
                      f"{holo_pdb}:{holo_ch}")

    pos = {r: i for i, r in enumerate(shared)}
    y = np.zeros(len(shared), np.int64)
    for r in t["label_residues"]:
        if r in pos:
            y[pos[r]] = 1
    anchor = np.array(sorted(pos[r] for r in t["active_site"] if r in pos), np.int64)
    if anchor.size == 0 or int(y.sum()) == 0:
        return None, (f"after intersection: {int(y.sum())} labels, "
                      f"{anchor.size} anchor residues")

    out = {}
    for side, table, pdb in (("apo", apo, apo_pdb), ("holo", holo, holo_pdb)):
        out[side] = dict(
            cb=np.stack([table[r][1] for r in shared]).astype(np.float32),
            resnums=np.asarray(shared, np.int32),
            anchor=anchor,
            y=y,
            chain_id=np.array([apo_ch if side == "apo" else holo_ch] * len(shared),
                              dtype="U4"),
            coord_type=np.asarray("CB"),
            label_rule=np.asarray("holo-derived-apo-paired"),
            polarity=np.asarray("positive"),
            conformation=np.asarray(side),
            source=np.asarray("quantum-allostery"),
            licence_tier=np.asarray("redistributable"),
            anchor_rule=np.asarray("ligand-derived"
                                   if "from_ligands" in spec["active_site"]
                                   else "motif-derived"),
            chain_id_source=np.asarray("recorded"),
            pdb=np.asarray(pdb),
            apo_pdb=np.asarray(apo_pdb),
            holo_pdb=np.asarray(holo_pdb),
            uniprot=np.asarray(str(spec.get("uniprot", ""))),
            modulator=np.asarray(str(spec["holo"].get("ligand", ""))),
            arm=np.asarray(arm),
            node_set=np.asarray("apo frozen node set intersected with holo-modelled "
                                "residues; identical on both members"),
            node_share=np.asarray(float(frac)),
        )
    return out, f"{len(shared)}/{len(ids)} shared ({frac:.0%}), {int(y.sum())} labels"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="checkout of quantum-allostery")
    ap.add_argument("--out", default=os.path.join(HERE, "sets", "apo_holo_paired"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, os.path.join(HERE, "spec"))
    from validate import check                                    # noqa: E402

    made, skipped = 0, []
    for tier in ("primary", "secondary"):
        base = os.path.join(a.src, "docs", "benchmark", tier)
        frozen = json.load(open(os.path.join(base, "frozen.json")))
        by_id = {t["id"]: t for t in
                 yaml.safe_load(open(os.path.join(base, "manifest.yaml")))["targets"]}
        for arm in frozen["targets"]:
            if arm not in by_id:
                skipped.append((arm, "no manifest entry"))
                continue
            pair, note = build(arm, frozen, by_id[arm],
                               os.path.join(a.src, "structures"))
            if pair is None:
                skipped.append((arm, note))
                continue
            bad = []
            for side, rec in pair.items():
                p = os.path.join(a.out, f"{arm}__{side}.npz")
                np.savez_compressed(p, **rec)
                bad += [f"{side}: {b}" for b in check(p)]
            if bad:
                for side in pair:
                    os.remove(os.path.join(a.out, f"{arm}__{side}.npz"))
                skipped.append((arm, bad[0]))
                continue
            made += 1
            print(f"  {arm:28s} {note}")
    print(f"\n{made} matched pairs -> {a.out}")
    for arm, why in skipped:
        print(f"  SKIPPED {arm}: {why}")


if __name__ == "__main__":
    main()
