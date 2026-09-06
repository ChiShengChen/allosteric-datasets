#!/usr/bin/env python
"""Convert the frozen apo->holo arms of George930502/quantum-allostery.

This is the only source here whose coordinates are **apo**. Both other sources label and
score the same holo structure, so a method there can be rewarded for finding a pocket the
ligand itself is holding open. These fifteen arms score an unliganded structure against a
label set transferred from a separate holo entry, which is the condition real use faces.

Fifteen arms is not a benchmark. It is a probe, and `sets/apo_holo/README.md` says so.

What the source froze, and what is therefore taken rather than re-derived:

* `residue_ids` -- the node set a method receives, chain A only
* `active_site` -- derived from a *rule* (`from_ligands` or `from_motifs`), never a
  pinned list, so their `allo benchmark verify` re-derives and checks it
* `label_residues` -- the holo ligand footprint at 4.5 A, transferred onto the apo
* `excluded_from_scoring` -- residues that score by construction and leave both classes
  (their ADR 0011). Carried across as positions, because a pool defined differently is a
  different experiment

Coordinates come from their tracked `structures/apo/*.cif.gz`, whose sha256 the freeze
pins against versioned wwPDB artifacts. mmCIF is parsed here directly: the environment
has no gemmi and no Biopython, and the parse needed is one loop over `_atom_site`.

**Numbering is verified, not assumed.** Their `residue_ids` are author numbers on some
arms and start at zero on others, so both `auth_seq_id` and `label_seq_id` are tried and
the one that actually covers the frozen node set is used. A mismatch is reported and the
arm is skipped rather than written with residues silently dropped.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def parse_cif_atoms(path, chain):
    """Return {auth_seq: (resname, cb_xyz)} and {label_seq: ...} for one chain.

    Model 1, CB where present and CA for glycine, first altloc only.

    HETATM is read as well as ATOM. Modified residues -- MSE, CSD, SEP, TPO, PTR --
    are deposited as HETATM and are genuine protein residues; 1RTJ:A/280 is CSD and
    dropping it cost the whole hiv_rt arm until this was fixed. Waters are excluded
    explicitly and everything else is filtered by the caller, which keeps only residues
    the source's freeze already placed in the node set.
    """
    cols, rows, in_loop = {}, [], False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            t = line.strip()
            if t.startswith("_atom_site."):
                cols[t.split(".", 1)[1]] = len(cols)
                in_loop = True
                continue
            if in_loop:
                if not t or t.startswith("#") or t.startswith("loop_"):
                    if rows:
                        break
                    continue
                rows.append(t.split())
    need = ("group_PDB", "label_atom_id", "label_comp_id", "auth_asym_id",
            "auth_seq_id", "label_seq_id", "Cartn_x", "Cartn_y", "Cartn_z",
            "pdbx_PDB_model_num", "label_alt_id")
    missing = [c for c in need if c not in cols]
    if missing:
        raise ValueError(f"{os.path.basename(path)}: _atom_site lacks {missing}")
    g = {c: cols[c] for c in need}

    by_auth, by_label, hetero = {}, {}, set()
    for r in rows:
        if len(r) <= max(g.values()) or r[g["group_PDB"]] not in ("ATOM", "HETATM"):
            continue
        if r[g["label_comp_id"]] == "HOH":
            continue
        if r[g["pdbx_PDB_model_num"]] != "1" or r[g["auth_asym_id"]] != chain:
            continue
        if r[g["label_alt_id"]] not in (".", "?", "A"):
            continue
        atom, comp = r[g["label_atom_id"]], r[g["label_comp_id"]]
        want = "CA" if comp == "GLY" else "CB"
        if atom != want:
            continue
        xyz = np.array([float(r[g[c]]) for c in ("Cartn_x", "Cartn_y", "Cartn_z")],
                       np.float32)
        try:
            k = int(r[g["auth_seq_id"]])
            by_auth.setdefault(k, (comp, xyz))
            if r[g["group_PDB"]] == "HETATM":
                hetero.add(k)
        except ValueError:
            pass
        try:
            by_label.setdefault(int(r[g["label_seq_id"]]), (comp, xyz))
        except ValueError:
            pass
    return by_auth, by_label, hetero


def convert(arm, frozen, manifest_by_id, structures, tier):
    t = frozen["targets"][arm]
    spec = manifest_by_id[arm]
    apo_pdb, chain = spec["apo"]["pdb"], spec["apo"]["chain"]
    path = os.path.join(structures, "apo", f"{apo_pdb}.cif.gz")
    if not os.path.exists(path):
        return None, f"no structure at {path}"

    ids = list(t["residue_ids"])
    by_auth, by_label, hetero = parse_cif_atoms(path, chain)
    cover = {k: len(set(ids) & set(v)) for k, v in (("auth", by_auth),
                                                    ("label", by_label))}
    which = max(cover, key=cover.get)
    table = by_auth if which == "auth" else by_label
    if cover[which] < len(ids):
        return None, (f"numbering mismatch: {which}_seq_id covers "
                      f"{cover[which]}/{len(ids)} frozen residues "
                      f"(auth {cover['auth']}, label {cover['label']})")

    pos = {r: i for i, r in enumerate(ids)}
    cb = np.stack([table[r][1] for r in ids]).astype(np.float32)
    y = np.zeros(len(ids), np.int64)
    for r in t["label_residues"]:
        y[pos[r]] = 1
    anchor = np.array(sorted(pos[r] for r in t["active_site"] if r in pos), np.int64)
    excluded = np.array(sorted(pos[r] for r in t["excluded_from_scoring"] if r in pos),
                        np.int64)

    rec = dict(
        cb=cb,
        resnums=np.asarray(ids, np.int32),
        anchor=anchor,
        y=y,
        chain_id=np.array([chain] * len(ids), dtype="U4"),
        coord_type=np.asarray("CB"),
        label_rule=np.asarray("holo-derived-apo-paired"),
        polarity=np.asarray("positive"),
        conformation=np.asarray("apo"),
        source=np.asarray("quantum-allostery"),
        licence_tier=np.asarray("redistributable"),
        anchor_rule=np.asarray("ligand-derived" if "from_ligands" in spec["active_site"]
                               else "motif-derived"),
        chain_id_source=np.asarray("recorded"),
        pdb=np.asarray(apo_pdb),
        apo_pdb=np.asarray(apo_pdb),
        holo_pdb=np.asarray(spec["holo"]["pdb"]),
        modulator=np.asarray(spec["holo"].get("ligand", "")),
        uniprot=np.asarray(spec.get("uniprot", "")),
        sha256=np.asarray(spec["apo"].get("sha256", "")),
        numbering=np.asarray(f"{which}_seq_id"),
        modified_residues=np.asarray(sorted(hetero & set(ids)), np.int32),
        tier=np.asarray(tier),
        arm_tier=np.asarray(str(t.get("tier", ""))),
        n_candidates=np.asarray(int(t["n_candidates"])),
        excluded_from_scoring=excluded,
        contact_cutoff=np.asarray(float(frozen["contact_cutoff_angstrom"])),
        protocol=np.asarray("quantum-allostery frozen input layer, "
                            f"frozen_on {frozen['frozen_on']}"),
    )
    return rec, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True, help="checkout of quantum-allostery")
    ap.add_argument("--out", default=os.path.join(HERE, "sets", "apo_holo"))
    a = ap.parse_args()

    os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, os.path.join(HERE, "spec"))
    from validate import check                                      # noqa: E402

    written, skipped = 0, []
    for tier in ("primary", "secondary"):
        base = os.path.join(a.src, "docs", "benchmark", tier)
        frozen = json.load(open(os.path.join(base, "frozen.json")))
        man = yaml.safe_load(open(os.path.join(base, "manifest.yaml")))
        by_id = {t["id"]: t for t in man["targets"]}
        for arm in frozen["targets"]:
            if arm not in by_id:
                skipped.append((arm, "no manifest entry"))
                continue
            rec, why = convert(arm, frozen, by_id, os.path.join(a.src, "structures"),
                               tier)
            if rec is None:
                skipped.append((arm, why))
                continue
            p = os.path.join(a.out, f"{arm}.npz")
            np.savez_compressed(p, **rec)
            bad = check(p)
            if bad:
                os.remove(p)
                skipped.append((arm, bad[0]))
                continue
            written += 1
            print(f"  {arm:28s} N={len(rec['y']):4d} anchor={len(rec['anchor']):3d} "
                  f"pos={int(rec['y'].sum()):3d} {str(rec['apo_pdb'])}->"
                  f"{str(rec['holo_pdb'])} [{str(rec['numbering'])}]")
    print(f"\n{written} arms written to {a.out}")
    for arm, why in skipped:
        print(f"  SKIPPED {arm}: {why}")


if __name__ == "__main__":
    main()
