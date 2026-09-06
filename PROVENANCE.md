# Provenance and terms

Per source: where it came from, what its terms allow, what this repository therefore
ships, and how to cite it. The short version is in the README; this is the part a lawyer
or a reviewer reads.

Coordinates are the easy half. Every structure in every set originates from the RCSB
PDB and is public domain. **Annotations are the hard half**, and they are what the two
tiers below are about.

---

## 1. AlloBench route — `sets/allobench/`, **not shipped**

**Upstream** [leo07010/allosteric-dataset-pipeline](https://github.com/leo07010/allosteric-dataset-pipeline),
MIT. Its `dataset/` directory holds 1,439 samples as `samples/*.npz` (Cα, residue
numbers, labels, active-site mask, metadata) with a parallel `cb/*.npz` (Cβ, Cα at
glycine, residue names), identically keyed.

**Terms on the annotations**, quoted from `dataset/README.md`:

> Coordinates originate from RCSB and are public domain. […] research use only; cite
> AlloBench and ASD, and do not redistribute them onward.

The labels derive from the Allosteric Site Database through AlloBench. Upstream
deliberately does not redistribute ASD-derived annotations, and neither does this
repository. `builders/from_allobench.py` is the only file here that touches them,
`sets/allobench/` is gitignored, and every file it writes carries
`licence_tier = build-only` and a `terms` string repeating the clause above.

**Cite** AlloBench and ASD when using anything this builder produces.

**What the builder changes** relative to how this route has been used before: it reads
`cb/` rather than `samples/`, so the output is Cβ. Every other set here is Cβ. The
earlier vendoring of this pipeline predates that directory and was necessarily Cα.

## 2. quantum-allostery — `sets/apo_holo/`, shipped

**Upstream** [George930502/quantum-allostery](https://github.com/George930502/quantum-allostery),
**MIT**. Fifteen frozen arms: six primary (three disease areas × mandated/corrected) and
nine secondary, each an apo structure paired with a holo entry whose ligand footprint at
4.5 Å becomes the label set.

Its own reuse condition is worth honouring literally: the freezes are versioned
artifacts and *"a number produced under one protocol version is not comparable to one
produced under another"*. Every file converted here records
`protocol = "quantum-allostery frozen input layer, frozen_on <date>"` so the version it
came from travels with it.

Taken from their freeze rather than re-derived: the node set (`residue_ids`), the active
site (derived from a rule, `from_ligands` or `from_motifs`, never a pinned list), the
label set, and `excluded_from_scoring` — residues that score by construction and leave
both classes under their ADR 0011.

**Three limitations of that source apply to every number drawn from these 15 arms**, and
they are theirs, stated plainly in their own `docs/benchmark/README.md`:

1. all fifteen arms use a synthetic small molecule as the effector, so no arm tests
   classical allosteric enzymology — cooperativity, metabolite feedback, a physiological
   effector;
2. the negative class has an unknown false-negative rate, so precision-style endpoints
   are more trustworthy here than recall-style ones;
3. ground truth is a **binding-site** label set, not a **coupling** label set. No
   structure pair can show that a method recovered coupling rather than a pocket.

**Cite** the repository and state the protocol version.

## 3. Curated, proxy, negative and matched sets — shipped

**Upstream** [ChiShengChen/allosteric-benchmark](https://github.com/ChiShengChen/allosteric-benchmark)
(`qasc_plus`), MIT.

`curated` and the positive half of `matched` derive from **Supplementary Table S2** of

> Wu, Amor, Schaub, Barahona et al., "Prediction of allosteric sites and signaling:
> Insights from benchmarking datasets", *Patterns* 2 (2021) 100408,
> doi:[10.1016/j.patter.2021.100408](https://doi.org/10.1016/j.patter.2021.100408)

which lists allosteric-site and active-site residues for 118 proteins drawn from ASBench
with sites from ASD Release 4.10. **Open access, and the supplement is redistributable**
— which is why these labels may ship while the AlloBench ones may not, even though both
trace back to ASD. The route matters: one is a published open-access table, the other is
a database extract under its own redistribution clause.

`proxy_tier_a`, `proxy_tier_b` and `multimer` use the field's operational convention
rather than curation: a residue contacting a drug-like ligand whose contact shell is at
least 8 Å from every anchor residue. `negatives` are structures carrying a cofactor and
no annotated allosteric site anywhere.

**Cite** Wu et al. 2021 for the curated labels, and the source repository for the sets.

---

## Anchor rules, and why they are recorded per file

`anchor_rule` is written into every file because the sets do not agree on it either:

| value | meaning | sets |
|---|---|---|
| `curated-table` | active-site residues as listed by Wu et al. Table S2 | `curated` |
| `cofactor-contact` | residues within a cutoff of a cofactor ligand | `proxy_*`, `multimer`, `negatives` |
| `uniform-matched` | one identical procedure applied to both classes | `matched` |
| `ligand-derived` / `motif-derived` | the frozen rule of that arm | `apo_holo` |
| `asd-annotated-active-site` | the mask AlloBench ships | `allobench` |

This is not bookkeeping. The source repository's §10.1 found that its protein-level
experiment was not measuring what it claimed, because positives took the anchor from the
curated table while negatives took it from a cofactor contact — and every method there is
seeded at the anchor, so every score inherited the difference. The `matched` set exists
to repair exactly that. A reader who cannot see which anchor rule produced a file cannot
see that class of error, so the tag travels with the data.

## This repository's own licence

MIT, covering the format specification, the validator, the builders and the
documentation — that is, everything except the data. Data carries the terms of its
source as set out above. `LICENSE` has the text.
