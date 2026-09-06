# allosteric-datasets

Three projects built allosteric-site datasets without comparing notes. This repository
puts them in one format so they can be used together — and, more importantly, records
the four ways they are **not** interchangeable, because that is what a naive merge would
destroy.

```bash
python3 builders/from_qasc.py --qasc ../qasc_plus          # 472 targets, shipped
python3 builders/from_quantum_allostery.py --src ../quantum-allostery   # 15, shipped
python3 builders/from_allobench.py --src ../allosteric-dataset-pipeline # 1,423, local only
python3 spec/validate.py sets/*                            # 0 only if everything conforms
python3 overlap/measure.py                                 # regenerate overlap/README.md
```

## What is here

| set | targets | label rule | conformation | ships? |
|---|---|---|---|---|
| `curated` | 94 | expert-curated | holo | yes |
| `proxy_tier_a` | 11 | proxy-distal-druglike | holo | yes |
| `proxy_tier_b` | 90 | proxy-distal-druglike | holo | yes |
| `multimer` | 87 | proxy-distal-druglike | holo | yes |
| `negatives` | 89 | none-annotated | holo | yes |
| `matched` | 101 | expert-curated + none-annotated | holo | yes |
| **`apo_holo`** | **15** | **holo-derived-apo-paired** | **apo** | yes |
| `allobench` | 1,423 | 4A-heavy-atom-to-modulator | holo | **no — build it** |

487 targets in the repository, 1,910 once the AlloBench builder has run.

One `.npz` per target, one schema, every rule tag mandatory. The schema and what each
rule means: [`spec/FORMAT.md`](spec/FORMAT.md). Nothing loads without saying which rule
produced it, which is the whole design.

## Four things that must not be smoothed over

**1. Five label rules, and they disagree.** "Allosteric site" here means a curator said
so, or a drug-like molecule crystallised far from the catalytic site, or any residue
within 4 Å of an annotated modulator, or a ligand footprint transferred from a different
structure — or, for the negatives, nothing at all. These are not degrees of the same
measurement. On the two sets where both rules were computed, 26% of the positives under
the 4 Å rule fall inside an 8 Å distal mask that the curated rule leaves empty.
**The sets are reported side by side and never pooled.**

**2. The sets share structures, heavily.** `allobench` and `curated` share 82 PDB
entries — **89% of `curated`**. `matched` and `negatives` share 92% of `negatives` by
construction. The full matrix is [`overlap/README.md`](overlap/README.md), regenerated
by a script rather than written by hand.

The consequence is specific and easy to walk into: *train on the large set, test on the
small curated one* is the obvious next experiment and it is **89% contaminated**. It
will produce a clean-looking generalisation number that means nothing.

**3. `apo_holo` is a probe, not a benchmark — and it is not independent either.** It is
the only set whose coordinates are apo, so it is the only place the question "does this
survive on an unliganded structure" can be asked at all. But n = 15, and **8 of its 12
UniProt accessions also appear in `allobench`**. Two thirds of the apo probe is
protein-level contaminated against the largest training set here. Use it to ask whether
a method degrades on apo input; do not use it to claim generalisation.

**4. Sequence identity is not measured anywhere.** Two proteins at 90% identity are two
accessions and two PDB ids, and every overlap number above will call them disjoint. All
of those numbers are lower bounds.

## Licensing is two-tier, and it is not a preference

| source | coordinates | annotations | may this repository ship it? |
|---|---|---|---|
| curated (Wu et al. 2021 Table S2) | RCSB, public domain | open-access supplement | **yes** |
| quantum-allostery | RCSB, public domain | MIT, derived from holo entries | **yes** |
| AlloBench / ASD | RCSB, public domain | *"research use only; cite AlloBench and ASD, and do not redistribute them onward"* | **no — builder only** |

So `sets/allobench/` is gitignored, every file it produces carries
`licence_tier = build-only`, and anyone who wants those 1,423 targets clones upstream and
runs one command. Full terms and citations: [`PROVENANCE.md`](PROVENANCE.md).

## What the integration actually changed

**The AlloBench route is now Cβ.** Upstream ships `dataset/cb/` — Cβ coordinates, Cα
substituted at glycine, keyed identically to the Cα samples. That directory did not
exist when the pipeline was first vendored into `qasc_plus`, so every published number on
that route was computed on Cα while every other set here is Cβ. One of the three stated
incompatibilities between those two sets is now removable. The other two — label rule and
holo-vs-holo — are not.

**Five files are quarantined rather than shipped.** They came out of the conversion with
residue numbers that are not a usable key: insertion codes collapsed onto duplicate
numbers, or hetero residues numbered 9001 sorted into the middle of a chain. See
[`quarantine/README.md`](quarantine/README.md). They are held rather than repaired,
because repairing means dropping residues from a node set, which shifts every `anchor`
index, and that belongs to the source repository.

## Sources

- [leo07010/allosteric-dataset-pipeline](https://github.com/leo07010/allosteric-dataset-pipeline) — the AlloBench route
- [George930502/quantum-allostery](https://github.com/George930502/quantum-allostery) — the apo→holo arms
- [ChiShengChen/allosteric-benchmark](https://github.com/ChiShengChen/allosteric-benchmark) — the curated, proxy, negative and matched sets
