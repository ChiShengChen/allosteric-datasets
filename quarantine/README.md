# Quarantine

Files that came out of the conversion without conforming to `spec/FORMAT.md`.

They are held rather than repaired: every repair here means dropping residues from the
node set, which shifts every `anchor` index, and that is a decision for the source
repository rather than for a converter. `builders/from_qasc.py` rewrites this list on
every run.

- `curated/3KGF_1` — resnums not strictly ascending within chain `B`
- `curated/3KGF_2` — resnums not strictly ascending within chain `B`
- `curated/3PJG` — resnums not strictly ascending within chain `A`
- `multimer/4HNB_multi` — resnums not strictly ascending within chain `B`
- `negatives/1A5Z_A` — resnums not strictly ascending within chain `A`

---

## What is actually wrong with each

The validator's message is the same for all five and the causes are not. Two distinct
defects, with different consequences.

### A. Non-amino-acid residues became graph nodes

| file | node indices | residue numbers | in anchor? | labelled? |
|---|---|---|---|---|
| `curated/3KGF_1` | 905–909 | 9003, 9004, 9001, 9004, 4188 | no | no |
| `curated/3KGF_2` | 905–909 | 9003, 9004, 9001, 9004, 4188 | no | no |
| `curated/3PJG` | 0 | 999 | no | no |

Residue numbers in the 9000s and 4188 are the conventional deposit numbering for
heteroatoms — ligands, ions, modified groups — not for amino acids. They carry no label
and are not anchors, so **the labels are unaffected**. But they are **nodes in the
contact graph**, and every spectral quantity in the source repository is a functional of
that graph: five extra nodes in 3KGF's 910 add edges to whatever they sit near, and every
eigenvalue, every ALPS perturbation response and every message-passing step sees them.

3KGF also has `9004` twice, so two distinct nodes share one identifier.

The magnitude is small — 5 nodes in 910, 1 in 389 — and the direction is unknown without
rerunning. These are curated-set members, so they are inside the 96 targets behind every
curated number.

### B. Insertion codes collapsed onto duplicate residue numbers

| file | duplicated `(chain, resnum)` |
|---|---|
| `negatives/1A5Z_A` | A/132, A/209, A/210, A/330 — six residues over four numbers |
| `multimer/4HNB_multi` | B/1 |

PDB insertion codes distinguish `209`, `209A`, `209B`; the builder keeps the number and
drops the code, so several residues end up sharing a key. `1A5Z_A` has residue 209 four
times.

This one **can** move labels. The curated annotations map to residues by *chain and
residue number*, so a duplicated key is ambiguous by construction: an annotation for
`A/209` matches four nodes. Neither of these two files is label-carrying — `1A5Z_A` is a
protein-level negative and `4HNB_multi` is a proxy-labelled multimer — so no label is
known to have landed wrongly here. It was not established whether any *other* file, one that
passed the ascending check because its duplicates happened to be adjacent and equal, had
the same defect without the symptom. **It was then checked, and six more do:**

| set | files with a duplicate `(chain, resnum)` |
|---|---|
| `targets_allobench` | **6** — `3BEU_A_NA_8dae1e` (9 duplicates), `3I78_A_NA_3f0b5d` (9), `3HO8_A_COA_8c3a68` (6), `3QH0_A_PLM_f35a37` (1), and two more |
| `targets_multimer` | 1 — `4HNB_multi` |
| `targets_negative` | 1 — `1A5Z_A` |
| every other set | 0 |

So the ordering check caught two of eight. The other six sit inside
`targets_allobench`, which is the 1,042-target set behind the AlloBench numbers, and
they are invisible to any check that looks at ordering rather than at uniqueness. This
is why the suggested check below is uniqueness and not ordering.

`1A5Z_A` is one of the 90 protein-level negatives, which are the evidence behind "nothing
exceeds AUC 0.57 against size-matched negatives".

## Suggested repair, for the source repository

Not applied here, because both repairs change the node set.

**For A:** drop residues whose numbering is outside the amino-acid range at build time,
and rebuild `anchor` afterwards. Better: filter on the record type in the source file
rather than on the number, so the rule does not depend on a deposit convention.

**For B:** keep the insertion code in the key — `resnums` as `int` plus an `icode`
array, or a string key `"209A"`. Any annotation mapping should then refuse to proceed
on an ambiguous key rather than take the first match.

**A check worth running regardless:** count duplicate `(chain, resnum)` pairs across
every target set. This directory holds the two files where the duplication also broke
the ordering, which is a symptom and not the defect.
