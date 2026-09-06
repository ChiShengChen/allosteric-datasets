# The format

One `.npz` per target. Every set in this repository converts into it, and the conversion
is lossy in one direction only: a field may be absent because a source does not carry it,
never present with a meaning that varies by source.

That last clause is the entire reason this repository exists. Three sources define an
"allosteric site" three different ways, on two different coordinate conventions, in two
different conformational states. Pooling them silently would be the most expensive
mistake available here, so the fields that say *which rule produced this file* are
mandatory and are checked mechanically by [`validate.py`](validate.py).

## Required arrays

| field | shape | dtype | meaning |
|---|---|---|---|
| `cb` | (N, 3) | float32 | Cβ coordinates, Cα for glycine. The node positions |
| `resnums` | (N,) | int32 | author residue numbers, ascending, the source's own numbering |
| `anchor` | (K,) | int64 | **indices into `cb`** of the active/orthosteric site |
| `y` | (N,) | int64 | 1 for an annotated allosteric-site residue, else 0 |
| `chain_id` | (N,) | U4 | chain per residue |

`anchor` is positional and `resnums` is nominal. Keeping both is what lets a file be
checked against its source without re-deriving the mapping, which is where the errors
live: two of the three sources use author numbering that does not start at 1, and one
has gaps.

## Required scalars — the rule tags

None of these has a default. A file without them does not load.

| field | values | why it must be explicit |
|---|---|---|
| `coord_type` | `CB` | every set is converted to Cβ; the tag stays so a future `CA` set cannot enter unlabelled |
| `label_rule` | `expert-curated`, `proxy-distal-druglike`, `4A-heavy-atom-to-modulator`, `holo-derived-apo-paired`, `none-annotated` | four incompatible definitions of a positive, plus the negatives which have none |
| `polarity` | `positive`, `negative` | a negative target is one with no annotated site at all; `y` is all zero **by construction**, not by failure |
| `conformation` | `holo`, `apo` | whether the effector was present in the structure the coordinates come from |
| `source` | `curated`, `allobench`, `quantum-allostery` | which repository and which terms apply |
| `licence_tier` | `redistributable`, `build-only` | whether this file may be committed anywhere |

## Optional scalars

`pdb`, `uniprot`, `modulator`, `fold`, `apo_pdb`, `holo_pdb`, `sha256`, `n_candidates`,
`excluded_from_scoring`. Present when the source carries them; never invented.

## What `label_rule` actually means

**`expert-curated`** — a human curator named these residues as the allosteric site.
Supplementary Table S2 of Wu, Amor, Schaub, Barahona et al., *Patterns* (2021),
doi:10.1016/j.patter.2021.100408, for 118 proteins drawn from ASBench with sites from
ASD Release 4.10. The strongest label rule here and the smallest set.

**`proxy-distal-druglike`** — a residue contacting any drug-like ligand (at least a
minimum heavy-atom count, not a cofactor, not a crystallisation additive) whose contact
shell is at least 8 Å from every anchor residue, with the anchor itself taken from
cofactor contacts. The field's operational convention, the same one ASBench and CASBench
use, and the rule behind most of the sets here. It is a *proxy*: "a drug-like molecule
crystallised there, far from the catalytic site" is not the same claim as "a curator
determined this is an allosteric site", and every conclusion drawn on these sets carries
that caveat.

**`4A-heavy-atom-to-modulator`** — any protein residue with a heavy atom within 4 Å of a
modulator heavy atom. Operational, reproducible, and systematically *closer to the active
site* than curated annotation: on the sets measured here, 26% of raw positives under this
rule fall inside an 8 Å distal mask that the curated rule leaves empty. Not a defect —
a different definition — but it is why the two are never pooled.

**`none-annotated`** — a protein-level negative. The structure carries a cofactor, so it
has an anchor and can be scored, but no allosteric site is annotated anywhere on it. The
question these exist to ask is the one a user actually faces: given a protein with no
site, does the method say so, or does it confidently name five residues anyway? `y` is
all zero and that is the correct value, so `validate.py` requires a positive only when
`polarity` is `positive`.

**`holo-derived-apo-paired`** — the label set is the ligand footprint in a **holo**
structure, transferred onto an **apo** structure of the same protein, which is what the
method actually receives. The only rule here where the coordinates a method sees and the
coordinates the labels come from are different files. Cutoff 4.5 Å, with 4.0 and 5.0
frozen alongside.

## Why `conformation` is a separate axis

The first two rules both label holo structures: the effector is present in the
coordinates being scored. A method can therefore be rewarded for finding a pocket that
the ligand itself is holding open. Only the third rule scores apo input, which is the
condition any real use faces, and it is the only place in this repository where "does
this survive on an unliganded structure" can be asked at all.

n = 15 there against n = 1,138 in the holo sets. It is a probe, not a benchmark, and
[`../sets/apo_holo/README.md`](../sets/apo_holo/README.md) says so at more length.

## Loading

```python
import numpy as np

z = np.load(path, allow_pickle=True)
assert str(z["coord_type"]) == "CB"
if str(z["label_rule"]) != "expert-curated":
    ...                      # do not compare across rules without saying so
cb, y, anchor = z["cb"], z["y"], z["anchor"]
```

`validate.py` refuses any file missing a required field, any `anchor` index out of range,
any `y` that is not 0/1, any `resnums` that is not ascending, and any unknown value of a
rule tag. Run it before committing a set and before trusting one.
