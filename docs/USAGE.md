# Attaching these energies to an evaluation set

This package answers one question: **for this `(pdb, mutation)`, what did FoldX compute?** It
ships no splits, no benchmark definitions, and no train/test partitions, and it never will —
those belong to whatever is being evaluated, not to the energies.

So every use follows the same shape: bring a list of rows, get their energies, and decide what to
do about the rows that have none.

## The join

```python
from skempi_foldx import MULTI_POINT, TERMS, load_bundled_store


def energies_for(rows, arms=("sp", "mp")):
    """rows: iterable of (pdb, mutation) using SKEMPI identifiers, multi-point comma-joined.

    Returns {(pdb, mutation): [12 floats in TERMS order]} for the rows that are covered.
    Rows with no FoldX result are absent rather than filled -- see "Uncovered rows" below.
    """
    index = {}
    if "sp" in arms:
        index.update({(pdb, mut): rec
                      for pdb, recs in load_bundled_store(key="skempi").items()
                      for mut, rec in recs.items()})
    if "mp" in arms:
        index.update({(pdb, mut): rec
                      for pdb, recs in load_bundled_store(MULTI_POINT).items()
                      for mut, rec in recs.items()})
    return {row: [float(index[row][t]) for t in TERMS] for row in rows if row in index}
```

`key="skempi"` on the single-point store is not optional. Without it about 29% of single-point
records are keyed in a form no SKEMPI-derived row list will match, and the join silently
under-covers — see [STORE.md](STORE.md).

## Single-point, multi-point, or both

```python
energies_for(rows, arms=("sp",))          # single-point only
energies_for(rows, arms=("mp",))          # multi-point only
energies_for(rows)                        # both -- 6003 entries indexed
```

**Using both is fine, and the warning elsewhere in these docs is about something else.**
`consolidate()` refuses to merge a single-point store with a multi-point one, and the two must
not share a directory — but that is about *building* stores. The two arms are separate BuildModel
campaigns whose per-complex files share a filename, so one directory cannot hold both without
each complex meaning whichever arm ran last.

Reading both and indexing them together, as above, is a different operation and is exactly right
for a model trained on single- and multi-point mutations together. The keys cannot collide: a
multi-point key contains a comma and a single-point key never does.

## Splits: by-complex, sequence-clustered, CATH-superfamily

The energies are independent of the split, which is the property that makes them useful — they
are computed once per `(complex, mutation)` from structure, then reused across every partition.
So there is nothing split-specific to do. Attach the energies to the rows, then partition; or
partition, then attach per fold. The result is identical either way.

The only thing a split changes is **which rows are uncovered**, and that is worth measuring per
fold rather than globally: a hold-out that happens to concentrate the uncovered rows will train
and test on different amounts of real signal.

```python
from skempi_foldx import coverage, load_bundled_store

covered, total, missing = coverage(load_bundled_store(key="skempi"), rows)
```

## Benchmark subsets (S1102, S1131, S2003, S4169)

Same shape: supply the subset's row list and join. The package cannot select a subset,
because it carries no benchmark definitions.

**Do not use `_source` to decide subset membership.** A record's `_source` names the *campaign
that computed the value* — `S4169`, `full_skempi`, or `full_skempi_multipoint` — which is
provenance, not membership. Campaigns computed overlapping and unioned mutation lists, and `consolidate()` picks
per mutation by source precedence, so a mutation's `_source` says nothing about which benchmarks
contain it. Take subset membership from the benchmark's own definition.

## Uncovered rows

`energies_for` omits rows it has no result for, rather than returning zeros. That is deliberate:
what a missing energy should become is a modelling decision with real consequences, and it is not
this package's to make.

Whatever is chosen, the two things worth doing:

- **Count them before training, not after.** Coverage is a property of the row list, so it is
  knowable up front.
- **Make "missing" distinguishable from "measured zero".** A fill value that also occurs
  naturally erases the distinction. Note that 296 single-point and 47 multi-point records are
  genuine all-zero vectors — non-interface mutations FoldX scores as having no effect — so zero
  is a real measurement here, not a safe sentinel.

## Interpretation

Every value is `IE(mutant) − IE(wild-type)` in kcal/mol, so negative means FoldX predicts tighter
binding. `Interaction Energy` is the summary term and the other eleven decompose it;
`skempi_foldx.TERMS` gives the canonical order for flattening a record into a feature vector.

Before comparing values across complexes, read [DETERMINISM.md](DETERMINISM.md) — a mutation's
energy depends on the whole mutation list it was computed in, and roughly three fifths of the
shipped single-point values came from subset lists rather than each complex's full list.
