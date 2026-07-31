# The result store

One JSON per complex, holding every mutation computed for it and the 12 decomposed energy terms
of each. This is the shipped artifact — the thing most users want and the thing that costs
a FoldX licence (free for academic and non-profit institutions, paid for commercial use)
and CPU-weeks to reproduce.

```python
from skempi_foldx import load_bundled_store, coverage, audit

store = load_bundled_store()              # {pdb: {mutation: {term: value}}}
```

---


FoldX outputs are per-complex JSON: `{"muts": {"<mutation>": {<12 terms>, "cleaned": …}}, "meta":
…}`
(multi-point files use `"variants"` and key on the raw comma-joined SKEMPI mutation string).

These accumulated across five separate single-point campaigns — the S1102 set, three benchmark
sets (S1131 / S2003 / S4169) and the full-SKEMPI remainder — which is why the working tree carries
several result directories. A sixth campaign produced the multi-point arm. **The redundancy is
almost total**: `S1131`, `S2003` and the S1102 set contribute *zero* complexes that `S4169` does
not already contain, and the full-SKEMPI remainder is disjoint from all of them by construction.
The union is 322 single-point complexes, plus 152 multi-point, for 344 distinct complexes overall.
`skempi_foldx/store.py` consolidates them, and reports what each source actually contributed:
**three of the five single-point directories contribute nothing** — not one complex, not one
mutation — because they are subsets of the largest benchmark campaign.

## Coverage

`coverage(store, wanted)` answers "how many of these `(pdb, mutation)` pairs does the store have",
against an explicit wanted-set rather than an inferred one. It performs no naming resolution, so
the store must be keyed the same way as the wanted-set: `worklist_single_point` yields SKEMPI's
author form, which is what `key="skempi"` selects. Reading the store with the default
`key="stored"` and asking for author forms scores 3012 of 4337 rather than 4238 — the 69.4% figure
this file warns about below. [FoldxLookup](USAGE.md) resolves the conventions instead of requiring
them to be matched by hand:

```python
from skempi_foldx import coverage, load_bundled_store, load_skempi, worklist_single_point

sk = load_skempi("skempi_v2.csv")
wanted = [(pdb, m) for pdb, muts in worklist_single_point(sk).items() for m in muts]
covered, total, missing = coverage(load_bundled_store(key="skempi"), wanted)
```

Against everything SKEMPI defines, the shipped stores are:

| | shipped | SKEMPI defines | |
|---|---:|---:|---|
| single-point | 4238 | 4337 | 97.7% |
| multi-point | 1765 | 1848 | 95.5% |

The shortfall has exactly two causes, and they are worth separating because only one of them is
in principle recoverable:

* **`1KBH` — 3 single-point and 83 multi-point entries.** FoldX cannot repair it, so no
  structure-based method has these. Not recoverable. See below.
* **96 single-point entries never computed.** The campaigns contributing those complexes ran
  against per-dataset mutation *subsets* rather than each complex's full list. Recoverable by
  recomputing with `worklist_single_point()`.

Excluding `1KBH`, that is 97.8% of the reachable single-point set and **100%** of the reachable
multi-point set.

## Two key conventions in the single-point store

**A lookup by SKEMPI's mutation string misses 1226 of the 4238 single-point records.** The store
was assembled from campaigns that used two different key conventions:

| keyed by | entries | example |
|---|---:|---|
| SKEMPI's author-chain form — what SKEMPI calls `cleaned` | 3012 | `1ACB` → `LI38D` |
| the role-chain form, where the partners are named `A`/`B` | 1226 | `1ACB` → `LB38D` |

Both name the same mutation. `store["1ACB"]["LI38D"]` raises `KeyError`; the value is under
`LB38D`. Every record carries the author-chain form in its `cleaned` field, and no complex has
two records resolving to the same `cleaned` value, so the mapping is unambiguous in both
directions.

To look up by SKEMPI identifier, ask for that convention at read time:

```python
store = load_bundled_store(key="skempi")
store["1ACB"]["LI38D"]["Interaction Energy"]     # 3.6068
```

`reindex_by_skempi_id(records)` does the same for one complex's records.

**The stored keys are deliberately not rewritten on disk.** A consumer that maps role chains onto
author chains itself needs the role-chain key to be present. Rewriting the store to a single
convention was measured against one such consumer: its coverage fell from 3300/3300 to 1092/3300,
with no error raised — the join simply stops matching and still writes well-formed output. Choose
the convention when reading, not on disk.

The multi-point store has no such split: all 1765 keys are SKEMPI variant strings, and its
records carry no `cleaned` field because the key already is it.

## Intractable complexes

One complex, **1KBH**, is excluded from FoldX compute outright: RepairPDB does not terminate on
it. It has been left running past 22.5 h and past 10 h, on two different machines, without
finishing.

Curation already dropped 1KBH — RDE's `block_list` does, so it appears in no split — but curation
is not what stops it being *computed*. A driver that enumerates its scope from the **result
store** rather than from the curated table picks it up regardless: a stale `1KBH.json` (3
mutations) can sit in a source results directory even when the canonical store has no such file.
Because a sweep driver is a loop, killing the hung job only postpones it: the next round
enumerates, sees 1KBH, and launches RepairPDB again.

So the exclusion is applied **before the worklist is built**, in
[`skempi_foldx/exclusions.py`](../skempi_foldx/exclusions.py), and enforced at four points:

| Where | What it prevents |
|---|---|
| `run_campaign` | the complex reaching a worker at all |
| `process_complex` | a driver that builds its own targets and skips `run_campaign` |
| `store.consolidate` | a stale JSON being carried into the canonical store |
| any script that materialises a unified results dir | that dir re-acquiring it from a source dir |

Each drop is announced twice — on the driver's progress stream and as a `RuntimeWarning`, which
survives a driver that redirects stdout.

**Re-enabling.** The entry records the FoldX version the behaviour was observed against (5.1),
because "does not terminate" is a property of an engine, not a law. To test a newer one:

```bash
SKEMPI_FOLDX_ALLOW_INTRACTABLE=1 python experiments/repair_ablation.py run \
    --complexes 1KBH --iterations 1
```

If it converges, delete the entry from `INTRACTABLE` and rebuild the consolidated store. The
environment bypass deliberately does **not** re-admit the complex to curation — a curation step
reads `INTRACTABLE` as a plain set, not through the bypass — so a complex can be tested without
silently changing which rows every split contains.

## Provenance of the shipped store

**Every value here was computed with a single `RepairPDB` pass.** That is not recorded in the
per-record `meta`, so it is stated here with its evidence — it determines whether these numbers
are comparable to a repair-count series, and it is expensive to reconstruct after the fact:

- All six campaign builders — the five single-point campaigns above plus the multi-point one —
  invoke `RepairPDB` exactly once per complex, with no iteration construct anywhere in them.
- Of the 717 per-complex working directories those campaigns left, **717 carry exactly one
  `<pdb>_Repair.pdb` and none carries a `repair_round_<i>.pdb`** — the file an iterated chain
  necessarily writes. No `_original.pdb` or `_chain.pdb` scaffolding either.
- Campaigns that reused structures *seeded* from an existing `<pdb>_Repair.pdb` rather than
  repairing again, so they inherited the same single repair. There is no path by which an
  iterated structure entered this store.

The on-disk check covers the single-point campaigns. The multi-point working directory is not
retained, so for that arm the evidence is its builder's code alone.

These values stay fixed and available even if a later store improves on them, because a result
computed from this store can only be reproduced against this store. A mutation's energy depends
on the whole list it was computed in, so a rebuilt store is a different measurement, not a
correction of the same one.

Two campaigns contributed the shipped single-point values, and `_source` on every record says
which. They are not independent measurements of the same thing: where they overlap, the
difference between them is the mutation-list effect described in
[DETERMINISM.md](DETERMINISM.md), not measurement noise. Treating agreement between them as
corroboration would be a mistake — and treating disagreement as a bug would be a different one.

## Two common misreadings

**Chain letters are positional, not authorial.** Mutations are keyed within the SKEMPI partner
group (A = first chain, B = second), which does not always match the author chain letters in the
PDB. `1ACB LB38D` and `1ACB LI38D` are the same physical residue under the two conventions.

**Coverage checks cannot validate values.** An audit that compares mutation *keys* verifies that
nothing was dropped and nothing about the numbers behind those keys. That blind spot is wide
enough for a real defect: a join that gained a second candidate key left coverage identical while
handing ~13 mutations another mutation's ΔΔG, and surfaced only because downstream accuracy
*fell*. Any change to the join should be checked by diffing the 12-term values per
`(pdb, mutation)`, not the key sets.

