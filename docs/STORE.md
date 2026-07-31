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

These accumulated across five separate campaigns — the S1102 set, three benchmark sets
(S1131 / S2003 / S4169) and the full-SKEMPI remainder — which is why the working tree carries
several result directories. **The redundancy is almost total**: `S1131`, `S2003` and the S1102 set
contribute *zero* complexes that `S4169` does not already contain, and the full-SKEMPI remainder
is
disjoint from all of them by construction. The union is 322 single-point complexes, plus 152
multi-point, for 344 distinct complexes overall. `skempi_foldx/store.py` consolidates them, and
reports what each
source actually contributed: **three of the five directories contribute nothing** — not one
complex, not one mutation — because they are subsets of the largest benchmark campaign.

## Coverage

Coverage is the fraction of split rows that receive a real FoldX score rather than the imputed
train-mean. It is not 100%, and the residual has two distinct causes worth separating:

* **genuinely absent** — the complex or mutation was never computed;
* **denied by the grouping guard** — the merger refuses to score a row whose SKEMPI chain grouping
  it cannot match unambiguously, rather than guessing. `3SE4.B.A` is the recurring case.

Current single-point coverage is ~99%, and it got there almost entirely by *plumbing* rather than
by new compute: the climb from 28% to 58% to 88% to 99% was mostly a matter of pointing the
mergers at result directories that already existed on disk. A low coverage number is worth
auditing against what is actually being read before any CPU is spent on it.

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

- All seven campaign builders invoke `RepairPDB` exactly once per complex, with no iteration
  construct anywhere in them.
- Of the 717 per-complex working directories those campaigns left, **717 carry exactly one
  `<pdb>_Repair.pdb` and none carries a `repair_round_<i>.pdb`** — the file an iterated chain
  necessarily writes. No `_original.pdb` or `_chain.pdb` scaffolding either.
- The campaigns that reuse structures (`_bench`, `_delta`) *seed* from an existing
  `<pdb>_Repair.pdb` rather than repairing again, so they inherit the same single repair. There
  is no path by which an iterated structure entered this store.

The on-disk check covers the single-point campaigns; the multi-point working directory has since
been pruned, so for that arm the evidence is its builder's code alone.

This store is also the **provenance record for a published set of ΔΔG model results**: those
splits were merged from exactly these values. It stays available and citable for that reason,
independent of any later store that improves on it.

Worth knowing before treating any two directories as independent measurements:

- `results_all` contains **no original compute**. It is 211 byte-identical copies of the S4169
  campaign plus 112 of the full-SKEMPI campaign.
- `results_remainder` is **not** an independent recompute, despite reading like one: it is
  **byte-identical to the full-SKEMPI campaign on all 65 complexes they share**. The
  often-quoted "40 of 2473 mutations differ" comparison is therefore roughly half
  self-against-self; its non-identical half compares two campaigns that shared repaired
  structures and agree to *r* = 0.9991.

So the 1.6% figure and the 66.7% figure are both correct and not in conflict: they measure
different things. **66.7% is the honest answer to "what changes if the mutation lists change";
1.6% describes two directories that are largely the same bytes.**

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

