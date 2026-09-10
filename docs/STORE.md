# The result store

One JSON per SKEMPI interface definition, holding every mutation computed against that pairing
and the 12 energy terms of each — the total and eleven of its components. This is the shipped
artifact, and the part that costs a FoldX licence (free for academic and non-profit institutions,
paid for commercial use) and CPU-weeks to reproduce.

```python
from skempi_foldx import load_bundled_store, coverage

store = load_bundled_store()              # {'1BRS_A_D': {'DA52A': {term: value}}, ...}
```

---

FoldX outputs are per-complex JSON: `{"muts": {"<mutation>": {<12 terms>, "cleaned": …}}, "meta":
…}`
(multi-point files use `"variants"` and key on the raw comma-joined SKEMPI mutation string).

Every value comes from one source, and there is nothing to arbitrate between sources. The repair
sweep's `round_1` supplies 4286 single-point and 1702 multi-point records — every code SKEMPI
defines once, computed against that code's union mutation list. The remaining 119 come from nine
per-definition campaigns covering `2C5D`, `3SE3` and `3SE4`, one pairing at a time.

That is a change of kind from 0.1.0, which consolidated five overlapping single-point campaigns and
had to record which one won each mutation and where the losers disagreed.
`skempi_foldx/data/CONSOLIDATION_REPORT.txt` still records provenance — a value is a property of
the list it was computed in — but it has no precedence order and no conflict section, because a
definition's mutations are computed together, in one list, against one interface.

The union is 322 single-point codes plus 152 multi-point, for 344 distinct codes overall, carried
as 323 and 154 interface definitions.

## Agreement with experiment

Coverage says how many rows carry a value; it says nothing about whether the values are right.
Against SKEMPI's own measured affinities, converted to ΔΔG as `RT·ln(Kd_mut / Kd_wt)` at each
row's recorded temperature:

| arm | n | Pearson r | Spearman ρ | MAE (kcal/mol) |
|---|---|---|---|---|
| single-point | 4081 | 0.418 | 0.437 | 1.14 |
| multi-point | 1580 | 0.518 | 0.587 | 1.87 |

**SKEMPI carries no ΔΔG column.** It records the two dissociation constants and a temperature —
`Affinity_mut (M)`, `Affinity_wt (M)`, `Temperature` — along with `dH` and `dS` for a subset of
rows. There is no free-energy field of any kind, so a ΔΔG compared against here is always derived
rather than read, and the derivation is a choice the comparing party makes. That is why the recipe
is stated rather than assumed, and why two papers quoting "SKEMPI ΔΔG" need not mean the same
numbers.

Reproducing these figures requires three such choices. Rows whose affinity is qualified rather than
numeric (`n.b`, `unf`, `>1E-04`) are dropped, and a mutation appearing on more than one SKEMPI row
is counted once, keeping the first — both move the result: keeping the qualified rows gives
n = 4190 and 1638, and keeping the last row instead of the first moves the single-point Spearman
to 0.439 and the multi-point one to 0.589. The third, 298 K standing in for a blank temperature,
changes no reported digit, because only one blank survives the other two filters.

The qualified rows are dropped rather than salvaged, which is a decision the table makes easy to
reverse and worth stating plainly. SKEMPI ships `Affinity_mut_parsed` and `Affinity_wt_parsed`
beside the raw fields, holding the number with its qualifier stripped: `>1.5E-04` becomes
`1.5E-04`. Where both forms are numeric they agree exactly, on all 6611 such rows, so the parsed
columns add nothing except those 187 rows — and what they add is an inequality recorded as a point
estimate. `>1.5E-04` means the complex binds more weakly than the assay resolves, so treating the
bound as the measurement understates how destabilising the mutation is — a bias with a known
sign, toward apparent agreement with any predictor evaluated on those rows. Either choice is
defensible; which one was made should be stated, since the two give different numbers.

A mutation is matched to its record by SKEMPI identifier, not by PDB code. That is what the third
decimal is sensitive to here, and it is why these are not the 0.1.0 figures: those were computed
against a store keyed by code, on values from a different set of mutation lists.

That is ordinary for FoldX on SKEMPI, and it sits alongside the published per-PPI figures quoted
in [DETERMINISM.md](DETERMINISM.md). The correlation is a property of FoldX rather than of this
store, and it is not offered as validation of the protocol. What it establishes is narrower, and
it is the check coverage cannot perform: the values behave like FoldX values, so a systematic
corruption between computation and packaging would have shown up here — an inverted sign above
all, which would leave every coverage figure unchanged.

## Coverage

`coverage(store, wanted)` answers "how many of these `(pdb, mutation)` pairs does the store have",
against an explicit wanted-set rather than an inferred one. It performs no naming resolution, so
the store must be keyed the same way as the wanted-set. Both are SKEMPI's author form here, so
they match directly — but a wanted-set of **role-chain** labels, which is what a split file gives,
scores zero against the store and raises nothing. [FoldxLookup](USAGE.md) resolves the conventions
instead of requiring them to be matched by hand:

```python
from skempi_foldx import coverage, load_bundled_store, load_skempi, worklist_single_point

sk = load_skempi("skempi_v2.csv")
wanted = [(pdb, m) for pdb, muts in worklist_single_point(sk).items() for m in muts]
covered, total, missing = coverage(load_bundled_store(key="skempi"), wanted)
```

Against everything SKEMPI defines, the shipped stores are:

| | shipped | SKEMPI defines | | excluding `1KBH` |
|---|---:|---:|---|---|
| single-point | 4340 | 4343 | 99.93% | **100%** |
| multi-point | 1767 | 1850 | 95.51% | **100%** |

The shortfall has exactly one cause:

* **`1KBH` — 3 single-point and 83 multi-point entries.** `RepairPDB` does not terminate on it, so
  no structure-based method has these. Not recoverable at any compute budget. See below.

The denominators count per interface definition — 324 single-point definitions and 155
multi-point, over 323 and 153 distinct PDB codes — because that is what SKEMPI names a row by.
Quote them with the `1KBH` exclusion stated; a bare "of 1850" implies a target that does not exist.

Two shortfalls that 0.1.0 carried are gone. **80 single-point entries were never computed**, since
the campaigns contributing those complexes ran against per-dataset mutation subsets; the rebuild
uses `worklist_single_point()` throughout. And **16 entries the chain mapping declined to
resolve** — eight complete pairs carrying the same substitution on two chains of the *same* group,
so the role letter could not separate them — are present now, because the store is keyed by the
author form the campaign computed rather than by a role form that had to be inferred.

How a role-chain label resolves onto these records is
[USAGE.md](USAGE.md#naming-conventions).

## Composition

Completeness is not the only shape these numbers have. The store is extremely uneven, and any
decision to split, weight or average over it is a decision about that unevenness:

![Single-point definitions are mostly pale and multi-point mostly red; a handful of tiles occupy a
third of the single-point panel and more of the multi-point one](assets/store_composition.png)

| | definitions | records | held by the largest 10 | holding a single record |
|---|---:|---:|---:|---:|
| single-point | 323 | 4340 | 34% | 105 |
| multi-point | 154 | 1767 | 42% | 46 |
| both arms, summed | 477 | 6107 | 25% | 151 |

**The last row sums the two arms; it does not deduplicate them.** 130 interface definitions appear
in both arms, so the store holds **347 distinct definitions**, not 477, and 477 counts the tiles
below — one per *(arm, definition)*. Records do not double-count, so 6107 is a true total; the
definition columns are per-arm sums. By distinct definition, **97** hold exactly one record rather
than 151.

**336 records — 292 single-point and 44 multi-point, 5.5% of the store — are all-zero on every
term.** These are mutations away from the interface, which FoldX scores as having no effect on
binding, so zero is a real measurement rather than a missing one. They enter the agreement figures
below unweighted, and anyone fitting on this column, ranking candidates, or choosing a fill value
should decide deliberately what to do with them — see
[USAGE.md](USAGE.md#coverage-of-an-evaluation-set).

Six single-point definitions hold at least 100 records each; a third of the 323 hold exactly one. A
definition holding one record supports no within-definition statistic at all — no correlation, no
rank, no spread — so the per-structure metrics common in this literature drop it. A figure
averaged over definitions and the same figure averaged over records are therefore answering
different questions, and on this store the gap is large.

The tile is an interface definition rather than a PDB code because that is the unit FoldX was
asked about, and because `2C5D`, `3SE3` and `3SE4` each carry two definitions. Keyed by code those
collapse into one tile or draw as several sharing a label; keyed by definition they are
`3SE4_B_A` and `3SE4_B_C`, and the picture says which pairing holds what.

Colour is the fraction of a definition's records that moved between 0.1.0 and 0.2.0, which is why
the single-point panel is pale and the multi-point panel is not. Tiles land at the ends rather
than in the middle — 360 of 476 paired tiles moved in no part, 57 moved entire, 59 split.

That is a statement about tiles, and it does not survive translation to definitions: a definition
present in both arms can be untouched in one and wholly moved in the other. Counted per distinct
definition the same data gives 239 still, 9 moved entire and 99 split, so the ends are far less
crowded. The tile is the right unit here because a mutation list belongs to an *(arm, definition)*
pair, which is the
visual form of the finding [VERSIONS.md](VERSIONS.md) measures: whether a record changed is a
visual form of the finding [VERSIONS.md](VERSIONS.md) measures: whether a record changed is a
property of its definition rather than of its mutation.

`experiments/store_composition.py` writes the figure, a hoverable version, a per-definition CSV,
and the full report at [`assets/STORE_COMPOSITION.txt`](assets/STORE_COMPOSITION.txt):

```
pip install -e ".[figures]"          # --plots needs matplotlib, --html needs plotly
git worktree add ../skempi-foldx-v0.1.0 v0.1.0
python experiments/store_composition.py --old ../skempi-foldx-v0.1.0 \
    --out docs/assets --plots --html --csv
```

`--old` is optional and only supplies the colour; without it the tiles are coloured by size class
and everything else on this page still reproduces. Nothing else is needed — no SKEMPI table, no
structures, no FoldX binary.

## One key convention, and the role name beside it

Every shipped key is SKEMPI's author-chain form — 4340 single-point and 1767 multi-point records
where `key == cleaned`, with no exceptions. 0.1.0's single-point arm split 3012/1226 between
that form and the
role-chain form, because it was assembled from campaigns run in different modes; the rebuild runs
every campaign in one mode, so the split is gone.

That alone would have stranded the labels a split file actually carries. Splits name the partners
`A`/`B` and renumber along each group, so the same substitution reads `LB38D` where the store says
`LI38D`, and a dictionary join on the role form now matches **only where the two coincide** —
silently, still writing
well-formed output. Rewriting a store to a single convention leaves such a consumer matching only
the subset where the two spellings coincide — a third of the rows in the split this was measured
on — with no error raised.

So a record whose role-chain name differs from its key ships that name in a **`role`** field:

```python
store = load_bundled_store()
rec = store["1ACB_E_I"]["LI38D"]
rec["cleaned"], rec["role"]        # ('LI38D', 'LB38D')
```

`role` is computed at build time from the SKEMPI `.mapping` files, by residue offset. 2763 of
4340 single-point records and 1124 of 1767 multi-point carry one — 3887 in all. Asked by
identifier, they recover all 1226 of the role-chain keys 0.1.0 shipped (by bare code, 1220: the
other 6 are `3SE4`, which is ambiguous by code however exactly it is named). The remaining 2220
need no field: their chain already sits first in its group, so its offset is 0 and the role form
*is* the key.

Because the field is built from chain *lengths*, a mapping file that parses short is a silent
corruption of it. That is not hypothetical: the files are fixed-width, so an author number of 1000
or more runs into the chain letter (`ALA C1001  1`) and a whitespace split yields three fields
instead of four. Requiring four dropped those residues — `2NYY`/`2NZ9` chain A measured 971 rather
than 1267, and `1S0W`/`1XXM` lost the whole of chain C — costing 16 records their role name.
`load_chain_mapping` now splits the field back apart. Recomputing all 6107 names with it
reproduces the shipped store exactly, so the store and its build script agree.

A field rather than a second key convention, deliberately. Deriving the role form at read time
needs `skempi_v2.csv` or the `.mapping` files, and this package redistributes neither — so a
consumer joining role-labelled rows would otherwise need the one input they cannot have.
`FoldxLookup` reads the field, so both names resolve with no configuration.

`reindex_by_skempi_id(records)` is an identity map on this store, and stays because it is not one
on a foreign store: a campaign run in role mode still keys its output that way.

## Intractable complexes

One complex, **1KBH**, is excluded from FoldX compute outright: RepairPDB does not terminate on
it. It has been left running past 22.5 h and past 10 h, on two different machines, without
finishing.

Curation already dropped 1KBH — RDE-Network's `block_list` does, so it appears in no split — but
curation
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

- The repair sweep saves every round it computes, and this store is assembled from `round_1` —
  the first pass — so the repair count is a property of which round was selected, not of a
  builder's flag.
- **The two inputs ran on different CPU architectures.** The sweep ran on Linux (4286
  single-point and 1702 multi-point records); the nine per-definition campaigns ran on Apple
  Silicon (the remaining 119). FoldX makes no statement about cross-platform bitwise agreement, so
  the two are not assumed interchangeable — see the measurement below, and
  [DETERMINISM.md](DETERMINISM.md) for what was and was not established.
- The nine per-definition campaigns *seeded* from an existing `<pdb>_Repair.pdb` rather than
  repairing again, so they inherited the same single repair. Each record's
  `meta.repair_seeded_from` names the campaign, complex and file it came from.
- Those seeds were verified byte-identical to the repairs the earlier campaigns built before any
  of them ran, so repair count is held fixed across both inputs rather than assumed to be.

**`meta.mutation_list_sha256` was added after the fact for the sweep's records.** The runner now
writes the hash of `individual_list.txt` as it runs, but the sweep's rounds were computed before it
did, so `experiments/repair_sweep/list_hashes.py` hashed the consumed lists afterwards and injected
the key. The 468 of 477 per-complex files it annotated say so in `meta.mutation_list_sha256_note`,
and that note names the commit that taught the runner to write the key ("predates dc5de42"). The
commit belongs to the history before this repository was published and is not resolvable here;
the hash it refers to is still the same quantity — the bytes BuildModel consumed — computed the
same way, which is why `experiments/repair_ablation.py compare` verifies these records natively.

That the two inputs are on one footing — the same repair count, and in this instance the same
answers across the architecture split — is measured rather than argued: `3SE3`'s `B,A` multi-point
arm exists in both, computed from different repair chains, and its 30 records agree to the last
decimal on an identical mutation-list hash.

**A record carries `_source`, `cleaned` and, where one exists, `role` — and nothing about how it
was computed.** The repair count, the FoldX version and the platform are properties of this
release, stated here, not fields you can carry with a subset. Anyone slicing this store into their
own table, or merging it with another campaign's output, has to carry that context themselves.

These values stay fixed and available even if a later store improves on them, because a result
computed from this store can only be reproduced against this store. A mutation's energy depends
on the whole list it was computed in, so a rebuilt store is a different measurement, not a
correction of the same one.

Ten source labels appear across the store — `sweep_4x_round_1` and the nine per-definition
campaigns — and `_source` on every record says which. They do not overlap: each record has exactly
one source, so there is no pair of values here to treat as corroborating or contradicting each
other. Where a value from *another* campaign disagrees with one of these, the difference is the
mutation-list effect described in [DETERMINISM.md](DETERMINISM.md), not measurement noise.

## Two common misreadings

**Chain letters are positional, not authorial.** Mutations are keyed within the SKEMPI partner
group (A = first chain, B = second), which does not always match the author chain letters in the
PDB. `1ACB LB38D` and `1ACB LI38D` are the same physical residue under the two conventions.

**Coverage checks cannot validate values.** An audit that compares mutation *keys* verifies that
nothing was dropped and nothing about the numbers behind those keys. That blind spot is wide
enough for a real class of defect: a join that gains a second candidate key can leave
coverage identical while handing a dozen or so mutations another mutation's ΔΔG, visible only as
a fall in downstream accuracy.
Any change to the join should be checked by diffing the 12-term values per
`(pdb, mutation)`, not the key sets.
