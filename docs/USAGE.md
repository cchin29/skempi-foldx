# Attaching these energies to an evaluation set

This package answers one question: **for this `(pdb, mutation)`, what did FoldX compute?** It
ships no splits, no benchmark definitions, and no train/test partitions, and it never will —
those belong to whatever is being evaluated, not to the energies.

So every use follows the same shape: bring a list of rows, get their energies, and decide what to
do about the rows that have none.

## The join

`FoldxLookup` is the entry point for all of it.

```python
from skempi_foldx import FoldxLookup

fx = FoldxLookup(skempi_csv="skempi_v2.csv", mapping_dir="PDBs/")

fx.get("1ACB", "LI38D")            # the 12-term record, or None
fx.require("1ACB", "LI38D")        # the same, but raises with the reason instead of None
fx.vector("1ACB", "LI38D")         # the record flattened into TERMS order
fx.arm_of("1ACB", "LI38D")         # "sp" or "mp", or None when uncovered
fx.energies_for(rows)              # {(pdb, mutation): [12 floats]} for the rows that resolve
fx.coverage(rows)                  # (covered, total, missing)
```

`rows` is any iterable of `(pdb, mutation)` pairs. `energies_for` keys its result by the row **as
given**, not by whatever key the store used internally, so the result joins straight back onto the
originating table.

`skempi_csv` is SKEMPI 2.0's `skempi_v2.csv`; `mapping_dir` is the directory of `<code>.mapping`
files distributed with SKEMPI's cleaned PDBs. Both are optional, and the next section is about
what each one buys.

## Naming conventions

**A mutation has more than one name, and that is where wrong answers come from.** The same
substitution appears as `LI38D` in SKEMPI's author-chain form and as `LB38D` in the role-chain
form a split file uses, and the shipped single-point store contains *both* — 3012 of its 4238
records under the first convention, 1226 under the second. Matching strings against one convention
therefore misses about 29% of the store without erroring.

`FoldxLookup` resolves the naming. The two optional constructor arguments widen what it can
resolve, and it reports which forms it will match rather than quietly matching fewer:

| constructed with | resolves |
|---|---|
| nothing | stored keys, and SKEMPI's author form via each record's `cleaned` field |
| `+ skempi_csv` | also role-chain forms, by substitution identity — exact for most complexes, unresolved where both partners admit the same substitution |
| `+ mapping_dir` | also role-chain forms **exactly**, by residue offset along each chain group |

```python
fx.conventions
# ['stored key', 'SKEMPI author form', 'role-chain form (exact, via residue offsets)']
```

The effect is large and it is measurable per row list. Over a 4159-row single-point split whose
labels are role-chain: stored keys and author forms alone resolve 65.4%, and adding both arguments
resolves 99.6%. Over a 1634-row multi-point split: 75.8% by identity matching alone against 98.8%
with the offset remap. **Supply both whenever row labels come from split files** — without them a
coverage figure over such labels is a lower bound, not a measurement.

Identity matching cannot reach the exact answer on its own: where both partners admit the same
substitution it has no basis for choosing, and `FoldxLookup` leaves such a row unresolved rather
than guessing. The mapping files supply each chain's length, which is what makes the role-chain
position computable rather than inferable.

`fx.require(pdb, mutation)` raises instead of returning `None`, and the error names how many
mutations the complex does hold, gives examples of their keys, and names the constructor argument
that would widen resolution:

```
KeyError: 1ACB has no record named 'LB999D'. It holds 6 mutations, e.g. ['LB38D',
'LB38E', 'LB38G']. This instance resolves only stored keys and SKEMPI author forms;
pass skempi_csv= and mapping_dir= to also resolve role-chain labels.
```

A miss is usually a naming mismatch rather than absent data, so `require` is the right call
wherever a miss should stop the run.

## Coverage of an evaluation set

Coverage is a property of the row list, so it is knowable before any training happens.
`experiments/coverage_report.py` reads splits or flat tables and reports it — the first thing
worth running against a new setup.

```bash
python experiments/coverage_report.py \
    --split path/to/splits_cath_kfold --label cath_kfold \
    --skempi-csv skempi_v2.csv --mapping-dir PDBs/
```

```
resolving: stored key, SKEMPI author form, role-chain form (exact, via residue offsets)

evaluation set                                 rows  covered            SP     MP
----------------------------------------------------------------------------------
cath_kfold                                     5770     5735  99.4%   4121   1614
                                             uncovered 35 over 6 complexes: 2C5D(19), 3HFM(6), 1DAN(4), 1DVF(2), 1DQJ(2)
```

Inputs are read as headerless TSV with the sequence id in column 0 and the mutation in column 2
(`--pdb-col` and `--mut-col` override that); `--split` walks a directory tree of `*.tsv` and
de-duplicates, `--table` takes one flat file, both are repeatable, and `--json` writes the same
report as JSON. The uncovered line names the complexes carrying the shortfall, which is what
distinguishes a diffuse gap from one complex missing entirely.

The same run without `--skempi-csv` and `--mapping-dir` reports 57.4% rather than 99.4% on that
split, and says so:

```
figures are a LOWER BOUND for role-chain labels -- pass --skempi-csv and
--mapping-dir to make them exact
```

The script holds no resolution logic of its own — it reads rows, asks `FoldxLookup`, and formats
the answer — so its figures and a program's own `fx.coverage(rows)` cannot diverge.

## Single-point, multi-point, and both

The `arms` **constructor** argument selects which stores are indexed. It takes the store directory
names, exported as `SINGLE_POINT` and `MULTI_POINT`; `energies_for` and the other lookup methods
take no such argument, because the instance already carries the answer.

```python
from skempi_foldx import FoldxLookup, MULTI_POINT, SINGLE_POINT

FoldxLookup(arms=(SINGLE_POINT,))                 # 4238 records
FoldxLookup(arms=(MULTI_POINT,))                  # 1765 records
FoldxLookup()                                     # both -- 6003 records
```

`arms` names the stores (`"results_sp"`, `"results_mp"`) while `arm_of()` reports the shorter
`"sp"` / `"mp"`. Passing the short form to the constructor raises `ValueError` rather than
silently indexing nothing.

An arm that is not loaded is not a miss with a hint, it is simply absent: over the 5770-row
combined split above, a single-point instance covers 4121 rows, a multi-point instance 1614, and
the default 5735.

**Loading both arms into one lookup is correct, and the prohibition in [STORE.md](STORE.md) is
about something else.** `consolidate()` refuses to merge a single-point store with a multi-point
one, and the two must never share a directory — but that is about *building* stores. The two arms
are separate `BuildModel` campaigns whose per-complex files share a filename, so one directory
cannot hold both without each complex meaning whichever arm ran last.

Reading both and indexing them together is a different operation, and it is exactly right for a
model trained on single- and multi-point mutations together. The keys cannot collide: every one of
the 1765 multi-point keys contains a comma and none of the 4238 single-point keys does.

## Splits: by-complex, sequence-clustered, CATH-superfamily

The energies are independent of the split, which is the property that makes them useful — they are
computed once per `(complex, mutation)` from structure, then reused across every partition. So
there is nothing split-specific to do. Attach the energies to the rows, then partition; or
partition, then attach per fold. The result is identical either way.

The only thing a split changes is **which rows are uncovered**, and that is worth measuring per
fold rather than globally: a hold-out that happens to concentrate the uncovered rows will train
and test on different amounts of real signal. One `FoldxLookup` serves every fold.

```python
fx = FoldxLookup(skempi_csv="skempi_v2.csv", mapping_dir="PDBs/")

for name, rows in folds.items():
    covered, total, missing = fx.coverage(rows)
```

The lower-level `coverage(store, wanted)` from `skempi_foldx.store` remains the right tool for one
case `FoldxLookup` does not serve: auditing a results directory that is not the bundled store — a
freshly consolidated campaign, say — since `FoldxLookup` only ever indexes the shipped stores. It
takes a `{pdb: {mutation: record}}` dict and performs no naming resolution, so its figures over
role-chain labels are the lower bound described above. [STORE.md](STORE.md) covers that path.

## Benchmark subsets (S1102, S1131, S2003, S4169)

Same shape: supply the subset's row list and join. The package cannot select a subset, because it
carries no benchmark definitions.

**A record's `_source` is not subset membership.** It names the *campaign that computed the
value* — the shipped stores carry exactly three, `S4169`, `full_skempi`, and
`full_skempi_multipoint` — which is provenance. Campaigns computed overlapping mutation lists and
`consolidate()` picks per mutation by source precedence, so a mutation's `_source` says nothing
about which benchmarks contain it. Three of the four named subsets contributed no record at all:
`skempi_foldx/data/CONSOLIDATION_REPORT.txt` records the `S1131`, `S2003` and S1102 campaigns as
contributing zero mutations to the union, since higher-precedence sources already held every one
of them. Take subset membership from the benchmark's own definition.

## Uncovered rows

`energies_for` omits rows it has no result for, rather than returning zeros. That is deliberate:
what a missing energy should become is a modelling decision with real consequences, and it is not
this package's to make.

Whatever is chosen, the two things worth doing:

- **Count them before training, not after.** Coverage is a property of the row list, so it is
  knowable up front — see the coverage report above.
- **Make "missing" distinguishable from "measured zero".** A fill value that also occurs naturally
  erases the distinction. 296 single-point and 47 multi-point records are all-zero vectors —
  mutations away from the interface, which FoldX scores as having no effect on binding — so zero
  is a real measurement here, not a safe sentinel.

## Complexes SKEMPI defines twice

Three PDB codes carry two different interface definitions in SKEMPI 2.0: `2C5D` (`A_C` and
`AB_CD`), `3SE3` and `3SE4` (both `B_A` and `B_C`). `load_skempi` keys on the code alone, so all of
a code's mutations are pooled under whichever definition its first row gave, and that is the
interface passed to `--analyseComplexChains`.

**31 shipped records therefore mutate a chain outside the interface they were scored against, and
20 of those are all-zero for that reason rather than as a measurement.** They are a subset of the
343 all-zero vectors above. `load_skempi` warns for each affected complex, and
`SkempiComplex.mutations_outside_interface()` lists the mutations concerned:

```python
from skempi_foldx import load_skempi
entry = load_skempi("skempi_v2.csv")["3SE4"]
entry.groups                          # 'B,C'  -- the definition that won
entry.alternate_groups                # {('B', 'A')}
entry.mutations_outside_interface()   # the 12 single-point mutations on chain A
```

The same collapse accounts for the largest single block of uncovered rows in any evaluation set
built from split files: `2C5D` contributes 19 uncovered multi-point rows because the splits number
it under `AB_CD` while the store computed it under `A_C`, so no label can match. Excluding these
three complexes is the conservative choice for work sensitive to interface definition.

`experiments/recompute_alternate_interfaces.py --plan` resolves the affected mutations and the
pairing each should be computed against; without `--plan` it runs those campaigns, one results
directory per interface. 3SE3 is defined twice as well and needs no recomputation: every mutation
of its second definition also appears in the first.

## Interpretation

Every value is `IE(mutant) − IE(wild-type)` in kcal/mol, so negative means FoldX predicts tighter
binding. `Interaction Energy` is the summary term and the other eleven decompose it;
`skempi_foldx.TERMS` gives the canonical order, and `vector()` flattens a record into it.

Before comparing values across complexes, read [DETERMINISM.md](DETERMINISM.md) — a mutation's
energy depends on the whole mutation list it was computed in, and 2494 of the 4238 shipped
single-point values came from a subset list rather than each complex's full list.
