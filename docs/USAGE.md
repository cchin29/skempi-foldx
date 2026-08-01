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
originating table. Being keyed, it also collapses a pair that appears twice, and `len()` of it is
therefore not a coverage count: on one 1389-row test fold, `coverage` reports 1385 covered where
`energies_for` returns 1379 entries. The six repeats there are `3SE4` rows carried once under each
of its two interface definitions — the collapse described below, not replicate measurements.

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
than guessing. The mapping files supply each chain's length, so the role-chain position is
computable rather than inferable.

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
worth running against a new setup. It lives in the repository rather than the installed package,
so it needs a checkout; `fx.coverage(rows)` gives the same figures from the package alone.

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

The blocks here show stdout. A `RuntimeWarning` about the interface collapse precedes it on
stderr on every run, because it is a property of the shipped store rather than of the rows given;
the section below on warnings says what each one means.

Inputs are read as headerless TSV with the sequence id in column 0 and the mutation in column 2
(`--pdb-col` and `--mut-col` override that); `--split` walks a directory tree of `*.tsv` and
de-duplicates, `--table` takes one flat file, both are repeatable, and `--json` writes the same
report as JSON. The uncovered line names the five complexes carrying most of the shortfall,
which distinguishes a diffuse gap from one complex missing entirely.

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

Rows of an arm that is not loaded simply count as uncovered: over the 5770-row combined split
above, a single-point instance covers 4121 rows, a multi-point instance 1614, and the default
5735. `require()` names that as the cause when the mutation's shape gives it away; for a complex
appearing only in the other arm it can say only that the complex is absent and which arms are
loaded. The coverage warning does not name it at all — it reports the shortfall and what would
widen the lookup, and leaves per-row causes to `require()`.

**Loading both arms into one lookup is correct. The prohibition on mixing them is about building
stores, not reading them.** `consolidate()` refuses to merge a single-point store with a
multi-point one, and the two must never share a results directory — see
[README](../README.md#running-the-pipeline). The two arms
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
partition rather than globally: a hold-out that happens to concentrate the uncovered rows will
train and test on different amounts of real signal. It does concentrate — over three folds of a
1634-row multi-point split, all 19 uncovered rows are `2C5D`, and they land wholly inside one
partition of each fold, taking a validation set to 88.4% covered against 100% for its siblings.
One `FoldxLookup` serves every fold.

`coverage_report.py --split` measures a whole split family, and pointing it at a single `fold_k/`
measures the same thing again: a fold's train, val and test files union back to the full dataset,
so each fold reports the same total. Per-partition figures come from `--table` on each file.

```python
fx = FoldxLookup(skempi_csv="skempi_v2.csv", mapping_dir="PDBs/")

for name, rows in folds.items():
    covered, total, missing = fx.coverage(rows)
```

## A store that did not ship with the package

The naming problem is the same for a campaign's own output, so `stores=` points the same
resolution at it — a results directory, an in-memory store, or several:

```python
fx = FoldxLookup(stores="my_campaign/results_sp")
fx = FoldxLookup(stores=["run_a/results_sp", "run_b/results_mp"])
```

The arm is read from the store rather than declared, because it is already recorded: on disk by
the key each file uses, in memory by whether any record name holds a comma. Nothing else changes —
`get`, `require`, `coverage` and the role-chain remap behave exactly as they do over the bundled
data, and `skempi_csv`/`mapping_dir` still widen what resolves.

The lower-level `coverage(store, wanted)` from `skempi_foldx.store` remains available for auditing
a raw `{pdb: {mutation: record}}` dict. It performs no naming resolution, so its figures over
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

## Warnings this package raises

Two conditions are announced as `RuntimeWarning` rather than left to be discovered in a result:
a lookup indexing records scored against the wrong interface, and a `coverage()` call where a
large share of rows failed to resolve and the instance could resolve more of them given more
inputs. Both are silent-wrong-answer conditions, so they interrupt rather than log.

**The coverage warning describes the lookup, not the rows.** It reports how many rows missed,
which conventions this instance matches, and which constructor argument would widen it — all
properties of the instance. It does not say why any particular row missed, because the label does
not carry the answer: a role-chain label and an author-form one that was never computed are the
same string for complexes whose author chains are themselves `A` and `B`, such as `1VFB` and
`1DQJ`. `require(pdb, mutation)` answers that question for a single row, where the record is
available to check against.

Three things silence it. Fewer than 100 rows; a shortfall under 5%; and `is_exact`, meaning every
mapping file this instance asked for has been read. That last is a statement about the mappings
asked for rather than about the rows — a `skempi_csv=` naming one complex is exact once that one
mapping is read, and an exact single-arm lookup is silent over every row of the arm it never
indexed.

A lookup built with `skempi_csv=` raises the interface warning **and** one per affected complex
from `load_skempi`. They are not duplicates. The per-complex ones describe the table supplied,
naming its pairings and how many mutations fall outside the kept one; the aggregate describes the
shipped store. A supplied table can collapse an interface the store's static registry knows
nothing about, and only the per-complex warning reports that.

Python's default filter shows each once per process. Under `-W error`, or a `filterwarnings =
error` setting in a test suite, constructing a `FoldxLookup` raises instead — deliberate, but not
always wanted in CI:

```python
warnings.filterwarnings("once", category=RuntimeWarning)
```

Filtering on `module="skempi_foldx"` does not work, and the reason is worth knowing: these
warnings are raised with `stacklevel=2` so the traceback points at the calling line rather than
at the package's own internals, and `module` is matched against that reported location. The
category is the reliable handle.

Suppressing them entirely restores the failure mode they exist to prevent, so narrowing to `once`
is the better trade.

## Complexes SKEMPI defines twice

Three PDB codes carry two different interface definitions in SKEMPI 2.0: `2C5D` (`A_C` and
`AB_CD`), `3SE3` and `3SE4` (both `B_A` and `B_C`). `load_skempi` keys on the code alone, so all of
a code's mutations are pooled under whichever definition its first row gave, and that is the
interface passed to `--analyseComplexChains`.

**31 shipped records therefore mutate a chain outside the interface they were scored against.**
Twenty are all-zero for that reason rather than as a measurement, and are a subset of the 343
(296 single-point and 47 multi-point)
all-zero vectors above. The other eleven are the dangerous ones: ten are multi-point variants of
`2C5D` with one component inside the analysed pair and one outside, so they carry a real energy
holding only the inside component's contribution — |IE| from 0.001 to 4.355 kcal/mol, one of
them negative — wrong without
looking wrong. The eleventh is `3SE4 YA55A`, a single-point record at 0.0002 — numerically an
artifact zero without being exactly zero. Discarding zeros is not sufficient protection; use
`interface_suspect()`. `load_skempi` warns for each affected complex, and
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
binding. `Interaction Energy` is the total. **The other eleven are components of it rather than a
decomposition of it**, and the difference matters: `AnalyseComplex` reports about twenty-one energy
columns, these are the eleven kept, and summing them does not reproduce the total — the residual
exceeds 0.1 kcal/mol on 47% of shipped records and reaches 4.6 on one. Dropping the total from a
feature vector therefore discards signal the components do not carry, and a consumer checking the
store by summation will wrongly conclude it is corrupt.

`skempi_foldx.TERMS` gives the canonical order, and `vector()` flattens a record into it.

Before comparing values across complexes, read [DETERMINISM.md](DETERMINISM.md) — a mutation's
energy depends on the whole mutation list it was computed in, and 2494 of the 4238 shipped
single-point values came from a campaign built around one benchmark rather than around each
complex's full list, 287 of them from a list that is a strict subset of it.
