# Attaching these energies to an evaluation set

This package answers one question: **for this `(pdb, mutation)`, what did FoldX compute?** It
ships no splits, no benchmark definitions, and no train/test partitions, and is not intended to —
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
originating table. Being keyed, it also collapses a pair that appears twice, so `len()` of it is
not a coverage count wherever a row list repeats a `(pdb, mutation)` pair. Use `coverage` for the
count and `energies_for` for the join.

`skempi_csv` is SKEMPI 2.0's `skempi_v2.csv`; `mapping_dir` is the directory of `<code>.mapping`
files distributed with SKEMPI's cleaned PDBs. Both are optional, and the next section is about
what each one buys.

## Naming conventions

**A mutation has more than one name, and that is where wrong answers come from.** The same
substitution appears as `LI38D` in SKEMPI's author-chain form and as `LB38D` in the role-chain
form a split file uses. Every shipped key is the author form, so a plain dictionary join keyed on
the role form matches **only the subset where the two coincide** — silently, still writing
well-formed output.

**For the shipped store, neither constructor argument is needed.** A record whose role-chain
name differs from its key carries that name in a `role` field, computed by residue offset when the
store was built; the rest need none, because their role form *is* the key. Either way
`FoldxLookup()` resolves both conventions out of the box:

```python
fx = FoldxLookup()
fx.conventions
# ['stored key', 'SKEMPI author form', 'role-chain form (exact, from the shipped `role` field)']
fx.is_exact        # True
```

Measured over the splits this store was built for, from the package alone:

| split | rows | covered | |
|---|---:|---:|---|
| single-point, by-complex | 4165 | 4165 | 100% |
| multi-point, by-complex | 1636 | 1636 | 100% |

**Ask by identifier, not by PDB code.** A split writes an interface as `1A22.A.B_A`, which is a
store file's name in different punctuation, so the conversion is
`field.split("_", 1)[0].replace(".", "_").upper()`. The `.upper()` is not decoration: store keys
are upper case and split files are not reliably so, and a whole-file case difference silently
halves the figure. Reducing the field to its code instead loses the 8 rows `3SE4` and `3SE3`
carry under **both** of their interface definitions — 6 and 2 — where a lookup by bare code has
two right answers and refuses to pick one. It also *conflates* them: the two pairings' rows
collapse onto one key, so the denominator falls too — which is why a by-code reading of these
two splits counts 4159 and 1634 rows rather than 4165 and 1636.

`experiments/coverage_report.py` does the conversion; read it there rather than reimplementing.

The two optional arguments still matter for a **foreign** store — one this package did not build,
whose records carry no `role` field:

| constructed with | resolves |
|---|---|
| nothing | stored keys and SKEMPI's author form via `cleaned`; no role forms, since a foreign store ships no `role` field |
| `+ skempi_csv` | also role-chain forms, by substitution identity — exact for most complexes, unresolved where both partners admit the same substitution |
| `+ mapping_dir` | also role-chain forms **exactly**, by residue offset along each chain group |

Identity matching cannot reach the exact answer on its own: where both partners admit the same
substitution it has no basis for choosing, and `FoldxLookup` leaves such a row unresolved rather
than guessing. The mapping files supply each chain's length, so the role-chain position is
computable rather than inferable — which is the computation the shipped `role` field already
records.

`fx.require(pdb, mutation)` raises instead of returning `None`, and the error names how many
mutations the complex does hold, gives examples of their keys, and says whether any input would
widen resolution:

```
KeyError: 1ACB has no record named 'LB999D'. It holds 6 mutations, e.g. ['LI38D',
'LI38E', 'LI38G']. This instance already resolves role-chain labels exactly, from
the shipped `role` field, so the name is absent rather than filed differently.
```

Over a store with no `role` field the last sentence is replaced by the constructor argument that
would widen it — `skempi_csv=` and `mapping_dir=`.

A miss is usually a naming mismatch rather than absent data, so `require` is the right call
wherever a miss should stop the run.

## Coverage of an evaluation set

Coverage is a property of the row list, so it is knowable before any training happens.
`experiments/coverage_report.py` reads splits or flat tables and reports it — the first thing
worth running against a new setup. It lives in the repository rather than the installed package,
so it needs a checkout; `fx.coverage(rows)` gives the same figures from the package alone.

```bash
python experiments/coverage_report.py \
    --split path/to/splits_bycomplex_all --label bycomplex_all
```

```
resolving: stored key, SKEMPI author form, role-chain form (exact, from the shipped `role` field)

evaluation set                                 rows  covered            SP     MP
----------------------------------------------------------------------------------
bycomplex_all                                  5801     5801 100.0%   4165   1636
```

A shortfall would print beneath the row, naming the complexes it falls on. Nothing prints here
because there is none.

The blocks here show stdout. A `RuntimeWarning` about the two malformed `2C5D` rows precedes it on
stderr on every run, because it is a property of the shipped store rather than of the rows given;
the section below on warnings says what each one means.

Inputs are read as headerless TSV with the sequence id in column 0 and the mutation in column 2
(`--pdb-col` and `--mut-col` override that); `--split` walks a directory tree of `*.tsv` and
de-duplicates, `--table` takes one flat file, both are repeatable, and `--json` writes the same
report as JSON. The uncovered line names the five complexes carrying most of the shortfall,
which distinguishes a diffuse gap from one complex missing entirely.

`--skempi-csv` and `--mapping-dir` change nothing here: the shipped `role` field already resolves
these labels exactly, and the run above omits both. They matter for a store this package did not
build, where their absence is announced rather than left to be inferred:

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

FoldxLookup(arms=(SINGLE_POINT,))                 # 4340 records
FoldxLookup(arms=(MULTI_POINT,))                  # 1767 records
FoldxLookup()                                     # both -- 6107 records
```

`arms` names the stores (`"results_sp"`, `"results_mp"`) while `arm_of()` reports the shorter
`"sp"` / `"mp"`. Passing the short form to the constructor raises `ValueError` rather than
silently indexing nothing.

Rows of an arm that is not loaded simply count as uncovered: over the 5801-row combined split
above, a single-point instance covers 4165 rows and a multi-point instance 1636, while the default
covers all 5801. `require()` names that as the cause when the mutation's shape gives it away; for
a complex
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
the 1767 multi-point keys contains a comma and none of the 4340 single-point keys does.

## Splits: by-complex, sequence-clustered, CATH-superfamily

The energies are independent of the split, which is the property that makes them useful — they are
computed once per `(complex, mutation)` from structure, then reused across every partition. So
there is nothing split-specific to do. Attach the energies to the rows, then partition; or
partition, then attach per fold. The result is identical either way.

The only thing a split changes is **which rows are uncovered**, and that is worth measuring per
partition rather than globally: a hold-out that happens to concentrate the uncovered rows will
train and test on different amounts of real signal. On the splits this store was built for there
is nothing left to concentrate — a 1636-row multi-point split resolves in full across all three
folds. Measure it anyway on a new split: the property is of the partition, not of the store, and
a store with a larger gap would put it somewhere. One `FoldxLookup` serves every fold.

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
fx = FoldxLookup(stores="campaign/results_sp")
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

These are the named row subsets of SKEMPI that this literature evaluates on, each defined by a
different paper and each identified by its row count; [REFERENCES.md](REFERENCES.md) gives the
sources. Same shape as above: supply the subset's row list and join. The package cannot select a
subset, because it carries no benchmark definitions.

**A record's `_source` is not subset membership.** It names the *campaign that computed the
value*, which is provenance. The shipped store carries ten labels — `sweep_4x_round_1` for
everything whose code SKEMPI defines once, and one per per-definition campaign for `2C5D`, `3SE3`
and `3SE4` — and **none of them is a benchmark name**. No `S1102`, `S1131`, `S2003` or `S4169`
value survives in 0.2.0: every one was recomputed against its complex's union mutation list, which
is the point of the release. Take subset membership from the benchmark's own definition.

## Uncovered rows

`energies_for` omits rows it has no result for, rather than returning zeros. That is deliberate:
what a missing energy should become is a modelling decision with real consequences, and it is not
this package's to make.

Whatever is chosen, the two things worth doing:

- **Count them before training, not after.** Coverage is a property of the row list, so it is
  knowable up front — see the coverage report above.
- **Make "missing" distinguishable from "measured zero".** A fill value that also occurs naturally
  erases the distinction. 292 single-point and 44 multi-point records are all-zero vectors —
  mutations away from the interface, which FoldX scores as having no effect on binding — so zero
  is a real measurement here, not a safe sentinel.

## Warnings this package raises

`FoldxLookup` announces five conditions as `RuntimeWarning` rather than leaving them to be
discovered in a result (the compute path raises its own, for an excluded complex and for a
malformed row reaching a worklist):

- a lookup indexing records whose SKEMPI row is malformed upstream;
- a `coverage()` call where a large share of rows failed to resolve and the instance could
  resolve more of them given more inputs;
- a name resolving to more than one record — across the stores supplied, or between two records
  of one complex claiming the same alias, which a single store can produce on its own. The first
  is kept, matching `consolidate()`. `FoldxLookup.conflicts` holds them, one entry per additional
  claimant, for a consumer that has narrowed the warning away;
- a PDB code indexed both as a record key and as interface definitions, which is what mixing a
  code-keyed store with a definition-keyed one produces. `FoldxLookup.pooled_codes` holds them.
  A bare-code lookup then returns the pooled record rather than the pairing-specific one, and the
  ambiguity refusal cannot fire because the code resolves exactly — so this one is reported
  rather than resolved, and the remedy is to load one store convention at a time;
- a `mapping_dir=` whose residue numbering disagrees with the `role` name a record already ships.
  The shipped name wins and the disagreements are listed in `FoldxLookup.role_mismatches`, because
  a recomputed role name that differs is a name belonging to some *other* record, and installing
  it as an alias is how a lookup returns the wrong energies while reporting `is_exact`.

All five are silent-wrong-answer conditions, so they interrupt rather than log.

**`coverage()` reports an ambiguous row as missing.** A `(code, mutation)` naming two interface
definitions cannot be resolved without the caller saying which pairing it wants, so it is counted
as uncovered and returned in `missing` — the same shape as a row that is genuinely absent. Eight
rows in the shipped store are ambiguous this way. `require()` distinguishes them, raising a
`KeyError` that names both identifiers; `definitions_of(code)` lists them. If a coverage shortfall
is small and unexplained, run `require()` over the missing rows before concluding the records do
not exist.

**The coverage warning describes the lookup, not the rows.** It reports how many rows missed,
which conventions this instance matches, and which constructor argument would widen it — all
properties of the instance. It does not say why any particular row missed, because the label does
not carry the answer: a role-chain label and an author-form one that was never computed are the
same string for complexes whose author chains are themselves `A` and `B`, such as `1VFB` and
`1DQJ`. `require(pdb, mutation)` answers that question for a single row, where the record is
available to check against.

Three things silence it. Fewer than 100 rows; a shortfall under 5%; and `is_exact`, meaning this
instance can resolve every role name it might be asked for. Two situations satisfy that. The
bundled store does, with no configuration and no file read, because every record whose role form
differs from its key ships its own
`role` field — which is why `FoldxLookup()` alone reports `is_exact` True. A lookup given
`mapping_dir=` instead satisfies it once every mapping file it asked for has been read, and that is
a statement about the mappings asked for rather than about the rows: naming one complex is exact
once that one mapping is read, and an exact single-arm lookup is silent over every row of the arm
it never indexed.

`load_skempi` does not warn about a second interface definition, because it discards neither:
each definition is its own entry, keyed by SKEMPI identifier. The warning that remains on
construction is about the two malformed `2C5D` rows, which is a property of SKEMPI's table rather
than of this store.

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
`AB_CD`), `3SE3` and `3SE4` (both `B_A` and `B_C`). Each definition is its own record, named by
the SKEMPI identifier, and each holds values computed against its own pairing:

```python
from skempi_foldx import FoldxLookup
fx = FoldxLookup()
fx.definitions_of("3SE4")        # ['3SE4_B_A', '3SE4_B_C']
fx.definitions_of("1BRS")        # ['1BRS_A_D']
fx.get("3SE4_B_A", "DA117A")     # scored against B,A -- the pairing that mutation belongs to
```

A lookup by bare code still works wherever the code is defined once, which is every code but
these three. Where it is defined twice, a mutation both definitions carry has two right answers:
`get` returns its default and `require` raises naming both identifiers, rather than choosing.

**Consequence of pooling a code's definitions.** Pooling a code's mutations under one pairing
scores the other definition's mutations against an interface they are not part of, and the result
does not look wrong. A mutation lying wholly outside the analysed pair moves neither side of
`IE(mutant) − IE(wild-type)` and returns a clean zero — indistinguishable from a measurement of no
effect. A multi-point variant with one component inside and one outside returns a real energy
carrying only the inside component's contribution. Across the three codes, 26 of 0.1.0's records
were zero on all twelve terms; recomputed against their own pairings, **15 of those 26 carry real
signal**, up to 6.85 kcal/mol:

| record | pooled | own pairing |
|---|---|---|
| `3SE4_B_A` `RA226A` | 0.0000 | 6.8535 |
| `2C5D_AB_CD` `LC105E,LD105E` | 0.0000 | 5.1289 |
| `2C5D_AB_CD` `KC178E,KD178E` | 0.0000 | 3.4407 |

`experiments/recompute_alternate_interfaces.py --plan` resolves the campaigns; `--all-definitions`
computes every definition of every multiply-defined code, which is what this store is built from.
Each pairing needs its own results directory: a campaign writes one file per record, and two
pairings of one code would otherwise overwrite each other.

## Interpretation

Every value is `IE(mutant) − IE(wild-type)` in kcal/mol, so negative means FoldX predicts tighter
binding. `Interaction Energy` is the total. **The other eleven are components of it rather than a
decomposition of it**, and the difference matters: `AnalyseComplex` reports about twenty-one energy
columns; these are the eleven kept, and summing them does not reproduce the total — the residual
exceeds 0.1 kcal/mol on 48% of shipped records and reaches 4.56 on one. Dropping the total from a
feature vector therefore discards signal the components do not carry, and a consumer checking the
store by summation will wrongly conclude it is corrupt.

`skempi_foldx.TERMS` gives the canonical order, and `vector()` flattens a record into it.

Before comparing values across complexes, read [DETERMINISM.md](DETERMINISM.md) — a mutation's
energy depends on the whole mutation list it was computed in, and 2494 of 0.1.0's 4238
single-point values came from a campaign built around one benchmark rather than around each
complex's full list, 287 of them from a list that is a strict subset of it.
