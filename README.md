# skempi-foldx

**Precomputed FoldX binding energies for SKEMPI 2.0 — a fully specified FoldX baseline for
mutation-effect prediction, computed once and shipped with the settings it was computed under.**

## Background

Mutating a residue at the interface between two proteins makes their binding stronger or weaker.
The size of that change is **ΔΔG of binding**, measured in kcal/mol, and predicting it is
central to antibody engineering, protein design, and interpreting disease variants.

Two things anchor that field:

- **[SKEMPI 2.0](https://life.bsc.es/pid/skempi2)** — the benchmark database of experimentally
  measured binding changes for thousands of mutations in protein–protein complexes. Nearly every
  ΔΔG method is evaluated on it.
- **[FoldX](https://foldxsuite.crg.eu/)** — a widely used empirical force field that estimates the
  same quantity from structure. Its number is the physics baseline that learned methods are
  expected to beat, and it is also a useful input feature for them.

Putting those together — running FoldX across all of SKEMPI — is something many groups need and
each one redoes. It takes a FoldX licence, correct handling of SKEMPI's chain conventions, and
CPU-weeks. **This repository is that computation, done once and shipped**, plus the pipeline that
produced it.

## Applications

- **Benchmarking a ΔΔG predictor.** A FoldX column is the conventional baseline to report against.
  This one is fully specified — protocol, repair count, mutation lists and known defects are all
  written down, so a reader can reproduce it or diverge from it deliberately.
  [docs/PROTOCOL.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/PROTOCOL.md) surveys
  what the comparator papers state about their own FoldX settings. It is not the same computation
  as any published FoldX row and should not be quoted as one.
- **Using FoldX as a feature.** The 12 energy terms per mutation — the total and eleven
  components — are a ready-made structural feature, independent of the encoder and of the
  evaluation split.
- **Running FoldX at dataset scale.** The pipeline handles SKEMPI parsing, chain-group mapping,
  resumption, and parallelism — and documents the failure modes that are expensive to discover
  independently.

**This is not a FoldX wrapper.** It does not wrap, bundle, or replace FoldX — for Python bindings
see [pyfoldx](https://github.com/leandroradusky/pyfoldx). This is a dataset-scale campaign runner
and the campaign's results.

## Quickstart

```bash
pip install skempi-foldx==0.2.2
```

**Pin the version.** Energies are tied to a release — [Scope and stability](#scope-and-stability)
says why, and
[docs/VERSIONS.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/VERSIONS.md) says what
moved since 0.1.0 and whether it affects a model already trained on it. An unpinned `pip install
skempi-foldx` resolves to whatever is newest, which for this package means a different measurement,
not just newer code.

0.2.2 is the first release on PyPI, and its energies are byte-identical to 0.2.0's — it is a
change of distribution, not of data. Every earlier release installs from its git tag instead:

```bash
pip install "git+https://github.com/cchin29/skempi-foldx@v0.1.0"
```

Or, from a clone of this repository:

```bash
pip install .
```

There are no dependencies outside the standard library, and the energies ship inside the wheel —
no separate data download follows. `pip install '.[dev]' && pytest` runs the suite, which pins the
store's shape, its energy values and the quickstart figures against the shipped data. It ends with
a summary of warnings — the malformed-row and excluded-complex notices below, raised deliberately
by the tests that assert on them — and that is expected, not a failure. Figures derived from SKEMPI's
table or from evaluation splits cannot be pinned here, because neither ships with this package.

A mutation is named as in SKEMPI: `LI38D` is Leu→Asp at position 38 of chain I — original
residue, chain, position, new residue. Energies are always mutant minus the unmutated
protein, the *wild type*.

```python
from skempi_foldx import FoldxLookup

fx = FoldxLookup()                            # both arms; len(fx) is 6107 records
fx.get("1BRS", "DA52A")["Interaction Energy"] # -0.6806, mutant − wild-type
fx.vector("1BRS", "DA52A")                    # 12 terms, in TERMS order

rows = [("1BRS", "DA52A"), ("1ACB", "LI38D")] # any (pdb, mutation) pairs
fx.coverage(rows)                             # (covered, total, missing)
fx.energies_for(rows)                         # {(pdb, mut): [12 floats]}
```

Constructing `FoldxLookup` prints a `RuntimeWarning` about two malformed rows. That is expected and
is not an installation problem: two SKEMPI rows name three substitutions while labelling themselves
as four, and they are computed and shipped rather than dropped, because substituting the intended
string would produce a record matching no SKEMPI row. `FoldxLookup.is_malformed()` tests for them
and `MALFORMED_ROWS` lists them; see
[docs/USAGE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/USAGE.md#warnings-this-package-raises).

## Naming, terms and provenance

`FoldxLookup` is the entry point rather than the raw dictionary because **a mutation carries more
than one name**. The same substitution is `LI38D` in SKEMPI's author-chain form and `LB38D` in the
role-chain form used by the benchmark *split files* — the train/test partition files this
literature ships, which rename each complex's two partners `A` and `B` regardless of their PDB
chain letters. Every shipped key is the author form, and a record whose role
name differs from its key carries that name in a `role` field — 3887 of 6107 records. For all
2220 of the rest the role form *is* the key, because the chain already sits first in its group, so
both names resolve with no configuration:

```python
# The author form and the role form reach the same record.
fx.get("1BRS", "DD35A") is fx.get("1BRS", "DB35A")
```

`DD35A` is one of the 3887: its key is the author form and its `role` field carries `DB35A`, which
is what makes the second lookup resolve. `DB35A` is not itself a key in the store.

A plain dictionary join keyed on the role form silently matches only the 2220 records whose
role form is already their key, and misses the other 3887. `fx.require()` raises
instead of returning `None`, and its error names the fix.

A bare PDB code resolves wherever it names one record. `2C5D`, `3SE3` and `3SE4` are defined
twice, but that only makes a lookup ambiguous for a mutation *both* definitions carry — **8 of the
111** mutations under those three codes. `2C5D`'s two definitions share none at all, so every one
of its records answers by code. For the 8, `get` returns its default and `require` raises naming
both identifiers rather than choosing an interface on the caller's behalf; ask by identifier
instead. `fx.definitions_of(code)` lists them.

The raw store stays available as `load_bundled_store()` for auditing and for consumers doing their
own chain mapping — see
[docs/STORE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/STORE.md), and
[docs/USAGE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/USAGE.md) for attaching
these energies to an evaluation set.

Each record carries all twelve terms from `AnalyseComplex`, FoldX's interface-energy command, as
mutant minus wild-type, plus bookkeeping — `_source`, `cleaned`, and `role` where the role-chain
name differs from the key. The record below is one of the 2220 where the two names coincide, so it
carries no `role`:

```python
from skempi_foldx import load_bundled_store

store = load_bundled_store()
store["1BRS_A_D"]["DA52A"]
# {'Interaction Energy': -0.6806,   'Backbone Hbond': 0.0,
#  'Sidechain Hbond': 0.0,          'Van der Waals': 0.0283,
#  'Electrostatics': -0.2432,       'Solvation Polar': -0.1739,
#  'Solvation Hydrophobic': 0.0086, 'Van der Waals clashes': -0.0001,
#  'entropy sidechain': -0.0381,    'entropy mainchain': 0.0,
#  'torsional clash': -0.0,         'backbone clash': 0.0,
#  'cleaned': 'DA52A',              '_source': 'sweep_4x_round_1'}
```

A negative interaction energy means FoldX predicts the mutant binds *more* tightly. For a single
number rather than twelve, use `Interaction Energy`. **The other eleven are components of it, not a
decomposition of it**: they do not sum to it, so a model given the components must keep
`Interaction Energy` as well, and a consumer checking the store by summation will wrongly conclude
it is corrupt. [docs/USAGE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/USAGE.md)
gives the size of the residual. `skempi_foldx.TERMS` gives the canonical order for flattening
records into a feature vector.

`_source` names the campaign that produced the value, which matters because a value is a property
of the mutation list it was computed in — see [running the pipeline](#running-the-pipeline).
**Every value comes from a union mutation list**: 4286 single-point and 1702 multi-point records
from the repair sweep's `round_1`, and the remainder from the nine per-definition campaigns that
computed `2C5D`, `3SE3` and `3SE4` one pairing at a time. 0.1.0 had 2494 of 4238 single-point
values from a per-benchmark campaign instead; that is what this release removed. How large the list
effect is, measured between campaign pairs, is in
[docs/DETERMINISM.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md).
`cleaned` is SKEMPI's own mutation string for the record.

## Contents and coverage

| | definitions | codes | entries | |
|---|---:|---:|---:|---|
| `load_bundled_store()` | 323 | 322 | 4340 | single-point mutations |
| `load_bundled_store("results_mp")` | 154 | 152 | 1767 | multi-point variants (two or more substitutions) |

Each row is one SKEMPI **interface definition**, named `<pdb>_<group1>_<group2>` — `1BRS_A_D`.
Three codes carry two definitions each, but not in both arms: only `3SE4` is doubled in the
single-point arm (323 − 322 = 1), and only `2C5D` and `3SE3` in the multi-point one (154 − 152 =
2).

Measured against what SKEMPI itself defines — 4343 single-point *(definition, mutation)* entries
over 324 definitions, 1850 multi-point over 155 — that is **4340 of 4343** and **1767 of 1850**.
Entries rather than mutations: a row under a doubled code is counted once per definition it
belongs to, so the shipped 4340 entries are 4334 distinct mutations and the 1767 are 1765. The gap does not need
decomposing, because it is one complex:

- **`1KBH`** accounts for all of it: 3 single-point and 83 multi-point. `RepairPDB`, FoldX's
  structure-preparation step, does not
  terminate on it under FoldX 5.1, so no amount of compute closes this and no method depending on
  a repaired structure will have these rows.

Excluding `1KBH`, **both arms are complete**. Quote the denominators with the exclusion stated; a
bare "of 1850" implies a target that does not exist. Completeness is about coverage, not about
every label being sound: two `2C5D` records carry an arity label that is wrong upstream in SKEMPI
and are shipped flagged rather than dropped — see `MALFORMED_ROWS`.

Against SKEMPI's measured affinities the shipped values reach Spearman **0.437** on the
single-point arm (4081 rows) and **0.587** on the multi-point one (1580) — ordinary for FoldX, and
the check coverage cannot perform. These are pooled over rows, so they are not comparable with the
per-PPI (per protein–protein complex) or per-structure correlations the comparator papers report —
those average over complexes, and a third of the single-point definitions here hold exactly one
record. For a figure that *is* comparable,
[docs/DETERMINISM.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md)
reports per-structure Spearman of **0.3984 [0.2759, 0.5193] on 13 complexes**, against the 0.4458
per-PPI published by CATH-ddG — a 2025 predictor whose FoldX baseline is the most fully specified
in this literature
([docs/PROTOCOL.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/PROTOCOL.md)). Below
it, on a sample far too small to separate the two.

Both figures are `RT·ln(Kd_mut/Kd_wt)` at each row's temperature — `Kd` is the dissociation
constant SKEMPI records for wild type and mutant, `R` the gas constant — over rows with numeric
affinities, deduplicated on `(identifier, mutation)` keeping the first row. The recipe matters in
both directions: keeping the last row instead moves the figures to 0.439 and 0.589, and keeping
the 187 rows whose affinities are qualified (`n.b`, no binding detected; `>1E-04`, weaker than the
assay resolves) would add 167 of them back after deduplication at a bias whose sign is known —
toward apparent agreement. Coverage, agreement and the key conventions
are detailed in [docs/STORE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/STORE.md).

## Scope and stability

Version 0.2.2. The energies are the artifact and they are tied to a version. A mutation's ΔΔG
depends on the whole mutation list it was computed in, so a rebuilt store is mostly a different
measurement rather than a correction of an earlier one — a published result quoting these numbers
should pin the version it installed — `skempi-foldx==0.2.2` — and say so. The git tag
`…/skempi-foldx@v0.2.2` names the same tree, and is the only way to name a release before this one.

**0.2.2 carries exactly 0.2.0's energies.** It is the first release distributed on PyPI and
changes no value; citing either is citing the same measurement.

**Two `v0.2.0` tags and a `v0.2.1` were published and withdrawn during August and early
September.** They behave differently, and the quiet one is the one to watch:

- **`v0.2.1` does not exist.** A pin naming it fails outright, which at least announces itself.
  0.2.1 is skipped, so that a surviving pin keeps failing rather than resolving to something
  unrelated.
- **`v0.2.0` exists and is a valid release**, but points at a tree the two withdrawn `v0.2.0`
  objects did not. A pin naming it still resolves — *silently, to different content than an
  August install had*. Check the commit below.

**This affects an install taken from this repository between 2026-08-02 and 2026-09-09** — from a
tag, a branch, or a clone. Nothing withdrawn ever reached PyPI, which had no release of this
package before 0.2.2. Every 0.1.0 install is unaffected.

A clone taken from one of the withdrawn trees shares no history with this one, which is a squashed
replacement. `git pull` will not reconcile them; re-clone, or
`git fetch origin && git reset --hard origin/main`.

No energy value differs between any of those trees and this one — 6107 records compared field by
field — so a result quoting the numbers stays correct. A reinstall adds the `role` name to 16
records that a mapping-file parsing bug had kept out of them, plus corrected documentation.

The commit an install was built from is what identifies it, and pip records it for a git install:

```bash
python -c "import importlib.metadata as m, json; \
           t = m.distribution('skempi-foldx').read_text('direct_url.json'); \
           d = json.loads(t) if t else {}; \
           print(d.get('vcs_info', {}).get('commit_id') or d.get('url') or 'installed from PyPI')"
```

Compare it against the commit this release's tag points at, which the repository will tell you
without a clone:

```bash
git ls-remote https://github.com/cchin29/skempi-foldx 'refs/tags/v0.2.2^{}'
```

Anything else predates this release. The command answers rather than raising whatever the install
is: `direct_url.json` is absent for a PyPI install, and for a wheel or sdist installed from a path
it records the archive it came from instead of a commit. A line that is not a commit hash is the
answer — this install did not come from a tag — not a failure to look.

*A release cannot state its own commit hash — the hash covers the file that would state it — which
is why this asks the remote instead of printing a literal.*

A second check, from the data rather than the metadata:

```python
FoldxLookup().get("1S0W", "AB142F")["role"]
```

`'AB142F'` on this release. A 0.1.0 install raises `KeyError: 'role'`, since 0.1.0 ships no `role`
field at all. A withdrawn tree carries `role` on 3871 records but lacks the name on the 16 this
one backfilled, of which `1S0W_A_C` is one — so the lookup returns `None` and this raises
`TypeError`. A typo raises `TypeError` too, so check the spelling before concluding anything from
a failure.

A withdrawn-tag install reports `0.2.0` or `0.2.1`, so installing 0.2.2 from PyPI moves off it
without needing `--force-reinstall` — the version differs, and pip will replace it:

```bash
pip install --upgrade skempi-foldx==0.2.2
```

To stay on a git install instead, pip needs telling, since `v0.2.0` still resolves and its version
string already matches:

```bash
pip install --force-reinstall "git+https://github.com/cchin29/skempi-foldx@v0.2.2"
```

The corrections to the tests, CI and contributor documentation are in the source tree, which a `pip
install` does not carry;
[CHANGELOG.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/CHANGELOG.md) itemises them. The
comparison below is against 0.1.0.

**1850 of 0.1.0's 6003 values changed, and 4145 did not**; the remaining eight could not be
compared, and
[docs/VERSIONS.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/VERSIONS.md) says what
they are and where the record-by-record list is. 26 of those changes are corrections: records under
`2C5D` and `3SE4`, which SKEMPI pairs two ways each, that 0.1.0 scored against an interface
excluding the chain they mutate. That failure is invisible in the output — a mutation outside the
analysed pair returns a clean zero, indistinguishable from a real measurement of no effect.
Separately, and confusingly the same count: across the three doubled codes, 26 of 0.1.0's records
were zero on all twelve terms. Call those the *clean-zero records*; they are not the same set as
the 26 corrections, though the two overlap in 15 members. Recomputing per pairing gives 15 of those
zeros real signal, up to 6.85 kcal/mol. **A 0.1.0 consumer whose rows include those records was
trained on values that are wrong.**

The other 1824 changes are not corrections. The store is rebuilt against union mutation lists
throughout, and a value computed in a different list is a different measurement, so records whose
0.1.0 list already matched are identical and the rest are not — neither release's value is the
more correct one. That makes the two *look* interchangeable across two thirds of the store while a
third of it disagrees, by a median of 0.208 kcal/mol and up to 9.96, which is why they must not be
mixed within one table.

Against experiment the two are equivalent on the records they share — Spearman 0.435 → 0.436
single-point and 0.574 → 0.587 multi-point, with paired bootstrap intervals spanning zero on the
records that changed — so **a model trained on 0.1.0 single-point features is not invalidated by
0.2.0**, and neither release is the better predictor. **A multi-point model is a different case**:
90% of that arm changed, at Spearman 0.944, MAE 0.469 and 112 sign flips, so retraining is what
establishes the size of that sensitivity rather than bounding it, and nothing here predicts which
way such a retrain moves. `experiments/compare_versions.py` writes one row per 0.1.0 record, which
settles the question for a particular subset before anyone commits to a retrain. It needs a
checkout of the older tree, and lives in the repository rather than in the distribution:

```bash
git worktree add ../skempi-foldx-v0.1.0 v0.1.0
python experiments/compare_versions.py --old ../skempi-foldx-v0.1.0 \
    --out scratch/version_compare --csv
```

As written that reproduces every section of
[`docs/assets/VERSION_COMPARE.txt`](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/assets/VERSION_COMPARE.txt)
except the two that score each release against measured affinities: those need `--skempi-csv
<skempi_v2.csv>`, which is not redistributable and so does not ship here. Without it the script
says so and omits them.

**v0.1.0 remains tagged, installable and citable.** Nothing in this release retracts it; a
published result built on it stays reproducible against it, provided the tag is the one cited.
[docs/VERSIONS.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/VERSIONS.md) gives the
record-by-record counts, which complexes moved and by how much, and the script that reproduces all
of it. See also [CHANGELOG.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/CHANGELOG.md).

`TERMS` and its order are stable, because a consumer's feature layout depends on them. Nothing
else in the API is settled yet.

Every value was computed with a **single** `RepairPDB` pass, matching CATH-ddG. FoldX's own
documentation says to repair before modelling without naming a count; reading that as a
recommendation of one is Usmanova et al.'s, and is attributed to them in
[docs/DETERMINISM.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md). It
is a protocol choice, not a calibrated optimum, and the evidence against treating it as settled is
stronger than a single pass suggests. A completed four-round sweep over all 344 complexes, with the
mutation list held byte-constant so repair count is the only variable, finds **no round after which
the store stops changing**: between rounds 3 and 4, 6.1% of single-point and 23.4% of multi-point
entries still move by more than FoldX's ~0.5 kcal/mol noise floor. The shipped store is round 1. A
smaller 13-complex pilot separately found 5× repair scored *better* against experiment
(per-structure Spearman 0.398 → 0.459), though it is too narrow to generalise from.
[docs/DETERMINISM.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md) has
both; read it before treating one pass as settled. The choice is not recorded per-record, so
[docs/STORE.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/STORE.md) states the
evidence for it — the store is assembled from the repair sweep's first round, the per-definition
campaigns were seeded from repairs verified byte-identical beforehand, and one complex computed
both ways agrees to the last decimal — and it determines what these numbers can legitimately be
compared against.

## Running the pipeline

This requires a FoldX 5 binary (`FOLDX_BIN`), SKEMPI's cleaned PDBs, and `skempi_v2.csv`. FoldX
is licensed software and is not redistributed here — free for academic and non-profit
institutions, paid for commercial use.

**Read
[docs/DETERMINISM.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md)
first.** The short version, because it is the failure mode most likely to cost a week:

FoldX `BuildModel`, the command that builds each mutant, is deterministic — same structure and
same mutation list, same answer to the
last decimal, over 1719 measured pairs. But a mutation's ΔΔG depends on the whole list it was
computed in, and not on the mutation alone. **Computing a subset of a complex's mutations gives
different numbers than computing all of them** — measured spreads
of median 0.067 over the pairs that differ — 0.008 over all pairs — and up to 11.03 kcal/mol
for the same mutation of the same structure.

A spread like that invites the reading "FoldX is stochastic, seed it". The measurements rule that
out, and acting on it means averaging noise that does not exist. The variable to control is the
mutation list:
`worklist_single_point()` takes the union by construction.

Two more that cost time to learn:

- **Keep single- and multi-point results in separate directories.** A result file keys its payload
  `muts` for single-point and `variants` for multi-point. Share a directory and each complex's
  file means whichever arm ran last, while the resume check reports "cached" — so an entire arm
  can vanish in silence. `consolidate()` refuses to mix them.
- **`1KBH` cannot be repaired.** `RepairPDB` has run past 22.5 h on it without finishing. It is
  excluded before any worklist is built, because killing a hung job is a per-round remedy and a
  sweep driver is a loop. `SKEMPI_FOLDX_ALLOW_INTRACTABLE=1` overrides that if a future FoldX
  handles it; the entry records that the observation was made against 5.1.

## Repository layout

| | |
|---|---|
| `skempi_foldx/` | `lookup` (resolving a mutation's name to its record), `run` (campaigns), `store` (results), `skempi` (parsing), `terms` (the 12-term contract), `exclusions`, `config` |
| `skempi_foldx/data/` | the computed energies, one JSON per interface definition, plus the provenance report |
| `experiments/coverage_report.py` | how much of *your own* evaluation set these energies cover — the first thing to run against a new setup; needs a `--split` or `--table` and exits with an error without one |
| `experiments/repair_ablation.py` | the repair-count ablation |
| `experiments/repair_sweep/` | the full-SKEMPI repair-count sweep, both arms |
| `docs/` | [determinism](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/DETERMINISM.md), [the store](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/STORE.md), [how this protocol compares with the literature](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/PROTOCOL.md), [references](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/REFERENCES.md), [attaching the energies to an evaluation set](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/USAGE.md) |
| upgrading | [what differs between 0.1.0 and 0.2.0](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/VERSIONS.md), record by record, and the [changelog](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/CHANGELOG.md) |

## Citation

Cite the version you installed, because a version names a measurement here rather than a code
state: values from 0.1.0 and 0.2.0 are different measurements of the same quantities and must not
be mixed in one table.
[CITATION.cff](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/CITATION.cff) carries the
machine-readable record, and GitHub's *Cite this repository* renders it. There is no DOI for the
software yet, so cite the repository and the tag.

Cite the upstream sources as well — SKEMPI 2.0 for the measurements and FoldX for the energies —
with DOIs in
[docs/REFERENCES.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/REFERENCES.md).
SKEMPI's CC BY 4.0 makes that attribution a licence term rather than a courtesy;
[NOTICE](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/NOTICE) says what it requires.

## Provenance and licence

The code is MIT. The data is derived from SKEMPI 2.0 (CC BY 4.0) and produced with licensed
software — read [NOTICE](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/NOTICE) before
redistributing it. Citations with DOIs for the data sources, the tools and the comparator methods
are in
[docs/REFERENCES.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/REFERENCES.md). Every
work named in the docs is cited there.

The code originated in a research fork of MuLAN, a sequence-based ΔΔG predictor cited in
[docs/REFERENCES.md](https://github.com/cchin29/skempi-foldx/blob/v0.2.2/docs/REFERENCES.md), where
these energies began as a score channel. None of that project's code is present here, which is what
leaves this package free to carry a permissive licence.

## Acknowledgements

This work was carried out during a 2026 summer research internship at the Laboratoire de Biologie
Computationnelle, Quantitative et Synthétique — the Laboratory of Computational, Quantitative and
Synthetic Biology ([CQSB](https://lcqb.fr/), UMR 7238, CNRS–Sorbonne Université), Paris.

The internship was supported by a fellowship from the France-Stanford Center for Interdisciplinary
Studies, Stanford Global Studies Division, Stanford University.
