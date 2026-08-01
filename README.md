# skempi-foldx

**Precomputed FoldX binding energies for SKEMPI 2.0 — the standard physics baseline for
mutation-effect prediction, without the licence or the compute.**

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

- **Benchmarking a ΔΔG predictor.** A FoldX column is the conventional baseline to report
  against. This one is fully specified — protocol, repair count, mutation lists and known defects
  are all written down, which is more than the surveyed comparator papers state
  ([docs/PROTOCOL.md](docs/PROTOCOL.md)). It is not the same computation as any published FoldX
  row and should not be quoted as one.
- **Using FoldX as a feature.** The 12 decomposed energy terms per mutation are a ready-made
  structural input for a learned model — computed once from structure, and independent of the
  encoder and the evaluation split.
- **Running FoldX at dataset scale.** The pipeline handles SKEMPI parsing, chain-group mapping,
  resumption, and parallelism — and documents the failure modes that are expensive to discover
  independently.

**This is not a FoldX wrapper.** It does not wrap, bundle, or replace FoldX — for Python bindings
see [pyfoldx](https://github.com/leandroradusky/pyfoldx). This is a dataset-scale campaign runner
and the campaign's results.

## Quickstart

```bash
pip install git+https://github.com/cchin29/skempi-foldx
```

Or, from a clone of this repository:

```bash
pip install .
```

This package is not on PyPI, so `pip install skempi-foldx` resolves to nothing. There are no
dependencies outside the standard library, and the energies ship inside the wheel — no separate
data download follows. `pip install '.[dev]' && pytest` runs the suite, which pins the store's
shape and the quickstart values against the shipped data. Figures derived from SKEMPI's table or
from evaluation splits cannot be pinned here, because neither ships with this package.

```python
from skempi_foldx import FoldxLookup

fx = FoldxLookup()                                # both arms, 6003 records
fx.get("1BRS", "DA52A")["Interaction Energy"]     # -0.6806  (mutant − wild-type)
fx.vector("1BRS", "DA52A")                        # the 12 terms in TERMS order

rows = [("1BRS", "DA52A"), ("1ACB", "LI38D")]     # any (pdb, mutation) pairs
fx.coverage(rows)                                 # (covered, total, missing)
fx.energies_for(rows)                             # {(pdb, mutation): [12 floats]}
```

## Naming, terms and provenance

`FoldxLookup` is the entry point rather than the raw dictionary because **a mutation carries more
than one name**. The same substitution is `LI38D` in SKEMPI's author-chain form and `LB38D` in the
role-chain form split files use, and the store holds both — 3012 single-point records under the
first convention and 1226 under the second. A plain dictionary join keyed on the author form
misses about 29% of the store without raising; keyed on the role form it misses 71%. Measured on a
4159-row split whose labels are role-chain: 65.4% resolved by dictionary lookup, **99.6%** by
`FoldxLookup(skempi_csv=..., mapping_dir=...)`. `fx.require()` raises instead of returning `None`,
and its error names the fix.

The raw store stays available as `load_bundled_store()` for auditing and for consumers doing their
own chain mapping — see [docs/STORE.md](docs/STORE.md), and
[docs/USAGE.md](docs/USAGE.md) for attaching these energies to an evaluation set.

`DA52A` is Asp→Ala at position 52 of chain A. Each record
carries all twelve `AnalyseComplex` terms as mutant minus wild-type, plus bookkeeping — `_source`
and, for single-point records, `cleaned`; a multi-point key already is SKEMPI's variant string:

```python
from skempi_foldx import load_bundled_store

store = load_bundled_store()
store["1BRS"]["DA52A"]
# {'Interaction Energy': -0.6806,   'Backbone Hbond': 0.0,
#  'Sidechain Hbond': 0.0,          'Van der Waals': 0.0283,
#  'Electrostatics': -0.2432,       'Solvation Polar': -0.1739,
#  'Solvation Hydrophobic': 0.0086, 'Van der Waals clashes': -0.0001,
#  'entropy sidechain': -0.0381,    'entropy mainchain': 0.0,
#  'torsional clash': -0.0,         'backbone clash': 0.0,
#  'cleaned': 'DA52A',              '_source': 'S4169'}
```

A negative interaction energy means FoldX predicts the mutant binds *more* tightly. For a single
number rather than twelve, use `Interaction Energy`. **The other eleven are components of it, not
a decomposition of it**: they do not sum to it, so a model given the components must keep
`Interaction Energy` as well, and a consumer checking the store by summation will wrongly conclude
it is corrupt. [docs/USAGE.md](docs/USAGE.md) gives the size of the residual.
`skempi_foldx.TERMS` gives the canonical order for flattening records into a feature vector.

`_source` names the campaign that produced the value, and it changes how the value should be read:
values from different campaigns were computed against different mutation lists, which changes them
for the reason described under [running the pipeline](#running-the-pipeline). **2494 of the 4238
single-point values (58.8%) are labelled `S4169`, meaning they were produced by a campaign built
around that benchmark rather than by one built around each complex's full mutation list** — which
matters because a value is a property of the list it was computed in. For 27 of those complexes,
carrying 287 records, that list is a strict subset of the complex's full one; for the other 183,
carrying 2207, it holds the same mutations. How large the list effect is, measured between
campaign pairs, is in [docs/DETERMINISM.md](docs/DETERMINISM.md). `cleaned` is SKEMPI's own
mutation string for the record.

## Contents and coverage

| | complexes | entries | |
|---|---:|---:|---|
| `load_bundled_store()` | 322 | 4238 | single-point mutations |
| `load_bundled_store("results_mp")` | 152 | 1765 | multi-point variants (two or more substitutions) |

Measured against what SKEMPI itself defines — 4337 single-point mutations over 323 complexes,
1848 multi-point variants over 153 — that is **97.7%** and **95.5%**. The gap decomposes exactly:

- **`1KBH`** accounts for 3 single-point and 83 multi-point entries. FoldX cannot repair it (see
  below), so no method that depends on a repaired structure will have these.
- **80 single-point** entries were never computed, because the campaigns that contributed those
  complexes ran against per-dataset mutation *subsets* rather than each complex's full list.
- The remaining **16** are eight complete pairs, both members absent: each carries the same
  substitution on two chains of one group, which the role letter cannot tell apart, so the chain
  mapping refuses both rather than guess. See [docs/STORE.md](docs/STORE.md).

Excluding `1KBH`, whose `RepairPDB` does not terminate under FoldX 5.1, that puts single-point at
**97.8%** and multi-point
at **100%** of what is reachable, counting by PDB code.

Those denominators count SKEMPI rows keyed by PDB code, which is how this package's parser counts
them, and three codes carry **two interface definitions each** — `2C5D`, `3SE3`, `3SE4`. The
parser keys on the code alone, so one definition wins and every mutation is pooled under it:
**31 shipped records are scored against an interface that excludes the chain they mutate.** Twenty
are all-zero as an artifact of the pairing rather than as a measurement. Ten of the remaining
eleven are multi-point variants carrying only one component's contribution, |IE| from 0.001 to
4.355 kcal/mol, which does not look wrong; the eleventh is a single-point record at 0.0002. Keying
per interface definition instead moves both denominators and both coverage figures downward, so
the 100% above is a property of the by-code parse rather than of the store —
[docs/STORE.md](docs/STORE.md) gives that parse for both arms, and 0.2.0 adopts it
([CHANGELOG.md](CHANGELOG.md)).

`FoldxLookup.is_interface_suspect(pdb, mutation)` identifies the affected records without needing
SKEMPI, and accepts any convention the lookup resolves; the module-level
`skempi_foldx.interface_suspect(pdb, cleaned)` behind it takes SKEMPI's author form only.
`FoldxLookup` also warns on construction. The mutation lists and the recomputation
plan are in [docs/USAGE.md](docs/USAGE.md).

Against SKEMPI's measured affinities the shipped values reach Spearman **0.435** on the
single-point arm and **0.574** on the multi-point one — ordinary for FoldX, and the check coverage
cannot perform. Coverage, agreement and the missing rows are detailed in
[docs/STORE.md](docs/STORE.md).

## Scope and stability

Version 0.1.0. The energies are the artifact and they are tied to this version. A mutation's ΔΔG
depends on the whole mutation list it was computed in, so a rebuilt store is a different
measurement rather than a correction of this one — a published result quoting these numbers should
pin `skempi-foldx==0.1.*` and say so.

Two changes are planned for 0.2.0, and both change published energy values, including for
mutations whose current value is not wrong: the store will be rebuilt against union mutation lists
throughout, and each SKEMPI interface definition will be computed as its own pairing. The second
also moves the coverage denominators above. Values from 0.1.0 and 0.2.0 are not interchangeable
and must not be mixed within one table.

`TERMS` and its order are stable, because a consumer's feature layout depends on them. Nothing
else in the API is settled yet.

Every value was computed with a **single** `RepairPDB` pass. That is not recorded per-record, so
[docs/STORE.md](docs/STORE.md) states the evidence for it — a 717-directory on-disk check for the
single-point arm, the builder's code for the multi-point one — and it determines what these
numbers can legitimately be compared against.

## Running the pipeline

This requires a FoldX 5 binary (`FOLDX_BIN`), SKEMPI's cleaned PDBs, and `skempi_v2.csv`. FoldX
is licensed software and is not redistributed here — free for academic and non-profit
institutions, paid for commercial use.

**Read [docs/DETERMINISM.md](docs/DETERMINISM.md) first.** The short version, because it is the
failure mode most likely to cost a week:

FoldX `BuildModel` is deterministic — same structure and same mutation list, same answer to the
last decimal, over 1722 measured pairs. But a mutation's ΔΔG depends on the whole list it was
computed in, and not on the mutation alone. **Computing a subset of a complex's mutations gives
different numbers than computing all of them** — measured spreads
of median 0.067 and up to 11.03 kcal/mol for the same mutation of the same structure.

The intuitive reading of a spread like that is "FoldX is stochastic, seed it". That is wrong, and
it leads to averaging noise that does not exist. The variable to control is the mutation list:
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
| `skempi_foldx/data/` | the computed energies, one JSON per complex, plus the consolidation report |
| `experiments/repair_ablation.py` | the repair-count ablation |
| `experiments/repair_sweep/` | the full-SKEMPI repair-count sweep, both arms |
| `docs/` | [determinism](docs/DETERMINISM.md), [the store](docs/STORE.md), [how this protocol compares with the literature](docs/PROTOCOL.md), [references](docs/REFERENCES.md), [attaching the energies to an evaluation set](docs/USAGE.md) |

## Provenance and licence

The code is MIT. The data is derived from SKEMPI 2.0 (CC BY 4.0) and produced with licensed
software — read [NOTICE](NOTICE) before redistributing it. Citations with DOIs for the data
sources, the tools and the comparator methods are in [docs/REFERENCES.md](docs/REFERENCES.md).
Every work named in the docs is cited there.

The code originated in a research fork of MuLAN, where these energies began as a score channel for
a ΔΔG predictor. None of that project's code is present here, which is what leaves this package
free to carry a permissive licence.

## Acknowledgements

This work was carried out during a 2026 summer research internship at the Laboratoire de Biologie
Computationnelle, Quantitative et Synthétique — the Laboratory of Computational, Quantitative and
Synthetic Biology ([CQSB](https://lcqb.fr/), UMR 7238, CNRS–Sorbonne Université), Paris.

The internship was supported by a fellowship from the France-Stanford Center for Interdisciplinary
Studies, Stanford Global Studies Division, Stanford University.
