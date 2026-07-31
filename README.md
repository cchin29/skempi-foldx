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
  against, and ideally the same one other papers quote. Read it here rather than reproduce it.
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
pip install skempi-foldx
```

No dependencies outside the standard library, and the results are installed with the package.

```python
from skempi_foldx import load_bundled_store

store = load_bundled_store()                      # {pdb: {mutation: {term: kcal/mol}}}
store["1BRS"]["DA52A"]["Interaction Energy"]      # -0.6806  (mutant − wild-type)
```

Mutation keys are SKEMPI's own form: `DA52A` is Asp→Ala at position 52 of chain A. Each record
carries all twelve `AnalyseComplex` terms as mutant minus wild-type, plus two fields of
bookkeeping:

```python
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
number rather than twelve, use `Interaction Energy` — the other eleven decompose it.
`skempi_foldx.TERMS` gives the canonical order for flattening records into a feature vector.

`_source` names the campaign that produced the value, and it is worth attention rather than being
ignorable bookkeeping: values from different campaigns were computed against different mutation
lists, which changes them for the reason described under
[running the pipeline](#running-the-pipeline). `cleaned` is SKEMPI's own mutation string for the
record.

## Contents and coverage

| | complexes | entries | |
|---|---:|---:|---|
| `load_bundled_store()` | 322 | 4238 | single-point mutations |
| `load_bundled_store("results_mp")` | 152 | 1765 | multi-point variants (two or more substitutions) |

Measured against what SKEMPI itself defines — 4337 single-point mutations over 323 complexes,
1848 multi-point variants over 153 — that is **97.7%** and **95.5%**. The gap decomposes exactly:

- **`1KBH`** accounts for 3 single-point and 83 multi-point entries. FoldX cannot repair it (see
  below), so no method that depends on a repaired structure will have these.
- The remaining **96 single-point** entries were never computed, because the campaigns that
  contributed those complexes ran against per-dataset mutation *subsets* rather than each
  complex's full list. Excluding `1KBH`, that puts single-point at 97.8% and **multi-point at
  100%** of what is reachable.

Coverage and the missing rows are detailed in [docs/STORE.md](docs/STORE.md).

Every value was computed with a **single** `RepairPDB` pass. That is not recorded per-record, so
[docs/STORE.md](docs/STORE.md) states the evidence for it — it determines what these numbers can
legitimately be compared against.

## Running the pipeline

This requires a FoldX 5 binary (`FOLDX_BIN`), SKEMPI's cleaned PDBs, and `skempi_v2.csv`. FoldX
is licensed software and is not redistributed here — free for academic and non-profit
institutions, paid for commercial use.

**Read [docs/DETERMINISM.md](docs/DETERMINISM.md) first.** The short version, because it is the
failure mode most likely to cost a week:

FoldX `BuildModel` is deterministic — same structure and same mutation list, same answer to the
last decimal, over 1722 measured pairs. But it walks the mutation list *sequentially in one
process*, so a mutation's ΔΔG depends on every entry before it in that list. **Computing a subset
of a complex's mutations gives different numbers than computing all of them** — measured spreads
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
| `skempi_foldx/` | `run` (campaigns), `store` (results), `skempi` (parsing), `terms` (the 12-term contract), `exclusions`, `config` |
| `skempi_foldx/data/` | the computed energies, one JSON per complex, plus the consolidation report |
| `experiments/repair_ablation.py` | does repairing more than once change the answer? |
| `experiments/repair_sweep/` | the full-SKEMPI repair-count sweep, both arms |
| `docs/` | [determinism](docs/DETERMINISM.md), [the store](docs/STORE.md), [how this protocol compares with the literature](docs/PROTOCOL.md), [references](docs/REFERENCES.md) |

## Provenance and licence

The code is MIT. The data is derived from SKEMPI 2.0 (CC BY 4.0) and produced with licensed
software — read [NOTICE](NOTICE) before redistributing it. Citations with DOIs for the data
sources, the tools and the comparator methods are in [docs/REFERENCES.md](docs/REFERENCES.md).
That list is not exhaustive: seven works named in the docs carry no citation because no DOI for
them could be confirmed, and REFERENCES.md names each one rather than leaving the gap silent.

Extracted from a research fork of MuLAN, where this began as a score channel for a ΔΔG predictor.
None of that project's code is here, which is what lets this one be permissively licensed.

## Acknowledgements

This work was carried out during a 2026 summer research internship at the Laboratoire de Biologie
Computationnelle, Quantitative et Synthétique (UMR 7238, CNRS–Sorbonne Université), Institut de
Biologie Paris-Seine, Paris.

The internship was supported by a fellowship from the France-Stanford Center for Interdisciplinary
Studies, Stanford Global Studies Division, Stanford University.
