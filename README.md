# skempi-foldx

**FoldX binding ΔΔG for SKEMPI 2.0 — the pipeline, and the computed results.**

```python
from skempi_foldx import load_bundled_store

store = load_bundled_store()          # the shipped single-point results
store["1BRS"]["DA52A"]["Interaction Energy"]     # -0.6806 kcal/mol, mutant − wild-type
```

Two things ship here, and for most people the second is the one that matters.

**The results.** 322 complexes / 4238 single-point mutations, and 152 complexes / 1765
multi-point variants, decomposed into FoldX's 12 energy terms — 99.2% of the rows SKEMPI
defines. Reproducing them takes a FoldX licence (free for academic and non-profit institutions,
paid for commercial use) and CPU-weeks. Reading them takes
neither, and this package has **no dependencies outside the standard library**.

**The pipeline** that produced them: `RepairPDB → BuildModel → AnalyseComplex` over SKEMPI's
chain groups, yielding ΔΔG_int = IE(mutant) − IE(wild-type). Resumable per complex, parallel
across complexes, and explicit about what it could not compute rather than silently filling in.

**This is not a FoldX wrapper.** It does not wrap, bundle, or replace FoldX — for Python
bindings see [pyfoldx](https://github.com/leandroradusky/pyfoldx). What is here is a dataset-scale
campaign runner plus the campaign's results.

## Install

```bash
pip install skempi-foldx        # to read the results
```

Running the pipeline additionally needs a FoldX 5 binary (`FOLDX_BIN`), SKEMPI's cleaned PDBs,
and `skempi_v2.csv`. FoldX is licensed software and is not redistributed here — free for
academic and non-profit use, paid for commercial. See [NOTICE](NOTICE).

## Read this before running a campaign

FoldX 5.1 `BuildModel` **is deterministic**. Measured over 1722 pairs computed from an identical
repaired structure *and* an identical `individual_list.txt` prefix: zero differ, to the last
decimal.

But a mutation's ΔΔG is a function of `(repaired structure, every entry preceding it in
individual_list.txt)` — not of the mutation alone. BuildModel walks the list sequentially in one
process, and each entry's side-chain optimisation inherits the state the previous entries left.
It also re-optimises a fresh wild-type reference per entry, so *both* terms of the subtraction
move.

The consequence is severe and easy to hit: **computing a subset of a complex's mutations gives
different numbers than computing the full set.** Across five campaigns that each covered a
different SKEMPI subset, this produced spreads of median 0.067 and up to 11.03 kcal/mol on the
same mutation of the same structure — with no randomness involved at all.

So: always compute the union of a complex's mutations. `worklist_single_point()` does this by
construction. `worklist_from_table()` does not, and exists only for reproducing campaigns that
were driven from a per-dataset mutation table.

This is worth stating loudly because the natural reading of a spread like that is "FoldX is
stochastic, seed it" — which is wrong, and leads you to average noise that isn't there instead
of fixing the mutation list. See [docs/DETERMINISM.md](docs/DETERMINISM.md).

## Layout

| | |
|---|---|
| `skempi_foldx/` | the package: `run` (campaigns), `store` (results), `skempi` (parsing), `terms` (the 12-term contract), `exclusions` (structures FoldX cannot handle), `config` |
| `skempi_foldx/data/results_sp/`, `.../results_mp/` | the computed energies, one JSON per complex — inside the package, so `pip install` ships them |
| `skempi_foldx/data/CONSOLIDATION_REPORT.txt` | how the five source campaigns were unioned, and where they disagree |
| `experiments/repair_ablation.py` | does repairing more than once change the answer? |
| `experiments/repair_sweep/` | the full-SKEMPI repair-count sweep, both arms |

## Two operational facts that cost time to learn

**Single- and multi-point results must live in separate directories.** A result JSON keys its
payload `muts` for single-point and `variants` for multi-point. Share one directory and each
complex's file means whichever arm ran last, while the resume check happily reports "cached" —
so an entire arm can be skipped in silence. `consolidate()` refuses to mix them.

**`1KBH` cannot be repaired.** RepairPDB has run past 22.5 h on it without finishing, on two
machines. It is in `skempi_foldx/exclusions.py` and is dropped before any worklist is built —
because killing the hung job is a per-round remedy and a sweep driver is a loop. Set
`SKEMPI_FOLDX_ALLOW_INTRACTABLE=1` to try it anyway on a future FoldX; the entry records that
the observation was made against 5.1.

## Provenance and attribution

The code is MIT. The data is derived from SKEMPI 2.0 and produced by licensed software —
read [NOTICE](NOTICE) before redistributing it.

Extracted from a research fork of MuLAN, where this began as a score channel for a ΔΔG
predictor. None of that project's code is here, which is what lets this one be permissively
licensed.
