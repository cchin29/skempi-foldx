# FoldX repair-count sweep — full SKEMPI

Repairs every complex **4 times**, saving each round, then runs BuildModel + AnalyseComplex from
**each** round. That gives the ΔΔG convergence curve, not just the endpoints. The repair chain is
paid once per complex and shared across rounds, and across the SP and MP arms.

Why: on the 13 CATH T≥10 complexes, 5× repair moved the FoldX-alone per-structure Spearman from
0.398 to 0.459 (paired Δ +0.060, 95% CI [+0.020, +0.100]) — past the published single-point FoldX row of 0.4458.
Structure converges by round 4 (median RMSD 0.013 Å from round 4→5). This widens that measurement
to all of SKEMPI.

## Components

| | |
|---|---|
| `../repair_ablation.py` | SP arm: `probe`, `run`, `sweep`, `compare` |
| `sweep_mp.py` | MP arm — the half `sweep` omits. Reuses the SP repair chains verbatim. |
| `residue_check.py` | Does any repair round change *which* mutations pass WT-residue validation? |
| `list_hashes.py` | Hashes the `individual_list.txt` files BuildModel actually consumed |
| `launch_sweep.sh` | SP launcher |
| `chain_finalize.sh` | verify → MP build → verify → record provenance, in that order |

## Invocation

```bash
export FOLDX_BIN=/path/to/foldx
JOBS=36 ./launch_sweep.sh                  # SP arm
JOBS=36 SWEEP_EXCLUDE="" ./chain_finalize.sh   # everything after it
```

Paths default to `<repo>/scratch/...` and are overridable via `SKEMPI_FOLDX_ROOT`, `SWEEP_RUN`,
`SWEEP_STORE`, `SKEMPI_PDBS`, `SKEMPI_CSV`. Resumable at complex granularity.

## Cost model

The two phases are billed differently, and a single "per mutation" rate predicts neither.
The last point is about provenance rather than cost, but it is the one that is expensive to
learn the hard way.

**1. Repair cost tracks atom count, not mutation count.** The plausible reasoning is that
BuildModel dominates, so the wall-clock floor belongs to the complex with the most mutations —
`1CHO` at 275 mutations and only 2150 atoms, or `3BT1` at 240. It does not. The floor is
**repair**, on `1KBH`: three single-point mutations, 33 060 atoms, RepairPDB running 10 h 24 m to
reach residue A36 of 2120. Budget the repair phase by atom count.

`1KBH` is in `skempi_foldx/exclusions.py` and is dropped before any worklist is built.

**2. BuildModel cost tracks entry composition, not entry count.** A multi-point entry repacks two
or more positions, so it costs roughly **3× a single-point entry**. The two arms of a full-SKEMPI
sweep therefore invert against intuition:

| arm | entries | complexes | measured per round |
|---|---:|---:|---|
| single-point | 4334 | 322 | faster, despite 2.5× the entries |
| multi-point | 1765 | 152 | **~2 h** |

Measured shape of one multi-point round, 152 complexes on 36 workers:

| window | complexes finished |
|---|---:|
| first 40 min | 97 |
| next 60 min | 45 |
| remainder | 10 (long tail) |

Two thirds land in the first third of the wall clock; the rest is tail. Four rounds run
sequentially inside one `sweep_mp.py build --rounds 1 2 3 4`, so budget ~8 h for the multi-point
arm and expect the estimate to be dominated by the tail rather than the mean.

**3. Each round's finish is set by one complex.** The tail is a single high-mutation-count
complex that the whole round waits on: `3BT1` (240 mutations) for single-point, `1JTG` (136) for
multi-point. Adding workers past that point buys nothing — the floor is one complex's serial
BuildModel walk, because entries within a complex are evaluated sequentially in one process.
That serialisation is the same property that makes a mutation's ΔΔG depend on the entries before
it; here it shows up as a scheduling floor.

**4. Do not delete `work*/` as "regenerable".** The `chain_work/` and per-round `work/` trees are
gigabytes, and it is tempting to bring back only `round_*/results/`. The energies are regenerable;
**the record of what produced them is not.** `individual_list.txt` lives in those trees, and a
mutation's ΔΔG depends on every entry preceding it in that list — so the list is what makes a
number interpretable at all.

`chain_finalize.sh` runs `list_hashes.py` **before** any cleanup, which hashes those bytes into
`list_manifest.json` and injects `meta.mutation_list_sha256` into each result JSON. After that the
work trees are genuinely disposable. Before it, they are the only copy.

## `residue_check.json`

`clean` and `full_coverage` are separate fields and both matter. `clean: true` with
`full_coverage: false` means *clean over the complexes checked so far* — a mid-run snapshot reads
`chains_complete: 78, chains_missing: 267, clean: true`. Gate on both: **a check that matches
nothing reports clean, which is indistinguishable from actually clean**, and that is how a whole
stale tier goes unnoticed.

## Output

```
scratch/foldx_repair_ablation/sweep_4x/chain_work/<PDB>/repair_round_<r>.pdb
scratch/foldx_repair_ablation/sweep_4x/round_<r>/results/<PDB>.json      # SP
scratch/foldx_repair_ablation/sweep_4x/round_<r>/results_mp/<PDB>.json   # MP
scratch/foldx_repair_ablation/sweep_4x/{list_manifest,residue_check}.json
```

SP and MP results live in **separate** directories on purpose: `process_complex` caches at
`results_dir/<pdb>.json` and keys the payload `muts` for `MODE_AUTHOR` but `variants` for
`MODE_VARIANT`. Sharing one directory makes each complex's JSON mean whichever arm ran last, and
the resume check then skips the other arm in silence. An SP and an MP energy for one complex are
**not** two entries of one list.

## Comparing rounds

```bash
python3 ../repair_ablation.py compare \
    --store    scratch/foldx_repair_ablation/sweep_4x/round_1/results \
    --ablation scratch/foldx_repair_ablation/sweep_4x/round_4/results
```

Repeat for rounds 2 and 3 to see where the energy converges — it may settle earlier than the
structure does, which would be the useful result.
