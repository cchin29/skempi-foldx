# Determinism, and the mutation-list effect

**Read this before running a campaign.** It is the single most important operational fact about
this pipeline, and the intuitive reading of the evidence is the wrong one.

The short version: FoldX 5.1 `BuildModel` is deterministic. Same repaired structure and same
`individual_list.txt` prefix gives the same number, to the last decimal, always. But a mutation's
ΔΔG depends on *every entry preceding it in that list* — so computing a subset of a complex's
mutations does not give you the same answers as computing all of them.

The inference a spread like this invites — "FoldX is stochastic, so seed it or average over
runs" — is the wrong one. There is nothing to average. Seeding would do nothing. The fix is to
control the mutation list, and `worklist_single_point()` does that by construction.

---

## The evidence

The five campaign result directories disagree with each other on two thirds of the mutations they
share, by up to 11 kcal/mol. That looks like stochastic side-chain optimisation and is not: it is
a setup effect, and it is fully explained.

`BuildModel` walks `individual_list.txt` sequentially in **one process**. Each entry's side-chain
optimisation inherits the state left by the entries before it, and a fresh wild-type reference is
re-optimised per entry — so both terms of `IE(mutant) − IE(wildtype)` move together. A mutation's
result is therefore a function of

    (repaired structure, every entry preceding it in individual_list.txt)

not of the mutation alone.

Measured across the five campaigns, on 5873 pairwise comparisons of the same mutation of the same
complex:

| | pairs | differ | median \|Δ\| | max \|Δ\| |
|---|---|---|---|---|
| identical repaired structure **and** identical list prefix | 1722 | **0** | 0 | **0** |
| different list prefix | 4151 | 3621 | 0.067 | 11.03 |

Zero counterexamples. Same inputs, same output, to the last decimal.

The repaired structures were never the cause: 252 complexes appear in more than one campaign and
**every one has a byte-identical `*_Repair.pdb`** — each was repaired exactly once and copied.
The chain groups were identical too. The only variable was which mutations each campaign asked
for, and a campaign covering a SKEMPI *subset* puts a given mutation at a different line number.

Worked example — `1PPF GB32Y`, same structure (`146d64b021…`) in all four campaigns:

| campaign | list length | index | list prefix | Interaction Energy |
|---|---|---|---|---|
| S1131 | 171 | 56 | `1d0364266aef` | 8.5354 |
| S4169 | 190 | 56 | `1d0364266aef` | **8.5354** |
| S2003 | 181 | 54 | `7d2f95aea5ae` | 9.8185 |
| S1102 | 177 | 53 | `82b0cf61579f` | 19.5701 |

Two campaigns with the same prefix agree exactly; the others do not.

**What follows from this.** Re-running any campaign with its original mutation list reproduces it
exactly, so the published store is reproducible. But a *subset* computation is not
interchangeable with a full one. `skempi_foldx/run.py` therefore computes the **union** of a
complex's mutations (`worklist_single_point`), and `worklist_from_table` — which takes a
per-dataset subset — is retained only for reproducing campaigns that were driven from a
per-dataset mutation table, and is flagged as such.

**This is not inert for the shipped store.** Counting the `_source` field of every shipped
single-point record:

| source campaign | records | complexes | mutation list |
|---|---:|---:|---|
| `foldx_s1102_results_S4169` | 2494 (58.8%) | 210 | a per-dataset **subset** |
| `foldx_skempi_full_results` | 1744 (41.2%) | 112 | the **union** for those complexes |

So roughly three fifths of the shipped single-point values were computed under a subset list and
are **not** the union-list answer for that mutation. `skempi_foldx/data/CONSOLIDATION_REPORT.txt`
corroborates
it directly: 1970 mutations appear in more than one source with differing values, and S4169 is
the kept source — including `1PPF GB32Y` at 8.5354, which is exactly the S4169 row of the worked
example above, where the same mutation reads 9.8185 and 19.5701 under other lists.

What that does and does not mean:

* Every shipped value is **reproducible** — re-run its campaign's list and you get it back exactly.
* Values are **internally consistent per complex**: a complex's records all come from one campaign,
  so within a complex the comparison is sound.
* But a value from an S4169 complex is **not interchangeable** with one computed over the union,
  and the 210 complexes above are not directly comparable to the 112 on that axis.

If you need union-list answers throughout, recompute with `worklist_single_point()`, which takes
the union by construction. The precedence order in `consolidate()` prefers the most complete list
available for each complex; it cannot manufacture a union that was never computed.

## How reproducible is a FoldX number here?

Because the spread is a function of the mutation list rather than of chance, "reproducibility"
has two different answers depending on what you are asking.

**Re-running a campaign as it was run: exact.** Same structure, same list, same numbers.

**Computing the same mutation under a different subset: not exact.** Quantified over the four
independently-listed campaigns (1528 mutations shared, 10631 campaign-pair comparisons), on
`Interaction Energy`:

| | |
|---|---|
| pooled RMSD between campaigns | **0.45 kcal/mol** |
| median absolute difference | 0.12 |
| 95% limits of agreement | ±0.89 |
| implied single-choice SD | 0.32, against a population SD of 1.63 |
| intraclass correlation | **0.96** (Pearson *r* between any two campaigns 0.956–0.987) |
| 99th percentile \|Δ\| | 2.6 |
| maximum | 11.0 (6.8 σ) |

Per-term ICC ranges from 0.99 (Van der Waals, Solvation Hydrophobic) down to **0.90 (torsional
clash)**, 0.93 (Van der Waals clashes) and 0.92 (entropy mainchain) — the clash and entropy terms
are the least stable, and they are terms the scalar arm never sees.

**It is concentrated, not diffuse.** 185 of 252 shared complexes agree on *every* mutation; 10
disagree on every one. Three complexes — `1PPF`, `1R0R`, `3SGB`, the OMTKY3 protease–inhibitor
saturation-mutagenesis sets — account for **half** of all absolute disagreement, and the top ten
for 80%. That bimodality is exactly what the list-position mechanism predicts: two campaigns
whose lists match for a complex agree on all of it, and two whose lists diverge early disagree on
most of it. Note these are not rare rows — the three account for ~14% of split rows.

By interface position, mean \|Δ\| is 0.154 kcal/mol at interface-core positions against
0.003–0.007 at non-interface ones (Kruskal–Wallis p < 1e-97) — expected, since a non-interface
mutation's *binding* ΔΔG is ≈0 under any list ordering.

**Sign stability.** Two campaigns disagree on the sign of `Interaction Energy` for 2.4% of shared
mutations — but almost all are noise around zero (median magnitude 0.22 kcal/mol). Flips where
both sides exceed 0.5 kcal/mol: **4 of 3096**. At ±2.0 kcal/mol: **none**.

**Downstream exposure.** Substituting one campaign's values for another's and re-standardizing
per fold moves, of FoldX-covered rows: ~1.3% by more than 0.5 standardized units on the **scalar**
arm, and ~5% on the **12-term** arm — the decomposed arm is roughly twice as exposed, because it
reads twelve chances at the noise including the least stable terms. Under an adversarial
worst-case swap those become 2.6% and 9.3%.

None of this affects the shipped store, which is a fixed set of byte-copies (below). The exposure
is to *regeneration*.

## Why it is deterministic, and what the literature says

Three independent lines agree, which matters because the vendor never states determinism either
way — no FoldX document claims it, and none denies it.

**1. The parameter semantics.** The entire suite exposes exactly one RNG control:
`--timeSeedRotabase` (bool, **default `false`**), documented as "set the seed of the rotabase to
the time (useful for MC sampling)". There is no `--seed` and no `--rngSeed`. An opt-in flag to
seed from the wall clock, off by default, only makes sense if the default seed is a fixed
constant. **Leave it false.** Setting it true would make runs genuinely irreproducible.

**2. The vendor's description of `numberOfRuns`.** From the Switch Lab's own manual: with more
than one run "the algorithm will do the same mutations but changing the rotamer set used and the
order of the rotamer moves". So variation between runs *within* one invocation is deliberate and
indexed — it does not imply variation between invocations. Every campaign here used
`--numberOfRuns=1`, i.e. run 1 every time.

**3. The direct measurement** in the table above: 1722 pairs with identical structure and
identical list prefix, zero differing. That is the evidence the literature does not contain, and
it is cheap to re-run if you want to confirm it on another binary.

**The published magnitudes independently rule out the stochastic explanation.** FoldX's stated
error margin is ~0.5 kcal/mol; its calibration SD against experiment is 0.46; its *structural*
sensitivity — the spread across different PDB structures of the same protein — is 0.61
(Caldararu, Blundell & Kepp, *BMC Bioinformatics* 2021). The widest published figure of any kind
is a ±3.5 kcal/mol prediction *interval* for binding ΔΔG, which is an accuracy bound rather than
run-to-run scatter. **Nothing in the literature is within an order of magnitude of the 11
kcal/mol observed here**, which is itself an argument that the spread is not FoldX noise.

## Is a single `RepairPDB` defensible? Measured, and the answer is *no*

FoldX's manual recommends one repair, CATH-ddG uses one, and the only controlled test in the
literature (Usmanova et al. 2018) found iterating made no difference. That predicts the repair
count should not matter. **It does.**

Measured on the 13 CATH T≥10 complexes, 5× versus 1× repair with identical mutation lists.
**Single-point only** — the CATH test set is 39% multi-point, so this is the CATH-*single*
number; the full-SKEMPI run covers both arms:

| | |
|---|---|
| structure movement, round 1 → 2 | median **0.25 Å** (converged by round 4) |
| mutations whose `Interaction Energy` changes | **97.6%** of 334 |
| median \|Δ\| / p90 / max | 0.070 / 0.62 / **8.51** kcal/mol |
| FoldX-alone per-structure Spearman, 1× | 0.3984 [0.2759, 0.5193] |
| FoldX-alone per-structure Spearman, **5×** | **0.4585** [0.3455, 0.5695] |
| paired Δ (cluster bootstrap over complexes) | **+0.060, 95% CI [+0.020, +0.100]** |

The published FoldX baseline on this split is 0.4303 — between the two. That is the likely
explanation for a single-repair FoldX-alone number sitting below the published one.

This does not contradict Usmanova: they measured *folding* ΔΔG and self-consistency bias; this
measures *binding* ΔΔG correlation with experiment.

Caveats worth keeping attached to the number: 11 complexes on one tier, and this is FoldX used
*alone* — whether a downstream model consuming these energies improves on better inputs requires
a retrain. Reproduce or widen with `experiments/repair_ablation.py`.

## `--numberOfRuns=1` is likewise defensible

It is the vendor default and the vendor's own documented recommendation ("normally it should be
set to 1"). The FoldX authors' 2025 revision paper uses 5 with the median, but three things
qualify that:

1. the reported gain is small (R 0.705 → 0.711, RMSE 1.250 → 1.238) **and was measured on folding
   benchmarks, not on SKEMPI binding data**;
2. it addresses *within-invocation* rotamer exploration, which is **not** the effect measured
   above — it would not remove the mutation-list dependence; computing the union list once does;
3. every published baseline row these results sit beside was produced with undocumented
   single-run defaults.

What actually distinguishes a defensible pipeline here is *reporting*, not the setting. The best
template found in this literature is a single sentence from Meli et al. (*IJMS* 2024): "The macros
'PositionScan' and 'BuildModel' were run by using default settings (i.e., number of runs: 1, pH 7,
temperature 298 K, and ionic strength 0.05 M)." That is more protocol detail than any of the
comparator papers give.

**Platform note.** The multi-point campaign ran on Linux while the single-point campaigns ran on
Apple Silicon. There is no vendor statement on cross-platform bitwise agreement, so the two should
not be assumed interchangeable. In practice it cannot matter here: the multi-point store shares
zero mutation keys with the single-point ones.

