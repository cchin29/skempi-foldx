# Determinism, and the mutation-list effect

**Read this before running a campaign.** It is the operational fact this pipeline turns on, and
the intuitive reading of the evidence does not hold.

**On the evidence below.** The determinism table, the agreement statistics and the convergence
tables are computed from campaign output — the `scratch/` directories a run leaves behind — which
is not distributed with this package. The figures here are therefore recorded rather than
reproducible from a checkout. Given that output, `experiments/agreement_stats.py` regenerates the
agreement statistics and `experiments/repair_convergence.py` the convergence tables; the shipped
store is what a reader can check independently.

The short version: FoldX 5.1 `BuildModel` is deterministic. Same repaired structure and same
`individual_list.txt` prefix gives the same number, to the last decimal, always. But a mutation's
ΔΔG depends on *every entry preceding it in that list* — so computing a subset of a complex's
mutations does not give the same answers as computing all of them. The size of that effect,
measured here: over the pairs that differ, the same mutation of the same complex moves by a
median of 0.067 kcal/mol between two campaigns' mutation lists, and by as much as 11.03 — 0.008
over all pairs, most of which are identical.

The inference a spread of that size invites — "FoldX is stochastic, so seed it or average over
runs" — is the wrong one. There is nothing to average. Seeding does not address it: the variable
to control is the mutation list, and `worklist_single_point()` does that by construction.

---

## The evidence

The five single-point campaign result directories disagree on 3621 of the 5870 campaign-pair
comparisons they share, 61.7%, by up to 11 kcal/mol. That looks like stochastic
side-chain optimisation and is not: it is a setup effect, and the dependence is fully
characterised below.

`BuildModel` walks `individual_list.txt` sequentially in **one process**, and a fresh wild-type
reference is rebuilt per entry — which is required, since each mutation repacks different
neighbours and so needs its own reference. Empirically, a mutation's result is a function of

    (repaired structure, the list it was computed in)

not of the mutation alone.

**The mechanism behind that is not established here.** FoldX documents each line of
`individual_list.txt` as an independent mutant built from the input structure, so the coupling is
not the mutant structures accumulating. The likeliest remaining candidate is process-level state
carried between entries — a rotamer search order or random stream not reset per line, which is
consistent with what the vendor does say about `--numberOfRuns` changing the rotamer set and the
order of rotamer moves. That is a hypothesis; what is measured below is the dependence itself.

Nor do the measurements separate *which* preceding entries matter from *how many*: list content
and list position covary in every campaign pair compared here, so "depends on the preceding
entries" and "depends on the index" are not distinguished. The operational consequence is the same
either way — a subset computation is not interchangeable with a union one — which is why the
guidance does not rest on the mechanism.

Measured over all five single-point campaign result directories — every pair of campaigns that
computed the same mutation of the same complex, 5870 such pairs, which is the population the two
rows below partition. `1KBH` is excluded here as it is everywhere else; counting its stale
three-record JSON would give 5873 and 1722:

| | pairs | differ | median \|Δ\| over the differing | max \|Δ\| |
|---|---|---|---|---|
| identical repaired structure **and** identical list prefix | 1719 | **0** | — | **0** |
| different list prefix | 4151 | 3621 | 0.067 | 11.03 |

Zero counterexamples. Same inputs, same output, to the last decimal.

The repaired structures were never the cause: **182 complexes appear in more than one campaign**
(181 excluding the excluded `1KBH`) and every one was repaired exactly once and copied, so the
`*_Repair.pdb` is shared rather than recomputed. Measured directly over the 187 complexes whose
`*_Repair.pdb` exists in more than one campaign work directory — a superset of those 182, since a
complex can be repaired for a campaign without sharing a mutation with another — the file is
byte-identical in every case, 187 of 187 with none differing.
The chain groups were identical too. The only variable was which mutations each campaign asked
for, and a campaign covering a SKEMPI *subset* puts a given mutation at a different line number.

Worked example — `1PPF GB32Y`, same structure (`146d64b021…`) in all four campaigns. The name is
the campaigns' role-chain form; the shipped record is `results_sp/1PPF_E_I.json` under the author
key `GI32Y`, with `GB32Y` in its `role` field, so a raw dictionary join needs the author form:

| campaign | list length | index | list prefix | Interaction Energy |
|---|---|---|---|---|
| S1131 | 171 | 56 | `1d0364266aef` | 8.5354 |
| S4169 | 190 | 56 | `1d0364266aef` | **8.5354** |
| S2003 | 181 | 54 | `7d2f95aea5ae` | 9.8185 |
| S1102 | 177 | 53 | `82b0cf61579f` | 19.5701 |

Two campaigns with the same prefix agree exactly; the others do not.

**Consequences.** Re-running any campaign with its original mutation list reproduces it
exactly, so the published store is reproducible. But a *subset* computation is not
interchangeable with a full one. `skempi_foldx/run.py` therefore computes the **union** of a
complex's mutations (`worklist_single_point`), and `worklist_from_table` — which takes a
per-dataset subset — is retained only for reproducing campaigns that were driven from a
per-dataset mutation table, and is flagged as such.

**This was not inert for the 0.1.0 store, and closing it is what 0.2.0 did.** Counting the
`_source` field of every 0.1.0 single-point record:

| source campaign | records | complexes | mutation list |
|---|---:|---:|---|
| `S4169` | 2494 (58.8%) | 210 | that benchmark's per-complex list |
| `full_skempi` | 1744 (41.2%) | 112 | the **union** for those complexes |

So roughly three fifths of 0.1.0's single-point values came from a campaign built around one
benchmark rather than around each complex's full list. For most of them that made no difference to
the *membership* of the list — 183 complexes carrying 2207 records held the complete list, and 27
carrying 287 held a strict subset — but membership is not the whole of it. The value depends on
the list a mutation was computed in, which is a property of the campaign rather than of the
mutation, so a value from an `S4169` complex was not interchangeable with one computed over the
union, and those 210 complexes were not directly comparable to the other 112.

**The shipped store no longer has that shape.** Every value comes from a union mutation list:
4286 single-point and 1702 multi-point records from the repair sweep's `round_1`, and 119 from the
nine per-definition campaigns. There is one source per record and no precedence order, so there is
no pair of populations here to keep apart. `skempi_foldx/data/CONSOLIDATION_REPORT.txt` records
provenance and nothing to arbitrate.

What survives the rebuild is the *mechanism* — everything above this section. A subset computation
is still not interchangeable with a union one, which is why `worklist_single_point()` takes the
union by construction and `worklist_from_table()` is retained only for reproducing campaigns that
were driven from a per-dataset table.

## Reproducibility

Because the spread is a function of the mutation list rather than of chance, "reproducibility"
has two different answers depending on what is being asked.

**Re-running a campaign as it was run: exact.** Same structure, same list, same numbers.

**Computing the same mutation under a different subset: not exact.** Quantified on
`Interaction Energy` over the campaigns driven from independently constructed mutation lists:
**1527 mutations shared, 5870 campaign-pair comparisons**, of which 3621 differ (61.7%).
Recomputed from the campaign output by
[`experiments/agreement_stats.py`](../experiments/agreement_stats.py); this is the same
population as the determinism table above, not a second one.

| | |
|---|---|
| pooled RMSD between campaigns | **0.524 kcal/mol** |
| median absolute difference | 0.008 over all pairs; 0.067 over the pairs that differ |
| 95% limits of agreement | ±1.03 |
| implied single-choice SD | 0.37, against a population SD of 2.12 |
| intraclass correlation | **0.972** |
| 99th percentile \|Δ\| | 2.09 |
| maximum | 11.03 |

Per-term ICC ranges from 0.991 (Solvation Hydrophobic) and 0.984 (Van der Waals) down to
**0.923 (entropy mainchain)**, 0.935 (Sidechain Hbond) and 0.941 (torsional clash). The entropy,
H-bond and clash terms are the least stable, and a consumer reading only the summary
`Interaction Energy` — itself at 0.972 — never sees them. The revision of this document shipped
at `v0.1.0` named torsional clash as the floor at 0.90; measured, it is third from the bottom at
0.941, and the floor is entropy mainchain.

**Populations.** There is one population, and every statistic in this section is computed over it:
**5870 comparisons across 1527 shared mutations** (5873 over 1530 if the stale three-record
`1KBH` JSON is counted; it is excluded everywhere else, so it is excluded here). The determinism
table above and the agreement statistics below describe the same set, which the differing-pair
count of 3621 confirms independently — the two agree on it exactly.

Three denominators quoted in the `v0.1.0` revision — 10631 comparisons, 1528 shared mutations,
and 3096 for sign stability — are not reproducible from the campaign directories and are retired.
Sign stability is reported below against the same 5870. Every quantity in the table above is
recomputed by `experiments/agreement_stats.py`, so a figure quoted from that earlier revision
should be re-read from this one.

**It is concentrated, not diffuse.** 115 of the 181 complexes holding a shared mutation agree on
*every* one of them; 12 disagree on every one. Three complexes — `1PPF`, `1R0R`, `3SGB`, three
structures of the OMTKY3 inhibitor bound to different proteases, each covered by saturation
mutagenesis (every amino acid tried at every interface position, which is why they hold so many
rows) — account for **68%** of all absolute
disagreement, and the top ten for 87%. That bimodality is exactly what the list-position
mechanism predicts: two campaigns whose lists match for a complex agree on all of it, and two
whose lists diverge early disagree on
most of it. Note these are not rare rows — the three hold 572 of the 4340 shipped single-point
records, 13.2%.

By interface position, mean \|Δ\| is far larger at interface-core positions than at non-interface
ones — expected, since a non-interface mutation's *binding* ΔΔG is ≈0 under any list ordering.
*(Unlike everything else in this section, this split is not recomputed by `agreement_stats.py`,
so the direction is stated and no figures are quoted for it.)*

**Sign stability.** Two campaigns disagree on the sign of `Interaction Energy` for **175 of 5870
comparisons, 3.0%** — but almost all are noise around zero, with a median smaller-side magnitude
of 0.054 kcal/mol. Flips where *both* sides exceed 0.5 kcal/mol: **9 of 5870**. At ±2.0 kcal/mol:
**none**.

**Downstream exposure.** Substituting one campaign's values for another's and re-standardizing per
fold moves a small but non-zero share of FoldX-covered rows, and the **12-term arm is roughly twice
as exposed** as the scalar one, because it reads twelve chances at the noise including the least
stable terms. *(The size of that share depends on the consumer's folds and standardisation,
neither of which this package ships, so the direction is stated and no figures are quoted.)*

None of this affects the shipped store, which is a fixed set of byte-copies — see
[STORE.md](STORE.md), *Provenance of the shipped store*. The exposure is to *regeneration*.

## Mechanism, and the literature

Three independent lines agree, which matters because the vendor never states determinism either
way — no FoldX document claims it, and none denies it.

**1. The parameter semantics.** The entire suite exposes exactly one RNG control:
`--timeSeedRotabase` (bool, **default `false`**) — the rotabase being FoldX's side-chain rotamer
library — documented as "set the seed of the rotabase to
the time (usefull for MC sampling)" [sic]. There is no `--seed` and no `--rngSeed`. An opt-in
flag to seed from the wall clock, off by default, only makes sense if the default seed is a fixed
constant. **Leave it false.** Setting it true would make runs genuinely irreproducible.

**2. The vendor's description of `numberOfRuns`.** From the Switch Lab's own manual: with more
than one run "the algorithm will do the same mutations but changing the rotamer set used and the
order of the rotamer moves". So variation between runs *within* one invocation is deliberate and
indexed — it does not imply variation between invocations. Every campaign here used
`--numberOfRuns=1`, i.e. run 1 every time. The two vendor pages differ on what to set it to:
`foldxsuite.crg.eu` says "Normally it should be set to 1", while the Switch Lab page quoted above
recommends 3. One is followed here; the point that matters for reproducibility is that a run index
is not a random seed, which both pages support.

**3. The direct measurement** in the table above: 1719 pairs with identical structure and
identical list prefix, zero differing. That is the evidence the literature does not contain, and
it is cheap to re-run to confirm on another binary.

**The published magnitudes independently rule out the stochastic explanation.** FoldX's stated
error margin is ~0.5 kcal/mol; its calibration SD against experiment is 0.46 (Guerois, Nielsen &
Serrano 2002, the calibration the force field was fitted on); its *structural* sensitivity — the
spread across different PDB structures of the same protein — is 0.61 (Caldararu, Blundell & Kepp,
*BMC Bioinformatics* 2021). The widest published figure of any kind is a ±3.5 kcal/mol prediction
*interval* for binding ΔΔG (Sapozhnikov, Patel, Ytreberg & Miller, *BMC Bioinformatics* 2023),
which is an accuracy bound rather than run-to-run scatter. **The largest published magnitude of
any kind is more than an order of magnitude below the 11 kcal/mol observed here**, which is itself
an argument that the spread is not FoldX noise.

## Repair count

The FoldX documentation says structures should be repaired before modelling without naming a
count; Usmanova et al. 2018 read that as a recommendation of one, and it is the reading followed
here. CATH-ddG uses one. Usmanova et al.'s study — the only controlled test of repair count in
the literature — found iterating made no difference to *folding* ΔΔG. For *binding* ΔΔG measured here, **the repair
count does move the numbers.**

Measured on 13 complexes, 5× versus 1× repair with identical mutation lists. The 13 are the
CATH-ddG test PPIs carrying at least ten rows — the filter that paper applies for its per-PPI
metrics — and they hold 699 of the split's 813 rows, unevenly: `1JTG_A_B` alone contributes 275.
Rows, not distinct mutations: SKEMPI records several measurements of the same substitution, and
those 275 rows are 194 distinct mutations across both arms, of which 58 are single-point —
which is what this store holds for `1JTG_A_B` in the arm this pilot measures.
**Single-point only**; 270 of the 813 rows, or 33%, are multi-point, so this is the CATH-*single*
number and the full-SKEMPI run covers both arms:

| | |
|---|---|
| structure movement, round 1 → 2 | median **0.25 Å** (converged by round 4) |
| mutations whose `Interaction Energy` changes | **97.6%** of 334 |
| median \|Δ\| / p90 / max | 0.070 / 0.62 / **8.51** kcal/mol |
| FoldX-alone per-structure Spearman, 1× | 0.3984 [0.2759, 0.5193] |
| FoldX-alone per-structure Spearman, **5×** | **0.4585** [0.3455, 0.5695] |
| paired Δ (cluster bootstrap over complexes) | **+0.060, 95% CI [+0.020, +0.100]** |

The published FoldX baseline for single-point mutations on this split is 0.4458 (CATH-ddG
Table 2, FoldX row, per-PPI SpearmanR; the All-mutation figure is 0.4303 and the multi-point
0.5479) — between the two. Repair count is **not** the explanation: CATH-ddG reports one repair
too ([PROTOCOL.md](PROTOCOL.md)), so both sit at ×1. What differs is the engine version (5.0 there,
5.1 here), CATH-ddG's additional `Optimize` step, and the mutation lists each value was computed
in. Which of those accounts for the gap is not determined here.

This does not contradict Usmanova: they measured *folding* ΔΔG and self-consistency bias; this
measures *binding* ΔΔG correlation with experiment.

Caveats worth keeping attached to the number: 13 complexes on one tier, and this is FoldX used
*alone* — whether a downstream model consuming these energies improves on better inputs requires
a retrain. Reproduce or widen with `experiments/repair_ablation.py`.

### Energy convergence against structural convergence

The pilot above compares two endpoints, 1× against 5×. A completed four-round sweep now measures
every step, over 322 single-point and 152 multi-point complexes with the mutation list held
byte-constant, so repair count is the only variable. Measured by
[`experiments/repair_convergence.py`](../experiments/repair_convergence.py) on `Interaction Energy`:

| transition | entries moved (sp) | median (sp) | p99 (sp) | entries moved (mp) | median (mp) | p99 (mp) |
|---|---|---|---|---|---|---|
| 1 → 2 | 92% | 0.0446 | 2.13 | 97% | 0.2227 | 2.95 |
| 2 → 3 | 91% | 0.0201 | 1.70 | 96% | 0.1577 | 2.86 |
| 3 → 4 | 87% | 0.0072 | 1.65 | 96% | 0.1094 | 2.86 |

The structure curve is now measured on the same 344 complexes — the union of the 322
single-point and 152 multi-point sets above — rather than the 13-complex pilot,
from the repair chain that produced these very energies (`repair_round_<r>.pdb`, all-atom RMSD):

| step | structure Å | SP energy mean | SP energy median | MP energy mean |
|---|---|---|---|---|
| 1 → 2 | 0.2663 | 0.2067 | 0.0446 | 0.4672 |
| 2 → 3 | 0.1184 | 0.1479 | 0.0201 | 0.4086 |
| 3 → 4 | 0.0609 | 0.1177 | 0.0072 | 0.3515 |
| **ratio 3→4 / 1→2** | **0.23** | **0.57** | **0.16** | **0.75** |

**`RepairPDB` moves side chains only.** Cα RMSD is *exactly* 0.0000 at every step for all 344
complexes, and Kabsch-superposed RMSD agrees with direct RMSD to four decimals, so there is no
backbone motion and no rigid-body drift. Whatever the repair count changes, it is packing.

**The answer depends on which statistic, and that is the finding.** Against the structure's 0.23,
the *typical* single-point mutation settles faster (median 0.16) while the *average* settles
slower (mean 0.57). Multi-point is slower on both (0.49 / 0.75). Mean and median disagree because
the distribution is heavy-tailed: **6.1% of single-point and 23.4% of multi-point entries still
move by more than FoldX's ~0.5 kcal/mol noise floor between rounds 3 and 4**, and the maximum is
not even monotone — the 3→4 step moves one entry by **11.3 kcal/mol**, more than 2→3's 8.5.

*(The median comparison is the one that reverses depending on which structural curve it is read
against, while the mean comparison holds either way — which is why the claim above rests on the
means.)*

**For ranking, it barely matters:** round 1 against round 4 gives Spearman **0.9295** single-point
and **0.9492** multi-point. So stopping early is defensible for ranking and for the typical
single-point entry, and not defensible for multi-point or for absolute ΔΔG.

**There is no round after which the store stops changing.** Between 87% and 97% of entries differ
from the neighbouring round at every step, including the last one measured.

**Multi-point converges more slowly on every statistic**, which is expected: an entry repacking
two or more positions has more to be sensitive about. Its median at 3→4 (0.109) is still larger
than single-point's at 1→2 (0.045).

The movement concentrates in the same complexes that dominate the campaign-pair disagreement —
at the final step, five complexes hold 37% of all single-point movement, led by `1PPF` and
`3SGB`. That is the third independent measurement pointing at the OMTKY3 saturation sets.

**What this supports is fixing and recording the repair count, not tuning it** — which is what
the store does. It still cannot speak to the 5× endpoint the pilot's +0.060 claim rests on, since
the sweep ran four rounds.

## `--numberOfRuns`

The shipped store uses `--numberOfRuns=1`, the vendor default and the vendor's own documented
recommendation ("normally it should be set to 1"). The FoldX authors' 2025 revision paper uses 5
with the median, but three things qualify that:

1. the reported gain is small (R 0.705 → 0.711, RMSE 1.250 → 1.238) **and was measured on folding
   benchmarks, not on SKEMPI binding data**;
2. it addresses *within-invocation* rotamer exploration, which is **not** the effect measured
   above — it would not remove the mutation-list dependence; computing the union list once does;
3. none of the papers whose baseline rows these results sit beside records a run count at all —
   `numberOfRuns` appears in none of them and in none of their supplements — so whether those
   rows used the single-run default is an inference from the default, not something stated.

What makes a run count comparable across papers is that it is reported, not which value it takes.
A compact template is a pair of sentences from Vincenzi *et al.* (*IJMS* 2024): "As the in
silico protocol started from a Haddock optimized structure, no 'RepairPDB' cycles were performed.
The macros 'PositionScan' and 'BuildModel' were run by using default settings (i.e., number of
runs: 1, pH 7, temperature 298 K, and ionic strength 0.05 M)." In two sentences it states the
repair count, the run count and the environment parameters — including a repair count of zero,
with the reason for it. [PROTOCOL.md](PROTOCOL.md) records what each comparator paper states.

**Platform note.** Both arms of the shipped store come overwhelmingly from the repair sweep,
which ran on Linux — 4286 single-point and 1702 multi-point records. The remaining 119, the
per-definition campaigns for `2C5D`, `3SE3` and `3SE4`, ran on Apple Silicon. There is no vendor
statement on cross-platform bitwise agreement, so the two should not be assumed interchangeable.

What bounds the exposure is that no record is a mixture: each comes from one campaign on one
machine, and the one complex computed on both — `3SE3`'s `B,A` multi-point arm, recomputed on
Apple Silicon from a Linux-built repair — agreed with the sweep on all 30 records to the last
decimal. That is one complex, not a cross-platform guarantee, but it is the only direct evidence
available and it points the right way.
