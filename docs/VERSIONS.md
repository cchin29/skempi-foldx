# Differences between 0.1.0 and 0.2.0

Energies are tied to a release. A rebuilt store is a different measurement rather than a correction
of the earlier one, so values from two releases are not interchangeable within one table.

That statement, on its own, does not tell a reader whose results rest on 0.1.0 how much of their
input moved or whether it moved anywhere that matters. This page quantifies both, and the last
section is the one that answers it.

Almost every figure here is produced by `experiments/compare_versions.py`, and the full report it
writes is [`assets/VERSION_COMPARE.txt`](assets/VERSION_COMPARE.txt). The exceptions, computed
from that report's CSV rather than printed by it, are the pooled cross-arm median shift, the
list-length correlation, and the corrections table below, which also needs 0.1.0's defect
registry.

## Origin of the disagreement

`BuildModel` is deterministic. The same repaired structure and the same `individual_list.txt`
prefix return the same number, to the last decimal, every time. What moves a number is the list:

    a record is a function of (repaired structure, mutation list), not of the mutation

so the same mutation computed in a per-benchmark subset and in its complex's union list are two
different measurements, and neither is the more correct one.

The dependence is measured; the mechanism is not established. FoldX documents each line of
`individual_list.txt` as an independent mutant built from the input structure, which rules out the
mutant structures accumulating, and [DETERMINISM.md](DETERMINISM.md) sets out both the evidence and
what it does not show.

0.2.0 computes every record in its complex's full union list, one list per interface definition.
Records whose 0.1.0 list already matched that come out identical; the rest move.

## Extent of the change

Matched on the SKEMPI identifier, resolving 0.2.0's `role` field as an alias so that 0.1.0's
role-chain keys pair with their author-chain records. "Identical" means all twelve terms are equal,
not merely `Interaction Energy`.

| arm | paired | identical | changed | |
|---|---:|---:|---:|---:|
| single-point | 4232 | 3970 | 262 | 6% |
| multi-point | 1763 | 175 | 1588 | 90% |
| **all** | **5995** | **4145** | **1850** | **31%** |

The asymmetry is the whole story, and it follows from where 0.1.0 got its values. Its multi-point
arm came entirely from one per-benchmark campaign whose lists are not the complexes' union lists,
so almost all of it moved. Its single-point arm was already mostly union-list work: the 0.1.0
release notes identified 287 records, over 27 complexes, whose benchmark list was a strict subset
of the full one, and the 262 single-point records that moved are consistent with that count.

Eight of 0.1.0's 6003 records are not comparable: the eight doubled-code mutations that no longer
resolve by code alone. All eight are listed in the full report. `1S0W AB142F` was a ninth until
the fixed chain-mapping parser gave its record a `role` name; it now pairs, and is identical
across the releases.

## Magnitude of the changes

Shift in `Interaction Energy`, over the records that changed, in kcal/mol.

| arm | n | median | mean | p90 | p99 | max | under 0.5 | sign flips |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| single-point | 262 | 0.052 | 0.251 | 0.738 | 1.875 | 6.854 | 84% | 17 |
| multi-point | 1588 | 0.242 | 0.521 | 1.343 | 3.691 | 9.960 | 67% | 112 |

The distribution is heavy-tailed, so the median understates the risk to any individual row: 69%
of the changed records move by less than 0.5 kcal/mol, and the largest single shift is 9.96 —
`4K71_A_BC GA503R,AA545V`, from 11.05 to 1.09.

A **sign flip** — a record whose energy changes between stabilising and destabilising — is the
change most likely to alter a conclusion drawn from one row. There are 129, 7.0% of the records
that changed and 2.2% of the paired store.

![Single-point values lie on the diagonal; multi-point values scatter around
it](assets/version_compare.png)

Three of the largest shifts are not list effects at all but the interface-definition fix:
`3SE4_B_A RA226A` (0.0000 → 6.8535) and `2C5D_AB_CD LC105E,LD105E` (0.0000 → 5.1289) were scored
in 0.1.0 against an interface excluding the chain they mutate, which returns a clean zero
indistinguishable from a real measurement of no effect. `2C5D_AB_CD RA32E,RB32E`
(0.0977 → 6.9109) is the other shape: only one of its two substitutions lies outside the pairing,
so it returned a plausible partial energy rather than a zero.
Those records were wrong in 0.1.0. They are three of the 26 corrections — the one part of this
comparison that is a correction rather than a different measurement.

### The 26 corrections, listed

A 0.1.0 consumer can check directly whether their rows contain a provably wrong value. These are
the records the 0.1.0 registry flagged as scored against a pairing excluding the chain they
mutate, and whose value moved on at least one of the twelve terms. Five further flagged records
are zero in both releases and are not listed, since nothing about them changed. Fifteen of the 26
were a clean zero in 0.1.0 — indistinguishable from a real measurement of no effect, which is what
made the defect invisible.

**Look these up by the 0.1.0 key, not the 0.2.0 one.** 11 of the 26 are keyed differently in
the two releases: 0.1.0 keys them under the bare code with the other pairing's chain letter, so
`3SE4_B_A DA117A` here is `3SE4 DB117A` in a 0.1.0 store. Matching on the 0.2.0 name alone finds
nothing for those rows and reads as "my data is unaffected" — which is the same silent-miss this
release exists to remove, and it hits the records with the largest corrections.

The `0.1.0 key` column is *same* where the two agree. Energies shown are `Interaction Energy`.

| 0.2.0 record | 0.2.0 mutation | 0.1.0 key | 0.1.0 | 0.2.0 | clean zero |
|---|---|---|---:|---:|:--:|
| `2C5D_AB_CD` | `AA126V,AB126V` | *same* | 0.0000 | -1.5421 | yes |
| `2C5D_AB_CD` | `EC30R,EC33R,ED30R,EC33R` ᵇ | *same* | 4.3550 | 6.2238 |  |
| `2C5D_AB_CD` | `EC30R,ED30R` | *same* | 3.2145 | 4.9774 |  |
| `2C5D_AB_CD` | `EC33R,ED33R` | *same* | 1.4313 | 1.5859 |  |
| `2C5D_AB_CD` | `KA34E,KB34E` | *same* | 0.8310 | 1.9457 |  |
| `2C5D_AB_CD` | `KC178E,KD178E` | *same* | 0.0000 | 3.4407 | yes |
| `2C5D_AB_CD` | `LC105E,LD105E` | *same* | 0.0000 | 5.1289 | yes |
| `2C5D_AB_CD` | `QC28R,QD28R` | *same* | 0.1133 | 0.2273 |  |
| `2C5D_AB_CD` | `QC96R,QD96R` | *same* | 0.0000 | -0.0128 | yes |
| `2C5D_AB_CD` | `RA32E,KA34E,RB32E,KA34E` ᵇ | *same* | 0.6868 | 2.3297 |  |
| `2C5D_AB_CD` | `RA32E,RB32E` | *same* | 0.0977 | 6.9109 |  |
| `2C5D_AB_CD` | `TA26D,TB26D` | *same* | 2.9207 | 6.7259 |  |
| `2C5D_AB_CD` | `TC182E,TD182E` | *same* | 0.0000 | -1.7535 | yes |
| `2C5D_AB_CD` | `TC51R,TD51R` | *same* | -0.6418 | -3.5586 |  |
| `2C5D_AB_CD` | `YA364A,YB364A` ᵃ | *same* | 0.0010 | 0.0010 |  |
| `3SE4_B_A` | `DA117A` | `3SE4` `DB117A` | 0.0000 | 1.4968 | yes |
| `3SE4_B_A` | `FA223A` | `3SE4` `FB223A` | 0.0000 | 1.7884 | yes |
| `3SE4_B_A` | `LA116A` | `3SE4` `LB116A` | 0.0000 | 0.6182 | yes |
| `3SE4_B_A` | `NA140T` | `3SE4` `NB140T` | 0.0000 | -0.0345 | yes |
| `3SE4_B_A` | `NA227A` | `3SE4` `NB227A` | 0.0000 | -0.0896 | yes |
| `3SE4_B_A` | `NA230A` | `3SE4` `NB230A` | 0.0000 | 0.0084 | yes |
| `3SE4_B_A` | `RA226A` | `3SE4` `RB226A` | 0.0000 | 6.8535 | yes |
| `3SE4_B_A` | `SA120A` | `3SE4` `SB120A` | 0.0000 | -0.0279 | yes |
| `3SE4_B_A` | `TA166A` | `3SE4` `TB166A` | 0.0000 | 1.2947 | yes |
| `3SE4_B_A` | `YA148A` | `3SE4` `YB148A` | 0.0000 | 0.0140 | yes |
| `3SE4_B_A` | `YA55A` | `3SE4` `YB55A` | 0.0002 | 1.2967 |  |

ᵃ `Interaction Energy` is unchanged for this record; it moved only in `Solvation Polar`
(0.0011 → 0.0010).

ᵇ These two rows are the entries of `skempi_foldx.MALFORMED_ROWS`. Their 0.2.0 value corrects the
interface they were scored against, but the SKEMPI row itself repeats one substitution where a
symmetric partner belongs — so the energy is right for the three substitutions named and wrong for
the four the row claims, in **both** releases. Correcting the pairing did not make these rows
sound; test for them with `FoldxLookup.is_malformed()`.

The same 31 records — the 26 above plus the 5 that are zero in both — ship as
`skempi_foldx/data/corrections_from_0_1_0.csv`, with both releases' keys, both `Interaction
Energy` values, and a flag for whether any of the twelve terms moved, so this can be joined
against a 0.1.0 store without retyping it. Which *term* moved is not in the file — for the one
record where that matters, footnote ᵃ above says so:

```python
import csv, importlib.resources as res
with res.files("skempi_foldx").joinpath("data/corrections_from_0_1_0.csv").open() as fh:
    corrections = {(r["pdb_010"], r["mutation_010"]) for r in csv.DictReader(fh)
                    if r["moved"] == "True"}
```

Movement is spread across all twelve terms rather than concentrated in one. `Van der Waals clashes`
carries almost the same maximum as `Interaction Energy` itself (9.85 against 9.96), which is what a
clash-dominated rebuild looks like.

The scatterplot matrix puts the quantities the comparison turns on into one field of view, split by
arm, so a relationship that holds in one and not the other is visible rather than argued:

![Scatterplot matrix over both releases' energies, experimental ΔΔG, the shift, clash and
definition size](assets/version_compare_matrix.png)

Three things read straight off it. The `0.1.0 IE` against `0.2.0 IE` panel is the same tight
diagonal with a red cloud around it. The `shift` row is a spike at zero for single-point and a
spread for multi-point, with its widest excursions at high `0.1.0 clash`. And both releases scatter
against `experimental ddG` in the same broad cloud — the panels are hard to tell apart, which is
the visual form of the bootstrap result below.

The complexes carrying the most changed records are `1JTG_A_B` (136 of 194), `2B2X_HL_A` (91 of
94), `3SGB_E_I` (88 of 279), `1CHO_EFG_I` (84 of 275) and `1DAN_HL_UT` (81 of 114). `3L5X_A_HL` is
the most thoroughly rebuilt — 34 of 35 records, at a median shift of 0.932.

## Characteristics of the records that moved

Two separate questions, with two different answers: whether a record changed at all is a property
of its **interface definition**, and how far it moved is a property of the **mutation**.

### The determinant of change: interface definition

A record changed if and only if the list it was computed in changed composition. Since a list
belongs to a definition, the changes arrive by whole definitions rather than scattered across the
store — which is why `3SE4`'s two single-point definitions moved differently from each other.

Drawn one tile per *(arm, definition)* pair, that is the whole picture: 360 of 476 paired tiles
moved in no part, 57 moved entire, and only 59 split — [STORE.md](STORE.md#composition) shows it
against how large each definition is. Tiles rather than definitions, because 130 definitions
appear in both arms and can move in one and not the other; counted per distinct definition the
split is 239 / 9 / 99 over 347.

Throughout this section an "N of M" over records uses M = the records that **paired** with 0.1.0,
not the definition's size in this release; the two differ wherever 0.2.0 holds records 0.1.0 never
had. `4RS1_A_B` is the widest gap below — 36 paired against 49 held. The same distinction is why
the two shipped assets disagree by one on how many single-point definitions hold exactly one
record: `VERSION_COMPARE.txt` says 106, counting paired records, and `STORE_COMPOSITION.txt` says
105, counting what the store holds.

In the single-point arm, all 262 changed records fall in **29 of the 323 paired interface
definitions**. The
other 294 are untouched, which is why that arm's overall change rate is 6%. The largest are
`1DAN_HL_UT` (55 of 85 records), `1IAR_A_B` (29 of 36), `1GC1_G_C` (28 of 54) and `4RS1_A_B`
(28 of 36).

In the multi-point arm the rule is starker still. Every definition that came through unchanged has
**four or fewer variants**, and 46 of the 66 have exactly one — a list of one entry is the same
list in both releases, by construction. Above that, list length predicts the fraction of a
definition that moved with ρ = +0.811:

| variants in the definition | definitions | mean fraction changed |
|---|---:|---:|
| 1–5 | 86 | 0.19 |
| 6–20 | 45 | 0.89 |
| 21–60 | 16 | 0.96 |
| 61+ | 6 | 0.99 |

The practical consequence for anything built on cross-validation folds: the affected complexes are
concentrated, not spread, so how much a fold moves depends on whether it holds one of them. The 29
single-point definitions above and any multi-point definition with more than five variants are the
set to check a split against.

### The determinant of magnitude: side-chain repacking

Among the records that did change, the strongest single predictor of the size of the shift is how
much clash energy the record already carried. Sorting by the 0.1.0 `Van der Waals clashes`
magnitude and quartering:

| quartile of 0.1.0 clash | single-point median \|shift\| | multi-point median \|shift\| |
|---|---:|---:|
| Q1 (least clash) | 0.010 | 0.095 |
| Q2 | 0.039 | 0.188 |
| Q3 | 0.055 | 0.317 |
| Q4 (most clash) | 0.183 | 0.478 |

Two further gradients point the same way. In the multi-point arm the shift grows with the number
of substitutions — median 0.164 kcal/mol at two, 0.263 at four, 0.537 at eight or more, though not
monotonically: the seven-substitution bucket sits at 0.355 against 0.373 at six, on 33 records. In
the single-point arm it grows with the rotamer freedom of the residue being introduced:
substitutions to a side chain with three or more χ angles — the rotatable bonds in a side
chain, so a proxy for how many shapes it can take — move by a median of 0.217 over 34
records, against 0.039 over 228 for two or fewer, and mutations **to alanine or glycine — 170 of
the 262 — move by 0.028**. The per-bucket counts and medians are in
[`docs/assets/VERSION_COMPARE.txt`](assets/VERSION_COMPARE.txt), which
[`experiments/compare_versions.py`](../experiments/compare_versions.py) regenerates; no
significance test is quoted, because no shipped script computes one and a figure a reader cannot
reproduce is worth less here than the gradient itself.

So the records that move are the ones whose energy depends most on how neighbouring side chains
are packed: crowded sites, many simultaneous substitutions, and bulky flexible replacements.
Alanine scans, which have no rotamer freedom to resolve, barely move at all.

That is consistent with the hypothesis [DETERMINISM.md](DETERMINISM.md) offers for the list effect
— rotamer search state carried between entries rather than reset per line — and it is worth saying
plainly that consistency is not confirmation. The mechanism remains unestablished; what this adds
is that the dependence concentrates exactly where side-chain repacking has the most freedom, which
is what that hypothesis would predict.

## Consequences for a model trained on 0.1.0

Disagreement between two stores is not the same question as whether either is a worse input. The
comparison that settles it is each release against experiment, on the records the two share, using
the recipe [STORE.md](STORE.md) states for the shipped agreement figures.

The records identical in both releases score identically by construction, so quoting only the full
population dilutes whatever difference exists. Splitting three ways puts the question where it can
be answered — and leaves a control that must show nothing.

| arm | subset | n | 0.1.0 ρ | 0.2.0 ρ | 0.1.0 r | 0.2.0 r | 0.1.0 MAE | 0.2.0 MAE |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| single-point | all shared | 4053 | 0.435 | 0.436 | 0.419 | 0.418 | 1.138 | 1.138 |
| single-point | **changed** | 257 | 0.374 | 0.392 | 0.416 | 0.393 | 0.904 | 0.909 |
| single-point | identical | 3796 | 0.438 | 0.438 | 0.420 | 0.420 | 1.154 | 1.154 |
| multi-point | all shared | 1576 | 0.574 | 0.587 | 0.515 | 0.518 | 1.885 | 1.871 |
| multi-point | **changed** | 1418 | 0.557 | 0.569 | 0.499 | 0.501 | 1.946 | 1.931 |
| multi-point | identical | 158 | 0.430 | 0.430 | 0.585 | 0.585 | 1.332 | 1.332 |

The identical rows move by exactly nothing, which is the control passing rather than a finding.

On the records that did change, a paired bootstrap over 1000 resamples — records drawn together so
each keeps its pair — puts every interval across zero:

| arm | subset | Δρ | 95% interval |
|---|---|---:|---|
| single-point | all shared | +0.001 | [−0.002, +0.004] |
| single-point | changed | +0.018 | [−0.036, +0.075] |
| multi-point | all shared | +0.012 | [−0.005, +0.029] |
| multi-point | changed | +0.012 | [−0.006, +0.031] |

**The data does not separate the two releases on agreement with experiment.** The multi-point
Spearman rises in the direction a union-list value should, and the interval is not narrow enough
to call it. The single-point changed subset moves both ways at once — Spearman up, Pearson *down*
from 0.416 to 0.393 on 257 records — which is what noise at this sample size looks like rather
than a signal with a sign.

So the case for 0.2.0 rests on what is demonstrable: the interface-definition fix, which corrects
records that were provably scored against the wrong pairing, and one list composition per
definition, which makes a value mean what the store says it means. Better agreement with experiment
is not part of that case, and quoting the 0.574 → 0.587 as an improvement would overstate it.

As a feature column, the two releases track each other closely:

| arm | n | Pearson r | Spearman ρ | MAE |
|---|---:|---:|---:|---:|
| single-point | 4232 | 0.996 | 0.994 | 0.016 |
| multi-point | 1763 | 0.945 | 0.944 | 0.469 |

**A single-point model trained on 0.1.0 features is not invalidated by 0.2.0.** Rank correlation of
0.994 at a mean absolute difference of 0.016 kcal/mol is below the resolution at which a learned
model distinguishes inputs, and the agreement with experiment does not move.

The multi-point arm needs the weaker statement. At ρ 0.944 and MAE 0.469, with 112 sign flips, a
result that leans on multi-point FoldX features is measurably sensitive to which release produced
them, and a retrain is what establishes the size of that sensitivity rather than bounding it from
these numbers. Nothing here predicts which way such a retrain moves: on the records that changed,
the paired bootstrap puts the change in agreement with experiment across zero in both arms.

Either way, what a published result needs is the release it used, stated. Both are tagged.

## Reproducing this

Everything except the comparison against experiment comes from the repository itself. The older
release is a `git worktree` at its tag, so nothing needs downloading:

```
git worktree add ../skempi-foldx-v0.1.0 v0.1.0
python experiments/compare_versions.py --old ../skempi-foldx-v0.1.0 \
    --out scratch/version_compare --csv --plots --html
```

`--skempi-csv` is the one input the repository cannot supply, SKEMPI's table not being
redistributable from here. It adds the two sections that compare each release against measured
ΔΔG, and the CSV's `ddg_experimental` column, and it supplies the `experimental ΔΔG` variable in
the scatterplot matrix; without it every other figure and table on this page still reproduces, and
the report says so where the table would have been.

That column is derived, not read: SKEMPI records dissociation constants and a temperature and has
no ΔΔG field, so the value here is `RT·ln(Kd_mut / Kd_wt)` under the recipe
[STORE.md](STORE.md#agreement-with-experiment) states. A different recipe gives different numbers
from the same table.

No FoldX binary is needed either way, since this compares shipped values rather than recomputing
them. `--plots` needs matplotlib and `--html` needs plotly, which the `figures` extra installs
(`pip install -e ".[figures]"`); both flags are optional and skip silently when the library is
absent, and the report and the CSV are written regardless.

`--csv` writes one row per 0.1.0 record — 6003 of them — carrying every covariate derived here:
both energies and all twelve terms for each release, the shift, whether the record changed or
flipped sign, its definition's size, its substitution count, its source in each release, and the
experimental ΔΔG where SKEMPI has one. Every figure on this page is recomputable from that file
alone, which is the point of it: a reader's question is usually about their own subset, and
answering it should not require re-deriving the pairing, the role aliasing or the experimental
join. The eight unpairable records are included and marked, rather than dropped into a file that
would otherwise read as complete.

```python
import pandas as pd
d = pd.read_csv("version_compare.csv").query("status == 'paired'")
rows = d[d.identifier.isin(split)]             # which of these rows moved, and by how much
rows.groupby("arm").agg(n=("changed", "size"), changed=("changed", "sum"))
```

`--html` writes an interactive version of the figure above, in which every point carries its
identifier, mutation, both energies, the shift, whether it flips sign, and the two variables that
predict the shift — the size of its definition's list and the number of substitutions. That
answers the question a static scatter cannot: not how much of the store moved, but whether a
particular row did. The file is self-contained, works offline, and is around 5 MB, which is why it
is generated on demand rather than committed.
