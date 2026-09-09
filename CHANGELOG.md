# Changelog

Energy values are tied to a version. A mutation's ΔΔG depends on the whole mutation list it was
computed in, so a rebuilt store is mostly a different measurement rather than a correction of the
same one — values from two versions are not interchangeable within one table. Where a release does
correct a value rather than recompute it, the entry says so and counts them. Entries below say
explicitly whether a release changes them.

## 0.2.2 — 2026-09-09

**First release on PyPI.** `pip install skempi-foldx==0.2.2` is now the documented install.
Earlier releases are installable only from their git tags, which is how every release before this
one was distributed.

**No data changed.** The shipped store is byte-identical to 0.2.0 — same 4340 single-point and
1767 multi-point records, same values on all twelve terms, verified by the digest the suite pins
over every shipped `(arm, identifier, mutation, term, value)`. A version names a measurement in
this package, so a version that moves no measurement says so explicitly: **0.2.0 and 0.2.2 are the
same numbers**, and a result computed against either stays correct without qualification. Nothing
here supersedes 0.2.0 as a measurement; it supersedes it only as a distribution.

### Why 0.2.2 and not 0.2.1

A `v0.2.1` tag was published in August and withdrawn, and `CHANGELOG.md` records that a pin naming
it now fails outright. Reusing the number would have made that true statement false and would have
silently re-pointed any surviving `@v0.2.1` pin at unrelated content — the same defect this
release series exists to correct. 0.2.1 is skipped.

### Changed

- The install instructions name PyPI.
- **Every link in `README.md` is absolute and pinned to this tag.** The README is the long
  description, so it is also the PyPI project page, and PyPI's renderer leaves a relative path
  alone: `](docs/STORE.md)` resolves against the project page URL and 404s. Thirty-two links --
  every `docs/` reference, `CHANGELOG.md`, `NOTICE`, `CITATION.cff` — were relative, correct on
  GitHub and dead for anyone arriving from PyPI. They name `blob/v0.2.2/` rather than `blob/main/`
  so an old release's page keeps showing that release's claims. In-page anchors were already
  safe; the renderer rewrites those. A test now fails on any relative link, on a pinned link
  naming another version, and on one naming a path not in the tree.
- **The install-identity command no longer raises on an install from a file.** It read
  `direct_url.json['vcs_info']['commit_id']`, guarded only against the file's absence — but a
  wheel or sdist installed from a path writes `archive_info` instead, so the command raised
  `KeyError` for a case the surrounding paragraph named as handled. It now reports the commit for
  a git install, the archive for a file install, and `installed from PyPI` when there is no
  `direct_url.json` at all.
- `CONTRIBUTING.md` describes publishing to an existing project rather than creating one, since
  the trusted publisher is now registered and the pending form no longer applies.
- `tests/test_version_consistency.py` now requires a `skempi-foldx==<version>` pin in the README
  and checks it tracks `pyproject.toml`. Until this release it asserted the opposite, because a
  PyPI pin was a command that failed for every reader.

## 0.2.0 — 2026-09-08

**An earlier `v0.2.0` tag, and a `v0.2.1` that existed briefly, have been withdrawn and replaced
by this release.** The two behave differently. `v0.2.1` no longer exists, so a pin naming it fails
outright; an install from it also reports `__version__ == "0.2.1"`, a version this changelog does
not record, and moving off it lowers the version string while raising the content. `v0.2.0` still
exists but now points at this tree, so a pin naming it still resolves — silently, to different
content than before, and `__version__` does not distinguish the two. Neither tag reached PyPI,
which had no release of this package until 0.2.2. No energy value differs
between those trees and this one,
so a result quoting the numbers stays correct; what they lacked was the `role` field on 16
records, described under [chain
mappings](#chain-mappings-and-16-role-names-lost-to-a-fused-column) below, and corrected metadata.
`__version__` does not distinguish them — for a git install the `direct_url.json` in an installed
distribution records the commit, and does.

The `v0.2.0` tag has been moved more than once: it was re-pointed on 2026-09-09 after this
release's documentation was corrected, and the earlier object is unreachable. **No energy value
differs across any of those moves** — the shipped store is byte-identical, and the suite now pins
a digest over every value so that a silent change is impossible within a version. What changed was
documentation, tests and CI. An install taken from `v0.2.0` at any point holds the right numbers;
reinstall only to get the corrected documentation and the packaging fixes.

To move an install off any of them, pip needs telling, since the version already matches:

```bash
pip install --force-reinstall \
    "git+https://github.com/cchin29/skempi-foldx@v0.2.0"
```

**Energy values: 1850 of 6003 changed, and 4145 did not**, and eight could not be compared.
26 of those changes are corrections — the interface-definition fix below, where 0.1.0 scored a
record against a pairing excluding the chain it mutates. Five more such records were registered as
defective but are zero in both releases, so they are not among the 1850: the registry counts 31
defective records, of which 26 moved. All 31 are under `2C5D` and `3SE4`: `3SE3` is the third code
SKEMPI defines twice, but its second definition contributes no mutation the first lacks, so
0.1.0's pooling under the first stranded nothing there. `3SE3`'s own values do change between the
releases, but as ordinary recomputation rather than as corrections. The other 1824 changes are not
corrections: the two releases disagree there because a FoldX energy is a property of the run that
produced it, not of the mutation alone, and neither run's value is the more correct one.

`BuildModel` is deterministic: the same repaired structure and the same `individual_list.txt`
prefix return the same number, to the last decimal, every time. What moves a number is the list it
was computed in —

    a record is a function of (repaired structure, mutation list), not of the mutation

— so the same mutation computed in a per-benchmark subset and in its complex's union list are two
different measurements, and neither is the more correct one. The dependence is measured; the
mechanism is not established, and FoldX's own documentation describes each line of the list as an
independent mutant built from the input structure. [docs/DETERMINISM.md](docs/DETERMINISM.md) sets
out the evidence and is explicit about what it does not show.

0.2.0 computes every record in its complex's full union list, one list per interface definition.
Records whose 0.1.0 list already matched came out identical across all twelve terms; the rest
moved. The split is very uneven between the arms, because 0.1.0's single-point values already came
mostly from union lists while its multi-point arm came from a per-benchmark campaign:

| arm | unchanged | changed | |
|---|---:|---:|---:|
| single-point | 3970 | 262 | 6% |
| multi-point | 175 | 1588 | 90% |

Eight of the 6003 are not comparable: the eight doubled-code mutations that no longer resolve by
code alone. `1S0W AB142F` was a ninth until the fixed chain-mapping parser gave its record a
`role` name; it now pairs, and is identical across the releases.

Among the 1850 that moved, the median shift in `Interaction Energy` is 0.208 kcal/mol and 69% are
under 0.5 — but the tail reaches 9.96, at `4K71_A_BC GA503R,AA545V` (11.05 → 1.09), and 129
records change sign between stabilising and destabilising.

The changes are concentrated rather than scattered: all 262 single-point ones fall in 29 of the
323 paired interface definitions, and every multi-point definition that came through unchanged
holds four or fewer variants. Within a changed list, the size of the shift tracks side-chain
repacking freedom — crowded sites, many simultaneous substitutions, bulky flexible replacements —
while alanine substitutions barely move.

Against experiment the two releases are equivalent on the records they share: single-point Spearman
0.435 → 0.436, multi-point 0.574 → 0.587, on 4053 and 1576 rows. Restricted to the records that
actually changed, a paired bootstrap puts every interval across zero, so **the data does not
separate the releases on agreement with experiment** — the multi-point figure moves in the
direction a union-list value should, and not by enough to claim. As a feature column they track at
ρ 0.994 (single-point, MAE 0.016) and ρ 0.944 (multi-point, MAE 0.469). **A single-point model
trained on 0.1.0 is not invalidated by 0.2.0**; a multi-point one is measurably sensitive to which
release produced its features, and a retrain is what settles the size.

The case for 0.2.0 is the interface-definition fix and one list composition per definition, not
better correlation with experiment.

**Do not mix 0.1.0 and 0.2.0 values in one table.** Two thirds of the store is unchanged, so a
mixed table looks self-consistent while a third of its rows are a different measurement. Pin a tag.

Record-by-record counts, magnitudes, the scatter and the reproduction recipe are in
[docs/VERSIONS.md](docs/VERSIONS.md).

### The record unit — the SKEMPI identifier

A record is now named `<pdb>_<group1>_<group2>` — `1BRS_A_D`, `2C5D_A_C`, `2C5D_AB_CD` — because
that is what SKEMPI names a row by, and three codes appear under two such names. 0.1.0 keyed on the
code, kept whichever pairing came first, and pooled the rest under it, so 31 records were scored
against an interface excluding the chain they mutate. That failure was invisible in the output: a
mutation outside the analysed pair moves neither side of `IE(mut) − IE(wt)` and returns a clean
zero, indistinguishable from a measurement of no effect.

Across the three codes, 26 of 0.1.0's records were zero on all twelve terms — a different set of
26 from the corrections counted above, overlapping it in 15 members. Recomputed against their own
pairings, **15 of those 26 carry real signal**, up to 6.85 kcal/mol — `3SE4_B_A RA226A` goes
0.0000 → 6.8535.

`load_skempi` is keyed by identifier and `SkempiComplex` carries `.id`; `by_pdb()` regroups for the
steps that are per structure. `FoldxLookup` resolves a bare code wherever it names one record —
ambiguity is per `(code, mutation)`, not per code, and affects **8 of the 111** mutations under
those three (`2C5D`'s two definitions share none at all). For the 8, `get` returns its default and
`require` raises naming both identifiers rather than choosing. `definitions_of(code)` lists them.

**Removed**: `COLLAPSED_INTERFACES`, `CollapsedInterface`, `N_INTERFACE_SUSPECT`,
`interface_suspect()`, `FoldxLookup.is_interface_suspect()`, `SkempiComplex.alternate_groups` and
`mutations_outside_interface()`. They described a defect this layout does not have, and a guard
that can no longer fire is worse than no guard — it reads as a check that passed.

**Added**: `SkempiComplex.id`, `by_pdb()`, `pool_by_code()`, `code_of()`,
`FoldxLookup.definitions_of()`, `FoldxLookup.is_malformed()`, `MALFORMED_ROWS`, `MalformedRow`,
`is_malformed()` and `announce_malformed()`.

**Added — `skempi_foldx/data/corrections_from_0_1_0.csv`.** The 31 records the 0.1.0 registry
flagged as scored against a pairing excluding the chain they mutate: the 26 whose values moved and
the 5 that are zero in both releases. It carries *both* releases' keys, because 11 of the 26 are
keyed differently — 0.1.0 names them under the bare code with the other pairing's chain letter, so
`3SE4_B_A DA117A` here is `3SE4 DB117A` there — and a join on the 0.2.0 name alone finds nothing
for exactly the records with the largest corrections, which reads as "unaffected". It ships in
both the sdist and the wheel, so the question "is one of my rows provably wrong?" is answerable
from an install, without a checkout:

```python
import csv, importlib.resources as res
with res.files("skempi_foldx").joinpath("data/corrections_from_0_1_0.csv").open() as fh:
    corrections = {(r["pdb_010"], r["mutation_010"]) for r in csv.DictReader(fh)
                   if r["moved"] == "True"}
```

Also added, all for asking a question that previously only a warning answered:
`check_identifier()` and `check_code()`, which validate a SKEMPI identifier and a PDB code where
either becomes a filename; `FoldxLookup.conflicts`, the names that resolved to more than one
record, one entry per additional claimant as `ConsolidationReport.summary` counts them; and
`FoldxLookup.pooled_codes`, the codes indexed both as a record key and as interface definitions.

That last one names the single route by which a bare doubled code can still answer. A store
keyed by PDB code — every store this package shipped before 0.2.0 — pools a code's definitions
under one name, and an exact key hit resolves before the ambiguity refusal is consulted, so
mixing that store with an identifier-keyed one returns the pooled record rather than the
pairing-specific one. Which record the caller wanted is not decidable here, so the condition is
reported rather than resolved. Loading one store convention at a time avoids it entirely.

**Changed — the agreement check no longer fails open.** `audit()` and `consolidate()` decide
whether two sources agree from `max|Δ|` across the shared terms, and two comparisons that cannot
be made were reading as comparisons that succeeded. A non-finite term returned `nan`, and every
comparison against `nan` is False, so a source pair disagreeing `1.0` against `NaN` was reported
as agreeing; `json.loads` accepts bare `NaN`, so a foreign results directory can supply one. Two
records sharing no term at all returned `0.0`, reported as identical. Both now return infinity.
`term_vector` refuses a non-finite term for the same reason, which `FoldxLookup.vector` already
did — the two are documented as the same operation and must not disagree about what a valid
record is.

**Changed — a `mapping_dir=` that contradicts the store is refused rather than believed.** A
record that ships a `role` already carries this computation's known-good answer. When the
recompute from a supplied mapping directory disagrees, the mapping files describe different
residue numbering than the ones that built the store, and installing the recomputed name would
alias a label belonging to another mutation — answered confidently and counted as covered. The
shipped name now wins, `FoldxLookup.role_mismatches` lists the disagreements and a warning names
them. With correct mappings the recompute is exactly idempotent, so this cannot fire on good
input. `load_chain_mapping` also now requires a single-character chain and a numeric residue
number, since a header row previously parsed as a chain of its own and shifted every later
offset — one way a mapping set comes to disagree in the first place.

**Changed — `reindex_by_skempi_id` announces a collision instead of dropping a record.** Two
records re-keying onto one name is a collision, not a merge; the first is kept, as `consolidate()`
does, and the rest are reported. It stays an identity map on the shipped store. The store-level
`coverage()` docstring now says plainly that it resolves no naming — on the bundled store it
reports 19.6% where `FoldxLookup.coverage()` reports 100% over the same rows.

**Changed — identity-based role inference no longer guesses a residue number.** With
`skempi_csv=` but no `mapping_dir=`, `FoldxLookup` infers a record's role-chain name by matching
substitution identity. The role convention renumbers each chain by the cumulative residue count of
the chains before it in its group; that inference has no residue counts, so the author position is
the role position only for the first chain of a group. For any other chain it produced a label in
neither convention — the role chain letter with the author number — which named no record and so
installed an alias onto whichever record carried that substitution. On the bundled store that was
863 aliases, each answering a role label that belongs to a different mutation, counted as covered,
with `is_exact` still reporting `True`. Such labels are now refused, leaving them unresolved;
`mapping_dir=` resolves them exactly, as it always did.

**Changed — every name that becomes a path is validated where it becomes one.** A record
identifier and a PDB code are both parsed from column 0 of a caller-supplied `skempi_v2.csv`,
which nothing upstream constrains, and both reach the filesystem: `write_store` and
`FoldxConfig.result_path`/`complex_work_dir`/`find_repaired` build filenames from an identifier,
`run.py` and `load_chain_mapping` build them from a code. Each validates at the join rather than
at the parse, because the join is the one place every route passes through.
`experiments/build_per_definition_store.py` checks the whole store before it removes the previous
one, since a removal partway through would cost the data the check exists to protect.

### Contents

| | definitions | codes | records | of SKEMPI | excluding `1KBH` |
|---|---:|---:|---:|---|---|
| single-point | 323 | 322 | **4340** | 4340 / 4343 | **100%** |
| multi-point | 154 | 152 | **1767** | 1767 / 1850 | **100%** |

The entire shortfall is `1KBH` — 3 single-point and 83 multi-point — whose `RepairPDB` does not
terminate under FoldX 5.1. No amount of compute closes it. Quote the denominators with the
exclusion stated; a bare "of 1850" implies a target that does not exist.

The denominators count per definition, which is why they differ from 0.1.0's by-code 4337 and 1848.

### One key convention, and a shipped role name

Every key is SKEMPI's author-chain form: 4340 and 1767 records with `key == cleaned`, no
exceptions. 0.1.0's single-point arm split 3012/1226 between that form and the role-chain form;
its multi-point arm was already all author form, so the mixture was 4777/1226 across the store.

Because one convention alone would strand role-labelled split files, a record whose role name
differs from its key now carries that name in a **`role`** field, computed from the chain mappings
at build time — 3887 of 6107 records. For all 2220 of the rest the role form *is* the key, because
the chain sits first in its group and so carries offset 0. Asked by identifier, the two together
recover every one of 0.1.0's 1226 role keys.

A field rather than a second key convention, deliberately: `skempi_v2.csv` is not redistributable,
so a consumer resolving role labels at read time would need the one input they cannot have.

### Chain mappings, and 16 role names lost to a fused column

`load_chain_mapping` reads each chain's residue count from a SKEMPI `.mapping` file, and
`to_role_form` accumulates those counts into the offsets a role-chain name is built from. The
files are fixed-width, so **an author number of 1000 or more runs into the chain letter** and the
line arrives as three whitespace fields rather than four — `ALA C1001  1`. Requiring four fields
dropped every such residue. Eleven of 345 mapping files carry them.

Two consequences, both silent:

* `2NYY` and `2NZ9` chain A measured **971 residues instead of 1267**, so `to_role_form` rejected
  14 in-range mutations as out of range.
* `1S0W` and `1XXM` lost the **whole** of chain C, whose residues are numbered 1001–1165, so the
  group loop bailed before any arithmetic and 2 more records got no name.

Neither raises. A role-keyed join simply missed those 16 records, since `None` is also how
`to_role_form` refuses an insertion code. Worse in principle than in fact: for a *later* chain in
a group, a short count does not decline anything — it shifts every subsequent offset, aliasing
records onto residues they do not describe. The field is now split back apart in
`load_chain_mapping`, where the damage starts, rather than guarded against downstream.

The 16 names were added to the shipped store by `experiments/backfill_role_names.py`, which
derives them from the interface identifier alone and so needs neither the mapping files nor
`skempi_v2.csv`. It is an audit script, not a CI gate; what CI holds is the store, through
`test_shipped_store_matches_the_backfill_rule`.

Verified two ways rather than argued. Recomputing the role name of **all 6107** shipped records
with the fixed parser reproduces every one — 0 additions, 0 losses, 0 disagreements — so the store
and its build script agree again. And against the withdrawn `v0.2.0` tag the store gains 16 role
fields and **changes no value**: 6107 records compared field by field, 0 non-`role` differences.
The same holds against the tree the withdrawn `v0.2.1` tag pointed at, whose store is identical to
the first's. Neither tree is reachable from any tag now, so the comparison is stated here rather
than reproducible.

### Coverage figures and the interface identifier

`experiments/coverage_report.py` reduced a split's sequence id to its bare PDB code, discarding
the chain groups the id carries. That made `3SE4` and `3SE3` — the codes SKEMPI defines under two
pairings — unresolvable, since a bare code then has two right answers, and it *conflated* the two
pairings' rows in the denominator as well. It now reads `1A22.A.B_A` as `1A22_A_B`.

**That fix alone accounts for the whole change**: the splits this store was built for go from
4153/4159 and 1632/1634 to full coverage at 4165 and 1636, with the role names making no
difference to those two figures. The 16 recovered names move a different split — a 1100-row
S1102 partition goes from 1099 to 1100, on `1S0W`.

### Provenance

Sources are the repair sweep's `round_1` for the 468 definition-arms whose code SKEMPI defines
once, and nine per-definition campaigns for `2C5D`, `3SE3` and `3SE4`, from
`recompute_alternate_interfaces.py --all-definitions`. `round_1` cannot serve those three: its
worklist pools both definitions into one list, so its values belong to no single pairing's list.

**Mixing the two is measured, not assumed.** `3SE3`'s `B,A` multi-point arm was recomputed as a
per-definition campaign seeded from the *original* campaign's repaired structure, while `round_1`
came from the *sweep's own* repair chain: 30 of 30 records agree to the last decimal, on an
identical mutation-list hash.

Per-record `_source` names the campaign. `CONSOLIDATION_REPORT.txt` is rewritten to match — each
record has exactly one source by construction, so there is no precedence order and no conflicts to
arbitrate. The report ships as written by `experiments/build_per_definition_store.py`; there is no
regeneration script, deliberately, since one pointed at this store would overwrite the shipped
report with a description of a store that no longer exists.

Against SKEMPI's measured ΔΔG the shipped values reach Spearman **0.437** single-point (4081 rows)
and **0.587** multi-point (1580), deduplicating on `(identifier, mutation)` and keeping the first
row. Keeping the last gives 0.439 and 0.589 — worth stating, since the two differ in the third
decimal and nothing in the data says which convention a reader assumed.

### Substitution order in a multi-point label

Six multi-point *mutants* are written two ways — SKEMPI lists both orderings, and both are
computed and shipped, so these six mutants account for twelve of the 1767 records:

| record | spellings | \|Δ Interaction Energy\| |
|---|---|---:|
| `1JTG_A_B` | `DB49A,KA208A` / `KA208A,DB49A` | 0.499 |
| `1JTG_A_B` | `DB49A,RA217A` / `RA217A,DB49A` | 0.647 |
| `1JTG_A_B` | `DB49A,SA105A` / `SA105A,DB49A` | 0.481 |
| `2WPT_A_B` | `DA30A,FB79A` / `FB79A,DA30A` | 0.639 |
| `3S9D_A_B` | `EB66A,RA120A` / `RA120A,EB66A` | 0.014 |
| `4L3E_ABC_DE` | `RA65A,YD26A` / `YD26A,RA65A` | 0.078 |

They are not duplicates of each other: order is part of the mutation list, so the two spellings
are two positions in that list and the energies differ accordingly — by up to 0.647 kcal/mol here.
A consumer canonicalising a label by sorting its substitutions collapses each pair and keeps an
arbitrary one of the two. Neither is more correct; if you canonicalise, do it consistently and say
which you kept.

### Malformed SKEMPI rows

`MALFORMED_ROWS` is a trust guard on computed values, not a compute guard. Nothing about the two
`2C5D` rows resists computation: the string names three substitutions, FoldX returns the energy of
those three, and every comparator pipeline whose source was read computes the same mutant. What is
defective is the arity the label claims — 3-substitution mutants recorded as 4.

So they are computed, shipped and flagged rather than dropped. Dropping them would leave
`is_malformed()` and the load-time warning guarding an empty set, and would move the other 17
values in that complex's list, since a record is a property of the list it was computed in.
`announce_malformed` names the rows and returns them without touching the list.

The registry is new in 0.2.0, so no 0.1.0 consumer has a symbol to migrate. It is recorded here
because the two rows are shipped in both releases and a reader needs to know what they are.

### Findings from the 0.1.0-era campaign directories

Four results computed from the pipeline's own campaign output — the `scratch/` directories a run
leaves behind, which are not shipped. The first two answer questions the 0.1.0 documentation left
open; the last two came out of the same material:

- `CONSOLIDATION_REPORT.txt` became regenerable again, which established one correction to the
  0.1.0 report: its conflict line counted source-pair rows and called them mutations, so "1970
  mutations present in more than one source" should have read "1027 mutations … giving 1970
  source-pair disagreements". The overstatement was ~1.9×. The report this release ships is
  written by `experiments/build_per_definition_store.py` and describes a store with one source
  per record, so it has no conflict line at all.
- The repair-count question is **answered on both halves**. Energy by
  `experiments/repair_convergence.py`; structure from the sweep's per-round `repair_round_<r>.pdb`
  chain. `RepairPDB` moves side chains only (Cα RMSD exactly 0.0000 at every step, all 344
  complexes). Against the structure's 0.23 convergence ratio, the *typical* single-point mutation
  settles faster (median 0.16) and the *average* settles slower (mean 0.57); multi-point is
  slower on both. Ranking is robust regardless — round 1 vs round 4 Spearman 0.9295 / 0.9492.
- **Two SKEMPI rows are malformed upstream**, both in `2C5D_AB_CD`, both shipped in 0.1.0:
`RA32E,KA34E,RB32E,KA34E` and `EC30R,EC33R,ED30R,EC33R` each repeat a substitution instead of
applying the symmetric partner (`KB34E`, `ED33R`). Verified against `skempi_v2.csv` — 2 of 7085
rows, the only such rows in the file. They are 3-substitution mutants labelled as 4-substitution.
**Registered** in `skempi_foldx.MALFORMED_ROWS` at mutation granularity and deliberately not
corrected. Reading the published source of the comparator pipelines shows no pipeline drops them
and all compute the same 3-substitution mutant, so the values here are not anomalous — the
discrepancy is confined to the arity label, which the surveyed pipelines record inconsistently (3
in USP-ddG, 4 in RDE-Network and MINT) while all compute the same three-substitution mutant.
- `docs/DETERMINISM.md`'s agreement statistics are recomputed and no longer marked unverified,
by `experiments/agreement_stats.py`. **A comparison count of 10631, quoted in 0.1.0's
`docs/DETERMINISM.md`, is not reproducible from the campaign directories and has been retired**;
the single correct population is 1527 shared mutations over 5870 comparisons, which the
determinism table already described. Every quantity in that table moved — pooled RMSD 0.45 →
0.524, ICC 0.96 → 0.972, sign flips 2.4% → 3.0% — and the least-stable term is entropy mainchain
(0.923), not torsional clash. The qualitative conclusions are unchanged.

### Fixed — claims about other groups' work, published in 0.1.0

These were measured or attributed wrongly in 0.1.0's documentation, which remains public. Each is
restated as what can be checked against the cited source:

- The CATH-ddG test split was described as 39% multi-point. It is **33%** — 270 of 813 rows,
  recomputed from that project's own released split. The 13 complexes the repair-count measurement
  runs on were labelled a "T≥10 tier"; that paper defines no such tier, and the threshold is its
  per-PPI metric filter of ten rows per PPI. They hold 699 of the 813, `1JTG_A_B` alone
  contributing 275 — rows rather than distinct mutations, a distinction the text now draws, since
  those 275 rows are the 194 mutations this store holds for that definition.
- The repair sweep's rationale claimed 5× repair moved the per-structure Spearman "past the
  published single-point FoldX row of 0.4458". The 5× interval is [0.3455, 0.5695] and spans that
  value rather than resolving against it — which is what `docs/DETERMINISM.md` already said, and
  what the sweep's own summary now says.
- The repair-count discussion asserted that two papers reporting gains on binding data "both
  confound repair count with other pipeline changes", naming and citing neither. Restated as what
  can be checked: no report isolates repair count from the other changes made alongside it.
- The comparator survey counted three works reporting a FoldX baseline that appears elsewhere.
  There is a fourth outside the ΔΔG tables: PPIformer reports no FoldX ΔΔG row, but its
  antibody-retrieval table carries one and states that it is reproduced from RDE-Network. The
  survey's own thesis is that a reported FoldX baseline is not always a fresh run, so an em-dashed
  row was hiding an instance of it. Of the three already counted, two were given an innocent
  reading and the third now carries the caveat the section's own sourcing note implies: a value
  may be attributed in a section not read.
- A comparator that dedupes without a comment was described as doing so "silently" beside one that
  does it "deliberately". The evidence in the same paragraph shows the behaviour is load-bearing
  either way, so the adverbs were doing work the data does not support.
- GearBind states "we use FoldX 4 for mutant structure generation" — the only work surveyed to
  name a major version other than 5, and a comparability axis the survey raises for CATH-ddG's 5.0
  against this store's 5.1. It went unrecorded.
- CATH-ddG's `RepairPDB` count was attributed to supplementary section 1.4, which prints the
  command lines but states no count. The quoted string is in section 1.11, on computational
  efficiency.
- RDE-Network's widely-quoted 0.3789 is a per-structure **Pearson**. Every other correlation in
  these documents is a Spearman, so leaving the metric unnamed invited the wrong reading.
- The FoldX documentation says structures should be repaired before modelling without naming a
count; "FoldX's manual recommends one repair" is Usmanova et al.'s reading of it, and is now
attributed to them. The `--numberOfRuns` discussion likewise now records that the vendor's two
pages differ — one says "normally it should be set to 1", the other recommends 3.
- Two figures carrying the argument that the observed spread is not FoldX noise — the ±3.5 kcal/mol
binding interval and the 0.46 kcal/mol calibration SD — were quoted without sources. Both are now
in `docs/REFERENCES.md`, which claims to cover every work named under `docs/`.
- A claim that every published baseline row "was produced at the single-run default" was an
  inference from the default presented as knowledge. What is checkable, and now stated, is that
  none of those papers records a run count at all.

### Fixed — citation metadata

- `CITATION.cff` carried `license: MIT AND CC-BY-4.0`. That is an SPDX *expression*, and CFF 1.2's
`license` field takes an identifier or an array of them — an array the schema reads as **OR**, so
neither spelling states the licence this work is under. The field is omitted rather than guessed
at: `LICENSE`, `NOTICE` and `pyproject.toml`'s `MIT AND CC-BY-4.0` are the record, and CFF cannot
express it. Dropping it is deliberate, not an oversight — an array here would publish that a
consumer may take the whole work under MIT alone, discarding the CC BY attribution that `NOTICE`
states as a licence term rather than a courtesy, on data derived from SKEMPI 2.0.
- Two works cited in the prose were absent from `references` — the Protein Data Bank, and
  RDE-Network among the surveyed comparator pipelines.

### Added

- Release-metadata tests, `tests/test_version_consistency.py`. The version is written in
  `__version__`, `pyproject.toml`, `CITATION.cff`, this file's newest heading and three times in the
  README; the supported Python range is written in the `requires-python` specifier, the
  classifiers and the CI matrix. Nothing checked either group, and an unchecked version
  declaration is how a tag comes to describe a tree it does not match.
- Among them, `test_readme_install_pin_matches_pyproject`: the install command is a version
  declaration too, and the one with the shortest path to a reader. For a package distributed by
  git tag, a stale `@vX.Y.Z` there does not merely misdocument, it installs the previous release.
- `test_dev_extra_names_the_parsers_this_file_needs`. Without `tomli` and `pyyaml` most of that
  suite is import-guarded into silence — a green bar reporting that nothing ran. Read as raw text,
  since using the TOML parser to check that the TOML parser is declared would skip in exactly the
  case worth catching.
- `test_requires_python_excludes_nothing_it_advertises`, for the `!=3.x.*` form. A floor and a
  ceiling cannot express it, so an exclusion would leave three files disagreeing with no bound out
  of place to notice.
- A CI step refusing a tag whose name disagrees with the version in the tree it points at. The tag
  is what publish.yml uploads, so it is the one version declaration that lives outside the
  tree and could not be checked from inside it.
- The bundled multi-point definition count and the role-name alias are asserted in the suite as
  well as in the wheel-build step of both workflows. There they read as duplicates of assertions
  made elsewhere and are not: the wheel check answers whether the store survived packaging, which
  no test of the source tree can.
- Continuous integration, `.github/workflows/tests.yml` — 0.1.0 shipped no `.github/` at all.
  Six interpreters, and it runs on version tags as well as `main` and pull requests: a tag can
  point at a commit the branch never ran, and for a git-tag-installed package that tag is what a
  consumer resolves. Every action in both workflows is pinned by
  commit SHA, with the version in a trailing comment: a tag is a mutable ref, and an action that
  moves under a workflow holding `id-token: write` moves what can mint a PyPI token.
- `publish.yml` also accepts `workflow_dispatch`, against a tag. A tag push is still the release,
  but a first upload to a project PyPI has never seen needs a trusted publisher registered first,
  and that is a browser action outside this repository — so every gate here can pass and the
  upload still fail on configuration alone. Without a manual trigger the only remedies are
  re-running a failed run or moving a published tag. The version assert reads `GITHUB_REF_NAME`,
  so dispatching against a branch is refused, which is the intended behaviour.
- The local-path scan, `.github/leakscan.sh`, a script rather than a command for each caller to
  copy, so that CI and the pre-release checklist in `CONTRIBUTING.md` cannot drift apart — two
  copies of a pattern this fiddly will not stay in step. The pattern does not match its own source
  line, which is what lets it scan `.github/` at all. It reaches capitalised home directories, the
  JSON-escaped spelling of a path separator, and remote-host forms that name no command, and it
  scans `NOTICE` and `LICENSE` — while not firing on a PEP 508 direct reference, a registry
  digest, or an email address followed by a colon.
- `CONTRIBUTING.md`. The licence is `MIT AND CC-BY-4.0`, not a plain OSI licence, so the terms a
  contribution is offered under are not inferable from the file layout — and a data contribution
  is not the same kind of thing as a code one here. It also carries the pre-release checks that
  only hold for the tree a release is built from rather than for the commit.

### Corrected before publication

A full review of every published surface — numeric recomputation against the shipped store,
cross-document consistency, standalone reading, three audiences, and release mechanics with
mutation testing — ran against this tree before the tag was moved. No energy value changed. What
it found:

- **The install instructions named PyPI, where this package does not exist.** The Quickstart's
  first command was `pip install skempi-foldx==0.2.0`, and five further passages asserted a PyPI
  presence. Every release installs from its git tag; the documents now say so, and a test refuses
  a `skempi-foldx==` pin in the README.
- **"Neither withdrawn tag resolves" was half wrong, in the dangerous direction.** `v0.2.1` fails
  outright, but `v0.2.0` resolves to this tree — silently, with different content than the
  withdrawn one. Both cases are now stated separately.
- **The version-detection one-liner crashed on the install it was written for.**
  `direct_url.json` exists only for a git or local install, so `json.loads` was handed `None`. It
  now reports rather than raising, and the release's commit can be obtained from the remote, since
  a release cannot state its own hash.
- **Counts that summed the two arms were labelled as distinct entities.** 477 is `(arm,
  definition)` pairs, not interface definitions — there are 347, because 130 appear in both arms.
  The same applied to 476 paired tiles, to 151 definitions holding one record (97 distinct), and
  to "4343 single-point mutations", which counts entries: the shipped 4340 are 4334 distinct
  mutations. Each is now labelled, with the distinct figure alongside.
- **`docs/DETERMINISM.md` gave a both-arms count in a single-point-only paragraph** (194 where the
  arm holds 58), and the changelog called six multi-point *mutants* six *records*; they are twelve.
- **The repair-count claim this release retracted survived in two code files.**
  `config.py` and `experiments/repair_ablation.py` still called one pass "FoldX's own documented
  recommendation" after the prose had re-attributed that reading to Usmanova et al.
  `docs/REFERENCES.md` stated their finding without the qualifier that it concerns *folding* ΔΔG.
- **The README under-stated the repair-count evidence its own docs carry**, quoting only a
  13-complex pilot and dismissing it, where the completed four-round sweep over 344 complexes
  finds no round after which the store stops changing.
- **Platform provenance was in the wrong document.** The architecture split between the two inputs
  is now in `docs/STORE.md`'s provenance section, alongside the repair count.
- **A statistic no shipped script computes has been removed** rather than restated: the χ-angle
  Mann-Whitney p-value could not be reproduced from the tree. The counts and medians it summarised
  remain, and are regenerable.
- **`experiments/repair_ablation.py` described the 0.1.0 store shape** and a scope size that the
  rebuild changed (98 complexes, now 99).
- **Both `MALFORMED_ROWS` sat unmarked in the corrections table**, which a 0.1.0 consumer is told
  to check against; they are now footnoted as wrong in *both* releases.
- **`CONTRIBUTING.md` mis-stated the release mechanics** — that the recommended push forces (it
  does not), that pushing a tag is the whole release (the upload can still fail on configuration),
  and that a *version* rather than a *project* PyPI has never seen is the case needing care.

### Hardened

- **`publish.yml` refused a branch by name, not by kind.** The version assert compared
  `GITHUB_REF_NAME` against `pyproject.toml`, so a `workflow_dispatch` against a branch named
  `0.2.0` — a release-prep branch named after its version is an ordinary thing to make — passed
  it and would have published an untagged tree. The job now refuses any ref that is not a tag,
  and does so as a failure rather than a skipped step.
- `actions/checkout` no longer persists the job token into the one job that executes the tagged
  tree.
- **Five invariants had no test that noticed them being broken**, found by mutating each and
  re-running the suite: the refusal to infer an author form from a name SKEMPI does not list
  (a third route to the fabricated-alias class this release exists to close); building a role name
  from the record's own interface definition rather than another pairing of its code; first-store-
  wins on a conflict; the shipped store containing no excluded complex; and the energy values
  themselves, now pinned by digest.

### Changed

- Python classifiers name **3.9 through 3.14** individually, and CI builds all six. 0.1.0's
  metadata said only `Programming Language :: Python :: 3`; there was never a resolution reason to
  advertise less than the package supports — it has no dependencies — only an untested one. The
  3.9 leg is pinned to `ubuntu-24.04` pre-emptively: `actions/python-versions` publishes no
  `linux-26.04` build for 3.9 at any patch level, 3.9 having reached end of life before that image
  existed. `ubuntu-latest` still resolves to 24.04 today, so that pin is
  currently a no-op; 3.10 reaches end of life in October 2026 and its 26.04 coverage begins only
  at 3.10.20, so a second entry is expected on the same grounds.
- `MANIFEST.in` names files by extension instead of grafting `docs/`, `tests/` and `experiments/`.
  A graft walks the filesystem rather than git, so it shipped whatever running the pipeline left
  behind — including `rotabase.txt`, which is FoldX's own licensed support data and not this
  project's to redistribute, and a generated CSV withheld because it carries per-row SKEMPI
  affinity data. The protection is defeated by a stale `egg-info/SOURCES.txt`, so `CONTRIBUTING.md`
  gives the clean-build check that catches it.
- **This release was not distributed on PyPI**, and installs from its git tag. `publish.yml`
  shipped ready — it uploads on a version tag via Trusted Publishing, gated on the ref being a
  tag, on the tag matching `pyproject.toml`, on the suite, and on the shipped store being present
  in the built wheel — but no upload was made from it. A PyPI upload cannot be undone, and here a
  version names a measurement rather than a code state, so the first upload was made its own
  release rather than a side effect of this one; that is 0.2.2, which carries these same energies.
- `Changelog` and `Issues` join the project URLs.
- The `dev` extra names `build`, which the pre-release checklist in `CONTRIBUTING.md` invokes.
- Entries here are newest-first, per Keep a Changelog.
- README code blocks are kept narrow, so they are not clipped where the rendered description
  column is narrower than GitHub's.

## 0.1.0 — 2026-08-01

First public release. **Energy values: initial.**

- 6003 records over 344 complexes: 4238 single-point mutations across 322, and 1765 multi-point
  variants across 152. Twelve `AnalyseComplex` terms per record, as mutant minus wild-type.
- `FoldxLookup` resolves the naming conventions a mutation appears under, which is the
  difference between 65.4% and 99.6% coverage on a split whose labels are role-chain.
- Every value computed with a single `RepairPDB` pass, at FoldX 5.1.

Known limitations, each documented where it applies:

- **Three PDB codes carry two SKEMPI interface definitions each** (`2C5D`, `3SE3`, `3SE4`). The
  parser keys on the code, so 31 records are scored against an interface excluding the chain they
  mutate. `interface_suspect()` identifies them; `experiments/recompute_alternate_interfaces.py`
  recomputes them. Fixing the parser moves the coverage denominators, so it is held for 0.2.0
  where the store is rebuilt to match.
- **2494 of 4238 single-point values came from a per-benchmark campaign** rather than from one
  built around each complex's full mutation list. For 27 of those complexes, 287 records, that
  list is a strict subset of the complex's full one, and a value is a property of the list it was
  computed in.
- **`1KBH` is excluded**: `RepairPDB` does not terminate on it under FoldX 5.1.
