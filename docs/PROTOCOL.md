# FoldX protocols in the published ΔΔG literature

For matching another group's FoldX numbers, or justifying these settings in a paper, this is the
survey. The short answer is that the field has no single convention, and that FoldX settings are
usually not recoverable from a published methods section.

That is a reasonable state of affairs rather than an oversight. In almost every paper here FoldX
is a baseline or a structure-generation step, not the contribution, and methods sections are
written to the length the contribution warrants. The consequence is still worth knowing: a FoldX
column cannot be reproduced from most of these papers, so this one writes its own down. FoldX 5.1,
`RepairPDB` once per structure, then:

```
--command=BuildModel     --pdb=<code>_Repair.pdb --mutant-file=individual_list.txt --numberOfRuns=1
--command=AnalyseComplex --pdb-list=pdblist.txt  --analyseComplexChains=<group1>,<group2>
```

No `--pH`, `--temperature`, `--ionStrength`, `--water` or `--vdwDesign` is passed on any
invocation, so the engine's defaults apply throughout. The repair count is evidenced in
[STORE.md](STORE.md) and the run count reasoned through in [DETERMINISM.md](DETERMINISM.md).

> **Sourcing.** The characterisations below were made by reading each paper's methods section and
> supplement. The rows are not individually keyed to table numbers — treat this as a working
> survey to check against the primary sources, not as a citable review. Any specific row being
> relied on should be verified against the paper itself.
>
> Every method in the table is cited in [REFERENCES.md](REFERENCES.md). Several are conference
> papers whose venues register no DOI, so an arXiv identifier is the citable form. Where a row says
> something is not stated, that means it was not found in what was read — a main text, and a
> supplement where one was available — which is not the same as it being absent from the work.

---

There is **no community convention** to follow. Of the ten works surveyed below, four name the
commands — CATH-ddG, ProtBFF, RDE-Network and GearBind. Two state a repair count: ProtBFF in its
methods text ("twice"), CATH-ddG in its supplement ("only once for each protein complex
structure"). CATH-ddG is the only one to print the command lines; the other three describe the
sequence in prose, and RDE-Network and ProtBFF name two commands rather than three.
`numberOfRuns`, `--vdwDesign`, `--pH`, `--temperature`, `--water` and `--ionStrength` appear in
none of them, and none discusses FoldX run-to-run behaviour. Several report a FoldX baseline
without an accompanying protocol.

| Paper | FoldX role | RepairPDB | `numberOfRuns` | FoldX baseline |
|---|---|---|---|---|
| **CATH-ddG** | mutant structures + AnalyseComplex terms → MLP; baseline | **×1** | not stated (output `_1` implies default) | own run |
| **USP-ddG** | `ΔG_FoldX` term + its own loss branch; baseline | not stated | not stated | its FoldX row matches CATH-ddG's to the reported digits |
| **ProtBFF** | mutant structures → dihedral/lDDT features | **×2** | not stated | own run |
| **RDE-Network** | baseline only | sequence stated; count not given | not stated | own run — the origin of the widely-quoted 0.3789 |
| **GearBind** | mutant structures for training; baseline | sequence stated; count not given | not stated | own run |
| Prompt-DDG | baseline only | — | — | quoted from RDE-Network, stated |
| BA-DDG | baseline only | — | — | matches RDE-Network's row; provenance not stated |
| GraphPPI | mutant structures | not stated | not stated | quoted |
| PPIformer | not used as a ΔΔG baseline | — | — | — |
| MINT | FoldX-derived ΔΔG as task labels, via a cited dataset; not a baseline | — | — | — |

Two things follow that are worth knowing before quoting any FoldX number.

**Fewer of these FoldX baselines are original runs than the count of them suggests.** Eight of the
ten works report a FoldX baseline. Five state an original run — CATH-ddG, ProtBFF, RDE-Network,
GearBind, and GraphPPI in its published form, though the survey row here is read from its preprint,
which does not. Three carry a value from elsewhere: Prompt-DDG, which names the row it took from
RDE-Network, and USP-ddG and BA-DDG, which do not say where theirs came from.
In the CATH-superfamily table, where a new
method's result would sit beside it, the FoldX row carries the same value in CATH-ddG and in
USP-ddG, and only CATH-ddG says how that value was produced. Agreeing with that row is therefore
agreeing with one measurement rather than with a consensus of several.

**This protocol matches the most fully specified one, step for step.** CATH-ddG's supplementary
section 1.4 gives `RepairPDB` (once) → `BuildModel` → `AnalyseComplex` on both structures, FoldX
5.0, defaults otherwise, with the command lines printed — the detail is in the supplement rather
than the article PDF. RDE-Network and GearBind give the same sequence in prose; printing the
commands and the version is what makes CATH-ddG reproducible from the paper alone.

That sequence is this pipeline, with two differences. CATH-ddG ran FoldX 5.0 and these values were
produced with 5.1, so the commands match step for step while the force field does not. And
CATH-ddG runs `Optimize` on the mutant structure, which this pipeline does not — that step feeds
their *structural* input, while only the `AnalyseComplex` energies are taken here, so it does not
enter the numbers being compared.
