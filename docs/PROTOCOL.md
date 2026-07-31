# How this protocol compares with the published ΔΔG literature

If you are trying to match someone else's FoldX numbers, or to justify these settings in a paper,
this is the survey. The short answer is that the field has no single convention, and several
widely-cited papers do not specify theirs at all.

> **Sourcing.** The characterisations below were made by reading each paper's methods section and
> supplement; this document does not carry per-claim citations with DOIs and table numbers. Treat
> it as a working survey to check against the primary sources, not as a citable review. If you
> are relying on a specific row, verify it against the paper itself.

---


There is **no community convention** to follow. Reading the comparator papers directly: of eleven
in the reference library, only one states how many times it ran `RepairPDB`, **none** states
`numberOfRuns`, none states `--vdwDesign`, `--pH`, `--temperature`, `--water` or `--ionStrength`,
and **none acknowledges FoldX run-to-run behaviour at all**. Several report a FoldX baseline with
no protocol whatsoever.

| Paper | FoldX role | RepairPDB | `numberOfRuns` | FoldX baseline |
|---|---|---|---|---|
| **CATH-ddG** | mutant structures + AnalyseComplex terms → MLP; baseline | **×1** | not stated (output `_1` implies default) | own run |
| **USP-ddG** | `ΔG_FoldX` term + its own loss branch; baseline | not stated | not stated | its FoldX row matches CATH-ddG's to the reported digits; no independent protocol is given |
| **ProtBFF** | mutant structures → dihedral/lDDT features | **×2** | not stated | own run |
| **RDE-Network** | baseline only | not stated | not stated | own run — the origin of the widely-copied 0.3789 |
| **GearBind** | mutant structures for training; baseline | not stated | not stated | own run |
| Prompt-DDG, BA-DDG | baseline only | — | — | copied from RDE-Network |
| GraphPPI | mutant structures | **never mentioned** | not stated | copied |
| PPIformer, MINT | not used as a ΔΔG baseline | — | — | — |

Two things follow that are worth knowing before quoting any FoldX number.

**The FoldX row a new method is compared against is very likely a single computation.** Across
the papers surveyed here, only a handful report an original FoldX run; the rest quote a FoldX
baseline without stating a protocol. In particular the FoldX row in the CATH-superfamily table —
the row a new method's result would be placed next to — appears with the same value in CATH-ddG
and in USP-ddG, and only CATH-ddG specifies how it was produced.

That is an observation about what the papers state, not a claim about attribution practice: the
value may well be cited in text not checked here. The practical point stands either way — the
number is a single measurement, so agreeing with it is weaker evidence than it looks.

**This protocol matches the only fully specified one, step for step.** CATH-ddG's supplement gives
`RepairPDB` (once) → `BuildModel` → `AnalyseComplex` on both structures, FoldX 5.0, defaults
otherwise. That is exactly this pipeline, and the 12-term arm here is the same construction as
their FoldX energy branch: `AnalyseComplex` interaction terms fed to a small MLP. The one step
they take that this pipeline does not is `Optimize` on the mutant structure — but that feeds their
*structural* input, and only energies are consumed here.

USP-ddG publishes **no** FoldX protocol at all — its appendices are not in the released PDF. So
"matching USP-ddG" is not a thing that can be done; matching CATH-ddG and saying so is.

