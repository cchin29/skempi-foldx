# References

Every DOI below was checked: resolved against the CrossRef API, or — for preprints, which
CrossRef indexes differently — confirmed to resolve at `doi.org`. None was inferred from a
citation pattern. (One near-miss worth recording: a plausible-looking DOI for the Usmanova paper,
differing by one digit, belongs to a different paper on the same topic in the same journal and
year.)

**This list is not exhaustive of the works named under `docs/`.** It covers the data sources, the
tools, the comparator methods this protocol is measured against, and the method background. Seven
further works are named in the docs without a citation because no DOI for them could be confirmed;
each is listed below at the point where it would otherwise appear, so the gap is visible rather
than silent.

## The data and the tool

**SKEMPI 2.0** — the source of every mutation, chain grouping and cleaned mutation string here.
Released under CC BY 4.0; see [NOTICE](../NOTICE) for what is and is not redistributed.

> Jankauskaitė, J., Jiménez-García, B., Dapkūnas, J., Fernández-Recio, J. & Moal, I. H.
> SKEMPI 2.0: an updated benchmark of changes in protein–protein binding energy, kinetics and
> thermodynamics upon mutation. *Bioinformatics* **35**(3):462–469 (2019).
> [10.1093/bioinformatics/bty635](https://doi.org/10.1093/bioinformatics/bty635) ·
> <https://life.bsc.es/pid/skempi2>

**FoldX** — produced every energy in this repository. Not redistributed; free for academic and
non-profit institutions, paid for commercial use.

> Delgado, J., Radusky, L. G., Cianferoni, D. & Serrano, L. FoldX 5.0: working with RNA, small
> molecules and a new graphical interface. *Bioinformatics* **35**(20):4168–4169 (2019).
> [10.1093/bioinformatics/btz184](https://doi.org/10.1093/bioinformatics/btz184) ·
> <https://foldxsuite.crg.eu/>

**FoldX, revised force field** — the FoldX authors' own 2025 revision, the source of the
`--numberOfRuns=5`-with-the-median recommendation weighed in [DETERMINISM.md](DETERMINISM.md).
Its reported gain was measured on folding benchmarks, not on SKEMPI binding data.

> Delgado, J., Reche, R., Cianferoni, D., Orlando, G., van der Kant, R., Rousseau, F.,
> Schymkowitz, J. & Serrano, L. FoldX force field revisited, an improved version.
> *Bioinformatics* **41**(2) (2025).
> [10.1093/bioinformatics/btaf064](https://doi.org/10.1093/bioinformatics/btaf064)

**pyFoldX** — Python bindings for FoldX. A different thing from this package, and the right
answer when bindings are what is needed.

> Radusky, L. G. & Serrano, L. pyFoldX: enabling biomolecular analysis and engineering along
> structural ensembles. *Bioinformatics* **38**(8):2353–2355 (2022).
> [10.1093/bioinformatics/btac072](https://doi.org/10.1093/bioinformatics/btac072)

## Methods this protocol is measured against

Characterised in [PROTOCOL.md](PROTOCOL.md).

**CATH-ddG** — the only comparator with a fully specified FoldX protocol, and the source of the
CATH-superfamily partition used as an evaluation split.

> Yu, G., Bi, X., Ma, T., Li, Y. & Wang, J. CATH-ddG: towards robust mutation effect prediction
> on protein–protein interactions out of CATH homologous superfamily. *Bioinformatics*
> **41**(Supplement 1):i362–i372 (2025).
> [10.1093/bioinformatics/btaf228](https://doi.org/10.1093/bioinformatics/btaf228)

**USP-ddG** — consumes FoldX empirical energy terms as one of its channels. **Shares its first,
second and last authors with CATH-ddG above**, which matters for any inference drawn from the two
reporting the same value.

> Yu, G., Bi, X., Zhao, Q. & Wang, J. USP-ddG: a unified structural paradigm with data efficacy
> and mixture-of-experts for predicting mutational effects on protein–protein interactions.
> *Bioinformatics* **42**(Supplement 1) (2026).
> [10.1093/bioinformatics/btag249](https://doi.org/10.1093/bioinformatics/btag249) ·
> preprint [10.1101/2025.11.09.687124](https://doi.org/10.1101/2025.11.09.687124)

**RDE-Network** — the source of the `block_list` that this package's `1KBH` exclusion follows.

> Luo, S., Su, Y., Wu, Z., Su, C., Peng, J. & Ma, J. Rotamer Density Estimator is an Unsupervised
> Learner of the Effect of Mutations on Protein–Protein Interaction. *ICLR* (2023);
> *bioRxiv* 2023.02.28.530137.
> [10.1101/2023.02.28.530137](https://doi.org/10.1101/2023.02.28.530137) ·
> <https://github.com/luost26/RDE-PPI>

**GearBind** — builds mutant structures with FoldX for training and reports a FoldX baseline,
without stating a repair count.

> Cai, H., Zhang, Z., Wang, M., Zhong, B., Li, Q., Zhong, Y., Wu, Y., Ying, T. *et al.*
> Pretrainable geometric graph neural network for antibody affinity maturation.
> *Nature Communications* **15**:7785 (2024).
> [10.1038/s41467-024-51563-8](https://doi.org/10.1038/s41467-024-51563-8)

**Named in [PROTOCOL.md](PROTOCOL.md) without a citation.** ProtBFF, Prompt-DDG, BA-DDG, GraphPPI,
PPIformer and MINT appear in that file's survey table. No DOI for any of them could be confirmed
against the CrossRef API, so none is given rather than one being inferred — including for
GraphPPI, whose row carries the strongest characterisation in the table. PROTOCOL.md states this
at the point of use.

**The benchmark subsets.** `S1102`, `S1131`, `S2003` and `S4169` are named throughout these docs
as evaluation sets and as campaign labels. They are conventional partitions of SKEMPI reused
across this literature rather than works in their own right; no citation is given because the
usage here is to SKEMPI 2.0 above, whose rows they select. A consumer taking subset membership
should take it from the benchmark definition it means to follow, not from this package.

## Method background

**Repair-count practice.** The only controlled study of iterating `RepairPDB` — its finding that
iteration does not improve ΔΔG or reduce bias is the reason a single repair is defensible, and
the reason the repair-count question here is measured rather than assumed. See
[DETERMINISM.md](DETERMINISM.md).

> Usmanova, D. R., Bogatyreva, N. S., Ariño Bernad, J. *et al.* Self-consistency test reveals
> systematic bias in programs for prediction change of stability upon mutation. *Bioinformatics*
> **34**(21):3653–3658 (2018).
> [10.1093/bioinformatics/bty340](https://doi.org/10.1093/bioinformatics/bty340)

**Structural sensitivity of FoldX.** The 0.61 kcal/mol spread across different PDB structures of
the same protein, quoted in [DETERMINISM.md](DETERMINISM.md) as one of the published magnitudes
that rule out a stochastic explanation for the 11 kcal/mol observed here.

> Caldararu, O., Blundell, T. L. & Kepp, K. P. A base measure of precision for protein stability
> predictors: structural sensitivity. *BMC Bioinformatics* **22**:88 (2021).
> [10.1186/s12859-021-04030-w](https://doi.org/10.1186/s12859-021-04030-w)

**Meli et al. 2024** — quoted in [DETERMINISM.md](DETERMINISM.md) as the best protocol-reporting
template found in this literature. **No citation is given here**: the paper could not be
identified against the CrossRef API from the *IJMS* 2024 attribution recorded with the quotation,
and inferring a DOI would be worse than leaving the gap visible. The quoted sentence should be
treated as unsourced until the paper is located.

**CD-HIT** — used to build the ≤60% sequence-identity clustered evaluation splits.

> Fu, L., Niu, B., Zhu, Z., Wu, S. & Li, W. CD-HIT: accelerated for clustering the
> next-generation sequencing data. *Bioinformatics* **28**(23):3150–3152 (2012).
> [10.1093/bioinformatics/bts565](https://doi.org/10.1093/bioinformatics/bts565)

**CATH** — the structural classification underlying the superfamily-level evaluation split.

> Sillitoe, I., Bordin, N., Dawson, N. *et al.* CATH: increased structural coverage of functional
> space. *Nucleic Acids Research* **49**(D1):D266–D273 (2021).
> [10.1093/nar/gkaa1079](https://doi.org/10.1093/nar/gkaa1079) · <https://www.cathdb.info/>

## Origin of this package

**MuLAN** — this code began as a FoldX score channel inside a research fork of MuLAN. None of
that project's code is present here, which is what allows this one to be MIT licensed; see
[NOTICE](../NOTICE).

> Lombardi, G. & Carbone, A. MuLAN: Mutation-driven Light Attention Networks for investigating
> protein-protein interactions from sequences. *bioRxiv* 2024.08.24.609515 (2024).
> [10.1101/2024.08.24.609515](https://doi.org/10.1101/2024.08.24.609515)
