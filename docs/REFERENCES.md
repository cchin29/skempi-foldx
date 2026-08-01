# References

Every DOI below was checked: resolved against the CrossRef API, or — for preprints, which
CrossRef indexes differently — confirmed to resolve at `doi.org`. None was inferred from a
citation pattern. (One near-miss worth recording: a plausible-looking DOI for the Usmanova paper,
differing by one digit, belongs to a different paper on the same topic in the same journal and
year.)

This covers every work named under `docs/`: the data sources, the tools, the comparator methods
this protocol is measured against, and the method background.

## The data and the tool

**SKEMPI 2.0** — the source of every mutation, chain grouping and cleaned mutation string here.
Released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); see
[NOTICE](../NOTICE) for what is and is not redistributed.

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

**CATH-ddG** — the comparator whose FoldX protocol is specified in most detail, and the source of the
CATH-superfamily partition used as an evaluation split.

> Yu, G., Bi, X., Ma, T., Li, Y. & Wang, J. CATH-ddG: towards robust mutation effect prediction
> on protein–protein interactions out of CATH homologous superfamily. *Bioinformatics*
> **41**(Supplement 1):i362–i372 (2025).
> [10.1093/bioinformatics/btaf228](https://doi.org/10.1093/bioinformatics/btaf228)

**USP-ddG** — consumes FoldX empirical energy terms as one of its channels, and reports the same
FoldX baseline row as CATH-ddG above.

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

**GearBind** — builds mutant structures with FoldX for training and reports a FoldX baseline; the
command sequence is given in prose, the repair count is not.

> Cai, H., Zhang, Z., Wang, M., Zhong, B., Li, Q., Zhong, Y., Wu, Y., Ying, T. *et al.*
> Pretrainable geometric graph neural network for antibody affinity maturation.
> *Nature Communications* **15**:7785 (2024).
> [10.1038/s41467-024-51563-8](https://doi.org/10.1038/s41467-024-51563-8)

**Also characterised in [PROTOCOL.md](PROTOCOL.md).** Six further methods appear in that file's
survey table. Most are conference papers whose venues do not register DOIs, so an arXiv identifier
is given where that is the citable form.

> Feldman, J., Maechler, A., Wang, D. & Shakhnovich, E. I. A General Framework for Injecting
> Biophysical Priors into Protein Embeddings. *bioRxiv* 2025.12.23.696257 (2025). — **ProtBFF**
> [10.64898/2025.12.23.696257](https://doi.org/10.64898/2025.12.23.696257)

> Wu, L., Tian, Y., Lin, H., Huang, Y., Li, S., Chawla, N. V. & Li, S. Z. Learning to Predict
> Mutational Effects of Protein-Protein Interactions by Microenvironment-aware Hierarchical Prompt
> Learning. *ICML* (2024), PMLR 235:53847–53859. — **Prompt-DDG**
> [arXiv:2405.10348](https://arxiv.org/abs/2405.10348)

> Jiao, X., Mao, W., Jin, W., Yang, P., Chen, H. & Shen, C. Boltzmann-Aligned Inverse Folding
> Model as a Predictor of Mutational Effects on Protein–Protein Interactions. *ICLR* (2025). —
> **BA-DDG** [arXiv:2410.09543](https://arxiv.org/abs/2410.09543)

> Liu, X., Luo, Y., Song, S. & Peng, J. Pre-training of Graph Neural Network for Modeling Effects
> of Mutations on Protein-Protein Binding Affinity (2020). — **GraphPPI**
> [arXiv:2008.12473](https://arxiv.org/abs/2008.12473)
>
> The published form of this preprint appears to be Liu, Luo, Li, Song & Peng, *PLoS Comput Biol*
> **17**(8):e1009284 (2021),
> [10.1371/journal.pcbi.1009284](https://doi.org/10.1371/journal.pcbi.1009284), where the method is
> renamed **GeoPPI**: arXiv's record for the preprint gives that DOI as its related publication,
> and the two share benchmarks, head architecture and case study. Neither document states the link,
> so both identifiers are given — the arXiv one for the name *GraphPPI*, which appears only in the
> preprint, and the DOI for the peer-reviewed form. The survey row in
> [PROTOCOL.md](PROTOCOL.md) is read from the preprint; the published version reports its FoldX
> baseline as its own run rather than a quoted one.

> Bushuiev, A., Bushuiev, R., Kouba, P. *et al.* Learning to Design Protein-Protein Interactions
> with Enhanced Generalization. *ICLR* (2024). — **PPIformer**
> [arXiv:2310.18515](https://arxiv.org/abs/2310.18515)

> Ullanat, V., Jing, B., Sledzieski, S. & Berger, B. Learning the language of protein-protein
> interactions. *Nature Communications* **17**:1199 (2026). — **MINT**
> [10.1038/s41467-025-67971-3](https://doi.org/10.1038/s41467-025-67971-3)

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

**Protocol reporting.** Quoted in [DETERMINISM.md](DETERMINISM.md) as the best protocol-reporting
template found in this literature: it states the repair count, the run count and the environment
parameters in two sentences, which no comparator paper does.

> Vincenzi, M., Mercurio, F. A., La Manna, S., Palumbo, R., Pirone, L., Marasco, D., Pedone, E. M.
> & Leone, M. Exploring a Potential Optimization Route for Peptide Ligands of the Sam Domain from
> the Lipid Phosphatase Ship2. *International Journal of Molecular Sciences* **25**(19):10616
> (2024). [10.3390/ijms251910616](https://doi.org/10.3390/ijms251910616)

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
that project's code is present here; see
[NOTICE](../NOTICE).

> Lombardi, G. & Carbone, A. MuLAN: Mutation-driven Light Attention Networks for investigating
> protein-protein interactions from sequences. *bioRxiv* 2024.08.24.609515 (2024).
> [10.1101/2024.08.24.609515](https://doi.org/10.1101/2024.08.24.609515)
