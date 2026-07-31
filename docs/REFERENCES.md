# References

Every DOI below was checked: resolved against the CrossRef API, or — for preprints, which
CrossRef indexes differently — confirmed to resolve at `doi.org`. None was inferred from a
citation pattern. (One near-miss worth recording: a plausible-looking DOI for the Usmanova paper,
differing by one digit, belongs to a different paper on the same topic in the same journal and
year.)

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

**USP-ddG** — consumes FoldX empirical energy terms as one of its channels.

> USP-ddG: a unified structural paradigm with data efficacy and mixture-of-experts for predicting
> mutational effects on protein–protein interactions. *bioRxiv* 2025.11.09.687124 (2025).
> [10.1101/2025.11.09.687124](https://doi.org/10.1101/2025.11.09.687124)

**RDE-Network** — the source of the `block_list` that this package's `1KBH` exclusion follows.

> Luo, S., Su, Y., Wu, Z., Su, C., Peng, J. & Ma, J. Rotamer Density Estimator is an Unsupervised
> Learner of the Effect of Mutations on Protein–Protein Interaction. *ICLR* (2023);
> *bioRxiv* 2023.02.28.530137.
> [10.1101/2023.02.28.530137](https://doi.org/10.1101/2023.02.28.530137) ·
> <https://github.com/luost26/RDE-PPI>

## Method background

**Repair-count practice.** The only controlled study of iterating `RepairPDB` — its finding that
iteration does not improve ΔΔG or reduce bias is the reason a single repair is defensible, and
the reason the repair-count question here is measured rather than assumed. See
[DETERMINISM.md](DETERMINISM.md).

> Usmanova, D. R., Bogatyreva, N. S., Ariño Bernad, J. *et al.* Self-consistency test reveals
> systematic bias in programs for prediction change of stability upon mutation. *Bioinformatics*
> **34**(21):3653–3658 (2018).
> [10.1093/bioinformatics/bty340](https://doi.org/10.1093/bioinformatics/bty340)

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
