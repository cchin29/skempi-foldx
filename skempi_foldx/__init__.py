"""FoldX binding ΔΔG for SKEMPI 2.0 — compute it, or just use the results.

    from skempi_foldx import load_bundled_store
    store = load_bundled_store()      # {'1BRS_A_D': {'DA52A': {12 energy terms}}, ...}

Two things ship here, and the second is the one most people want:

1. **The pipeline.** ``RepairPDB -> BuildModel -> AnalyseComplex`` over SKEMPI's chain groups,
   yielding ΔΔG_int = IE(mutant) - IE(wild-type), reported as the total and eleven of its
   component terms.
   Resumable per complex, parallel across complexes, and recording every complex it could not
   compute rather than omitting it.

2. **The results.** 4340 single-point mutations over 323 interface definitions, and 1767
   multi-point variants over 154 — 4340 of the 4343 and 1767 of the 1850 SKEMPI defines, which is
   all of both once the one unrepairable structure is set aside. Reproducing these needs a FoldX
   licence and CPU-weeks; reading them needs neither, and this package has no dependencies beyond
   the standard library.

   Two units share the word "record" in this literature, so this package fixes both. The **store
   unit** — what a results file holds, and what the store is keyed by — is one SKEMPI
   **interface definition**, ``<pdb>_<group1>_<group2>``, because three codes are defined under two
   chain pairings each and the pairing is half of what a value means. A **record** is one mutation
   within such a definition: one set of twelve terms. Counts in the documentation are records in
   that second sense — 6107 of them, over 477 definition-arms and 347 distinct definitions.

Read ``docs/DETERMINISM.md`` before running a campaign. The short version: FoldX 5.1 is
deterministic, but a mutation's ΔΔG depends on *every entry preceding it* in
``individual_list.txt`` — so computing a subset of a complex's mutations gives different numbers
than computing all of them. That is not a caveat, it is the single most important operational
fact about this pipeline.
"""

from .config import FoldxConfig, default_config
from .lookup import FoldxLookup, load_chain_mapping, to_role_form
from .exclusions import (
    INTRACTABLE,
    code_of,
    Exclusion,
    filter_complexes,
    announce_malformed,
    filter_worklist,
    is_excluded,
    is_malformed,
    MALFORMED_ROWS,
    MalformedRow,
)
from .run import (
    MODE_AUTHOR,
    MODE_ROLE,
    MODE_VARIANT,
    ComplexResult,
    exclude_already_computed,
    format_individual_list,
    process_complex,
    run_campaign,
    worklist_from_table,
    pool_by_code,
    worklist_multi_point,
    worklist_single_point,
)
from .skempi import (Mutation, SkempiComplex, by_pdb, load_skempi,
                     map_role_to_author)
from .store import (
    check_code,
    check_identifier,
    DATA_DIR,
    MULTI_POINT,
    SINGLE_POINT,
    ConsolidationReport,
    bundled_path,
    load_bundled_store,
    audit,
    consolidate,
    coverage,
    infer_kind,
    load_complex,
    load_complex_kind,
    load_store,
    reindex_by_skempi_id,
    source_label,
    store_kind,
    write_store,
)
from .terms import N_TERMS, SCALAR_TERM, TERMS, term_vector

__version__ = "0.2.0"

__all__ = [
    "N_TERMS", "SCALAR_TERM", "TERMS", "term_vector",
    "FoldxConfig", "default_config",
    "FoldxLookup", "load_chain_mapping", "to_role_form",
    "INTRACTABLE", "Exclusion", "is_excluded", "filter_complexes", "filter_worklist",
    "code_of",
    "MALFORMED_ROWS", "MalformedRow", "is_malformed", "announce_malformed",
    "MODE_AUTHOR", "MODE_ROLE", "MODE_VARIANT",
    "ComplexResult", "process_complex", "run_campaign", "format_individual_list",
    "worklist_from_table", "worklist_single_point", "worklist_multi_point",
    "pool_by_code",
    "exclude_already_computed",
    "Mutation", "SkempiComplex", "by_pdb", "load_skempi", "map_role_to_author",
    "ConsolidationReport", "audit", "consolidate", "coverage", "source_label",
    "DATA_DIR", "SINGLE_POINT", "MULTI_POINT", "bundled_path", "load_bundled_store",
    "load_complex", "load_complex_kind", "load_store", "store_kind",
    "reindex_by_skempi_id",
    "infer_kind", "write_store", "check_code", "check_identifier",
]
