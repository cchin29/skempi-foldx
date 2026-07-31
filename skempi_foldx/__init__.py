"""FoldX binding ΔΔG for SKEMPI 2.0 — compute it, or just use the results.

    from skempi_foldx import load_bundled_store
    store = load_bundled_store()               # {pdb: {mutation: {12 energy terms}}}

Two things ship here, and the second is the one most people want:

1. **The pipeline.** ``RepairPDB -> BuildModel -> AnalyseComplex`` over SKEMPI's chain groups,
   yielding ΔΔG_int = IE(mutant) - IE(wild-type) decomposed into FoldX's 12 energy terms.
   Resumable per complex, parallel across complexes, and honest about what it could not compute.

2. **The results.** 322 complexes / 4238 single-point mutations and 152 / 1765 multi-point
   variants -- 97.7% and 95.5% of what SKEMPI defines, and 97.8%/100% of what is reachable
   once the one unrepairable structure is set aside. Reproducing these needs a
   FoldX licence and CPU-weeks; reading them needs neither, and this package has no dependencies
   beyond the standard library.

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
    Exclusion,
    filter_complexes,
    filter_worklist,
    is_excluded,
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
    worklist_multi_point,
    worklist_single_point,
)
from .skempi import Mutation, SkempiComplex, load_skempi, map_role_to_author
from .store import (
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

__version__ = "0.1.0"

__all__ = [
    "N_TERMS", "SCALAR_TERM", "TERMS", "term_vector",
    "FoldxConfig", "default_config",
    "FoldxLookup", "load_chain_mapping", "to_role_form",
    "INTRACTABLE", "Exclusion", "is_excluded", "filter_complexes", "filter_worklist",
    "MODE_AUTHOR", "MODE_ROLE", "MODE_VARIANT",
    "ComplexResult", "process_complex", "run_campaign", "format_individual_list",
    "worklist_from_table", "worklist_single_point", "worklist_multi_point",
    "exclude_already_computed",
    "Mutation", "SkempiComplex", "load_skempi", "map_role_to_author",
    "ConsolidationReport", "audit", "consolidate", "coverage", "source_label",
    "DATA_DIR", "SINGLE_POINT", "MULTI_POINT", "bundled_path", "load_bundled_store",
    "load_complex", "load_complex_kind", "load_store", "store_kind",
    "reindex_by_skempi_id",
    "infer_kind", "write_store",
]
