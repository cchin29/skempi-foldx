"""The FoldX energy-term contract.

These twelve terms, in this order, are what a result record contains: ``Interaction Energy``,
the total, followed by eleven of its components. ``AnalyseComplex`` reports about twenty-one
energy columns, so the eleven do not sum to the total and the total is not redundant with them.
The order is load-bearing for any consumer that flattens a record into a feature vector. Written under one order and read under another, it yields permuted features and
a plausible-looking, wrong result.

This module is the single definition, and consumers should import it rather than keep a copy: two
copies of an ordering can agree for a long time with nothing enforcing that they continue to.

Consumer-side constants — the clip bound, the fill value for an uncovered row, the column
layouts of a standardized split file — deliberately do **not** live here. Those describe one
particular model's feature encoding, not FoldX. Keeping them out is what lets this package stay
dependency-free and useful to a consumer that encodes its features differently.
"""

from __future__ import annotations

from typing import List, Sequence

#: The AnalyseComplex terms captured as the decomposed feature vector, as ``mutant - wildtype``.
TERMS: List[str] = [
    "Interaction Energy",
    "Backbone Hbond",
    "Sidechain Hbond",
    "Van der Waals",
    "Electrostatics",
    "Solvation Polar",
    "Solvation Hydrophobic",
    "Van der Waals clashes",
    "entropy sidechain",
    "entropy mainchain",
    "torsional clash",
    "backbone clash",
]

#: The summary term. A consumer wanting one number rather than twelve wants this one.
SCALAR_TERM = "Interaction Energy"

N_TERMS = len(TERMS)

assert TERMS[0] == SCALAR_TERM, "the scalar term must be first, so it is column 0 of any fit"


def term_vector(entry: dict, terms: Sequence[str] = TERMS) -> List[float]:
    """Pull the term vector out of a per-mutation FoldX record, in canonical order.

    Raises rather than silently zero-filling a missing term: a record that lacks a term is a
    parsing failure upstream, and imputing it here would hide that.
    """
    missing = [t for t in terms if t not in entry]
    if missing:
        raise KeyError(f"FoldX record is missing term(s): {missing}")
    return [float(entry[t]) for t in terms]
