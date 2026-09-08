"""What must not be handed to FoldX, and what must not be trusted after it.

Two registries, at two granularities, for two different defects:

* :data:`INTRACTABLE` -- whole complexes FoldX cannot compute at all.
* :data:`MALFORMED_ROWS` -- rows that are wrong in ``skempi_v2.csv`` before FoldX sees them, per
  mutation. A newer engine cannot fix these, which is why they ignore the bypass below.

The first is a compute guard; the second is a trust guard on values already computed. There is
no registry of records scored against the wrong chain pairing, because the store is keyed by
SKEMPI identifier: each pairing carries its own values, so the condition such a registry would
name cannot arise.

This is not a curation filter. Curation already drops these complexes from every split — the
whole point is that they are *absent* from the downstream tables, which is exactly why nothing
downstream stops a compute driver from picking them up.

The failure this prevents
-------------------------
``1KBH`` cannot be repaired. RepairPDB on it has been left running past 22.5 h on one machine and
past 10 h on another, finishing on neither. Killing the job works, but killing is a per-round
remedy and the driver is a loop: the next round enumerates its scope, sees 1KBH again, and
launches RepairPDB again. An intractable complex has to be dropped *before* the worklist is
built, not reaped after.

Route by which an excluded complex re-enters scope
--------------------------------------------------
A stale ``1KBH.json`` — 3 mutations, computed before the complex was known to be intractable —
can sit in a source results directory even when the canonical store has no such file. Any scope
enumerated from the *store* rather than from the *curated table* therefore picks 1KBH up, and
pointing a merger at extra source directories to lift coverage is precisely what puts such a
directory in scope. The stale artifact, not the curation, is the entry point. So the guard lives
here, and any step that materialises a unified results directory must refuse to copy an excluded
complex into it.

Re-enabling
-----------
Each entry records the **FoldX version it was observed against**, because "does not terminate" is
a property of a particular engine, not a law. A future FoldX may well repair 1KBH fine. To find
out::

    SKEMPI_FOLDX_ALLOW_INTRACTABLE=1 python experiments/repair_ablation.py run \\
        --complexes 1KBH --iterations 1

The bypass restores the complex everywhere in the compute path at once. If it now converges,
delete the entry and rebuild the consolidated store, which will then carry it. What the bypass
does *not* do is re-admit the complex to curation: :data:`INTRACTABLE` is consulted there as a
plain set, so a consumer's curation keeps blocking it until the entry is actually removed.
That asymmetry is intended: a complex can be tested without silently changing which
rows every split contains.
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Optional

#: Set to ``1`` to run an excluded complex anyway.
ALLOW_ENV = "SKEMPI_FOLDX_ALLOW_INTRACTABLE"


@dataclass(frozen=True)
class Exclusion:
    pdb: str
    reason: str
    evidence: str
    observed_with: str
    """The FoldX version the behaviour was observed against — what to re-test on an upgrade."""

    def __str__(self) -> str:
        return (f"{self.pdb}: {self.reason} [observed with {self.observed_with}] "
                f"({self.evidence})")


INTRACTABLE: Dict[str, Exclusion] = {
    "1KBH": Exclusion(
        pdb="1KBH",
        reason="RepairPDB does not terminate; the complex has no usable FoldX result",
        evidence="observed >22.5 h and >10 h without finishing on two machines; already dropped "
                 "from curation by RDE's block_list, so it appears in no split; the only "
                 "artifact is a stale 3-mutation JSON in a source results directory",
        observed_with="FoldX 5.1",
    ),
}


def allowed() -> bool:
    """Whether the exclusion list is currently bypassed.

    The test is positive, and deliberately so. A guard whose failure mode is a multi-day stall
    must not disarm on an unrecognised value: under a negative test, ``ALLOW=no`` and ``ALLOW=off``
    both read as a bypass and re-arm ``RepairPDB`` on the very complex this module exists to keep
    out of the worklist.
    """
    return os.environ.get(ALLOW_ENV, "").strip().lower() in ("1", "true", "yes", "on")


def code_of(name: str) -> str:
    """The PDB code in a SKEMPI identifier, or the argument unchanged if it is already one.

    ``1KBH_A_B`` and ``1KBH`` name the same structure, and intractability is a property of the
    structure. A registry keyed by code that is queried with an identifier answers ``None`` for
    every entry it holds — which reads exactly like "nothing is excluded".
    """
    return name.split("_", 1)[0]


def is_excluded(pdb: str) -> Optional[Exclusion]:
    """The :class:`Exclusion` for ``pdb``, or ``None``. Always ``None`` when bypassed.

    Accepts a PDB code or a SKEMPI identifier; both resolve to the structure.
    """
    if allowed():
        return None
    return INTRACTABLE.get(code_of(pdb))


def filter_complexes(pdbs: Iterable[str], on_note=print, context: str = "") -> List[str]:
    """Drop excluded ids from a list of complexes, loudly.

    Silence would be the wrong behaviour twice over: a scope that shrinks by one without comment
    is how a driver ends up reporting a coverage number nobody can account for, and a reader who
    never sees the notice has no way to learn the bypass exists.

    So it is announced twice — through ``on_note`` (the driver's own progress stream, where a
    human watching a run will see it) *and* as a :class:`RuntimeWarning` (which survives a driver
    that redirects or swallows stdout, and which pytest and CI surface on their own).
    """
    kept, dropped = [], []
    for pdb in pdbs:
        (dropped if is_excluded(pdb) else kept).append(pdb)
    where = f" from {context}" if context else ""
    for pdb in dropped:
        message = (f"EXCLUDED{where}: {INTRACTABLE[code_of(pdb)]}. "
                   f"Set {ALLOW_ENV}=1 to run it anyway.")
        if on_note:
            on_note("[foldx] " + "!" * 3 + " " + message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    return kept


def filter_worklist(worklist: Mapping[str, Iterable[str]], on_note=print,
                    context: str = "") -> Dict[str, List[str]]:
    """:func:`filter_complexes` over a ``{pdb: mutations}`` worklist, announcing malformed rows.

    Only the compute guard removes anything. An intractable complex is dropped whole, because
    FoldX cannot compute it; a malformed row is **kept and announced**, because FoldX computes it
    perfectly well. What is wrong with such a row is its arity label, not its chemistry: the
    string names three substitutions and FoldX returns the energy of exactly those three. Every
    comparator pipeline whose source was read computes the same mutant.

    Dropping them here would push the defect somewhere it cannot be seen. The store would carry no
    record, so :func:`is_malformed`, ``FoldxLookup.is_malformed`` and the load-time warning could
    never fire on anything -- a trust guard over an empty set -- and a consumer joining the store
    onto SKEMPI would lose two rows with nothing to explain the shortfall. Worse, removing an
    entry from a mutation list changes the values of every entry after it, so dropping two rows
    would silently move the other 17 in that complex's list.

    **An empty mutation list is passed through, not dropped.** Nothing here can empty one -- the
    malformed announce removes no entries -- so an empty list is a caller's instruction, and the
    one caller that gives it is a repair-only campaign, which has no mutations to build by
    definition. Dropping those left ``repair_ablation probe`` and both sweeps' repair phases
    enqueuing nothing and exiting 0, and the sweep then advised re-running to finish a chain that
    re-running could never start.
    """
    keep = set(filter_complexes(list(worklist), on_note=on_note, context=context))
    out = {}
    for pdb, muts in worklist.items():
        if pdb not in keep:
            continue
        muts = list(muts)
        announce_malformed(code_of(pdb), muts, on_note=on_note, context=context)
        out[pdb] = muts
    return out


@dataclass(frozen=True)
class MalformedRow:
    """A SKEMPI row whose ``Mutation(s)_cleaned`` field is internally inconsistent."""

    pdb: str
    cleaned: str
    """The mutation string exactly as SKEMPI ships it."""
    intended: str
    """What the row almost certainly meant, from the surrounding rows. NOT applied."""
    reason: str
    evidence: str

    def __str__(self) -> str:
        return (f"{self.pdb} {self.cleaned}: {self.reason} "
                f"(intended {self.intended}; {self.evidence})")


#: SKEMPI rows that are malformed **upstream**, keyed by ``(pdb, cleaned)``.
#:
#: Distinct from :data:`INTRACTABLE`, which is a whole complex FoldX cannot compute, and from
#: a correct row scored against the wrong pairing, which the identifier-keyed store removes. Here the
#: *row itself* is wrong before FoldX sees it, so the computed value is faithful to a mutant nobody
#: intended.
#:
#: Granularity is the mutation, not the complex, and that is the point. Both rows sit in
#: ``2C5D_AB_CD``, which holds 19 multi-point variants; excluding that record to drop two would
#: discard 17 sound ones. Excluding the whole code would discard 33, since ``2C5D``'s other
#: definition carries 12 more multi-point variants and 4 single-point mutations.
#:
#: **These are not corrected.** Substituting the intended string would put a value in the store
#: that corresponds to no SKEMPI row, which is worse than a documented defect: a consumer joining
#: on SKEMPI identifiers would silently match a mutant that was never requested.
#:
#: **How the field handles them, surveyed against the published source of each comparator
#: pipeline.** No pipeline drops these rows; every block list found is ``{'1KBH'}`` and nothing
#: more. All
#: of them end up computing the same **3-substitution** mutant computed here, by two different
#: routes:
#:
#: * ``PPIformer`` dedupes and documents it -- a code comment names ``RA32E,KA34E,RB32E,KA34E``
#:   in ``2C5D_AB_CD`` verbatim, the only such note found in any source read.
#: * ``USP-ddG`` dedupes via ``set(mut_str.split(','))``, without a comment. That is load-bearing
#:   rather than incidental: a later assertion compares the mutation count against the number of
#:   differing residues, and a list would make it fire ``4 != 3``. Its redistributed table records
#:   ``num_muts = 3`` for both rows -- exactly 2 of its 6706 rows disagree with their own comma
#:   count, and they are these two.
#: * ``RDE-Network`` and ``MINT`` keep a list and apply the substitution twice. The second write
#:   is idempotent, so the mutant is identical; they record ``num_muts = 4``.
#:
#: **So the energies here are not anomalous -- the label is.** Every pipeline scores the same
#: 3-substitution mutant, and the field disagrees only on whether to call it 3 or 4. Cross-paper
#: comparisons stratified by mutation arity are therefore comparing tables that disagree about
#: these two rows. A 2-row effect and numerically negligible, but a real definitional mismatch if
#: per-arity breakdowns are quoted.
MALFORMED_ROWS: Dict[tuple, MalformedRow] = {
    ("2C5D", "RA32E,KA34E,RB32E,KA34E"): MalformedRow(
        pdb="2C5D", cleaned="RA32E,KA34E,RB32E,KA34E", intended="RA32E,KA34E,RB32E,KB34E",
        reason="a symmetric-dimer variant repeats chain A's substitution instead of chain B's, "
               "so it is a 3-substitution mutant labelled as a 4-substitution one",
        evidence="scanning all 7085 rows of skempi_v2.csv for a repeated substitution returns "
                 "exactly two, both this complex; the neighbouring rows KA34E,KB34E and "
                 "RA32E,RB32E are cleanly symmetric, and no row anywhere mutates one position "
                 "to two different residues",
    ),
    ("2C5D", "EC30R,EC33R,ED30R,EC33R"): MalformedRow(
        pdb="2C5D", cleaned="EC30R,EC33R,ED30R,EC33R", intended="EC30R,EC33R,ED30R,ED33R",
        reason="same defect on the C/D chain pair: ED33R is replaced by a repeat of EC33R",
        evidence="neighbouring rows EC30R,ED30R and EC33R,ED33R are cleanly symmetric",
    ),
}


def is_malformed(pdb: str, cleaned: str) -> Optional[MalformedRow]:
    """The :class:`MalformedRow` for a ``(pdb, cleaned)`` pair, or ``None``.

    Unlike :func:`is_excluded` this is **not** bypassed by :data:`ALLOW_ENV`. That variable exists
    to re-test a complex against a newer FoldX, and a newer engine cannot repair a defect in the
    source table -- the row will still say what it says.
    """
    return MALFORMED_ROWS.get((code_of(pdb), cleaned))


def announce_malformed(pdb: str, mutations: Iterable[str], on_note=print,
                       context: str = "") -> List[str]:
    """Name the malformed rows in one complex's list, loudly, and return them.

    Nothing is removed. Announced through both channels for the same reason as
    :func:`filter_complexes`: a defect nobody is told about is how a coverage number becomes
    unaccountable -- and here the number is right while a label on two of its rows is not.
    """
    code = code_of(pdb)
    found = [m for m in mutations if is_malformed(code, m)]
    where = f" in {context}" if context else ""
    for mut in found:
        message = (f"MALFORMED{where}: {MALFORMED_ROWS[(code, mut)]}. "
                   f"Computed as written and shipped; the energies are right for the three "
                   f"substitutions named and wrong for the four the row claims. Not correctable "
                   f"without changing what SKEMPI says.")
        if on_note:
            on_note("[foldx] " + "!" * 3 + " " + message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    return found
