"""Complexes that must never be handed to FoldX.

This is not a curation filter. Curation already drops these complexes from every split — the
whole point is that they are *absent* from the downstream tables, which is exactly why nothing
downstream stops a compute driver from picking them up.

The failure this prevents
-------------------------
``1KBH`` cannot be repaired. RepairPDB on it has been left running past 22.5 h without finishing,
twice, on two different machines. Killing the job works, but killing is a per-round remedy and
the driver is a loop: the next round enumerates its scope, sees 1KBH again, and launches RepairPDB
again. An intractable complex has to be dropped *before* the worklist is built, not reaped after.

Why a driver sees it at all, given curation drops it
----------------------------------------------------
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
That asymmetry is intended — you should be able to test a complex without silently changing which
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
    """Whether the exclusion list is currently bypassed."""
    return os.environ.get(ALLOW_ENV, "") not in ("", "0", "false", "False")


def is_excluded(pdb: str) -> Optional[Exclusion]:
    """The :class:`Exclusion` for ``pdb``, or ``None``. Always ``None`` when bypassed."""
    if allowed():
        return None
    return INTRACTABLE.get(pdb)


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
        message = (f"EXCLUDED{where}: {INTRACTABLE[pdb]}. "
                   f"Set {ALLOW_ENV}=1 to run it anyway.")
        if on_note:
            on_note("[foldx] " + "!" * 3 + " " + message)
        warnings.warn(message, RuntimeWarning, stacklevel=2)
    return kept


def filter_worklist(worklist: Mapping[str, Iterable[str]], on_note=print,
                    context: str = "") -> Dict[str, List[str]]:
    """:func:`filter_complexes` over a ``{pdb: mutations}`` worklist."""
    keep = set(filter_complexes(list(worklist), on_note=on_note, context=context))
    return {pdb: list(muts) for pdb, muts in worklist.items() if pdb in keep}
