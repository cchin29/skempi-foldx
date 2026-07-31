"""Reading SKEMPI and mapping dataset mutations onto real PDB chains.

The mapping is the subtle part of the whole pipeline. A dataset row names its mutation by
**role** — ``A`` is the first partner, ``B`` the second — while FoldX needs the **author** chain
letter as it appears in the structure. SKEMPI's ``Mutation(s)_cleaned`` column carries the author
form and uses the same residue numbering, so the two are matched on ``(wild-type, position,
mutant)`` and disambiguated by role when a complex has the same substitution on both partners.

Every mapping is then validated against the wild-type residue actually present in the repaired
structure, so a numbering mismatch fails loudly instead of scoring the wrong residue.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, Tuple

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


@dataclass(frozen=True)
class Mutation:
    """One substitution: ``LI38S`` -> wt ``L``, chain ``I``, position 38, mutant ``S``.

    ``chain`` is an author chain letter for a cleaned SKEMPI mutation, and a role letter
    (``A``/``B``) for a dataset mutation. Which one it is depends on where the string came
    from — that ambiguity is the reason :func:`map_role_to_author` exists.
    """

    wt: str
    chain: str
    position: int
    mutant: str

    @classmethod
    def parse(cls, text: str) -> "Mutation":
        return cls(wt=text[0], chain=text[1], position=int(text[2:-1]), mutant=text[-1])

    def __str__(self) -> str:
        return f"{self.wt}{self.chain}{self.position}{self.mutant}"

    @property
    def identity(self) -> Tuple[str, int, str]:
        """The chain-independent part, used to match a role mutation to an author one."""
        return (self.wt, self.position, self.mutant)


@dataclass
class SkempiComplex:
    """One PDB's SKEMPI record: its two interacting chain groups and its mutations."""

    pdb: str
    group1: str
    group2: str
    single: Set[str] = field(default_factory=set)
    """Single-point cleaned mutation strings."""
    multi: Set[str] = field(default_factory=set)
    """Multi-point variants, as the raw comma-joined cleaned string."""

    @property
    def groups(self) -> str:
        """The ``--analyseComplexChains`` argument."""
        return f"{self.group1},{self.group2}"


def load_skempi(path: Path) -> Dict[str, SkempiComplex]:
    """Parse SKEMPI 2.0 into ``{pdb: SkempiComplex}``.

    Single- and multi-point rows are kept separately rather than one being dropped: the
    single-point pipeline and the multi-point one both read this, so filtering either out here
    would force the other arm to re-parse the table for itself.
    """
    complexes: Dict[str, SkempiComplex] = {}
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter=";")
        next(reader, None)  # header
        for row in reader:
            if not row or not row[0]:
                continue
            pdb_field, cleaned = row[0], row[2]      # '<PDB>_<group1>_<group2>'
            parts = pdb_field.split("_")
            if len(parts) < 3:
                continue
            pdb, group1, group2 = parts[0], parts[1], parts[2]
            entry = complexes.setdefault(pdb, SkempiComplex(pdb, group1, group2))
            (entry.multi if "," in cleaned else entry.single).add(cleaned)
    return complexes


def map_role_to_author(
    entry: SkempiComplex, role_mutations: Sequence[str]
) -> Tuple[Dict[str, str], List[str]]:
    """Map role-chain mutations onto SKEMPI's author-chain (``cleaned``) strings.

    Returns ``(mapping, unresolved)``. A mutation resolves when exactly one cleaned mutation
    shares its ``(wt, position, mutant)``. When several do — the same substitution present on
    both partners — the role letter selects the chain group: role ``A`` must land in group 1,
    role ``B`` in group 2. If that still leaves more than one candidate, the mutation is left
    unresolved rather than guessed at.
    """
    by_identity: Dict[Tuple[str, int, str], List[str]] = {}
    for cleaned in entry.single:
        by_identity.setdefault(Mutation.parse(cleaned).identity, []).append(cleaned)

    mapping: Dict[str, str] = {}
    unresolved: List[str] = []
    for role_text in role_mutations:
        role_mut = Mutation.parse(role_text)
        candidates = by_identity.get(role_mut.identity, [])
        if len(candidates) == 1:
            mapping[role_text] = candidates[0]
            continue
        if len(candidates) > 1:
            group = entry.group1 if role_mut.chain == "A" else entry.group2
            narrowed = [c for c in candidates if Mutation.parse(c).chain in group]
            if len(narrowed) == 1:
                mapping[role_text] = narrowed[0]
                continue
        unresolved.append(role_text)
    return mapping, unresolved


def repaired_wt_residues(pdb_path: Path) -> Dict[Tuple[str, int], str]:
    """``(chain, position) -> one-letter residue`` from a structure's CA records.

    Used to confirm that the residue FoldX is about to mutate is the one the dataset says is
    there. Cheap, and it catches numbering-convention mismatches that would otherwise produce
    confidently wrong energies.
    """
    residues: Dict[Tuple[str, int], str] = {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            try:
                position = int(line[22:26])
            except ValueError:
                continue
            residues[(line[21], position)] = AA3TO1.get(line[17:20].strip(), "X")
    return residues


def validate_against_structure(
    mapping: Dict[str, str], residues: Dict[Tuple[str, int], str]
) -> Tuple[List[Tuple[str, str]], List[str]]:
    """Split a mapping into ``(validated pairs, rejected descriptions)``."""
    validated: List[Tuple[str, str]] = []
    rejected: List[str] = []
    for role_text, cleaned in mapping.items():
        mut = Mutation.parse(cleaned)
        found = residues.get((mut.chain, mut.position))
        if found == mut.wt:
            validated.append((role_text, cleaned))
        else:
            rejected.append(f"{role_text}->{cleaned}(structure has {found})")
    return validated, rejected


def parse_interaction_file(path: Path, terms: Sequence[str]) -> Dict[str, float]:
    """Read an ``AnalyseComplex`` Interaction ``.fxout`` into ``{term: value}``."""
    with open(path) as fh:
        lines = fh.readlines()
    header = next(l for l in lines if l.startswith("Pdb\t")).rstrip("\n").split("\t")
    data = next(l for l in lines if l.split("\t")[0].endswith(".pdb")).rstrip("\n").split("\t")
    row = dict(zip(header, data))
    return {t: float(row[t]) for t in terms if t in row}
