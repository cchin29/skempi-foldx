"""One object that answers "what is the FoldX energy for this mutation", whatever it is called.

A mutation in this domain has more than one name, and that is the source of nearly every way a
consumer gets a wrong answer here:

    LI38D          SKEMPI's author-chain form -- chain I, residue 38
    LB38D          the role-chain form -- partner B, renumbered along its group
    LI38D,GI32Y    a multi-point variant, comma-joined, in either form

The shipped single-point store contains **both** conventions: 3012 records keyed the first way,
1226 the second. A lookup that matches strings therefore misses about 29% of it, silently: asked
for SKEMPI's 4337 single-point mutations by their author form, such a lookup reports 3012 covered,
69.4%, where the store in fact holds 4238 of them, 97.7%.

:class:`FoldxLookup` resolves the naming so a caller does not have to::

    fx = FoldxLookup()                                    # store only
    fx = FoldxLookup(skempi_csv=..., mapping_dir=...)      # + exact role-chain resolution

    fx.get("1ACB", "LI38D")        # author form
    fx.get("1ACB", "LB38D")        # role form -- same record
    fx.coverage(rows)              # (covered, total, missing)

What it can resolve depends on what it was given, and it says so rather than failing quietly:
:attr:`conventions` lists the forms it will match, and a miss that exists under another name
raises an error naming that form and what to pass to resolve it automatically.
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .skempi import Mutation, load_skempi
from .store import MULTI_POINT, SINGLE_POINT, load_bundled_store
from .terms import TERMS


def load_chain_mapping(mapping_dir, code: str) -> Optional[dict]:
    """``chain -> {"seq": [...]}`` from a SKEMPI ``<code>.mapping`` file.

    Lines are ``RESNAME CHAIN AUTHOR_NUMBER SEQUENCE_INDEX``. Only the chain and the number of
    residues per chain matter here: the role-chain remap needs each chain's length to accumulate
    the offsets a split file numbers by.
    """
    path = Path(mapping_dir) / f"{code}.mapping"
    if not path.exists():
        return None
    chains: Dict[str, dict] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 4:
            chains.setdefault(parts[1], {"seq": []})["seq"].append(parts[2])
    return chains or None


def to_role_form(mutation: str, groups: Tuple[str, str], chains: dict) -> Optional[str]:
    """An author-chain mutation in the role numbering a split file uses, or ``None``.

    Splits name the partners ``A`` and ``B`` and number residues consecutively along each group,
    so ``LI38D`` becomes ``LB38D`` when chain ``I`` is the first chain of group 2.

    Identity matching cannot do this reliably: where both partners admit the same substitution it
    has no way to choose. The offsets do, which is why the mapping files turn a lower bound into
    an exact answer. Returns ``None`` rather than guessing when the mutation carries an insertion
    code or falls outside the chain.
    """
    if not groups or chains is None or len(mutation) < 4:
        return None
    chain = mutation[1]
    g1, g2 = groups
    if chain in g1:
        side, group = "A", g1
    elif chain in g2:
        side, group = "B", g2
    else:
        return None
    residue = mutation[2:-1]
    if not residue.isdigit():
        return None                  # insertion codes are dropped by the split builders too
    offsets, running = {}, 0
    for c in group:
        if c not in chains:
            return None
        offsets[c] = running
        running += len(chains[c]["seq"])
    n = int(residue)
    if not 1 <= n <= len(chains[chain]["seq"]):
        return None
    return f"{mutation[0]}{side}{offsets[chain] + n}{mutation[-1]}"


class FoldxLookup:
    """Resolve a mutation to its FoldX record regardless of naming convention.

    ``skempi_csv`` and ``mapping_dir`` are optional and additive:

    =========================== ==========================================================
    given                       resolves
    =========================== ==========================================================
    nothing                     stored keys, and SKEMPI's author form via each record's
                                ``cleaned`` field
    ``+ skempi_csv``            also role-chain forms, by matching ``(wt, position,
                                mutant)`` identity -- exact for most complexes, ambiguous
                                where both partners admit the same substitution
    ``+ mapping_dir``           also role-chain forms exactly, by residue offset
    =========================== ==========================================================

    Supplying both is worth it whenever labels come from split files. Without them a coverage
    figure over such labels is a lower bound: measured on a 1634-row multi-point split, 37.6% from
    the store alone and 75.8% by identity matching, against an exact 98.8%.
    """

    def __init__(self, arms: Sequence[str] = (SINGLE_POINT, MULTI_POINT),
                 skempi_csv=None, mapping_dir=None):
        self._records: Dict[Tuple[str, str], dict] = {}
        self._arm: Dict[Tuple[str, str], str] = {}
        self._alias: Dict[Tuple[str, str], Tuple[str, str]] = {}
        self._mappings_read = 0
        self.arms = tuple(arms)
        self.skempi = load_skempi(Path(skempi_csv)) if skempi_csv else None
        self.mapping_dir = Path(mapping_dir) if mapping_dir else None

        for which in self.arms:
            arm = "sp" if which == SINGLE_POINT else "mp"
            for pdb, recs in load_bundled_store(which).items():
                for key, rec in recs.items():
                    self._records[(pdb, key)] = rec
                    self._arm[(pdb, key)] = arm
                    cleaned = rec.get("cleaned")
                    if cleaned and cleaned != key:
                        self._alias.setdefault((pdb, cleaned), (pdb, key))

        if self.skempi:
            self._add_role_forms()

    # ------------------------------------------------------------------ construction helpers
    def _add_role_forms(self):
        cache: Dict[str, Optional[dict]] = {}
        for (pdb, key), rec in list(self._records.items()):
            entry = self.skempi.get(pdb)
            if entry is None:
                continue
            author = rec.get("cleaned", key)
            groups = (entry.group1, entry.group2)
            role = None
            if self.mapping_dir is not None:
                if pdb not in cache:
                    cache[pdb] = load_chain_mapping(self.mapping_dir, pdb)
                    if cache[pdb]:
                        self._mappings_read += 1
                parts = [to_role_form(a, groups, cache[pdb]) for a in author.split(",")]
                role = ",".join(parts) if all(parts) else None
            if role is None:
                role = self._role_by_identity(entry, author)
            if role and (pdb, role) not in self._records:
                self._alias.setdefault((pdb, role), (pdb, key))

    def _role_by_identity(self, entry, author: str) -> Optional[str]:
        """Fallback when no mapping file is available: match on substitution identity."""
        pool: Dict[tuple, List[str]] = {}
        for cleaned in set(entry.single) | {c for v in entry.multi for c in v.split(",")}:
            try:
                mut = Mutation.parse(cleaned)
            except ValueError:
                continue     # insertion code or malformed row: it can name nothing here anyway
            pool.setdefault(mut.identity, []).append(cleaned)
        out = []
        for part in author.split(","):
            try:
                mut = Mutation.parse(part)
            except Exception:
                return None
            if mut.chain in entry.group1:
                side = "A"
            elif mut.chain in entry.group2:
                side = "B"
            else:
                return None      # chain outside both groups; "B" here would be a guess
            if len(pool.get(mut.identity, [])) != 1:
                return None          # ambiguous across partners; do not guess
            out.append(f"{mut.wt}{side}{mut.position}{mut.mutant}")
        return ",".join(out)

    # ------------------------------------------------------------------------------ lookups
    @property
    def conventions(self) -> List[str]:
        """Which naming forms this instance will match, in the order they were added.

        The exact form is claimed only when mapping files were actually read, not merely when
        ``mapping_dir`` was supplied: a path that is empty or misspelt degrades silently to
        identity matching, and a report that still called itself exact would be a lower bound
        wearing the wrong label.
        """
        forms = ["stored key", "SKEMPI author form"]
        if self.skempi and self._mappings_read:
            forms.append("role-chain form (exact, via residue offsets)")
        elif self.skempi:
            forms.append("role-chain form (identity matching; ambiguous cases unresolved)")
        return forms

    @property
    def is_exact(self) -> bool:
        """Whether role-chain labels resolve by residue offset rather than by identity."""
        return bool(self.skempi and self._mappings_read)

    def _resolve(self, pdb: str, mutation: str) -> Optional[Tuple[str, str]]:
        key = (pdb, mutation)
        if key in self._records:
            return key
        return self._alias.get(key)

    def __contains__(self, row) -> bool:
        return self._resolve(*row) is not None

    def get(self, pdb: str, mutation: str, default=None):
        """The 12-term record, or ``default``. Accepts any convention this instance resolves."""
        found = self._resolve(pdb, mutation)
        return self._records[found] if found else default

    def vector(self, pdb: str, mutation: str) -> Optional[List[float]]:
        """The record flattened into ``TERMS`` order, for use as a feature vector."""
        rec = self.get(pdb, mutation)
        return [float(rec[t]) for t in TERMS] if rec else None

    def arm_of(self, pdb: str, mutation: str) -> Optional[str]:
        """``"sp"`` or ``"mp"``, or ``None`` when the mutation is not covered."""
        found = self._resolve(pdb, mutation)
        return self._arm[found] if found else None

    def require(self, pdb: str, mutation: str) -> dict:
        """Like :meth:`get`, but raises with the reason -- and the fix -- instead of returning
        ``None``.

        A miss is usually a naming mismatch rather than absent data, so the error reports how many
        mutations the complex does hold, gives examples of the form they are filed under, and names
        the constructor argument that would widen resolution. Matching one convention and calling
        the remainder absent is what turns a 97.7%-covered store into a 69.4% figure.
        """
        rec = self.get(pdb, mutation)
        if rec is not None:
            return rec
        others = [k for (p, k) in self._records if p == pdb]
        if not others:
            raise KeyError(f"{pdb} is not in the shipped store "
                           f"(arms loaded: {', '.join(self.arms)})")
        hint = ""
        if self.skempi is None:
            hint = (" This instance resolves only stored keys and SKEMPI author forms; pass "
                    "skempi_csv= and mapping_dir= to also resolve role-chain labels.")
        elif self.mapping_dir is None:
            hint = (" Role-chain labels are being resolved by identity matching, which cannot "
                    "disambiguate when both partners admit the same substitution; pass "
                    "mapping_dir= for the exact residue-offset remap.")
        raise KeyError(
            f"{pdb} has no record named {mutation!r}. It holds {len(others)} mutations, e.g. "
            f"{sorted(others)[:3]}.{hint}"
        )

    def coverage(self, rows: Iterable[Tuple[str, str]]):
        """``(covered, total, missing)`` over ``(pdb, mutation)`` pairs."""
        rows = list(rows)
        missing = [r for r in rows if r not in self]
        return len(rows) - len(missing), len(rows), missing

    def energies_for(self, rows: Iterable[Tuple[str, str]]) -> Dict[Tuple[str, str], List[float]]:
        """``{(pdb, mutation): [12 floats]}`` for the rows that resolve, keyed as given."""
        out = {}
        for row in rows:
            vec = self.vector(*row)
            if vec is not None:
                out[row] = vec
        return out

    def __len__(self) -> int:
        return len(self._records)

    def __repr__(self) -> str:
        return (f"FoldxLookup({len(self._records)} records, arms={list(self.arms)}, "
                f"resolves={self.conventions})")
