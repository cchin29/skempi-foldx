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

import math
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .skempi import Mutation, load_skempi
from .exclusions import COLLAPSED_INTERFACES, N_INTERFACE_SUSPECT, interface_suspect
from .store import (MULTI_KEY, MULTI_POINT, SINGLE_POINT, infer_kind, load_bundled_store,
                    load_store, store_kind)
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
                 skempi_csv=None, mapping_dir=None, stores=None):
        self._records: Dict[Tuple[str, str], dict] = {}
        self._arm: Dict[Tuple[str, str], str] = {}
        self._alias: Dict[Tuple[str, str], Tuple[str, str]] = {}
        self._mappings_read = 0
        self._mappings_wanted = 0
        self.arms = tuple(arms)
        # load_skempi's per-complex warning and the aggregate raised below are not the same fact,
        # so neither is suppressed: the first describes the table the caller supplied, naming its
        # pairings and per-complex counts, and the second describes the shipped store, which a
        # caller-supplied table can collapse in ways the store's static registry does not know.
        self.skempi = load_skempi(Path(skempi_csv)) if skempi_csv else None
        self.mapping_dir = Path(mapping_dir) if mapping_dir else None

        conflicts: List[Tuple[str, str]] = []
        sources = self._sources(stores, self.arms)
        if stores is not None:
            # Otherwise `arms` still describes the bundled default, and every message derived from
            # it -- the arm diagnosis in require(), the "arms loaded" list -- describes a store
            # this instance never opened.
            seen = {arm for arm, _ in sources}
            self.arms = tuple(w for w in (SINGLE_POINT, MULTI_POINT)
                              if ("sp" if w == SINGLE_POINT else "mp") in seen)
        for arm, source in sources:
            for pdb, recs in source.items():
                for key, rec in recs.items():
                    if (pdb, key) in self._records or (pdb, key) in self._alias:
                        # Alias as well as key: a store keyed by author form collides with an
                        # earlier store's `cleaned` alias rather than with its key, and that
                        # shadowing is invisible at key level.
                        conflicts.append((pdb, key))
                        continue          # first store wins, as consolidate() does
                    self._records[(pdb, key)] = rec
                    self._arm[(pdb, key)] = arm
                    cleaned = rec.get("cleaned")
                    if cleaned and cleaned != key:
                        if (pdb, cleaned) in self._records:
                            # An earlier record already owns that name as a key, so this record's
                            # author form cannot alias onto it. The record itself is kept; only
                            # the alias is refused, and silence here would hide two stores
                            # disagreeing about one mutation.
                            conflicts.append((pdb, cleaned))
                        else:
                            self._alias.setdefault((pdb, cleaned), (pdb, key))

        if conflicts:
            warnings.warn(
                f"{len(conflicts)} mutations resolve to more than one record across the supplied "
                f"stores, e.g. {conflicts[:3]}. The first is kept, matching consolidate(); the "
                f"others are discarded without being compared.",
                RuntimeWarning, stacklevel=2)

        if self.skempi:
            self._add_role_forms()

        self.interface_suspects = frozenset(
            (p, k) for (p, k) in self._records
            if interface_suspect(p, self._records[(p, k)].get("cleaned", k)))
        if self.interface_suspects:
            codes = sorted({p for p, _ in self.interface_suspects})
            warnings.warn(
                f"{len(self.interface_suspects)} indexed records were scored against an interface "
                f"SKEMPI does not pair them with ({', '.join(codes)}): SKEMPI defines these codes "
                f"under more than one pairing, and the store keys on the code alone. Most are "
                f"all-zero as an artifact of the pairing, and the rest carry a partial energy "
                f"that does not look wrong. Test a record with "
                f"FoldxLookup.is_interface_suspect(pdb, mutation).",
                RuntimeWarning, stacklevel=2)

    # ------------------------------------------------------------------ construction helpers
    @staticmethod
    def _sources(stores, arms):
        """``(arm, {pdb: {mutation: record}})`` for each store to index.

        ``stores`` is ``None`` for the bundled data, or a results directory, an in-memory store,
        or a sequence of either. A campaign's own output is the case that matters: the naming
        conventions are the same problem there as here, so a store this class cannot index is one
        its owner has to join by hand -- which is where the wrong answers come from.

        The arm is read from the store rather than declared, because it is already recorded: on
        disk by the key each file uses, in memory by whether any record name holds a comma.
        """
        if stores is None:
            return [("sp" if which == SINGLE_POINT else "mp", load_bundled_store(which))
                    for which in arms]
        if isinstance(stores, (str, Path, dict)):
            stores = [stores]
        out = []
        for item in stores:
            if isinstance(item, (str, Path)):
                kind = store_kind(Path(item))
                out.append(("mp" if kind == MULTI_KEY else "sp", load_store(Path(item))))
            else:
                out.append(("mp" if infer_kind(item) == MULTI_KEY else "sp", item))
        return out

    def _add_role_forms(self):
        cache: Dict[str, Optional[dict]] = {}
        for (pdb, key), rec in list(self._records.items()):
            entry = self.skempi.get(pdb)
            if entry is None:
                continue
            # Counted before any refusal below, so a complex this method declines to resolve
            # still registers as one whose mapping was wanted. Otherwise exactness is claimed over
            # precisely the complexes that got no role resolution at all.
            if self.mapping_dir is not None and pdb not in cache:
                cache[pdb] = load_chain_mapping(self.mapping_dir, pdb)
                self._mappings_wanted += 1
                if cache[pdb]:
                    self._mappings_read += 1

            author = rec.get("cleaned")
            if author is None:
                # No `cleaned`, so the key is being *inferred* to be an author form. Only a
                # mutation SKEMPI lists for this complex may be inferred that way: otherwise a
                # role-keyed store has its key taken for an author form and remapped again,
                # aliasing a real record onto a residue it does not describe. The bundled
                # multi-point store carries no `cleaned`, so the shape is not hypothetical.
                author = key
                if author not in entry.single and author not in entry.multi:
                    continue
            # A record carrying `cleaned` declares its own author form, so nothing is inferred and
            # SKEMPI membership is irrelevant to it. Requiring membership here would discard
            # correct aliases whenever the supplied table is a subset of the one that built the
            # store -- silently, and by 9 points of coverage on a real split.
            groups = (entry.group1, entry.group2)
            role = None
            if self.mapping_dir is not None:
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
        if self.is_exact:
            forms.append("role-chain form (exact, via residue offsets)")
        elif self.skempi and self._mappings_read:
            forms.append(f"role-chain form (partly exact: {self._mappings_read} of "
                         f"{self._mappings_wanted} mappings read, remainder by identity)")
        elif self.skempi:
            forms.append("role-chain form (identity matching; ambiguous cases unresolved)")
        return forms

    def is_interface_suspect(self, pdb: str, mutation: str) -> bool:
        """Whether this record was scored against an interface SKEMPI does not pair it with.

        Accepts any convention this instance resolves, so a split-file label answers the same as
        a stored key. See :data:`skempi_foldx.COLLAPSED_INTERFACES`.
        """
        found = self._resolve(pdb, mutation)
        return found in self.interface_suspects if found else False

    @property
    def is_exact(self) -> bool:
        """Whether role-chain labels resolve by residue offset rather than by identity.

        Requires a mapping for *every* complex that needed one. A directory holding some of them
        resolves the rest by identity, and reporting that as exact would suppress the lower-bound
        warning over precisely the rows still estimated.
        """
        return bool(self.skempi and self._mappings_wanted
                    and self._mappings_read == self._mappings_wanted)

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
        if rec is None:
            return None
        out = [float(rec[t]) for t in TERMS]
        # json.loads accepts bare NaN and Infinity, so a results directory can carry them all the
        # way into a feature vector, where they are far harder to trace back.
        bad = [t for t, v in zip(TERMS, out) if not math.isfinite(v)]
        if bad:
            raise ValueError(f"{pdb} {mutation!r} has non-finite values for {bad}; "
                             f"a FoldX term is always a finite number")
        return out

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
            # Diagnose before declaring absence: a case-only mismatch and an unloaded arm both
            # present as "not in the store", and both send a reader looking in the wrong place.
            cased = [p for (p, _) in self._records if p.upper() == pdb.upper()]
            if cased:
                raise KeyError(f"{pdb} is not in the store, but {cased[0]} is -- complex codes "
                               f"are upper case here.")
            raise KeyError(f"{pdb} is not in the store "
                           f"(arms loaded: {', '.join(self.arms) or 'none'})")
        if "," in mutation and MULTI_POINT not in self.arms:
            raise KeyError(
                f"{mutation!r} is a multi-point variant and only the single-point arm is loaded. "
                f"Construct with arms=(SINGLE_POINT, MULTI_POINT) -- the default -- or "
                f"arms=(MULTI_POINT,).")
        if "," not in mutation and SINGLE_POINT not in self.arms:
            raise KeyError(
                f"{mutation!r} is a single-point mutation and only the multi-point arm is loaded. "
                f"Construct with arms=(SINGLE_POINT, MULTI_POINT) -- the default.")
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
        """``(covered, total, missing)`` over ``(pdb, mutation)`` pairs.

        Warns when a large share fails to resolve and this instance could resolve more of them
        given more inputs, since a shortfall that size is worth interrupting for while a handful
        of genuinely uncomputed rows is not: on a 4159-row split the same rows read 65.4% here and
        99.6% with ``skempi_csv`` and ``mapping_dir`` supplied.

        **The warning describes this lookup, not the rows.** Five successive attempts to have it
        classify why a row missed -- guessing from chain letters, from the arm's shape, from the
        code's case -- each shipped a statement that was false of rows a consumer really had. Every
        such test is a property of the label, and the label does not carry the answer. What is
        knowable without inspecting a row is what this instance can match and what would widen it,
        so that is all this says; :meth:`require` diagnoses a single row, where the record itself
        is available to check against.

        ``is_exact`` suppresses it, which is a statement about the mappings this instance was
        asked for rather than about the rows: a ``skempi_csv`` covering one complex reports exact
        once that one mapping is read, and an exact single-arm lookup is silent over every row of
        the arm it never indexed.
        """
        rows = list(rows)
        missing = [r for r in rows if r not in self]
        if missing and not self.is_exact and len(rows) >= 100 and len(missing) >= 0.05 * len(rows):
            if self.skempi is None:
                widen = (" skempi_csv= and mapping_dir= widen it to the role-chain labels split "
                         "files use.")
            elif self.mapping_dir is None:
                widen = " mapping_dir= makes its role-chain resolution exact rather than partial."
            elif self._mappings_read < self._mappings_wanted:
                widen = (f" Its role-chain resolution is partial: {self._mappings_read} of "
                         f"{self._mappings_wanted} mapping files were read.")
            else:
                # skempi_csv and mapping_dir are both set and everything asked for was read, yet
                # is_exact is false -- no mapping was wanted, because no indexed complex is in the
                # table. Naming an input to widen would be a guess about which one is narrow.
                widen = (" No complex it indexes appears in the supplied skempi_csv=, so no "
                         "mapping was read.")
            warnings.warn(
                f"{len(missing)} of {len(rows)} rows ({100 * len(missing) / len(rows):.1f}%) did "
                f"not resolve. This lookup matches: {', '.join(self.conventions)}.{widen} "
                f"require(pdb, mutation) gives the reason for any one row.",
                RuntimeWarning, stacklevel=2)
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
