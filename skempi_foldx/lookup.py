"""One object that answers "what is the FoldX energy for this mutation", whatever it is called.

A mutation in this domain has more than one name, and that is the source of nearly every way a
consumer gets a wrong answer here:

    LI38D          SKEMPI's author-chain form -- chain I, residue 38
    LB38D          the role-chain form -- partner B, renumbered along its group
    LI38D,GI32Y    a multi-point variant, comma-joined, in either form

Every shipped key is the author form. A record whose role name *differs* from its key carries
that name in a ``role`` field; the rest need none, because their role form is the key. So a plain
dictionary join keyed on the author form works, and one keyed on the role form — which is what a
split file carries — matches only the subset where the two coincide, silently, and still writes
well-formed output. That is the failure this class exists to prevent.

There is a second name in play. The store is keyed by SKEMPI *interface definition*,
``<pdb>_<group1>_<group2>``, because three codes are defined under two pairings each and the
pairing is half of what a value means; a *record* is one mutation inside such a definition. A bare
code resolves wherever it is unambiguous, which is every code but those three.

:class:`FoldxLookup` resolves all of it so a caller does not have to::

    fx = FoldxLookup()                                    # store only
    fx = FoldxLookup(skempi_csv=..., mapping_dir=...)      # + role forms for a foreign store

    fx.get("1ACB", "LI38D")        # author form
    fx.get("1ACB", "LB38D")        # role form -- same record, from the shipped `role` field
    fx.get("1ACB_E_I", "LI38D")    # by identifier
    fx.definitions_of("2C5D")      # ['2C5D_AB_CD', '2C5D_A_C'] -- ask by identifier
    fx.coverage(rows)              # (covered, total, missing)

What it can resolve depends on what it was given, and it says so rather than failing quietly:
:attr:`conventions` lists the forms it will match, and a miss that exists under another name
raises an error naming that form and what to pass to resolve it automatically.
"""

from __future__ import annotations

import math
import re
import warnings
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .skempi import Mutation, load_skempi
from .exclusions import code_of, is_malformed
from .store import (MULTI_KEY, MULTI_POINT, SINGLE_POINT, check_code, infer_kind,
                    load_bundled_store, load_store, store_kind)
from .terms import TERMS


#: A ``<chain><author number>`` field that lost its separating space. See
#: :func:`load_chain_mapping`.
_FUSED_CHAIN_NUMBER = re.compile(r"[A-Za-z]\d+\Z")


def load_chain_mapping(mapping_dir, code: str) -> Optional[dict]:
    """``chain -> {"seq": [...]}`` from a SKEMPI ``<code>.mapping`` file.

    Lines are ``RESNAME CHAIN AUTHOR_NUMBER SEQUENCE_INDEX``. Only the chain and the number of
    residues per chain matter here: the role-chain remap needs each chain's length to accumulate
    the offsets a split file numbers by.

    **The columns are fixed-width, so a four-digit author number runs into the chain letter** and
    the line arrives as three fields, not four -- ``ALA C1001  1``. Splitting on whitespace and
    requiring four fields therefore drops every residue numbered 1000 or above, which is not a
    rare shape: 11 of 345 mapping files carry such lines. In six of them -- ``1S0W``, ``1XXM``,
    ``2NOJ``, ``2REX``, ``3Q3J``, ``5E6P`` -- an entire chain is numbered above 999, so that chain
    vanishes from the result; in the other five a chain is merely read short, ``2GOX`` chain A
    worst at 9 residues of 297.

    Both failures are silent and both corrupt the offsets rather than raising. A short chain makes
    :func:`to_role_form` reject residues that are really in range -- 14 records in ``2NYY``/
    ``2NZ9``, whose chain A parses as 971 residues instead of 1267. A missing chain makes it
    return ``None`` outright, and for a *later* chain in a group it would silently shift every
    subsequent offset, aliasing records onto residues they do not describe. So the field is split
    back apart here rather than guarded against downstream.

    Verified: with this parsing, recomputing the role name of all 6107 shipped records reproduces
    every one of them -- 0 additions, 0 losses, 0 disagreements.
    """
    path = Path(mapping_dir) / f"{check_code(code)}.mapping"
    if not path.exists():
        return None
    chains: Dict[str, dict] = {}
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) == 3 and _FUSED_CHAIN_NUMBER.match(parts[1]):
            parts = [parts[0], parts[1][0], parts[1][1:], parts[2]]
        # The residue number must be a number and a chain is one character. Without both, a
        # header row parses as a chain of its own and every offset after it shifts -- which is
        # how a mapping file silently comes to describe different numbering than the one that
        # built the store.
        if len(parts) >= 4 and len(parts[1]) == 1 and parts[2].isdigit():
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
    # A count, not a membership test, and deliberately. The mapping file's third column is the
    # author number and chains do not start at 1 -- `1ACB` chain I runs 8..70 over 63 residues --
    # so this bound is not the obvious one and looks like an off-by-one waiting to be fixed.
    # It is not: split files number by the author value, so `offsets + n` is the right arithmetic,
    # and replacing the bound with "is `n` one of this chain's author numbers" was measured to
    # lose the role name on 899 of the 6107 shipped records while recovering none -- counted per
    # record, since a multi-point record's role name is all-or-nothing across its parts.
    # The bound declines nothing in the shipped store: the only records it could turn away were
    # the `2NYY`/`2NZ9` ones whose chain A parsed short, and splitting the fused column above
    # removes that category. So the loose bound costs nothing here and a tight one costs 899.
    if not 1 <= n <= len(chains[chain]["seq"]):
        return None
    return f"{mutation[0]}{side}{offsets[chain] + n}{mutation[-1]}"


class FoldxLookup:
    """Resolve a mutation to its FoldX record regardless of naming convention.

    ``skempi_csv`` and ``mapping_dir`` are optional and additive:

    =========================== ==========================================================
    given                       resolves
    =========================== ==========================================================
    nothing                     stored keys, SKEMPI's author form via each record's
                                ``cleaned`` field, and role-chain forms via its ``role``
                                field -- exact, and all the bundled store needs
    ``+ skempi_csv``            also role-chain forms, by matching ``(wt, position,
                                mutant)`` identity -- exact for most complexes, ambiguous
                                where both partners admit the same substitution
    ``+ mapping_dir``           also role-chain forms exactly, by residue offset
    =========================== ==========================================================

    Neither is needed for the bundled store, whose records carry their own role names: measured
    from the package alone, a 1636-row multi-point split resolves all 1636 and a 4165-row
    single-point split all 4165, and supplying both arguments changes neither figure.

    They matter for a store this package did not build. Without a ``role`` field a coverage figure
    over role-chain labels is a lower bound, and :meth:`coverage` says so rather than presenting
    it as the answer.
    """

    def __init__(self, arms: Sequence[str] = (SINGLE_POINT, MULTI_POINT),
                 skempi_csv=None, mapping_dir=None, stores=None):
        self._records: Dict[Tuple[str, str], dict] = {}
        self._arm: Dict[Tuple[str, str], str] = {}
        self._alias: Dict[Tuple[str, str], Tuple[str, str]] = {}
        self._mappings_read = 0
        self._mappings_wanted = 0
        self._shipped_roles = 0
        self._unbacked_records = 0
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
        alias_source: Dict[Tuple[str, str], int] = {}
        for index, (arm, source) in enumerate(sources):
            # Whether this source was built by something that computes role names.
            #
            # The test is per source, and it has to be: a record with no `role` is not evidence of
            # a gap, because 2220 of the bundled store's 6107 legitimately have none -- their role
            # form *is* their key. Distinguishing "no role because none differs" from "no role
            # because none was computed" needs the chain mappings, which is exactly the input this
            # field exists to spare a consumer. So the signal available is the source, not the
            # record: a store carrying role names was built by a tool that computes them.
            #
            # What that cannot catch is a hand-assembled store with one role field among many,
            # which will claim exactness it has not earned. Documented on :attr:`is_exact` rather
            # than guessed at, since no threshold here is defensible.
            n_records = sum(len(recs) for recs in source.values())
            if n_records and not any(r.get("role") for recs in source.values()
                                     for r in recs.values()):
                self._unbacked_records += n_records
            for pdb, recs in source.items():
                for key, rec in recs.items():
                    if (pdb, key) in self._records:
                        conflicts.append((pdb, key))
                        continue          # first store wins, as consolidate() does
                    if (pdb, key) in self._alias:
                        conflicts.append((pdb, key))
                        if alias_source.get((pdb, key)) != index:
                            continue      # an earlier store owns the name; first store wins
                        # Same store: a real record outranks an alias standing on its name.
                        # Discarding the record instead loses it entirely and leaves the alias
                        # serving a sibling record's energies under it, decided by JSON key order.
                        del self._alias[(pdb, key)]
                    self._records[(pdb, key)] = rec
                    self._arm[(pdb, key)] = arm
                    # The shipped role-chain name, where the store carries one. It is the alias
                    # split-file labels actually hit, and it is a field rather than a second key
                    # convention so the store stays single-keyed while the join still works with
                    # no skempi_v2.csv -- which this package does not redistribute.
                    role = rec.get("role")
                    if role:
                        self._shipped_roles += 1
                        if role != key and (pdb, role) not in self._records:
                            if self._alias.get((pdb, role), (pdb, key)) != (pdb, key):
                                # Two records of one complex claiming one role name. Whichever
                                # arrived first keeps it, so the other's energies would answer to
                                # a name that is not theirs -- decided by JSON key order, and
                                # silent until this was counted as a conflict like any other.
                                conflicts.append((pdb, role))
                            else:
                                alias_source[(pdb, role)] = index
                            self._alias.setdefault((pdb, role), (pdb, key))
                    cleaned = rec.get("cleaned")
                    if cleaned and cleaned != key:
                        if (pdb, cleaned) in self._records:
                            # An earlier record already owns that name as a key, so this record's
                            # author form cannot alias onto it. The record itself is kept; only
                            # the alias is refused, and silence here would hide two stores
                            # disagreeing about one mutation.
                            conflicts.append((pdb, cleaned))
                        elif self._alias.get((pdb, cleaned), (pdb, key)) != (pdb, key):
                            # Same collision as the role branch above, and it was silent here for
                            # as long as it was silent there: two records claiming one author
                            # form, resolved by whichever JSON key came first.
                            #
                            # Compared against the target rather than tested for presence. A
                            # record whose `role` and `cleaned` are the same string -- both
                            # differing from its key -- registers the alias in the branch above
                            # and then meets its own entry here. One record, one alias, correct
                            # resolution: reporting that as two records claiming one name sends a
                            # consumer looking for a collision that does not exist.
                            conflicts.append((pdb, cleaned))
                        else:
                            alias_source[(pdb, cleaned)] = index
                            self._alias[(pdb, cleaned)] = (pdb, key)

        #: Records whose shipped ``role`` disagrees with what ``mapping_dir=`` recomputes. Empty
        #: unless both a mapping directory and a role-bearing store are supplied.
        self.role_mismatches: Tuple[Tuple[str, str, str, str], ...] = ()

        #: One entry per *additional* claimant on a name, in the order found -- so a mutation
        #: present in three stores contributes two, matching ConsolidationReport.summary. The
        #: warning below is the loud form; this is the one a caller can test, for a consumer that
        #: narrows warnings to ``once`` as docs/USAGE.md recommends.
        self.conflicts: Tuple[Tuple[str, str], ...] = tuple(conflicts)

        if conflicts:
            warnings.warn(
                f"Names resolving to more than one record: {len(set(conflicts))} over "
                f"{len(conflicts)} claims, e.g. "
                f"{conflicts[:3]} -- across the supplied stores, or between two records of one "
                f"complex claiming the same alias. The first is kept, matching consolidate(); the "
                f"others are discarded without being compared.",
                RuntimeWarning, stacklevel=2)

        # Built before role forms are added, because that pass resolves each record's identifier
        # against the SKEMPI table and needs the code to find a chain mapping file.
        self._by_code: Dict[str, List[str]] = {}
        for ident, _ in self._records:
            bucket = self._by_code.setdefault(code_of(ident), [])
            if ident not in bucket:
                bucket.append(ident)
        for bucket in self._by_code.values():
            bucket.sort()

        #: Codes that are themselves a record key while also naming interface definitions --
        #: what mixing a pre-0.2.0 code-keyed store with an identifier-keyed one produces.
        #:
        #: This is the one route by which a bare doubled code still answers. `_resolve` returns
        #: an exact `_records` hit before it consults `_candidates`, so the ambiguity refusal
        #: never runs, and `("2C5D", m)` and `("2C5D_A_C", m)` are different keys, so nothing is
        #: recorded as a conflict either. The exact hit is the right answer when the legacy store
        #: is used alone -- the caller asked for that store's record -- and the wrong one when it
        #: shadows a definition-keyed record for the same mutation. Which of those it is depends
        #: on what the caller meant, so this reports rather than decides.
        self.pooled_codes: Tuple[str, ...] = tuple(sorted(
            code for code, idents in self._by_code.items()
            if len(idents) > 1 and code in idents))
        if self.pooled_codes:
            warnings.warn(
                f"{len(self.pooled_codes)} PDB code(s) are indexed both as a record key and as "
                f"interface definitions: {list(self.pooled_codes)[:3]}. A store keyed by code "
                f"pools every definition of that code under one name, so a bare-code lookup "
                f"returns the pooled record and not the pairing-specific one -- and the "
                f"ambiguity refusal cannot fire, because the code resolves exactly. Ask by "
                f"identifier (definitions_of() lists them), or load one store convention at a "
                f"time.",
                RuntimeWarning, stacklevel=2)

        if self.skempi:
            self._add_role_forms()

        self.malformed = frozenset(
            (p, k) for (p, k) in self._records
            if is_malformed(p, self._records[(p, k)].get("cleaned", k)))
        if self.malformed:
            warnings.warn(
                f"{len(self.malformed)} indexed records come from SKEMPI rows that are malformed "
                f"upstream: the mutation string repeats one substitution where a symmetric "
                f"partner belongs, so they are 3-substitution mutants labelled as "
                f"4-substitution. The energies are correct for what was requested and wrong for "
                f"what was meant. They are NOT corrected here, because substituting the intended "
                f"string would produce a record matching no SKEMPI row. Test with "
                f"FoldxLookup.is_malformed(pdb, mutation); see skempi_foldx.MALFORMED_ROWS.",
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
        mismatched: List[Tuple[str, str, str, str]] = []
        for (pdb, key), rec in list(self._records.items()):
            # `pdb` is a SKEMPI identifier, which is what the table is keyed by. A store written
            # under bare codes still resolves, because a code that SKEMPI defines once is the
            # identifier's prefix and the fallback below finds it.
            entry = self.skempi.get(pdb) or self._sole_definition(pdb)
            if entry is None:
                continue
            code = code_of(pdb)
            # Counted before any refusal below, so a complex this method declines to resolve
            # still registers as one whose mapping was wanted. Otherwise exactness is claimed over
            # precisely the complexes that got no role resolution at all.
            if self.mapping_dir is not None and code not in cache:
                cache[code] = load_chain_mapping(self.mapping_dir, code)
                self._mappings_wanted += 1
                if cache[code]:
                    self._mappings_read += 1

            author = rec.get("cleaned")
            if author is None:
                # No `cleaned`, so the key is being *inferred* to be an author form. Only a
                # mutation SKEMPI lists for this complex may be inferred that way: otherwise a
                # role-keyed store has its key taken for an author form and remapped again,
                # aliasing a real record onto a residue it does not describe. Every bundled record
                # carries `cleaned`, so this path is reached only by a foreign store.
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
                parts = [to_role_form(a, groups, cache[code]) for a in author.split(",")]
                role = ",".join(parts) if all(parts) else None
            if role is None:
                role = self._role_by_identity(entry, author)
            # A record that already carries a `role` has the same computation's known-good answer
            # attached to it. If the recompute disagrees, the mapping files this instance was
            # given describe different residue numbering than the ones that built the store --
            # shifted offsets, a mis-parsed file, an older mapping set -- and installing the
            # recomputed name would put an alias on a label belonging to a different mutation,
            # answered confidently and counted as covered. That is the failure `_role_by_identity`
            # was fixed for, reached through the exact path instead of the inferred one. Keep the
            # shipped name, and say so once; with correct mappings the recompute is exactly
            # idempotent, so this cannot fire on good input.
            shipped = rec.get("role")
            if shipped and role and role != shipped:
                mismatched.append((pdb, key, shipped, role))
                role = shipped
            if role and (pdb, role) not in self._records:
                self._alias.setdefault((pdb, role), (pdb, key))

        self.role_mismatches = tuple(mismatched)
        if mismatched:
            warnings.warn(
                f"{len(mismatched)} records carry a `role` that disagrees with what mapping_dir= "
                f"recomputes, e.g. {mismatched[:2]}. The mapping files describe different residue "
                f"numbering than the ones that built this store, so the recomputed names would "
                f"alias labels belonging to other mutations. The shipped `role` is kept. Check "
                f"that mapping_dir= holds the SKEMPI mapping files for this release.",
                RuntimeWarning, stacklevel=2)

    def _sole_definition(self, name: str):
        """The SKEMPI entry for a bare PDB code, when the table defines it exactly once.

        Lets a store written under bare codes -- every store this package shipped before 0.2.0 --
        still resolve role forms against an identifier-keyed table. ``None`` where the code
        carries two definitions, since the pairing decides the role numbering and picking one
        would remap residues against an interface the record may not belong to.
        """
        if self.skempi is None:
            return None
        found = [e for e in self.skempi.values() if e.pdb == name]
        return found[0] if len(found) == 1 else None

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
            # The role convention renumbers each chain by the cumulative residue count of the
            # chains before it in its group (`to_role_form`: offsets[chain] + n). This fallback
            # has no residue counts -- that is what `mapping_dir` supplies -- so the author
            # position is the role position only for the first chain of a group, where the offset
            # is zero. Emitting it for any other chain produces a label in neither convention:
            # the role chain letter with the author number. Such a label names no record, so it
            # would install an alias onto whichever record happened to carry the substitution,
            # answer confidently, and count as covered. Refuse instead; `mapping_dir` resolves it.
            group = entry.group1 if mut.chain in entry.group1 else entry.group2
            if mut.chain != group[0]:
                return None
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
        if self._shipped_roles and not self._unbacked_records:
            forms.append("role-chain form (exact, from the shipped `role` field)")
            return forms
        if self.is_exact:
            forms.append("role-chain form (exact, via residue offsets)")
        elif self.skempi and self._mappings_read:
            forms.append(f"role-chain form (partly exact: {self._mappings_read} of "
                         f"{self._mappings_wanted} mappings read, remainder by identity)")
        elif self.skempi:
            forms.append("role-chain form (identity matching; ambiguous cases unresolved)")
        return forms

    def definitions_of(self, pdb: str) -> List[str]:
        """Every SKEMPI identifier this store holds for a PDB code, sorted.

        One entry for almost every structure. Two for the codes SKEMPI defines under two chain
        pairings, where a lookup by code alone is genuinely ambiguous and this says so.
        """
        return list(self._by_code.get(pdb, ()))

    def is_malformed(self, pdb: str, mutation: str) -> bool:
        """Whether this record's SKEMPI row is malformed upstream.

        The row was wrong before FoldX saw it: the mutation string repeats a substitution where a
        symmetric partner belongs, so the energies are right for the three substitutions named and
        wrong for the four the label claims. See :data:`skempi_foldx.MALFORMED_ROWS`.
        """
        found = self._resolve(pdb, mutation)
        return found in self.malformed if found else False

    @property
    def is_exact(self) -> bool:
        """Whether role-chain labels resolve by residue offset rather than by identity.

        Two ways to be exact. A store whose records ship a ``role`` field already is: that name was
        computed by residue offset when the store was built, which is the same computation
        ``mapping_dir`` would repeat. Otherwise it takes a mapping for *every* complex that needed
        one -- a directory holding some of them resolves the rest by identity, and reporting that
        as exact would suppress the lower-bound warning over precisely the rows still estimated.

        A record with no ``role`` needs none: for all 2220 of them the role form *is* the key,
        because the chain already sits first in its group. Since the fixed-width split in
        :func:`load_chain_mapping` was corrected, no record in the bundled store is left without an
        exact role name -- the 16 that once had none were all chains parsed short or lost. So
        exactness stays a claim about the method rather than a promise that every row of an
        arbitrary store resolves, but the bundled store now meets it row by row as well.

        **The role-backed test is per store, not per record**, and cannot be otherwise: an absent
        ``role`` is not evidence of a gap, for the reason just given. A store that carries role
        names is taken to have been built by something that computes them. A hand-assembled store
        with a single role field among thousands will therefore be reported exact when it is not
        -- the one shape this cannot detect, and no threshold that would catch it is defensible.
        Indexing such a store alongside the bundled one is the case to avoid.
        """
        if self._shipped_roles and not self._unbacked_records:
            return True
        return bool(self.skempi and self._mappings_wanted
                    and self._mappings_read == self._mappings_wanted)

    def _candidates(self, name: str, mutation: str) -> List[Tuple[str, str]]:
        """Every record ``(name, mutation)`` could mean, where ``name`` may be a bare PDB code.

        The store is keyed by SKEMPI identifier, and a caller holding only a code is asking a
        question with one answer for almost every structure and two for the three SKEMPI defines
        twice. Enumerating rather than picking is what keeps that second answer from being chosen
        silently -- which is the defect this store layout exists to remove.
        """
        out = []
        for ident in self._by_code.get(name, ()):
            key = (ident, mutation)
            hit = key if key in self._records else self._alias.get(key)
            if hit and hit not in out:
                out.append(hit)
        return out

    def _resolve(self, pdb: str, mutation: str) -> Optional[Tuple[str, str]]:
        """The one record named, or ``None`` -- including when more than one is named.

        An identifier resolves directly. A bare code resolves only when its definitions agree on a
        single record; two candidates mean the caller has not said which interface it wants, and
        returning either would answer a question that was never asked.
        """
        key = (pdb, mutation)
        if key in self._records:
            return key
        hit = self._alias.get(key)
        if hit:
            return hit
        found = self._candidates(pdb, mutation)
        return found[0] if len(found) == 1 else None

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
        # way into a feature vector, where they are far harder to trace back. term_vector() makes
        # the same refusal; this one names the record, which is what a caller needs to find it.
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
        the remainder absent is what turns a nearly complete store into a two-thirds figure.
        """
        rec = self.get(pdb, mutation)
        if rec is not None:
            return rec
        ambiguous = self._candidates(pdb, mutation)
        if len(ambiguous) > 1:
            raise KeyError(
                f"{pdb} names {len(ambiguous)} interface definitions holding {mutation!r} "
                f"({', '.join(i for i, _ in ambiguous)}). SKEMPI defines this code under more "
                f"than one chain pairing, and the two are different measurements of the same "
                f"substitution. Ask for the identifier, not the code.")
        others = [k for (p, k) in self._records if p == pdb or code_of(p) == pdb]
        if not others:
            # Diagnose before declaring absence: a case-only mismatch and an unloaded arm both
            # present as "not in the store", and both send a reader looking in the wrong place.
            # Compared on the code, since `pdb` may be either a code or an identifier and the
            # records are keyed by identifier: an identifier-to-code comparison never matches, so
            # a lower-cased code would be reported as absent rather than as mis-cased.
            cased = [p for p in self._by_code if p.upper() == code_of(pdb).upper()]
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
        if self._shipped_roles and not self._unbacked_records:
            hint = (" This instance already resolves role-chain labels exactly, from the shipped "
                    "`role` field, so the name is absent rather than filed differently.")
        elif self.skempi is None:
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
        of genuinely uncomputed rows is not. The bundled store resolves exactly on its own, so it
        never raises this; a store with no ``role`` field is what it is for.

        **The warning describes this lookup, not the rows.** Why a particular row missed is not
        recoverable from the label: any test available here -- the chain letters, the arm's shape,
        the code's case -- is a property of the label, and the label does not carry the answer. So
        the warning states only what this instance can match and what would widen it.
        :meth:`require` diagnoses a single row, where the record itself is available to check
        against.

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
