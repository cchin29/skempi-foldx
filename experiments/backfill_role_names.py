#!/usr/bin/env python3
"""Add the 16 ``role`` names a mapping-file parsing bug kept out of the shipped store.

``build_per_definition_store.py`` writes each record's role-chain name through
``lookup.to_role_form``, which reads chain lengths from the SKEMPI ``.mapping`` files. Those files
are fixed-width, so an author number of 1000 or more runs into the chain letter and the line
arrives as three whitespace fields instead of four. ``load_chain_mapping`` dropped such lines, so
``2NYY``/``2NZ9`` chain A measured 971 residues instead of 1267 and its high-numbered mutations
fell outside the range check, while ``1S0W``/``1XXM`` lost the whole of chain C and produced no
role name at all. ``load_chain_mapping`` now splits the fused field; this script brings the store
that was built before it into line.

Rebuilding the store from the campaigns would also fix it, and is the wrong instrument: the inputs
are CPU-weeks of FoldX and no value may move. This adds the missing field and nothing else.

**It needs no chain mappings**, which is why it can run from a bare checkout -- the package
redistributes neither them nor ``skempi_v2.csv``. Every one of the 16 is a chain sitting FIRST in
its group, where the offset is 0 and the role form follows from the identifier alone::

    1S0W_A_C   AC142F  ->  AB142F        chain C is group 2, and its first chain

That is a narrower rule than ``to_role_form``, and the two are held to each other by
``test_backfilled_roles_match_to_role_form`` wherever both apply.

Idempotent: a second run reports 0 changes. ``--check`` exits non-zero if any record is missing a
name this rule would give it *or* carries one that disagrees with it, and writes nothing; so does a
run against a directory holding no store, which is otherwise indistinguishable from a clean one.
Nothing runs it automatically -- it is a manual audit, not a CI gate; what CI holds is the store
itself, via ``test_shipped_store_matches_the_backfill_rule``.

    python experiments/backfill_role_names.py --check
    python experiments/backfill_role_names.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from skempi_foldx.store import DATA_DIR, MULTI_KEY, SINGLE_KEY  # noqa: E402


def role_of(mutation: str, groups) -> str | None:
    """The role form of a mutation whose chain is first in its group, else ``None``.

    Deliberately NOT ``to_role_form``: that one takes chain mappings this script does not have,
    and the case it covers here is exactly the one needing none -- a first chain has offset 0, so
    no length enters the arithmetic. It is the narrower function: it declines every later chain,
    and it applies no upper bound because it has nothing to bound against.
    """
    if len(mutation) < 4:
        return None
    chain, g1, g2 = mutation[1], *groups
    if chain in g1:
        side, group = "A", g1
    elif chain in g2:
        side, group = "B", g2
    else:
        return None
    if group[0] != chain or not mutation[2:-1].isdigit():
        return None
    return f"{mutation[0]}{side}{mutation[2:-1]}{mutation[-1]}"


def role_for_record(key: str, record: dict, groups) -> str | None:
    """A whole record's role name, or ``None`` if any of its substitutions declines.

    All or nothing, matching ``build_per_definition_store.role_form``: a multi-point name with one
    substitution left in author form names no residue set at all.
    """
    parts = [role_of(m, groups) for m in record.get("cleaned", key).split(",")]
    return ",".join(parts) if all(parts) else None


def scan(data_dir: Path):
    """``(pending, wrong)`` -- names this rule would add, and stored ones it contradicts.

    Absence is not the only way a role name can be defective, and checking only for absence is how
    a bad one would ship: the audit passed on this store while a role read `DB9999A`. So every
    record the rule can decide is compared, not just the empty ones. `wrong` is never written --
    a stored name that disagrees with the rule is a fault to look at, not one to overwrite.
    """
    store = sorted(data_dir.glob("results_*/*.json"))
    if not store:
        # An audit that cannot fail on an unreadable target reports success for a typo in
        # `--data-dir`, which is the one input most likely to be wrong.
        raise SystemExit(f"{data_dir}: no results_*/*.json under here -- nothing was checked.")
    pending, wrong = [], []
    for path in store:
        ident = path.stem.split("_")
        if len(ident) != 3:
            # A store keyed by bare PDB code carries no pairing, and the role side is undecidable
            # without one. Skipping is right: this script only ever adds what is already implied.
            continue
        groups = (ident[1], ident[2])
        payload = json.loads(path.read_text())
        records = payload.get(SINGLE_KEY) or payload.get(MULTI_KEY) or {}
        for key, record in records.items():
            role = role_for_record(key, record, groups)
            if role is None:
                continue                 # a later chain: this rule has no opinion either way
            stored = record.get("role")
            if stored is None:
                if role != key:
                    pending.append((path, key, role))
            elif stored != role:
                wrong.append((path, key, stored, role))
    return pending, wrong


def apply(pending) -> None:
    """Append ``role`` to each named record, preserving every other byte's order and value."""
    by_file: dict[Path, dict[str, str]] = {}
    for path, key, role in pending:
        by_file.setdefault(path, {})[key] = role
    for path, roles in by_file.items():
        payload = json.loads(path.read_text())
        which = SINGLE_KEY if payload.get(SINGLE_KEY) else MULTI_KEY
        for key, role in roles.items():
            # Appended last, which is where build_per_definition_store puts it: that script sets
            # `_source` and then `role` on a record that already carries `cleaned`, so every one
            # of the 3871 records built with a role name orders them `cleaned, _source, role`.
            # Inserting after `cleaned` instead would leave these 16 the only records in the
            # store with a different field order -- invisible in a value diff, but it would make
            # a later rebuild move the field and so produce a diff that is not only role lines.
            payload[which][key]["role"] = role
        # `indent=1`, no trailing newline: byte-for-byte what write_store and
        # build_per_definition_store emit. A stray "\n" here would show as a change on every
        # touched file, in a diff whose whole point is that only `role` lines move.
        path.write_text(json.dumps(payload, indent=1))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data-dir", type=Path, default=DATA_DIR,
                    help="store root holding results_sp/ and results_mp/ (default: the shipped)")
    ap.add_argument("--check", action="store_true",
                    help="report and exit 1 if anything is missing or disagrees; write nothing")
    args = ap.parse_args()

    pending, wrong = scan(args.data_dir)
    for path, key, stored, role in wrong:
        print(f"  {path.parent.name}/{path.stem:16s} {key:34s} stored {stored} != {role}")
    for path, key, role in pending:
        print(f"  {path.parent.name}/{path.stem:16s} {key:38s} -> {role}")

    if wrong:
        # Never written over, whichever mode: the rule is narrow enough that a disagreement is at
        # least as likely to mean the rule is wrong here as that the store is.
        print(f"\n{len(wrong)} stored role names disagree with this rule. Not written; look at "
              f"them.", file=sys.stderr)
        return 1
    if not pending:
        print(f"{args.data_dir}: every first-chain record carries the role name this rule gives "
              f"it.")
        return 0
    if args.check:
        print(f"\n{len(pending)} records are missing a role name this rule would give them.",
              file=sys.stderr)
        return 1
    apply(pending)
    print(f"\nAdded {len(pending)} role names.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
