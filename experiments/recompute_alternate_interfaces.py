#!/usr/bin/env python3
"""Recomputation of the complexes SKEMPI defines under two interfaces.

SKEMPI keys a row by ``<pdb>_<group1>_<group2>``, and three codes appear with two different
partner pairings. :func:`~skempi_foldx.load_skempi` keys by that identifier, so each pairing is
its own entry — but a *campaign* still computes one pairing at a time, and this resolves which
campaigns are needed and what each must submit.

Before 0.2.0 the parser keyed on the code alone: one pairing won and every mutation was pooled
under it, including the ones belonging to the other, which were then scored against an interface
they are not part of. This script is what recomputed them.

Two campaigns close the *mis-pairing*. Each recomputes a discarded definition's mutations against
their own pairing:

    2C5D   AB,CD   19 multi-point variants   (kept pairing: A,C)
    3SE4   B,A     18 single-point mutations (kept pairing: B,C)

3SE3 is defined twice as well and contributes nothing to that repair: every mutation of its second
definition also appears in the first, so none was scored against a pairing SKEMPI does not
associate with it.

``--all-definitions`` asks a wider question. Repairing the mis-paired records leaves the store one
record per *mutation*; computing every definition as its own pairing leaves it one record per
*(definition, mutation)*, which is what the SKEMPI identifier ``<pdb>_<group1>_<group2>`` denotes
and what the coverage denominators count. The two differ on 3SE3's ``B,C`` pair: its 2 multi-point
variants hold a value computed against ``B,A``, and ``B,C`` is a different interface. Nothing marks
that record as belonging to another pairing, because under the narrower question it does not.

Two properties of the pipeline shape how this must run.

**The whole list, not the missing part.** A mutation's energy depends on every entry preceding it
in ``individual_list.txt``, so the 12 wrong records in 3SE4 cannot be repaired by recomputing 12
mutations: the run submits all 18 that the ``B_A`` definition carries, and the 6 it shares with
``B_C`` are recomputed too, under the pairing being rebuilt.

**One file per complex.** A campaign writes ``<pdb>.json`` keyed by mutation, so a second pairing
for the same code cannot share a results directory with the first — for 3SE4 six mutation keys
occur in both, and the later write would silently replace them. Each pairing therefore gets its
own directory, named for the interface that produced it.

    python experiments/recompute_alternate_interfaces.py --skempi-csv skempi_v2.csv --plan
    python experiments/recompute_alternate_interfaces.py --skempi-csv skempi_v2.csv \\
        --out results_alt_interfaces --pdb-dir PDBs/ --jobs 4

``--plan`` resolves and prints the campaigns without invoking FoldX, which is the cheap way to
confirm the worklist before committing a run.

**`--seed-repairs` saves the repair, but the seed must be a 1x structure.** Measured 2026-08-01:
`RepairPDB` is deterministic given the same binary and input, so a fresh repair of `2C5D` and
`3SE4` came out byte-identical to the repair sweep's `repair_round_1.pdb`. Seeding is therefore
safe and saves roughly 20 minutes of the ~37 these two complexes took.

The trap is which file to seed from. The sweep's `chain_work/<PDB>/<PDB>_Repair.pdb` is the
**round-4** structure, not round 1 -- seeding from it computes these records at 4x repair while
the shipped store is 1x, and nothing downstream would flag the mismatch. Use
`round_1/.../repair_round_1.pdb`, or omit the flag and pay the repair.
"""

from __future__ import annotations

import argparse
import csv
import collections
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skempi_foldx import FoldxConfig, SkempiComplex, run_campaign  # noqa: E402
from skempi_foldx.run import MODE_AUTHOR, MODE_VARIANT            # noqa: E402
from skempi_foldx.store import check_identifier                  # noqa: E402


def definitions(path: Path):
    """``{pdb: {(group1, group2): {cleaned, ...}}}``, preserving the pairing of every row.

    Read straight from the table rather than through :func:`~skempi_foldx.load_skempi`, because
    this needs the *order* the rows appear in: the first pairing seen per code is the one 0.1.0's
    collapse kept, which is what makes a record "already computed against its own interface" or
    not. The package parser keys by identifier and no longer records that ordering.
    """
    rows = collections.defaultdict(lambda: collections.defaultdict(set))
    first = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=";")
        next(reader, None)
        for row in reader:
            if not row or not row[0]:
                continue
            parts = row[0].split("_")
            if len(parts) < 3:
                continue
            pdb, groups = parts[0], (parts[1], parts[2])
            first.setdefault(pdb, groups)
            rows[pdb][groups].add(row[2])
    return rows, first


def campaigns(path: Path, include_kept: bool = False, all_definitions: bool = False):
    """One entry per definition to compute, for each complex carrying more than one.

    A mutation appearing under both pairings does not by itself trigger a campaign — it was
    already computed against a pairing SKEMPI associates with it — but it is still submitted as
    part of the list, because the list is what determines the values.

    ``include_kept`` adds the surviving definition's own list. That is the difference between
    repairing the mis-paired records and making the complex internally consistent: the kept
    definition's values were themselves computed in a list holding both definitions' mutations,
    so they belong to no single pairing's list either. Recomputing only the discarded definition
    leaves one complex spanning two list compositions.

    ``all_definitions`` computes every definition of every multiply-defined complex, whether or
    not it holds a mutation the kept one lacks. It implies ``include_kept``, and it is the mode
    that matches a per-definition store: the unit of the store becomes the SKEMPI identifier
    rather than the PDB code, so a definition whose mutations are all shared still needs its own
    values, computed against its own interface.
    """
    rows, first = definitions(path)
    out = []
    for pdb, defs in sorted(rows.items()):
        if len(defs) < 2:
            continue
        kept = defs[first[pdb]]
        if not all_definitions and not any(muts - kept for g, muts in defs.items()
                                           if g != first[pdb]):
            continue          # every alternate mutation is already in the kept list
        for groups, muts in sorted(defs.items()):
            if groups == first[pdb]:
                if not (include_kept or all_definitions):
                    continue
            elif not all_definitions and not (muts - kept):
                continue
            multi = {m for m in muts if "," in m}
            single = muts - multi
            for label, subset, mode in (("mp", multi, MODE_VARIANT),
                                        ("sp", single, MODE_AUTHOR)):
                if subset:
                    out.append({
                        "pdb": pdb, "groups": groups, "arm": label, "mode": mode,
                        "mutations": sorted(subset),
                        "corrects": len({m for m in subset if m not in kept}),
                    })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skempi-csv", required=True, help="SKEMPI 2.0 table")
    ap.add_argument("--out", help="parent directory for the per-interface results directories")
    ap.add_argument("--plan", action="store_true", help="resolve and print, without running FoldX")
    ap.add_argument("--include-kept", action="store_true",
                    help="also recompute the surviving definition's own list, so each pairing's "
                         "values come from one list rather than from the pooled one. An arm the "
                         "other definition never touches is unaffected and recomputes identically")
    ap.add_argument("--all-definitions", action="store_true",
                    help="compute every definition of every multiply-defined complex, including "
                         "one whose mutations are all shared with the kept definition. Implies "
                         "--include-kept. Use when the store is keyed by SKEMPI identifier rather "
                         "than by PDB code, so each definition needs values of its own")
    ap.add_argument("--seed-repairs", action="append", default=[], metavar="DIR",
                    help="campaign work directory holding an existing <pdb>/<pdb>_Repair.pdb "
                         "(repeatable); reuses the repair the shipped values were built on. "
                         "MUST be a 1x repair -- see the note in the module docstring")
    ap.add_argument("--pdb-dir", help="directory of SKEMPI cleaned PDB files")
    ap.add_argument("--binary", help="FoldX 5 executable")
    ap.add_argument("--jobs", type=int, default=1)
    args = ap.parse_args()

    plan = campaigns(Path(args.skempi_csv), include_kept=args.include_kept,
                     all_definitions=args.all_definitions)
    if not plan:
        print("no complex carries a discarded definition with mutations of its own")
        return

    for c in plan:
        g = ",".join(c["groups"])
        print(f"{c['pdb']}  --analyseComplexChains {g:6s} [{c['arm']}]  "
              f"{len(c['mutations']):3d} mutations, {c['corrects']} of them currently "
              f"computed against the other interface")
    if args.plan:
        return
    if not args.out:
        ap.error("--out is required unless --plan is given")
    if not args.pdb_dir:
        ap.error("--pdb-dir is required to run; FoldX needs the structures")

    for c in plan:
        # Named for the interface, not the complex: two pairings of one code write the same
        # filename, and sharing a directory would let the later run replace the earlier one.
        results = Path(args.out) / check_identifier(
            f"{c['pdb']}_{c['groups'][0]}_{c['groups'][1]}_{c['arm']}")
        config = FoldxConfig(
            results_dir=results,
            # One work directory per pairing as well: the mutant structures BuildModel writes are
            # named by index within the list, so two pairings of one complex would overwrite each
            # other's intermediates exactly as they would overwrite each other's results.
            work_dir=Path(args.out) / "work" / results.name,
            pdb_dir=Path(args.pdb_dir),
            skempi_csv=Path(args.skempi_csv),
            binary=Path(args.binary) if args.binary else None,
            repair_seed_dirs=tuple(Path(s) for s in args.seed_repairs),
        )
        entry = SkempiComplex(c["pdb"], c["groups"][0], c["groups"][1])
        print(f"\n[run] {c['pdb']} {entry.groups} -> {results}")
        status = run_campaign({c["pdb"]: c["mutations"]}, {c["pdb"]: entry},
                              config, mode=c["mode"], jobs=args.jobs)
        print(f"[run] {status}")


if __name__ == "__main__":
    main()
