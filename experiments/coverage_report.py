#!/usr/bin/env python3
"""FoldX coverage of an evaluation set.

Answers the question worth asking before committing to a setup: of the rows in this split or
benchmark subset, how many carry a real FoldX energy, and what is the shortfall made of.

The package ships no splits, so paths are supplied. Two input shapes are accepted:

    --split  <dir>    a directory of fold_<k>/<base>_{train,val,test}.tsv
    --table  <file>   a flat TSV

Both are read as headerless TSV where column 0 is a sequence id starting with the PDB code
(``1AK4.A.D_A`` and ``1A22_A`` are both understood) and column 2 is the mutation. ``--pdb-col``
and ``--mut-col`` override that.

    python experiments/coverage_report.py --split path/to/splits_cath_kfold \\
        --skempi-csv skempi_v2.csv --mapping-dir PDBs/

Split files name mutations by role chain -- the partners are ``A``/``B`` -- while every store
record is keyed by author chain, and neither string matches the other.

**The shipped store needs neither ``--skempi-csv`` nor ``--mapping-dir``**: each record carries its
own role-chain name, so the figures are already exact. Measured on a 1636-row multi-point split,
every row resolves with no arguments at all, and supplying both changes nothing.

**Pass them for a store this package did not build**, whose records carry no ``role`` field. There
the figures are a LOWER BOUND without them, and say so.

All naming resolution lives in :class:`skempi_foldx.FoldxLookup`. This script reads rows, asks it,
and formats the answer; it deliberately holds no resolution logic of its own, because two copies
of that logic is how the figures diverged in the first place.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skempi_foldx import FoldxLookup  # noqa: E402


def sequence_id_to_name(field: str) -> str:
    """The store name a split's column 0 refers to: ``1A22.A.B_A`` -> ``1A22_A_B``.

    Splits write an interface as ``<code>.<group1>.<group2>``, then a partner suffix. That is the
    same thing a store file is named, in different punctuation, so keeping it costs one
    substitution and buys the disambiguation the 0.2.0 layout exists for: reducing the field to
    its PDB code instead makes ``3SE4`` and ``3SE3`` -- the codes SKEMPI defines under two
    pairings -- unresolvable, and their rows are then reported uncovered when the split said
    exactly which interface it meant.

    Everything before the FIRST underscore, not the last: an augmented split's reverse-mutation
    labels carry a further underscore-joined suffix (``1CSE.E.I_E_LI38S-GI32Y``), and cutting at
    the last would leave ``1CSE.E.I_E``.

    A label with no groups (``1A22_A``, the S1102 form) yields the bare code, which
    :class:`FoldxLookup` resolves wherever it is unambiguous.
    """
    # Store keys are upper case; split files are not reliably so, and a whole-file case
    # difference would otherwise halve the figure this script exists to compute.
    return field.split("_", 1)[0].replace(".", "_").upper()


def read_rows(path: Path, pdb_col: int, mut_col: int):
    """Yield ``(name, mutation)`` from a headerless TSV."""
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) <= max(pdb_col, mut_col) or not row[pdb_col]:
                continue
            yield sequence_id_to_name(row[pdb_col]), row[mut_col]


def collect(target: Path, pdb_col: int, mut_col: int):
    """Distinct rows under a split directory, or in a flat table, in first-seen order."""
    files = sorted(target.rglob("*.tsv")) if target.is_dir() else [target]
    rows, seen = [], set()
    for f in files:
        for row in read_rows(f, pdb_col, mut_col):
            if row not in seen:
                seen.add(row)
                rows.append(row)
    return rows, len(files)


def report(label: str, rows, lookup: FoldxLookup, exact: bool) -> dict:
    covered, total, missing = lookup.coverage(rows)
    arms = {"sp": 0, "mp": 0}
    for row in rows:
        arm = lookup.arm_of(*row)
        if arm:
            arms[arm] += 1
    by_complex: dict = {}
    for pdb, _ in missing:
        by_complex[pdb] = by_complex.get(pdb, 0) + 1
    worst = sorted(by_complex.items(), key=lambda kv: -kv[1])[:5]
    return {
        "label": label,
        "rows": total,
        "covered": covered,
        "pct": round(100.0 * covered / total, 1) if total else 0.0,
        "single_point": arms["sp"],
        "multi_point": arms["mp"],
        "uncovered": len(missing),
        "uncovered_complexes": len(by_complex),
        "worst_complexes": [{"pdb": p, "rows": n} for p, n in worst],
        "exact": exact,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

    # --split and --table share one ordered destination. Collecting them separately and
    # concatenating reorders the inputs relative to --label, which silently mislabels every row
    # once the two are interleaved on the command line.
    class _Input(argparse.Action):
        def __call__(self, parser, ns, value, option_string=None):
            ns.inputs = (getattr(ns, "inputs", None) or []) + [(option_string.lstrip("-"), value)]

    ap.add_argument("--split", action=_Input, dest="inputs", help="split directory (repeatable)")
    ap.add_argument("--table", action=_Input, dest="inputs", help="flat TSV (repeatable)")
    ap.add_argument("--label", action="append", default=[], help="label per input, in order")
    ap.add_argument("--pdb-col", type=int, default=0)
    ap.add_argument("--mut-col", type=int, default=2)
    ap.add_argument("--skempi-csv", help="SKEMPI 2.0 table; needed to resolve role-chain labels")
    ap.add_argument("--mapping-dir", help="directory of SKEMPI <code>.mapping files; makes the "
                                          "role-chain resolution exact rather than a lower bound")
    ap.add_argument("--json", help="also write the report as JSON")
    args = ap.parse_args()

    pairs = getattr(args, "inputs", None) or []
    if not pairs:
        ap.error("give at least one --split or --table")
    targets = [Path(v) for _, v in pairs]
    if args.label and len(args.label) != len(targets):
        ap.error(f"{len(args.label)} --label values for {len(targets)} inputs; "
                 f"give one per input or none")
    labels = list(args.label) or [t.name for t in targets]

    lookup = FoldxLookup(skempi_csv=args.skempi_csv, mapping_dir=args.mapping_dir)
    # Asked of the lookup, not inferred from the flags: a --mapping-dir that is empty or misspelt
    # resolves nothing, and the report must not call that exact.
    exact = lookup.is_exact
    print(f"resolving: {', '.join(lookup.conventions)}")
    if not exact:
        if lookup.skempi is None:
            print("figures are a LOWER BOUND for role-chain labels -- pass --skempi-csv and "
                  "--mapping-dir to make them exact")
        elif lookup._mappings_read:
            print(f"figures are a LOWER BOUND for role-chain labels in the complexes whose "
                  f"mapping was missing: {lookup._mappings_read} of "
                  f"{lookup._mappings_wanted} mappings were read")
        else:
            print("figures are a LOWER BOUND for role-chain labels -- no mapping file was read, "
                  "so role forms fall back to identity matching; check --mapping-dir")
    print()

    print(f"{'evaluation set':44s} {'rows':>6s} {'covered':>8s} {'':>6s} {'SP':>6s} {'MP':>6s}")
    print("-" * 82)
    out = []
    for target, label in zip(targets, labels):
        if not target.exists():
            print(f"{label:44s} {'MISSING':>6s}")
            continue
        rows, n_files = collect(target, args.pdb_col, args.mut_col)
        r = report(label, rows, lookup, exact)
        r["files_read"] = n_files
        out.append(r)
        print(f"{label:44s} {r['rows']:6d} {r['covered']:8d} {r['pct']:5.1f}% "
              f"{r['single_point']:6d} {r['multi_point']:6d}")
        if r["uncovered"]:
            worst = ", ".join(f"{w['pdb']}({w['rows']})" for w in r["worst_complexes"])
            print(f"{'':44s} uncovered {r['uncovered']} over "
                  f"{r['uncovered_complexes']} complexes: {worst}")
    if args.json:
        Path(args.json).write_text(json.dumps(out, indent=1) + "\n")
        print(f"\nwrote {args.json}")


if __name__ == "__main__":
    main()
