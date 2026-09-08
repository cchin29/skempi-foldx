#!/usr/bin/env python3
"""Assembly of the shipped store, one record per SKEMPI interface definition.

SKEMPI names a row ``<pdb>_<group1>_<group2>``. Three codes appear under two such names, so the
identifier and the PDB code are not the same thing, and a store keyed by code has to discard one
pairing per multiply-defined complex. This builds the store keyed by the identifier instead, which
is what the coverage denominators count and what makes a value explainable: the pairing is half of
what ``AnalyseComplex`` was asked.

Two kinds of input, and the distinction is the whole of the logic.

**A complex with one definition** is read from the repair sweep's ``round_1``, whose worklist is
each complex's *union* mutation list. That is the property the release is for: a value computed in
a per-benchmark subset is not interchangeable with one computed in the union, and 2494 of the
0.1.0 single-point values came from a subset campaign.

**A complex with more than one** is read from a per-definition campaign, one per
``(pdb, groups, arm)``, produced by ``recompute_alternate_interfaces.py --all-definitions``.
``round_1`` cannot serve these: its worklist pools both definitions' mutations into one list, so
its values belong to no single pairing's list, and for the mutations of a discarded definition it
also scored them against the wrong interface.

Both inputs are one ``RepairPDB`` pass at FoldX 5.1, which is what the shipped store has always
been. The per-definition campaigns were seeded from the repaired structures the earlier campaigns
built, verified byte-identical beforehand, so repair count is held fixed across the two inputs and
is not a variable in the assembled store.

The assembled store is checked against ``skempi_v2.csv`` rather than against the inputs: every
definition's key set must equal the mutation set SKEMPI lists under that identifier, less anything
the pipeline rejected on validation. A shortfall is reported per definition rather than summed, so
a missing complex cannot hide inside a total.

    python experiments/build_per_definition_store.py --skempi-csv skempi_v2.csv \\
        --round1 <sweep>/round_1 --definitions <defs> --out skempi_foldx/data --dry-run
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import shutil
import sys
from pathlib import Path
from typing import Optional

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skempi_foldx import is_excluded  # noqa: E402
from skempi_foldx.lookup import load_chain_mapping, to_role_form  # noqa: E402
from skempi_foldx.store import check_identifier  # noqa: E402
from skempi_foldx.store import sanitise_meta  # noqa: E402

SINGLE_KEY = "muts"
MULTI_KEY = "variants"


def role_form(cleaned: str, groups, chains) -> str:
    """A record's role-chain name (``LB38D``), or ``""`` when it cannot be derived exactly.

    Shipped as a field because the store is keyed by one convention and split files use the other.
    Deriving it at read time needs ``skempi_v2.csv`` or the ``.mapping`` files, and this package
    redistributes neither -- so a consumer joining role-labelled rows would be left with the one
    lookup they cannot perform. Computing it once, here, is what keeps that join working with
    nothing but the installed package.

    Empty rather than approximate: an insertion code or an out-of-range residue has no exact role
    name, and a guess would alias a record onto a residue it does not describe.
    """
    if not chains:
        return ""
    parts = [to_role_form(part, groups, chains) for part in cleaned.split(",")]
    return ",".join(parts) if all(parts) else ""


def skempi_definitions(path: Path):
    """``{(pdb, groups): {"sp": {...}, "mp": {...}}}`` straight from the table.

    The authority for what the store should contain. Read here rather than through
    :func:`~skempi_foldx.load_skempi` on purpose: that function answers "which mutations belong to
    this complex", and the question here is "which belong to this *identifier*", which is the
    distinction the release exists to restore.
    """
    out = collections.defaultdict(lambda: {"sp": set(), "mp": set()})
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=";")
        next(reader, None)
        for row in reader:
            if not row or not row[0]:
                continue
            parts = row[0].split("_")
            if len(parts) < 3:
                continue
            pdb, groups, cleaned = parts[0], (parts[1], parts[2]), row[2]
            out[(pdb, groups)]["mp" if "," in cleaned else "sp"].add(cleaned)
    return dict(out)


def identifier(pdb: str, groups) -> str:
    return f"{pdb}_{groups[0]}_{groups[1]}"


def read_payload(path: Path):
    """``(records, meta, arm)`` from a campaign's ``<pdb>.json``, whichever key it uses."""
    payload = json.loads(path.read_text())
    if SINGLE_KEY in payload:
        return payload[SINGLE_KEY], payload.get("meta", {}), "sp"
    return payload.get(MULTI_KEY, {}), payload.get("meta", {}), "mp"


def index_definition_campaigns(root: Path):
    """``{(pdb, groups, arm): path}`` over the per-definition campaign directories.

    Directory names are ``<pdb>_<g1>_<g2>_<arm>``, written that way because two pairings of one
    code produce the same ``<pdb>.json`` filename and would otherwise overwrite each other. The
    pairing is re-read from each file's ``meta.groups`` rather than parsed back out of the
    directory name, so a mislabelled directory is caught instead of trusted.
    """
    found = {}
    for results in sorted(p for p in root.iterdir() if p.is_dir() and p.name != "work"):
        arm = results.name.rsplit("_", 1)[-1]
        for path in sorted(results.glob("*.json")):
            records, meta, kind = read_payload(path)
            if not records:
                continue
            groups = tuple(meta["groups"].split(","))
            if arm != kind:
                raise ValueError(f"{results.name} holds {kind} records; its name says {arm}")
            key = (meta["pdb"], groups, kind)
            if key in found:
                raise ValueError(f"two campaigns for {identifier(meta['pdb'], groups)} [{kind}]")
            found[key] = path
    return found


def build(skempi_csv: Path, round1: Path, definitions_root: Path,
          mapping_dir: Optional[Path] = None, source_label: str = "sweep_4x_round_1"):
    """``(store, report)`` — ``store`` maps ``(identifier, arm)`` to ``(records, meta)``."""
    wanted = skempi_definitions(skempi_csv)
    per_code = collections.Counter(pdb for pdb, _ in wanted)
    campaigns = index_definition_campaigns(definitions_root)

    store, report = {}, []
    for (pdb, groups), arms in sorted(wanted.items()):
        excluded = is_excluded(pdb)
        for arm in ("sp", "mp"):
            expected = arms[arm]
            if not expected:
                continue
            ident = identifier(pdb, groups)
            if excluded:
                report.append({"id": ident, "arm": arm, "expected": len(expected), "got": 0,
                               "source": None, "note": f"excluded: {excluded.reason}"})
                continue
            path = campaigns.get((pdb, groups, arm))
            source = "per-definition campaign"
            if path is None:
                if per_code[pdb] > 1:
                    raise ValueError(
                        f"{ident} [{arm}] has no per-definition campaign, and {pdb} carries "
                        f"{per_code[pdb]} definitions -- round_1 pools them, so its values belong "
                        f"to no single pairing. Run recompute_alternate_interfaces.py "
                        f"--all-definitions.")
                path = round1 / ("results" if arm == "sp" else "results_mp") / f"{pdb}.json"
                source = "round_1 union list"
                if not path.exists():
                    report.append({"id": ident, "arm": arm, "expected": len(expected), "got": 0,
                                   "source": None, "note": "no input"})
                    continue
            records, meta, _ = read_payload(path)
            records = {k: dict(v) for k, v in records.items()}   # never mutate the input
            if meta.get("groups") and tuple(meta["groups"].split(",")) != groups:
                raise ValueError(f"{ident} [{arm}]: {path} was computed against "
                                 f"{meta['groups']}, not {','.join(groups)}")
            rejected = {r.split("->")[0] for r in meta.get("validation_failed") or []}
            missing = {m for m in expected
                       if m not in records and not any(m == r.get("cleaned")
                                                       for r in records.values())}
            # Provenance per record, not only per file: a consumer holding one value needs to
            # know which campaign's mutation list produced it, and that is what makes the value
            # explainable at all.
            label = source_label if source == "round_1 union list" else path.parent.name
            chains = load_chain_mapping(mapping_dir, pdb) if mapping_dir else None
            for key, rec in records.items():
                rec["_source"] = label
                role = role_form(rec.get("cleaned", key), groups, chains)
                if role and role != key:
                    rec["role"] = role
            store[(ident, arm)] = (records, sanitise_meta(
                dict(meta, pdb=pdb, groups=",".join(groups), skempi_id=ident, source=source)))
            report.append({"id": ident, "arm": arm, "expected": len(expected),
                           "got": len(records), "source": source,
                           "note": (f"{len(missing)} not in store" if missing else "")
                                   or (f"{len(rejected)} rejected on validation"
                                       if rejected else "")})
    return store, report


def write(store, destination: Path):
    """One ``<identifier>.json`` per definition and arm, into ``results_sp`` / ``results_mp``."""
    # Every name is checked before anything is removed. The rmtree below is not recoverable, so
    # a rejected identifier partway through would leave the previous store deleted and the new
    # one half-written -- the validation would have cost the data it exists to protect.
    for ident, _ in store:
        check_identifier(ident)
    for arm, sub in (("sp", "results_sp"), ("mp", "results_mp")):
        target = destination / sub
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
    for (ident, arm), (records, meta) in sorted(store.items()):
        sub = "results_sp" if arm == "sp" else "results_mp"
        key = SINGLE_KEY if arm == "sp" else MULTI_KEY
        meta = dict(meta, n_mutations=len(records))
        (destination / sub / f"{ident}.json").write_text(
            json.dumps({key: records, "meta": meta}, indent=1))


def write_report(store, report, destination: Path) -> str:
    """``CONSOLIDATION_REPORT.txt`` — what produced each value, and what is absent.

    The 0.1.0 report answered a question this store does not raise: five campaigns had overlapping
    scopes, so it had to say which one won each mutation and where the losers disagreed. Here each
    definition has exactly one source by construction, so there is nothing to arbitrate. What is
    still worth recording is provenance -- a value is a property of the mutation list it was
    computed in, so the campaign that computed it is part of being able to explain it -- and the
    shortfall, stated per definition so it cannot hide inside a total.
    """
    lines = [
        "Source labels name the campaign that computed a value. They record PROVENANCE, not",
        "membership: a value's label says which mutation list produced it, which is what makes",
        "the number explainable, and says nothing about which benchmarks contain the mutation.",
        "",
        "Each record has exactly one source. A definition's mutations are computed together, in",
        "one list, against one interface -- so unlike 0.1.0 there is no precedence order here and",
        "no conflicts to arbitrate. The unit is the SKEMPI identifier <pdb>_<group1>_<group2>.",
        "",
    ]
    for arm, title in (("sp", "Single-point store"), ("mp", "Multi-point store")):
        entries = {i: (r, m) for (i, a), (r, m) in store.items() if a == arm}
        n = sum(len(r) for r, _ in entries.values())
        codes = {i.split("_")[0] for i in entries}
        lines += [title, "=" * 40,
                  f"{len(entries)} interface definitions over {len(codes)} PDB codes / "
                  f"{n} mutations", "", "contribution by source:"]
        by_src = collections.Counter()
        for recs, _ in entries.values():
            for rec in recs.values():
                by_src[rec.get("_source", "?")] += 1
        for label, count in sorted(by_src.items(), key=lambda kv: (-kv[1], kv[0])):
            lines.append(f"  {label:26s} {count:5d} mutations")
        lines.append("")

    absent = [r for r in report if r["got"] != r["expected"]]
    lines += ["Absent from the store", "=" * 40]
    if not absent:
        lines.append("  nothing: every mutation SKEMPI lists under every definition is present.")
    for r in absent:
        lines.append(f"  {r['id']:14s} [{r['arm']}]  {r['expected'] - r['got']:3d} of "
                     f"{r['expected']:4d} absent -- {r['note']}")
    lines.append("")
    text = "\n".join(lines)
    (destination / "CONSOLIDATION_REPORT.txt").write_text(text)
    return text


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--skempi-csv", required=True)
    ap.add_argument("--round1", required=True, help="sweep round_1 dir holding results/ results_mp/")
    ap.add_argument("--definitions", required=True, help="per-definition campaign parent directory")
    ap.add_argument("--out", required=True, help="package data directory to write")
    ap.add_argument("--mapping-dir", help="SKEMPI .mapping files, for the shipped role-chain name")
    ap.add_argument("--dry-run", action="store_true", help="report without writing")
    args = ap.parse_args()

    store, report = build(Path(args.skempi_csv), Path(args.round1), Path(args.definitions),
                          Path(args.mapping_dir) if args.mapping_dir else None)

    shortfalls = [r for r in report if r["got"] != r["expected"]]
    by_arm = collections.Counter()
    for r in report:
        by_arm[(r["arm"], "expected")] += r["expected"]
        by_arm[(r["arm"], "got")] += r["got"]
        if r["got"]:
            by_arm[(r["arm"], "definitions")] += 1

    for arm in ("sp", "mp"):
        print(f"{arm}: {by_arm[(arm,'got')]} of {by_arm[(arm,'expected')]} records over "
              f"{by_arm[(arm,'definitions')]} definitions")
    print(f"\nshortfalls ({len(shortfalls)}):")
    for r in shortfalls:
        print(f"  {r['id']:14s} [{r['arm']}]  {r['got']:4d} of {r['expected']:4d}  {r['note']}")

    by_source = collections.Counter(r["source"] for r in report if r["source"])
    print("\nby source:")
    for s, n in by_source.most_common():
        print(f"  {s:26s} {n} definition-arms")

    if not args.dry_run:
        write(store, Path(args.out))
        write_report(store, report, Path(args.out))
        print(f"\nwritten to {args.out}")


if __name__ == "__main__":
    main()
