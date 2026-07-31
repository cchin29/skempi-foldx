#!/usr/bin/env python3
"""Hash the `individual_list.txt` each round actually fed to BuildModel, and verify it never drifts.

Covers all 345 complexes of a sweep, including runs whose result JSONs carry no recorded
`meta["mutation_list_sha256"]` and would otherwise be reported "not verifiable".

Why this can be done after the fact without weakening the evidence
------------------------------------------------------------------
`process_complex` writes `work/<pdb>/individual_list.txt` and leaves it there. So the exact bytes
BuildModel consumed are still on disk, per complex, per round. Hashing them is *stronger* evidence
than re-deriving the accepted list from `chain_work/` — a re-derivation proves that the list which
would be computed today matches; this proves the list FoldX actually used matched. The only
requirement is to run before the `work*/` trees are deleted, which is why this is wired into the
chain script ahead of any cleanup.

What "drift" would mean
-----------------------
A complex whose list hash differs between rounds was built from a different mutation set in each,
so its per-round ΔΔG difference mixes the repair-count effect with a list-composition effect and
the complex must be dropped from the comparison. This is the failure that leaves no trace in the
energies themselves.

Hash convention
---------------
  file_sha256  SHA-256 of the raw `individual_list.txt` bytes: `"<author>;\\n"` per accepted entry,
               concatenated, trailing newline included. **This is byte-identical to the
               `meta["mutation_list_sha256"]` written by `skempi_foldx.run`** — checked against
               the worked example `'LI38S;\\nGI32Y;\\n'` -> `ba5f85e1...07952`. It is emitted under
               that key name, so these stores verify natively rather than through an alias.

  list_sha256  Entries joined by "\\n", `;` stripped. **DIAGNOSTIC ONLY — never emitted as a peer
               of the canonical key.** The comparison tooling deliberately does not accept this
               normalisation: treating it as equivalent would manufacture mismatches
               indistinguishable from real list drift.

Why hashing the FILE is immune to the author-chain caveat
---------------------------------------------------------
A hash taken over role-chain keys differs from one taken over the file bytes on every SP complex
of a MODE_ROLE campaign, despite identical content. That cannot arise here: this hashes the bytes
`process_complex` wrote, and it writes `f"{author};\\n" for _, author in ordered` — always the
author-chain (SKEMPI cleaned) string. Doubly so for this sweep, which runs MODE_AUTHOR and
MODE_VARIANT, where role == author. The role-chain hazard is real only for a hash rebuilt from
`mapping` keys under MODE_ROLE.

Usage
-----
    python3 list_hashes.py verify                    # report only
    python3 list_hashes.py verify --write-manifest   # + sweep_4x/list_manifest.json (ships w/ results)
    python3 list_hashes.py verify --inject           # + add meta.* keys into each result JSON
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "scratch" / "foldx_repair_ablation" / "sweep_4x"

# arm -> (work dir name, results dir name)
ARMS = {"sp": ("work", "results"), "mp": ("work_mp", "results_mp")}


def hash_list(path: Path):
    raw = path.read_bytes()
    entries = [ln.strip().rstrip(";") for ln in raw.decode().splitlines() if ln.strip()]
    return {
        "file_sha256": hashlib.sha256(raw).hexdigest(),
        "list_sha256": hashlib.sha256("\n".join(entries).encode()).hexdigest(),
        "n_entries": len(entries),
    }


def collect(rounds):
    """{arm: {round: {pdb: hashes}}} for every individual_list.txt on disk."""
    out = {a: defaultdict(dict) for a in ARMS}
    for arm, (workname, _) in ARMS.items():
        for r in rounds:
            wd = BASE / f"round_{r}" / workname
            if not wd.is_dir():
                continue
            for d in sorted(wd.iterdir()):
                f = d / "individual_list.txt"
                if f.is_file():
                    out[arm][r][d.name] = hash_list(f)
    return out


def verify(data, rounds, key="list_sha256"):
    """Report per-arm drift. Returns the number of drifting complexes."""
    total_drift = 0
    for arm in ARMS:
        per_round = data[arm]
        present = [r for r in rounds if per_round.get(r)]
        if not present:
            print(f"[{arm}] no individual_list.txt on disk yet — BuildModel has not run")
            continue

        allc = sorted({p for r in present for p in per_round[r]})
        complete = [p for p in allc if all(p in per_round[r] for r in present)]
        partial = [p for p in allc if p not in complete]

        drift = []
        for p in complete:
            hashes = {per_round[r][p][key] for r in present}
            if len(hashes) > 1:
                drift.append(p)

        print(f"[{arm}] rounds on disk {present} · {len(allc)} complexes "
              f"({len(complete)} in every round, {len(partial)} partial)")
        if drift:
            print(f"[{arm}] ⚠ LIST DRIFT across rounds in {len(drift)} complexes — these must be")
            print(f"[{arm}]   dropped from the repair-count comparison: {', '.join(drift[:20])}")
        elif complete:
            n = sum(per_round[present[0]][p]["n_entries"] for p in complete)
            print(f"[{arm}] ✅ identical across all rounds on disk: {len(complete)}/{len(complete)} "
                  f"complexes, {n} entries. Repair count is the only variable.")
        if partial:
            print(f"[{arm}]   not yet in every round (run still in flight): {len(partial)}")
        total_drift += len(drift)
    return total_drift


def inject(data, rounds):
    """Add the hashes into each result JSON's meta, so a store carries its own provenance."""
    n = 0
    for arm, (_, resname) in ARMS.items():
        for r in rounds:
            rd = BASE / f"round_{r}" / resname
            if not rd.is_dir():
                continue
            for pdb, h in data[arm].get(r, {}).items():
                jf = rd / f"{pdb}.json"
                if not jf.is_file():
                    continue
                payload = json.loads(jf.read_text())
                meta = payload.setdefault("meta", {})
                # The canonical keys, in the byte-identical convention `skempi_foldx.run` uses.
                # Emitting these means the store verifies under `repair_ablation.py compare`
                # natively instead of relying on its file_sha256 alias.
                meta["mutation_list_sha256"] = h["file_sha256"]
                meta["n_accepted"] = h["n_entries"]
                # Deliberately NOT emitting list_sha256 here: different normalisation, and a peer
                # key would be indistinguishable from real drift. It stays in the manifest only.
                meta["mutation_list_sha256_note"] = (
                    "sha256 of individual_list.txt bytes ('<author>;\\n' per accepted entry). "
                    "Added by list_hashes.py from the file BuildModel actually consumed, for a "
                    "run whose results were written without this key."
                )
                jf.write_text(json.dumps(payload, indent=1))
                n += 1
    print(f"[inject] annotated {n} result JSONs")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["verify"])
    ap.add_argument("--rounds", nargs="*", type=int, default=[1, 2, 3, 4])
    # Default to the canonical convention, so the drift verdict is computed over exactly the bytes
    # `meta["mutation_list_sha256"]` hashes rather than an equivalent-but-different normalisation.
    ap.add_argument("--key", default="file_sha256", choices=["file_sha256", "list_sha256"])
    ap.add_argument("--write-manifest", action="store_true")
    ap.add_argument("--inject", action="store_true",
                    help="only safe once the arm has finished writing its result JSONs")
    args = ap.parse_args()

    data = collect(args.rounds)
    drift = verify(data, args.rounds, key=args.key)

    if args.write_manifest:
        out = BASE / "list_manifest.json"
        out.write_text(json.dumps(
            {"rounds": args.rounds, "key_checked": args.key, "drifting_complexes": drift,
             "arms": {a: {str(r): v for r, v in data[a].items()} for a in ARMS}}, indent=1))
        print(f"[manifest] {out}  ({out.stat().st_size/1024:.0f} KB)")
    if args.inject:
        inject(data, args.rounds)

    sys.exit(1 if drift else 0)


if __name__ == "__main__":
    main()
