#!/usr/bin/env python3
"""WT-residue validation across repair rounds, over all 344 complexes.

The failure this guards against
-------------------------------
`process_complex` derives `repaired_wt_residues` from **that round's** repaired structure and drops
mutations whose WT residue does not match. So the accepted mutation list is a function of the
repaired structure — which is the variable under test. If a repair round altered even one residue
identity, round 1 and round 4 would build from different `individual_list.txt` files, and the
per-round ΔΔG difference would mix the repair-count effect with a list-composition effect.

It leaves no trace in the energies, and it applies to the SP arm as much as the MP one.

Relationship to `list_hashes.py`
--------------------------------
`list_hashes.py` is the stronger check but only covers complexes BuildModel has already reached.
This one needs `chain_work/` alone, so it can run the moment the repair phase finishes — before
any BuildModel — and covers all 344. Run both: this one early and broadly, that one as the record
of what FoldX actually consumed.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from skempi_foldx import load_skempi, worklist_single_point  # noqa: E402
from skempi_foldx.run import (  # noqa: E402
    repaired_wt_residues,
    validate_against_structure,
    worklist_multi_point,
)
from skempi_foldx.skempi import Mutation  # noqa: E402
from skempi_foldx.exclusions import is_excluded  # noqa: E402

ITERATIONS = 4
BASE = ROOT / "scratch" / "foldx_repair_ablation" / f"sweep_{ITERATIONS}x"


def accepted_mp(variants, residues):
    keep = []
    for v in variants:
        subs = [Mutation.parse(s) for s in v.split(",")]
        if not [s for s in subs if residues.get((s.chain, s.position)) != s.wt]:
            keep.append(v)
    return keep


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--skempi-csv", default=str(ROOT / "scratch/skempi_v2.csv"))
    ap.add_argument("--json-out", default=str(BASE / "residue_check.json"))
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="PDB ids not in the sweep. Excluding them here keeps the verdict "
                         "honest: an excluded complex is not a missing chain, and counting it "
                         "as one would report full_coverage=false forever.")
    args = ap.parse_args()

    chain = BASE / "chain_work"
    skempi = load_skempi(Path(args.skempi_csv))
    sp, mp = worklist_single_point(skempi), worklist_multi_point(skempi)
    # Both mechanisms subtract from the universe. If either is missed here, an excluded
    # complex is counted as a MISSING CHAIN and full_coverage is pinned false forever.
    drop = set(args.exclude or []) | {c for c in set(sp) | set(mp) if is_excluded(c)}
    sp = {p: v for p, v in sp.items() if p not in drop}
    mp = {p: v for p, v in mp.items() if p not in drop}
    universe = sorted(set(sp) | set(mp))
    if drop:
        print(f"excluded from the universe: {', '.join(sorted(drop))}")

    done, incomplete = [], []
    for pdb in universe:
        d = chain / pdb
        (done if all((d / f"repair_round_{i}.pdb").exists()
                     for i in range(1, ITERATIONS + 1)) else incomplete).append(pdb)

    print(f"universe (SP u MP): {len(universe)} complexes")
    print(f"  chains complete : {len(done)}")
    print(f"  chains missing  : {len(incomplete)}")

    bad_res, bad_sp, bad_mp = [], [], []
    for pdb in done:
        res = [repaired_wt_residues(chain / pdb / f"repair_round_{i}.pdb")
               for i in range(1, ITERATIONS + 1)]
        if any(r != res[0] for r in res[1:]):
            bad_res.append(pdb)
        if pdb in sp:
            accs = [[a for _, a in validate_against_structure({m: m for m in sp[pdb]}, r)[0]]
                    for r in res]
            if any(a != accs[0] for a in accs[1:]):
                bad_sp.append(pdb)
        if pdb in mp:
            accs = [accepted_mp(mp[pdb], r) for r in res]
            if any(a != accs[0] for a in accs[1:]):
                bad_mp.append(pdb)

    print()
    print(f"residue-set differs across rounds : {len(bad_res)}  {bad_res[:20]}")
    print(f"SP accepted list differs          : {len(bad_sp)}  {bad_sp[:20]}")
    print(f"MP accepted list differs          : {len(bad_mp)}  {bad_mp[:20]}")
    print()
    clean = not (bad_res or bad_sp or bad_mp)
    if clean and not incomplete:
        print(f"✅ CLEAN across all {len(done)}/{len(universe)} complexes — the mutation list is "
              f"identical in every round, so repair count is the only variable.")
    elif clean:
        print(f"✅ clean on {len(done)}/{len(universe)} — but {len(incomplete)} chains are still "
              f"incomplete, so this is NOT yet the full-344 result. Re-run when the sweep ends.")
    else:
        print("⚠ CONFOUNDED for the complexes listed above: their list composition co-varies with "
              "repair round. Drop them from the comparison, or the deltas are not interpretable.")

    Path(args.json_out).write_text(json.dumps(
        {"universe": len(universe), "chains_complete": len(done),
         "chains_missing": len(incomplete), "residue_set_differs": bad_res,
         "sp_list_differs": bad_sp, "mp_list_differs": bad_mp,
         "full_coverage": not incomplete, "clean": clean}, indent=1))
    print(f"\nwrote {args.json_out}")
    sys.exit(0 if clean else 1)


if __name__ == "__main__":
    main()
