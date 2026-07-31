#!/usr/bin/env python3
"""Multi-point arm of the repair-count sweep — the half `repair_ablation.py sweep` omits.

Why this file exists
--------------------
`cmd_sweep` builds its worklist with `worklist_single_point`, so it sweeps single-point mutations
only. The sweep's cost model — 345 complexes, 6006 mutations+variants — covers both arms: 345 is
exactly |SP complexes (323) u MP complexes (153)|. This file drives the multi-point half, without
touching the SP driver or the SP results.

Everything expensive is shared: the repair chain is per COMPLEX, not per mutation, so the four
rounds computed for the SP sweep are reused verbatim here. Only BuildModel is paid again.

Why MP cannot go in the SP results directory
--------------------------------------------
`process_complex` caches at `results_dir/<pdb>.json` and keys the payload "muts" for MODE_AUTHOR
but "variants" for MODE_VARIANT. Writing both arms into one directory would make each complex's
JSON mean whichever arm ran last, and the `out.exists()` resume check would then skip the other
arm entirely. The canonical store keeps them apart for the same reason (results_sp vs results_mp),
so this mirrors that: `round_<r>/results_mp/` and `round_<r>/work_mp/`.

The union rule still holds
--------------------------
A variant's ΔΔG depends on every entry before it in individual_list.txt, so `worklist_multi_point`
takes the complex's whole variant list. Do not subset it. Note that SP and MP are separate
BuildModel invocations with separate lists — that is how the canonical store is built too, so
these values are comparable to it, but an SP and an MP energy for one complex are NOT two entries
of a single list.

Usage
-----
    export FOLDX_BIN=/path/to/foldx
    export PYTHONPATH=<repo root>

    # phase 1 - repair chains. Restrict to complexes NOT already in flight elsewhere.
    python3 sweep_mp.py repair --only-missing --jobs 4

    # phase 2 - BuildModel from each saved round. Refuses to start if a chain is incomplete.
    python3 sweep_mp.py build --rounds 1 2 3 4 --jobs 36
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from skempi_foldx import FoldxConfig, load_skempi, run_campaign  # noqa: E402
from skempi_foldx.run import MODE_VARIANT, worklist_multi_point  # noqa: E402
from skempi_foldx.exclusions import filter_worklist  # noqa: E402

ITERATIONS = 4
BASE = ROOT / "scratch" / "foldx_repair_ablation" / f"sweep_{ITERATIONS}x"


def chain_config(args) -> FoldxConfig:
    """Identical to the SP sweep's chain config, so `complex_work_dir` resolves to the same
    per-complex directory and an already-repaired complex is reused rather than redone."""
    return FoldxConfig(
        results_dir=BASE / "chain_results",
        work_dir=BASE / "chain_work",
        pdb_dir=Path(args.pdb_dir),
        skempi_csv=Path(args.skempi_csv),
        repair_iterations=ITERATIONS,
    )


def mp_worklist(args, skempi):
    pdb_dir = Path(args.pdb_dir)
    work = {p: v for p, v in worklist_multi_point(skempi).items()
            if (pdb_dir / f"{p}.pdb").exists()}
    if args.complexes:
        work = {p: v for p, v in work.items() if p in set(args.complexes)}
    # The registry first -- it is the exclusion you cannot forget to type. `--exclude` is
    # additive on top, for a structure that stalls before there is evidence to make it permanent.
    work = filter_worklist(work, context="multi-point worklist")
    drop = set(args.exclude or [])
    hit = sorted(drop.intersection(work))
    if hit:
        print(f"[mp] excluding {len(hit)}: {', '.join(hit)} "
              f"({sum(len(work[p]) for p in hit)} variants)")
    return {p: v for p, v in work.items() if p not in drop}


def chain_complete(chain: FoldxConfig, pdb: str) -> bool:
    """All ITERATIONS rounds present. `<pdb>_Repair.pdb` alone is not proof: it is moved aside
    between rounds and only survives after the last one, but a process killed mid-write can leave
    a partial file behind. The per-round copies are the durable record."""
    d = chain.complex_work_dir(pdb)
    return all((d / f"repair_round_{i}.pdb").exists() for i in range(1, ITERATIONS + 1))


# ------------------------------------------------------------------------------------- repair
def cmd_repair(args):
    chain = chain_config(args)
    chain.ensure_dirs()
    chain.resolve_binary()
    skempi = load_skempi(chain.skempi_csv)
    work = mp_worklist(args, skempi)

    targets = sorted(work)
    if args.only_missing:
        skipped = [p for p in targets if chain_complete(chain, p)]
        targets = [p for p in targets if not chain_complete(chain, p)]
        print(f"[mp-repair] {len(skipped)} complexes already have all {ITERATIONS} rounds; "
              f"{len(targets)} to do")

    # A complex being repaired by the SP sweep right now would be a genuine collision: two
    # processes writing the same work dir. In-flight complexes have a work dir but no round 4.
    inflight = [p for p in targets if chain.complex_work_dir(p).exists()
                and not chain_complete(chain, p)]
    if inflight and not args.force:
        print(f"[mp-repair] ⚠ {len(inflight)} target(s) already have a work dir without a "
              f"complete chain — the SP sweep may be repairing them RIGHT NOW:")
        print(f"[mp-repair]   {', '.join(inflight[:12])}{'...' if len(inflight) > 12 else ''}")
        print("[mp-repair] Refusing to race it. Re-run once the SP sweep's phase 1 is done, "
              "or pass --force if you know it is not.")
        return

    if not targets:
        print("[mp-repair] nothing to do")
        return
    print(f"[mp-repair] repairing {len(targets)} complexes x{ITERATIONS}, jobs={args.jobs}")
    run_campaign({p: [] for p in targets}, skempi, chain, mode=MODE_VARIANT,
                 jobs=args.jobs, repair_only=True)


# -------------------------------------------------------------------------------------- build
def cmd_build(args):
    chain = chain_config(args)
    chain.resolve_binary()
    skempi = load_skempi(chain.skempi_csv)
    work = mp_worklist(args, skempi)
    rounds = sorted(set(args.rounds or list(range(1, ITERATIONS + 1))))

    incomplete = [p for p in sorted(work) if not chain_complete(chain, p)]
    if incomplete:
        sys.exit(f"[mp-build] {len(incomplete)} MP complexes have no complete repair chain: "
                 f"{', '.join(incomplete[:12])}{'...' if len(incomplete) > 12 else ''}\n"
                 f"[mp-build] Run `sweep_mp.py repair --only-missing` first. Building from a "
                 f"partial chain would silently mix repair counts across complexes.")

    n_var = sum(len(v) for v in work.values())
    print(f"[mp-build] {len(work)} complexes, {n_var} variants, rounds={rounds}")

    for r in rounds:
        print(f"\n[mp-build] BuildModel from repair round {r}")
        cfg = FoldxConfig(results_dir=BASE / f"round_{r}" / "results_mp",
                          work_dir=BASE / f"round_{r}" / "work_mp",
                          pdb_dir=chain.pdb_dir, skempi_csv=chain.skempi_csv)
        cfg.ensure_dirs()
        seeded = 0
        for pdb in sorted(work):
            src = chain.complex_work_dir(pdb) / f"repair_round_{r}.pdb"
            dst = cfg.complex_work_dir(pdb)
            dst.mkdir(parents=True, exist_ok=True)
            # Pre-placing the repaired structure makes process_complex skip repair entirely.
            if not (dst / f"{pdb}_Repair.pdb").exists():
                shutil.copy(src, dst / f"{pdb}_Repair.pdb")
            if not (dst / f"{pdb}.pdb").exists():
                shutil.copy(cfg.pdb_dir / f"{pdb}.pdb", dst / f"{pdb}.pdb")
            seeded += 1
        print(f"[mp-build]   seeded {seeded}/{len(work)} repaired structures")
        run_campaign(work, skempi, cfg, mode=MODE_VARIANT, jobs=args.jobs)

    print(f"\n[mp-build] done. MP stores under {BASE}/round_<r>/results_mp")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["repair", "build"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--rounds", nargs="*", type=int)
    ap.add_argument("--complexes", nargs="*")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="PDB ids to drop. Without this, `build` refuses to start (by design) "
                         "when any MP complex has no complete repair chain — so one structure "
                         "FoldX cannot repair blocks the whole MP arm.")
    ap.add_argument("--only-missing", action="store_true",
                    help="repair: skip complexes whose 4-round chain already exists")
    ap.add_argument("--force", action="store_true",
                    help="repair: proceed even if a target looks in-flight elsewhere")
    ap.add_argument("--pdb-dir", default=str(ROOT / "scratch/skempi2/PDBs"))
    ap.add_argument("--skempi-csv", default=str(ROOT / "scratch/skempi_v2.csv"))
    args = ap.parse_args()
    {"repair": cmd_repair, "build": cmd_build}[args.command](args)


if __name__ == "__main__":
    main()
