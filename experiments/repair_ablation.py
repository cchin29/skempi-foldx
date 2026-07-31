#!/usr/bin/env python3
"""Does iterating `RepairPDB` change the FoldX ΔΔG channel? Measure it.

The pipeline runs `RepairPDB` **once**. That is FoldX's own documented recommendation, and it
matches CATH-ddG — the only comparator paper with a fully specified protocol, and the source of
the FoldX baseline row these results are placed beside. But some groups iterate 5–10 times "to
convergence", and a reviewer may ask why this does not.

The evidence for iterating is thinner than its popularity suggests: the only controlled test
(Usmanova et al., *Bioinformatics* 34(21):3653, 2018) ran ten rounds instead of one and found
stability change and bias unchanged (r = 0.99). That paper is also the origin of the widely-cited
"plateau after 7–10 rounds" line, which is routinely quoted *without* the null result in the
following sentence. Two papers do report gains on binding data, but both confound repair count
with other pipeline changes.

So this is a question worth answering with a measurement rather than a citation — especially since
it is cheap. The per-structure metric only scores complexes with at least ten mutations, and on
the frontier-comparable CATH tier that is **13 complexes**.

What this does
--------------
  probe    RepairPDB N times per complex, no mutations. Reports the heavy-atom RMSD between
           successive rounds, which answers "does it converge, and by which round?" for a few
           minutes of compute.
  run      the full pipeline at --iterations N into a SEPARATE results directory, leaving the
           canonical store untouched.
  sweep    repair N times, then BuildModel from EACH round -- the energy convergence curve.
           The repair chain is paid once and shared across rounds.
  compare  diff two stores per term. It warns when they were built from different mutation
           lists, because that alone changes energies and would confound the measurement.

A clean repair-count measurement needs a CONTROL at --iterations 1 on the same union lists.
Diffing the 5x run against the canonical store is NOT clean: the canonical store was built from
per-campaign mutation subsets, so such a diff mixes repair count with list composition.

Nothing here writes to the canonical store or to any split directory.

Usage
-----
    export FOLDX_BIN=/path/to/foldx
    python experiments/repair_ablation.py probe   --iterations 5 --jobs 8
    python experiments/repair_ablation.py run     --iterations 5 --jobs 8
    python experiments/repair_ablation.py run     --iterations 1 --jobs 8   # control
    python experiments/repair_ablation.py compare \\
        --store    scratch/foldx_repair_ablation/repair_1x/results \\
        --ablation scratch/foldx_repair_ablation/repair_5x/results

`--complexes` defaults to the CATH tier's T>=10 set. Use `--scope t10` for every complex the
per-structure metric can score (98 complexes, ~1 h on 40 cores), or `--scope all` for the whole
single-point store (~2.5 h on 40 cores).

⚠ FoldX is CPU-bound. Do not run this on a machine that is training.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import statistics
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent

# The package is dependency-free and sits at the repo root; no layout disambiguation needed.
sys.path.insert(0, str(ROOT))

from skempi_foldx import (  # noqa: E402
    FoldxConfig,
    MODE_AUTHOR,
    filter_complexes,
    load_skempi,
    load_store,
    process_complex,
    run_campaign,
    worklist_single_point,
)
from skempi_foldx.terms import TERMS  # noqa: E402

#: The CATH tier's T>=10 complexes — the set that actually drives the frontier-comparable
#: per-structure Spearman. Thirteen complexes; about ten minutes on sixteen cores.
CATH_T10 = ["1AK4", "1BRS", "1EMV", "1FFW", "1JTD", "1JTG", "2G2U",
            "2J0T", "2KSO", "2PCC", "2WPT", "3QHY", "3SZK"]


# ---------------------------------------------------------------------------- structure diff
def heavy_atoms(path: Path):
    """``{(chain, resseq, atom): (x, y, z)}`` for non-hydrogen ATOM records."""
    out = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith("ATOM") or len(line) < 54:
            continue
        name = line[12:16].strip()
        if name.startswith("H") or line[76:78].strip() == "H":
            continue
        try:
            out[(line[21], line[22:27].strip(), name)] = (
                float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
    return out


def rmsd(a: Path, b: Path):
    """Heavy-atom RMSD over shared atoms, without superposition.

    No superposition on purpose: `RepairPDB` does not move the backbone frame, so the quantity
    of interest is how far atoms moved in place. Superposing would mask exactly that.
    """
    x, y = heavy_atoms(a), heavy_atoms(b)
    shared = set(x) & set(y)
    if not shared:
        return None, 0
    total = sum((x[k][0] - y[k][0]) ** 2 + (x[k][1] - y[k][1]) ** 2 + (x[k][2] - y[k][2]) ** 2
                for k in shared)
    return math.sqrt(total / len(shared)), len(shared)


# ---------------------------------------------------------------------------------- scoping
def resolve_complexes(args, store_dir: Path, skempi=None, pdb_dir=None):
    """Which complexes to run.

    `all` is derived from SKEMPI intersected with the available PDBs, NOT from a results store:
    reading the store means that on a machine without that particular directory the scope
    resolves to the empty set and the run completes successfully having done nothing. `t10` does
    need a store — "has at least ten scored mutations" is a property of results — so it says so
    rather than silently returning nothing.
    """
    if args.complexes:
        targets = sorted(set(args.complexes))
    elif args.scope == "cath":
        targets = list(CATH_T10)
    elif args.scope == "t10":
        store = load_store(store_dir)
        if not store:
            sys.exit(f"--scope t10 needs a results store to count mutations, but {store_dir} "
                     f"is empty or missing. Pass --store, or use --scope all.")
        targets = sorted(p for p, m in store.items() if len(m) >= 10)
    elif skempi is not None:
        targets = sorted(p for p, e in skempi.items() if e.single or e.multi)
    else:
        store = load_store(store_dir)
        if not store:
            sys.exit(f"--scope all could not enumerate complexes: {store_dir} is empty or "
                     f"missing and no SKEMPI table was supplied.")
        targets = sorted(store)

    # Before the PDB check, so the note reads as an exclusion rather than a missing file.
    targets = filter_complexes(targets, context=f"--scope {args.scope}")

    # --exclude is additive on top of the registry, never a replacement for it. The registry is
    # what cannot be forgotten; the flag is the escape hatch for a complex that is not yet
    # knowledge -- a structure that stalls mid-sweep, before there is enough evidence to justify
    # writing it into INTRACTABLE. Applied after --complexes as well, because an exclusion is a
    # statement about the structure itself, not about how it got selected.
    drop = set(args.exclude or [])
    hit = sorted(drop.intersection(targets))
    if hit:
        print(f"[scope] --exclude dropping {len(hit)}: {', '.join(hit)}")
        print("[scope] an id passed twice belongs in skempi_foldx/exclusions.py "
              "with its evidence -- a flag typed twice is a registry entry waiting to happen.")
        targets = [p for p in targets if p not in drop]

    if pdb_dir is not None:
        missing = [p for p in targets if not (Path(pdb_dir) / f"{p}.pdb").exists()]
        if missing:
            print(f"[scope] {len(missing)} complexes have no PDB, skipped: "
                  f"{', '.join(missing[:6])}{'...' if len(missing) > 6 else ''}")
        targets = [p for p in targets if (Path(pdb_dir) / f"{p}.pdb").exists()]

    if not targets:
        sys.exit("Resolved 0 complexes — refusing to run, since this would exit successfully "
                 "having done nothing. Check --scope / --store / --pdb-dir.")
    return targets



def make_config(args, suffix: str) -> FoldxConfig:
    """A config writing into its own directory tree — never the canonical store."""
    base = ROOT / "scratch" / "foldx_repair_ablation" / suffix
    return FoldxConfig(
        results_dir=base / "results",
        work_dir=base / "work",
        pdb_dir=Path(args.pdb_dir),
        skempi_csv=Path(args.skempi_csv),
        repair_iterations=args.iterations,
    )


# ------------------------------------------------------------------------------------ probe
def cmd_probe(args):
    """Repair each complex N times and report how far it moves between rounds."""
    config = make_config(args, f"probe_{args.iterations}x")
    config.ensure_dirs()
    binary = config.resolve_binary()
    skempi = load_skempi(config.skempi_csv)
    targets = resolve_complexes(args, Path(args.store), skempi, config.pdb_dir)
    print(f"[probe] {len(targets)} complexes x {args.iterations} repair rounds")
    print(f"[probe] FoldX: {binary}")

    run_campaign({p: [] for p in targets if p in skempi}, skempi, config,
                 mode=MODE_AUTHOR, jobs=args.jobs, repair_only=True)

    print(f"\n{'complex':<10}" + "".join(f"  r{i}->r{i+1}" for i in range(1, args.iterations)))
    rows = defaultdict(list)
    for pdb in targets:
        work = config.complex_work_dir(pdb)
        rounds = [work / f"repair_round_{i}.pdb" for i in range(1, args.iterations + 1)]
        if not all(r.exists() for r in rounds):
            print(f"{pdb:<10}  (incomplete)")
            continue
        deltas = []
        for i in range(len(rounds) - 1):
            value, _ = rmsd(rounds[i], rounds[i + 1])
            deltas.append(value)
            rows[i].append(value)
        print(f"{pdb:<10}" + "".join(f"  {d:7.4f}" for d in deltas if d is not None))

    print()
    for i in sorted(rows):
        vals = [v for v in rows[i] if v is not None]
        if vals:
            print(f"round {i+1} -> {i+2}: median RMSD {statistics.median(vals):.4f} A   "
                  f"max {max(vals):.4f} A   (n={len(vals)})")
    print("\nConvergence would show as these shrinking toward zero. If round 1->2 is already "
          "small, a single repair is not leaving anything material on the table.")


# -------------------------------------------------------------------------------------- run
def cmd_run(args):
    config = make_config(args, f"repair_{args.iterations}x")
    config.ensure_dirs()
    config.resolve_binary()
    skempi = load_skempi(config.skempi_csv)
    targets = resolve_complexes(args, Path(args.store), skempi, config.pdb_dir)

    # Union mutation lists, deliberately: a mutation's energy depends on the whole list, so a
    # per-dataset subset would confound the repair-count effect with a list-composition effect.
    worklist = {p: m for p, m in worklist_single_point(skempi, only=targets).items()}
    n_mut = sum(len(v) for v in worklist.values())
    print(f"[run] {len(worklist)} complexes, {n_mut} mutations, "
          f"repair x{args.iterations}, jobs={args.jobs}")
    print(f"[run] writing to {config.results_dir}  (canonical store untouched)")
    run_campaign(worklist, skempi, config, mode=MODE_AUTHOR, jobs=args.jobs)
    print(f"\nNow: python {Path(__file__).name} compare "
          f"--ablation {config.results_dir} --store {args.store}")


def list_hashes(directory):
    """``{pdb: mutation_list_sha256}`` for the complexes that recorded one.

    The accepted mutation list is derived from that round's repaired structure — process_complex
    drops any mutation whose wild-type residue does not match — so two stores can differ in list
    composition without anything obvious going wrong, and a diff across that difference is not
    measuring what it claims to.
    """
    out = {}
    for path in sorted(Path(directory).glob("*.json")):
        meta = json.loads(path.read_text()).get("meta", {})
        # `file_sha256` is an accepted alias: same normalization, raw individual_list.txt bytes.
        # `list_sha256` (';' stripped, '\n'-joined, no trailing newline) is deliberately NOT
        # accepted -- it is a different normalization, so treating it as equivalent would
        # manufacture mismatches that look like real list drift.
        h = meta.get("mutation_list_sha256") or meta.get("file_sha256")
        if h:
            out[path.stem] = h
    return out


# ---------------------------------------------------------------------------------- compare
def cmd_compare(args):
    canonical = load_store(Path(args.store))
    ablation = load_store(Path(args.ablation))
    if not ablation:
        sys.exit(f"No results under {args.ablation} — run the `run` subcommand first.")

    per_term = defaultdict(list)
    pairs = 0
    only_in_ablation = 0
    for pdb, muts in ablation.items():
        for mut, record in muts.items():
            base = canonical.get(pdb, {}).get(mut)
            if not base:
                only_in_ablation += 1
                continue
            pairs += 1
            for term in TERMS:
                if term in base and term in record:
                    per_term[term].append(abs(float(base[term]) - float(record[term])))

    h_ref, h_abl = list_hashes(Path(args.store)), list_hashes(Path(args.ablation))
    shared_h = set(h_ref) & set(h_abl)
    mismatched = sorted(p for p in shared_h if h_ref[p] != h_abl[p])
    unverifiable = len(set(canonical) & set(ablation)) - len(shared_h)

    print(f"  reference : {args.store}")
    print(f"  ablation  : {args.ablation}")
    print(f"compared {pairs} mutations over {len(ablation)} complexes")
    if mismatched:
        print(f"\n⚠ MUTATION LIST DIFFERS on {len(mismatched)} complexes: "
              f"{', '.join(mismatched[:6])}{'...' if len(mismatched) > 6 else ''}")
        print("  A mutation's energy depends on every entry preceding it in the list, so these")
        print("  complexes are NOT a clean comparison — the list changed as well as the")
        print("  variable under test. The list is derived from the repaired structure, so a")
        print("  round altering a residue identity is one way this happens.")
    elif shared_h:
        print(f"  mutation lists verified identical on {len(shared_h)} complexes ✓")
    if unverifiable:
        print(f"  ({unverifiable} complexes predate list-hash recording — not verifiable)")
    print()

    # A mutation-count mismatch means the two stores were built from DIFFERENT mutation lists.
    # Because a mutation's energy depends on every entry preceding it in individual_list.txt,
    # that difference alone changes values -- so a diff across it is NOT a clean measurement of
    # whatever else changed. This is easy to do by accident and impossible to see in the output.
    if only_in_ablation:
        print(f"⚠ CONFOUNDED COMPARISON — {only_in_ablation} mutations exist only in the")
        print("  ablation store, so the two runs used different mutation lists. A mutation's")
        print("  energy depends on the whole preceding list, so these deltas mix the effect")
        print("  under test with a list-composition effect.")
        print("  For a clean repair-count measurement, generate a control with the SAME lists:")
        print("      python experiments/repair_ablation.py run --iterations 1 ...")
        print("      python experiments/repair_ablation.py compare \\")
        print("          --store  scratch/foldx_repair_ablation/repair_1x/results \\")
        print("          --ablation scratch/foldx_repair_ablation/repair_5x/results\n")

    print(f"{'term':<24} {'n diff':>7} {'%':>6} {'median':>9} {'p90':>9} {'max':>9}")
    for term in TERMS:
        vals = per_term.get(term, [])
        diff = [v for v in vals if v > 1e-6]
        if not vals:
            continue
        p90 = sorted(diff)[int(0.9 * len(diff)) - 1] if diff else 0.0
        print(f"{term:<24} {len(diff):>7} {100*len(diff)/len(vals):>5.1f}% "
              f"{statistics.median(diff) if diff else 0:>9.4f} {p90:>9.4f} "
              f"{max(diff) if diff else 0:>9.4f}")

    print("\nWhat to conclude (for a CLEAN comparison — same lists, differing only in repair):")
    print("  - all-zero  -> repair count does not matter here; report it and keep one repair.")
    print("  - small     -> keep one repair, quote these numbers as the sensitivity.")
    print("  - large     -> the choice is load-bearing; rescore the FoldX-alone baseline")
    print("                 against the published FoldX row (0.383 here vs 0.430 published)")
    print("                 before deciding which protocol to publish.")



# ------------------------------------------------------------------------------------ sweep
def cmd_sweep(args):
    """Repair N times, then run BuildModel from EACH round — the energy convergence curve.

    Two phases, because the repair chain is shared. Phase 1 produces `repair_round_<r>.pdb` for
    every round. Phase 2 seeds each round's structure into its own work tree as
    `<pdb>_Repair.pdb`, which makes `process_complex` skip repair and build straight from it.

    So N repair rounds are paid once and reused across every round's BuildModel, rather than
    re-repairing per arm.
    """
    rounds = sorted(set(args.rounds or list(range(1, args.iterations + 1))))
    if max(rounds) > args.iterations:
        sys.exit(f"--rounds asks for round {max(rounds)} but --iterations is {args.iterations}")

    base = ROOT / "scratch" / "foldx_repair_ablation" / f"sweep_{args.iterations}x"
    chain = FoldxConfig(results_dir=base / "chain_results", work_dir=base / "chain_work",
                        pdb_dir=Path(args.pdb_dir), skempi_csv=Path(args.skempi_csv),
                        repair_iterations=args.iterations)
    chain.ensure_dirs()
    chain.resolve_binary()
    skempi = load_skempi(chain.skempi_csv)
    targets = resolve_complexes(args, Path(args.store), skempi, chain.pdb_dir)
    worklist = worklist_single_point(skempi, only=targets)
    n_mut = sum(len(v) for v in worklist.values())

    print(f"[sweep] {len(targets)} complexes, {n_mut} single-point mutations")
    print(f"[sweep] ⚠ SINGLE-POINT ONLY — multi-point variants are not covered by this "
          f"version. The CATH test set is ~39% multi-point, so a CATH-*all* number needs both.")
    print(f"[sweep] phase 1/2: {args.iterations} repair rounds (shared)")
    run_campaign({p: [] for p in targets}, skempi, chain, mode=MODE_AUTHOR,
                 jobs=args.jobs, repair_only=True)

    for r in rounds:
        print(f"\n[sweep] phase 2/2: BuildModel from repair round {r}")
        cfg = FoldxConfig(results_dir=base / f"round_{r}" / "results",
                          work_dir=base / f"round_{r}" / "work",
                          pdb_dir=chain.pdb_dir, skempi_csv=chain.skempi_csv)
        cfg.ensure_dirs()
        seeded = 0
        for pdb in targets:
            src = chain.complex_work_dir(pdb) / f"repair_round_{r}.pdb"
            if not src.exists():
                continue
            work = cfg.complex_work_dir(pdb)
            work.mkdir(parents=True, exist_ok=True)
            # Pre-placing the repaired structure makes process_complex skip repair entirely.
            if not (work / f"{pdb}_Repair.pdb").exists():
                shutil.copy(src, work / f"{pdb}_Repair.pdb")
            if not (work / f"{pdb}.pdb").exists():
                shutil.copy(cfg.pdb_dir / f"{pdb}.pdb", work / f"{pdb}.pdb")
            seeded += 1
        print(f"[sweep]   seeded {seeded}/{len(targets)} repaired structures")
        run_campaign({p: v for p, v in worklist.items()}, skempi, cfg,
                     mode=MODE_AUTHOR, jobs=args.jobs)

    print(f"\n[sweep] done. Stores under {base}/round_<r>/results")
    print(f"[sweep] compare adjacent rounds, e.g.:")
    print(f"    python {Path(__file__).name} compare \\")
    print(f"        --store    {base}/round_1/results \\")
    print(f"        --ablation {base}/round_{max(rounds)}/results")


# ------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["probe", "run", "sweep", "compare"])
    ap.add_argument("--iterations", type=int, default=5,
                    help="RepairPDB rounds (default 5; the literature's claimed plateau is 7-10)")
    ap.add_argument("--jobs", type=int, default=4, help="parallel complexes")
    ap.add_argument("--scope", choices=["cath", "t10", "all"], default="cath",
                    help="cath = the 13 CATH T>=10 complexes (default, ~10 min on 16 cores); "
                         "t10 = all 98 metric-scoring complexes; all = the full store")
    ap.add_argument("--complexes", nargs="*", help="explicit PDB ids, overriding --scope")
    ap.add_argument("--exclude", nargs="*", default=[],
                    help="PDB ids to drop, on top of skempi_foldx.exclusions.INTRACTABLE. For a "
                         "structure that stalls before there is evidence to make it permanent.")
    ap.add_argument("--rounds", nargs="*", type=int,
                    help="sweep only: which repair rounds to BuildModel from "
                         "(default: every round 1..iterations)")
    ap.add_argument("--store", default=str(ROOT / "skempi_foldx/data/results_sp"),
                    help="canonical store to compare against")
    ap.add_argument("--ablation", default=None, help="ablation store (compare only)")
    ap.add_argument("--pdb-dir", default=str(ROOT / "scratch/skempi2/PDBs"))
    ap.add_argument("--skempi-csv", default=str(ROOT / "scratch/skempi_v2.csv"))
    args = ap.parse_args()

    if args.command == "compare" and not args.ablation:
        args.ablation = str(ROOT / "scratch/foldx_repair_ablation" /
                            f"repair_{args.iterations}x" / "results")
    {"probe": cmd_probe, "run": cmd_run, "sweep": cmd_sweep,
     "compare": cmd_compare}[args.command](args)


if __name__ == "__main__":
    main()
