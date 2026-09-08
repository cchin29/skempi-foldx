#!/usr/bin/env python3
"""Energy convergence against structural convergence, across repair rounds.

The repair-count sweep repairs each complex four times and runs ``BuildModel`` from every round,
holding the mutation list byte-constant so repair count is the only variable. It produces the
*energy* curve. The *structure* curve comes from ``repair_ablation.py probe``, which measures
heavy-atom RMSD between consecutive repair rounds.

The question matters because repairing to round 4 costs roughly four times repairing once. If
the energy were stable by round 2 while the structure kept moving, the cheaper protocol would be
defensible on the quantity that is actually consumed.

**The answer depends on which statistic, and the interesting part is why.** Run against the
completed 4-round sweep:

* The **median** does converge, at about the same rate as the structure — single-point medians
  fall 0.045 → 0.020 → 0.007 kcal/mol across 1→2, 2→3, 3→4, a factor of ~2.5 per step. Measured
  against the full 344-complex structure curve rather than the 13-complex pilot, the typical
  single-point mutation settles *faster* than the structure — see ``docs/DETERMINISM.md``,
  *Repair count*, which carries the reconciled figures.
* The **tail does not converge at all.** The 99th percentile of |Δ| goes 2.13 → 1.70 → 1.65
  kcal/mol single-point, and 2.95 → 2.86 → 2.86 multi-point — essentially flat. The maximum is
  not even monotone: single-point 3→4 moves one entry by 11.3 kcal/mol, *more* than 2→3's 8.5.
* The **fraction of entries that move at all never settles**: 92% / 91% / 87% single-point, and
  97% / 96% / 96% multi-point. There is no round after which the store stops changing.
* **Multi-point converges more slowly than single-point** on every statistic, which is expected —
  a multi-point entry repacks two or more positions, so it has more to be sensitive about.

The operational reading: a repair-count choice cannot be justified from the median, because the
median is not where the risk is. Whatever round is chosen, ~90% of entries differ from the
neighbouring round, and a small number differ by more than a kcal/mol. That is an argument for
*fixing and recording* the repair count, which the store does, rather than for tuning it.

**On comparing this with the structure curve.** Both are now measured over the same 344
complexes from the repair chain that produced these energies, so the rates and the endpoints are
comparable. Neither speaks to the 5x endpoint: the sweep ran four rounds.

Usage — the sweep tree is not shipped, so point at an extracted one::

    python experiments/repair_convergence.py --sweep-root <...>/sweep_4x
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

TERM = "Interaction Energy"
ARMS = [("results", "muts", "single-point"), ("results_mp", "variants", "multi-point")]


def load_round(root: Path, rnd: int, arm: str, key: str, term: str):
    out = {}
    d = root / f"round_{rnd}" / arm
    if not d.is_dir():
        raise SystemExit(f"missing {d}")
    for f in sorted(d.glob("*.json")):
        for mut, rec in (json.loads(f.read_text()).get(key) or {}).items():
            v = rec.get(term)
            if isinstance(v, (int, float)):
                out[(f.stem, mut)] = v
    return out


def quantile(sorted_vals, q):
    if not sorted_vals:
        return float("nan")
    return sorted_vals[min(len(sorted_vals) - 1, int(q * len(sorted_vals)))]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sweep-root", type=Path, required=True,
                    help="the sweep_4x directory holding round_1..round_4")
    ap.add_argument("--rounds", type=int, default=4)
    ap.add_argument("--term", default=TERM)
    args = ap.parse_args()

    rounds = list(range(1, args.rounds + 1))
    for arm, key, label in ARMS:
        data = {r: load_round(args.sweep_root, r, arm, key, args.term) for r in rounds}
        n = len(data[rounds[0]])
        print(f"\n=== {label}: {n} entries, {args.term} ===")
        print(f"{'transition':<12}{'n':>7}{'differ':>9}{'':>6}"
              f"{'median':>10}{'p90':>9}{'p99':>9}{'max':>9}")
        steps = [(a, b) for a, b in zip(rounds, rounds[1:])] + [(rounds[0], rounds[-1])]
        for a, b in steps:
            shared = data[a].keys() & data[b].keys()
            deltas = sorted(abs(data[a][k] - data[b][k]) for k in shared)
            moved = [d for d in deltas if d > 1e-9]
            print(f"{f'{a}->{b}':<12}{len(deltas):>7}{len(moved):>9}"
                  f"{f'{100*len(moved)/len(deltas):.0f}%':>6}"
                  f"{statistics.median(deltas):>10.4f}{quantile(deltas,0.90):>9.3f}"
                  f"{quantile(deltas,0.99):>9.3f}{deltas[-1]:>9.2f}")

        # Where the tail lives: a protocol decision hinges on whether the movers are a
        # handful of complexes (bounded risk, nameable) or spread across the store.
        a, b = rounds[-2], rounds[-1]
        per = {}
        for (pdb, mut) in data[a].keys() & data[b].keys():
            per.setdefault(pdb, []).append(abs(data[a][(pdb, mut)] - data[b][(pdb, mut)]))
        worst = sorted(per.items(), key=lambda kv: -max(kv[1]))[:5]
        tot = sum(sum(v) for v in per.values())
        top = sum(sum(v) for _, v in worst)
        print(f"  final step {a}->{b}: top 5 complexes hold "
              f"{100*top/tot:.0f}% of all absolute movement — "
              + ", ".join(f"{p}({max(v):.2f})" for p, v in worst))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
