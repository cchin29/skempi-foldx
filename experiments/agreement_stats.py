#!/usr/bin/env python3
"""Recompute the cross-campaign agreement statistics quoted in ``docs/DETERMINISM.md``.

Those statistics were recorded once, by hand, against a population that could not later be
reconciled -- the document carried "1528 shared mutations, 10631 campaign-pair comparisons"
alongside a determinism table over 5873 pairs, and a third denominator of 3096 for sign
stability. Three populations, no reconciliation, and the campaign directories were not to hand,
so the numbers sat marked **unverified**.

This script recomputes every one of them from the campaign directories, so the document can
quote measured values with a stated denominator.

**The campaign directories are not shipped here** -- point ``--campaign-root`` at an extracted
snapshot of the ``scratch/`` directory the campaigns wrote to.

Two findings worth stating up front, both established by running this:

* **10631 is not reproducible from any subset of the campaign directories.** The correct single
  population is 5870 comparisons over 1527 shared mutations (5873 / 1530 if the stale 1KBH JSON
  is left in). 10631 and 5873 were describing the same thing.
* **The four- vs five-campaign question was moot.** ``full_skempi`` shares *zero* mutations with
  the four benchmark campaigns, so both readings give the same population -- which is why the
  document's C(4,2)-vs-C(5,2) arithmetic never resolved anything.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from skempi_foldx.terms import TERMS  # noqa: E402

SCALAR = "Interaction Energy"

CAMPAIGNS = [
    ("foldx_skempi_full/results", "full_skempi"),
    ("foldx_s1102/results_S4169", "S4169"),
    ("foldx_s1102/results_S2003", "S2003"),
    ("foldx_s1102/results_S1131", "S1131"),
    ("foldx_s1102/results", "S1102"),
]


def load(root: Path, drop_intractable: bool):
    """-> {campaign: {(pdb, mutation): record}}"""
    from skempi_foldx.exclusions import INTRACTABLE
    out = {}
    for rel, name in CAMPAIGNS:
        d = root / rel
        if not d.is_dir():
            raise SystemExit(f"missing campaign directory: {d}\n"
                             "Point --campaign-root at the directory the campaigns wrote to.")
        recs = {}
        for f in sorted(d.glob("*.json")):
            pdb = f.stem
            if drop_intractable and pdb in INTRACTABLE:
                continue
            payload = json.loads(f.read_text())
            for mut, rec in (payload.get("muts") or {}).items():
                recs[(pdb, mut)] = rec
        out[name] = recs
    return out


def icc_one_way(groups):
    """ICC(1,1) by one-way random effects, unequal group sizes.

    groups: list of lists of measurements of the same target.
    """
    groups = [g for g in groups if len(g) >= 2]
    if not groups:
        return float("nan")
    k = len(groups)
    n_i = [len(g) for g in groups]
    N = sum(n_i)
    grand = sum(sum(g) for g in groups) / N
    ms_b = sum(len(g) * (statistics.fmean(g) - grand) ** 2 for g in groups) / (k - 1)
    ms_w = sum((x - statistics.fmean(g)) ** 2 for g in groups for x in g) / (N - k)
    # n0: average group size correction for unequal sizes
    n0 = (N - sum(n * n for n in n_i) / N) / (k - 1)
    denom = ms_b + (n0 - 1) * ms_w
    return (ms_b - ms_w) / denom if denom else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign-root", type=Path, required=True)
    ap.add_argument("--keep-intractable", action="store_true",
                    help="include 1KBH's stale JSON (reproduces the 1530/5873 variant)")
    args = ap.parse_args()

    data = load(args.campaign_root, drop_intractable=not args.keep_intractable)

    # Which campaigns hold each mutation.
    holders = defaultdict(list)
    for name, recs in data.items():
        for key in recs:
            holders[key].append(name)
    shared = {k: v for k, v in holders.items() if len(v) >= 2}

    print(f"campaigns loaded      : {', '.join(f'{n}={len(r)}' for n, r in data.items())}")
    print(f"shared mutations      : {len(shared)}")

    # Does full_skempi overlap the benchmark campaigns at all?
    bench = {"S4169", "S2003", "S1131", "S1102"}
    cross = sum(1 for v in shared.values() if "full_skempi" in v and bench & set(v))
    print(f"full_skempi x benchmark overlap: {cross} mutations "
          f"({'four- and five-campaign populations are identical' if cross == 0 else 'they differ'})")

    # Pairwise comparisons on the scalar term.
    deltas, pairs = [], 0
    per_campaign_vals = defaultdict(dict)
    sign_flips, flip_mags, flip_big, flip_huge = 0, [], 0, 0
    by_complex = defaultdict(lambda: [0, 0, 0.0])   # pdb -> [pairs, differing, sum|Δ|]
    for key, names in shared.items():
        vals = {n: data[n][key].get(SCALAR) for n in names}
        vals = {n: v for n, v in vals.items() if isinstance(v, (int, float))}
        for n, v in vals.items():
            per_campaign_vals[n][key] = v
        for a, b in itertools.combinations(sorted(vals), 2):
            d = vals[a] - vals[b]
            deltas.append(d); pairs += 1
            by_complex[key[0]][0] += 1
            by_complex[key[0]][2] += abs(d)
            if abs(d) > 1e-9:
                by_complex[key[0]][1] += 1
            if vals[a] * vals[b] < 0:
                sign_flips += 1
                m = min(abs(vals[a]), abs(vals[b]))
                flip_mags.append(m)
                if m > 0.5:
                    flip_big += 1
                if m > 2.0:
                    flip_huge += 1

    ab = sorted(abs(d) for d in deltas)
    allv = [v for m in per_campaign_vals.values() for v in m.values()]
    rmsd = math.sqrt(sum(d * d for d in deltas) / len(deltas))
    sd_d = statistics.pstdev(deltas)
    differing = [d for d in deltas if abs(d) > 1e-9]

    print(f"\ncampaign-pair comparisons : {pairs}")
    print(f"differing (>1e-9)         : {len(differing)}  ({100*len(differing)/pairs:.1f}%)")
    print(f"pooled RMSD               : {rmsd:.3f}")
    print(f"median |d| all pairs      : {statistics.median(ab):.4f}")
    print(f"median |d| differing only : {statistics.median([abs(d) for d in differing]):.4f}")
    print(f"95% limits of agreement   : +/-{1.96*sd_d:.3f}")
    print(f"implied single-choice SD  : {rmsd/math.sqrt(2):.3f}")
    print(f"population SD of values   : {statistics.pstdev(allv):.3f}")
    print(f"p99 |d|                   : {ab[int(0.99*len(ab))]:.3f}")
    print(f"max |d|                   : {ab[-1]:.3f}")

    groups = [[data[n][k].get(SCALAR) for n in names
               if isinstance(data[n][k].get(SCALAR), (int, float))]
              for k, names in shared.items()]
    print(f"ICC (scalar)              : {icc_one_way(groups):.3f}")

    print("\nper-term ICC:")
    rows = []
    for term in TERMS:
        g = [[data[n][k].get(term) for n in names
              if isinstance(data[n][k].get(term), (int, float))]
             for k, names in shared.items()]
        rows.append((icc_one_way(g), term))
    for v, term in sorted(rows):
        print(f"  {v:.3f}  {term}")

    print(f"\nsign flips                : {sign_flips}/{pairs} = {100*sign_flips/pairs:.1f}%")
    if flip_mags:
        print(f"  median smaller magnitude: {statistics.median(flip_mags):.3f}")
    print(f"  both sides >0.5 kcal/mol: {flip_big}")
    print(f"  both sides >2.0 kcal/mol: {flip_huge}")

    tot = sum(v[2] for v in by_complex.values())
    ranked = sorted(by_complex.items(), key=lambda kv: -kv[1][2])
    agree_all = sum(1 for _, v in by_complex.items() if v[1] == 0)
    disagree_all = sum(1 for _, v in by_complex.items() if v[1] == v[0] and v[0] > 0)
    print(f"\ncomplexes with shared muts: {len(by_complex)}")
    print(f"  agree on every mutation : {agree_all}")
    print(f"  disagree on every one   : {disagree_all}")
    print(f"  top 3 share of |d|      : {100*sum(v[2] for _, v in ranked[:3])/tot:.0f}%"
          f"  ({', '.join(k for k, _ in ranked[:3])})")
    print(f"  top 10 share of |d|     : {100*sum(v[2] for _, v in ranked[:10])/tot:.0f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
