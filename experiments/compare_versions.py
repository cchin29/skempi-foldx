#!/usr/bin/env python3
"""Comparison of two released stores, record by record.

A rebuilt store is a different measurement rather than a correction of the earlier one, because a
FoldX energy is a property of ``(repaired structure, mutation list)`` rather than of the mutation
alone -- see ``docs/DETERMINISM.md``. That statement alone does not tell a reader whose results
rest on the older release how much of their input moved, or whether it moved anywhere that
matters. This quantifies both.

Three questions, in the order a downstream consumer asks them:

**How much changed?** Records are matched across releases and compared on all twelve terms.
0.1.0 keys single-point records under a mixture of author- and role-chain names, so the match
resolves 0.2.0's ``role`` field as an alias; the report names every record that stays unpaired.

**How large are the changes?** The distribution of the shift in ``Interaction Energy``, plus which
complexes carry it. A median is not enough here: the distribution is heavy-tailed, and a consumer
needs to know whether their rows are in the tail.

**Does the change matter for a model trained on these values?** The decisive comparison is not
between the two releases but between each release and experiment, on the rows they share. Two
stores can disagree substantially, record by record, and still carry the same signal -- and if
they do, a model trained on the older one is not invalidated by the newer.

    python experiments/compare_versions.py --old <checkout of the older release> \\
        --out scratch/version_compare --csv --plots --html

``--old`` is the root of a checkout of the earlier release (a ``git worktree`` at its tag is the
cheapest way to get one). Only its ``skempi_foldx/data`` tree is read.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import random
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skempi_foldx import TERMS  # noqa: E402

ARMS = (("results_sp", "muts", "single-point"), ("results_mp", "variants", "multi-point"))
GAS_CONSTANT = 0.0019872041  # kcal/(mol*K)
DEFAULT_TEMPERATURE = 298.0


# -- loading -------------------------------------------------------------------------------------
def load_store(root: Path):
    """``{(file stem, mutation, arm): record}``, and ``{(file stem, arm): definition size}``.

    The size is how many mutations the definition's list holds, which is the variable that
    predicts whether a record moved between releases at all.
    """
    out, sizes = {}, {}
    for sub, key, _ in ARMS:
        for path in sorted((root / "skempi_foldx" / "data" / sub).glob("*.json")):
            payload = json.load(open(path))
            field = key if key in payload else next(k for k in payload if k != "meta")
            sizes[(path.stem, sub)] = len(payload[field])
            for mut, rec in payload[field].items():
                out[(path.stem, mut, sub)] = rec
    return out, sizes


def index_by_alias(new):
    """``{(code, name, arm): [(stem, record)]}``, indexing each record under its key and its role.

    The older release keys some records by role-chain name, so a code-and-name match has to accept
    either spelling. Where a code carries two interface definitions the same name can reach two
    distinct records; those are reported rather than resolved, since nothing in the older store
    says which pairing its value belonged to.
    """
    idx = defaultdict(list)
    for (stem, mut, arm), rec in new.items():
        code = stem.split("_")[0]
        idx[(code, mut, arm)].append((stem, rec))
        role = rec.get("role")
        if role and role != mut:
            idx[(code, role, arm)].append((stem, rec))
    return idx


# -- experiment ----------------------------------------------------------------------------------
def experimental_ddg(path: Path):
    """``{(code, cleaned): ddG}`` from SKEMPI's affinities, by ``RT*ln(Kd_mut/Kd_wt)``.

    Follows the recipe ``docs/STORE.md`` states for the shipped agreement figures: qualified
    affinities are dropped rather than coerced, a mutation appearing on several rows is counted
    once keeping the first, and a blank temperature stands in as 298 K.
    """
    out = {}
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.reader(fh, delimiter=";")
        next(reader, None)
        for row in reader:
            if len(row) < 14 or not row[0]:
                continue
            code = row[0].split("_")[0]
            cleaned = row[2]
            if (code, cleaned) in out:
                continue                       # keep the first row for a repeated mutation
            try:
                kd_mut, kd_wt = float(row[6]), float(row[8])
            except ValueError:
                continue                       # 'n.b', 'unf', '>1E-04' -- qualified, not numeric
            if kd_mut <= 0 or kd_wt <= 0:
                continue
            try:
                temp = float((row[13] or "").strip().rstrip("()").split("(")[0]) or DEFAULT_TEMPERATURE
            except ValueError:
                temp = DEFAULT_TEMPERATURE
            out[(code, cleaned)] = GAS_CONSTANT * temp * math.log(kd_mut / kd_wt)
    return out


def correlations(xs, ys):
    """Pearson r, Spearman rho and MAE, without pulling in scipy."""
    n = len(xs)
    if n < 3:
        return float("nan"), float("nan"), float("nan")

    def pearson(a, b):
        ma, mb = statistics.fmean(a), statistics.fmean(b)
        da, db = [x - ma for x in a], [y - mb for y in b]
        num = sum(x * y for x, y in zip(da, db))
        den = math.sqrt(sum(x * x for x in da) * sum(y * y for y in db))
        return num / den if den else float("nan")

    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        ranks = [0.0] * len(v)
        i = 0
        while i < len(order):                  # average ranks over ties, as Spearman requires
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            shared = (i + j) / 2 + 1
            for k in range(i, j + 1):
                ranks[order[k]] = shared
            i = j + 1
        return ranks

    mae = statistics.fmean(abs(x - y) for x, y in zip(xs, ys))
    return pearson(xs, ys), pearson(rank(xs), rank(ys)), mae


def bootstrap_delta_rho(old_vals, new_vals, truth, draws: int, seed: int = 0):
    """95% interval for the change in Spearman between releases, resampling records together.

    The two releases are scored on the *same* records, so the comparison is paired and the
    resample has to keep each record's pair intact -- drawing the two independently would compare
    populations rather than releases and would widen the interval for the wrong reason.

    Reported because a difference in the third decimal invites being read as an improvement. An
    interval spanning zero says the data does not separate the releases, which is a different and
    more honest claim than either "better" or "the same".
    """
    rnd = random.Random(seed)
    n = len(truth)
    deltas = []
    for _ in range(draws):
        idx = [rnd.randrange(n) for _ in range(n)]
        t = [truth[i] for i in idx]
        deltas.append(correlations([new_vals[i] for i in idx], t)[1]
                      - correlations([old_vals[i] for i in idx], t)[1])
    deltas.sort()
    return deltas[int(0.025 * draws)], deltas[int(0.975 * draws)]


# -- comparison ----------------------------------------------------------------------------------
def compare(old, new, sizes=None):
    """Pair records across releases; return the pairs and whatever could not be paired."""
    sizes = sizes or {}
    idx = index_by_alias(new)
    paired, unmatched, ambiguous = [], [], []
    for (stem, mut, arm), orec in old.items():
        code = stem.split("_")[0]
        cands = idx.get((code, mut, arm))
        if not cands:
            unmatched.append((stem, mut, arm))
            continue
        if len({id(rec) for _, rec in cands}) > 1:
            ambiguous.append((stem, mut, arm))
            continue
        nstem, nrec = cands[0]
        if any(orec.get(t) is None or nrec.get(t) is None for t in TERMS):
            unmatched.append((stem, mut, arm))
            continue
        changed = any(orec[t] != nrec[t] for t in TERMS)
        paired.append({"code": code, "identifier": nstem, "mutation": mut, "arm": arm,
                       "old": orec, "new": nrec, "changed": changed,
                       "listlen": sizes.get((nstem, arm)),
                       "nsub": mut.count(",") + 1,
                       "delta": nrec["Interaction Energy"] - orec["Interaction Energy"]})
    return paired, unmatched, ambiguous


def quantile(sorted_values, q):
    if not sorted_values:
        return float("nan")
    return sorted_values[min(len(sorted_values) - 1, int(q * len(sorted_values)))]


def write_report(paired, unmatched, ambiguous, exp, out: Path, draws: int = 0):
    lines = []
    w = lines.append
    w("Comparison of two skempi-foldx stores, record by record.")
    w("")
    arm_label = {sub: label for sub, _, label in ARMS}

    w("== How much changed ==")
    w("")
    w(f"{'arm':<14}{'paired':>8}{'identical':>11}{'changed':>9}{'changed %':>11}")
    totals = Counter()
    for sub, _, label in ARMS:
        rows = [p for p in paired if p["arm"] == sub]
        ch = sum(p["changed"] for p in rows)
        totals["paired"] += len(rows)
        totals["changed"] += ch
        pct = 100 * ch / len(rows) if rows else 0
        w(f"{label:<14}{len(rows):>8}{len(rows) - ch:>11}{ch:>9}{pct:>10.1f}%")
    pct = 100 * totals["changed"] / totals["paired"] if totals["paired"] else 0
    w(f"{'all':<14}{totals['paired']:>8}{totals['paired'] - totals['changed']:>11}"
      f"{totals['changed']:>9}{pct:>10.1f}%")
    w("")
    w("'identical' means all twelve terms equal, not merely Interaction Energy.")
    w(f"Unpaired: {len(ambiguous)} ambiguous (a name reaching two interface definitions), "
      f"{len(unmatched)} unmatched.")
    for stem, mut, arm in ambiguous + unmatched:
        w(f"    {stem:<12}{mut:<28}{arm_label[arm]}")
    w("")

    w("== How large the changes are ==")
    w("")
    w("Shift in Interaction Energy over the records that changed, in kcal/mol.")
    w("")
    w(f"{'arm':<14}{'n':>7}{'median':>9}{'mean':>9}{'p90':>9}{'p99':>9}{'max':>9}"
      f"{'<0.5':>8}{'sign flip':>11}")
    for sub, _, label in ARMS:
        rows = [p for p in paired if p["arm"] == sub and p["changed"]]
        if not rows:
            continue
        mags = sorted(abs(p["delta"]) for p in rows)
        flips = sum(1 for p in rows
                    if p["old"]["Interaction Energy"] * p["new"]["Interaction Energy"] < 0)
        under = 100 * sum(1 for m in mags if m < 0.5) / len(mags)
        w(f"{label:<14}{len(mags):>7}{statistics.median(mags):>9.3f}"
          f"{statistics.fmean(mags):>9.3f}{quantile(mags, 0.90):>9.3f}"
          f"{quantile(mags, 0.99):>9.3f}{mags[-1]:>9.3f}{under:>7.0f}%"
          f"{flips:>11}")
    w("")
    w("A sign flip is a record whose Interaction Energy changes sign -- stabilising to")
    w("destabilising or the reverse -- which is the change most likely to alter a conclusion.")
    w("")

    w("Largest individual shifts:")
    w("")
    for p in sorted((p for p in paired if p["changed"]),
                    key=lambda p: -abs(p["delta"]))[:15]:
        # Name the record the way the shipped store does. `mutation` is the 0.1.0 key, which for
        # the doubled codes is the other pairing's spelling -- printing it beside a 0.2.0
        # identifier produced a label that matches no record in either release.
        w(f"    {p['identifier']:<12}{p['new'].get('cleaned', p['mutation']):<30}"
          f"{p['old']['Interaction Energy']:>10.4f} -> {p['new']['Interaction Energy']:>9.4f}"
          f"   ({p['delta']:+.4f})")
    w("")

    w("Complexes where the most records moved:")
    w("")
    per = defaultdict(list)
    for p in paired:
        per[p["identifier"]].append(p)
    ranked = sorted(per.items(),
                    key=lambda kv: (-sum(x["changed"] for x in kv[1]),
                                    -max(abs(x["delta"]) for x in kv[1])))
    w(f"    {'identifier':<14}{'records':>8}{'changed':>9}{'median |d|':>12}{'max |d|':>10}")
    for ident, rows in ranked[:15]:
        ch = [r for r in rows if r["changed"]]
        if not ch:
            continue
        mags = sorted(abs(r["delta"]) for r in ch)
        w(f"    {ident:<14}{len(rows):>8}{len(ch):>9}"
          f"{statistics.median(mags):>12.3f}{mags[-1]:>10.3f}")
    w("")

    w("Movement by term, over the records that changed:")
    w("")
    w(f"    {'term':<26}{'moved':>8}{'median |d|':>12}{'max |d|':>10}")
    changed_rows = [p for p in paired if p["changed"]]
    for term in TERMS:
        deltas = [abs(p["new"][term] - p["old"][term]) for p in changed_rows]
        moved = [d for d in deltas if d != 0]
        if not moved:
            w(f"    {term:<26}{0:>8}{'-':>12}{'-':>10}")
            continue
        w(f"    {term:<26}{len(moved):>8}{statistics.median(moved):>12.4f}{max(moved):>10.4f}")
    w("")

    w("== What characterises the records that moved ==")
    w("")
    w("Whether a record changed is a property of its interface definition: a record moves exactly")
    w("when the list it was computed in changed composition, and a list belongs to a definition.")
    w("")
    for sub, _, label in ARMS:
        defs = defaultdict(list)
        for p in paired:
            if p["arm"] == sub:
                defs[p["identifier"]].append(p)
        touched = {k: v for k, v in defs.items() if any(x["changed"] for x in v)}
        whole = sum(1 for v in defs.values() if all(x["changed"] for x in v))
        w(f"    {label}: {len(touched)} of {len(defs)} paired definitions hold at least one changed "
          f"record, and between them hold all of them ({whole} changed entirely)")
        untouched = [len(v) for k, v in defs.items() if not any(x["changed"] for x in v)]
        if untouched:
            w(f"        unchanged definitions hold at most {max(untouched)} mutations "
              f"({sum(1 for n in untouched if n == 1)} hold exactly one -- a one-entry list is")
            w(f"        the same list in both releases, by construction)")
    w("")
    w("    list length vs the fraction of a definition that moved:")
    w("")
    w(f"        {'arm':<14}{'length':<10}{'definitions':>12}{'mean fraction changed':>24}")
    for sub, _, label in ARMS:
        defs = defaultdict(list)
        for p in paired:
            if p["arm"] == sub:
                defs[p["identifier"]].append(p)
        fracs = [(len(v), sum(x["changed"] for x in v) / len(v)) for v in defs.values()]
        for lo, hi, name in ((1, 5, "1-5"), (6, 20, "6-20"), (21, 60, "21-60"), (61, 10**9, "61+")):
            chunk = [f for n, f in fracs if lo <= n <= hi]
            if chunk:
                w(f"        {label:<14}{name:<10}{len(chunk):>12}{statistics.fmean(chunk):>24.2f}")
    w("")
    w("How far a record moved tracks side-chain repacking freedom. Quartiles of the clash")
    w("energy the record already carried in the older release:")
    w("")
    w(f"        {'arm':<14}{'quartile':<20}{'n':>7}{'median |shift|':>16}")
    for sub, _, label in ARMS:
        rows = sorted((p for p in paired if p["arm"] == sub and p["changed"]),
                      key=lambda p: abs(p["old"]["Van der Waals clashes"]))
        if len(rows) < 8:
            continue
        q = len(rows) // 4
        for i, name in enumerate(("Q1 (least clash)", "Q2", "Q3", "Q4 (most clash)")):
            chunk = rows[i * q:(i + 1) * q] if i < 3 else rows[3 * q:]
            med = statistics.median(sorted(abs(p["delta"]) for p in chunk))
            w(f"        {label:<14}{name:<20}{len(chunk):>7}{med:>16.3f}")
    w("")
    w("    multi-point: shift by number of substitutions")
    w("")
    by_n = defaultdict(list)
    for p in paired:
        if p["arm"] == "results_mp" and p["changed"]:
            by_n[min(p["mutation"].count(",") + 1, 8)].append(abs(p["delta"]))
    w(f"        {'substitutions':<16}{'n':>7}{'median |shift|':>16}")
    for k in sorted(by_n):
        w(f"        {(str(k) + '+' if k == 8 else str(k)):<16}{len(by_n[k]):>7}"
          f"{statistics.median(sorted(by_n[k])):>16.3f}")
    w("")
    w("    single-point: shift by the rotamer freedom of the residue introduced")
    w("")
    chi = {"A": 0, "G": 0, "C": 1, "S": 1, "T": 1, "V": 1, "P": 1, "I": 2, "L": 2, "D": 2,
           "N": 2, "F": 2, "Y": 2, "H": 2, "W": 2, "M": 3, "E": 3, "Q": 3, "K": 4, "R": 4}
    by_chi = defaultdict(list)
    for p in paired:
        if p["arm"] != "results_sp" or not p["changed"]:
            continue
        subs = [s for s in p["mutation"].split(",") if len(s) >= 4]
        if len(subs) == 1 and subs[0][-1] in chi:
            by_chi[chi[subs[0][-1]]].append(abs(p["delta"]))
    w(f"        {'chi angles':<16}{'n':>7}{'median |shift|':>16}")
    for k in sorted(by_chi):
        w(f"        {k:<16}{len(by_chi[k]):>7}{statistics.median(sorted(by_chi[k])):>16.3f}")
    w("")

    w("== Whether the change matters ==")
    w("")
    if not exp:
        w("The comparison against experiment needs --skempi-csv, which was not given. SKEMPI's")
        w("table is the one input this repository does not ship, not being redistributable here;")
        w("every other section above is derived from the two stores alone. The cross-release")
        w("figures below stand without it.")
    else:
        w("Each release against SKEMPI's measured ddG, on the records the two share. A model")
        w("trained on the older values is affected by the difference only to the extent that")
        w("these disagree.")
        w("")
        w("Split three ways, because the records identical in both releases score identically by")
        w("construction and dilute the comparison: the middle block is where any real difference")
        w("has to show up, and the last is the control that must show none.")
        w("")
        w(f"{'arm':<14}{'subset':<22}{'release':<9}{'n':>6}{'Pearson r':>11}"
          f"{'Spearman rho':>14}{'MAE':>8}")
        for sub, _, label in ARMS:
            shared = [p for p in paired
                      if p["arm"] == sub and (p["code"], p["new"].get("cleaned")) in exp]
            if not shared:
                continue
            for name, rows in (("all shared records", shared),
                               ("changed in 0.2.0", [p for p in shared if p["changed"]]),
                               ("identical in both", [p for p in shared if not p["changed"]])):
                if len(rows) < 3:
                    continue
                truth = [exp[(p["code"], p["new"]["cleaned"])] for p in rows]
                for release, key in (("0.1.0", "old"), ("0.2.0", "new")):
                    vals = [p[key]["Interaction Energy"] for p in rows]
                    r, rho, mae = correlations(vals, truth)
                    w(f"{label:<14}{name:<22}{release:<9}{len(rows):>6}{r:>11.3f}"
                      f"{rho:>14.3f}{mae:>8.3f}")
        if draws:
            w("")
            w(f"Paired bootstrap on the change in Spearman, {draws} resamples, records drawn")
            w("together so each keeps its pair:")
            w("")
            w(f"{'arm':<14}{'subset':<22}{'delta rho':>11}{'95% interval':>22}")
            for sub, _, label in ARMS:
                shared = [p for p in paired
                          if p["arm"] == sub and (p["code"], p["new"].get("cleaned")) in exp]
                for name, rows in (("all shared records", shared),
                                   ("changed in 0.2.0", [p for p in shared if p["changed"]])):
                    if len(rows) < 3:
                        continue
                    truth = [exp[(p["code"], p["new"]["cleaned"])] for p in rows]
                    a = [p["old"]["Interaction Energy"] for p in rows]
                    b = [p["new"]["Interaction Energy"] for p in rows]
                    delta = correlations(b, truth)[1] - correlations(a, truth)[1]
                    lo, hi = bootstrap_delta_rho(a, b, truth, draws)
                    w(f"{label:<14}{name:<22}{delta:>+11.3f}"
                      f"{f'[{lo:+.3f}, {hi:+.3f}]':>22}")
            w("")
            w("An interval spanning zero says these records do not separate the two releases.")
    w("")
    w("And the two releases against each other, which is how much a feature column moves:")
    w("")
    w(f"{'arm':<14}{'n':>7}{'Pearson r':>12}{'Spearman rho':>15}{'MAE':>9}")
    for sub, _, label in ARMS:
        rows = [p for p in paired if p["arm"] == sub]
        if not rows:
            continue
        a = [p["old"]["Interaction Energy"] for p in rows]
        b = [p["new"]["Interaction Energy"] for p in rows]
        r, rho, mae = correlations(a, b)
        w(f"{label:<14}{len(rows):>7}{r:>12.3f}{rho:>15.3f}{mae:>9.3f}")
    w("")

    text = "\n".join(lines) + "\n"
    out.write_text(text)
    return text


def write_plots(paired, out_dir: Path):
    """Scatter of the two releases, per arm, plus the distribution of the shift."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[plots] matplotlib not installed; skipping", file=sys.stderr)
        return []

    written = []
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    for col, (sub, _, label) in enumerate(ARMS):
        rows = [p for p in paired if p["arm"] == sub]
        same = [p for p in rows if not p["changed"]]
        diff = [p for p in rows if p["changed"]]

        ax = axes[0][col]
        lo = min(p["old"]["Interaction Energy"] for p in rows)
        hi = max(p["old"]["Interaction Energy"] for p in rows)
        lo = min(lo, min(p["new"]["Interaction Energy"] for p in rows))
        hi = max(hi, max(p["new"]["Interaction Energy"] for p in rows))
        pad = 0.03 * (hi - lo)
        ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], lw=1, color="0.6", zorder=1)
        ax.scatter([p["old"]["Interaction Energy"] for p in same],
                   [p["new"]["Interaction Energy"] for p in same],
                   s=6, alpha=0.35, color="#4c72b0", zorder=2,
                   label=f"unchanged ({len(same)})")
        ax.scatter([p["old"]["Interaction Energy"] for p in diff],
                   [p["new"]["Interaction Energy"] for p in diff],
                   s=8, alpha=0.55, color="#c44e52", zorder=3,
                   label=f"changed ({len(diff)})")
        ax.set_xlabel("0.1.0  Interaction Energy (kcal/mol)")
        ax.set_ylabel("0.2.0  Interaction Energy (kcal/mol)")
        ax.set_title(f"{label}: on the diagonal is unchanged")
        ax.legend(loc="upper left", fontsize=8, frameon=False)

        ax = axes[1][col]
        mags = [abs(p["delta"]) for p in diff]
        if mags:
            ax.hist(mags, bins=60, color="#c44e52", alpha=0.85)
            ax.set_yscale("log")
            ax.axvline(statistics.median(mags), color="0.2", lw=1, ls="--",
                       label=f"median {statistics.median(mags):.3f}")
            ax.legend(fontsize=8, frameon=False)
        ax.set_xlabel("|shift| in Interaction Energy (kcal/mol)")
        ax.set_ylabel("records (log scale)")
        ax.set_title(f"{label}: size of the shift, where it changed")

    fig.suptitle("skempi-foldx 0.1.0 vs 0.2.0, matched records", fontsize=13)
    fig.tight_layout()
    path = out_dir / "version_compare.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    written.append(path)
    return written


def write_csv(paired, unmatched, ambiguous, old, exp, out: Path):
    """One row per 0.1.0 record, with every covariate this analysis derives.

    Written so that a reader can ask the question this page cannot anticipate -- which of *their*
    rows moved, and whether the movement matters for their subset -- without re-deriving the
    pairing, the role aliasing or the experimental join. All twelve terms travel for both releases,
    so the per-term breakdown is reproducible too.

    Every one of 0.1.0's records appears, including the nine that could not be paired: a file that
    silently held only the pairable rows would read as complete while omitting exactly the records
    whose absence is the interesting part. ``status`` says which is which.
    """
    columns = (["status", "identifier", "code", "arm", "mutation", "cleaned", "role",
                "substitutions", "definition_size", "ddg_experimental",
                "ie_0_1_0", "ie_0_2_0", "shift", "changed", "sign_flip",
                "source_0_1_0", "source_0_2_0"]
               + [f"old_{t}" for t in TERMS] + [f"new_{t}" for t in TERMS])

    arm_of = {"results_sp": "single-point", "results_mp": "multi-point"}

    def row_for(p):
        o, n = p["old"], p["new"]
        flip = (o["Interaction Energy"] * n["Interaction Energy"]) < 0
        return {
            "status": "paired", "identifier": p["identifier"], "code": p["code"],
            "arm": arm_of[p["arm"]], "mutation": p["mutation"], "cleaned": n.get("cleaned", ""),
            "role": n.get("role", ""), "substitutions": p["nsub"],
            "definition_size": p["listlen"],
            "ddg_experimental": exp.get((p["code"], n.get("cleaned")), ""),
            "ie_0_1_0": o["Interaction Energy"], "ie_0_2_0": n["Interaction Energy"],
            "shift": round(n["Interaction Energy"] - o["Interaction Energy"], 6),
            "changed": int(p["changed"]), "sign_flip": int(flip),
            "source_0_1_0": o.get("_source", ""), "source_0_2_0": n.get("_source", ""),
            **{f"old_{t}": o[t] for t in TERMS},
            **{f"new_{t}": n[t] for t in TERMS},
        }

    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, restval="")
        writer.writeheader()
        for p in sorted(paired, key=lambda p: (p["arm"], p["identifier"], p["mutation"])):
            writer.writerow(row_for(p))
        for status, keys in (("ambiguous", ambiguous), ("unmatched", unmatched)):
            for stem, mut, arm in sorted(keys):
                rec = old.get((stem, mut, arm), {})
                writer.writerow({
                    "status": status, "identifier": "", "code": stem.split("_")[0],
                    "arm": arm_of[arm], "mutation": mut,
                    "cleaned": rec.get("cleaned", ""),
                    "substitutions": mut.count(",") + 1,
                    "ie_0_1_0": rec.get("Interaction Energy", ""),
                    "source_0_1_0": rec.get("_source", ""),
                    **{f"old_{t}": rec.get(t, "") for t in TERMS},
                })
    return out


def matrix_variables(paired, exp):
    """The columns of the scatterplot matrix: what each record is, and what predicts its shift.

    ``experimental ddG`` is included only when the SKEMPI table was supplied, since it is the one
    quantity that cannot be derived from the two stores.
    """
    cols = [
        ("0.1.0 IE", lambda p: p["old"]["Interaction Energy"]),
        ("0.2.0 IE", lambda p: p["new"]["Interaction Energy"]),
        ("shift", lambda p: p["delta"]),
        ("0.1.0 clash", lambda p: p["old"]["Van der Waals clashes"]),
        ("definition size", lambda p: p["listlen"]),
    ]
    if exp:
        cols.insert(2, ("experimental ddG",
                        lambda p: exp.get((p["code"], p["new"].get("cleaned")))))
    return cols


def write_matrix(paired, exp, out_dir: Path):
    """Scatterplot matrix over the quantities the comparison turns on, split by arm.

    The pairwise panels are the point: agreement with experiment, the two releases against each
    other, and the two variables that predict the shift are all in one field of view, so a
    relationship that holds in one arm and not the other is visible rather than argued.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("[matrix] matplotlib not installed; skipping", file=sys.stderr)
        return []

    cols = matrix_variables(paired, exp)
    arms = (("results_sp", "single-point", "#4c72b0"), ("results_mp", "multi-point", "#c44e52"))
    rows = {sub: [p for p in paired if p["arm"] == sub
                  and all(f(p) is not None for _, f in cols)] for sub, _, _ in arms}

    n = len(cols)
    # Shared per column and per row: without it every panel autoscales alone, and two panels in a
    # column can then show the same variable on different ranges, which is the one thing a matrix
    # is read for. The diagonal's histogram counts go on a twin so they do not fight the row scale.
    fig, axes = plt.subplots(n, n, figsize=(2.35 * n, 2.35 * n), sharex="col", sharey="row")
    for i, (yname, yf) in enumerate(cols):
        for j, (xname, xf) in enumerate(cols):
            ax = axes[i][j]
            if i == j:
                twin = ax.twinx()
                for sub, label, colour in arms:
                    vals = [xf(p) for p in rows[sub]]
                    if vals:
                        twin.hist(vals, bins=40, color=colour, alpha=0.55, label=label)
                twin.set_yticks([])
            else:
                for sub, label, colour in arms:
                    ax.scatter([xf(p) for p in rows[sub]], [yf(p) for p in rows[sub]],
                               s=2.5, alpha=0.30, color=colour, linewidths=0, label=label)
            if i == n - 1:
                ax.set_xlabel(xname, fontsize=9)
            if j == 0:
                ax.set_ylabel(yname, fontsize=9)
            ax.tick_params(labelsize=7)

    handles = [plt.Line2D([], [], marker="o", ls="", color=c, label=f"{lab} ({len(rows[s])})")
               for s, lab, c in arms]
    fig.legend(handles=handles, loc="upper right", frameon=False, fontsize=10)
    fig.suptitle("skempi-foldx 0.1.0 vs 0.2.0 — scatterplot matrix of matched records",
                 fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    path = out_dir / "version_compare_matrix.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    return [path]


def write_interactive_matrix(paired, exp, out_dir: Path):
    """The same matrix with a hover card per record, so an outlying panel names its own points."""
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("[matrix-html] plotly not installed; skipping", file=sys.stderr)
        return []

    cols = matrix_variables(paired, exp)
    fig = go.Figure()
    for sub, label, colour in (("results_sp", "single-point", "#4c72b0"),
                               ("results_mp", "multi-point", "#c44e52")):
        rows = [p for p in paired if p["arm"] == sub
                and all(f(p) is not None for _, f in cols)]
        if not rows:
            continue
        fig.add_trace(go.Splom(
            name=f"{label} ({len(rows)})",
            dimensions=[dict(label=name, values=[f(p) for p in rows]) for name, f in cols],
            marker=dict(size=3, color=colour, opacity=0.4, line_width=0),
            showupperhalf=False, diagonal_visible=False,
            customdata=[[p["identifier"], p["mutation"], p["nsub"], p["listlen"]] for p in rows],
            hovertemplate=("<b>%{customdata[0]}</b>  %{customdata[1]}<br>"
                           "%{xaxis.title.text}: %{x:.4f}<br>%{yaxis.title.text}: %{y:.4f}<br>"
                           "substitutions %{customdata[2]}  ·  definition holds "
                           "%{customdata[3]}<extra></extra>"),
        ))
    size = 250 * len(cols)
    fig.update_layout(
        title=dict(text="skempi-foldx 0.1.0 vs 0.2.0 — scatterplot matrix, hover any point",
                   x=0.5, xanchor="center"),
        width=size, height=size, dragmode="select", hovermode="closest",
        legend=dict(orientation="h", yanchor="bottom", y=-0.08, xanchor="center", x=0.5),
    )
    path = out_dir / "version_compare_matrix.html"
    path.write_text(fig.to_html(full_html=True, include_plotlyjs="inline",
                                config={"displaylogo": False}))
    return [path]


def write_interactive(paired, out_dir: Path, fragment: bool = False):
    """The same four panels as the static figure, with a hover card per record.

    A static scatter answers "how much moved"; it cannot answer "which of *my* rows moved", which
    is the question a reader with their own subset actually has. Each point carries its identifier,
    mutation, both energies, the shift, and the two variables that predict it -- definition size
    and substitution count -- so that question is answered by pointing at the outlier.
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("[html] plotly not installed (pip install plotly); skipping", file=sys.stderr)
        return []

    labels = {sub: label for sub, _, label in ARMS}
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=[f"{labels[s]}: on the diagonal is unchanged" for s, _, _ in ARMS]
                       + [f"{labels[s]}: size of the shift, where it changed" for s, _, _ in ARMS],
        vertical_spacing=0.13, horizontal_spacing=0.09,
    )

    # customdata plus one shared template, rather than a formatted string per point: the strings
    # are what dominate the file size, and 6000 of them make the page several megabytes heavier.
    HOVER = (
        "<b>%{customdata[0]}</b>  %{customdata[1]}<br>"
        "0.1.0 <b>%{x:+.4f}</b>  →  0.2.0 <b>%{y:+.4f}</b> kcal/mol<br>"
        "shift <b>%{customdata[2]:+.4f}</b>%{customdata[6]}<br>"
        "substitutions %{customdata[3]}   ·   definition holds %{customdata[4]} mutations<br>"
        "0.1.0 clash %{customdata[7]:+.4f}   ·   source %{customdata[5]}"
        "<extra></extra>"
    )

    def custom(p):
        flip = (p["old"]["Interaction Energy"] * p["new"]["Interaction Energy"]) < 0
        return [p["identifier"], p["mutation"], round(p["delta"], 4), p["nsub"], p["listlen"],
                p["new"].get("_source", "?"), "  <b>(sign flip)</b>" if flip else "",
                round(p["old"]["Van der Waals clashes"], 4)]

    for col, (sub, _, label) in enumerate(ARMS, start=1):
        rows = [p for p in paired if p["arm"] == sub]
        vals = [p["old"]["Interaction Energy"] for p in rows] + \
               [p["new"]["Interaction Energy"] for p in rows]
        lo, hi = min(vals), max(vals)
        pad = 0.03 * (hi - lo)
        fig.add_trace(go.Scatter(
            x=[lo - pad, hi + pad], y=[lo - pad, hi + pad], mode="lines",
            line=dict(color="rgba(120,120,120,0.8)", width=1),
            hoverinfo="skip", showlegend=False), row=1, col=col)

        for changed, colour, name in ((False, "#4c72b0", "unchanged"), (True, "#c44e52", "changed")):
            pts = [p for p in rows if p["changed"] is changed]
            if not pts:
                continue
            fig.add_trace(go.Scattergl(
                x=[p["old"]["Interaction Energy"] for p in pts],
                y=[p["new"]["Interaction Energy"] for p in pts],
                mode="markers", name=f"{name} ({len(pts)})",
                marker=dict(size=5, color=colour, opacity=0.55),
                customdata=[custom(p) for p in pts], hovertemplate=HOVER,
                legendgroup=f"{label}-{name}",
                showlegend=True, legendgrouptitle_text=label if not changed else None,
            ), row=1, col=col)

        changed_pts = [p for p in rows if p["changed"]]
        if changed_pts:
            mags = [abs(p["delta"]) for p in changed_pts]
            fig.add_trace(go.Histogram(
                x=mags, nbinsx=60, marker_color="#c44e52", opacity=0.85,
                showlegend=False,
                hovertemplate="|shift| %{x} kcal/mol<br><b>%{y}</b> records<extra></extra>",
            ), row=2, col=col)
            fig.add_vline(x=statistics.median(mags), line_dash="dash", line_width=1,
                          line_color="#333", row=2, col=col,
                          annotation_text=f"median {statistics.median(mags):.3f}",
                          annotation_position="top right")

        fig.update_xaxes(title_text="0.1.0  Interaction Energy (kcal/mol)", row=1, col=col)
        fig.update_yaxes(title_text="0.2.0  Interaction Energy (kcal/mol)", row=1, col=col)
        fig.update_xaxes(title_text="|shift| in Interaction Energy (kcal/mol)", row=2, col=col)
        fig.update_yaxes(title_text="records (log scale)", type="log", row=2, col=col)

    fig.update_layout(
        title=dict(text="skempi-foldx 0.1.0 vs 0.2.0, matched records — hover any point",
                   x=0.5, xanchor="center"),
        height=900, width=1250, hovermode="closest",
        hoverlabel=dict(align="left", font_size=12),
        legend=dict(orientation="h", yanchor="bottom", y=-0.13, xanchor="center", x=0.5),
        margin=dict(t=90, b=110),
    )

    path = out_dir / ("version_compare_fragment.html" if fragment else "version_compare.html")
    html = fig.to_html(full_html=not fragment, include_plotlyjs="inline",
                       config={"displaylogo": False,
                               "modeBarButtonsToRemove": ["select2d", "lasso2d"]})
    path.write_text(html)
    return [path]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True, help="checkout root of the earlier release")
    ap.add_argument("--new", default=str(ROOT), help="checkout root of the later release")
    ap.add_argument("--skempi-csv",
                    help="SKEMPI 2.0 table. Optional, and the only input not in the repository: "
                         "it supplies the comparison against experiment and the CSV's "
                         "ddg_experimental column. Everything else derives from the two stores")
    ap.add_argument("--out", required=True, help="directory for the report and plots")
    ap.add_argument("--plots", action="store_true", help="also write the scatter and histogram")
    ap.add_argument("--bootstrap", type=int, default=1000, metavar="N",
                    help="resamples for the paired interval on the change in Spearman "
                         "(default 1000; 0 disables). Needs --skempi-csv to have any effect")
    ap.add_argument("--csv", action="store_true",
                    help="also write one row per 0.1.0 record with every covariate derived here, "
                         "for analysis elsewhere; needs no third-party library")
    ap.add_argument("--html", action="store_true",
                    help="also write an interactive version with a hover card per record "
                         "(needs plotly; the file is self-contained and works offline)")
    ap.add_argument("--html-fragment", action="store_true",
                    help="with --html, emit the body content alone rather than a whole document, "
                         "for embedding in a page that supplies its own <head>")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    old, _ = load_store(Path(args.old))
    new, sizes = load_store(Path(args.new))
    print(f"[load] {len(old)} records in --old, {len(new)} in --new")

    paired, unmatched, ambiguous = compare(old, new, sizes)
    if args.skempi_csv:
        exp = experimental_ddg(Path(args.skempi_csv))
        print(f"[load] {len(exp)} experimental ddG values")
    else:
        exp = {}
        print("[load] no --skempi-csv; skipping the comparison against experiment")

    text = write_report(paired, unmatched, ambiguous, exp, out / "VERSION_COMPARE.txt",
                        draws=args.bootstrap)
    print(text)
    if args.csv:
        path = write_csv(paired, unmatched, ambiguous, old, exp, out / "version_compare.csv")
        print(f"[csv] {path}  ({len(paired) + len(unmatched) + len(ambiguous)} rows)")
    if args.plots:
        for path in write_plots(paired, out) + write_matrix(paired, exp, out):
            print(f"[plot] {path}")
    if args.html:
        for path in (write_interactive(paired, out, fragment=args.html_fragment)
                     + write_interactive_matrix(paired, exp, out)):
            print(f"[html] {path}  ({path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
