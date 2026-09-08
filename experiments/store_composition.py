#!/usr/bin/env python3
"""Composition of the shipped store, one tile per interface definition.

``docs/STORE.md`` says how *complete* the store is -- 4340 of 4343 single-point entries, 1767 of
1850 multi-point -- and completeness is the wrong shape to read those numbers in. The store is
extremely uneven: a handful of protease and antibody definitions hold a third of the records while
a third of definitions hold exactly one. Anyone splitting these records, weighting them, or
averaging a metric over them is making a decision about that unevenness whether or not they know
it, and a table of totals does not show it.

The figure is a treemap: one tile per interface definition, area proportional to how many records
that definition holds.

Two properties make this the right unit rather than the PDB code:

**A definition is what FoldX was actually asked.** ``docs/DETERMINISM.md`` establishes that a
record is a property of ``(repaired structure, mutation list)``. The list belongs to the
definition, so the definition -- not the code -- is the thing whose size predicts a record's
behaviour.

**Three codes carry two definitions each.** ``2C5D``, ``3SE3`` and ``3SE4`` are paired two ways by
SKEMPI. Keyed by code they collapse into one tile, or draw as two tiles with the same label and no
way to tell which is which; keyed by definition they are ``3SE4_B_A`` and ``3SE4_B_C``, and the
picture says which pairing holds what.

Given ``--old``, tiles are coloured by the fraction of their records that moved between that
release and this one, which is the visual form of the finding in ``docs/VERSIONS.md``: whether a
record changed is a property of its definition rather than of the mutation. Without ``--old`` they
are coloured by size class.

    python experiments/store_composition.py --out scratch/composition --plots --html
    python experiments/store_composition.py --old ../skempi-foldx-v0.1.0 \\
        --out scratch/composition --plots --html --csv

``--old`` is the root of a checkout of the earlier release; only its ``skempi_foldx/data`` tree is
read. It needs nothing else -- no SKEMPI table, no structures.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from compare_versions import ARMS, compare, load_store  # noqa: E402

# Size classes, and the colour each gets when no older release is supplied. The thresholds are the
# ones the distribution actually has joints at, not round numbers: >=100 isolates the six
# definitions that dominate, and 1 isolates the third of the store that cannot support any
# within-definition statistic at all.
CLASSES = (
    (100, "holds >=100 records", "#0e7c66"),
    (10, "holds 10-99", "#2a9d7f"),
    (2, "holds 2-9", "#6b7b8c"),
    (1, "holds 1", "#3d4652"),
)


def size_class(n):
    for threshold, label, colour in CLASSES:
        if n >= threshold:
            return label, colour
    return CLASSES[-1][1], CLASSES[-1][2]


# -- layout --------------------------------------------------------------------------------------
def squarify(sizes, x, y, dx, dy):
    """Squarified treemap layout: ``[(x, y, dx, dy)]``, one per size, in the order given.

    Bruls, Huizing & van Wijk (2000). Sizes must be sorted descending and already scaled so that
    they sum to ``dx * dy``. Written out rather than taken from the ``squarify`` package because
    this file's only hard dependency should be the standard library -- matplotlib and plotly are
    already optional here, and a layout algorithm is not worth a third.

    Iterative rather than recursive: the store has 323 single-point definitions, and one frame per
    row placed is a needless way to approach the recursion limit.
    """
    sizes = [float(s) for s in sizes]
    rects = []

    def row(batch, x, y, dx, dy):
        """Place a batch along the shorter side, so tiles stay near-square."""
        covered = sum(batch)
        out = []
        if dx >= dy:
            width = covered / dy if dy else 0.0
            for s in batch:
                out.append((x, y, width, s / width if width else 0.0))
                y += s / width if width else 0.0
        else:
            height = covered / dx if dx else 0.0
            for s in batch:
                out.append((x, y, s / height if height else 0.0, height))
                x += s / height if height else 0.0
        return out

    def worst(batch, x, y, dx, dy):
        placed = row(batch, x, y, dx, dy)
        ratios = [max(w / h, h / w) for _, _, w, h in placed if w > 0 and h > 0]
        return max(ratios) if ratios else float("inf")

    while sizes:
        if dx <= 0 or dy <= 0:
            rects.extend((x, y, 0.0, 0.0) for _ in sizes)
            break
        # Grow the current row while adding one more tile keeps it squarer, then commit it.
        i = 1
        while i < len(sizes) and worst(sizes[:i], x, y, dx, dy) >= worst(sizes[:i + 1], x, y, dx, dy):
            i += 1
        batch, sizes = sizes[:i], sizes[i:]
        rects.extend(row(batch, x, y, dx, dy))
        covered = sum(batch)
        if dx >= dy:
            width = covered / dy if dy else 0.0
            x, dx = x + width, dx - width
        else:
            height = covered / dx if dx else 0.0
            y, dy = y + height, dy - height
    return rects


# -- measurement ---------------------------------------------------------------------------------
def composition(sizes, changed_by_definition=None):
    """``{arm: [tile]}``, tiles sorted descending by record count.

    ``changed_by_definition`` maps ``(stem, arm)`` to ``(changed, paired)`` counts; where it is
    absent or a definition is unpaired, ``fraction`` is ``None``, which ``tile_colour`` renders as
    its own colour rather than as the pale end of the ramp.
    """
    per_arm = {}
    for sub, _, label in ARMS:
        tiles = []
        for (stem, arm), n in sizes.items():
            if arm != sub:
                continue
            changed = paired = None
            fraction = None
            if changed_by_definition is not None:
                changed, paired = changed_by_definition.get((stem, arm), (0, 0))
                fraction = changed / paired if paired else None
            tiles.append({
                "identifier": stem,
                "code": stem.split("_")[0],
                "arm": sub,
                "arm_label": label,
                "records": n,
                "changed": changed,
                "paired": paired,
                "fraction": fraction,
            })
        per_arm[sub] = sorted(tiles, key=lambda t: (-t["records"], t["identifier"]))
    per_arm["combined"] = sorted(
        [t for sub, _, _ in ARMS for t in per_arm[sub]],
        key=lambda t: (-t["records"], t["identifier"]),
    )
    return per_arm


def changed_by_definition(paired):
    """``{(stem, arm): (changed, paired)}`` from ``compare_versions.compare`` output."""
    out = {}
    for p in paired:
        key = (p["identifier"], p["arm"])
        changed, total = out.get(key, (0, 0))
        out[key] = (changed + (1 if p["changed"] else 0), total + 1)
    return out


PANELS = [(sub, label) for sub, _, label in ARMS] + [("combined", "single- and multi-point")]


def concentration(tiles):
    """The numbers the figure is meant to make unmissable, so they can also be quoted."""
    counts = [t["records"] for t in tiles]
    total = sum(counts)
    top10 = sum(counts[:10])
    return {
        "definitions": len(tiles),
        "records": total,
        "top10_records": top10,
        "top10_share": top10 / total if total else 0.0,
        "giant": sum(1 for n in counts if n >= 100),
        "large": sum(1 for n in counts if 10 <= n < 100),
        "medium": sum(1 for n in counts if 2 <= n < 10),
        "singleton": sum(1 for n in counts if n == 1),
        "singleton_records": sum(n for n in counts if n == 1),
        "median": sorted(counts)[len(counts) // 2] if counts else 0,
    }


# -- report --------------------------------------------------------------------------------------
def write_report(per_arm, out: Path, with_change: bool):
    lines = []
    w = lines.append
    w("Composition of the shipped store, by interface definition.")
    w("")
    w("Each definition is one interface pairing of one PDB code, which is the unit FoldX was")
    w("asked about and therefore the unit a record's behaviour is a property of.")
    w("")

    w("Concentration")
    w("-" * 76)
    w(f"{'':<24}{'defs':>7}{'records':>9}{'top 10':>9}{'share':>8}{'median':>8}")
    for key, label in PANELS:
        c = concentration(per_arm[key])
        w(f"{label:<24}{c['definitions']:>7}{c['records']:>9}{c['top10_records']:>9}"
          f"{c['top10_share'] * 100:>7.0f}%{c['median']:>8}")
    w("")
    w("The share column is what the treemap shows: the ten largest definitions hold that much of")
    w("the arm. A metric averaged over definitions and a metric averaged over records are")
    w("therefore answering different questions on this store, and the gap is not small.")
    w("")

    w("Size classes")
    w("-" * 76)
    w(f"{'':<24}{'>=100':>8}{'10-99':>8}{'2-9':>8}{'1':>8}{'recs in singletons':>21}")
    for key, label in PANELS:
        c = concentration(per_arm[key])
        w(f"{label:<24}{c['giant']:>8}{c['large']:>8}{c['medium']:>8}{c['singleton']:>8}"
          f"{c['singleton_records']:>21}")
    w("")
    w("A definition holding one record supports no within-definition statistic at all -- no")
    w("correlation, no rank, no spread. Under the per-structure metrics common in this")
    w("literature those definitions are dropped, and they are a third of the store by count.")
    w("")

    w("The largest definitions")
    w("-" * 76)
    for key, label in PANELS[:2]:
        w(f"{label}:")
        for t in per_arm[key][:10]:
            share = t["records"] / concentration(per_arm[key])["records"] * 100
            tail = ""
            if with_change and t["paired"]:
                tail = f"   {t['changed']}/{t['paired']} changed"
            w(f"    {t['identifier']:<16}{t['records']:>6} records  {share:>5.1f}%{tail}")
        w("")

    doubled = {}
    for key, _ in PANELS[:2]:
        for t in per_arm[key]:
            doubled.setdefault(t["code"], set()).add(t["identifier"])
    twins = {c: sorted(v) for c, v in doubled.items() if len(v) > 1}
    w("Codes carrying more than one definition")
    w("-" * 76)
    if twins:
        for code, ids in sorted(twins.items()):
            w(f"    {code}: {', '.join(ids)}")
        w("")
        w("Keyed by PDB code these would draw as one tile, or as several tiles sharing a label.")
        w("They are the records the 0.2.0 interface fix is about; see docs/VERSIONS.md.")
    else:
        w("    none")
    w("")

    if with_change:
        w("What moved between releases")
        w("-" * 76)
        w(f"{'':<24}{'defs paired':>13}{'wholly still':>14}{'partly moved':>14}{'wholly moved':>14}")
        for key, label in PANELS:
            rows = [t for t in per_arm[key] if t["paired"]]
            still = sum(1 for t in rows if t["fraction"] == 0)
            whole = sum(1 for t in rows if t["fraction"] == 1)
            part = len(rows) - still - whole
            w(f"{label:<24}{len(rows):>13}{still:>14}{part:>14}{whole:>14}")
        w("")
        w("Definitions land at the ends rather than in the middle: a definition's records mostly")
        w("moved together or not at all, which is what it means to say that changing is a")
        w("property of the definition. docs/VERSIONS.md measures this.")
        w("")

    path = out / "STORE_COMPOSITION.txt"
    path.write_text("\n".join(lines) + "\n")
    return path


def write_csv(per_arm, out: Path):
    """One row per definition, so a reader can ask which of theirs are in the tail."""
    path = out / "store_composition.csv"
    with open(path, "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["identifier", "code", "arm", "records", "size_class",
                     "records_paired", "records_changed", "fraction_changed"])
        for sub, _, _ in ARMS:
            for t in per_arm[sub]:
                wr.writerow([
                    t["identifier"], t["code"], t["arm"], t["records"],
                    size_class(t["records"])[0],
                    "" if t["paired"] is None else t["paired"],
                    "" if t["changed"] is None else t["changed"],
                    "" if t["fraction"] is None else f"{t['fraction']:.4f}",
                ])
    return path


# -- figures -------------------------------------------------------------------------------------
UNPAIRED = "#3f6f8f"  # off the red ramp entirely, so it cannot be read as "barely moved"


def tile_colour(t, with_change):
    """Fraction changed on a pale-to-red ramp, or the size class when there is nothing to compare.

    Unpaired definitions get their own colour under ``--old`` rather than the pale end of the
    ramp: no comparison was possible, which is a different statement from no movement, and the
    two should not be adjacent shades.
    """
    if not with_change:
        return size_class(t["records"])[1]
    if t["fraction"] is None:
        return UNPAIRED
    lo, hi = (0.94, 0.94, 0.90), (0.72, 0.11, 0.16)
    return tuple(lo[i] + (hi[i] - lo[i]) * t["fraction"] for i in range(3))


def label_colour(colour):
    """Black or white, whichever the tile can actually be read against.

    The changed-fraction ramp starts near white, so a fixed white label disappears on exactly the
    tiles that carry most of the store.
    """
    from matplotlib.colors import to_rgb

    r, g, b = to_rgb(colour)
    return "black" if (0.299 * r + 0.587 * g + 0.114 * b) > 0.6 else "white"


_WIDTHS = {}


def label_width(text, fig):
    """Width of ``text`` in points per point of font size, measured by the renderer.

    Two cheaper estimates were tried and both overflow. A per-character average is wrong because
    in DejaVu Sans a digit is 0.65 em while ``M`` is 0.86 and ``W`` is 1.01. ``TextPath`` is wrong
    because it returns the ink bounding box, 3-10% narrower than the advance width the renderer
    lays out with -- and since a width-limited label is sized to exactly fill its tile, a 5%
    underestimate makes every one of them spill. Asking the renderer costs one probe per distinct
    label and is exact.
    """
    if text not in _WIDTHS:
        probe = fig.text(0, 0, text, fontsize=10)
        extent = probe.get_window_extent(fig.canvas.get_renderer())
        probe.remove()
        _WIDTHS[text] = extent.width * 72 / fig.dpi / 10
    return _WIDTHS[text]


def write_plots(per_arm, out_dir: Path, with_change: bool, old_label: str = "the older release"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.patches import Patch, Rectangle
    except ImportError:
        print("[plots] matplotlib not installed; skipping", file=sys.stderr)
        return []

    W = H = 100.0
    fig, axes = plt.subplots(1, 3, figsize=(19, 7))
    for ax, (key, label) in zip(axes, PANELS):
        c = concentration(per_arm[key])
        ax.set_xlim(0, W)
        ax.set_ylim(0, H)
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"{label}\n{c['definitions']} definitions, {c['records']} records "
                     f"- the largest 10 hold {c['top10_share'] * 100:.0f}%", fontsize=10)

    if with_change:
        handles = [Patch(facecolor=tile_colour({"fraction": f, "records": 50}, True),
                         edgecolor="white", label=lab)
                   for f, lab in ((0.0, "none of its records moved"), (0.5, "half moved"),
                                  (1.0, "all moved"))]
        handles.append(Patch(facecolor=UNPAIRED, edgecolor="white",
                             label="not paired across releases"))
        subtitle = ("Tile = one interface definition, area proportional to records it holds; "
                    f"colour = the fraction that moved since {old_label}")
    else:
        handles = [Patch(facecolor=col, edgecolor="white", label=lab)
                   for _, lab, col in CLASSES]
        subtitle = ("Tile = one interface definition, area proportional to the number of records "
                    "it holds")
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=9)
    fig.suptitle(f"skempi-foldx store composition\n{subtitle}", fontsize=12)

    # Settle the layout before drawing a single tile, so that what fits inside a tile is measured
    # off the final axes rather than estimated from the figure size. An estimate overflows the
    # narrow tiles, and an identifier spilling into its neighbour is worse than no identifier.
    fig.tight_layout(rect=(0, 0.045, 1, 1))
    fig.canvas.draw()

    for ax, (key, _) in zip(axes, PANELS):
        box = ax.get_window_extent()
        px = box.width * 72 / fig.dpi / W   # points per layout unit, horizontally
        py = box.height * 72 / fig.dpi / H  # and vertically -- the panels are not square on paper
        tiles = per_arm[key]
        total = sum(t["records"] for t in tiles)
        scaled = [t["records"] * W * H / total for t in tiles]
        for t, (x, y, dx, dy) in zip(tiles, squarify(scaled, 0.0, 0.0, W, H)):
            colour = tile_colour(t, with_change)
            ax.add_patch(Rectangle((x, y), dx, dy, facecolor=colour,
                                   edgecolor="white", linewidth=0.5))
            # Fit the label to the tile rather than guessing a threshold: an identifier is 8-12
            # characters, so a fixed size either overflows the narrow tiles or wastes the wide
            # ones. Below 4.2pt it stops being legible, and the tail is texture anyway -- the
            # interactive page is where those get read.
            # 0.92 of the tile, not all of it: the probe measures at 10pt and these labels render
            # at 4-7pt, where hinting rounds advance widths non-linearly, so a label sized to fill
            # its tile exactly still spills by a couple of percent.
            size = min(6.6, 0.92 * dx * px / label_width(t["identifier"], fig))
            if size >= 4.2 and dy * py >= 2.4 * size:
                ax.text(x + dx / 2, y + dy / 2, f"{t['identifier']}\n{t['records']}",
                        ha="center", va="center", fontsize=size, color=label_colour(colour),
                        linespacing=1.15)

    path = out_dir / "store_composition.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return [path]


def write_interactive(per_arm, out_dir: Path, with_change: bool, fragment: bool = False,
                      old_label: str = "the older release"):
    """The same three panels, with every tile hoverable.

    The static figure can only label the tiles big enough to hold text, which is roughly the top
    fifth. The tail is where a reader's own complex usually is, so it needs to be readable by
    pointing at it.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("[html] plotly not installed (pip install plotly); skipping", file=sys.stderr)
        return []

    def rgb(t):
        c = tile_colour(t, with_change)
        if isinstance(c, str):
            return c
        return "rgb({},{},{})".format(*[int(round(255 * v)) for v in c])

    buttons, traces = [], []
    for i, (key, label) in enumerate(PANELS):
        tiles = per_arm[key]
        total = sum(t["records"] for t in tiles)
        custom = [[t["identifier"], t["records"], t["records"] / total * 100, t["arm_label"],
                   size_class(t["records"])[0],
                   "n/a" if t["paired"] is None else f"{t['changed']} of {t['paired']}",
                   "n/a" if t["fraction"] is None else f"{t['fraction'] * 100:.0f}%"]
                  for t in tiles]
        traces.append(go.Treemap(
            # Explicit ids, because a definition present in both arms carries the same identifier
            # twice in the combined panel. Without ids plotly keys its data-join on the label and
            # silently drops the duplicates -- 130 of 477 tiles, a fifth of the record mass, would
            # render as holes in a panel whose title still counted them.
            ids=[f"{t['arm']}:{t['identifier']}" for t in tiles],
            labels=[t["identifier"] for t in tiles],
            parents=[""] * len(tiles),
            values=[t["records"] for t in tiles],
            marker=dict(colors=[rgb(t) for t in tiles], line=dict(color="white", width=1)),
            customdata=custom,
            texttemplate="%{customdata[0]}<br>%{customdata[1]}",
            hovertemplate=(
                "<b>%{customdata[0]}</b>   %{customdata[3]}<br>"
                "%{customdata[1]} records   ·   %{customdata[2]:.1f}% of this set<br>"
                "%{customdata[4]}<br>"
                + (f"moved since {old_label}: %{{customdata[5]}}   (%{{customdata[6]}})"
                   if with_change else "")
                + "<extra></extra>"
            ),
            visible=(i == 0),
            tiling=dict(packing="squarify"),
        ))
        c = concentration(tiles)
        buttons.append(dict(
            label=label, method="update",
            args=[{"visible": [j == i for j in range(len(PANELS))]},
                  {"title": {"text": title_for(label, c, with_change, old_label)}}],
        ))

    c0 = concentration(per_arm[PANELS[0][0]])
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=title_for(PANELS[0][1], c0, with_change, old_label), x=0.5, xanchor="center"),
        updatemenus=[dict(type="buttons", direction="right", buttons=buttons,
                          x=0.5, xanchor="center", y=1.10, yanchor="bottom", showactive=True)],
        margin=dict(t=140, l=20, r=20, b=20), height=760,
    )
    path = out_dir / "store_composition.html"
    fig.write_html(path, include_plotlyjs="cdn", full_html=not fragment)
    return [path]


def title_for(label, c, with_change, old_label="the older release"):
    colour = (f"colour = fraction of its records that moved since {old_label}"
              if with_change else "colour = size class")
    return (f"<b>skempi-foldx store composition — {label}</b><br>"
            f"<span style='font-size:12px'>{c['definitions']} definitions, {c['records']} "
            f"records; the largest 10 hold {c['top10_share'] * 100:.0f}%, and "
            f"{c['singleton']} hold a single record<br>"
            f"tile = one interface definition, area ∝ records; {colour}</span>")


# -- entry point ---------------------------------------------------------------------------------
def check(root: Path, flag: str):
    """``load_store`` a repository root, refusing a path that holds no store.

    ``load_store`` globs, so a wrong path returns an empty dict rather than raising, and every
    figure downstream renders a confident picture of nothing.
    """
    if not (root / "skempi_foldx" / "data").is_dir():
        sys.exit(f"{flag} {root} has no skempi_foldx/data tree")
    records, sizes = load_store(root)
    if not records:
        sys.exit(f"{flag} {root} holds no records")
    return records, sizes


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", type=Path, default=ROOT,
                    help="repository root whose skempi_foldx/data is read (default: this one)")
    ap.add_argument("--old", type=Path,
                    help="checkout of an earlier release; colours tiles by what moved")
    ap.add_argument("--old-label", help="how to name that release in captions "
                                        "(default: the --old directory's trailing version)")
    ap.add_argument("--out", type=Path, required=True, help="output directory")
    ap.add_argument("--csv", action="store_true", help="write one row per definition")
    ap.add_argument("--plots", action="store_true", help="write the static PNG (matplotlib)")
    ap.add_argument("--html", action="store_true", help="write the hoverable page (plotly)")
    ap.add_argument("--html-fragment", action="store_true",
                    help="emit the page as an embeddable div rather than a whole document")
    args = ap.parse_args()

    new, sizes = check(args.store, "--store")
    changes, old_label = None, "the older release"
    if args.old:
        # A mistyped path globs to nothing rather than failing, and the figure that comes out of
        # an empty comparison is not blank -- it is every tile drawn as unpaired, which reads as
        # a finding. Fail instead.
        old, _ = check(args.old, "--old")
        paired, _, _ = compare(old, new, sizes)
        if not paired:
            sys.exit(f"--old {args.old} pairs against none of {args.store}: not the same store")
        changes = changed_by_definition(paired)
        old_label = args.old_label or args.old.resolve().name.split("-")[-1]
        print(f"[compare] {len(paired)} records paired against {args.old}", file=sys.stderr)

    per_arm = composition(sizes, changes)
    args.out.mkdir(parents=True, exist_ok=True)
    written = [write_report(per_arm, args.out, changes is not None)]
    if args.csv:
        written.append(write_csv(per_arm, args.out))
    if args.plots:
        written += write_plots(per_arm, args.out, changes is not None, old_label)
    if args.html or args.html_fragment:
        written += write_interactive(per_arm, args.out, changes is not None,
                                    args.html_fragment, old_label)
    for path in written:
        print(f"wrote {path}")


if __name__ == "__main__":
    main()
