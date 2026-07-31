"""The FoldX result store: many campaign directories -> one canonical set.

FoldX results accumulate across separate campaigns — the S1102 set, three benchmark sets
(S1131 / S2003 / S4169) and the full-SKEMPI remainder — each writing its own ``results*/``
directory. A consumer told to read only one of them sees only part of what has been computed, and
the shortfall can be large: single-point coverage of 28% against results for another 30% of rows
already sitting on disk unread, a gap closed with no new compute at all. This module folds the
directories into one store.

Two properties of that store are deliberate, because both are the kind of thing that looks fine
in review:

**Real files, not symlinks.** A store that symlinks into its source directories is not
self-contained: pruning a source silently breaks it, which defeats the point of consolidating.

**Value audits, not key-set audits.** Comparing *mutation keys* between directories and reporting
agreement says nothing about the numbers behind those keys, and that blind spot is wide enough to
hide a join bug: coverage identical before and after, while ~13 mutations had been handed another
mutation's energies. :func:`audit` compares values.
"""

from __future__ import annotations

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .exclusions import filter_complexes
from .terms import TERMS

#: Per-mutation records live under this key for single-point results and under ``variants``
#: for multi-point ones, because a PDB may appear in both and must not collide.
SINGLE_KEY = "muts"
MULTI_KEY = "variants"


def source_label(path: Path) -> str:
    """A campaign name that survives two campaigns having the same directory name.

    The canonical layout puts the campaign under its tree -- ``foldx_s1102/results`` and
    ``foldx_skempi_full/results`` are different campaigns with the same basename. Keying on the
    basename alone collides them: contribution counts merge into one row, the "safe to retire"
    list names a directory that could be either, and -- worst -- every record gets a ``_source``
    that no longer identifies where the number came from, which is the one job that field has.
    Qualify with the parent, matching the ``foldx_s1102_results_S4169`` form the shipped store
    already uses.
    """
    path = Path(path)
    parent = path.parent.name
    return f"{parent}_{path.name}" if parent else path.name


#: Where the shipped results live inside the installed package.
DATA_DIR = Path(__file__).resolve().parent / "data"

SINGLE_POINT = "results_sp"
MULTI_POINT = "results_mp"


def bundled_path(which: str = SINGLE_POINT) -> Path:
    """Absolute path to a shipped store, wherever the package is installed.

    The results are the point of this package, so they travel inside it. Use this rather than a
    relative ``"data/results_sp"``: that only resolves when the working directory happens to be a
    git clone, so it breaks for anyone reading the results from an installed wheel.
    """
    if which not in (SINGLE_POINT, MULTI_POINT):
        raise ValueError(f"which must be {SINGLE_POINT!r} or {MULTI_POINT!r}, not {which!r}")
    path = DATA_DIR / which
    if not path.is_dir():
        raise FileNotFoundError(
            f"{path} is missing. The shipped results should be installed with the package; "
            f"when running from a source tree, check that skempi_foldx/data/ exists."
        )
    return path


def load_bundled_store(which: str = SINGLE_POINT) -> Dict[str, Dict[str, dict]]:
    """The shipped results, loaded by name: ``results_sp`` or ``results_mp``."""
    return load_store(bundled_path(which))


def load_complex(path: Path) -> Tuple[Dict[str, dict], dict]:
    """Read one ``<pdb>.json``, returning ``(records, meta)`` whichever key it uses."""
    records, _, meta = load_complex_kind(path)
    return records, meta


def load_complex_kind(path: Path) -> Tuple[Dict[str, dict], str, dict]:
    """As :func:`load_complex`, but also reports which payload key held the records.

    The key is the only thing distinguishing a single-point store from a multi-point one, so it
    has to survive a consolidate/write round trip. A writer that always emits ``muts`` silently
    relabels multi-point variants as single-point, leaving the SP/MP guard unable to fire on the
    very stores this module produces -- which is why :func:`write_store` takes the key it was
    read with.
    """
    payload = json.loads(Path(path).read_text())
    if payload.get(SINGLE_KEY):
        return payload[SINGLE_KEY], SINGLE_KEY, payload.get("meta", {})
    if payload.get(MULTI_KEY):
        return payload[MULTI_KEY], MULTI_KEY, payload.get("meta", {})
    # An empty payload is a recorded failure; fall back to whichever key is present.
    kind = MULTI_KEY if MULTI_KEY in payload else SINGLE_KEY
    return {}, kind, payload.get("meta", {})


def load_store(directory: Path) -> Dict[str, Dict[str, dict]]:
    """Read a results directory into ``{pdb: {mutation: record}}``.

    Raises if the directory does not exist, rather than globbing a missing path and returning an
    empty store -- otherwise a wrong path surfaces as a ``KeyError`` on the first lookup, several
    frames from the actual mistake.
    """
    if not Path(directory).is_dir():
        raise FileNotFoundError(f"No results directory at {directory}")
    store: Dict[str, Dict[str, dict]] = {}
    for path in sorted(Path(directory).glob("*.json")):
        records, _ = load_complex(path)
        if records:
            store[path.stem] = records
    return store


def store_kind(directory: Path) -> Optional[str]:
    """Whether a results directory holds single-point (``muts``) or multi-point (``variants``)
    records. ``None`` if it is empty or holds neither."""
    for path in sorted(Path(directory).glob("*.json")):
        records, kind, _ = load_complex_kind(path)
        if records:
            return kind
    return None


@dataclass
class ConsolidationReport:
    """What :func:`consolidate` did, in enough detail to justify the output."""

    complexes: int = 0
    mutations: int = 0
    per_source_unique_complexes: Dict[str, int] = field(default_factory=dict)
    per_source_contributed_mutations: Dict[str, int] = field(default_factory=dict)
    redundant_sources: List[str] = field(default_factory=list)
    conflicts: List[dict] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"{self.complexes} complexes / {self.mutations} mutations",
            "",
            "contribution by source (in precedence order):",
        ]
        for name, count in self.per_source_contributed_mutations.items():
            unique = self.per_source_unique_complexes.get(name, 0)
            lines.append(f"  {name:<28} {count:>6} mutations   {unique:>4} complexes not "
                         f"present in any higher-precedence source")
        if self.redundant_sources:
            lines += ["", "sources contributing nothing unique (safe to retire):"]
            lines += [f"  {name}" for name in self.redundant_sources]
        if self.conflicts:
            lines += ["", f"{len(self.conflicts)} mutations present in more than one source "
                          f"with differing values. FoldX is deterministic -- this means the "
                          f"sources passed different mutation lists to BuildModel. "
                          f"See docs/DETERMINISM.md:"]
            worst = sorted(self.conflicts, key=lambda c: -c["max_abs_delta"])[:10]
            for c in worst:
                lines.append(f"  {c['pdb']:<8} {c['mutation']:<12} "
                             f"max|Δ| {c['max_abs_delta']:.3f}  "
                             f"({c['kept_from']} kept over {c['also_in']})")
        return "\n".join(lines)


def consolidate(
    sources: Sequence[Path],
    destination: Optional[Path] = None,
    conflict_tolerance: float = 1e-6,
) -> Tuple[Dict[str, Dict[str, dict]], ConsolidationReport]:
    """Union every source into one store, earlier sources winning.

    The union is **per mutation**, not per file. Picking, for each complex, the single source
    file with the most mutations is lossless only when that file is a superset of the others —
    a property nothing enforces and a key-set comparison cannot verify.

    ``sources`` is in precedence order: when the same ``(pdb, mutation)`` appears more than
    once, the earliest source's values are kept.

    Two sources can genuinely hold different numbers for the same mutation — but **not because
    FoldX is random**. It is deterministic. The difference arises when the campaigns passed
    different mutation lists to BuildModel, which shifts each entry's position in a sequentially
    evaluated list and changes its result. Prefer the source whose mutation list was most
    complete: it is the one closest to the campaign-invariant answer.
    """
    # Single- and multi-point results come from SEPARATE BuildModel invocations with separate
    # individual_list.txt files, and a mutation's energy depends on every entry preceding it in
    # its own list. Unioning the two into one store would imply they are entries of one list,
    # which they are not -- and would write them back under a single key, erasing which arm each
    # came from. The canonical store keeps them apart (results_sp vs results_mp).
    kinds = {}
    for source in sources:
        kind = store_kind(Path(source))
        if kind:
            kinds.setdefault(kind, []).append(source_label(source))
    if len(kinds) > 1:
        raise ValueError(
            "Refusing to consolidate single-point and multi-point stores together: "
            f"{kinds[SINGLE_KEY]} hold {SINGLE_KEY!r}, {kinds[MULTI_KEY]} hold {MULTI_KEY!r}. "
            "They are separate BuildModel campaigns -- a mutation's energy depends on the whole "
            "preceding mutation list, so these are not two halves of one computation. "
            "Consolidate each arm separately, as the canonical store does."
        )

    merged: Dict[str, Dict[str, dict]] = defaultdict(dict)
    kind = next(iter(kinds), SINGLE_KEY)
    metas: Dict[str, dict] = {}
    report = ConsolidationReport()
    seen_complexes: set = set()

    for source in sources:
        name = source_label(source)
        store = load_store(Path(source))
        for path in sorted(Path(source).glob("*.json")):
            _, _, m = load_complex_kind(path)
            if m:
                metas.setdefault(path.stem, {}).update(
                    {k: v for k, v in m.items() if k not in metas.get(path.stem, {})})
        # A stale JSON for an intractable complex, carried into the canonical store, is picked
        # back up by any scope enumerated from results rather than from the curated table. Drop
        # it here so the store cannot re-arm the compute loop. Any other step that materialises a
        # unified results directory needs the same guard.
        store = {pdb: recs for pdb, recs in store.items()
                 if pdb in filter_complexes(list(store), on_note=None, context=name)}
        contributed = 0
        unique_complexes = 0
        for pdb, records in store.items():
            if pdb not in seen_complexes:
                unique_complexes += 1
            for mutation, record in records.items():
                existing = merged[pdb].get(mutation)
                if existing is None:
                    merged[pdb][mutation] = dict(record, _source=name)
                    contributed += 1
                    continue
                delta = _max_abs_delta(existing, record)
                if delta > conflict_tolerance:
                    report.conflicts.append({
                        "pdb": pdb, "mutation": mutation,
                        "kept_from": existing.get("_source", "?"), "also_in": name,
                        "max_abs_delta": delta,
                    })
        seen_complexes |= set(store)
        report.per_source_contributed_mutations[name] = contributed
        report.per_source_unique_complexes[name] = unique_complexes
        if contributed == 0:
            report.redundant_sources.append(name)

    report.complexes = len(merged)
    report.mutations = sum(len(v) for v in merged.values())

    if destination is not None:
        write_store(merged, Path(destination), kind=kind, metas=metas)
    return dict(merged), report


def _max_abs_delta(a: dict, b: dict) -> float:
    """Largest absolute disagreement across the shared numeric terms."""
    deltas = [abs(float(a[t]) - float(b[t])) for t in TERMS if t in a and t in b]
    return max(deltas) if deltas else 0.0


def infer_kind(store: Dict[str, Dict[str, dict]]) -> str:
    """``variants`` if any record name holds more than one substitution, else ``muts``.

    A multi-point record is keyed by a comma-joined variant (``"VD36I,ED5Q"``); a single-point
    one never contains a comma. Kind is a property of the whole store, not of one complex, so a
    single multi-point record anywhere makes it a multi-point store.

    :func:`load_complex_kind` is the more authoritative answer — it reads the key the file
    actually used rather than inferring from the shape of a name — but it only works when the
    store came off disk. This one works for a store assembled in memory, which is what
    :func:`write_store` needs when its caller does not say.
    """
    for records in store.values():
        if any("," in name for name in records):
            return MULTI_KEY
    return SINGLE_KEY


def write_store(store: Dict[str, Dict[str, dict]], destination: Path,
                kind: Optional[str] = None, metas: Optional[Dict[str, dict]] = None) -> None:
    """Materialize a store as one real ``<pdb>.json`` per complex.

    Real files, not symlinks: the store must stay valid if a source directory is pruned, which
    is the whole point of consolidating. Each record keeps a ``_source`` field so provenance
    survives — since a mutation's ΔΔG depends on the whole ``individual_list.txt`` it was
    computed in, knowing which campaign produced a number is part of being able to explain it.

    ``kind`` defaults to :func:`infer_kind` rather than to ``muts``. This function is public, and
    a default of ``muts`` would mislabel every store written by a caller that omits the argument:
    multi-point variants under the single-point key, which leaves the SP/MP guard unable to fire
    on the stores this module produces. :func:`consolidate` always passes an explicit ``kind``
    read off disk, so the authoritative answer still wins wherever there is one.
    """
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    kind = kind or infer_kind(store)
    metas = metas or {}
    for pdb, records in store.items():
        # Start from the upstream meta so provenance survives — notably
        # `mutation_list_sha256`, which the repair-count comparison relies on to detect a
        # mutation list drifting between rounds. Regenerating meta from scratch would drop it.
        meta = dict(metas.get(pdb, {}))
        meta.update({
            "pdb": pdb,
            "n_mutations": len(records),
            "sources": sorted({r.get("_source", "?") for r in records.values()}),
        })
        (destination / f"{pdb}.json").write_text(
            json.dumps({kind: records, "meta": meta}, indent=1))


def audit(sources: Sequence[Path], conflict_tolerance: float = 1e-6) -> str:
    """Report overlap and **value** agreement across sources, without writing anything."""
    _, report = consolidate(sources, destination=None,
                            conflict_tolerance=conflict_tolerance)
    header = "FoldX result-store audit\n" + "=" * 40
    return f"{header}\n{report.summary()}"


def coverage(
    store: Dict[str, Dict[str, dict]], wanted: Iterable[Tuple[str, str]]
) -> Tuple[int, int, List[Tuple[str, str]]]:
    """How many ``(pdb, mutation)`` pairs the store covers.

    Returns ``(covered, total, missing)``. Coverage is the honest headline number for this
    channel, and it is worth computing against an explicit wanted-set rather than inferring it.
    """
    wanted = list(wanted)
    missing = [(pdb, mut) for pdb, mut in wanted if mut not in store.get(pdb, {})]
    return len(wanted) - len(missing), len(wanted), missing
