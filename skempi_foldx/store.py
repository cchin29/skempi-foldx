"""Reading and writing a FoldX result store, and folding several into one.

A store is a directory of ``<skempi_id>.json``, one per SKEMPI interface definition. The shipped
one has a single source per record and needs no folding — but the tools here are what built it,
and what a consumer needs for a store of their own.

:func:`consolidate` exists because FoldX results accumulate across campaigns, each writing its own
``results*/`` directory, and a consumer told to read only one of them sees part of what has been
computed. The shortfall can be large: single-point coverage of 28% against results for another 30%
of rows already sitting on disk unread, a gap closed with no new compute at all.

Two properties of that store are deliberate, and both are the kind of defect that produces a
store which loads without complaint and is wrong:

**Real files, not symlinks.** A store that symlinks into its source directories is not
self-contained: pruning a source silently breaks it, which defeats the point of consolidating.

**Value audits, not key-set audits.** Comparing *mutation keys* between directories and reporting
agreement says nothing about the numbers behind those keys, and that blind spot is wide enough to
hide a join bug: coverage stays identical while a dozen or so mutations carry another mutation's
energies. :func:`audit` compares values.
"""

from __future__ import annotations

import math
import re
import warnings

import json
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .exclusions import code_of, filter_complexes
from .terms import TERMS

#: Per-mutation records live under this key for single-point results and under ``variants``
#: for multi-point ones, because a PDB may appear in both and must not collide.
SINGLE_KEY = "muts"
MULTI_KEY = "variants"


#: ``meta`` fields whose value is a filesystem path on the machine that ran the campaign.
PATH_FIELDS = ("repair_seeded_from",)

#: Directory names after which the *next* component is a username.
_USER_PARENTS = {"Users", "home"}

#: Directory names that are themselves a home, with no username after them.
_BARE_HOMES = {"root"}


def shorten_path(value: str) -> str:
    """A path with the machine-specific prefix removed, keeping enough to explain a record.

    A campaign records where it seeded its repaired structure from, as an absolute path on the
    machine that ran it. That travels into ``meta``, and ``meta`` travels into shipped package
    data — where it names somebody's home directory, permanently.

    Keeps the last three components: campaign, complex, file. A short path can put a username in
    that window -- a home directory two levels above the file -- so a leading home marker and the
    name after it are dropped first. Truncating without that check leaves the username at the
    front of the result, where no absolute-path grep will catch it, because the leading slash is
    gone.
    """
    raw = str(value).replace("\\", "/")
    parts = [p for p in raw.split("/") if p not in ("", ".", "..")]
    # Scan the whole path, not just its front. A home directory can sit at any depth --
    # `/data/home/<user>/...`, `//fileserver/home/<user>/...` -- and anchoring at position 0 leaves
    # the username at the head of the result, where an absolute-path grep cannot see it because
    # the leading slash is gone. The last match wins, so a nested home is cut at the deepest one.
    cut, names = 0, set()
    for i, part in enumerate(parts):
        if part in _USER_PARENTS and i + 1 < len(parts):
            cut = i + 2                      # the marker and the username after it
            names.add(parts[i + 1])
        elif part in _BARE_HOMES:
            cut = i + 1                      # `/root` is the home; nothing follows to drop
    # A username often reappears deeper as an ordinary directory -- a scratch tree named after its
    # owner. Nothing marks those, but the name is known by now, so drop every later repeat of it.
    parts = [p for p in parts[cut:] if p not in names]
    # Never fall back to the raw basename: where stripping consumed everything, the basename is
    # the username itself.
    return "/".join(parts[-3:])


def sanitise_meta(meta: dict) -> dict:
    """``meta`` with every :data:`PATH_FIELDS` entry shortened by :func:`shorten_path`."""
    out = dict(meta)
    for field in PATH_FIELDS:
        value = out.get(field)
        if isinstance(value, str) and value:
            out[field] = shorten_path(value)
    return out


def source_label(path: Path) -> str:
    """A campaign name that survives two campaigns having the same directory name.

    The canonical layout puts the campaign under its tree -- ``foldx_s1102/results`` and
    ``foldx_skempi_full/results`` are different campaigns with the same basename. Keying on the
    basename alone collides them: contribution counts merge into one row, the "safe to retire"
    list names a directory that could be either, and -- worst -- every record gets a ``_source``
    that no longer identifies where the number came from, which is the one job that field has.
    Qualify with the parent, matching the campaign-qualified form the shipped store uses --
    ``sweep_4x_round_1`` for the sweep, ``<identifier>_<arm>`` for the per-definition runs.
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


def reindex_by_skempi_id(records: Dict[str, dict]) -> Dict[str, dict]:
    """Re-key one complex's records by SKEMPI's own mutation string.

    Every shipped record is already keyed by SKEMPI's author-chain form, so this is an identity
    map on the bundled store. It stays because it is not one on a *foreign* store: a campaign run
    in role mode keys its output ``LB38D`` where SKEMPI says ``LI38D``, and this re-indexes it in
    memory. 0.1.0's own store was such a mixture, 3012 records under one convention and 1226 under
    the other.

    The stored keys are deliberately NOT rewritten on disk. A consumer that joins these energies
    onto its own labels may map role chains onto author chains itself, in which case it needs the
    role-chain key to be present; rewriting the store to one convention drops such a consumer's
    coverage without any error, because the join simply stops matching. Choose the convention at
    read time instead.

    Two records re-keying onto one name is a collision, not a merge: the loser is dropped, chosen
    by dict order. That is likeliest on exactly the mixed-convention store this exists for, so it
    warns rather than losing a record quietly. ``FoldxLookup`` records the same class of clash as
    a conflict; this is the standalone path.
    """
    out: Dict[str, dict] = {}
    collided = []
    for key, rec in records.items():
        name = rec.get("cleaned", key)
        if name in out:
            collided.append(name)
            continue                      # first wins, as consolidate() does
        out[name] = rec
    if collided:
        warnings.warn(
            f"{len(collided)} records re-key onto a name another record already holds, e.g. "
            f"{collided[:3]}. The first is kept and the rest dropped, by dict order. A store "
            f"mixing naming conventions can do this; re-index it per convention instead.",
            RuntimeWarning, stacklevel=2)
    return out


def load_bundled_store(which: str = SINGLE_POINT,
                       key: str = "stored") -> Dict[str, Dict[str, dict]]:
    """The shipped results, loaded by name: ``results_sp`` or ``results_mp``.

    ``key="skempi"`` re-indexes each record set by SKEMPI's mutation string. The outer key is the
    SKEMPI identifier either way, so the lookup is ``store["1ACB_E_I"]["LI38D"]``; only the inner
    key changes. On the shipped store that is an identity map -- see :func:`reindex_by_skempi_id`
    for the foreign store it is not one on.
    """
    if key not in ("stored", "skempi"):
        raise ValueError(f'key must be "stored" or "skempi", not {key!r}')
    store = load_store(bundled_path(which))
    if key == "skempi":
        store = {pdb: reindex_by_skempi_id(recs) for pdb, recs in store.items()}
    return store


def load_complex(path: Path) -> Tuple[Dict[str, dict], dict]:
    """Read one ``<skempi_id>.json``, returning ``(records, meta)`` whichever payload key it uses."""
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
    """Read a results directory into ``{skempi_id: {mutation: record}}``.

    The outer key is each file's stem, so a directory written by this package is keyed by SKEMPI
    identifier and one written before 0.2.0 is keyed by PDB code. Both load; :class:`FoldxLookup`
    resolves either.

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
    records. ``None`` if it is empty or holds neither.

    Raises if a single directory holds both. That is not a hypothetical: the two arms are separate
    BuildModel campaigns whose per-complex files share a filename, so pointing both at one
    directory leaves each complex holding whichever arm ran last.

    Returning on the first non-empty file would report one kind for such a directory, and the
    cross-source check in :func:`consolidate` would then see a single consistent kind and pass --
    the guard failing open on exactly the input it exists to reject. Every file is inspected.
    """
    found: Dict[str, str] = {}
    for path in sorted(Path(directory).glob("*.json")):
        records, kind, _ = load_complex_kind(path)
        if records:
            found.setdefault(kind, path.name)
    if len(found) > 1:
        raise ValueError(
            f"{directory} holds both {SINGLE_KEY!r} and {MULTI_KEY!r} records "
            f"(e.g. {found[SINGLE_KEY]} and {found[MULTI_KEY]}). Single- and multi-point results "
            f"are separate campaigns and share a per-complex filename, so one directory cannot "
            f"hold both without each complex meaning whichever arm ran last. Keep them apart."
        )
    return next(iter(found), None)


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
            # One row per (mutation, additional source holding it), NOT one per mutation: a
            # mutation in three sources yields two rows. Reporting the row count as a mutation
            # count overstated the disagreement by ~1.9x in the shipped store (1970 rows over
            # 1027 mutations), so state both and say which is which.
            distinct = len({(c["pdb"], c["mutation"]) for c in self.conflicts})
            lines += ["", f"{distinct} mutations are held by more than one source with differing "
                          f"values, giving {len(self.conflicts)} source-pair disagreements "
                          f"(one row per additional source holding the mutation, so a mutation "
                          f"in three sources contributes two). FoldX is deterministic -- this "
                          f"means the sources passed different mutation lists to BuildModel. "
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
        keep = set(filter_complexes(list(store), on_note=None, context=name))
        store = {pdb: recs for pdb, recs in store.items() if pdb in keep}
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
    """Largest absolute disagreement across the shared numeric terms.

    A non-finite term returns ``inf`` rather than ``nan``. Every comparison against ``nan`` is
    False, so returning it would make ``delta > conflict_tolerance`` false and report two sources
    that disagree as agreeing -- the audit failing open on exactly the values it exists to catch.
    ``json.loads`` accepts bare ``NaN`` and ``Infinity``, so a foreign results directory can
    supply one.
    """
    deltas = []
    for term in TERMS:
        if term not in a or term not in b:
            continue
        x, y = float(a[term]), float(b[term])
        if not (math.isfinite(x) and math.isfinite(y)):
            return math.inf
        deltas.append(abs(x - y))
    # No shared term is not agreement. Returning 0.0 would report two records with nothing in
    # common as identical, which is the same failure as the nan case above reached by a
    # different route: a comparison that cannot be made must not read as one that succeeded.
    return max(deltas) if deltas else math.inf


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


#: A SKEMPI identifier as this package writes one: a PDB code, optionally followed by the two
#: chain groups. Anchored, so nothing here can contain a separator or a parent reference.
_IDENTIFIER_RE = re.compile(r"[A-Za-z0-9]{4}(_[A-Za-z0-9]+){0,2}\Z")


#: A PDB code as SKEMPI writes one. Separate from the identifier pattern because a code names a
#: *structure* -- a file under ``pdb_dir`` and a repair directory -- while an identifier names a
#: record.
_CODE_RE = re.compile(r"[A-Za-z0-9]{4}\Z")


def check_code(code: str) -> str:
    """``code`` unchanged, or ``ValueError``.

    ``SkempiComplex.pdb`` reaches the filesystem twice -- the input structure under ``pdb_dir``
    and the repair directory searched by :meth:`FoldxConfig.find_repaired` -- and it is parsed
    from the same unconstrained column as the identifier. Validated for the same reason and in
    the same place: where the name becomes a path.
    """
    if not isinstance(code, str) or not _CODE_RE.match(code):
        raise ValueError(f"{code!r} is not a PDB code")
    return code


def check_identifier(ident: str) -> str:
    """``ident`` unchanged, or ``ValueError``.

    Every identifier this package turns into a filename passes through here. A store's keys come
    from column 0 of a caller-supplied ``skempi_v2.csv``, and nothing upstream constrains them, so
    ``../../x`` or an absolute path would otherwise be written wherever it pointed -- validated
    where the name becomes a path rather than where it is parsed, because that is the only place
    every route passes through.
    """
    if not isinstance(ident, str) or not _IDENTIFIER_RE.match(ident):
        raise ValueError(f"{ident!r} is not a SKEMPI identifier")
    return ident


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
    for ident, records in store.items():
        check_identifier(ident)
        # Start from the upstream meta so provenance survives — notably
        # `mutation_list_sha256`, which the repair-count comparison relies on to detect a
        # mutation list drifting between rounds. Regenerating meta from scratch would drop it.
        # Sanitised here rather than only in the build script: this is the public, exported way
        # to materialise a store, so a consumer running it over campaign output that still holds
        # absolute paths would otherwise write them into shipped package data.
        meta = sanitise_meta(dict(metas.get(ident, {})))
        meta.update({
            # Both, because they answer different questions and are equal for most entries. The
            # identifier names the record — a pairing's worth of energies. The code names the
            # structure those energies were computed on, which two identifiers can share.
            "skempi_id": ident,
            # code_of on both sides: a campaign written before the record/structure split put
            # the identifier in `pdb`, and keeping it would leave this field naming a record in
            # the one place whose whole job is to name the structure.
            "pdb": code_of(meta.get("pdb") or ident),
            "n_mutations": len(records),
            "sources": sorted({r.get("_source", "?") for r in records.values()}),
        })
        (destination / f"{ident}.json").write_text(
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

    Returns ``(covered, total, missing)``. Computed against an explicit wanted-set rather than
    inferred from the store, since a store cannot report rows nobody asked it for.

    **This is a plain dict join and resolves no naming.** A row labelled in the role-chain
    convention misses a record stored under its author-chain name, and is reported as uncovered
    rather than as unresolvable — on the bundled store that is the difference between 19.6% and
    100% over the same rows, silently. Use :meth:`FoldxLookup.coverage`, which resolves the
    conventions and warns when a shortfall looks like a naming problem; this function is for a
    caller that has already settled naming and wants the arithmetic.
    """
    wanted = list(wanted)
    missing = [(pdb, mut) for pdb, mut in wanted if mut not in store.get(pdb, {})]
    return len(wanted) - len(missing), len(wanted), missing
