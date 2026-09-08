"""The FoldX compute engine: structure + mutation list -> per-complex ΔΔG JSON.

Paths travel in a :class:`~skempi_foldx.config.FoldxConfig` passed to each task rather than in
module-level globals. Globals would make the engine's behaviour a function of whichever driver
assigned to them last, and every driver would need a process-pool initializer to re-apply them in
each worker. Passing the config explicitly keeps campaigns independent and needs no initializer.

One complex, once
-----------------
    RepairPDB           slow, and the result is reusable across campaigns
    BuildModel          one mutant + one wild-type structure per mutation
    AnalyseComplex      interaction energies for both, on the SKEMPI chain groups
    ΔΔG_int = IE(mutant) - IE(wild-type), per term

Resumable at complex granularity: a complex whose JSON already exists is skipped.

Determinism and mutation-list dependence
----------------------------------------
FoldX 5.1 BuildModel **is deterministic**. Measured over 1719 pairs of results computed from an
identical repaired structure *and* an identical ``individual_list.txt`` prefix: zero differ, to
the last decimal.

But a mutation's ΔΔG is a function of ``(repaired structure, every entry preceding it in
individual_list.txt)`` — not of the mutation alone. BuildModel walks the list sequentially in one
process, and each entry's side-chain optimisation inherits the state left by the entries before
it. It also re-optimises a fresh wild-type reference per entry, so *both* terms of the
subtraction move.

The practical consequence is severe: computing a **subset** of a complex's mutations gives
different numbers than computing the full set, because the subset shifts every entry's position.
Across five campaigns that each covered a different SKEMPI subset, this produced spreads of
median 0.067 and up to 11.03 kcal/mol on the same mutation of the same structure — with no
randomness involved at all.

**So: always compute the union of a complex's mutations, never a per-dataset subset.**
:func:`worklist_single_point` does this by construction (it takes everything SKEMPI lists for the
complex); :func:`worklist_from_table` does **not**, and is retained only for reproducing campaigns
that were driven from a per-dataset mutation table. Results are cached per complex, so computing
the union once serves every downstream dataset.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .config import FoldxConfig
from .exclusions import ALLOW_ENV, code_of, filter_worklist, is_excluded
from .skempi import (
    Mutation,
    SkempiComplex,
    map_role_to_author,
    parse_interaction_file,
    repaired_wt_residues,
    validate_against_structure,
)
from .store import check_code, shorten_path
from .terms import TERMS

# How the mutation strings handed to a campaign should be interpreted.
MODE_ROLE = "role"
"""Dataset mutations naming chains by role (A/B); need mapping onto author chains."""
MODE_AUTHOR = "author"
"""Cleaned SKEMPI mutations, already on author chains; self-map."""
MODE_VARIANT = "variant"
"""Comma-joined multi-point variants, already on author chains."""


@dataclass
class ComplexResult:
    pdb: str
    mutations: Dict[str, dict]
    meta: dict
    status: str

    def to_json(self) -> dict:
        key = "variants" if self.meta.get("mode") == MODE_VARIANT else "muts"
        return {key: self.mutations, "meta": self.meta}


def format_individual_list(ordered: Sequence[Tuple[str, str]]) -> str:
    """The exact bytes handed to BuildModel: one ``<author>;`` per line.

    ``ordered`` is ``(key, author)`` pairs. The **author** string is written — the SKEMPI
    ``cleaned`` form on the real PDB chain — never the key, which under ``MODE_ROLE`` is a
    role-chain string (``LB38S`` where the author form is ``LI38S``).

    That distinction matters beyond FoldX correctness: ``meta["mutation_list_sha256"]`` is a hash
    of *these bytes*, and it is what detects a mutation list drifting between repair rounds.
    Rebuilding that hash from the role→author ``mapping`` keys instead of from ``ordered`` would
    hash role-chain strings, so every single-point complex in a ``MODE_ROLE`` campaign would
    mismatch a store hashed from the file — a difference indistinguishable from real drift.
    Hash the file, or build from ``ordered``; never from the mapping keys.
    """
    return "".join(f"{author};\n" for _, author in ordered)


def _foldx(argv: Sequence[str], cwd: Path, log: Path) -> None:
    """Run FoldX with the argument vector as given, no shell.

    ``pdb`` and ``entry.groups`` reach these arguments from column 0 of a caller-supplied
    ``skempi_v2.csv``, and a shell would make that column executable: an identifier containing
    ``$(...)`` runs it. Nothing here needs a shell -- there is no pipe, glob or redirection, and
    the log is captured by file handle -- so the vector is passed straight to ``execve``. It also
    removes the quoting the string form needed around ``binary`` for a path with a space in it.
    """
    with open(log, "a") as fh:
        subprocess.run([str(a) for a in argv], shell=False, check=False,
                       cwd=cwd, stdout=fh, stderr=fh)


def _finish(config: FoldxConfig, pdb: str, meta: dict, status: str) -> ComplexResult:
    """Write an empty result so a failure is recorded and not silently retried forever."""
    meta["error"] = status
    result = ComplexResult(pdb, {}, meta, status)
    config.result_path(pdb).write_text(json.dumps(result.to_json(), indent=1))
    return result


def process_complex(
    pdb: str,
    mutations: Sequence[str],
    entry: SkempiComplex,
    config: FoldxConfig,
    mode: str = MODE_ROLE,
    terms: Sequence[str] = TERMS,
    repair_only: bool = False,
) -> ComplexResult:
    """Run the full pipeline for one interface definition and cache the result.

    ``pdb`` names the *record*: it is the results filename and the work directory. ``entry.pdb``
    names the *structure*. They are the same string for a complex SKEMPI defines once, and they
    differ when a caller computes each definition of a multiply-defined code separately, where
    the record wants the full identifier and the structure is still one file on disk.

    ``repair_only`` stops after the RepairPDB rounds — used by the repair-count convergence
    probe, which has no mutations to build."""
    excluded = is_excluded(entry.pdb)
    if excluded is not None:
        # Second guard, for a driver that builds its own targets and never calls run_campaign.
        # Nothing is written: a results JSON here would land in the store, be picked up by the
        # next scope enumeration, and recreate the situation this exists to prevent.
        warnings.warn(f"EXCLUDED: {excluded}. Set {ALLOW_ENV}=1 to run it anyway.",
                      RuntimeWarning, stacklevel=2)
        return ComplexResult(pdb, {}, {"skempi_id": pdb, "pdb": entry.pdb,
                                       "excluded": excluded.reason}, "excluded")

    siblings = tuple(s for s in (config.sibling_records or {}).get(entry.pdb, ()) if s != pdb)

    out = config.result_path(pdb)
    if out.exists():
        payload = json.loads(out.read_text())
        key = "variants" if mode == MODE_VARIANT else "muts"
        other = "muts" if key == "variants" else "variants"
        # Single- and multi-point results share a filename but not a payload key. Reading a
        # single-point file in multi-point mode finds no records, so a bare `.get(key, {})` would
        # report "cached" and skip an entire arm in silence. They must live in separate results
        # dirs; if they do not, say so rather than quietly producing nothing.
        if not payload.get(key) and payload.get(other):
            raise RuntimeError(
                f"{out} holds {other!r} results but {mode!r} mode wants {key!r}. Single- and "
                f"multi-point campaigns must use separate results_dir values -- sharing one "
                f"makes each complex's JSON mean whichever arm ran last, and the resume check "
                f"then skips the other arm entirely."
            )
        return ComplexResult(pdb, payload.get(key, {}), payload.get("meta", {}), "cached")

    binary = config.resolve_binary()
    work = config.complex_work_dir(pdb)
    work.mkdir(parents=True, exist_ok=True)
    log = work / "foldx.log"
    log.write_text("")
    # Both, because they differ whenever a code is computed one definition at a time, and only
    # `pdb` (the structure) can be recovered from neither the other field nor the filename.
    meta = {"skempi_id": pdb, "pdb": entry.pdb, "mode": mode, "n_requested": len(mutations)}

    # The structure is named by the code; everything inside the work directory is named by the
    # record. Copying it in under the record's name is what lets two definitions of one code run
    # side by side without their intermediates -- which FoldX names after the input stem --
    # overwriting each other.
    source = config.pdb_dir / f"{check_code(entry.pdb)}.pdb"
    if not source.exists():
        return _finish(config, pdb, meta, "no_pdb")
    if not (work / f"{pdb}.pdb").exists():
        shutil.copy(source, work / f"{pdb}.pdb")

    # --- RepairPDB (or reuse one from a sibling campaign) ------------------------------------
    repaired = work / f"{pdb}_Repair.pdb"
    # A partially-written repair chain is an inconsistent state, and it must be detected even
    # when `<pdb>_Repair.pdb` happens to exist — otherwise the resume short-circuits and the
    # missing rounds stay missing while the run reports success.
    if config.repair_iterations > 1:
        _rounds = [work / f"repair_round_{i}.pdb"
                   for i in range(1, config.repair_iterations + 1)]
        if any(r.exists() for r in _rounds) and not all(r.exists() for r in _rounds):
            for r in _rounds:
                r.unlink(missing_ok=True)
            repaired.unlink(missing_ok=True)
            meta["restarted_incomplete_repair_chain"] = True
    if not repaired.exists():
        # By structure, then by any sibling record of the same structure. A repair belongs to
        # the code, but a work directory is named after the *record*, so a sibling definition's
        # repair sits under its own identifier. Passing this record's own name instead finds
        # nothing by construction: that path is `repaired`, which the enclosing branch has
        # already established does not exist.
        seeded = (config.find_repaired(entry.pdb, records=siblings)
                  if config.repair_iterations == 1 else None)
        if seeded and seeded != repaired:
            # Written aside and renamed, because a sibling definition of the same code may be
            # searching this directory concurrently. `find_repaired` tests existence only, so a
            # half-copied structure is indistinguishable from a finished one -- and it would go
            # straight into BuildModel. Unreachable before siblings were searched at all.
            staged = repaired.with_suffix(".pdb.partial")
            shutil.copy(seeded, staged)
            staged.replace(repaired)
            # Shortened at the point of writing: this lands in the shipped store, so an
            # absolute path here names the machine that ran the campaign and is permanent once
            # published.
            meta["repair_seeded_from"] = shorten_path(str(seeded))
        else:
            # Iterated repair feeds each output back in as the next input. The input
            # `<pdb>.pdb` is NEVER mutated: the chain runs in `<pdb>_chain.pdb`, and FoldX names
            # its output after the input stem, so each round reads `_chain` and writes
            # `_chain_Repair`. Every round is preserved for convergence analysis.
            #
            # Running the chain over `<pdb>.pdb` itself and restoring the original after the loop
            # would compound on a crash: a crash between rounds leaves `<pdb>.pdb` holding a
            # partially repaired structure, the next run copies that over `_original.pdb` —
            # destroying the record of the true input — and repairs N more times from it. A
            # 5-round resume after a round-3 crash then yields an 8x-repaired structure reported
            # as "repaired x5", with the saved rounds holding 4x-8x.
            chain = work / f"{pdb}_chain.pdb"
            chain_out = work / f"{pdb}_chain_Repair.pdb"
            rounds = [work / f"repair_round_{i}.pdb"
                      for i in range(1, config.repair_iterations + 1)]
            shutil.copy(work / f"{pdb}.pdb", chain)
            for i, round_path in enumerate(rounds, start=1):
                _foldx([binary, "--command=RepairPDB", f"--pdb={chain.name}"], work, log)
                if not chain_out.exists():
                    break
                shutil.copy(chain_out, round_path)
                shutil.move(str(chain_out), chain)
            if rounds[-1].exists():
                shutil.copy(rounds[-1], repaired)
            chain.unlink(missing_ok=True)
            meta["repair_iterations"] = config.repair_iterations
    if not repaired.exists():
        return _finish(config, pdb, meta, "repair_failed")

    if repair_only:
        # Convergence probe: the repair rounds are the experiment. No results JSON is written,
        # so this cannot be mistaken for a cached complex result later.
        meta["n_ok"] = 0
        return ComplexResult(pdb, {}, meta, f"repaired x{config.repair_iterations}")

    # --- resolve the mutation list to author-chain form -------------------------------------
    if mode == MODE_ROLE:
        mapping, unresolved = map_role_to_author(entry, mutations)
        meta["unresolved"] = unresolved
    else:
        mapping = {m: m for m in mutations}
    if not mapping:
        return _finish(config, pdb, meta, "no_mapping")

    # --- validate the wild-type residues against the structure -------------------------------
    residues = repaired_wt_residues(repaired)
    if mode == MODE_VARIANT:
        ordered, rejected = [], []
        for variant in mapping:
            subs = [Mutation.parse(s) for s in variant.split(",")]
            bad = [s for s in subs if residues.get((s.chain, s.position)) != s.wt]
            (rejected if bad else ordered).append(
                f"{variant}({[str(b) for b in bad]})" if bad else (variant, variant)
            )
    else:
        ordered, rejected = validate_against_structure(mapping, residues)
    meta["validation_failed"] = rejected
    if not ordered:
        return _finish(config, pdb, meta, "all_validation_failed")

    # --- BuildModel -------------------------------------------------------------------------
    # One line per mutation, semicolon-terminated. A multi-point variant is already
    # comma-joined, which is exactly the form BuildModel expects for a combined mutant.
    individual_list = format_individual_list(ordered)
    list_path = work / "individual_list.txt"
    # BuildModel outputs are numbered by POSITION in this list, and nothing cleans the work dir.
    # If the list changed since the last run -- the exact reason to recompute a complex -- the
    # stale numbered .pdb/.fxout files still exist and would be read as this run's results,
    # reporting ok(n/n) with the previous list's energies. Clear them whenever the list differs.
    if list_path.exists() and list_path.read_text() != individual_list:
        for pattern in (f"{pdb}_Repair_*.pdb", f"WT_{pdb}_Repair_*.pdb",
                        f"Interaction_*_AC.fxout", "*_AC.fxout"):
            for stale in work.glob(pattern):
                stale.unlink(missing_ok=True)
        meta["cleared_stale_buildmodel_outputs"] = True
    list_path.write_text(individual_list)
    # Record the accepted list. A mutation's energy depends on every entry preceding it, and the
    # list is derived from THIS round's repaired structure -- process_complex drops any mutation
    # whose wild-type residue does not match. So if a repair round ever altered a residue
    # identity, two rounds would silently build from different lists and a between-round
    # comparison would stop isolating repair count. Storing the hash lets `compare` detect that
    # instead of leaving it to be noticed.
    meta["mutation_list_sha256"] = hashlib.sha256(individual_list.encode()).hexdigest()
    meta["n_accepted"] = len(ordered)
    _foldx(
        [binary, "--command=BuildModel", f"--pdb={pdb}_Repair.pdb",
         "--mutant-file=individual_list.txt", "--numberOfRuns=1"],
        work, log,
    )
    n = len(ordered)
    mutant_pdbs = [f"{pdb}_Repair_{i + 1}.pdb" for i in range(n)]
    wildtype_pdbs = [f"WT_{pdb}_Repair_{i + 1}.pdb" for i in range(n)]
    if not all((work / p).exists() for p in mutant_pdbs + wildtype_pdbs):
        return _finish(config, pdb, meta, "build_failed")

    # --- AnalyseComplex on mutants and their paired wild types --------------------------------
    (work / "pdblist.txt").write_text("\n".join(mutant_pdbs + wildtype_pdbs) + "\n")
    _foldx(
        [binary, "--command=AnalyseComplex", "--pdb-list=pdblist.txt",
         f"--analyseComplexChains={entry.groups}"],
        work, log,
    )

    results: Dict[str, dict] = {}
    for i, (key, author) in enumerate(ordered):
        mutant_file = work / f"Interaction_{pdb}_Repair_{i + 1}_AC.fxout"
        wildtype_file = work / f"Interaction_WT_{pdb}_Repair_{i + 1}_AC.fxout"
        if not (mutant_file.exists() and wildtype_file.exists()):
            continue
        mutant = parse_interaction_file(mutant_file, terms)
        wildtype = parse_interaction_file(wildtype_file, terms)
        record = {
            t: round(mutant[t] - wildtype[t], 4)
            for t in terms
            if t in mutant and t in wildtype
        }
        record["cleaned"] = author
        results[key] = record

    meta["n_ok"] = len(results)
    meta["groups"] = entry.groups
    result = ComplexResult(pdb, results, meta, f"ok({len(results)}/{len(mutations)})")
    out.write_text(json.dumps(result.to_json(), indent=1))
    return result


def _task(args):
    pdb, mutations, entry, config, mode, repair_only = args
    try:
        return process_complex(pdb, mutations, entry, config, mode,
                               repair_only=repair_only).status
    except Exception as exc:  # a single bad complex must not abandon the campaign
        return f"EXC {exc}"


def run_campaign(
    worklist: Dict[str, Sequence[str]],
    skempi: Dict[str, SkempiComplex],
    config: FoldxConfig,
    mode: str = MODE_ROLE,
    jobs: int = 1,
    on_result=print,
    repair_only: bool = False,
) -> Dict[str, str]:
    """Run a whole worklist, optionally across processes.

    ``worklist`` maps a record name to the mutation strings to compute for it. A campaign is
    exactly a choice of worklist and ``mode``; nothing else distinguishes one from another.

    ``skempi`` may be keyed by SKEMPI identifier -- as :func:`~skempi_foldx.load_skempi` returns
    -- or by bare PDB code, and a worklist key is matched against either. Both shapes are real:
    :func:`worklist_single_point` yields identifiers, while :func:`worklist_from_table` reads a
    dataset table whose row labels only carry the code. Requiring the two to agree made a
    mismatched pair resolve nothing and the campaign report "0 complexes, 0 mutations" on its way
    to exiting successfully.
    """
    config.ensure_dirs()
    config.resolve_binary()   # fail fast, before spawning workers

    # Before anything is enqueued. An intractable complex must never reach a worker: killing the
    # hung job is a per-round remedy, and this function is called once per round.
    worklist = filter_worklist(worklist, on_note=on_result, context="worklist")

    # Both spellings of every entry, so a worklist keyed either way resolves. A code that SKEMPI
    # defines twice is deliberately absent from the code-level index: there is no single entry it
    # could mean, and picking one would score half the list against an interface it is not part
    # of -- the defect the identifier keying exists to remove.
    by_code: Dict[str, list] = {}
    for entry in skempi.values():
        by_code.setdefault(entry.pdb, []).append(entry)
    resolvable = dict(skempi)
    for code, entries in by_code.items():
        if len(entries) == 1 and code not in resolvable:
            resolvable[code] = entries[0]
    # And the other direction, for a caller whose `skempi` is keyed by code: an identifier-keyed
    # worklist against it would otherwise resolve nothing and report every complex absent.
    for entry in skempi.values():
        resolvable.setdefault(entry.id, entry)

    # One record, one run. A worklist naming the same definition twice -- once by code and once by
    # identifier -- would otherwise process it twice into the same work directory and results file.
    #
    # The two lists are POOLED, not one dropped. Dropping is the worse failure by far: this
    # module's whole subject is that a mutation's energy depends on the list it was computed in,
    # so silently computing a subset of what was asked for produces values that are wrong by up
    # to 11 kcal/mol and look fine. The identifier wins the name, because the code is the spelling
    # that cannot say which pairing it means.
    merged: Dict[str, List[str]] = {}
    winner: Dict[str, str] = {}
    collapsed = []
    for key in sorted(worklist, key=lambda k: (k == code_of(k), k)):
        entry = resolvable.get(key)
        if entry is None:
            merged.setdefault(key, list(worklist[key]))
            continue
        name = winner.setdefault(entry.id, key)
        if name != key:
            collapsed.append((name, key))
        pooled = merged.setdefault(name, [])
        pooled.extend(m for m in worklist[key] if m not in pooled)
    if collapsed:
        on_result(f"[foldx] {len(collapsed)} worklist keys name a record already queued under "
                  f"another name, e.g. {collapsed[:3]}; their mutations are pooled into one run.")
    worklist = {k: sorted(v) if v else [] for k, v in merged.items()}

    # Which record names share a structure, so the second definition of a code reuses the first's
    # repair instead of paying for it again. Derived here rather than asked of the caller: this is
    # the one place that knows both the worklist's spelling and the table behind it, and a facility
    # no driver populates is a facility that never runs.
    if config.sibling_records is None:
        shared: Dict[str, List[str]] = {}
        for key in worklist:
            entry = resolvable.get(key)
            if entry is not None:
                shared.setdefault(entry.pdb, []).append(key)
        # A copy, not an assignment: writing it back would leave a caller's config carrying one
        # campaign's sibling map into the next, and `{}` is falsy but not None, so a first
        # campaign with no siblings would disable reuse for every later one.
        config = replace(config, sibling_records={k: v for k, v in shared.items() if len(v) > 1})

    # An empty list means "repair this and stop", which is only meaningful under repair_only.
    # Anywhere else it would pay a full RepairPDB and then write an error JSON that marks the
    # complex done for every later resume and scope enumeration.
    if not repair_only:
        empty = sorted(k for k, v in worklist.items() if not v)
        if empty:
            on_result(f"[foldx] {len(empty)} complexes have no mutations to build and are not a "
                      f"repair-only run, skipped: {', '.join(empty[:8])}")
            worklist = {k: v for k, v in worklist.items() if v}

    targets = [(pdb, list(muts), resolvable[pdb], config, mode, repair_only)
               for pdb, muts in sorted(worklist.items()) if pdb in resolvable]
    skipped = sorted(set(worklist) - set(resolvable))
    if skipped:
        doubled = sorted(k for k in skipped if len(by_code.get(k, ())) > 1)
        on_result(f"[foldx] {len(skipped)} complexes not resolved, skipped: "
                  f"{', '.join(skipped[:8])}{'...' if len(skipped) > 8 else ''}")
        if doubled:
            on_result(f"[foldx] of those, {len(doubled)} are in SKEMPI but named by a code it "
                      f"defines more than once ({', '.join(doubled)}); there is no single entry "
                      f"they could mean. Key the worklist by identifier, or pass entries keyed "
                      f"by code -- see skempi_foldx.by_pdb.")

    on_result(f"[foldx] {len(targets)} complexes, "
              f"{sum(len(t[1]) for t in targets)} mutations, jobs={jobs}")

    statuses: Dict[str, str] = {}
    if jobs == 1:
        for args in targets:
            statuses[args[0]] = _task(args)
            on_result(f"[foldx] {args[0]}: {statuses[args[0]]}")
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_task, a): a[0] for a in targets}
            for future in as_completed(futures):
                pdb = futures[future]
                statuses[pdb] = future.result()
                on_result(f"[foldx] {pdb}: {statuses[pdb]}")
    on_result("[foldx] done")
    return statuses


# --------------------------------------------------------------------------------------------
# Worklist builders — one per campaign shape.
# --------------------------------------------------------------------------------------------

def worklist_from_table(table_path: Path) -> Dict[str, List[str]]:
    """Single-point role mutations from a tab-separated dataset mutation table (``MODE_ROLE``).

    **This produces a per-dataset SUBSET of each complex's mutations, which makes the computed
    energies dataset-dependent** — see "Determinism" in the module docstring. It is kept for
    backward compatibility with campaigns that were driven from such a table. For new compute use
    :func:`worklist_single_point`, which takes the union, so a mutation gets the same value
    regardless of which dataset wanted it.

    Column 1 is the row label and column 3 the role-chain mutation; the PDB id is the part of the
    label before ``_``. Rows whose mutation field contains a comma are multi-point and skipped.
    """
    by_pdb: Dict[str, set] = {}
    with open(table_path) as fh:
        for line in fh:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or "," in fields[2]:
                continue
            by_pdb.setdefault(fields[0].split("_")[0], set()).add(fields[2])
    return {pdb: sorted(muts) for pdb, muts in by_pdb.items()}


def _selected(ident: str, entry: SkempiComplex, wanted: Optional[set]) -> bool:
    """Whether ``only=`` names this definition, by identifier **or** by PDB code.

    A caller filtering a worklist usually holds codes — a list of structures to run — while the
    worklist is keyed by identifier. Matching only on the key made ``only=["1KBH"]`` select
    nothing and the campaign exit successfully having done nothing, which is the failure mode the
    scope guard exists to prevent and could not see.
    """
    return wanted is None or ident in wanted or entry.pdb in wanted


def pool_by_code(worklist: Dict[str, Sequence[str]]) -> Dict[str, List[str]]:
    """Re-key a per-definition worklist by PDB code, pooling the definitions of each code.

    For the campaigns that are genuinely per *structure* rather than per pairing — the repair
    sweep, the repair-count ablation, the residue check — where the unit of work is one repaired
    PDB and one results file named after it.

    **This reintroduces the pooling the store exists to avoid, and that is the point of naming it
    explicitly.** A code SKEMPI defines twice gets one list holding both definitions' mutations,
    computed against whichever pairing the entry carries, so its values belong to no single
    pairing's list. That is acceptable for a convergence study, which compares a complex against
    itself across repair rounds and never joins the result onto a SKEMPI row. It is not acceptable
    for the shipped store, which is why the three affected codes are computed by
    ``recompute_alternate_interfaces.py --all-definitions`` instead.
    """
    out: Dict[str, List[str]] = {}
    for key, muts in worklist.items():
        # A bare string satisfies `Sequence[str]` and would be extended character by character,
        # turning one mutation into a list of letters. Cheap to refuse, expensive to debug.
        if isinstance(muts, str):
            raise TypeError(f"worklist[{key!r}] is a string, not a sequence of mutation strings")
        out.setdefault(code_of(key), []).extend(muts)
    return {pdb: sorted(set(muts)) for pdb, muts in out.items()}


def worklist_single_point(
    skempi: Dict[str, SkempiComplex], only: Optional[Iterable[str]] = None
) -> Dict[str, List[str]]:
    """Every single-point cleaned mutation SKEMPI knows about, per definition (``MODE_AUTHOR``).

    Keyed by SKEMPI identifier, following :func:`~skempi_foldx.load_skempi`. Each key is one
    campaign: a complex SKEMPI defines twice yields two, and they must be run into separate
    results directories, since a campaign writes one file per key and the pairing is what
    distinguishes them.

    The union, deliberately: because BuildModel's result for an entry depends on the entries
    before it, computing the complete list once is what makes a value independent of which
    dataset asked for it. Sorting makes the list — and therefore every energy — reproducible.
    """
    wanted = set(only) if only is not None else None
    return {
        ident: sorted(entry.single)
        for ident, entry in skempi.items()
        if entry.single and _selected(ident, entry, wanted)
    }


def worklist_multi_point(
    skempi: Dict[str, SkempiComplex], only: Optional[Iterable[str]] = None
) -> Dict[str, List[str]]:
    """Every multi-point variant (``MODE_VARIANT``)."""
    wanted = set(only) if only is not None else None
    return {
        ident: sorted(entry.multi)
        for ident, entry in skempi.items()
        if entry.multi and _selected(ident, entry, wanted)
    }


def exclude_already_computed(
    worklist: Dict[str, Sequence[str]], results_dirs: Iterable[Path]
) -> Dict[str, Sequence[str]]:
    """Drop complexes already covered by any of ``results_dirs``."""
    done = {p.stem for d in results_dirs for p in Path(d).glob("*.json")}
    return {pdb: muts for pdb, muts in worklist.items() if pdb not in done}
