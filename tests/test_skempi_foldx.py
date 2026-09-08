"""The contracts that, if broken, fail silently.

FoldX failures mostly do not raise. A wrong subtraction direction, a mutation paired with another
mutation's energy file, a misread SKEMPI column or a swapped chain group all produce a
correctly-shaped result with the wrong numbers in it, and every count downstream still reports
success. So the suite pins the arithmetic and the pairing by value, not just the shapes.

The numeric core is covered end to end from fixture `.fxout` files, asserted per mutation by name;
that block is verified by mutation testing -- flipping the sign, shifting the index, doubling the
parse, reading the neighbouring SKEMPI column, swapping the chain groups, or reversing the term
order each make it fail.

No FoldX binary, network, or working directory outside tmp_path is required.
"""

from __future__ import annotations

import json
import pathlib
import warnings

import pytest

from skempi_foldx import consolidate, load_store, source_label
from skempi_foldx.skempi import Mutation, SkempiComplex, map_role_to_author
from skempi_foldx.terms import N_TERMS, SCALAR_TERM, TERMS


# ============================================================================================
# The term contract
# ============================================================================================

def test_term_order_is_pinned():
    """Column order is consumed positionally by any consumer that flattens a record into a
    feature vector. Reordering silently permutes the features, and what is fit on them looks
    plausible and is wrong."""
    assert TERMS == (
        "Interaction Energy", "Backbone Hbond", "Sidechain Hbond", "Van der Waals",
        "Electrostatics", "Solvation Polar", "Solvation Hydrophobic",
        "Van der Waals clashes", "entropy sidechain", "entropy mainchain",
        "torsional clash", "backbone clash",
    )
    # A tuple rather than a list, and that is the point: TERMS is public and is the default
    # argument of term_vector(), so a list would let any importer permute every feature vector
    # produced afterwards -- process-wide, silently, and exactly the failure this module exists
    # to prevent.
    assert isinstance(TERMS, tuple)
    with pytest.raises((AttributeError, TypeError)):
        TERMS.reverse()          # type: ignore[attr-defined]


# ============================================================================================
# Standardization
# ============================================================================================

def test_role_maps_to_author_chain():
    entry = SkempiComplex("1XXX", "E", "I", single={"LI38S"})
    mapping, unresolved = map_role_to_author(entry, ["LB38S"])
    assert mapping == {"LB38S": "LI38S"} and unresolved == []


def test_ambiguous_substitution_is_disambiguated_by_role():
    """The same substitution on both partners resolves via the chain group."""
    entry = SkempiComplex("1XXX", "E", "I", single={"LE38S", "LI38S"})
    mapping, _ = map_role_to_author(entry, ["LA38S", "LB38S"])
    assert mapping == {"LA38S": "LE38S", "LB38S": "LI38S"}


def test_unresolvable_mutation_is_reported_not_guessed():
    entry = SkempiComplex("1XXX", "EF", "I", single={"LE38S", "LF38S"})
    mapping, unresolved = map_role_to_author(entry, ["LA38S"])
    assert mapping == {} and unresolved == ["LA38S"]   # both candidates in group 1


def test_mutation_roundtrip():
    m = Mutation.parse("LI38S")
    assert (m.wt, m.chain, m.position, m.mutant) == ("L", "I", 38, "S")
    assert str(m) == "LI38S" and m.identity == ("L", 38, "S")


# ============================================================================================
# Store consolidation
# ============================================================================================

def _write_store(directory, records):
    directory.mkdir(parents=True, exist_ok=True)
    for pdb, muts in records.items():
        (directory / f"{pdb}.json").write_text(json.dumps({"muts": muts, "meta": {}}))


def test_consolidate_unions_per_mutation_not_per_file(tmp_path):
    """Picking one whole file per complex is lossless only if that file is a superset of the
    others -- a property a key-set audit cannot check."""
    a, b = tmp_path / "a", tmp_path / "b"
    _write_store(a, {"1XXX": {"M1": {t: 1.0 for t in TERMS}}})
    _write_store(b, {"1XXX": {"M2": {t: 2.0 for t in TERMS}}})
    store, report = consolidate([a, b])
    assert set(store["1XXX"]) == {"M1", "M2"}
    assert report.mutations == 2


def test_consolidate_precedence_and_conflict_reporting(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    _write_store(a, {"1XXX": {"M1": {t: 1.0 for t in TERMS}}})
    _write_store(b, {"1XXX": {"M1": {t: 9.0 for t in TERMS}}})
    store, report = consolidate([a, b])
    assert store["1XXX"]["M1"]["Interaction Energy"] == 1.0     # first source wins
    assert len(report.conflicts) == 1
    assert report.conflicts[0]["max_abs_delta"] == pytest.approx(8.0)
    # Labels are parent-qualified so two campaigns named `results` cannot collide.
    assert report.redundant_sources == [source_label(b)]        # contributed nothing new


def test_written_store_is_self_contained(tmp_path):
    """Real files, not symlinks: pruning a source must not break the canonical store."""
    a, dest = tmp_path / "a", tmp_path / "dest"
    _write_store(a, {"1XXX": {"M1": {t: 1.0 for t in TERMS}}})
    consolidate([a], destination=dest)
    import shutil
    shutil.rmtree(a)
    assert load_store(dest)["1XXX"]["M1"]["Interaction Energy"] == 1.0


# ============================================================================================
# Merge, end to end on a fixture
# ============================================================================================

def test_repair_iterations_defaults_to_one_and_is_guarded(tmp_path):
    """One repair is FoldX's own recommendation and matches CATH-ddG. It is configurable so the
    question can be measured, not so it can drift."""
    from skempi_foldx import FoldxConfig
    kw = dict(results_dir=tmp_path/"r", work_dir=tmp_path/"w",
              pdb_dir=tmp_path, skempi_csv=tmp_path/"s.csv")
    assert FoldxConfig(**kw).repair_iterations == 1
    assert FoldxConfig(**kw, repair_iterations=5).repair_iterations == 5
    with pytest.raises(ValueError, match=">= 1"):
        FoldxConfig(**kw, repair_iterations=0)


def test_iterated_repair_does_not_reuse_a_seeded_structure(tmp_path):
    """Seeding a sibling campaign's single-repair structure would silently defeat the ablation:
    the run would report N iterations while using a once-repaired input."""
    from skempi_foldx import FoldxConfig
    seed = tmp_path/"seed"/"1ABC"
    seed.mkdir(parents=True)
    (seed/"1ABC_Repair.pdb").write_text("ATOM\n")
    cfg = FoldxConfig(results_dir=tmp_path/"r", work_dir=tmp_path/"w", pdb_dir=tmp_path,
                      skempi_csv=tmp_path/"s.csv", repair_seed_dirs=[tmp_path/"seed"])
    assert cfg.find_repaired("1ABC") is not None      # available for the 1x path

    # Drive the guard rather than the config: at repair_iterations > 1 the seed must be ignored,
    # or the run reports N iterations while building from a once-repaired structure.
    from skempi_foldx import process_complex
    from skempi_foldx.skempi import SkempiComplex

    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    (tmp_path / "1ABC.pdb").write_text(
        "ATOM      1  CA  LEU A   1       0.0   0.0   0.0  1.00  0.00           C\n")
    cfg5 = FoldxConfig(results_dir=tmp_path/"r5", work_dir=tmp_path/"w5", pdb_dir=tmp_path,
                       skempi_csv=tmp_path/"s.csv", repair_seed_dirs=[tmp_path/"seed"],
                       repair_iterations=5, binary=binary)
    cfg5.ensure_dirs()
    res = process_complex("1ABC", [], SkempiComplex("1ABC", "A", "B"), cfg5, repair_only=True)
    assert "repair_seeded_from" not in res.meta, \
        "a 5x run reused a 1x seed -- the ablation would compare a structure against itself"

    # and the mirror: at 1x the seed IS the point, so it must be used
    cfg1 = FoldxConfig(results_dir=tmp_path/"r1", work_dir=tmp_path/"w1", pdb_dir=tmp_path,
                       skempi_csv=tmp_path/"s.csv", repair_seed_dirs=[tmp_path/"seed"],
                       repair_iterations=1, binary=binary)
    cfg1.ensure_dirs()
    res1 = process_complex("1ABC", [], SkempiComplex("1ABC", "A", "B"), cfg1, repair_only=True)
    assert "repair_seeded_from" in res1.meta, "a 1x run failed to reuse an available seed"


def test_ablation_rmsd_is_unsuperposed_and_correct(tmp_path):
    """RepairPDB does not move the backbone frame, so the quantity of interest is in-place
    displacement. Superposing would mask exactly that."""
    import sys
    from pathlib import Path as _P
    # absolute: a relative insert makes these tests pass only from the repo root
    sys.path.insert(0, str(_P(__file__).resolve().parent.parent / "experiments"))
    from repair_ablation import heavy_atoms, rmsd

    a, b = tmp_path/"a.pdb", tmp_path/"b.pdb"
    a.write_text("ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n")
    b.write_text("ATOM      1  CA  ALA A   1       3.000   4.000   0.000  1.00  0.00           C\n")
    assert rmsd(a, b)[0] == pytest.approx(5.0)
    assert rmsd(a, a)[0] == pytest.approx(0.0)
    assert len(heavy_atoms(a)) == 1


def test_ablation_rmsd_ignores_hydrogens(tmp_path):
    import sys
    from pathlib import Path as _P
    # absolute: a relative insert makes these tests pass only from the repo root
    sys.path.insert(0, str(_P(__file__).resolve().parent.parent / "experiments"))
    from repair_ablation import heavy_atoms
    p = tmp_path/"h.pdb"
    p.write_text(
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00  0.00           C\n"
        "ATOM      2  HA  ALA A   1       1.000   0.000   0.000  1.00  0.00           H\n")
    assert len(heavy_atoms(p)) == 1


def test_single_and_multi_point_results_must_not_share_a_directory(tmp_path):
    """A single-point results file read in multi-point mode holds no records, so a resume check
    that only looks for its own key reports "cached" and skips a whole arm in silence."""
    import json
    from skempi_foldx import FoldxConfig, MODE_VARIANT
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = FoldxConfig(results_dir=tmp_path/"res", work_dir=tmp_path/"work",
                      pdb_dir=tmp_path, skempi_csv=tmp_path/"s.csv")
    cfg.ensure_dirs()
    (cfg.results_dir/"1XXX.json").write_text(json.dumps(
        {"muts": {"LA1A": {"Interaction Energy": 1.0}}, "meta": {}}))
    entry = SkempiComplex("1XXX", "A", "B", single={"LA1A"}, multi={"LA1A,LB2B"})
    with pytest.raises(RuntimeError, match="separate results_dir"):
        process_complex("1XXX", ["LA1A,LB2B"], entry, cfg, mode=MODE_VARIANT)


def test_matching_mode_still_resumes_from_cache(tmp_path):
    import json
    from skempi_foldx import FoldxConfig, MODE_AUTHOR
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = FoldxConfig(results_dir=tmp_path/"res", work_dir=tmp_path/"work",
                      pdb_dir=tmp_path, skempi_csv=tmp_path/"s.csv")
    cfg.ensure_dirs()
    (cfg.results_dir/"1XXX.json").write_text(json.dumps(
        {"muts": {"LA1A": {"Interaction Energy": 1.0}}, "meta": {}}))
    entry = SkempiComplex("1XXX", "A", "B", single={"LA1A"})
    r = process_complex("1XXX", ["LA1A"], entry, cfg, mode=MODE_AUTHOR)
    assert r.status == "cached" and "LA1A" in r.mutations


def test_recorded_failure_still_resumes(tmp_path):
    """An empty payload is a recorded failure, not a mode collision — must not raise."""
    import json
    from skempi_foldx import FoldxConfig, MODE_VARIANT
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = FoldxConfig(results_dir=tmp_path/"res", work_dir=tmp_path/"work",
                      pdb_dir=tmp_path, skempi_csv=tmp_path/"s.csv")
    cfg.ensure_dirs()
    (cfg.results_dir/"1XXX.json").write_text(json.dumps(
        {"variants": {}, "meta": {"error": "no_pdb"}}))
    entry = SkempiComplex("1XXX", "A", "B", multi={"LA1A,LB2B"})
    r = process_complex("1XXX", ["LA1A,LB2B"], entry, cfg, mode=MODE_VARIANT)
    assert r.status == "cached" and r.mutations == {}


def test_consolidate_refuses_to_mix_single_and_multi_point(tmp_path):
    """SP and MP come from separate BuildModel invocations with separate mutation lists. Merging
    the stores would imply they are entries of one list and erase which arm each came from."""
    import json
    sp, mp = tmp_path/"results_sp", tmp_path/"results_mp"
    sp.mkdir(); mp.mkdir()
    (sp/"1XXX.json").write_text(json.dumps({"muts": {"LA1A": {t: 1.0 for t in TERMS}}, "meta": {}}))
    (mp/"1XXX.json").write_text(json.dumps({"variants": {"LA1A,LB2B": {t: 9.0 for t in TERMS}}, "meta": {}}))
    with pytest.raises(ValueError, match="separate BuildModel campaigns"):
        consolidate([sp, mp])
    # each arm on its own is fine
    assert consolidate([sp])[1].mutations == 1
    assert consolidate([mp])[1].mutations == 1


def test_mutation_list_hash_detects_a_changed_list(tmp_path):
    """The accepted mutation list is derived from each round's repaired structure, so a repair
    round that altered a residue identity would silently change it — and a between-round
    comparison would stop isolating repair count."""
    import json, sys
    from pathlib import Path as _P
    sys.path.insert(0, str(_P(__file__).resolve().parent.parent / "experiments"))
    from repair_ablation import list_hashes

    a, b = tmp_path/"a", tmp_path/"b"
    a.mkdir(); b.mkdir()
    (a/"1XXX.json").write_text(json.dumps({"muts": {}, "meta": {"mutation_list_sha256": "aaa"}}))
    (b/"1XXX.json").write_text(json.dumps({"muts": {}, "meta": {"mutation_list_sha256": "bbb"}}))
    (a/"2YYY.json").write_text(json.dumps({"muts": {}, "meta": {"mutation_list_sha256": "ccc"}}))
    (b/"2YYY.json").write_text(json.dumps({"muts": {}, "meta": {"mutation_list_sha256": "ccc"}}))
    (a/"3ZZZ.json").write_text(json.dumps({"muts": {}, "meta": {}}))   # predates the hash

    ha, hb = list_hashes(a), list_hashes(b)
    shared = set(ha) & set(hb)
    assert sorted(p for p in shared if ha[p] != hb[p]) == ["1XXX"]
    assert "3ZZZ" not in ha, "a complex without a recorded hash must not be reported as verified"


def test_file_sha256_is_accepted_but_list_sha256_is_not(tmp_path):
    """The sweep tooling records two hashes per complex under different normalizations.
    `file_sha256` is the canonical normalization (raw individual_list.txt bytes) and is an
    accepted alias. `list_sha256` strips ';' and joins with '\\n' — treating it as equivalent
    would manufacture mismatches indistinguishable from real list drift."""
    import json, sys
    from pathlib import Path as _P
    sys.path.insert(0, str(_P(__file__).resolve().parent.parent / "experiments"))
    from repair_ablation import list_hashes

    d = tmp_path/"s"; d.mkdir()
    (d/"1XXX.json").write_text(json.dumps({"muts": {}, "meta": {"file_sha256": "abc"}}))
    (d/"2YYY.json").write_text(json.dumps({"muts": {}, "meta": {"list_sha256": "def"}}))
    (d/"3ZZZ.json").write_text(json.dumps(
        {"muts": {}, "meta": {"mutation_list_sha256": "ours", "file_sha256": "abc"}}))
    h = list_hashes(d)
    assert h["1XXX"] == "abc", "file_sha256 should be accepted"
    assert "2YYY" not in h, "list_sha256 uses a different normalization and must be ignored"
    assert h["3ZZZ"] == "ours", "mutation_list_sha256 takes precedence when both are present"


def test_individual_list_uses_author_chains_not_role_keys():
    """meta["mutation_list_sha256"] hashes these bytes. Rebuilding that hash from the role->author
    mapping KEYS instead would hash role-chain strings, so every single-point complex in a
    MODE_ROLE campaign would mismatch a store hashed from the file — indistinguishable from real
    list drift."""
    import hashlib
    from skempi_foldx.run import format_individual_list

    ordered = [("LB38S", "LI38S"), ("GB32Y", "GI32Y")]   # (role key, author) — they differ
    text = format_individual_list(ordered)
    assert text == "LI38S;\nGI32Y;\n", "must write the AUTHOR form, one per line, ';'-terminated"
    assert "LB38S" not in text and "GB32Y" not in text, "role keys must never reach the file"
    # the pinned digest of the canonical convention
    assert hashlib.sha256(text.encode()).hexdigest().startswith("ba5f85e1")


def test_individual_list_is_order_preserving():
    """Order is the whole point: an entry's energy depends on every entry before it."""
    from skempi_foldx.run import format_individual_list
    a = format_individual_list([("x", "AAA"), ("y", "BBB")])
    b = format_individual_list([("y", "BBB"), ("x", "AAA")])
    assert a != b and a == "AAA;\nBBB;\n"


def _fake_foldx(tmp_path):
    """A stand-in FoldX that appends one mark per repair, so rounds are countable."""
    import stat
    fx = tmp_path / "foldx"
    fx.write_text('#!/bin/sh\nfor a in "$@"; do case "$a" in --pdb=*) P=${a#--pdb=};; esac; done\n'
                  'B=$(basename "$P" .pdb)\ncp "$P" "${B}_Repair.pdb"\n'
                  'echo MARK >> "${B}_Repair.pdb"\n')
    fx.chmod(fx.stat().st_mode | stat.S_IEXEC)
    return fx


def _repair_setup(tmp_path, iterations):
    from skempi_foldx import FoldxConfig
    pdbs = tmp_path / "pdbs"; pdbs.mkdir()
    # A real CA record for LEU A1, so `repaired_wt_residues` can validate mutation LA1A.
    (pdbs / "1ABC.pdb").write_text(
        "ATOM      1  CA  LEU A   1       0.000   0.000   0.000  1.00  0.00           C\n")
    cfg = FoldxConfig(results_dir=tmp_path/"res", work_dir=tmp_path/"work", pdb_dir=pdbs,
                      skempi_csv=tmp_path/"s.csv", binary=_fake_foldx(tmp_path),
                      repair_iterations=iterations)
    cfg.ensure_dirs()
    cfg.complex_work_dir("1ABC").mkdir(parents=True)
    return cfg


def _marks(p):
    return p.read_text().count("MARK") if p.exists() else None


def test_iterated_repair_never_mutates_the_input_structure(tmp_path):
    from skempi_foldx import MODE_AUTHOR
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = _repair_setup(tmp_path, 5)
    w = cfg.complex_work_dir("1ABC")
    process_complex("1ABC", ["LA1A"], SkempiComplex("1ABC", "A", "B", single={"LA1A"}),
                    cfg, mode=MODE_AUTHOR, repair_only=True)
    assert [_marks(w / f"repair_round_{i}.pdb") for i in range(1, 6)] == [1, 2, 3, 4, 5]
    assert _marks(w / "1ABC.pdb") == 0, "the input structure must never be repaired in place"


def test_interrupted_repair_chain_restarts_instead_of_compounding(tmp_path):
    """The failure this guards: a crash between rounds leaves a partially-repaired structure in
    `<pdb>.pdb`, the next run copies that over `_original.pdb` — destroying the record of the true
    input — and repairs N more times from it. A 5-round resume after a round-3 crash then yields
    an 8x structure reported as "repaired x5", with the saved rounds holding 4x-8x."""
    from skempi_foldx import MODE_AUTHOR
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = _repair_setup(tmp_path, 5)
    w = cfg.complex_work_dir("1ABC")
    entry = SkempiComplex("1ABC", "A", "B", single={"LA1A"})
    process_complex("1ABC", ["LA1A"], entry, cfg, mode=MODE_AUTHOR, repair_only=True)

    # partial chain + a stale mutable working file, i.e. killed mid-loop
    for i in (4, 5):
        (w / f"repair_round_{i}.pdb").unlink()
    (w / "1ABC_chain.pdb").write_text("ATOM\n" + "MARK\n" * 3)

    result = process_complex("1ABC", ["LA1A"], entry, cfg, mode=MODE_AUTHOR, repair_only=True)
    assert result.meta.get("restarted_incomplete_repair_chain") is True
    assert [_marks(w / f"repair_round_{i}.pdb") for i in range(1, 6)] == [1, 2, 3, 4, 5], \
        "rounds must be 1x..5x, not a continuation of the interrupted chain"
    assert _marks(w / "1ABC.pdb") == 0


def test_write_store_preserves_the_multi_point_key(tmp_path):
    """A write_store that always emits `muts` relabels a consolidated multi-point store's variants
    as single-point — and the SP/MP guard then cannot fire on the stores this module itself
    produces."""
    import json
    from skempi_foldx.store import consolidate, store_kind

    src, dest = tmp_path / "mp", tmp_path / "out"
    src.mkdir()
    (src / "1XXX.json").write_text(json.dumps(
        {"variants": {"LA1A,LB2B": {t: 1.0 for t in TERMS}},
         "meta": {"mutation_list_sha256": "abc123"}}))
    consolidate([src], destination=dest)
    assert store_kind(dest) == "variants"
    payload = json.loads((dest / "1XXX.json").read_text())
    assert "muts" not in payload and "LA1A,LB2B" in payload["variants"]
    # provenance must survive: the repair-count comparison depends on this hash
    assert payload["meta"].get("mutation_list_sha256") == "abc123"


def test_consolidating_the_two_arms_still_refuses_after_a_round_trip(tmp_path):
    import json
    from skempi_foldx.store import consolidate

    sp, mp = tmp_path / "sp", tmp_path / "mp"
    for d, key, k in ((sp, "muts", "LA1A"), (mp, "variants", "LA1A,LB2B")):
        d.mkdir()
        (d / "1XXX.json").write_text(json.dumps({key: {k: {t: 1.0 for t in TERMS}}, "meta": {}}))
    consolidate([sp], destination=tmp_path / "sp2")
    consolidate([mp], destination=tmp_path / "mp2")
    with pytest.raises(ValueError, match="separate BuildModel campaigns"):
        consolidate([tmp_path / "sp2", tmp_path / "mp2"])


def test_stale_buildmodel_outputs_are_cleared_when_the_list_changes(tmp_path):
    """BuildModel outputs are numbered by POSITION in the mutation list and nothing cleans the
    work dir, so recomputing a complex with a corrected list reads the previous list's energies
    and reports ok(n/n) unless the stale outputs are cleared."""
    import stat, json
    from skempi_foldx import MODE_AUTHOR
    from skempi_foldx.run import process_complex
    from skempi_foldx.skempi import SkempiComplex

    cfg = _repair_setup(tmp_path, 1)
    w = cfg.complex_work_dir("1ABC")
    (w / "individual_list.txt").write_text("OLD1;\n")
    stale = w / "1ABC_Repair_1.pdb"; stale.write_text("stale")
    (w / "Interaction_1ABC_Repair_1_AC.fxout").write_text("stale")
    entry = SkempiComplex("1ABC", "A", "B", single={"LA1A"})
    result = process_complex("1ABC", ["LA1A"], entry, cfg, mode=MODE_AUTHOR)
    assert result.meta.get("cleared_stale_buildmodel_outputs") is True
    assert not stale.exists(), "outputs from the previous mutation list must not be reused"


# ------------------------------------------------------------------- intractable exclusions
def test_run_campaign_never_enqueues_an_intractable_complex(tmp_path, monkeypatch, capsys):
    """The failure this guards: RepairPDB on 1KBH does not terminate, and the driver is a loop,
    so killing the hung job only postpones it to the next round."""
    from skempi_foldx import FoldxConfig, run_campaign
    from skempi_foldx.exclusions import INTRACTABLE
    from skempi_foldx.skempi import SkempiComplex

    pdb = next(iter(INTRACTABLE))
    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 1\n")
    binary.chmod(0o755)
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv", binary=binary)
    (tmp_path / "p").mkdir()
    (tmp_path / "p" / f"{pdb}.pdb").write_text("ATOM      1  CA  LEU A   1       0.0   0.0   0.0\n")
    entry = SkempiComplex(pdb, "A", "B", single={"LA930W"})

    statuses = run_campaign({pdb: ["LA930W"]}, {pdb: entry}, config, jobs=1)
    assert statuses == {}, "an excluded complex reached a worker"
    assert not (config.results_dir / f"{pdb}.json").exists(), \
        "an excluded complex must leave no store entry -- that is how it gets re-enumerated"
    assert pdb in capsys.readouterr().out, "the exclusion must be reported, not silent"


def test_process_complex_refuses_directly_and_writes_nothing(tmp_path):
    """Second guard, for a driver that builds its own targets (the `sweep_mp.py` shape)."""
    from skempi_foldx import FoldxConfig, process_complex
    from skempi_foldx.exclusions import INTRACTABLE
    from skempi_foldx.skempi import SkempiComplex

    pdb = next(iter(INTRACTABLE))
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv")
    config.ensure_dirs()
    result = process_complex(pdb, ["LA930W"],
                             SkempiComplex(pdb, "A", "B", single={"LA930W"}),
                             config)
    assert result.status == "excluded"
    assert result.mutations == {}
    assert not (config.results_dir / f"{pdb}.json").exists()


def test_the_exclusion_is_announced_as_a_warning_not_only_a_print():
    """A print alone is not loud enough: drivers redirect stdout, and a reader who never sees
    the notice never learns the bypass exists."""
    from skempi_foldx.exclusions import ALLOW_ENV, INTRACTABLE, filter_complexes

    pdb = next(iter(INTRACTABLE))
    with pytest.warns(RuntimeWarning, match="EXCLUDED") as record:
        filter_complexes([pdb], on_note=None)
    assert ALLOW_ENV in str(record[0].message), "the notice must say how to re-enable it"


def test_exclusion_is_bypassable_but_only_deliberately(monkeypatch):
    from skempi_foldx.exclusions import ALLOW_ENV, INTRACTABLE, filter_complexes, is_excluded

    pdb = next(iter(INTRACTABLE))
    assert is_excluded(pdb) is not None
    assert filter_complexes([pdb, "1BRS"], on_note=None) == ["1BRS"]
    monkeypatch.setenv(ALLOW_ENV, "1")
    assert is_excluded(pdb) is None
    assert filter_complexes([pdb, "1BRS"], on_note=None) == [pdb, "1BRS"]


def test_shipped_store_holds_no_excluded_complex():
    """A stale JSON in the store is the entry point: scope enumerated from results rather than
    from the curated table picks the complex back up."""
    from pathlib import Path

    from skempi_foldx.exclusions import INTRACTABLE

    root = Path(__file__).resolve().parents[1]
    for store in ("skempi_foldx/data/results_sp", "skempi_foldx/data/results_mp"):
        for pdb in INTRACTABLE:
            assert not (root / store / f"{pdb}.json").exists(), \
                f"{store}/{pdb}.json re-arms the compute loop for an intractable complex"


def test_write_store_without_a_kind_infers_it_rather_than_assuming_single_point():
    """`write_store` is public. A default of `muts` would mislabel every store written by a caller
    that omits the argument -- multi-point variants relabelled as single-point, which leaves the
    SP/MP guard unable to fire on the stores this module produces. `consolidate` always passes an
    explicit kind, so this only bites direct callers -- the latent case."""
    import tempfile
    from pathlib import Path

    from skempi_foldx import load_store, write_store
    from skempi_foldx.store import infer_kind

    mp = {"1XXX": {"LI38S,GI32Y": {t: 0.0 for t in TERMS}}}
    sp = {"1XXX": {"LI38S": {t: 0.0 for t in TERMS}}}
    assert infer_kind(mp) == "variants"
    assert infer_kind(sp) == "muts"

    with tempfile.TemporaryDirectory() as td:
        write_store(mp, Path(td))                      # no kind= on purpose
        payload = json.loads((Path(td) / "1XXX.json").read_text())
        assert "variants" in payload and "muts" not in payload, \
            "an MP store written without an explicit kind was relabelled single-point"
        assert load_store(Path(td)) == mp


def test_a_complete_repair_chain_is_not_redone(tmp_path):
    """Resume must stay cheap. The 4-round chains for 322 complexes are hours of RepairPDB; if a
    relaunch redid them the sweep could never be restarted."""
    from skempi_foldx import FoldxConfig, process_complex
    from skempi_foldx.skempi import SkempiComplex

    pdb = "1ZZZ"
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv",
                         binary=tmp_path / "foldx", repair_iterations=4)
    (tmp_path / "p").mkdir()
    (tmp_path / "p" / f"{pdb}.pdb").write_text(
        "ATOM      1  CA  LEU A   1       0.0   0.0   0.0  1.00  0.00           C\n")
    config.ensure_dirs()
    work = config.complex_work_dir(pdb)
    work.mkdir(parents=True, exist_ok=True)
    for i in range(1, 5):
        (work / f"repair_round_{i}.pdb").write_text("done\n")
    repaired = work / f"{pdb}_Repair.pdb"
    repaired.write_text("done\n")
    before = [(work / f"repair_round_{i}.pdb").stat().st_mtime_ns for i in range(1, 5)]

    # A binary that would fail loudly if invoked at all.
    (tmp_path / "foldx").write_text("#!/bin/sh\nexit 3\n")
    (tmp_path / "foldx").chmod(0o755)
    process_complex(pdb, [], SkempiComplex(pdb, "A", "B"), config, repair_only=True)

    after = [(work / f"repair_round_{i}.pdb").stat().st_mtime_ns for i in range(1, 5)]
    assert before == after, "a complete 4-round chain was redone on resume"
    assert repaired.read_text() == "done\n"


def test_two_campaigns_with_the_same_dirname_stay_distinct(tmp_path):
    """`_source` exists to say which campaign produced a number. Keying it on the directory
    BASENAME collides foldx_s1102/results with foldx_skempi_full/results -- same name, different
    campaigns, different mutation lists, therefore different energies."""
    from skempi_foldx import consolidate, source_label

    assert source_label("scratch/foldx_s1102/results") == "foldx_s1102_results"
    assert source_label("scratch/foldx_skempi_full/results") == "foldx_skempi_full_results"
    assert source_label("scratch/foldx_s1102/results") != source_label(
        "scratch/foldx_skempi_full/results")

    a, b = tmp_path / "campaign_a" / "results", tmp_path / "campaign_b" / "results"
    for d, val in ((a, 1.0), (b, 2.0)):
        d.mkdir(parents=True)
        (d / "1XXX.json").write_text(json.dumps(
            {"muts": {"LI38S": {t: val for t in TERMS}}, "meta": {"pdb": "1XXX"}}))

    _, report = consolidate([a, b], destination=None)
    assert len(report.per_source_contributed_mutations) == 2, \
        "two same-named campaign dirs collapsed into one report row"
    assert set(report.per_source_contributed_mutations) == {
        "campaign_a_results", "campaign_b_results"}


def test_default_config_stays_inside_the_package_tree():
    """`default_config()` infers the root from this file's location, so it is sensitive to how
    deep in a tree the package sits. Resolving one level too high puts the default
    results_dir/work_dir OUTSIDE the repo, and ensure_dirs() then creates them there -- from an
    installed wheel, next to site-packages."""
    import pathlib

    import skempi_foldx
    from skempi_foldx import default_config

    root = pathlib.Path(skempi_foldx.__file__).resolve().parents[1]
    cfg = default_config()
    for name, path in (("results_dir", cfg.results_dir), ("work_dir", cfg.work_dir),
                       ("pdb_dir", cfg.pdb_dir)):
        assert root in path.parents, f"{name}={path} escapes the package tree at {root}"


def test_bundled_results_are_reachable_without_a_clone():
    """The results are the deliverable, so they must resolve from the installed package rather
    than from a clone-relative 'data/results_sp', which only works when the cwd is a git clone."""
    from skempi_foldx import MULTI_POINT, bundled_path, load_bundled_store

    assert bundled_path().is_dir() and bundled_path(MULTI_POINT).is_dir()
    sp = load_bundled_store()
    assert len(sp) == 323 and sum(len(v) for v in sp.values()) == 4340
    assert sp["1BRS_A_D"]["DA52A"]["Interaction Energy"] == -0.6806

    # The multi-point definition count belongs beside the single-point one. It had been asserted
    # only by the workflow's wheel-check step, so trimming that step -- which looks like it merely
    # duplicates the counts above -- would have removed the sole check on it.
    mp = load_bundled_store(MULTI_POINT)
    assert len(mp) == 154 and sum(len(v) for v in mp.values()) == 1767


def test_load_store_raises_on_a_missing_directory(tmp_path):
    """Globbing a missing path and returning {} would surface a wrong path as a KeyError on the
    first lookup, several frames from the actual mistake."""
    import pytest as _pytest

    from skempi_foldx import load_store

    with _pytest.raises(FileNotFoundError):
        load_store(tmp_path / "nope")


def test_a_one_round_agreement_check_does_not_report_success():
    """A single round cannot disagree with itself, so comparing it against nothing yields an
    empty drift set — which reads exactly like a genuine pass. Agreement is only evidence when
    there are at least two rounds available to disagree.

    This matters because the verdict is used as a gate: the finalize chain treats a clean list
    check as licence to attribute a per-round energy difference to repair count alone."""
    import io
    from contextlib import redirect_stdout

    from list_hashes import verify

    # real manifest records carry both digests; verify() defaults to file_sha256
    rec = {"file_sha256": "x", "list_sha256": "x", "n_entries": 5}
    one = {"sp": {1: {"1ABC": dict(rec)}}, "mp": {}}
    two = {"sp": {1: {"1ABC": dict(rec)}, 2: {"1ABC": dict(rec)}}, "mp": {}}

    out = io.StringIO()
    with redirect_stdout(out):
        assert verify(one, [1, 2, 3, 4]) == 0
    text = out.getvalue()
    assert "identical across rounds" not in text, "one round reported as a successful comparison"
    assert "nothing to compare" in text

    out = io.StringIO()
    with redirect_stdout(out):
        assert verify(two, [1, 2, 3, 4]) == 0
    assert "identical across rounds" in out.getvalue(), "two agreeing rounds should pass"


# ----------------------------------------------------------------- the numeric core
# These pin the arithmetic and the index->mutation pairing. Without them, a sign flip on the
# subtraction, an off-by-one between a mutation and its energy file, a doubled parse, a wrong
# SKEMPI column, or swapped chain groups all leave the suite green -- which was the case.

def _ac_fxout(path, values):
    """An AnalyseComplex Interaction .fxout, in the shape parse_interaction_file expects."""
    path.write_text(
        "Interaction Residues Analysis\n\n"
        + "Pdb\t" + "\t".join(TERMS) + "\n"
        + "some.pdb\t" + "\t".join(f"{values[t]:.4f}" for t in TERMS) + "\n"
    )


def _ca(chain, pos, res3):
    return (f"ATOM  {pos:5d}  CA  {res3} {chain}{pos:4d}"
            f"       0.000   0.000   0.000  1.00  0.00           C\n")


def test_energies_are_mutant_minus_wildtype_and_paired_by_position(tmp_path):
    """End-to-end over the numeric core: two mutations with deliberately distinct energies,
    asserted BY MUTATION NAME.

    This pins three things at once. The subtraction direction: swapping it negates every term.
    The index->mutation pairing: BuildModel numbers its outputs by position in individual_list.txt,
    so an off-by-one hands a mutation another mutation's energies while still reporting ok(n/n).
    And the parse: a scaled or mis-keyed read changes the values.
    """
    from skempi_foldx import MODE_AUTHOR, FoldxConfig, process_complex
    from skempi_foldx.skempi import SkempiComplex

    pdb = "1TST"
    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 0\n")   # every FoldX output is pre-placed below
    binary.chmod(0o755)
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv", binary=binary)
    (tmp_path / "p").mkdir()
    structure = _ca("A", 38, "LEU") + _ca("B", 12, "GLY")
    (tmp_path / "p" / f"{pdb}.pdb").write_text(structure)
    config.ensure_dirs()
    work = config.complex_work_dir(pdb)
    work.mkdir(parents=True, exist_ok=True)
    (work / f"{pdb}.pdb").write_text(structure)
    (work / f"{pdb}_Repair.pdb").write_text(structure)      # skip RepairPDB

    # sorted() puts GB12A before LA38S, so entry 1 is GB12A and entry 2 is LA38S
    base = {t: 0.0 for t in TERMS}
    energies = {
        1: {"mut": dict(base, **{"Interaction Energy": 5.0, "Van der Waals": 1.5}),
            "wt":  dict(base, **{"Interaction Energy": 2.0, "Van der Waals": 1.0})},
        2: {"mut": dict(base, **{"Interaction Energy": -1.0, "Electrostatics": 0.25}),
            "wt":  dict(base, **{"Interaction Energy": 3.0, "Electrostatics": 0.75})},
    }
    for i, e in energies.items():
        (work / f"{pdb}_Repair_{i}.pdb").write_text(structure)
        (work / f"WT_{pdb}_Repair_{i}.pdb").write_text(structure)
        _ac_fxout(work / f"Interaction_{pdb}_Repair_{i}_AC.fxout", e["mut"])
        _ac_fxout(work / f"Interaction_WT_{pdb}_Repair_{i}_AC.fxout", e["wt"])

    entry = SkempiComplex(pdb, "A", "B", single={"LA38S", "GB12A"})
    result = process_complex(pdb, sorted(entry.single), entry, config, mode=MODE_AUTHOR)

    assert result.status == "ok(2/2)", result.status
    assert set(result.mutations) == {"GB12A", "LA38S"}

    # entry 1 == GB12A: 5.0 - 2.0 = +3.0, and 1.5 - 1.0 = +0.5
    assert result.mutations["GB12A"]["Interaction Energy"] == 3.0
    assert result.mutations["GB12A"]["Van der Waals"] == 0.5
    # entry 2 == LA38S: -1.0 - 3.0 = -4.0, and 0.25 - 0.75 = -0.5
    assert result.mutations["LA38S"]["Interaction Energy"] == -4.0
    assert result.mutations["LA38S"]["Electrostatics"] == -0.5
    # a sign flip would make these +/-; a pairing swap would exchange the two rows
    assert result.mutations["GB12A"]["Interaction Energy"] > 0
    assert result.mutations["LA38S"]["Interaction Energy"] < 0


def test_parse_interaction_file_reads_the_named_terms(tmp_path):
    """Directly: every term keyed by name, values not scaled or shifted."""
    from skempi_foldx.skempi import parse_interaction_file

    values = {t: round(0.5 + i, 4) for i, t in enumerate(TERMS)}
    path = tmp_path / "Interaction_X_AC.fxout"
    _ac_fxout(path, values)
    got = parse_interaction_file(path, TERMS)
    assert got == values, {k: (values[k], got.get(k)) for k in TERMS if got.get(k) != values[k]}


def test_load_skempi_reads_the_right_columns_and_splits_single_from_multi(tmp_path):
    """Chain groups feed --analyseComplexChains, so an orientation swap computes the wrong
    interface. The single/multi split decides which arm a mutation lands in."""
    from skempi_foldx import load_skempi

    csv = tmp_path / "skempi.csv"
    # The PDB-numbered and cleaned columns DIFFER, as they do in real SKEMPI. Identical values
    # would make reading the wrong column undetectable -- the pipeline consumes `cleaned`.
    csv.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        "1AAA_E_I;LI138S;LI38S;1e-9;1e-10\n"
        "1AAA_E_I;GI140A;GI40A;1e-9;1e-10\n"
        "1AAA_E_I;LI138S,GI140A;LI38S,GI40A;1e-9;1e-10\n"
        "2BBB_AB_CD;YC105F;YC5F;1e-9;1e-10\n"
    )
    sk = load_skempi(csv)
    # Keyed by SKEMPI identifier, not by PDB code: the chain pairing is half of what a value
    # means, so it belongs in the name of the record rather than in a field beside it.
    assert set(sk) == {"1AAA_E_I", "2BBB_AB_CD"}
    one, two = sk["1AAA_E_I"], sk["2BBB_AB_CD"]
    assert (one.group1, one.group2) == ("E", "I"), "chain groups swapped or misread"
    assert (two.group1, two.group2) == ("AB", "CD")
    assert one.groups == "E,I", "groups string feeds --analyseComplexChains"
    assert one.id == "1AAA_E_I" and one.pdb == "1AAA"
    assert one.single == {"LI38S", "GI40A"}
    assert one.multi == {"LI38S,GI40A"}, "multi-point row landed in the wrong arm"
    assert two.single == {"YC5F"} and not two.multi
    # the cleaned column, not the PDB-numbered one
    assert "LI138S" not in one.single, "read Mutation(s)_PDB instead of Mutation(s)_cleaned"


def test_term_vector_preserves_the_canonical_order():
    """Consumers flatten records into feature vectors; a reordering here permutes every feature
    downstream while every value remains individually correct."""
    from skempi_foldx.terms import term_vector

    rec = {t: float(i) for i, t in enumerate(TERMS)}
    assert term_vector(rec) == [float(i) for i in range(len(TERMS))]
    assert term_vector(rec)[0] == rec[SCALAR_TERM], "column 0 must be the scalar term"
    with pytest.raises(KeyError):
        term_vector({t: 0.0 for t in TERMS[:-1]})


def test_exclude_already_computed_actually_excludes(tmp_path):
    """It decides what a resumed campaign recomputes. Returning nothing would silently redo every
    complex -- expensive, and with the union-list rule it can change values on a partial rerun."""
    from skempi_foldx import exclude_already_computed

    done = tmp_path / "done"
    done.mkdir()
    (done / "1AAA.json").write_text(json.dumps(
        {"muts": {"LI38S": {t: 0.0 for t in TERMS}}, "meta": {"pdb": "1AAA"}}))

    worklist = {"1AAA": ["LI38S"], "2BBB": ["YC5F"]}
    remaining = exclude_already_computed(worklist, [done])
    assert "1AAA" not in remaining, "an already-computed complex was not excluded"
    assert remaining.get("2BBB") == ["YC5F"], "an uncomputed complex was dropped"


def test_skempi_reindex_is_lossless_and_the_default_is_untouched():
    """Every shipped key is now SKEMPI's own author form, so re-indexing is an identity map on
    this store. It stays because it is not an identity map on a *foreign* store: a campaign run
    in role mode still keys its output that way, and this is what re-indexes it.

    The default must stay byte-for-byte the stored keys either way -- rewriting a store to one
    convention was measured to drop a consumer's coverage from 3300/3300 to 1092/3300, silently,
    since the join simply stops matching and still produces well-formed output."""
    from skempi_foldx import load_bundled_store, reindex_by_skempi_id

    stored = load_bundled_store()
    skempi = load_bundled_store(key="skempi")

    assert sum(map(len, stored.values())) == sum(map(len, skempi.values())), "records lost"
    assert "LI38D" in stored["1ACB_E_I"], "the shipped key is SKEMPI's author form"
    assert stored["1ACB_E_I"] == skempi["1ACB_E_I"], "re-index is not an identity on this store"

    # a foreign store keyed the other way is what the re-index is for
    role_keyed = {"9ZZZ": {"LB38D": {"cleaned": "LI38D", "Interaction Energy": 1.0}}}
    assert "LI38D" in reindex_by_skempi_id(role_keyed["9ZZZ"])

    # no complex loses a record to a key collision
    for pdb, recs in stored.items():
        assert len(reindex_by_skempi_id(recs)) == len(recs), f"{pdb}: collision on re-index"

    with pytest.raises(ValueError):
        load_bundled_store(key="nonsense")


def test_a_directory_holding_both_arms_is_refused(tmp_path):
    """The two arms share a per-complex filename, so one directory cannot hold both without each
    complex meaning whichever arm ran last.

    store_kind used to return on the first non-empty file, so a mixed directory reported a single
    kind and consolidate's cross-source check then saw one consistent kind and passed -- the guard
    failing open on precisely the input it exists to reject. Reproduced before the fix: a mixed
    directory consolidated with no error and wrote a multi-point variant under the 'muts' key."""
    from skempi_foldx import consolidate, store_kind

    d = tmp_path / "mixed"
    d.mkdir()
    terms = {t: 0.0 for t in TERMS}
    (d / "1AAA.json").write_text(json.dumps({"muts": {"LI38S": terms}, "meta": {"pdb": "1AAA"}}))
    (d / "2BBB.json").write_text(json.dumps(
        {"variants": {"LI38S,GI40A": terms}, "meta": {"pdb": "2BBB"}}))

    with pytest.raises(ValueError, match="both"):
        store_kind(d)
    with pytest.raises(ValueError, match="both"):
        consolidate([d], destination=tmp_path / "out")
    assert not (tmp_path / "out").exists() or not list((tmp_path / "out").glob("*.json")), \
        "a mixed directory produced output"

    # a directory holding one arm is still fine
    single = tmp_path / "single"
    single.mkdir()
    (single / "1AAA.json").write_text(json.dumps({"muts": {"LI38S": terms}, "meta": {}}))
    assert store_kind(single) == "muts"
    # A missing directory has no kind. (load_store raises for one; store_kind reports None.)
    assert store_kind(tmp_path / "empty_missing") is None


def test_source_labels_are_the_documented_set():
    """`_source` is provenance and nothing joins on it, but it is read by humans, so the labels
    must say what they mean. The originals were directory paths -- `foldx_s1102_results_S4169`
    reads as two different benchmark subsets at once, because the S4169 campaign happened to run
    under a directory created for the S1102 work."""
    from skempi_foldx import MULTI_POINT, load_bundled_store

    expected = {"sweep_4x_round_1", "2C5D_A_C_sp", "2C5D_A_C_mp", "2C5D_AB_CD_mp",
                "3SE3_B_A_sp", "3SE3_B_A_mp", "3SE3_B_C_mp",
                "3SE4_B_A_sp", "3SE4_B_C_sp", "3SE4_B_C_mp"}
    seen = {r.get("_source") for which in ("results_sp", MULTI_POINT)
            for recs in load_bundled_store(which).values() for r in recs.values()}
    seen.discard(None)
    assert seen == expected, f"unexpected source labels: {seen ^ expected}"


# --------------------------------------------------------------------- naming resolution
# A mutation has more than one name, and every wrong answer this package has produced traced
# back to that. These pin the resolution rather than the plumbing around it.

def test_lookup_resolves_both_stored_conventions_without_extra_inputs():
    """Every shipped record is keyed by SKEMPI's author form and carries its role-chain name in
    a `role` field. A consumer holding either name must get the same record, with no
    configuration."""
    from skempi_foldx import FoldxLookup

    fx = FoldxLookup()
    author = fx.get("1ACB", "LI38D")      # SKEMPI's own form
    role = fx.get("1ACB", "LB38D")        # the form this record is filed under
    assert author is not None and role is not None
    assert author is role, "the two names resolved to different records"
    assert fx.arm_of("1ACB", "LI38D") == "sp"
    assert fx.arm_of("1A4Y", "KB40G,DA435A") == "mp", "multi-point variant did not resolve"


def test_lookup_reports_which_conventions_it_can_resolve():
    """What it can resolve depends on what it was given, and silently resolving less is how a
    coverage figure reads 69% for data that is 99% covered."""
    from skempi_foldx import FoldxLookup

    plain = FoldxLookup().conventions
    assert "stored key" in plain[0] and len(plain) == 3
    assert "shipped `role` field" in plain[2], \
        "the shipped role name is what resolves role-chain labels with no configuration"


def test_a_miss_names_the_other_form_and_the_fix():
    """A miss is usually a naming mismatch, not absent data. The error has to say so, or the
    caller concludes the coverage is worse than it is -- which is exactly what happened."""
    from skempi_foldx import FoldxLookup

    fx = FoldxLookup()
    with pytest.raises(KeyError, match="is not in the store"):
        fx.require("9ZZZ", "AA1G")
    with pytest.raises(KeyError) as caught:
        fx.require("1ACB", "NOPE1A")
    msg = str(caught.value)
    assert "It holds" in msg, "did not say what the complex does contain"
    # The bundled store resolves role-chain labels from its own `role` field, so there is no
    # input left to ask for. Advising skempi_csv= here would send a reader after a file that
    # cannot change the answer -- and one this package does not redistribute.
    assert "already resolves role-chain labels exactly" in msg, \
        "did not say the miss is an absent name rather than a form it could not match"
    assert "skempi_csv=" not in msg, "asked for an input that would change nothing"


def test_lookup_and_vector_agree_with_the_store():
    from skempi_foldx import TERMS, FoldxLookup, load_bundled_store

    fx = FoldxLookup()
    rec = load_bundled_store()["1BRS_A_D"]["DA52A"]
    assert fx.get("1BRS", "DA52A") == rec
    assert fx.vector("1BRS", "DA52A") == [float(rec[t]) for t in TERMS]
    assert fx.vector("9ZZZ", "AA1G") is None
    assert len(fx) == 6107, "both arms should be loaded by default"

    cov, total, missing = fx.coverage([("1BRS", "DA52A"), ("9ZZZ", "AA1G")])
    assert (cov, total, missing) == (1, 2, [("9ZZZ", "AA1G")])
    assert fx.energies_for([("1BRS", "DA52A"), ("9ZZZ", "AA1G")]) == {
        ("1BRS", "DA52A"): [float(rec[t]) for t in TERMS]}


def test_role_form_remap_is_exact_given_chain_mappings(tmp_path):
    """The offset remap is what turns a lower bound into an answer: identity matching cannot
    choose when both partners admit the same substitution, and offsets can."""
    from skempi_foldx import to_role_form

    chains = {"E": {"seq": ["x"] * 245}, "I": {"seq": ["x"] * 70}}
    # chain I is the whole of group 2, so its residues keep their numbering under side B
    assert to_role_form("LI38D", ("E", "I"), chains) == "LB38D"
    # chain E is group 1, likewise
    assert to_role_form("YE20A", ("E", "I"), chains) == "YA20A"
    # a second chain in a group is offset by the first chain's length
    two = {"A": {"seq": ["x"] * 10}, "B": {"seq": ["x"] * 10}, "C": {"seq": ["x"] * 5}}
    assert to_role_form("KB3G", ("AB", "C"), two) == "KA13G", "group offset not applied"
    # refuse rather than guess
    assert to_role_form("LI38D", ("E", "I"), None) is None
    assert to_role_form("LI38aD", ("E", "I"), chains) is None, "insertion code should not resolve"
    assert to_role_form("LZ38D", ("E", "I"), chains) is None, "chain outside both groups"
    assert to_role_form("LI900D", ("E", "I"), chains) is None, "residue beyond the chain"


# --------------------------------------------------------------- naming resolution, continued
# `to_role_form` and the aliasing built on it are what turn a 69% coverage figure into a 99% one,
# so the boundaries of the remap and the refusals around it are pinned by value here.

def _skempi_table(tmp_path, rows):
    """A minimal SKEMPI 2.0 table: `<pdb>_<group1>_<group2>;<pdb form>;<cleaned>;...`."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    path = tmp_path / "skempi.csv"
    path.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        + "".join(f"{ident};{cleaned};{cleaned};1e-9;1e-10\n" for ident, cleaned in rows))
    return path


def _mapping_dir(tmp_path, code, lengths):
    """A SKEMPI `<code>.mapping`: `RESNAME CHAIN AUTHOR_NUMBER SEQUENCE_INDEX`, one line per
    residue. Only the per-chain residue count is read, so the identities are filler."""
    directory = tmp_path / "mappings"
    directory.mkdir(exist_ok=True)
    lines = [f"ALA {chain} {i} {i}"
             for chain, n in lengths.items() for i in range(1, n + 1)]
    (directory / f"{code}.mapping").write_text("\n".join(lines) + "\n")
    return directory


def test_role_form_remap_includes_the_last_residue_of_a_chain():
    """The bound is inclusive at both ends. An off-by-one on the upper one returns None for every
    C-terminal mutation -- and None is also how this function refuses an insertion code, so the
    loss reads as "unresolvable label" rather than as a bug and coverage just quietly drops."""
    from skempi_foldx import to_role_form

    chains = {"E": {"seq": ["x"] * 245}, "I": {"seq": ["x"] * 70}}
    assert to_role_form("LI70D", ("E", "I"), chains) == "LB70D", "last residue must resolve"
    assert to_role_form("LI1D", ("E", "I"), chains) == "LB1D", "first residue must resolve"
    assert to_role_form("LI71D", ("E", "I"), chains) is None, "one past the end must not"
    assert to_role_form("LI0D", ("E", "I"), chains) is None, "position 0 must not"

    # the same at a group's internal boundary, where the offset is what is being pinned
    two = {"A": {"seq": ["x"] * 10}, "B": {"seq": ["x"] * 10}, "C": {"seq": ["x"] * 5}}
    assert to_role_form("KA10G", ("AB", "C"), two) == "KA10G", "first chain must not be offset"
    assert to_role_form("KB10G", ("AB", "C"), two) == "KA20G", "offset lost at the boundary"
    assert to_role_form("KB11G", ("AB", "C"), two) is None, "past the end of the second chain"


def test_chain_mapping_splits_a_fused_chain_and_author_number(tmp_path):
    """A four-digit author number runs into the chain letter, and dropping those lines is silent.

    The mapping files are fixed-width, so `ALA C1001  1` arrives as three whitespace fields. A
    parser requiring four drops every residue numbered 1000 or above: `2NYY` chain A then measures
    971 residues instead of 1267 and `to_role_form` rejects in-range mutations, while `1S0W` loses
    the whole of chain C and yields no role name at all. Neither raises, and for a LATER chain in
    a group a short length would shift every subsequent offset instead -- aliasing records onto
    residues they do not describe.
    """
    from skempi_foldx.lookup import load_chain_mapping, to_role_form

    directory = tmp_path / "mappings"
    directory.mkdir()
    (directory / "1ZZZ.mapping").write_text(
        "ALA A  26  1\n"
        "GLY A  27  2\n"
        "ALA C1001  1\n"        # fused: chain C, author 1001
        "GLY C1002  2\n"
        "VAL C1003  3\n")

    chains = load_chain_mapping(directory, "1ZZZ")
    assert set(chains) == {"A", "C"}, "the fused chain must not vanish"
    assert len(chains["C"]["seq"]) == 3
    assert chains["C"]["seq"] == ["1001", "1002", "1003"], "author numbers, separated"

    # With chain C present, a mutation on it resolves; without it the group loop bails to None.
    assert to_role_form("AC3F", ("A", "C"), chains) == "AB3F"


def test_shipped_store_matches_the_backfill_rule():
    """Every role-less shipped record is one whose role form IS its key.

    This is the gate on `experiments/backfill_role_names.py`: it re-derives, from the identifier
    alone, the same names that script writes, and fails if the store is missing any of them. The
    16 records the mapping-file bug once suppressed are covered by construction rather than by
    being listed, so the test does not go stale if another turns up.
    """
    import json

    from skempi_foldx.store import DATA_DIR, MULTI_KEY, SINGLE_KEY

    def role_of(mutation, groups):
        """A first chain has offset 0, so its role form needs no chain length."""
        if len(mutation) < 4:
            return None
        chain, (g1, g2) = mutation[1], groups
        if chain in g1:
            side, group = "A", g1
        elif chain in g2:
            side, group = "B", g2
        else:
            return None
        if group[0] != chain or not mutation[2:-1].isdigit():
            return None
        return f"{mutation[0]}{side}{mutation[2:-1]}{mutation[-1]}"

    missing = []
    for path in sorted(DATA_DIR.glob("results_*/*.json")):
        parts = path.stem.split("_")
        if len(parts) != 3:
            continue
        groups = (parts[1], parts[2])
        payload = json.loads(path.read_text())
        for key, record in (payload.get(SINGLE_KEY) or payload.get(MULTI_KEY) or {}).items():
            derived = [role_of(m, groups) for m in record.get("cleaned", key).split(",")]
            derived = ",".join(derived) if all(derived) else None
            if derived is None:
                continue                 # a later chain: this rule has no opinion either way
            stored = record.get("role")
            if stored is None:
                if derived != key:
                    missing.append((path.stem, key, "absent", derived))
            elif stored != derived:
                missing.append((path.stem, key, stored, derived))
    assert not missing, (
        "first-chain records whose stored `role` is absent or disagrees with the role form "
        f"derived from the identifier -- (file, key, stored, expected): {missing}")


def test_backfilled_roles_match_to_role_form():
    """The backfill's narrow rule and `to_role_form` agree wherever both apply.

    Two derivations of the same name exist -- one from chain mappings, one from the identifier
    alone -- and the store is written by the first while the audit script uses the second. Pinned
    here so they cannot drift: `to_role_form` is given a mapping long enough that its bound does
    not fire, which is the case where both are defined. `role_of` is imported from the script
    rather than restated, so that sabotaging the script's arithmetic fails this test.
    """
    from backfill_role_names import role_of
    from skempi_foldx.lookup import to_role_form

    chains = {"A": {"seq": ["x"] * 2000}, "C": {"seq": ["x"] * 2000},
              "D": {"seq": ["x"] * 2000}}
    for mutation, groups in (("AC142F", ("A", "C")), ("KA1028A", ("DC", "A")),
                             ("RA1266A", ("DC", "A")), ("YA1D", ("A", "C"))):
        assert to_role_form(mutation, groups, chains) == role_of(mutation, groups), mutation


def test_chain_mapping_reports_each_chain_length_and_nothing_at_all_for_a_missing_file(tmp_path):
    """The offsets a split file numbers by are accumulated from these lengths, so a chain read
    short shifts every later chain onto the wrong numbers. A missing file must read as "no
    mapping" rather than as an empty one, or the remap silently degrades to a guess."""
    from skempi_foldx import load_chain_mapping, to_role_form

    directory = _mapping_dir(tmp_path, "1CBW", {"F": 223, "G": 5, "H": 20, "I": 58})
    chains = load_chain_mapping(directory, "1CBW")
    assert {c: len(v["seq"]) for c, v in chains.items()} == {"F": 223, "G": 5, "H": 20, "I": 58}
    assert load_chain_mapping(directory, "9ZZZ") is None, "a missing mapping read as empty"

    # a group naming a chain the mapping does not carry cannot be offset, so it must refuse
    assert to_role_form("YF20A", ("FGH", "I"), {"F": {"seq": ["x"] * 223},
                                                "I": {"seq": ["x"] * 58}}) is None


def test_role_chain_labels_resolve_exactly_once_the_mappings_are_given(tmp_path):
    """1CBW is stored under its author form (`GI12A` -- chain I of `1CBW_FGH_I`). A split file
    calls that same mutation `GB12A`, and neither string matches the other; that mismatch is the
    whole of the ~29% a string-matching lookup misses."""
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A"), ("1CBW_FGH_I", "FI33A")])
    mappings = _mapping_dir(tmp_path, "1CBW", {"F": 223, "G": 5, "H": 20, "I": 58})

    fx = FoldxLookup(skempi_csv=csv, mapping_dir=mappings)
    assert fx.get("1CBW", "GB12A") is fx.get("1CBW", "GI12A") is not None, \
        "the role-chain label did not resolve to the stored record"
    assert fx.get("1CBW", "FB33A") is fx.get("1CBW", "FI33A") is not None
    assert fx.arm_of("1CBW", "GB12A") == "sp"
    # the side is read off the chain group, not assumed: chain I is group 2, so B and never A
    assert fx.get("1CBW", "GA12A") is None, "resolved a role form on the wrong partner"
    # and a mutation SKEMPI does not list for this complex still misses
    assert fx.get("1CBW", "GB99A") is None


def test_identity_matching_refuses_to_guess_where_both_partners_share_a_substitution(tmp_path):
    """Without the mapping files the role form is inferred from `(wt, position, mutant)` identity,
    which cannot choose when the same substitution exists on both partners. Guessing there is
    exactly how a record gets filed under another mutation's name; the miss must stand."""
    from skempi_foldx import FoldxLookup

    # G12A is present on chain I (group 2) AND on chain F (group 1) -- ambiguous.
    # F33A is on chain I only -- unambiguous.
    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A"), ("1CBW_FGH_I", "GF12A"),
                                   ("1CBW_FGH_I", "FI33A")])

    # A store with no shipped `role` field, so resolution has to come from identity matching.
    # The bundled store carries one, computed from the chain mappings at build time, and is
    # therefore exact by construction -- this invariant is about every other store.
    terms = {t: 1.0 for t in TERMS}
    foreign = tmp_path / "foreign_sp"
    foreign.mkdir()
    (foreign / "1CBW_FGH_I.json").write_text(json.dumps({"muts": {
        "GI12A": dict(terms, cleaned="GI12A"),
        "GF12A": dict(terms, cleaned="GF12A"),
        "FI33A": dict(terms, cleaned="FI33A")}, "meta": {}}))

    fx = FoldxLookup(stores=foreign, skempi_csv=csv)      # identity matching only
    assert fx.get("1CBW", "FB33A") is fx.get("1CBW", "FI33A") is not None, \
        "an unambiguous role form should still resolve without the mappings"
    assert fx.get("1CBW", "GB12A") is None, \
        "guessed a role form where both partners admit the substitution"

    # the mappings are what turn that refusal into an answer
    mappings = _mapping_dir(tmp_path, "1CBW", {"F": 223, "G": 5, "H": 20, "I": 58})
    exact = FoldxLookup(stores=foreign, skempi_csv=csv, mapping_dir=mappings)
    assert exact.get("1CBW", "GB12A") is exact.get("1CBW", "GI12A") is not None


def test_an_unresolved_role_form_registers_no_alias_at_all(tmp_path):
    """When neither the mappings nor identity matching can name a record's role form the answer is
    None, and None must not be stored as a name. An alias keyed on it makes `get(pdb, None)`
    return somebody's energies -- a lookup answering a question that was never asked."""
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A"), ("1CBW_FGH_I", "GF12A")])
    fx = FoldxLookup(skempi_csv=csv)                      # G12A is ambiguous, so role is None

    assert fx.get("1CBW", None) is None, "an unresolved role form was aliased under None"
    assert ("1CBW", None) not in fx
    assert fx.coverage([("1CBW", None)]) == (0, 1, [("1CBW", None)])


def test_conventions_report_which_resolution_is_actually_in_force(tmp_path):
    """The three states are not interchangeable and the difference is invisible in the numbers:
    identity matching leaves the ambiguous cases unresolved, so its coverage figure is a lower
    bound. Claiming the exact form without the mappings is how such a figure gets believed."""
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A")])
    mappings = _mapping_dir(tmp_path, "1CBW", {"F": 223, "G": 5, "H": 20, "I": 58})

    # A store with no shipped role names, so the three states are the ones actually in force.
    terms = {t: 1.0 for t in TERMS}
    foreign = tmp_path / "foreign_sp"
    foreign.mkdir()
    (foreign / "1CBW_FGH_I.json").write_text(json.dumps(
        {"muts": {"GI12A": dict(terms, cleaned="GI12A")}, "meta": {}}))

    plain = FoldxLookup(stores=foreign).conventions
    identity = FoldxLookup(stores=foreign, skempi_csv=csv).conventions
    exact = FoldxLookup(stores=foreign, skempi_csv=csv, mapping_dir=mappings).conventions

    assert len(plain) == 2 and not any("role-chain" in c for c in plain)
    assert len(identity) == 3 and "identity matching" in identity[2], identity
    assert len(exact) == 3 and "residue offsets" in exact[2], exact
    assert identity[2] != exact[2], "the two role-chain states read the same"
    assert "role-chain form (exact, via residue offsets)" in repr(
        FoldxLookup(stores=foreign, skempi_csv=csv, mapping_dir=mappings)), \
        "repr does not say what it resolves"

    # the shipped store needs neither argument, because the role name travels with the record
    assert "shipped `role` field" in FoldxLookup().conventions[2]


def test_the_arms_argument_restricts_what_is_loaded():
    """`arms=` is what a consumer of one arm passes to keep the other out of its coverage figures.
    Loading both regardless would report multi-point variants as covered by a single-point run."""
    from skempi_foldx import MULTI_POINT, SINGLE_POINT, FoldxLookup

    single = FoldxLookup(arms=(SINGLE_POINT,))
    multi = FoldxLookup(arms=(MULTI_POINT,))

    assert len(single) == 4340 and len(multi) == 1767
    assert len(FoldxLookup()) == len(single) + len(multi)

    assert single.get("1BRS", "DA52A") is not None
    assert single.get("1A4Y", "KB40G,DA435A") is None, "a multi-point record loaded under arms=sp"
    assert single.arm_of("1A4Y", "KB40G,DA435A") is None

    assert multi.get("1A4Y", "KB40G,DA435A") is not None
    assert multi.get("1BRS", "DA52A") is None, "a single-point record loaded under arms=mp"
    assert multi.arm_of("1A4Y", "KB40G,DA435A") == "mp"

    # and the miss says which arms were loaded, so the shortfall is attributable
    with pytest.raises(KeyError, match=MULTI_POINT):
        multi.require("9ZZZ", "AA1G")


def test_require_returns_the_record_and_otherwise_names_the_input_that_would_resolve_it(tmp_path):
    """Three different misses, three different fixes. A miss that reports "absent" when the record
    is filed under another name is what produced a 69% coverage figure for 99%-covered data, so
    each branch has to name the input that resolves the form it could not match."""
    from skempi_foldx import FoldxLookup, load_bundled_store

    plain = FoldxLookup()
    assert plain.require("1BRS", "DA52A") == load_bundled_store()["1BRS_A_D"]["DA52A"], \
        "require did not return the record on a hit"
    assert plain.require("1ACB", "LI38D") is plain.get("1ACB", "LB38D"), \
        "require and get disagreed on the same record under its two names"

    with pytest.raises(KeyError, match="is not in the store"):
        plain.require("9ZZZ", "AA1G")

    with pytest.raises(KeyError, match="already resolves role-chain labels exactly") as absent:
        plain.require("1BRS", "NOPE1A")
    assert "It holds" in str(absent.value), "did not say what the complex does contain"

    # A store with no role names is where the two constructor arguments still buy something.
    terms = {t: 1.0 for t in TERMS}
    foreign = tmp_path / "foreign_sp"
    foreign.mkdir()
    (foreign / "1BRS_A_D.json").write_text(json.dumps(
        {"muts": {"DA52A": dict(terms, cleaned="DA52A")}, "meta": {}}))
    (foreign / "1CBW_FGH_I.json").write_text(json.dumps(
        {"muts": {"GI12A": dict(terms, cleaned="GI12A")}, "meta": {}}))
    bare = FoldxLookup(stores=foreign)
    with pytest.raises(KeyError, match="skempi_csv=") as unwidened:
        bare.require("1BRS", "NOPE1A")
    assert "mapping_dir=" in str(unwidened.value), "did not name both missing inputs"

    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A")])
    identity = FoldxLookup(stores=foreign, skempi_csv=csv)
    with pytest.raises(KeyError, match="identity matching") as inexact:
        identity.require("1CBW", "NOPE1A")
    assert "mapping_dir=" in str(inexact.value), \
        "an identity-matching instance did not point at the exact remap"
    assert "skempi_csv=" not in str(inexact.value), "asked for an input it already has"


def test_vector_is_terms_ordered_floats_and_none_on_a_miss():
    """Consumers flatten a record positionally, so a vector in any other order is a silently
    permuted feature set. The fixture is a real record whose terms differ from one another --
    identical values would make a permutation undetectable."""
    from skempi_foldx import TERMS, FoldxLookup, load_bundled_store

    fx = FoldxLookup()
    rec = load_bundled_store()["1BRS_A_D"]["DA52A"]
    vector = fx.vector("1BRS", "DA52A")

    assert len(set(vector)) > 1, "a record with identical terms cannot detect a reordering"
    assert vector == [float(rec[t]) for t in TERMS]
    assert vector[0] == float(rec["Interaction Energy"]), "column 0 must be the scalar term"
    assert vector[TERMS.index("Van der Waals")] == float(rec["Van der Waals"])
    assert all(isinstance(v, float) for v in vector), "records may hold ints; the vector must not"
    assert fx.vector("1BRS", "NOPE1A") is None and fx.vector("9ZZZ", "AA1G") is None


def test_an_alias_is_kept_from_the_arm_that_registered_it_first(monkeypatch):
    """`cleaned` is not guaranteed unique across the two arms, and the alias table is filled in
    `arms` order. A later arm overwriting an existing alias would move an author-form name onto a
    different record without touching the stored key, so nothing downstream could notice."""
    from skempi_foldx import MULTI_POINT, SINGLE_POINT, FoldxLookup
    from skempi_foldx import lookup as lookup_module

    fake = {
        SINGLE_POINT: {"1XXX": {"LB38D": dict({t: 1.0 for t in TERMS}, cleaned="LI38D")}},
        MULTI_POINT: {"1XXX": {"LC38D": dict({t: 9.0 for t in TERMS}, cleaned="LI38D")}},
    }
    monkeypatch.setattr(lookup_module, "load_bundled_store", lambda which: fake[which])

    fx = FoldxLookup()
    assert fx.get("1XXX", "LI38D")["Interaction Energy"] == 1.0, \
        "a later arm overwrote an alias registered by an earlier one"
    assert fx.arm_of("1XXX", "LI38D") == "sp"
    assert fx.get("1XXX", "LC38D")["Interaction Energy"] == 9.0, "a stored key was shadowed"


# --------------------------------------------------------------------- the coverage reporter

def test_coverage_report_counts_what_resolves_and_names_what_does_not(tmp_path, monkeypatch,
                                                                     capsys):
    """The script exists to answer "how much FoldX signal does this evaluation set have", and it
    is read as a headline number. Two things must hold: the count is over rows it actually
    resolved, and a run without the mappings says its figure is a lower bound rather than
    presenting it as the answer."""
    import sys

    import coverage_report

    table = tmp_path / "rows.tsv"
    table.write_text("1BRS_A\tddg\tDA52A\n"          # stored key, single point
                     "1ACB.E.I_A\tddg\tLB38D\n"      # stored key, dotted id form
                     "1A4Y_A\tddg\tKB40G,DA435A\n"   # multi-point variant
                     "9ZZZ_A\tddg\tAA1G\n"           # absent
                     "9ZZZ_A\tddg\tAA1G\n")          # and the same row again: counted once
    out = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["coverage_report.py", "--table", str(table),
                                      "--label", "demo", "--json", str(out)])
    coverage_report.main()

    payload = json.loads(out.read_text())
    assert len(payload) == 1
    row = payload[0]
    assert row["label"] == "demo"
    assert row["rows"] == 4, "duplicate rows were not collapsed"
    assert (row["covered"], row["uncovered"]) == (3, 1)
    assert row["pct"] == 75.0
    assert (row["single_point"], row["multi_point"]) == (2, 1), "arms mis-attributed"
    assert row["worst_complexes"] == [{"pdb": "9ZZZ", "rows": 1}]
    assert row["exact"] is True, \
        "the shipped role names make role-chain resolution exact with no extra inputs"

    printed = capsys.readouterr().out
    assert "LOWER BOUND" not in printed, \
        "an exact figure was hedged as a lower bound; the shipped role names resolve these"
    assert "shipped `role` field" in printed, "did not say which resolution was in force"
    assert "9ZZZ(1)" in printed, "the uncovered complexes were not named"


def test_coverage_report_reads_a_split_directory_and_labels_each_input_in_order(tmp_path,
                                                                               monkeypatch):
    """--split and --table share one destination precisely so that labels stay paired with the
    inputs they name; collecting them separately mislabels every row once the two are
    interleaved."""
    import sys

    import coverage_report

    # the real layout: fold_<k>/<base>_{train,val,test}.tsv under one split directory
    split = tmp_path / "splits_kfold"
    (split / "fold_1").mkdir(parents=True)
    (split / "fold_2").mkdir(parents=True)
    (split / "fold_1" / "base_train.tsv").write_text("1BRS_A\tddg\tDA52A\n")
    (split / "fold_2" / "base_test.tsv").write_text("9ZZZ_A\tddg\tAA1G\n")
    table = tmp_path / "flat.tsv"
    table.write_text("1BRS_A\tddg\tDA52A\n")

    out = tmp_path / "report.json"
    monkeypatch.setattr(sys, "argv", ["coverage_report.py",
                                      "--table", str(table), "--split", str(split),
                                      "--label", "flat", "--label", "kfold",
                                      "--json", str(out)])
    coverage_report.main()

    flat, kfold = json.loads(out.read_text())
    assert (flat["label"], flat["rows"], flat["covered"]) == ("flat", 1, 1)
    assert (kfold["label"], kfold["rows"], kfold["covered"]) == ("kfold", 2, 1), \
        "labels were paired with the wrong inputs"
    assert kfold["files_read"] == 2, "the split directory was not read recursively"


# ------------------------------------------------------------------- the multi-point arm
# The single-point path is pinned end to end above; the variant path has its own validation, its
# own payload key and its own resume rule, and none of them is exercised by that test.

def test_a_variant_is_validated_substitution_by_substitution_and_rejected_whole(tmp_path):
    """MODE_VARIANT validates every substitution of a comma-joined variant against the repaired
    structure. Accepting a variant with one substitution whose wild-type residue is not there
    hands BuildModel a mutation of a residue that does not exist -- which it reports as a success,
    with a number attached."""
    from skempi_foldx import MODE_VARIANT, FoldxConfig, process_complex
    from skempi_foldx.skempi import SkempiComplex

    pdb = "1TST"
    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 0\n")     # every FoldX output is pre-placed below
    binary.chmod(0o755)
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv", binary=binary)
    (tmp_path / "p").mkdir()
    structure = _ca("A", 38, "LEU") + _ca("B", 12, "GLY")
    (tmp_path / "p" / f"{pdb}.pdb").write_text(structure)
    config.ensure_dirs()
    work = config.complex_work_dir(pdb)
    work.mkdir(parents=True, exist_ok=True)
    (work / f"{pdb}.pdb").write_text(structure)
    (work / f"{pdb}_Repair.pdb").write_text(structure)

    base = {t: 0.0 for t in TERMS}
    (work / f"{pdb}_Repair_1.pdb").write_text(structure)
    (work / f"WT_{pdb}_Repair_1.pdb").write_text(structure)
    _ac_fxout(work / f"Interaction_{pdb}_Repair_1_AC.fxout",
              dict(base, **{"Interaction Energy": 5.0}))
    _ac_fxout(work / f"Interaction_WT_{pdb}_Repair_1_AC.fxout",
              dict(base, **{"Interaction Energy": 2.0}))

    good, bad = "LA38S,GB12A", "LA38S,WB12A"     # structure has G at B12, not W
    entry = SkempiComplex(pdb, "A", "B", multi={good, bad})
    result = process_complex(pdb, [good, bad], entry, config, mode=MODE_VARIANT)

    assert result.status == "ok(1/2)", result.status
    assert set(result.mutations) == {good}, "a variant with an absent wild-type residue was built"
    assert result.mutations[good]["Interaction Energy"] == 3.0
    assert any("WB12A" in r for r in result.meta["validation_failed"]), \
        "the rejected substitution was not named"
    # the list handed to BuildModel holds the accepted variant only, comma-joined as FoldX wants
    assert (work / "individual_list.txt").read_text() == f"{good};\n"
    # and the payload keeps the multi-point key, so the SP/MP guard can still fire on it
    payload = json.loads(config.result_path(pdb).read_text())
    assert "variants" in payload and "muts" not in payload


def test_a_mutation_whose_wild_type_residue_is_absent_is_rejected_not_renamed(tmp_path):
    """The validation is what turns a numbering-convention mismatch into a loud rejection instead
    of a confidently wrong energy, so both halves matter: the match accepts, and the mismatch is
    reported with what the structure actually holds."""
    from skempi_foldx.skempi import repaired_wt_residues, validate_against_structure

    structure = tmp_path / "1TST_Repair.pdb"
    structure.write_text(_ca("A", 38, "LEU") + _ca("B", 12, "GLY"))
    residues = repaired_wt_residues(structure)
    assert residues == {("A", 38): "L", ("B", 12): "G"}

    validated, rejected = validate_against_structure(
        {"LB38S": "LA38S", "WB12A": "WB12A", "LA99S": "LA99S"}, residues)
    assert validated == [("LB38S", "LA38S")], "a validated pair lost its role-chain key"
    assert len(rejected) == 2
    assert "WB12A(structure has G)" in rejected[0], "the rejection must say what is really there"
    assert "structure has None" in rejected[1], "a position absent from the structure passed"


# ------------------------------------------------------------------------- worklist builders

def test_worklist_builders_take_the_union_per_complex_and_keep_the_arms_apart(tmp_path):
    """A mutation's energy depends on every entry preceding it in `individual_list.txt`, so a
    per-dataset SUBSET of a complex's mutations gives a different number for the same mutation --
    measured up to 11.03 kcal/mol apart. The union, sorted, is what makes a value reproducible;
    `worklist_from_table` is the one builder that does not do this, and it must at least keep
    multi-point rows out of the single-point arm."""
    from skempi_foldx import (load_skempi, worklist_from_table, worklist_multi_point,
                              worklist_single_point)

    csv = tmp_path / "skempi.csv"
    csv.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        "1AAA_E_I;LI138S;LI38S;1e-9;1e-10\n"
        "1AAA_E_I;GI140A;GI40A;1e-9;1e-10\n"
        "1AAA_E_I;LI138S,GI140A;LI38S,GI40A;1e-9;1e-10\n"
        "2BBB_AB_CD;YC105F;YC5F;1e-9;1e-10\n"
    )
    sk = load_skempi(csv)

    assert worklist_single_point(sk) == {"1AAA_E_I": ["GI40A", "LI38S"],
                                        "2BBB_AB_CD": ["YC5F"]}, \
        "the single-point worklist is not the sorted union of the definition's mutations"
    assert worklist_multi_point(sk) == {"1AAA_E_I": ["LI38S,GI40A"]}, "an arm leaked into the other"
    assert worklist_single_point(sk, only=["2BBB_AB_CD"]) == {"2BBB_AB_CD": ["YC5F"]}
    assert worklist_multi_point(sk, only=["2BBB_AB_CD"]) == {}

    table = tmp_path / "rows.tsv"
    table.write_text("1AAA_E\tddg\tLB38S\n1AAA_E\tddg\tGB40A\n1AAA_E\tddg\tLB38S,GB40A\n")
    assert worklist_from_table(table) == {"1AAA": ["GB40A", "LB38S"]}, \
        "a multi-point row reached the single-point worklist"


def test_store_coverage_counts_only_the_pairs_the_store_actually_holds(tmp_path):
    """`coverage` is the headline number for this channel and is exported for consumers to
    compute their own. Counting the wanted set instead of the found set reports 100% for a store
    that holds nothing."""
    from skempi_foldx import coverage

    store = {"1AAA": {"LI38S": {t: 0.0 for t in TERMS}}}
    covered, total, missing = coverage(store, [("1AAA", "LI38S"), ("1AAA", "GI40A"),
                                               ("9ZZZ", "AA1G")])
    assert (covered, total) == (1, 3)
    assert missing == [("1AAA", "GI40A"), ("9ZZZ", "AA1G")]
    assert coverage(store, []) == (0, 0, [])


# ------------------------------------------------------- interfaces, bypasses and honest labels
# Three guards whose failure mode is silence: a second interface definition quietly discarded, a
# stall guard quietly disarmed, and a lower bound quietly labelled exact.

def test_a_second_interface_definition_is_its_own_record(tmp_path):
    """SKEMPI records a few PDB codes under two interfaces. Keying by code pooled them under
    whichever row came first, so the other definition's mutations were scored against an
    interface they are not part of -- scoring all-zero for that reason, indistinguishable from a
    real 'no effect'. Keyed by identifier there is nothing to pool and nothing to announce."""
    from skempi_foldx import load_skempi
    from skempi_foldx.skempi import by_pdb

    csv = _skempi_table(tmp_path, [("3SE4_B_C", "EC69A"), ("3SE4_B_A", "DA117A")])
    with warnings.catch_warnings():
        warnings.simplefilter("error")       # nothing is discarded, so nothing is warned about
        sk = load_skempi(csv)

    assert set(sk) == {"3SE4_B_C", "3SE4_B_A"}, "a definition was dropped"
    assert sk["3SE4_B_C"].single == {"EC69A"}
    assert sk["3SE4_B_A"].single == {"DA117A"}, \
        "the second definition's mutation was pooled into the first"
    assert sk["3SE4_B_C"].groups == "B,C" and sk["3SE4_B_A"].groups == "B,A"

    grouped = by_pdb(sk)
    assert list(grouped) == ["3SE4"], "both definitions belong to one structure"
    assert [e.id for e in grouped["3SE4"]] == ["3SE4_B_A", "3SE4_B_C"], "not sorted by pairing"


def test_a_definition_whose_mutations_are_all_shared_is_still_its_own_record(tmp_path):
    """3SE3's second definition contributes no mutation the first lacks, so nothing was ever
    mis-scored for it -- and it still needs its own record, because its 2 variants hold values
    computed against the other pairing and only the identifier says so."""
    from skempi_foldx import load_skempi

    shared = _skempi_table(tmp_path / "shared", [("9CCC_B_A", "LA5D"), ("9CCC_B_C", "LA5D")])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        sk = load_skempi(shared)
    assert set(sk) == {"9CCC_B_A", "9CCC_B_C"}
    assert sk["9CCC_B_A"].single == sk["9CCC_B_C"].single == {"LA5D"}, \
        "a shared mutation belongs to both definitions, as two measurements"


@pytest.mark.parametrize("value,bypassed", [
    ("1", True), ("true", True), ("TRUE", True), ("yes", True), ("on", True), (" 1 ", True),
    ("", False), ("0", False), ("false", False), ("no", False), ("off", False), ("n", False),
])
def test_the_intractable_bypass_only_opens_on_an_affirmative_value(monkeypatch, value, bypassed):
    """The guard prevents a multi-day RepairPDB stall, so an unrecognised value must leave it
    armed. Under a negative test, ALLOW=no reads as a bypass and re-arms the very complex the
    module exists to keep out of the worklist."""
    from skempi_foldx.exclusions import ALLOW_ENV, allowed, is_excluded

    monkeypatch.setenv(ALLOW_ENV, value)
    assert allowed() is bypassed
    assert (is_excluded("1KBH") is None) is bypassed


def test_exactness_is_claimed_only_when_a_mapping_was_actually_read(tmp_path):
    """A misspelt or empty --mapping-dir resolves nothing and degrades to identity matching. A
    report that still called itself exact would be a lower bound wearing the wrong label."""
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("1CBW_FGH_I", "GI12A")])
    empty = tmp_path / "no_mappings"
    empty.mkdir()

    # A store carrying no shipped role names, so the mapping files are the only route to exact.
    terms = {t: 1.0 for t in TERMS}
    foreign = tmp_path / "foreign_sp"
    foreign.mkdir()
    (foreign / "1CBW_FGH_I.json").write_text(json.dumps(
        {"muts": {"GI12A": dict(terms, cleaned="GI12A")}, "meta": {}}))

    for mapping_dir in (empty, tmp_path / "does_not_exist"):
        fx = FoldxLookup(stores=foreign, skempi_csv=csv, mapping_dir=mapping_dir)
        assert fx.is_exact is False, f"claimed exactness with {mapping_dir.name}"
        assert not any("exact" in c for c in fx.conventions), \
            "conventions advertised a resolution that never ran"

    real = _mapping_dir(tmp_path, "1CBW", {"F": 223, "G": 5, "H": 20, "I": 58})
    fx = FoldxLookup(stores=foreign, skempi_csv=csv, mapping_dir=real)
    assert fx.is_exact is True
    assert any("exact" in c for c in fx.conventions)


def test_identity_matching_refuses_a_chain_that_is_in_neither_group(tmp_path):
    """A chain outside both groups has no role letter. Labelling it 'B' is a guess, in the one
    method whose contract is that it does not guess."""
    from skempi_foldx import load_skempi
    from skempi_foldx.lookup import FoldxLookup

    csv = _skempi_table(tmp_path, [("3SE4_B_C", "EC69A"), ("3SE4_B_A", "DA117A")])
    entry = load_skempi(csv)["3SE4_B_C"]
    fx = FoldxLookup.__new__(FoldxLookup)
    assert fx._role_by_identity(entry, "DA117A") is None, \
        "invented a role letter for a chain outside both groups"
    assert fx._role_by_identity(entry, "EC69A") == "EB69A"


# ------------------------------------------------- indexing a store that is not the shipped one
# The naming problem is identical for a campaign's own output, so a lookup that can only index the
# bundled data leaves its owner joining dictionaries by hand -- which is the failure this class
# exists to prevent.

def test_a_lookup_can_index_a_store_that_did_not_ship_with_the_package(tmp_path):
    terms = {t: 1.0 for t in TERMS}
    own = tmp_path / "results_sp"
    own.mkdir()
    (own / "1XYZ.json").write_text(json.dumps(
        {"muts": {"LB38D": dict(terms, cleaned="LI38D", _source="mine")}, "meta": {}}))

    from skempi_foldx import FoldxLookup

    fx = FoldxLookup(stores=own)
    assert len(fx) == 1, "did not index the supplied store"
    assert fx.get("1XYZ", "LB38D") is not None, "stored key did not resolve"
    assert fx.get("1XYZ", "LI38D") is fx.get("1XYZ", "LB38D"), \
        "the author form did not alias onto the record, which is the whole point"
    assert fx.get("1BRS", "DA52A") is None, "the bundled store leaked in alongside it"


def test_the_arm_of_a_supplied_store_is_read_from_it_rather_than_declared(tmp_path):
    """A caller who has to declare the arm can declare it wrongly. It is already recorded -- on
    disk by the key each file uses."""
    terms = {t: 1.0 for t in TERMS}
    mp = tmp_path / "mine_mp"
    mp.mkdir()
    (mp / "1XYZ.json").write_text(json.dumps(
        {"variants": {"LI38D,GI32Y": dict(terms, _source="mine")}, "meta": {}}))

    from skempi_foldx import FoldxLookup

    assert FoldxLookup(stores=mp).arm_of("1XYZ", "LI38D,GI32Y") == "mp"
    assert FoldxLookup(stores={"1XYZ": {"LI38D,GI32Y": terms}}).arm_of("1XYZ", "LI38D,GI32Y") == "mp"
    assert FoldxLookup(stores={"1XYZ": {"LI38D": terms}}).arm_of("1XYZ", "LI38D") == "sp"


def test_a_doubled_code_resolves_only_when_the_caller_names_the_pairing():
    """The registry of mis-paired records is gone because the records are no longer mis-paired.
    What replaces it is an ambiguity a caller cannot walk past: asking for a code that SKEMPI
    defines twice, for a mutation both definitions carry, has two right answers and no default.

    `skempi_v2.csv` is not redistributed, so this has to hold without it."""
    from skempi_foldx import FoldxLookup

    fx = FoldxLookup()
    doubled = [c for c in ("2C5D", "3SE3", "3SE4") if len(fx.definitions_of(c)) > 1]
    assert doubled, "no code in the shipped store carries two interface definitions"

    for code in doubled:
        idents = fx.definitions_of(code)
        assert all(i.startswith(code + "_") for i in idents)
        # Each identifier answers on its own.
        for ident in idents:
            held = [k for (p, k) in fx._records if p == ident]
            assert held, f"{ident} holds no records"
            assert fx.get(ident, held[0]) is not None, "an identifier failed to resolve"

    assert len(fx.definitions_of("1BRS")) == 1, "a singly-defined code gained a second definition"

    # The ambiguity itself, which is the part the docstring promises and the rest of this test
    # only sets up. Without these assertions `_resolve` could return the first candidate instead
    # of refusing, and every check above would still pass -- the store would answer a question
    # about an interface the caller never named, which is the 0.1.0 defect wearing a new shape.
    ambiguous = [
        (code, name)
        for code in doubled
        for name in {k for (ident, k) in fx._records if ident in fx.definitions_of(code)}
        if len(fx._candidates(code, name)) > 1
    ]
    assert len(ambiguous) == 8, "the shipped store no longer holds the 8 documented collisions"
    for code, name in ambiguous:
        assert fx.get(code, name) is None, f"{code} {name} resolved to one of two definitions"
        with pytest.raises(KeyError, match="interface definitions"):
            fx.require(code, name)


def test_identity_inference_refuses_a_chain_whose_offset_it_cannot_know(tmp_path):
    """`_role_by_identity` is the fallback used when no `mapping_dir` is supplied, and it must not
    invent a label.

    The role convention renumbers each chain by the cumulative residue count of the chains before
    it in its group, and this path has no residue counts. The author position is therefore the
    role position only for the first chain of a group. Emitting it for any other chain produces a
    label in neither convention — the role chain letter with the author number — which names no
    record, so it installs an alias onto whichever record carries that substitution and answers
    confidently. That is the silent-wrong-answer failure this package exists to prevent, so the
    fallback refuses instead."""
    from skempi_foldx import FoldxLookup, SkempiComplex

    fx = FoldxLookup()
    # A three-chain group: a mutation on B carries an unknown offset, one on A does not.
    entry = SkempiComplex("9XYZ", "ABC", "DE", single={"LA10D", "LB20D"})
    assert fx._role_by_identity(entry, "LA10D") == "LA10D", \
        "refused a first-in-group chain, whose offset is zero"
    assert fx._role_by_identity(entry, "LB20D") is None, \
        "invented a role label for a chain whose offset it cannot know"

    second = SkempiComplex("9XYZ", "A", "DE", single={"LE10D"})
    assert second.group2[0] == "D"
    assert fx._role_by_identity(second, "LE10D") is None, \
        "invented a role label for a non-first chain of the second group"


def test_a_role_alias_is_never_overwritten_by_a_later_claimant():
    """The documented resolution is that the first claimant keeps a contested alias, matching
    `consolidate()`. Only the warning was asserted before; this pins the answer itself, which is
    what a consumer actually gets."""
    import warnings as _warnings

    from skempi_foldx import FoldxLookup
    from skempi_foldx.terms import TERMS

    terms = {t: 1.0 for t in TERMS}
    store = {"1ABC_A_B": {"K1": dict(terms, role="AB9C"),
                          "K2": dict(terms, role="AB9C", **{TERMS[0]: -999.0})}}
    with _warnings.catch_warnings(record=True):
        _warnings.simplefilter("ignore")
        fx = FoldxLookup(stores=[store])
    got = fx.get("1ABC_A_B", "AB9C")
    assert got is not None and got[TERMS[0]] == 1.0, \
        "a later record took a contested alias from the one that claimed it first"
    assert fx.conflicts, "the contested alias was resolved without being reported"


def test_reindexing_announces_a_name_two_records_both_claim():
    """The function exists for stores that mix naming conventions, which is exactly where two
    records can re-key onto one name. Dropping the loser by dict order and saying nothing loses a
    record from a store the caller believes was only re-indexed."""
    import warnings as _warnings

    from skempi_foldx import load_bundled_store, reindex_by_skempi_id

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        out = reindex_by_skempi_id({"K1": {"cleaned": "SHARED", "Interaction Energy": 1.0},
                                    "K2": {"cleaned": "SHARED", "Interaction Energy": -9.0}})
    assert len(out) == 1 and out["SHARED"]["Interaction Energy"] == 1.0, "first no longer wins"
    assert [w for w in caught if "re-key onto a name" in str(w.message)], \
        "a record was dropped without saying so"

    # And it stays an identity map, without warning, on the shipped store.
    sp = load_bundled_store()
    with _warnings.catch_warnings(record=True) as clean:
        _warnings.simplefilter("always")
        for ident, records in sp.items():
            assert len(reindex_by_skempi_id(records)) == len(records), ident
    assert not [w for w in clean if "re-key onto a name" in str(w.message)]


def test_a_mapping_file_header_row_does_not_become_a_chain(tmp_path):
    """Mapping files are caller-supplied and are how residue offsets are computed. Any line with
    four fields used to become a residue, so a header row parsed as a chain of its own and every
    offset after it shifted — which is precisely how a mapping set comes to describe different
    numbering than the one that built the store."""
    from skempi_foldx import load_chain_mapping

    maps = tmp_path / "maps"
    maps.mkdir()
    (maps / "9XYZ.mapping").write_text(
        "RESNAME CHAIN AUTHOR SEQ\n"          # a header row: four fields, none of them residues
        "ALA A 1 1\n"
        "ALA A 2 2\n"
        "ALA B 1 3\n")
    chains = load_chain_mapping(maps, "9XYZ")
    assert set(chains) == {"A", "B"}, f"parsed a phantom chain: {sorted(chains)}"
    assert len(chains["A"]["seq"]) == 2 and len(chains["B"]["seq"]) == 1


def test_the_exclusion_registry_answers_an_identifier_not_only_a_code():
    """`INTRACTABLE` is keyed by PDB code, and 0.2.0 keys stores by interface identifier.

    `exclusions`' own module docstring names this as the failure to avoid: a registry keyed by
    code and queried with an identifier answers None for every entry it holds, which reads exactly
    like "nothing is excluded". `consolidate` calls `filter_complexes(list(store))` with
    identifiers, so the identifier form is the one that matters, and every existing test passed
    the bare code."""
    from skempi_foldx import filter_complexes, is_excluded

    for name in ("1KBH", "1KBH_A_B", "1KBH_AB_CD"):
        assert is_excluded(name), f"{name} is not recognised as excluded"
    assert not is_excluded("1BRS") and not is_excluded("1BRS_A_D")

    kept = filter_complexes(["1BRS_A_D", "1KBH_A_B", "1KBH", "1ACB_E_I"])
    assert "1KBH_A_B" not in kept and "1KBH" not in kept, \
        "an excluded complex survived filtering under its identifier"
    assert kept == ["1BRS_A_D", "1ACB_E_I"]


def test_a_mapping_dir_that_contradicts_the_store_is_refused_not_believed(tmp_path):
    """The same fabrication as `_role_by_identity`, reached through the exact path.

    A record that ships a `role` already carries this computation's known-good answer. If
    `mapping_dir=` recomputes a different name, the mapping files describe different residue
    numbering than the ones that built the store, and installing the recomputed name would alias
    a label belonging to another mutation — answered confidently and counted as covered. The
    shipped name wins and the disagreement is announced."""
    import warnings as _warnings

    from skempi_foldx import FoldxLookup

    maps = tmp_path / "maps"
    maps.mkdir()
    # 1AO7's chain A deliberately short, so every later chain's offset shifts.
    lines = ([f"ALA A {i} {i}" for i in range(1, 160)]
             + [f"ALA B {i} {i}" for i in range(1, 200)]
             + [f"ALA C {i} {i}" for i in range(1, 20)]
             + [f"ALA D {i} {i}" for i in range(1, 200)]
             + [f"ALA E {i} {i}" for i in range(1, 200)])
    (maps / "1AO7.mapping").write_text("\n".join(lines) + "\n")

    table = _skempi_table(tmp_path, [("1AO7_ABC_DE", "AE50P")])
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        fx = FoldxLookup(skempi_csv=table, mapping_dir=maps)

    assert fx.role_mismatches, "a mapping that contradicts the shipped role went unnoticed"
    assert [w for w in caught if "disagrees with what mapping_dir" in str(w.message)]

    # Nothing fabricated: every alias still points at a record whose own role name it is.
    for (pdb, name), target in fx._alias.items():
        shipped = fx._records[target].get("role")
        assert shipped in (None, name) or name == fx._records[target].get("cleaned"), \
            f"{pdb} {name} aliases a record whose role is {shipped}"

    # And the guard must not fire when nothing contradicts it.
    with _warnings.catch_warnings(record=True):
        _warnings.simplefilter("ignore")
        plain = FoldxLookup()
    assert plain.role_mismatches == (), "the shipped store contradicts itself"


def test_the_corrections_manifest_ships_and_carries_both_releases_keys():
    """The manifest exists so a 0.1.0 consumer can answer "is one of my rows provably wrong?"
    without a git checkout and without re-deriving the join.

    The join is the point. Eleven of the 26 are keyed differently across the releases — 0.1.0
    keys them under the bare code with the other pairing's chain letter — so a manifest carrying
    only the 0.2.0 name would find nothing for exactly the records with the largest corrections
    and read as "unaffected". Both keys ship, and every 0.1.0 key must resolve in a 0.1.0 store.
    """
    import csv
    import importlib.resources as res

    with res.files("skempi_foldx").joinpath("data/corrections_from_0_1_0.csv").open() as fh:
        rows = list(csv.DictReader(fh))

    assert len(rows) == 31, "the manifest no longer holds the whole 0.1.0 defect registry"
    moved = [r for r in rows if r["moved"] == "True"]
    assert len(moved) == 26
    assert sum(1 for r in moved if r["clean_zero_in_010"] == "True") == 15
    assert sum(1 for r in moved if r["interaction_energy_moved"] == "False") == 1

    differing = [r for r in moved if r["mutation_010"] != r["mutation_020"]]
    assert len(differing) == 11, "the cross-release key difference changed shape"
    assert all(r["identifier_020"].startswith("3SE4") for r in differing)

    # Every 0.2.0 side must name a record this release actually ships.
    from skempi_foldx import FoldxLookup
    fx = FoldxLookup()
    for r in moved:
        assert fx.get(r["identifier_020"], r["mutation_020"]) is not None, \
            f"{r['identifier_020']} {r['mutation_020']} is not in the shipped store"


def test_the_shipped_store_announces_its_malformed_rows_on_construction():
    """The third construction warning, and the only one the shipped store always raises.

    Two SKEMPI rows are malformed upstream and are computed and shipped rather than dropped, so
    every `FoldxLookup()` indexes them. A consumer who never sees this has values whose arity
    label is wrong and no reason to suspect it -- which is why the condition is a warning at
    construction rather than a note in a document."""
    import warnings as _warnings

    from skempi_foldx import MALFORMED_ROWS, FoldxLookup

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        fx = FoldxLookup()
    said = [w for w in caught
            if "malformed" in str(w.message) and w.category is RuntimeWarning]
    assert said, "the shipped store indexed malformed rows without saying so"
    assert "3-substitution mutants labelled as" in str(said[0].message), \
        "the warning no longer says what is wrong with the rows"
    assert len(fx.malformed) == len(MALFORMED_ROWS) == 2
    for pdb, mutation in fx.malformed:
        assert fx.is_malformed(pdb, mutation), "a row it warned about does not test as malformed"


def test_every_name_that_becomes_a_path_is_validated(tmp_path):
    """One test per join, because each is a separate door onto the same hazard.

    Identifiers and PDB codes both come from column 0 of a caller-supplied `skempi_v2.csv`, which
    nothing upstream constrains. `write_store` was closed first; these are the rest of the places
    a name reaches the filesystem, and they were closed later and separately."""
    from skempi_foldx import FoldxConfig, load_chain_mapping
    from skempi_foldx.store import check_code, check_identifier

    escape = "../" + "ESCAPED"
    cfg = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                      pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv")
    for bad in (escape, "..", "a/b", ""):
        with pytest.raises(ValueError, match="not a SKEMPI identifier"):
            cfg.result_path(bad)
        with pytest.raises(ValueError, match="not a SKEMPI identifier"):
            cfg.complex_work_dir(bad)
        with pytest.raises(ValueError, match="not a SKEMPI identifier"):
            cfg.find_repaired(bad)

    maps = tmp_path / "maps"
    maps.mkdir()
    (tmp_path / "secret.mapping").write_text("ALA A 1 1\n")
    with pytest.raises(ValueError, match="not a PDB code"):
        load_chain_mapping(maps, escape + "/secret")
    assert load_chain_mapping(maps, "1ABC") is None, "a valid code stopped resolving"

    # The anchor matters: without it the pattern matches a prefix and a trailing newline or path
    # separator rides along into the filename.
    for bad in ("1ABC\n", "1ABCD", "1AB", "1ABC/x", "../1"):
        with pytest.raises(ValueError, match="not a PDB code"):
            check_code(bad)
    assert check_code("2C5D") == "2C5D"
    # A bare code must still pass the identifier pattern: repairs are sought by code.
    assert check_identifier("2C5D") == "2C5D"


def test_a_hostile_pdb_code_never_reaches_the_input_directory(tmp_path):
    """`entry.pdb` is the other name that becomes a path -- the input structure under `pdb_dir`
    -- and it is parsed from the same unconstrained column as the identifier. `load_skempi`
    splits on `_` so it cannot produce one of these, but any caller constructing a
    `SkempiComplex` directly can, and `experiments/recompute_alternate_interfaces.py` does."""
    from skempi_foldx import FoldxConfig, SkempiComplex, process_complex

    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)
    (tmp_path / "secret.pdb").write_text("ATOM\n")
    cfg = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                      pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv", binary=binary)
    cfg.ensure_dirs()
    (tmp_path / "p").mkdir(exist_ok=True)

    hostile = SkempiComplex("../" + "secret", "A", "B")
    with pytest.raises(ValueError, match="not a PDB code"):
        process_complex("1ABC_A_B", ["LA1A"], hostile, cfg)


def test_a_relative_binary_name_is_resolved_before_it_is_run(tmp_path):
    """`_foldx` runs with `cwd=` set to a work directory, and `execvp` searches PATH rather than
    the current directory. A bare name that `exists()` relative to the caller's cwd would pass the
    check here and then run whatever PATH finds -- so the check and the execution would be about
    two different files."""
    import os

    from skempi_foldx import FoldxConfig

    binary = tmp_path / "foldx"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(0o755)

    cwd = os.getcwd()
    try:
        os.chdir(tmp_path)
        resolved = FoldxConfig(results_dir=tmp_path, work_dir=tmp_path,
                               pdb_dir=tmp_path, skempi_csv=tmp_path / "s.csv",
                               binary="foldx").resolve_binary()
    finally:
        os.chdir(cwd)
    assert resolved.is_absolute(), "a relative binary name survived resolution"
    assert resolved == binary.resolve()


def test_the_store_builder_validates_before_it_deletes(tmp_path):
    """`write()` removes the previous store before writing the new one, and that removal cannot
    be walked back. Validating partway through would cost the data the validation exists to
    protect, so every name is checked before anything is removed."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "_bpds", pathlib.Path(__file__).resolve().parents[1]
        / "experiments" / "build_per_definition_store.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    dest = tmp_path / "store"
    (dest / "results_sp").mkdir(parents=True)
    keep = dest / "results_sp" / "KEEP.json"
    keep.write_text("{}")

    terms = {t: 1.0 for t in TERMS}
    poisoned = {("1ABC_A_B", "sp"): ({"LA1A": dict(terms)}, {}),
                ("../" + "ESCAPED", "sp"): ({"LA1A": dict(terms)}, {})}
    with pytest.raises(ValueError, match="not a SKEMPI identifier"):
        module.write(poisoned, dest)
    assert keep.exists(), "the previous store was deleted before the new names were checked"


def test_a_code_keyed_store_alongside_this_one_is_announced(tmp_path):
    """The one route by which a bare doubled code still answers.

    `_resolve` returns an exact key hit before consulting `_candidates`, so a pre-0.2.0 store
    keyed by PDB code puts a pooled record under `2C5D` where the refusal cannot reach it, and
    `("2C5D", m)` versus `("2C5D_A_C", m)` are different keys so no alias conflict is recorded
    either. Which record the caller wanted depends on what they meant, so this reports rather
    than decides -- but it must not be silent."""
    import warnings as _warnings
    from skempi_foldx import MULTI_POINT, SINGLE_POINT, FoldxLookup, bundled_path

    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "2C5D.json").write_text(json.dumps(
        {"muts": {"DC61G": {t: -999.0 for t in TERMS}}, "meta": {}}, indent=1))

    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        mixed = FoldxLookup(stores=[str(legacy), str(bundled_path(SINGLE_POINT)),
                                    str(bundled_path(MULTI_POINT))])
    assert mixed.pooled_codes == ("2C5D",)
    assert [w for w in caught if "indexed both as a record key" in str(w.message)], \
        "a code-keyed store shadowed a definition-keyed one silently"

    with _warnings.catch_warnings(record=True) as clean_run:
        _warnings.simplefilter("always")
        shipped = FoldxLookup()
    assert shipped.pooled_codes == (), "the shipped store pools a code under its own name"
    assert not [w for w in clean_run if "indexed both as a record key" in str(w.message)], \
        "the warning fires on a store that mixes no conventions"


def test_a_non_finite_term_is_refused_by_both_flatteners_and_by_the_audit():
    """Three doors onto the same hazard, and they must agree.

    `json.loads` accepts bare `NaN` and `Infinity`, so a foreign results directory can carry one.
    `FoldxLookup.vector` and `term_vector` are both documented as the way to flatten a record, so
    a value one refuses and the other returns would mean the choice of flattener decides whether
    a feature matrix is poisoned. `_max_abs_delta` is the third: every comparison against `nan`
    is False, so returning it would make `audit()` report two sources that disagree as agreeing --
    failing open on exactly the values the other two refuse.
    """
    import math

    from skempi_foldx import term_vector
    from skempi_foldx.store import _max_abs_delta

    finite = {t: 1.0 for t in TERMS}
    for bad in (float("nan"), float("inf"), float("-inf")):
        poisoned = dict(finite, **{TERMS[0]: bad})
        with pytest.raises(ValueError, match="non-finite"):
            term_vector(poisoned)
        assert _max_abs_delta(finite, poisoned) == math.inf, \
            "a non-finite term compared as agreement"
        assert _max_abs_delta(finite, poisoned) > 1e-6, "the audit would not flag it"

    assert term_vector(finite) == [1.0] * len(TERMS)
    assert _max_abs_delta(finite, dict(finite, **{TERMS[0]: 2.0})) == 1.0

    # The other way the same guard fails open: a comparison that could not be made must not read
    # as one that succeeded. Two records sharing no term have not been shown to agree.
    assert _max_abs_delta({}, {}) == math.inf, "two empty records compared as identical"
    assert _max_abs_delta({TERMS[0]: 1.0}, {TERMS[1]: 99.0}) == math.inf, \
        "records with no term in common compared as identical"


def test_an_identifier_cannot_become_a_path(tmp_path):
    """Store keys come from column 0 of a caller-supplied `skempi_v2.csv`, which nothing upstream
    constrains, and `write_store` turns each of them into a filename. Validated where the name
    becomes a path, since that is the one place every route passes through.

    The read side is already closed -- `bundled_path` whitelists its argument -- so this is the
    write side, which a consumer running `consolidate()` over their own table reaches."""
    from skempi_foldx import write_store
    from skempi_foldx.store import check_identifier

    terms = {t: 1.0 for t in TERMS}
    outside = tmp_path / "ESCAPED.json"
    for bad in ("../ESCAPED", "/" + "tmp/ESCAPED", "..", "a/b", "", "1ABC_A_B/../../x"):
        with pytest.raises(ValueError, match="not a SKEMPI identifier"):
            write_store({bad: {"LA1A": dict(terms)}}, tmp_path / "out", kind="muts")
    assert not outside.exists(), "a crafted identifier wrote outside the destination"

    write_store({"1ABC_A_B": {"LA1A": dict(terms)}}, tmp_path / "out", kind="muts")
    assert (tmp_path / "out" / "1ABC_A_B.json").exists(), "a valid identifier stopped working"
    assert check_identifier("2C5D") == "2C5D"


def test_foldx_is_invoked_without_a_shell():
    """`pdb` and `entry.groups` reach the FoldX argument list from the SKEMPI table, so a shell
    would make that column executable. Nothing in the invocation needs one."""
    import inspect

    from skempi_foldx import run as run_module

    source = inspect.getsource(run_module)
    assert "shell=True" not in source, "a shell reappeared in the compute path"
    assert "shell=False" in source, "the no-shell choice is no longer explicit"
    argv = inspect.signature(run_module._foldx).parameters
    assert "argv" in argv, "_foldx no longer takes an argument vector"


def test_the_conflict_guard_is_queryable_and_does_not_fire_on_one_record():
    """Two facts about the guard a consumer cannot get from a warning.

    First, it has to be checkable. docs/USAGE.md recommends narrowing these warnings to `once`,
    and a consumer that does has no way to ask afterwards which names collided -- so the answer
    is an attribute, not only a message.

    Second, it must not cry wolf. A record whose `role` and `cleaned` are the same string, both
    differing from its key, registers one alias and then meets its own entry; testing for
    presence rather than for a different target reported that as two records claiming one name."""
    import warnings as _warnings
    from skempi_foldx import FoldxLookup
    from skempi_foldx.terms import TERMS

    terms = {t: 1.0 for t in TERMS}

    one = {"1ABC_A_B": {"KEYFORM": dict(terms, role="RB1A", cleaned="RB1A")}}
    with _warnings.catch_warnings(record=True) as caught:
        _warnings.simplefilter("always")
        fx = FoldxLookup(stores=[one])
    assert fx.conflicts == (), "one record with two spellings of one alias reported a collision"
    assert not [w for w in caught if "resolve to more than one" in str(w.message)]
    assert fx.get("1ABC_A_B", "RB1A") is not None and fx.get("1ABC_A_B", "KEYFORM") is not None

    two = {"1ABC_A_B": {"K1": dict(terms, cleaned="SHARED"),
                        "K2": dict(terms, cleaned="SHARED")}}
    with _warnings.catch_warnings(record=True):
        _warnings.simplefilter("ignore")
        real = FoldxLookup(stores=[two])
    assert real.conflicts == (("1ABC_A_B", "SHARED"),), \
        "two records claiming one author form went unrecorded"


def test_a_caller_supplied_table_refuses_a_doubled_code_too(tmp_path):
    """`_sole_definition` is the same refusal as `_resolve`, on the other input.

    It exists so a store written under bare codes -- every store this package shipped before
    0.2.0 -- can still resolve role forms against an identifier-keyed table. Returning the first
    of two definitions there would remap residues against an interface the record may not belong
    to, which is the 0.1.0 defect reached by a different route. The shipped store cannot exercise
    it, because the path needs a caller-supplied SKEMPI table; a synthetic one can."""
    from skempi_foldx import FoldxLookup

    table = _skempi_table(tmp_path, [
        ("9XYZ_A_B", "AB10C"),      # one code, two interface definitions
        ("9XYZ_A_C", "AC20D"),
        ("8ABC_A_B", "AB30E"),      # and one defined exactly once
    ])
    fx = FoldxLookup(skempi_csv=table)

    assert fx._sole_definition("9XYZ") is None, "picked one of two definitions for a bare code"
    sole = fx._sole_definition("8ABC")
    assert sole is not None and sole.pdb == "8ABC"
    assert fx._sole_definition("0000") is None, "invented a definition for an absent code"


def test_a_miss_says_which_of_the_three_causes_it_was():
    """'Not in the store' is the same message for an absent complex, a lower-cased code, and an
    arm that was never loaded. Only the first is what it sounds like."""
    from skempi_foldx import MULTI_POINT, SINGLE_POINT, FoldxLookup

    sp = FoldxLookup(arms=(SINGLE_POINT,))
    with pytest.raises(KeyError, match="multi-point variant and only the single-point arm"):
        sp.require("1JTG", "LI38D,GI32Y")

    mp = FoldxLookup(arms=(MULTI_POINT,))
    with pytest.raises(KeyError, match="single-point mutation and only the multi-point arm"):
        mp.require("1JTG", "LI38D")

    fx = FoldxLookup()
    with pytest.raises(KeyError, match="upper case"):
        fx.require("1brs", "DA52A")
    with pytest.raises(KeyError, match="is not in the store"):
        fx.require("9ZZZ", "DA52A")


def test_consolidate_refuses_to_carry_an_intractable_complex_into_the_store():
    """The store is what a later scope gets enumerated from, so a stale JSON copied into it
    re-arms the compute loop the exclusion exists to keep shut."""
    import tempfile
    from pathlib import Path
    from skempi_foldx import consolidate

    terms = {t: 1.0 for t in TERMS}
    with tempfile.TemporaryDirectory() as tmp:
        src = Path(tmp) / "src"
        src.mkdir()
        (src / "1BRS.json").write_text(json.dumps({"muts": {"DA52A": terms}, "meta": {}}))
        (src / "1KBH.json").write_text(json.dumps({"muts": {"LA930W": terms}, "meta": {}}))
        with pytest.warns(RuntimeWarning, match="1KBH"):
            store, _report = consolidate([src])

    assert "1BRS" in store, "consolidation dropped a complex it should have kept"
    assert "1KBH" not in store, \
        "an intractable complex was carried into the consolidated store"


def test_a_lower_cased_complex_code_still_finds_its_records():
    """Split files and hand-written tables are not reliably upper case, and a whole-file case
    difference would otherwise read as a coverage collapse rather than as a formatting mistake."""
    import subprocess
    import sys as _sys
    import tempfile
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    with tempfile.TemporaryDirectory() as tmp:
        table = Path(tmp) / "rows.tsv"
        table.write_text("1brs_A\tx\tDA52A\n1BRS_A\tx\tDB35A\n")
        out = subprocess.run(
            [_sys.executable, str(root / "experiments" / "coverage_report.py"),
             "--table", str(table), "--label", "mixed_case"],
            capture_output=True, text=True, cwd=str(root))
    assert out.returncode == 0, out.stderr
    line = [l for l in out.stdout.splitlines() if l.startswith("mixed_case")][0]
    assert " 2 " in line and "100.0%" in line, \
        f"a lower-cased code was counted as uncovered: {line!r}"


def test_the_recompute_plan_submits_each_pairings_whole_list(tmp_path):
    """It decides which shipped values get replaced, and a subset list would produce values from
    a mutation list that never existed."""
    import sys as _sys
    from pathlib import Path
    _sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
    from recompute_alternate_interfaces import campaigns

    csv_path = _skempi_table(tmp_path, [
        ("9AAA_B_C", "EC69A"), ("9AAA_B_C", "KB134A"),      # kept pairing: 2
        ("9AAA_B_A", "DA117A"), ("9AAA_B_A", "KB134A"),     # alternate: 2, one shared
        ("9BBB_A_C", "LC5A"), ("9BBB_AB_CD", "LC5A"),       # shares everything -> no campaign
    ])
    plan = campaigns(csv_path)
    assert [c["pdb"] for c in plan] == ["9AAA"], \
        "a complex whose definitions share every mutation should need no campaign"
    only = plan[0]
    assert only["groups"] == ("B", "A")
    assert only["mutations"] == ["DA117A", "KB134A"], \
        "submitted the corrections rather than the pairing's whole list"
    assert only["corrects"] == 1, "counted a shared mutation as a correction"

    both = campaigns(csv_path, include_kept=True)
    assert sorted(c["groups"] for c in both) == [("B", "A"), ("B", "C")], \
        "--include-kept did not add the surviving definition's own list"


def test_coverage_says_so_when_most_of_a_row_list_fails_to_resolve():
    """The naming gap is the larger hazard by two orders of magnitude -- 34.6% of a real split
    against the interface collapse's 0.5% of the store -- so it warrants the louder signal, not
    the quieter one. Below the threshold a handful of genuinely uncomputed rows stays silent."""
    from skempi_foldx import FoldxLookup, load_bundled_store

    # The warning fires for a lookup that could resolve more given more inputs. The shipped store
    # already resolves exactly, so it is silent by design -- this is about every other store.
    import tempfile
    from pathlib import Path as _Path
    with tempfile.TemporaryDirectory() as tmp:
        foreign = _Path(tmp) / "foreign_sp"
        foreign.mkdir()
        terms = {t: 1.0 for t in TERMS}
        (foreign / "1ACB_E_I.json").write_text(json.dumps(
            {"muts": {"LI38D": dict(terms, cleaned="LI38D")}, "meta": {}}))
        bare = FoldxLookup(stores=foreign)
        assert bare.is_exact is False, "a store with no role names cannot be exact"
        role_rows = [("1ACB", f"LB{n}D") for n in range(38, 188)]
        with pytest.warns(RuntimeWarning, match="did not resolve"):
            covered, total, missing = bare.coverage(role_rows)
        assert covered == 0 and total == 150, f"got {covered}/{total}"

    fx = FoldxLookup()
    assert fx.is_exact is True

    # a shortfall too small to be a naming problem must not interrupt: one miss in a hundred is
    # an absent mutation, not a convention mismatch.
    store = load_bundled_store()
    one = next(i for i in store if i == "1BRS" or i.startswith("1BRS_"))
    known = [(one, k) for k in sorted(store[one])]
    mostly_fine = (known * 20)[:120] + [(one, "ZZ999Z")]
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fx.coverage(mostly_fine)

    # and the shipped store resolves exactly, so its shortfalls are absent rows rather than
    # unmatched names -- which is what makes the silence above correct rather than lucky
    assert fx.is_exact


def test_coverage_stays_silent_below_the_threshold_and_when_nothing_more_could_resolve(tmp_path):
    """Three separate gates, each of which has let a real defect through when broken: the row
    floor, the missing-share threshold, and is_exact."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {f"LI{n}D": terms for n in range(1, 101)}, "meta": {}}))
    fx = FoldxLookup(stores=store)
    absent = [("1XYZ", f"LB{n}D") for n in range(500, 700)]
    present = [("1XYZ", f"LI{n}D") for n in range(1, 101)]

    # under 100 rows: silent however bad the shortfall
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fx.coverage(absent[:99])
    with pytest.warns(RuntimeWarning, match="did not resolve"):
        fx.coverage(absent[:100])

    # at least 5% missing, measured against the row count
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        fx.coverage(present * 2 + absent[:10])          # 210 rows, 10 missing = 4.8%
    with pytest.warns(RuntimeWarning, match="did not resolve"):
        fx.coverage(present * 2 + absent[:11])          # 211 rows, 11 missing = 5.2%


def test_the_coverage_warning_asserts_nothing_about_the_rows(tmp_path):
    """Five attempts to classify why a row missed each shipped a claim that was false of rows a
    consumer really had: guessed from chain letters, from the arm's shape, from the code's case.
    The message now describes the lookup only. Every retracted claim is barred by meaning."""
    from skempi_foldx import SINGLE_POINT, FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {f"LI{n}D": terms for n in range(1, 51)}, "meta": {}}))
    fx = FoldxLookup(stores=store, arms=(SINGLE_POINT,))

    # every shape that previously drew a different guess, in one list
    rows = ([("1XYZ", f"LI{n}D,GI{n}Y") for n in range(1, 61)]      # unloaded arm by shape
            + [("1xyz", f"LI{n}D") for n in range(1, 41)]           # lower-cased, held by the store
            + [("1XYZ", f"LB{n}D") for n in range(1, 41)])          # role-chain form
    with pytest.warns(RuntimeWarning) as caught:
        fx.coverage(rows)
    message = str(caught[0].message)

    for claim in ("already matches", "absent from the store", "no input will change",
                  "no campaign indexed here computed them", "arm this lookup did not load",
                  "resolve once the PDB code is upper-cased", "If these labels come from"):
        assert claim not in message, f"asserted a cause it cannot establish, as {claim!r}"
    assert "require(pdb, mutation) gives the reason for any one row" in message, message


def test_the_warning_names_the_input_that_would_widen_this_lookup(tmp_path):
    """Four states, each a fact about the instance rather than a guess about the rows."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {"LI38D": dict(terms, cleaned="LI38D")}, "meta": {}}))
    maps = tmp_path / "maps"
    maps.mkdir()
    rows = [("1XYZ", f"LB{n}D") for n in range(1, 151)]

    def message_for(**kw):
        with pytest.warns(RuntimeWarning) as caught:
            FoldxLookup(stores=store, **kw).coverage(rows)
        return str(caught[-1].message)

    assert "skempi_csv= and mapping_dir= widen it" in message_for(), "no table given"

    naming = tmp_path / "s.csv"
    naming.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        "1XYZ_I_E;LI38D;LI38D;1e-9;1e-9\n"
    )
    assert "mapping_dir= makes its role-chain resolution exact" in message_for(skempi_csv=naming)

    partial = message_for(skempi_csv=naming, mapping_dir=maps)
    assert "0 of 1 mapping files were read" in partial, partial

    # a table naming no complex this store holds: nothing was ever asked for, so naming an input
    # to widen would be a guess about which of the two is the narrow one.
    other = tmp_path / "other.csv"
    other.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        "9ZZZ_A_B;LA38D;LA38D;1e-9;1e-9\n"
    )
    none_wanted = message_for(skempi_csv=other, mapping_dir=maps)
    assert "No complex it indexes appears in the supplied skempi_csv=" in none_wanted, none_wanted
    assert "0 of 0" not in none_wanted, f"blamed a directory nothing was asked of: {none_wanted}"


def _exact_lookup(tmp_path):
    """A lookup that has read every mapping it asked for, so ``is_exact`` is true."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {f"LI{n}D": dict(terms, cleaned=f"LI{n}D") for n in range(1, 11)}, "meta": {}}))
    csv_path = tmp_path / "s.csv"
    csv_path.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        + "".join(f"1XYZ_I_E;LI{n}D;LI{n}D;1e-9;1e-9\n" for n in range(1, 11))
    )
    maps = tmp_path / "maps"
    maps.mkdir()
    (maps / "1XYZ.mapping").write_text(
        "".join(f"LEU I {n}  {n}\n" for n in range(1, 11))
        + "".join(f"GLU E {n}  {n}\n" for n in range(1, 11))
    )
    fx = FoldxLookup(stores=store, skempi_csv=csv_path, mapping_dir=maps)
    assert fx.is_exact, (
        f"fixture is not exact: {fx._mappings_read} of {fx._mappings_wanted} mappings read"
    )
    return fx


def test_coverage_is_silent_once_no_naming_input_could_resolve_more(tmp_path):
    """The is_exact gate. Dropping it entirely once left all tests green, so it is pinned here
    rather than only described."""
    fx = _exact_lookup(tmp_path)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        covered, total, _ = fx.coverage([("1XYZ", f"LI{n}D") for n in range(500, 700)])
    assert covered == 0 and total == 200, f"expected none to resolve, got {covered}/{total}"


def test_coverage_warns_at_the_threshold_itself(tmp_path):
    """`>=` against `>`: the boundary is where an off-by-one hides, and picking 4.8% and 5.2%
    never touches it."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {f"LI{n}D": terms for n in range(1, 191)}, "meta": {}}))
    fx = FoldxLookup(stores=store)
    present = [("1XYZ", f"LI{n}D") for n in range(1, 191)]
    absent = [("1XYZ", f"LB{n}D") for n in range(1, 11)]

    with pytest.warns(RuntimeWarning, match=r"10 of 200 rows \(5.0%\)"):
        fx.coverage(present[:190] + absent)          # exactly 5%


def test_coverage_reports_the_share_that_failed_not_the_share_that_worked(tmp_path):
    """The headline number is the one a reader acts on, and both halves are plausible in place."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1XYZ.json").write_text(json.dumps(
        {"muts": {f"LI{n}D": terms for n in range(1, 21)}, "meta": {}}))
    fx = FoldxLookup(stores=store)
    rows = [("1XYZ", f"LI{n}D") for n in range(1, 21)] + [("1XYZ", f"LB{n}D") for n in range(1, 81)]

    with pytest.warns(RuntimeWarning) as caught:
        fx.coverage(rows)
    message = str(caught[0].message)
    assert "80 of 100 rows (80.0%) did not resolve" in message, message


def test_a_store_written_under_bare_codes_still_resolves_against_an_identifier_table(tmp_path):
    """Every store shipped before 0.2.0 is keyed by PDB code. Reading one against a table keyed
    by identifier must still resolve, or the release strands its own predecessors -- but only
    where the code is unambiguous, since a doubled code has no single pairing to number roles by.
    """
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    store = tmp_path / "mine_sp"
    store.mkdir()
    (store / "1ABC.json").write_text(json.dumps(
        {"muts": {"LI38D": dict(terms, cleaned="LI38D")}, "meta": {}}))
    (store / "2XYZ.json").write_text(json.dumps(
        {"muts": {"LA38D": dict(terms, cleaned="LA38D")}, "meta": {}}))
    csv_path = tmp_path / "s.csv"
    csv_path.write_text(
        "#Pdb;Mutation(s)_PDB;Mutation(s)_cleaned;Affinity_mut (M);Affinity_wt (M)\n"
        "1ABC_E_I;LI38D;LI38D;1e-9;1e-9\n"
        "2XYZ_A_B;LA38D;LA38D;1e-9;1e-9\n"
        "2XYZ_C_D;LC38D;LC38D;1e-9;1e-9\n"
    )
    fx = FoldxLookup(stores=store, skempi_csv=csv_path)
    assert fx.get("1ABC", "LI38D") is not None, "a singly-defined code lost its record"
    assert fx.get("1ABC", "LB38D") is not None, \
        "the role form did not resolve through the sole definition"
    # 2XYZ is defined twice in the table, so the bare-code store entry gets no role alias: the
    # pairing decides the role numbering and there is no basis for choosing one.
    assert fx.get("2XYZ", "LA38D") is not None, "the stored key itself must still resolve"


def test_a_supplied_store_sets_the_arms_it_actually_holds(tmp_path):
    """Otherwise `arms` still describes the bundled default, and every message derived from it
    describes a store this instance never opened."""
    from skempi_foldx import MULTI_POINT, SINGLE_POINT, FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    sp = tmp_path / "mine_sp"
    sp.mkdir()
    (sp / "1XYZ.json").write_text(json.dumps({"muts": {"LI38D": terms}, "meta": {}}))

    fx = FoldxLookup(stores=sp)
    assert fx.arms == (SINGLE_POINT,), f"arms should describe what was loaded, got {fx.arms}"
    with pytest.raises(KeyError, match="multi-point variant and only the single-point arm"):
        fx.require("1XYZ", "LI38D,GI32Y")
    with pytest.raises(KeyError, match=r"arms loaded: results_sp\b"):
        fx.require("9ZZZ", "LI38D")

    mp = tmp_path / "mine_mp"
    mp.mkdir()
    (mp / "1XYZ.json").write_text(json.dumps({"variants": {"LI38D,GI32Y": terms}, "meta": {}}))
    assert FoldxLookup(stores=[sp, mp]).arms == (SINGLE_POINT, MULTI_POINT)


def test_a_record_that_declares_its_author_form_is_trusted_over_the_table(tmp_path):
    """A `cleaned` field IS the author form, so nothing is inferred and SKEMPI membership is
    irrelevant. Requiring membership discards correct aliases whenever the supplied table is a
    subset of the one that built the store -- silently, and by nine points of coverage."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    # The table lists LI38G but not LI38D; the store holds a record declaring LI38D via `cleaned`.
    csv_path = _skempi_table(tmp_path, [("1ACB_E_I", "LI38G")])
    store = {"1ACB": {"REC1": dict(terms, cleaned="LI38D")}}

    fx = FoldxLookup(stores=store, skempi_csv=csv_path)
    assert fx.get("1ACB", "LI38D") is not None, \
        "a record's own `cleaned` field stopped being honoured when the table lacked the mutation"


def test_a_key_that_is_only_inferred_to_be_an_author_form_must_be_listed(tmp_path):
    """Without `cleaned` the key is a guess at an author form. Remapping an unlisted guess is how
    a role-keyed store aliases a real record onto a residue it does not describe."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    csv_path = _skempi_table(tmp_path, [("1ACB_E_I", "LI38D")])
    mappings = _mapping_dir(tmp_path, "1ACB", {"E": 245, "I": 70})

    listed = FoldxLookup(stores={"1ACB": {"LI38D": dict(terms)}},
                         skempi_csv=csv_path, mapping_dir=mappings)
    assert listed.get("1ACB", "LI38D") is not None

    unlisted = FoldxLookup(stores={"1ACB": {"LB38D": dict(terms)}},
                           skempi_csv=csv_path, mapping_dir=mappings)
    assert unlisted.get("1ACB", "LB38D") is not None, "the stored key should still resolve"
    assert len(unlisted._alias) == 0, \
        "an unlisted key was inferred to be an author form and remapped anyway"


def test_a_complex_counts_as_wanting_a_mapping_even_when_no_role_form_is_built(tmp_path):
    """Counting downstream of the refusals lets exactness be claimed over precisely the complexes
    that got no role resolution at all."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    csv_path = _skempi_table(tmp_path, [("1ACB_E_I", "LI38D")])
    empty = tmp_path / "no_mappings"
    empty.mkdir()

    fx = FoldxLookup(stores={"1ACB": {"LB38D": dict(terms)}},
                     skempi_csv=csv_path, mapping_dir=empty)
    assert fx._mappings_wanted == 1, "a complex whose records were all refused was not counted"
    assert fx.is_exact is False


def test_a_later_store_cannot_shadow_an_earlier_alias(tmp_path):
    """First-wins has to hold across the key/alias boundary: a store keyed by author form
    collides with an earlier store's `cleaned` alias, not with its key."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    first = {"1ACB": {"LB38D": dict(terms, cleaned="LI38D")}}
    second = {"1ACB": {"LI38D": {t: -999.0 for t in TERMS}}}

    with pytest.warns(RuntimeWarning, match="resolving to more than one record"):
        fx = FoldxLookup(stores=[first, second])
    assert fx.get("1ACB", "LI38D")["Interaction Energy"] == 1.0, \
        "the second store shadowed the first store's alias for the same mutation"
    assert fx.get("1ACB", "LB38D") is fx.get("1ACB", "LI38D")


def test_a_non_finite_term_never_reaches_a_feature_vector():
    """json.loads accepts bare NaN and Infinity, so a results directory can carry them all the way
    into a model, where they are far harder to trace back."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    fx = FoldxLookup(stores={"1ZZZ": {"LI38D": dict(terms, **{TERMS[1]: float("inf")})}})
    with pytest.raises(ValueError, match=TERMS[1]):
        fx.vector("1ZZZ", "LI38D")


def test_a_collision_is_announced_whichever_order_the_stores_arrive_in(tmp_path):
    """Detection tested the incoming key against existing aliases but never the reverse, so the
    same disagreement was caught in one order and passed over in the other."""
    from skempi_foldx import FoldxLookup

    role_keyed = {"1ACB": {"LB38D": dict({t: 1.0 for t in TERMS}, cleaned="LI38D")}}
    author_keyed = {"1ACB": {"LI38D": {t: -999.0 for t in TERMS}}}

    for order, label in (([role_keyed, author_keyed], "role-keyed first"),
                         ([author_keyed, role_keyed], "author-keyed first")):
        with pytest.warns(RuntimeWarning, match="resolving to more than one record"):
            fx = FoldxLookup(stores=order)
        assert fx.get("1ACB", "LI38D") is not None, label


# ----------------------------------------------------------------- malformed SKEMPI rows
def test_malformed_rows_are_mutation_scoped_not_complex_scoped():
    """Excluding the complex to drop two rows would discard 33 sound 2C5D records."""
    from skempi_foldx.exclusions import MALFORMED_ROWS, INTRACTABLE, is_malformed
    assert all(pdb == "2C5D" for pdb, _ in MALFORMED_ROWS)
    # the complex itself must NOT be intractable -- the rest of it is fine
    assert "2C5D" not in INTRACTABLE
    assert is_malformed("2C5D", "RA32E,KA34E,RB32E,KA34E")
    assert not is_malformed("2C5D", "RA32E,RB32E")


def test_malformed_rows_repeat_a_substitution_and_intend_the_symmetric_partner():
    from skempi_foldx.exclusions import MALFORMED_ROWS
    for row in MALFORMED_ROWS.values():
        parts = row.cleaned.split(",")
        assert len(parts) != len(set(parts)), "shipped form must contain the repeat"
        intended = row.intended.split(",")
        assert len(intended) == len(set(intended)), "intended form must not"
        assert len(parts) == len(intended)
        # exactly one position differs, and it is the repeated one
        assert sum(a != b for a, b in zip(parts, intended)) == 1


def test_malformed_ignores_the_intractable_bypass(monkeypatch):
    """ALLOW_ENV re-tests a complex on a newer FoldX; no engine fixes a bad source row."""
    from skempi_foldx.exclusions import ALLOW_ENV, is_malformed, is_excluded
    monkeypatch.setenv(ALLOW_ENV, "1")
    assert is_excluded("1KBH") is None          # bypass works for the compute guard
    assert is_malformed("2C5D", "RA32E,KA34E,RB32E,KA34E")   # and not for this one


def test_filter_worklist_keeps_malformed_rows_and_drops_intractable_complexes():
    """The two registries act at different stages: one guards compute, the other guards trust.

    A malformed row survives the worklist because FoldX computes it correctly -- what is wrong is
    the arity its label claims, which no compute decision can repair. Removing it would also
    change the values of every entry after it in the list.
    """
    from skempi_foldx.exclusions import filter_worklist
    wl = {"2C5D": ["RA32E,RB32E", "RA32E,KA34E,RB32E,KA34E"], "1BRS": ["DA52A"], "1KBH": ["LA1A"]}
    out = filter_worklist(wl, on_note=None)
    assert out["2C5D"] == ["RA32E,RB32E", "RA32E,KA34E,RB32E,KA34E"]
    assert out["1BRS"] == ["DA52A"]
    assert "1KBH" not in out, "the intractable complex is still dropped whole"


def test_announce_malformed_names_the_rows_without_removing_them():
    from skempi_foldx.exclusions import announce_malformed
    seen = []
    muts = ["RA32E,RB32E", "RA32E,KA34E,RB32E,KA34E"]
    with pytest.warns(RuntimeWarning, match="MALFORMED"):
        found = announce_malformed("2C5D", muts, on_note=seen.append)
    assert found == ["RA32E,KA34E,RB32E,KA34E"]
    assert muts == ["RA32E,RB32E", "RA32E,KA34E,RB32E,KA34E"], "the input list is not mutated"
    assert len(seen) == 1 and "shipped" in seen[0]


# ---------------------------------------------------------------- the 0.2.0 re-key, and its edges
# Every test below pins a defect that shipped green: the suite passed while the behaviour was
# wrong, because nothing exercised the shape that breaks it.

def test_shorten_path_leaves_no_username_at_any_depth():
    """A campaign records where it seeded its repair from, and that lands in package data.

    The leak this prevents is not the obvious one. An absolute path is easy to grep for; what
    survives a naive truncation is a *relative* fragment with the username at its head, which no
    absolute-path pattern matches because the leading slash is gone.
    """
    from skempi_foldx.store import shorten_path

    user, u, h = "someuser", "Users", "home"
    for raw in (f"/{u}/{user}/scratch/camp/2C5D/f.pdb",
                f"/{h}/{user}/2C5D/f.pdb",
                f"/data/{h}/{user}/2C5D/f.pdb",           # home is not at position 0
                f"//fileserver/{h}/{user}/2C5D/f.pdb",
                f"../../{h}/{user}/2C5D/f.pdb",
                f"/{u}/{user}/scratch/{user}/2C5D/f.pdb",     # name recurs deeper
                f"C:\\{u}\\{user}\\camp\\2C5D\\f.pdb",         # backslashes
                f"/{u}/{user}"):                              # nothing but the home
        assert user not in shorten_path(raw), raw

    # `/root` IS the home, so only one component goes -- stripping two would eat the complex.
    # Assembled rather than written out, like the fixtures above: .github/leakscan.sh scans this
    # file, and a literal home-directory path here fails the release gate it exists to protect.
    bare_home = "/".join(["", "root", "camp", "2C5D", "f.pdb"])
    assert shorten_path(bare_home) == "camp/2C5D/f.pdb"
    assert shorten_path("") == ""
    already = "2C5D_A_C_mp/2C5D/2C5D_Repair.pdb"
    assert shorten_path(already) == already, "an already-safe path was mangled"
    assert shorten_path(shorten_path(f"/Users/{user}/a/b/c/f.pdb")).count("/") <= 2


def test_the_shipped_store_carries_no_path_at_all():
    """The data is the deliverable, so it is the last place a local path may survive -- and the
    place a source-only hygiene grep does not look."""
    import re
    from skempi_foldx import bundled_path, MULTI_POINT

    # Assembled rather than written out, so this file does not itself trip the CI hygiene grep
    # that greps *.py. A gate that has to skip its own tests is a gate with a blind spot.
    sep = "/"
    bad = re.compile("|".join([f"{sep}Users{sep}[a-z]", f"{sep}home{sep}[a-z]",
                               f"{sep}private{sep}tmp{sep}", r"[A-Za-z]:\\\\"]))
    for which in (bundled_path(), bundled_path(MULTI_POINT)):
        for path in which.glob("*.json"):
            assert not bad.search(path.read_text()), path.name


def test_pool_by_code_conserves_mutations_and_refuses_a_bare_string():
    from skempi_foldx import pool_by_code

    pooled = pool_by_code({"3SE4_B_A": ["DA117A", "SHARED"], "3SE4_B_C": ["EC69A", "SHARED"],
                           "1BRS_A_D": ["DA52A"]})
    assert pooled == {"3SE4": ["DA117A", "EC69A", "SHARED"], "1BRS": ["DA52A"]}
    assert pool_by_code({}) == {}
    assert pool_by_code({"1BRS": ["DA52A"]}) == {"1BRS": ["DA52A"]}, "not idempotent on codes"
    # A str satisfies Sequence[str] and would be pooled character by character.
    with pytest.raises(TypeError):
        pool_by_code({"1BRS_A_D": "DA52A"})


def test_two_spellings_of_one_record_pool_rather_than_one_being_dropped(tmp_path):
    """A worklist can name one definition by code and by identifier at once. Running it twice
    would collide; dropping one silently computes a *subset*, which this module's whole subject
    is about -- the same mutation in a shorter list is a different number, by up to 11 kcal/mol.
    """
    from skempi_foldx import FoldxConfig, MODE_AUTHOR, SkempiComplex, run_campaign

    entry = SkempiComplex("1ABC", "A", "B", single={"LA1A", "LA2A", "LA3A"})
    notes = []
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv",
                         binary="/bin/echo")
    (tmp_path / "p").mkdir()
    run_campaign({"1ABC": ["LA1A"], "1ABC_A_B": ["LA2A", "LA3A"]}, {"1ABC_A_B": entry},
                 config, mode=MODE_AUTHOR, jobs=1, on_result=notes.append)

    submitted = [n for n in notes if "complexes," in n][0]
    assert "1 complexes, 3 mutations" in submitted, submitted
    assert any("pooled into one run" in n for n in notes)
    # The identifier wins the name: a bare code cannot say which pairing it means.
    assert (tmp_path / "r" / "1ABC_A_B.json").exists() or not (tmp_path / "r" / "1ABC.json").exists()


def test_an_empty_mutation_list_only_runs_under_repair_only(tmp_path):
    """`{pdb: []}` means "repair and stop". Outside a repair-only run it would pay a full
    RepairPDB and then write a results JSON that marks the complex done for every later resume.
    """
    from skempi_foldx import FoldxConfig, MODE_AUTHOR, SkempiComplex, run_campaign

    entry = SkempiComplex("1ABC", "A", "B", single={"LA1A"})
    (tmp_path / "p").mkdir()
    config = FoldxConfig(results_dir=tmp_path / "r", work_dir=tmp_path / "w",
                         pdb_dir=tmp_path / "p", skempi_csv=tmp_path / "s.csv",
                         binary="/bin/echo")
    notes = []
    run_campaign({"1ABC_A_B": []}, {"1ABC_A_B": entry}, config, mode=MODE_AUTHOR, jobs=1,
                 on_result=notes.append)
    assert any("no mutations to build" in n for n in notes)
    assert not list((tmp_path / "r").glob("*.json")), "wrote a done-marker for an unbuilt complex"

    # ...and a repair-only campaign still gets its empty list through.
    notes.clear()
    run_campaign({"1ABC_A_B": []}, {"1ABC_A_B": entry}, config, mode=MODE_AUTHOR, jobs=1,
                 repair_only=True, on_result=notes.append)
    assert any("1 complexes, 0 mutations" in n for n in notes), notes


def test_only_selects_by_code_as_well_as_by_identifier(tmp_path):
    """Callers hold structures to run, not SKEMPI identifiers. Matching on the key alone made
    `only=["1KBH"]` select nothing and the campaign exit successfully having done nothing."""
    from skempi_foldx import load_skempi, worklist_single_point

    csv = _skempi_table(tmp_path, [("1AAA_E_I", "LI38S"), ("2BBB_A_B", "YC5F")])
    sk = load_skempi(csv)
    assert set(worklist_single_point(sk, only=["1AAA"])) == {"1AAA_E_I"}
    assert set(worklist_single_point(sk, only=["1AAA_E_I"])) == {"1AAA_E_I"}
    assert worklist_single_point(sk, only=["1AA"]) == {}, "matched on a prefix"


def test_two_records_of_one_complex_claiming_one_alias_are_announced(tmp_path):
    """Whichever key JSON iteration reached first won the name, and the other record's energies
    answered to it. Both alias branches are symmetric and both were silent."""
    from skempi_foldx import FoldxLookup

    terms = {t: 1.0 for t in TERMS}
    for label, recs in (
        ("role", {"LI38D": dict(terms, cleaned="LI38D", role="LB38D"),
                  "LJ38D": dict(terms, cleaned="LJ38D", role="LB38D")}),
        ("cleaned", {"LI38D": dict(terms, cleaned="ZZ9Z"),
                     "LJ38D": dict(terms, cleaned="ZZ9Z")}),
    ):
        with pytest.warns(RuntimeWarning, match="resolving to more than one record"):
            FoldxLookup(stores={"9ZZZ_A_B": recs})


def test_exactness_is_a_property_of_every_indexed_store(tmp_path):
    """A store with no role names is still an estimate, and indexing one beside the bundled store
    must not let the bundled store vouch for it. An *empty* store must not do the reverse."""
    from skempi_foldx import FoldxLookup, load_bundled_store, load_store

    terms = {t: 1.0 for t in TERMS}
    roleless = tmp_path / "roleless_sp"
    roleless.mkdir()
    (roleless / "9ZZZ_A_B.json").write_text(json.dumps(
        {"muts": {"LA38D": dict(terms, cleaned="LA38D")}, "meta": {}}))

    assert FoldxLookup().is_exact is True
    assert FoldxLookup(stores=[load_bundled_store(), {}]).is_exact is True, \
        "an empty store flipped exactness off with nothing in it"
    assert FoldxLookup(stores=[load_bundled_store(), load_store(roleless)]).is_exact is False, \
        "a roleless store borrowed the bundled store's exactness"


def test_write_store_records_the_structure_in_pdb_not_the_identifier(tmp_path):
    """`meta.pdb` names the structure; `skempi_id` names the record. A campaign written before the
    two were separated put the identifier in both, and consolidating it would carry that through
    into the field whose only job is to say which structure the energies came from."""
    from skempi_foldx import write_store

    terms = {t: 1.0 for t in TERMS}
    store = {"2C5D_AB_CD": {"LA1A": dict(terms)}}
    for supplied, expected in (({"pdb": "2C5D_AB_CD"}, "2C5D"),   # pre-split campaign
                               ({"pdb": "2C5D"}, "2C5D"),         # already a code
                               ({}, "2C5D")):                     # absent
        out = tmp_path / f"out_{expected}_{len(supplied)}"
        write_store(store, out, kind="muts", metas={"2C5D_AB_CD": supplied})
        meta = json.loads((out / "2C5D_AB_CD.json").read_text())["meta"]
        assert meta["pdb"] == expected, (supplied, meta)
        assert meta["skempi_id"] == "2C5D_AB_CD"


def test_write_store_shortens_a_path_a_campaign_left_in_meta(tmp_path):
    """consolidate()/write_store is the exported way to materialise a store, so a consumer running
    it over pre-fix campaign output would otherwise reproduce the leak this release removed."""
    from skempi_foldx import write_store

    terms = {t: 1.0 for t in TERMS}
    leaky = "/".join(["", "Users", "someuser", "camp", "1ABC", "1ABC_Repair.pdb"])
    write_store({"1ABC_A_B": {"LA1A": dict(terms)}}, tmp_path / "out", kind="muts",
                metas={"1ABC_A_B": {"repair_seeded_from": leaky}})
    meta = json.loads((tmp_path / "out" / "1ABC_A_B.json").read_text())["meta"]
    assert "someuser" not in meta["repair_seeded_from"]
    assert meta["repair_seeded_from"] == "camp/1ABC/1ABC_Repair.pdb"


def test_a_role_name_and_an_author_name_reach_one_bundled_record():
    """Both spellings of a mutation resolve, and to the same record rather than to two copies.

    The role form is what split files carry, and `1PPF`'s inhibitor chain is `I` by author name
    and `B` by role, so `GI32Y` and `GB32Y` are one substitution written two ways. Asserted on the
    bundled store specifically: the `role` field is computed at build time from the chain mappings,
    so it is a property of what shipped rather than of the loader.

    Identity, not equality -- two distinct records that happened to carry identical terms would
    satisfy `==` while meaning the alias had been lost. `FoldxLookup.get` returns the stored dict,
    which is what makes `is` the available check; a loader that started returning a copy would
    need this rewritten rather than deleted.
    """
    from skempi_foldx import FoldxLookup

    fx = FoldxLookup()
    author, role = fx.get("1PPF", "GI32Y"), fx.get("1PPF", "GB32Y")
    assert author is not None, "the author-chain spelling did not resolve"
    assert author is role, "the role-chain spelling reached a different record, or none"
    assert fx.is_exact, "the shipped role names did not survive packaging"


def test_an_inferred_author_form_must_be_one_skempi_lists(tmp_path):
    """A role-keyed foreign store must not have its keys taken for author forms.

    A record without `cleaned` does not declare which convention its key is in, so the key is
    *inferred* to be an author form and a role name is built from it. If the key was already a role
    name, that builds a second name out of a name that was never an author form, and installs it as
    an alias onto a record it does not describe.

    The conditions are specific, which is why nothing caught this: the role letter must also lie
    inside the author groups, so group1 needs more than one chain, and `mapping_dir=` must be
    supplied or no name can be built at all. `2C5D_AB_CD` satisfies the first — role chain `B` is
    also author chain `B` — so the role-shaped key `KB178E` is re-read as author chain B residue
    178 and remapped to `KA378E`, which names chain A residue 378 and belongs to nothing here.
    """
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("2C5D_AB_CD", "KC178E")])
    maps = _mapping_dir(tmp_path, "2C5D", {"A": 200, "B": 200, "C": 200, "D": 200})
    store = {"2C5D_AB_CD": {"KB178E": {t: 1.0 for t in TERMS}}}   # no `cleaned`

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fx = FoldxLookup(skempi_csv=str(csv), mapping_dir=str(maps), stores=[store])

    fabricated = sorted(name for pdb, name in fx._alias
                        if pdb == "2C5D_AB_CD" and name != "KB178E")
    assert not fabricated, (
        f"a role form was built from a key that is not an author form SKEMPI lists: {fabricated}"
    )
    assert fx.get("2C5D_AB_CD", "KB178E") is not None, "the literal key must still resolve"


def test_a_role_form_is_built_from_the_records_own_definition(tmp_path):
    """A role name is only meaningful under the chain grouping the record belongs to.

    Role letters are assigned by which group a chain sits in, and positions are renumbered by the
    cumulative offset along that group -- so the same author mutation has different role names
    under different pairings of the same PDB code. `_sole_definition` is therefore asked about the
    record's *identifier*. Asking about its bare code instead lets a record keyed to one pairing be
    resolved against another, and the name comes out describing a different residue.

    Here the table knows only `1ACB_I_AE`, which puts chain I in group1 (role letter `A`), while
    the record is keyed `1ACB_E_I`, where chain I is group2 (role letter `B`). Falling back by code
    installs `LA38D` on a record whose own pairing would never produce that name.
    """
    from skempi_foldx import FoldxLookup

    csv = _skempi_table(tmp_path, [("1ACB_I_AE", "LI38D")])
    maps = _mapping_dir(tmp_path, "1ACB", {"A": 100, "E": 245, "I": 70})
    store = {"1ACB_E_I": {"LI38D": dict({t: 1.0 for t in TERMS}, cleaned="LI38D")}}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fx = FoldxLookup(skempi_csv=str(csv), mapping_dir=str(maps), stores=[store])

    invented = sorted(name for pdb, name in fx._alias
                      if pdb == "1ACB_E_I" and name != "LI38D")
    assert not invented, (
        f"a role name was built from another pairing of the same code: {invented}"
    )


def test_the_first_store_supplied_wins_a_conflict():
    """Source precedence is the documented contract, shared with `consolidate()`.

    The conflict is announced and listed either way, so a reversal is invisible in the warning and
    in `conflicts`; only the returned energies change.
    """
    from skempi_foldx import FoldxLookup

    preferred = {"9ZZZ_A_B": {"LI38D": {t: 1.0 for t in TERMS} | {"cleaned": "LI38D"}}}
    fallback = {"9ZZZ_A_B": {"LI38D": {t: 2.0 for t in TERMS} | {"cleaned": "LI38D"}}}
    with pytest.warns(RuntimeWarning, match="resolving to more than one record"):
        fx = FoldxLookup(stores=[preferred, fallback])
    assert fx.get("9ZZZ_A_B", "LI38D")[SCALAR_TERM] == 1.0, "the first store must win"
    with pytest.warns(RuntimeWarning, match="resolving to more than one record"):
        rev = FoldxLookup(stores=[fallback, preferred])
    assert rev.get("9ZZZ_A_B", "LI38D")[SCALAR_TERM] == 2.0, "order must decide, not chance"


def test_the_shipped_store_holds_no_excluded_complex():
    """The deliverable is checked against the package's own exclusion registry.

    `consolidate()` refuses to build a store containing one, and that refusal is tested -- but
    nothing asserted it of the artifact that actually ships. A file named for an excluded complex
    would leave every count unchanged and every other test green.
    """
    from skempi_foldx import load_bundled_store
    from skempi_foldx.exclusions import is_excluded

    for arm in ("results_sp", "results_mp"):
        offenders = sorted(k for k in load_bundled_store(arm) if is_excluded(k))
        assert not offenders, f"{arm} ships excluded complexes: {offenders}"


#: sha256 over every shipped (arm, identifier, mutation, term, value), sorted. Recompute with
#: `experiments/` disabled and a clean tree; changing it is a deliberate act, not a fixup.
PINNED_ENERGY_DIGEST = "7d0c4a7e17a35c3d3c729d8829753bd18ddbc2681a07d9c60cb19d1e51a75c6d"


def test_the_shipped_energies_are_pinned():
    """A version names a measurement, so the values themselves are pinned, not only the shape.

    Counts, keys and field consistency are asserted exhaustively elsewhere. None of that notices a
    rebuild or a backfill that perturbs a number, which is the one change this package's premise
    says must never happen silently within a version.
    """
    import hashlib

    from skempi_foldx import load_bundled_store

    h = hashlib.sha256()
    for arm in ("results_sp", "results_mp"):
        store = load_bundled_store(arm)
        for ident in sorted(store):
            for key in sorted(store[ident]):
                rec = store[ident][key]
                h.update(f"{arm}\0{ident}\0{key}".encode())
                for term in TERMS:
                    h.update(f"\0{term}={rec[term]!r}".encode())
    assert h.hexdigest() == PINNED_ENERGY_DIGEST, (
        "the shipped energy values changed. If that is intended, this release needs a new version "
        "and a CHANGELOG entry saying what moved -- a version names a measurement here."
    )
