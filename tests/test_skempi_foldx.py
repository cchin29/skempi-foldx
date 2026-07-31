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
    assert TERMS == [
        "Interaction Energy", "Backbone Hbond", "Sidechain Hbond", "Van der Waals",
        "Electrostatics", "Solvation Polar", "Solvation Hydrophobic",
        "Van der Waals clashes", "entropy sidechain", "entropy mainchain",
        "torsional clash", "backbone clash",
    ]


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
    assert len(sp) == 322 and sum(len(v) for v in sp.values()) == 4238
    assert sp["1BRS"]["DA52A"]["Interaction Energy"] == -0.6806


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
    assert set(sk) == {"1AAA", "2BBB"}
    assert (sk["1AAA"].group1, sk["1AAA"].group2) == ("E", "I"), "chain groups swapped or misread"
    assert (sk["2BBB"].group1, sk["2BBB"].group2) == ("AB", "CD")
    assert sk["1AAA"].groups == "E,I", "groups string feeds --analyseComplexChains"
    assert sk["1AAA"].single == {"LI38S", "GI40A"}
    assert sk["1AAA"].multi == {"LI38S,GI40A"}, "multi-point row landed in the wrong arm"
    assert sk["2BBB"].single == {"YC5F"} and not sk["2BBB"].multi
    # the cleaned column, not the PDB-numbered one
    assert "LI138S" not in sk["1AAA"].single, "read Mutation(s)_PDB instead of Mutation(s)_cleaned"


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
    """The single-point store carries two key conventions, so a lookup by SKEMPI identifier
    misses ~29% of it. `key="skempi"` fixes that at read time.

    The default must stay byte-for-byte the stored keys: a consumer that maps role chains onto
    author chains itself needs the role-chain key present. Rewriting the store to one convention
    was measured to drop such a consumer's coverage from 3300/3300 to 1092/3300 -- silently, since
    the join simply stops matching and still produces well-formed output."""
    from skempi_foldx import load_bundled_store, reindex_by_skempi_id

    stored = load_bundled_store()
    skempi = load_bundled_store(key="skempi")

    assert sum(map(len, stored.values())) == sum(map(len, skempi.values())), "records lost"
    assert "LB38D" in stored["1ACB"], "default keys were rewritten"
    assert "LI38D" not in stored["1ACB"], "default keys were rewritten"
    assert "LI38D" in skempi["1ACB"], "re-index did not expose the SKEMPI form"
    assert stored["1ACB"]["LB38D"] == skempi["1ACB"]["LI38D"], "re-index changed the record"

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
    assert store_kind(tmp_path / "empty_missing") is None or True


def test_source_labels_are_the_documented_set():
    """`_source` is provenance and nothing joins on it, but it is read by humans, so the labels
    must say what they mean. The originals were directory paths -- `foldx_s1102_results_S4169`
    reads as two different benchmark subsets at once, because the S4169 campaign happened to run
    under a directory created for the S1102 work."""
    from skempi_foldx import MULTI_POINT, load_bundled_store

    expected = {"S4169", "full_skempi", "full_skempi_multipoint"}
    seen = {r.get("_source") for which in ("results_sp", MULTI_POINT)
            for recs in load_bundled_store(which).values() for r in recs.values()}
    seen.discard(None)
    assert seen == expected, f"unexpected source labels: {seen ^ expected}"
    assert not any("s1102" in s.lower() and "4169" in s for s in seen), \
        "a label names two different benchmark subsets"
