#!/usr/bin/env bash
# Full post-SP chain: verify -> MP build -> verify -> record provenance.
# Every step is idempotent and safe to re-run by hand.
#
# The point of the ordering: the integrity checks run BEFORE the work*/ trees are deleted,
# because individual_list.txt lives there and it is the only record of what BuildModel actually
# consumed. Delete work*/ first and the run becomes permanently unverifiable. Those trees look
# like disposable scratch and are not until list_hashes.py has run over them.
#
# EXCLUDE is extra ids on top of skempi_foldx.exclusions.INTRACTABLE, which both scripts apply on
# their own. residue_check.py needs the exclusions so an excluded complex is not counted as a
# missing chain (which would pin full_coverage=false forever); sweep_mp.py build needs them
# because its incomplete-chain guard is a hard sys.exit by design.
set -euo pipefail

# Repo root: two levels up from this script. Override any of these in the environment.
ROOT=${SKEMPI_FOLDX_ROOT:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"}
RUN=${SWEEP_RUN:-"$ROOT/scratch/foldx_repair_sweep"}
PDBS=${SKEMPI_PDBS:-"$ROOT/scratch/skempi2/PDBs"}
BASE="$RUN/scratch/foldx_repair_ablation/sweep_4x"
EXCLUDE=${SWEEP_EXCLUDE:-""}
JOBS=${JOBS:-36}
: "${FOLDX_BIN:?set FOLDX_BIN to the FoldX 5 binary}"
export PYTHONPATH="$ROOT"
cd "$(dirname "${BASH_SOURCE[0]}")"
say(){ echo "[chain $(date '+%F %T')] $*"; }

say "waiting for the SP sweep to exit (extra exclusions: ${EXCLUDE:-none})..."
while pgrep -f 'repair_ablation\.py sweep' >/dev/null; do sleep 60; done
say "SP arm finished."

# shellcheck disable=SC2086   # word-splitting is intended: EXCLUDE is a list of ids
say "residue check across the sweep universe (needs chain_work/ only)"
python3 residue_check.py --exclude $EXCLUDE 2>&1 | sed 's/^/    /'

say "hashing the SP individual_list.txt files BuildModel consumed"
python3 list_hashes.py verify --write-manifest --inject 2>&1 | sed 's/^/    /'

say "starting MP BuildModel"
python3 sweep_mp.py build --rounds 1 2 3 4 --jobs "$JOBS" --exclude $EXCLUDE 2>&1 | sed 's/^/    /'

say "final integrity pass over both arms"
python3 residue_check.py --exclude $EXCLUDE 2>&1 | sed 's/^/    /'
python3 list_hashes.py verify --write-manifest --inject 2>&1 | sed 's/^/    /'

say "chain-input audit: every repair input must still equal the pristine source"
BASE="$BASE" PDBS="$PDBS" python3 - <<'PY' 2>&1 | sed 's/^/    /'
import filecmp, os, sys
from pathlib import Path

# WHICH FILE CARRIES THE INVARIANT depends on which repair implementation wrote the tree.
#
# A repair implementation that repairs `<pdb>.pdb` in place keeps a write-once
# `<pdb>_original.pdb` backup, and that backup carries the invariant. skempi_foldx/run.py does
# something stronger: it never mutates the input at all -- the chain runs in `<pdb>_chain.pdb`, so
# `<pdb>.pdb` IS the pristine copy and `_original.pdb` is never created.
#
# Auditing only for `_original.pdb` against a tree this repo wrote finds zero files and reports
# clean by matching nothing, which is indistinguishable from actually clean. So: check whichever
# backup exists, name it, and fail loudly when a chain dir has neither or when nothing was
# audited at all.
B, P = Path(os.environ["BASE"]), Path(os.environ["PDBS"])
work = B / "chain_work"
if not work.is_dir():
    sys.exit(f"no chain_work under {B} -- nothing to audit (refusing to report clean)")

bad, checked, none = [], 0, []
for d in sorted(p for p in work.iterdir() if p.is_dir()):
    pristine = P / f"{d.name}.pdb"
    if not pristine.exists():
        continue
    original, working = d / f"{d.name}_original.pdb", d / f"{d.name}.pdb"
    ref = original if original.exists() else (working if working.exists() else None)
    if ref is None:
        none.append(d.name)
        continue
    checked += 1
    if not filecmp.cmp(ref, pristine, shallow=False):
        bad.append(f"{d.name}({ref.name})")

print(f"chain inputs audited: {checked}   corrupted: {len(bad)} {bad[:20]}")
if none:
    print(f"WARNING {len(none)} chain dirs have NEITHER _original.pdb nor <pdb>.pdb: {none[:20]}")
if not checked:
    sys.exit("audited 0 inputs -- treat this as FAILED, not clean")
sys.exit(1 if bad or none else 0)
PY

say "DONE. Ship round_*/results{,_mp}/ plus list_manifest.json and residue_check.json."
say "Do NOT delete work*/ before the manifest exists -- it is the only record of the lists."
