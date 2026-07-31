#!/usr/bin/env bash
# SP arm of the FoldX repair-count sweep.
#
# --store: NOT the driver's default (ROOT/data/foldx/results_sp). On a compute box that path may
# not exist, load_store returns {}, and `--scope all` would silently run zero complexes. Point it
# at the campaign store that actually has results.
#
# 1KBH is not named here: it lives in skempi_foldx/exclusions.py, which every entry point consults,
# so it cannot be reintroduced by forgetting a flag -- 2120 residues / 33k atoms, against 26k for
# 3VR6; its RepairPDB reached residue A36 of 2120 in 10 h 24 m with vdW clash energies over
# 1.2M kcal/mol. It reaches no split under data/, and its 3 SP + 83 MP entries are never scored.
# Use SWEEP_EXCLUDE for anything not yet in the registry.
set -uo pipefail

# Repo root: two levels up from this script. Override any of these in the environment.
ROOT=${SKEMPI_FOLDX_ROOT:-"$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"}
STORE=${SWEEP_STORE:-"$ROOT/scratch/foldx_skempi_full/results_all"}
PDBS=${SKEMPI_PDBS:-"$ROOT/scratch/skempi2/PDBs"}
SKEMPI=${SKEMPI_CSV:-"$ROOT/scratch/skempi_v2.csv"}
JOBS=${JOBS:-36}
EXCLUDE=${SWEEP_EXCLUDE:-""}
: "${FOLDX_BIN:?set FOLDX_BIN to the FoldX 5 binary}"
export PYTHONPATH="$ROOT"

[ -d "$STORE" ] || { echo "store $STORE does not exist -- --scope all would run zero complexes"; exit 1; }

# shellcheck disable=SC2086   # word-splitting is intended: EXCLUDE is a list of ids
exec python3 "$ROOT/experiments/repair_ablation.py" sweep \
    --iterations 4 --rounds 1 2 3 4 \
    --scope all --jobs "$JOBS" ${EXCLUDE:+--exclude $EXCLUDE} \
    --store      "$STORE" \
    --pdb-dir    "$PDBS" \
    --skempi-csv "$SKEMPI"
