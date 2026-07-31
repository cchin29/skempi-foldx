"""Where FoldX, its inputs, and its outputs live.

These are carried in a config object passed explicitly to each task, not in module-level globals.
Globals a driver assigns at import time (``bd.WORK = ...; bd.RESULTS = ...``) work, but they make
the engine's behaviour depend on which module last assigned to them, and every driver then needs
its own process-pool initializer to re-apply them inside each worker.

Passing a config means a worker process gets it as an argument, and two campaigns can run in the
same interpreter without interfering.
"""

from __future__ import annotations

import os
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence


@dataclass
class FoldxConfig:
    """Paths and binary location for a FoldX campaign."""

    results_dir: Path
    """Per-complex ``<pdb>.json`` output."""

    work_dir: Path
    """FoldX scratch: repaired PDBs, mutant PDBs, ``.fxout`` files. Large — gigabytes per
    campaign — and regenerable, so it is never published."""

    pdb_dir: Path
    """Input structures, one ``<pdb>.pdb`` per complex."""

    skempi_csv: Path
    """SKEMPI 2.0 table, semicolon-delimited. Supplies the chain groups and the cleaned
    mutation strings."""

    binary: Optional[Path] = None
    """The FoldX executable. Resolved from ``FOLDX_BIN`` when not given."""

    repair_iterations: int = 1
    """How many times to apply ``RepairPDB``, feeding each result into the next.

    One is FoldX's own documented recommendation and matches CATH-ddG, the only comparator with
    a fully specified protocol. Some groups iterate 5-10 times; the only controlled test of that
    (Usmanova et al. 2018) found no effect on folding ΔΔG or bias. Kept configurable so the
    question can be **measured** here rather than assumed — see
    ``experiments/repair_ablation.py``.

    Changing this invalidates every downstream energy for the complex, so an ablation must
    recompute BuildModel and AnalyseComplex too.
    """

    repair_seed_dirs: Sequence[Path] = field(default_factory=tuple)
    """Other campaigns' work directories to look in for an existing ``<pdb>_Repair.pdb``.

    RepairPDB dominates the runtime, and a complex repaired for one campaign is repaired for
    all of them. Seeding from a sibling campaign turns a full run into BuildModel-only."""

    def __post_init__(self):
        self.results_dir = Path(self.results_dir)
        self.work_dir = Path(self.work_dir)
        self.pdb_dir = Path(self.pdb_dir)
        self.skempi_csv = Path(self.skempi_csv)
        if self.binary is not None:
            self.binary = Path(self.binary)
        self.repair_seed_dirs = tuple(Path(p) for p in self.repair_seed_dirs)
        if self.repair_iterations < 1:
            raise ValueError("repair_iterations must be >= 1")

    # -- binary ------------------------------------------------------------------------------
    def resolve_binary(self) -> Path:
        """Locate the FoldX executable, or explain how to supply it.

        FoldX is licensed software -- free for academic and non-profit institutions, paid for
        commercial use -- and is deliberately not vendored. There is also deliberately no
        fallback to a hardcoded path: a path baked in from one machine makes a missing binary
        surface elsewhere as a confusing FoldX-not-found error instead of as a clear instruction.
        """
        candidate = self.binary or os.environ.get("FOLDX_BIN")
        if not candidate:
            raise RuntimeError(
                "FoldX binary not configured. Set FOLDX_BIN to your FoldX 5 executable, or "
                "pass binary=... . FoldX needs a licence (free for academic/non-profit); "
                "from https://foldxsuite.crg.eu/ . FoldX 5.1 needs no rotabase.txt; earlier "
                "major versions do."
            )
        path = Path(candidate)
        if not path.exists():
            resolved = shutil.which(str(candidate))
            if resolved:
                return Path(resolved)
            raise FileNotFoundError(f"FOLDX_BIN points at {path}, which does not exist.")
        return path

    # -- layout ------------------------------------------------------------------------------
    def ensure_dirs(self) -> None:
        self.results_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir.mkdir(parents=True, exist_ok=True)

    def result_path(self, pdb: str) -> Path:
        return self.results_dir / f"{pdb}.json"

    def complex_work_dir(self, pdb: str) -> Path:
        return self.work_dir / pdb

    def find_repaired(self, pdb: str) -> Optional[Path]:
        """An existing ``<pdb>_Repair.pdb`` from this or any seed campaign."""
        name = f"{pdb}_Repair.pdb"
        for root in (self.work_dir, *self.repair_seed_dirs):
            candidate = Path(root) / pdb / name
            if candidate.exists():
                return candidate
        return None


def default_config(
    root: Optional[Path] = None,
    results_dir: str = "scratch/foldx/results",
    work_dir: str = "scratch/foldx/work",
    **kw,
) -> FoldxConfig:
    """A config rooted at the repository, for the common layout.

    ``root`` defaults to the repository root inferred from this file's location, so the result
    does not depend on the current working directory — resolving ``scratch/...`` against the cwd
    instead makes a script silently write somewhere else when it is run from elsewhere.
    """
    root = Path(root) if root else Path(__file__).resolve().parents[1]
    return FoldxConfig(
        results_dir=root / results_dir,
        work_dir=root / work_dir,
        pdb_dir=Path(kw.pop("pdb_dir", root / "scratch" / "skempi2" / "PDBs")),
        skempi_csv=Path(kw.pop("skempi_csv", root / "scratch" / "skempi_v2.csv")),
        **kw,
    )
