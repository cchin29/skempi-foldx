"""Put `experiments/` on the path so the sweep drivers are importable by name.

They are scripts, not part of the installed package -- a driver is not something a consumer
imports -- but their pure-python helpers (RMSD, list hashing) carry real invariants and are
worth testing.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments" / "repair_sweep"))


import pytest


@pytest.fixture(autouse=True)
def _isolate_exclusion_bypass(monkeypatch):
    """The package documents SKEMPI_FOLDX_ALLOW_INTRACTABLE as a user-settable override, so it can
    be set in the shell that runs the suite. Four tests assert exclusion behaviour and fail if it
    is; that is the environment leaking in, not a real failure."""
    from skempi_foldx.exclusions import ALLOW_ENV

    monkeypatch.delenv(ALLOW_ENV, raising=False)
