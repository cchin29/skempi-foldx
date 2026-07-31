"""Put `experiments/` on the path so the sweep drivers are importable by name.

They are scripts, not part of the installed package -- a driver is not something a consumer
imports -- but their pure-python helpers (RMSD, list hashing) carry real invariants and are
worth testing.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "experiments"))
