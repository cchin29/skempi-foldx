# Changelog

Energy values are tied to a version. A mutation's ΔΔG depends on the whole mutation list it was
computed in, so a rebuilt store is a different measurement rather than a correction of the same
one — values from two versions are not interchangeable within one table. Entries below say
explicitly whether a release changes them.

## 0.1.0 — 2026-08-01

First public release. **Energy values: initial.**

- 6003 records over 344 complexes: 4238 single-point mutations across 322, and 1765 multi-point
  variants across 152. Twelve `AnalyseComplex` terms per record, as mutant minus wild-type.
- `FoldxLookup` resolves the naming conventions a mutation appears under, which is the
  difference between 65.4% and 99.6% coverage on a split whose labels are role-chain.
- Every value computed with a single `RepairPDB` pass, at FoldX 5.1.

Known limitations, each documented where it applies:

- **Three PDB codes carry two SKEMPI interface definitions each** (`2C5D`, `3SE3`, `3SE4`). The
  parser keys on the code, so 31 records are scored against an interface excluding the chain they
  mutate. `interface_suspect()` identifies them; `experiments/recompute_alternate_interfaces.py`
  recomputes them. Fixing the parser moves the coverage denominators, so it is held for 0.2.0
  where the store is rebuilt to match.
- **2494 of 4238 single-point values came from a per-benchmark campaign** rather than from one
  built around each complex's full mutation list. For 27 of those complexes, 287 records, that
  list is a strict subset of the complex's full one, and a value is a property of the list it was
  computed in.
- **`1KBH` is excluded**: `RepairPDB` does not terminate on it under FoldX 5.1.

## Unreleased (0.2.0)

Planned. **Energy values: will change**, including for mutations whose current value is not wrong.

- Rebuild against union mutation lists throughout.
- Compute each SKEMPI interface definition as its own pairing, which also moves the coverage
  denominators to 4343 single-point over 324 interface definitions and 1850 multi-point over
  155, of 323 and 153 distinct PDB codes. Against
  those, and counting only records scored against a pairing SKEMPI associates with them, the
  0.1.0 stores cover 97.3% and 94.4%.
