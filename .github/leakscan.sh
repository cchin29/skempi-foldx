#!/usr/bin/env bash
# Local paths, hosts and machine names must not reach a published commit. The first one is
# permanent: a hardcoded home directory or a hostname is cheap to catch before a push and
# impossible to remove from history afterwards.
#
# The shipped JSON is in scope, and it is the part most easily forgotten. A campaign records
# where it seeded its repaired structure from, in `meta.repair_seeded_from`, and a scan restricted
# to source files does not look there. Data is the thing this package exists to publish, so it is
# the last place to leave unchecked.
#
# Every text file is scanned. There is no include list: an allowlist of extensions silently drops
# whatever it does not name -- `.gitignore`, `.gitattributes`, a Makefile, a `.cfg` -- and the
# files it drops are exactly the ones nobody thinks to check. `-I` skips binaries, so the figures
# under `docs/assets/` cost nothing.
#
# The invariant this pattern has to hold is that it does not match its own source line. That is
# what lets the scan cover this file and the workflow directory at all: excluding `.github/` to
# stop the pattern matching its own source would leave everything under it unscanned.
# Most alternatives get there by bracketing one of their own characters (`tm[p]`, `[S]ynology`,
# `[s]sh`); the two that quantify a character class are self-immune because the character before
# their `@` here is not in the class they require. That second route is fragile under editing, so
# the check to run after any edit is the whole-pattern one: run this script and expect `clean`.
# The same applies to test fixtures: a home-directory path written literally in a test would fail
# this gate, so `tests/` assembles them with `"/".join([...])`.
#
# Details worth keeping when editing:
#   * The first path segment accepts capitals. A macOS home directory is routinely capitalised,
#     and a lowercase-only class walked straight past one. So does the user@host form: a
#     capitalised login or hostname is ordinary.
#   * The optional backslash before each home-directory slash covers the JSON-escaped spelling.
#     JSON permits it and the shipped store is JSON.
#   * A remote host matches both after a transfer command and bare. The bare form requires `:/`
#     followed by a non-slash, so it does not fire on a PEP 508 direct reference
#     (`name@https://...`), a registry digest (`image@sha256:...`), or an email before a colon.
#   * An IPv4 address matches only after `//` or `@`, not bare, so a four-component version
#     number does not fail the gate.
#   * Two patterns are deliberately absent, both because their false-positive rate on this tree
#     is total: a tilde-username (`~user/...`) cannot be told from Sphinx's `:class:`~module``
#     cross-reference syntax, and a UNC path (two backslashes then a name) cannot be told from
#     the `\\n` that every shipped JSON note contains. A username after a tilde and a Windows
#     share are therefore not caught.
#   * What this cannot catch either: `store.shorten_path()` drops a leading home marker and the
#     name after it, but a username sitting under any other parent -- a mounted volume, say,
#     where the path is `<mount>/<disk>/<user>/...` -- survives into a relative two-segment
#     path with no leading slash, and no pattern here distinguishes that from an ordinary
#     relative path. The store is held to its shape by
#     `test_the_shipped_store_carries_no_path_at_all` instead, which is the check that does
#     generalise.
#   * `--exclude-dir=.git` keeps the object store out of scope now that nothing else does.
#   * `--exclude-dir=.venv*` is the one concession to where this runs. CONTRIBUTING.md tells a
#     contributor to build the environment as `.venv` inside the repo and then run this script,
#     and what it holds matches on two counts: third-party sources carry other people's home
#     directories (pip and pytest each ship one in a comment), and an editable install records
#     *this* machine's checkout path in `direct_url.json` and the `__editable___*_finder.py`
#     shim. Those hits are not this tree's and cannot reach a commit
#     (`.gitignore` carries the same `.venv*/`), but they exit 1 and so fail the release gate
#     the moment anyone follows the documented setup. CI never saw it: there `pip install -e`
#     goes to the runner's own Python, so no venv exists in the checkout to walk.
#
#     Only the venv is excluded, and that is measured rather than assumed. Scanned with the full
#     documented dev state present -- venv, `.pytest_cache`, `skempi_foldx.egg-info`, the
#     `__pycache__` trees, `build/` and `dist/` from the pre-release checks, and the generated
#     `docs/assets` files -- every text hit was under `.venv`.
set -uo pipefail

# grep exits 0 on a match, 1 on none, and 2 on an error -- an unreadable file, a bad pattern. A
# bare `if grep ...` treats 2 like 1, so a scan that never ran would report clean. Check the code.
#
# The root is checked separately because the two greps disagree about it: GNU grep (what CI runs)
# exits 2 on a missing directory, BSD grep (what macOS ships) exits 1, which is indistinguishable
# from a clean scan. A scan of nothing must not pass.
root="${1:-.}"
if [ ! -e "$root" ]; then
  echo "::error::scan root does not exist: $root"
  exit 2
fi

rc=0
grep -rInE '/User[s]\\?/[A-Za-z]|/hom[e]\\?/[A-Za-z]|/roo[t]\\?/[A-Za-z]|/private/tm[p]/|/Volume[s]/|/mn[t]/[A-Za-z]|/ne[t]/[A-Za-z]|/sr[v]/[A-Za-z]|/medi[a]/[A-Za-z]|[~]/[A-Za-z]|[A-Za-z0-9-]+\.loca[l]|[A-Za-z]:\\Use[r]s|[s]mb://|(//|@)[0-9]{1,3}(\.[0-9]{1,3}){3}|[A-Za-z0-9_-]+@[A-Za-z0-9._-]+:/[^/]|([s]sh|[s]cp|[r]sync|[s]ftp)[[:space:]]+[A-Za-z0-9_-]+@|[S]ynologyDrive' \
     --exclude-dir=.git --exclude-dir='.venv*' "$root" || rc=$?

case "$rc" in
  0) echo "::error::local path, address or machine name found above"; exit 1 ;;
  1) echo clean ;;
  *) echo "::error::scan did not run (grep exit $rc)"; exit "$rc" ;;
esac
