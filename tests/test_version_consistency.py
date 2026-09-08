"""Release metadata: one version number and one supported range, each declared several times.

Two groups of hand-maintained facts drift independently, and both drifted here before this file
existed.

The version is written in `skempi_foldx.__version__`, `pyproject.toml`, `CITATION.cff`, the newest
heading of `CHANGELOG.md`, and three times in `README.md` -- once in the install command every
reader copies. A release that bumps some and not others ships a citation naming a different
version than the code reports, or hands its readers a command that installs the previous release.

The supported Python range is written in the `requires-python` specifier, the per-minor
classifiers, and the CI matrix. A version advertised but never built is a claim nothing checks; a
version built but never advertised is work whose result nobody can act on.

Most of these are file-to-file comparisons, and those run unconditionally: the files are located
relative to this file, so they are checkable from a checkout or an unpacked sdist whether or not
the installed `skempi_foldx` is this tree. Only the two that compare `__version__` against those
files need the imported package to *be* this tree, and only those skip. That split matters:
guarding all of them goes silent, green, on a plain `pip install .` followed by `pytest`, because
`skempi_foldx` then resolves to site-packages rather than to the tree under test.
"""
import itertools
import re
from pathlib import Path

import pytest

import skempi_foldx

_ROOT = Path(__file__).resolve().parents[1]
_PYPROJECT = _ROOT / "pyproject.toml"
_CITATION = _ROOT / "CITATION.cff"
_CHANGELOG = _ROOT / "CHANGELOG.md"
_README = _ROOT / "README.md"
_CI = _ROOT / ".github" / "workflows" / "tests.yml"

_SEMVER = r"\d+\.\d+\.\d+"


def _toml_load(path: Path) -> dict:
    """Parse TOML with the stdlib where it exists, else `tomli`, else skip.

    `tomllib` is stdlib only from 3.11 and `requires-python` still admits 3.9 and 3.10, so the
    dev extra names `tomli` there. A missing parser is an untestable environment rather than a
    failing invariant -- but an environment where it is missing is one where most of this file
    goes quiet, which is why `test_dev_extra_names_the_parsers_this_file_needs` asserts the extra
    still names it. That test deliberately does not use this function.
    """
    try:
        import tomllib
    except ImportError:  # pragma: no cover -- 3.9/3.10 path
        tomllib = pytest.importorskip("tomli", reason="no TOML parser (tomllib/tomli) available")
    with open(path, "rb") as f:
        return tomllib.load(f)


def _project() -> dict:
    return _toml_load(_PYPROJECT)["project"]


def _cff_scalar(path: Path, key: str) -> str:
    """A top-level scalar from CITATION.cff.

    Read line-wise rather than through a YAML parser: the keys needed here are all top-level
    scalars, and anchoring the match at column 0 cannot pick up a same-named key nested inside
    the file's `references` block, several of which carry their own `title` and `year`. YAML
    cannot place a nested key at column 0, so the anchor is sufficient rather than merely likely.

    Quoted values keep their contents and drop a trailing inline comment; a block scalar
    introducer is rejected rather than returned verbatim, since it would otherwise compare as the
    literal `>-` and fail with a message naming the wrong problem.
    """
    pat = re.compile(rf"^{re.escape(key)}:\s*(.+?)\s*$")
    for line in path.read_text().splitlines():
        m = pat.match(line)
        if not m:
            continue
        raw = m.group(1)
        assert not raw.startswith(("|", ">")), (
            f"{path.name}'s {key!r} is a block scalar; this reader handles plain and quoted ones"
        )
        q = re.match(r"""^(['"])(.*?)\1""", raw)
        if q:
            return q.group(2)
        return raw.split(" #", 1)[0].strip()
    raise AssertionError(f"{path.name} has no top-level {key!r} key")


def _source_version() -> str:
    """`__version__` as written in this tree's `__init__.py`, not as imported.

    An installed copy can shadow the tree -- `pip install .` without `-e`, then `pytest` from the
    checkout, resolves `skempi_foldx` to site-packages. Reading the literal keeps every
    file-to-file comparison a statement about this tree, so a stale installed copy cannot make a
    mismatch look like agreement.
    """
    src = (_ROOT / "skempi_foldx" / "__init__.py").read_text()
    m = re.search(r"""^__version__\s*=\s*["'](.+?)["']""", src, re.M)
    assert m, "skempi_foldx/__init__.py has no __version__ assignment"
    return m.group(1)


def _imported_from_repo() -> bool:
    """Whether the imported package is this source tree rather than an installed copy."""
    try:
        Path(skempi_foldx.__file__).resolve().relative_to(_ROOT)
    except ValueError:
        return False
    return True


_NOT_THIS_TREE = "the imported skempi_foldx is not this source tree; the comparison is meaningless"


def _release_headings(text: str):
    """Every `## x.y.z` release heading, in file order, as `(major, minor, patch)` tuples."""
    found = re.findall(rf"^## ({_SEMVER})\b", text, re.M)
    return [tuple(int(part) for part in v.split(".")) for v in found]


def _changelog_latest_version(path: Path) -> str:
    """The version of the topmost release heading.

    Entries are newest-first, so the top heading is the release being cut and everything below is
    history that must not move. Headings here are `## x.y.z -- date`, without the bracketed link
    form, because this changelog carries no link-reference section.

    The *first* `##` heading is required to be a parseable release heading rather than searched
    past. Skipping over an unparseable one -- `## 0.3.0rc1`, or a bracketed `## [0.3.0]` -- would
    report the release two entries down, so the failure would name the wrong line of the file.
    """
    first = re.search(r"^## (.+)$", path.read_text(), re.M)
    assert first, f"{path.name} has no `## ` heading at all"
    m = re.match(rf"^({_SEMVER})\b", first.group(1))
    assert m, f"{path.name}'s first `## ` heading is not `x.y.z`: {first.group(1)!r}"
    return m.group(1)


def _declared_python_minors(project: dict) -> set:
    """The `3.x` minors claimed by the `Programming Language :: Python :: 3.x` classifiers.

    `:: 3 :: Only` and `:: 3.14 :: Only` are trailing-qualifier classifiers, not version claims,
    and are excluded rather than crashed on -- `int("14 :: Only")` raises from inside three
    separate tests and names none of them usefully.
    """
    out = set()
    for c in project["classifiers"]:
        m = re.fullmatch(r"Programming Language :: Python :: 3\.(\d+)", c)
        if m:
            out.add(int(m.group(1)))
    return out


def _requires_python_bounds(spec: str):
    """`(floor_minor, ceiling_minor_or_None)` parsed from a `requires-python` string.

    Written to survive the ceiling being absent, which is the intended state. An upper bound in
    `requires-python` is a resolver gate rather than documentation: above it pip does not report
    that the interpreter is too new, it backtracks through this project's own tags. Both spellings
    of a ceiling are read -- `<3.14` bounds the advertised range exclusively and `<=3.13`
    inclusively -- and both are normalised to the exclusive form returned here.
    """
    lo = re.search(r">=\s*3\.(\d+)", spec)
    assert lo, f"no parseable >=3.x floor in requires-python: {spec!r}"
    hi_incl = re.search(r"<=\s*3\.(\d+)", spec)
    hi_excl = re.search(r"(?<!<)<\s*3\.(\d+)", spec)
    if hi_incl:
        return int(lo.group(1)), int(hi_incl.group(1)) + 1
    if hi_excl:
        return int(lo.group(1)), int(hi_excl.group(1))
    return int(lo.group(1)), None


def _requires_python_exclusions(spec: str) -> set:
    """The `3.x` minors excluded by a `!=3.x.*` clause, which the bounds above cannot express."""
    return {int(m) for m in re.findall(r"!=\s*3\.(\d+)\.\*", spec)}


# --- the version, in seven places -------------------------------------------------------------

def _optional_dependency_group(src: str, name: str):
    """The raw body of `name = [...]` under `[project.optional-dependencies]`, or None.

    Scanned as text rather than parsed, because the caller is checking that the TOML parser is
    declared -- using the parser to do it would skip in exactly the case worth catching.

    Table-aware on purpose. A same-named list under any other table is not what
    `pip install .[extra]` resolves, and neither is a PEP 735 `[dependency-groups]` entry, so
    finding one of those and reporting success would be the silent pass this test exists to
    prevent. Bracket depth is counted outside quotes, so a requirement carrying its own brackets
    (`"coverage[toml]"`) does not end the list early.
    """
    lines = src.splitlines()
    try:
        start = next(i for i, l in enumerate(lines)
                     if l.strip() == "[project.optional-dependencies]")
    except StopIteration:
        return None

    body, depth, collecting = [], 0, False
    for line in lines[start + 1:]:
        if not collecting:
            if line.lstrip().startswith("[") and line.rstrip().endswith("]") and "=" not in line:
                return None  # next table, without having found the group
            if not re.match(rf"^\s*{re.escape(name)}\s*=\s*\[", line):
                continue
            collecting, line = True, line.split("[", 1)[1]

        quote = None
        for ch in line:
            if quote:
                if ch == quote:
                    quote = None
                continue
            if ch in "\"'":
                quote = ch
            elif ch == "[":
                depth += 1
            elif ch == "]":
                if depth == 0:
                    return "".join(body) + line
                depth -= 1
        body.append(line)
    return None


def _matrix_python_versions(matrix: dict) -> set:
    """The Python versions a GitHub Actions matrix actually builds.

    Modelled in the order Actions applies, which is not the order the keys appear in: the base
    matrix is expanded to its full product, `exclude:` removes combinations matching *every* key
    it names, and `include:` adds afterwards -- so an `include`-added combination cannot be
    excluded, and an `exclude` naming the only value of a key wipes the base matrix entirely.

    Getting that order wrong is a silent pass in both directions: an `exclude` that names only a
    python-version drops it from every OS, and an `exclude` naming a value no combination has
    removes nothing.
    """
    base = {k: v for k, v in matrix.items() if k not in ("include", "exclude")}
    keys = sorted(base)
    combos = [dict(zip(keys, vals)) for vals in itertools.product(*(base[k] for k in keys))]

    for ex in matrix.get("exclude", []):
        combos = [c for c in combos
                  if not all(str(c.get(k)) == str(v) for k, v in ex.items())]

    built = {str(c["python-version"]) for c in combos if "python-version" in c}
    built |= {str(e["python-version"]) for e in matrix.get("include", []) if "python-version" in e}
    return built


def test_version_is_a_release_number():
    # Read from the source file: an installed copy can shadow the tree, and a pre-release left in
    # `__init__.py` would otherwise pass by being checked against something else entirely.
    assert re.fullmatch(_SEMVER, _source_version()), _source_version()


def test_pyproject_version_matches_package():
    if not _imported_from_repo():
        pytest.skip(_NOT_THIS_TREE)
    assert _project()["version"] == skempi_foldx.__version__


def test_all_declarations_agree():
    # The pairwise tests localize a mismatch; this one states the invariant itself, so a further
    # location can join the mapping without inventing a new pairing for it.
    if not _imported_from_repo():
        pytest.skip(_NOT_THIS_TREE)
    declared = {
        "skempi_foldx/__init__.py": _source_version(),
        "skempi_foldx.__version__ (imported)": skempi_foldx.__version__,
        "pyproject.toml": _project()["version"],
        "CITATION.cff": _cff_scalar(_CITATION, "version"),
        "CHANGELOG.md (newest heading)": _changelog_latest_version(_CHANGELOG),
    }
    assert len(set(declared.values())) == 1, declared


def test_citation_version_matches_pyproject():
    assert _cff_scalar(_CITATION, "version") == _project()["version"]


def test_changelog_documents_the_packaged_version():
    """The version being released has an entry, and it is the newest one.

    A bumped version with no entry, or an entry that is not at the top, is invisible until a
    reader goes looking after the fact -- which for a package installed by git tag is after they
    have already resolved it.
    """
    assert _changelog_latest_version(_CHANGELOG) == _project()["version"]


def test_readme_install_pin_matches_pyproject():
    """The version in the README's own install command is the version being released.

    This is the declaration with the shortest path to a reader. A release that bumps everything
    else and leaves this behind does not merely misdocument -- it hands every reader a command
    that installs the previous release.

    The package is not distributed on PyPI, so the install command is a git-tag URL and this
    checks that form. One git pin in the README is deliberately frozen -- the example showing how
    to install 0.1.0 -- so this asserts that the *current* version is pinned, rather than that
    every pin names it. `test_readme_git_pins_are_historical_only` covers the rest.

    Two occurrences, not one: the install command under Quickstart, and the pin quoted under
    *Scope and stability*. Anchoring on the command alone left the second unchecked, and it could
    be bumped, or not bumped, with nothing failing.
    """
    text = _README.read_text()
    version = _project()["version"]

    # The FIRST install command in the file is the Quickstart, and it is the one a reader runs.
    # Counting pins is not enough: the historical 0.1.0 example and the withdrawn-tag recovery
    # command are also `pip install ...@vX.Y.Z`, so a Quickstart left on the previous release
    # still leaves plenty of correct pins elsewhere in the file to satisfy a count.
    installs = re.findall(rf"pip install[^\n]*?skempi-foldx(?:\.git)?@v({_SEMVER})", text)
    assert installs, "README.md has no `pip install ...@vX.Y.Z` command"
    assert installs[0] == version, {
        "first install command in README.md names": installs[0],
        "pyproject.toml version": version,
        "note": "the Quickstart install must name the release being shipped",
    }

    pins = re.findall(rf"skempi-foldx(?:\.git)?@v({_SEMVER})", text)
    assert pins.count(version) >= 2, {
        "git-tag pins in README.md": sorted(set(pins)),
        "pyproject.toml version": version,
        "note": "expected the Quickstart install and the Scope-and-stability pin to name it",
    }

    # The prose version line is a declaration too, and the one a reader is most likely to quote.
    stated = set(re.findall(rf"^Version ({_SEMVER})\.", text, re.M))
    assert stated == {version}, {
        "stated in README.md prose": sorted(stated),
        "pyproject.toml version": version,
    }

    # The package is not on PyPI. A `skempi-foldx==X.Y.Z` pin would be a command that fails for
    # every reader, which is exactly the state this release corrected.
    assert not re.search(r"pip install\s+skempi-foldx==", text), (
        "README.md shows a PyPI install; this package is not distributed on PyPI"
    )


def test_readme_git_pins_are_historical_only():
    """Every git-tag URL in the README names a release the changelog records.

    Those URLs are the install path, so most of them track the current version. The exception is
    the example showing how to install the previous release, which is frozen by design. This pins
    the intent: a git pin may name any version this changelog records, and nothing else -- so a
    typo, or a bump to a version that was never released, fails.

    Any *released* version, deliberately, rather than `{previous, current}`: an allowance tied to
    the current version would stop permitting the historical example the moment the version moves
    past it, failing the next release on a string the release must not touch.
    """
    text = _README.read_text()
    git_pins = set(re.findall(rf"skempi-foldx(?:\.git)?@v({_SEMVER})", text))
    released = set(re.findall(rf"^## ({_SEMVER})", _CHANGELOG.read_text(), re.M))
    assert git_pins <= released, {
        "git-tag pins in README.md": sorted(git_pins),
        "versions CHANGELOG.md records": sorted(released),
        "note": "a git pin here must name a release this changelog documents",
    }
    # `<= released` alone would let the 0.1.0 example be bumped to the current version, which is
    # in the allowed set for a different reason. 0.1.0 is a historical snapshot a reader may still
    # need to name, so the example that shows how to install it has to survive every bump.
    assert "0.1.0" in git_pins, (
        "README.md no longer shows how to install 0.1.0"
    )


def test_changelog_entries_are_newest_first():
    """Release headings strictly descend.

    Keep a Changelog order, and the reason it matters here is specific: energy values are tied to
    a release, so a reader who takes the first entry as current takes the wrong store's counts.
    Strict rather than non-decreasing, so a bad merge that duplicates a section is caught too.
    """
    versions = _release_headings(_CHANGELOG.read_text())
    assert versions, "no `## x.y.z` release headings"
    assert all(a > b for a, b in zip(versions, versions[1:])), versions


def test_citation_release_date_is_iso():
    # Not pinned to a literal date, which would rot every release. The failure this catches is a
    # date left in a shape no citation manager can read.
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", _cff_scalar(_CITATION, "date-released"))


# --- the supported Python range, in three places --------------------------------------------

def test_python_classifiers_agree_with_requires_python():
    """The classifiers state the range that is built; `requires-python` is only the gate.

    Deliberately not an equality against a hardcoded range: with no ceiling there is no upper
    bound to compare against, and hardcoding one would have to be edited for every future CPython
    -- the maintenance this project avoids by not capping. What must hold is that the two never
    contradict: the floor is advertised and nothing below it is, and nothing sits above a ceiling
    if one is ever added.
    """
    project = _project()
    floor, ceiling = _requires_python_bounds(project["requires-python"])
    declared = _declared_python_minors(project)

    assert declared, "no `Programming Language :: Python :: 3.x` classifiers"
    assert min(declared) == floor, {
        "requires-python floor": f"3.{floor}",
        "lowest classifier": f"3.{min(declared)}",
        "why": "the floor must be advertised, and nothing below it may be",
    }
    if ceiling is not None:
        assert max(declared) < ceiling, {
            "requires-python ceiling (exclusive)": f"<3.{ceiling}",
            "highest classifier": f"3.{max(declared)}",
        }


def test_requires_python_excludes_nothing_it_advertises():
    """A `!=3.x.*` clause is invisible to the floor/ceiling pair, and silently contradicts.

    The resolver would refuse that interpreter while the classifiers advertise it and CI builds
    it -- three files disagreeing with no bound out of place to notice.
    """
    project = _project()
    excluded = _requires_python_exclusions(project["requires-python"])
    assert not (excluded & _declared_python_minors(project)), {
        "excluded by requires-python": sorted(f"3.{m}" for m in excluded),
        "but advertised by a classifier": sorted(
            f"3.{m}" for m in excluded & _declared_python_minors(project)
        ),
    }


def test_python_classifiers_have_no_gaps():
    """The claimed minors are contiguous.

    Catches the drift this guards in practice: adding a new CPython to the classifiers or the CI
    matrix and missing the one before it. A genuine hole in support is not expressible here, but
    for a package with no dependencies it is not a reachable state either.
    """
    declared = _declared_python_minors(_project())
    assert declared, "no `Programming Language :: Python :: 3.x` classifiers"
    expected = set(range(min(declared), max(declared) + 1))
    assert declared == expected, {"missing": sorted(expected - declared)}


def test_ci_matrix_matches_classifiers():
    """Every advertised Python is built, and nothing is built that is not advertised.

    Reads the matrix rather than restating it, so adding a version means editing one list and not
    two. The gap this closes: the metadata said `Programming Language :: Python :: 3` while CI ran
    3.9 and 3.12, so four of the six versions now claimed had never been executed.

    `.github/` is deliberately absent from the sdist, so this is the one check here that a
    self-verifying sdist cannot run; it skips there rather than failing.
    """
    if not _CI.is_file():
        pytest.skip("no CI workflow in this tree (expected in an sdist)")
    yaml = pytest.importorskip("yaml")

    matrix = yaml.safe_load(_CI.read_text())["jobs"]["test"]["strategy"]["matrix"]
    built = _matrix_python_versions(matrix)

    declared = {f"3.{m}" for m in _declared_python_minors(_project())}
    assert built == declared, {
        "advertised but never built": sorted(declared - built),
        "built but not advertised": sorted(built - declared),
    }


def test_dev_extra_names_the_parsers_this_file_needs():
    """`tomli` and `pyyaml` stay in the dev extra, so nothing here skips on a supported floor.

    Without them most of this file is import-guarded into silence -- a green bar reporting that
    nothing was checked. An extra that loses `tomli` produces exactly that on a 3.9 or 3.10 leg,
    which are also the two interpreters least likely to be run locally.

    Read as raw text on purpose. Using the TOML parser to check that the TOML parser is declared
    is circular -- on 3.9 without `tomli` this test would itself skip, which is the case it
    exists to catch.
    """
    body = _optional_dependency_group(_PYPROJECT.read_text(), "dev")
    assert body is not None, (
        "pyproject.toml has no `dev` list under [project.optional-dependencies]. A `dev` group "
        "declared anywhere else -- another table, or PEP 735 [dependency-groups] -- is not what "
        "`pip install .[dev]` resolves, which is what CI runs."
    )
    missing = [n for n in ("tomli", "pyyaml") if not re.search(rf"[\"']{n}\b", body)]
    assert not missing, {
        "missing from the dev extra": missing,
        "consequence": "release-metadata tests skip silently instead of running",
    }
