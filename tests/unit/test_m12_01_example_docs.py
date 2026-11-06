"""M12-01 PR 8a: the example page, the example tutorial and their links (section 4 "Docs", AC 15,
Q-K4, Q-L5). Spec: docs/specs/M12-01-example-dataset.md.

Already covered elsewhere, not repeated here:

- AC 16 (``## gippa example fetch`` in ``docs/reference/cli.md``, ``gen_reference.py --check``):
  ``tests/unit/test_m12_01_example_cli_reference.py``.
- The exact Get started and Tutorials nav lists: ``tests/unit/test_m11_07_mkdocs_config.py``
  (``get-started/example.md`` right after the quick start is already there; the Tutorials list is
  edited by the lead for the rename, Q-L5).
- The ``PEPTIDE_LIKE`` guard on every ``docs/tutorials/*.md``: ``tests/unit/test_m11_08_docs.py``.

AC 15: the fenced route commands of ``docs/get-started/example.md`` equal
``src/gippa/data/example_routes.yaml``'s, with ``{dir}`` written as ``gippa-example`` (the default
``--outdir``) and so ``{out}`` (``{dir}/runs``) as ``gippa-example/runs``.

Q-K4 (default yes): ``docs/tutorials/example-routes.md``, "Tutorial: one Nissle region, four
routes", replaces the empty ``docs/tutorials/ecoli-k12.md``; the K-12 tutorial goes to
``docs/backlog.md``.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml
from helpers.example_m12_01 import CLI_MD, ROUTES_YAML, cli_reference, load_routes, option_problems

REPO = Path(__file__).resolve().parents[2]
DOCS = REPO / "docs"
EXAMPLE = "get-started/example.md"
TUTORIAL = "tutorials/example-routes.md"
OLD_TUTORIAL = "tutorials/ecoli-k12.md"
TUTORIAL_TITLE = "# Tutorial: one Nissle region, four routes"
DIR = "gippa-example"
OUT = f"{DIR}/runs"
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # the M10-05 / M11-08 guard regex
FENCE = re.compile(r"^```[^\n]*\n(.*?)^```", re.S | re.M)
ROUTE_SUBS = ("convert", "preflight", "run")
TIERS = ("mayu", "protein_list", "two_peptide")


def _text(rel: str) -> str:
    path = DOCS / rel
    assert path.is_file(), f"docs/{rel} does not exist"
    return path.read_text(encoding="utf-8")


def _norm(cmd: str) -> str:
    return " ".join(cmd.split())


def _fenced_commands(text: str) -> list[str]:
    """``gippa ...`` lines of the fenced blocks: continuations joined, comments dropped."""
    out = []
    for block in FENCE.findall(text):
        for line in block.replace("\\\n", " ").splitlines():
            cmd = _norm(line.split(" #")[0])
            if cmd.startswith("gippa "):
                out.append(cmd)
    return out


def _route_commands(cmds: list[str]) -> list[str]:
    return [c for c in cmds if c.split()[1] in ROUTE_SUBS]


def _expected(cmd: str) -> str:
    return _norm(cmd.replace("{out}", OUT).replace("{dir}", DIR))


def _section(text: str, heading: str) -> str:
    """The body under ``heading`` (a full heading line) up to the next heading of its level or
    higher."""
    level = len(heading) - len(heading.lstrip("#"))
    m = re.search(rf"^{re.escape(heading)}[ \t]*$", text, re.M)
    assert m, f"heading {heading!r} not found"
    nxt = re.compile(rf"^#{{1,{level}}} ", re.M).search(text, m.end())
    return text[m.end() : nxt.start() if nxt else None]


def _nav() -> list[str]:
    class _Loader(yaml.SafeLoader):
        pass

    _Loader.add_multi_constructor("!", lambda loader, suffix, node: None)
    cfg = yaml.load((REPO / "mkdocs.yml").read_text(encoding="utf-8"), Loader=_Loader)

    def flat(node: object) -> list[str]:
        if isinstance(node, str):
            return [node]
        items = node.values() if isinstance(node, dict) else node
        return [p for item in items for p in flat(item)]  # type: ignore[union-attr]

    return flat(cfg["nav"])


@pytest.fixture(scope="module")
def routes() -> dict[str, dict[str, object]]:
    assert ROUTES_YAML.is_file(), "src/gippa/data/example_routes.yaml does not exist"
    return load_routes()


@pytest.fixture(scope="module")
def ref() -> dict[tuple[str, ...], set[str]]:
    return cli_reference(CLI_MD.read_text(encoding="utf-8"))


# ---------------------------------------------------------------- AC 15: example.md


def test_example_page_route_commands_equal_the_yaml(routes: dict[str, dict[str, object]]) -> None:
    """AC 15: every fenced convert/preflight/run line of the page, in order, is a route command of
    ``example_routes.yaml`` with ``{dir}`` -> ``gippa-example`` and ``{out}`` ->
    ``gippa-example/runs``; no command is missing or extra."""
    want = [_expected(c) for r in routes.values() for c in r["commands"]]  # type: ignore[attr-defined]
    assert _route_commands(_fenced_commands(_text(EXAMPLE))) == want


@pytest.mark.parametrize("rid", ["A", "B", "C", "D"])
def test_example_page_route_section_holds_its_commands(
    routes: dict[str, dict[str, object]], rid: str
) -> None:
    """AC 15, per route: the section headed by the route's YAML title holds exactly its commands."""
    route = routes[rid]
    body = _section(_text(EXAMPLE), f"### {route['title']}")
    want = [_expected(c) for c in route["commands"]]  # type: ignore[attr-defined]
    assert _route_commands(_fenced_commands(body)) == want


def test_example_page_has_no_placeholder_left() -> None:
    """The fenced commands are runnable after a default fetch: no ``{dir}``/``{out}`` and no
    ``<dir>``/``<out>`` placeholder in a route command."""
    for cmd in _route_commands(_fenced_commands(_text(EXAMPLE))):
        assert not re.search(r"\{(dir|out)\}|<(dir|out)>", cmd), cmd


def test_example_page_has_no_row_level_data() -> None:
    """AC 15: the M11-08 guard (``PEPTIDE_LIKE``) covers this page too."""
    text = _text(EXAMPLE)
    assert not PEPTIDE_LIKE.search(text), PEPTIDE_LIKE.search(text)


def test_example_page_options_exist_in_the_cli_reference(
    ref: dict[tuple[str, ...], set[str]],
) -> None:
    """Option drift on every fenced ``gippa`` line of the page (multi-word sections, Q-L4)."""
    problems = [p for c in _fenced_commands(_text(EXAMPLE)) for p in option_problems(c, ref)]
    assert problems == []


def test_example_page_content() -> None:
    """Section 4: fetch and ``--from-file``, the tiers per route, route C as a format demo, the
    comparison with ``expected/README.md``, and the pointer to From raw spectra."""
    text = _text(EXAMPLE)
    for needle in ("gippa example fetch", "--from-file", "PXD078459", "gippa setup --no-genomad",
                   "`mayu`", "format demo", "expected/README.md", "(from-spectra.md)"):  # fmt: skip
        assert needle in text, needle


# ---------------------------------------------------------------- the tutorial (Q-K4)


def test_k12_stub_is_replaced() -> None:
    """Q-K4: ``tutorials/ecoli-k12.md`` is gone from disk and nav; the example tutorial is in."""
    nav = _nav()
    assert not (DOCS / OLD_TUTORIAL).exists()
    assert OLD_TUTORIAL not in nav
    assert TUTORIAL in nav
    assert "tutorials/nissle-pilot.md" in nav  # stays a stub (out of scope)


def test_tutorial_title_and_not_a_stub() -> None:
    text = _text(TUTORIAL)
    assert text.splitlines()[0] == TUTORIAL_TITLE
    assert "Stub" not in text and "Filled by" not in text


def test_tutorial_reads_genes_and_report_per_route() -> None:
    """Run the four routes, read ``genes.tsv`` and the report; the tier changes with the input."""
    text = _text(TUTORIAL)
    for needle in ("genes.tsv", "report.html", "../get-started/example.md"):
        assert needle in text, needle
    for tier in TIERS:
        assert f"`{tier}`" in text, tier
    for rid in "ABCD":
        assert f"{OUT}/{rid}" in text, rid


def test_tutorial_compares_the_four_runs(ref: dict[tuple[str, ...], set[str]]) -> None:
    """``gippa compare`` of the four run directories in order A-D with ``--labels A,B,C,D``.

    ``--outdir`` is required by ``gippa compare`` (docs/reference/cli.md), so the fenced command
    carries it although the spec's one-liner omits it."""
    cmds = [c for c in _fenced_commands(_text(TUTORIAL)) if c.split()[1] == "compare"]
    assert len(cmds) == 1, cmds
    cmd = cmds[0]
    runs = " ".join(f"{OUT}/{rid}" for rid in "ABCD")
    assert cmd.startswith(f"gippa compare {runs} "), cmd
    assert "--labels A,B,C,D" in cmd and "--outdir " in cmd, cmd
    assert option_problems(cmd, ref) == []


def test_tutorial_commands_are_route_commands(
    routes: dict[str, dict[str, object]], ref: dict[tuple[str, ...], set[str]]
) -> None:
    """Every fenced convert/preflight/run line of the tutorial is one of the YAML route commands
    (written as on the example page), and every fenced line passes option drift."""
    allowed = {_expected(c) for r in routes.values() for c in r["commands"]}  # type: ignore[attr-defined]
    cmds = _fenced_commands(_text(TUTORIAL))
    assert [c for c in _route_commands(cmds) if c not in allowed] == []
    assert [p for c in cmds for p in option_problems(c, ref)] == []


def test_tutorial_holds_no_row_level_data() -> None:
    """Aggregates only (guarded for the new file even before the lead's TUTORIALS edit)."""
    assert not PEPTIDE_LIKE.search(_text(TUTORIAL))


def test_k12_tutorial_is_in_the_backlog() -> None:
    """Q-K4: the K-12 tutorial moves to ``docs/backlog.md``."""
    lines = _text("backlog.md").splitlines()
    assert any(re.search(r"K-12", ln) and re.search(r"[Tt]utorial", ln) for ln in lines)


# ---------------------------------------------------------------- links from other pages


def test_quickstart_points_to_the_example() -> None:
    """Section 1 of the quick start names ``gippa example fetch`` after ``gippa demo`` and links
    the page."""
    sec = _section(_text("get-started/quickstart.md"), "## 1. Install and check")
    assert "gippa demo" in sec and "gippa example fetch" in sec
    assert sec.index("gippa example fetch") > sec.index("gippa demo")
    assert "(example.md)" in sec


@pytest.mark.parametrize(
    ("rel", "link"),
    [("get-started/from-spectra.md", "(example.md)"),
     ("user-guide/convert.md", "(../get-started/example.md)")],
)  # fmt: skip
def test_page_links_the_example_for_routes_a_c_d(rel: str, link: str) -> None:
    """One sentence linking the example as a way to try routes A, C and D on real files."""
    text = _text(rel)
    paras = [p for p in re.split(r"\n\s*\n", text) if link in p]
    assert paras, f"{rel} does not link {link}"
    assert any(re.search(r"\bA\b.*\bC\b.*\bD\b", " ".join(p.split())) for p in paras), paras
