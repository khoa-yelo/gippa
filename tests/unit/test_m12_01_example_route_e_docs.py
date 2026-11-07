"""M12-01 route E (amendment 2026-11-06), work-split row 11a: the docs of route E (section 4,
AC 15 for E, AC 24, "five routes").

Spec docs/specs/M12-01-example-dataset.md section 4 > *Route E, 2026-11-06*: the example page's
``### Route E: DIA-NN report, six-frame database`` section holds exactly route E's three commands,
the ``gippa makedb --six-frame-orfs`` command in prose with a link to ``../user-guide/makedb.md``,
and the route E caveats; "four routes" becomes "five routes" on the page, in the quick-start line
and in the ``gippa example`` group help; the tutorial runs and compares five routes; release step
14 prints the five routes.
"""

from __future__ import annotations

import re

from helpers.example_m12_01 import CLI_MD, REPO, flat, load_routes, option_problems
from helpers.example_m12_01 import cli_reference as _cli_reference

DOCS = REPO / "docs"
DIR = "gippa-example"
OUT = f"{DIR}/runs"
TITLE_E = "Route E: DIA-NN report, six-frame database"
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
FENCE = re.compile(r"^```[^\n]*\n(.*?)^```", re.S | re.M)


def _text(rel: str) -> str:
    path = DOCS / rel
    assert path.is_file(), f"docs/{rel} does not exist"
    return path.read_text(encoding="utf-8")


def _fenced(text: str) -> list[str]:
    out = []
    for block in FENCE.findall(text):
        for line in block.replace("\\\n", " ").splitlines():
            cmd = " ".join(line.split(" #")[0].split())
            if cmd.startswith("gippa "):
                out.append(cmd)
    return out


def _section(text: str, heading: str) -> str:
    level = len(heading) - len(heading.lstrip("#"))
    m = re.search(rf"^{re.escape(heading)}[ \t]*$", text, re.M)
    assert m, f"heading {heading!r} not found"
    nxt = re.compile(rf"^#{{1,{level}}} ", re.M).search(text, m.end())
    return text[m.end() : nxt.start() if nxt else None]


def _page_cmd(cmd: str) -> str:
    return " ".join(cmd.replace("{out}", OUT).replace("{dir}", DIR).split())


# ---------------------------------------------------------------- the example page


def test_route_e_section_holds_exactly_its_commands() -> None:
    """AC 15 for route E."""
    body = _section(_text("get-started/example.md"), f"### {TITLE_E}")
    got = [c for c in _fenced(body) if c.split()[1] in ("convert", "preflight", "run")]
    assert got == [_page_cmd(c) for c in load_routes()["E"]["commands"]]


def test_route_e_section_explains_the_six_frame_database() -> None:
    body = flat(_section(_text("get-started/example.md"), f"### {TITLE_E}"))
    for needle in (
        "gippa makedb --six-frame-orfs",
        "(../user-guide/makedb.md)",
        "`mayu`",
        "floor",
        "1,400,000",
    ):  # fmt: skip  (route E region, 2026-11-06)
        assert needle in body, needle


def test_example_page_says_five_routes() -> None:
    text = flat(_text("get-started/example.md"))
    assert "five routes" in text
    assert "four routes" not in text


def test_example_page_makedb_command_is_valid() -> None:
    """A fenced ``gippa makedb`` line, if any, names only options of its CLI reference section."""
    ref = _cli_reference(CLI_MD.read_text(encoding="utf-8"))
    cmds = [c for c in _fenced(_text("get-started/example.md")) if c.split()[1] == "makedb"]
    assert [p for c in cmds for p in option_problems(c, ref)] == []


def test_quickstart_says_five_input_routes() -> None:
    text = flat(_text("get-started/quickstart.md"))
    assert "five input routes" in text and "four input routes" not in text


def test_cli_group_help_says_five_input_routes() -> None:
    from click.testing import CliRunner

    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["example", "--help"])
    assert res.exit_code == 0, res.output
    assert "five input routes" in flat(res.output)
    text = flat(CLI_MD.read_text(encoding="utf-8"))
    assert "five input routes" in text and "four input routes" not in text


def test_release_step_14_prints_the_five_routes() -> None:
    text = flat((DOCS / "release.md").read_text(encoding="utf-8"))
    assert "prints the five routes" in text and "prints the four routes" not in text


# ---------------------------------------------------------------- the tutorial (AC 24)


def test_tutorial_title_five_routes() -> None:
    lines = _text("tutorials/example-routes.md").splitlines()
    assert lines[0] == "# Tutorial: Nissle regions, five routes"  # route E region, 2026-11-06


def test_tutorial_reads_route_e_novel_orfs() -> None:
    text = _text("tutorials/example-routes.md")
    assert f"{OUT}/E" in text
    assert "novel_orfs.tsv" in text


def test_tutorial_compares_routes_a_to_d() -> None:
    """AC 24 (route E region, 2026-11-06): E runs on another genome, so only A-D are compared."""
    ref = _cli_reference(CLI_MD.read_text(encoding="utf-8"))
    cmds = [c for c in _fenced(_text("tutorials/example-routes.md")) if c.split()[1] == "compare"]
    assert len(cmds) == 1, cmds
    runs = " ".join(f"{OUT}/{r}" for r in "ABCD")
    assert cmds[0].startswith(f"gippa compare {runs} "), cmds[0]
    assert "--labels A,B,C,D " in cmds[0] + " " and "--outdir " in cmds[0], cmds[0]
    assert option_problems(cmds[0], ref) == []


def test_tutorial_holds_no_row_level_content() -> None:
    text = _text("tutorials/example-routes.md")
    assert not PEPTIDE_LIKE.search(text), PEPTIDE_LIKE.search(text)
    assert "GIPPAORF_" not in text


def test_tutorial_runs_route_e_commands() -> None:
    cmds = _fenced(_text("tutorials/example-routes.md"))
    for c in load_routes()["E"]["commands"]:
        assert _page_cmd(c) in cmds, c
