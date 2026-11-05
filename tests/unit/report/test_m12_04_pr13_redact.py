"""M12-04 PR 13, R-09 and AC 13.2: ``gippa report --redact`` makes a shareable report
(docs/specs/M12-04-outputs-revision.md § PR 13).

``--redact`` rewrites the provenance command line: every absolute path becomes
``<input>/<basename>``, the output directory ``<outdir>``, the ``--profile`` value ``<profile>``.
The default report is unchanged (absolute paths shown as file names, the outdir as ``OUT``,
M9-02 behaviour 6; tests/unit/report/test_command_line_redaction.py). A manifest without
``command_line`` has nothing to rewrite.

Interface pinned here (PR 13 test author, beyond the spec text):
- ``gippa.report.content.redact_command_line(args: Sequence[str], outdir: str) -> str``: pure;
  the manifest's ``command_line`` joined by single spaces, redacted. The path forms are those of
  the default rule (a whole token, after the last ``=``, an item of a comma list, attached to a
  short option ``-o/x``) and, review r1 M1, every ``/`` that starts a path inside a ``--set``
  YAML value: after ``=`` (the first one: ``--set=key=/home/x=y/z`` is the path
  ``/home/x=y/z``), ``[``, ``,``, ``{``, a space, a quote or ``:`` (``k=x:/a/b``); the path runs
  to the next ``,``, ``]``, ``}``, quote, whitespace or the end, and the brackets, braces and
  quotes around it are kept. A path under ``outdir`` becomes ``<outdir>/<relative>``
  (``<outdir>`` for the outdir itself), any other ``<input>/<last component>``; the profile
  value follows ``--profile`` as the next token or after ``--profile=``; relative arguments are
  unchanged. The default (non-redact) rule is unchanged;
- ``gippa.report.build.build_report(..., redact=True)`` shows that string in the provenance
  ``command line`` row (``<table data-table-provenance="run">``); ``redact`` defaults to False
  and the default page is byte-identical to a build without the keyword;
- ``gippa report OUT --redact`` writes ``OUT/report/report.html`` redacted;
- review r1: under ``redact`` no ``<p class="note">`` of the page (an unreadable table, e.g.
  ``detection_summary.tsv: could not be read: …``) holds an absolute path; it still starts with
  ``<file>: could not be read:``. Without ``redact`` the note is unchanged.

AC 13.2 reading: "a path component" is a directory component of an absolute input path (the
basename is kept on purpose, after ``<input>/``).
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from pathlib import Path, PurePosixPath
from types import ModuleType
from typing import Any

import pytest
from click.testing import CliRunner, Result
from helpers.report_content import manifest, write_manifest
from helpers.report_fixtures import Page, cfg, parse, run_out, swap_registry
from hypothesis import given, settings
from hypothesis import strategies as st

OUTDIR = "/zzscratch/zzlabgroup/zzruns/A"
GENOME = "/zzscratch/zzlabgroup/zzproject/zzgenome/nissle_r200k.fna.gz"
PEPTIDES = "/zzhome/zzuser/zzdata/report.pr_matrix.tsv"
PROFILE = "zzsitecluster,apptainer"
PROFILE_NAMES = ("zzsitecluster", "apptainer")

#: an absolute path in a rendered command line: one starting a token, following ``=``, ``,``,
#: ``[``, ``{``, a quote or ``:``, or attached to a short option (``-o/x``); ``<input>/x`` and
#: ``<outdir>/x`` are not
ABS = re.compile(r"(?:^|[\s=,\[{'\":])/|(?:^|\s)-[A-Za-z]/")


def redact(args: Sequence[str], outdir: str = OUTDIR) -> str:
    from gippa.report.content import redact_command_line

    return redact_command_line(["gippa", "run", *args], outdir)


def dir_components(*paths: str) -> set[str]:
    return {part for p in paths for part in PurePosixPath(p).parent.parts if part != "/"}


# ---------------------------------------------------------------- the pure rule


def test_spec_example() -> None:
    """Before: ``gippa run --profile <site profile>,apptainer --genome /…/genome/nissle….fna.gz``;
    after: ``gippa run --profile <profile> --genome <input>/nissle_r200k.fna.gz``."""
    got = redact(["--profile", PROFILE, "--genome", GENOME, "--outdir", OUTDIR])
    assert got == ("gippa run --profile <profile> --genome <input>/nissle_r200k.fna.gz "
                   "--outdir <outdir>")  # fmt: skip


@pytest.mark.parametrize(
    ("args", "shown"),
    [
        (["--profile=" + PROFILE], "--profile=<profile>"),
        (["--profile", "apptainer"], "--profile <profile>"),
        (["--outdir=" + OUTDIR], "--outdir=<outdir>"),
        (["--outdir", OUTDIR + "/"], "--outdir <outdir>"),
        (["--set", f"k={OUTDIR}/tables/x.tsv"], "--set k=<outdir>/tables/x.tsv"),
        (["--peptides", PEPTIDES], "--peptides <input>/report.pr_matrix.tsv"),
        (["--genome=" + GENOME], "--genome=<input>/nissle_r200k.fna.gz"),
        (["--set", "db=/zzhome/zzuser/db/"], "--set db=<input>/db"),
        (["--set=k=/zzhome/zzuser/c.yaml"], "--set=k=<input>/c.yaml"),
        (["--reads", "a.fq,/zzhome/zzuser/b.fq"], "--reads a.fq,<input>/b.fq"),
        (["-o/zzhome/zzuser/x"], "-o<input>/x"),
        # a sibling with the same string prefix is not under the outdir
        (["--genome", OUTDIR + "x/y.fa"], "--genome <input>/y.fa"),
        # review r1 M1: paths inside ``--set`` YAML values
        (["--set", "k=[/a/b,/c/d]"], "--set k=[<input>/b,<input>/d]"),
        (["--set", "k='/a/b'"], "--set k='<input>/b'"),
        (["--set", 'k="/a/b"'], '--set k="<input>/b"'),
        (["--set", "k={d: /a/b}"], "--set k={d: <input>/b}"),
        (["--set", "k=x:/a/b"], "--set k=x:<input>/b"),
        (["--set=key=/home/x=y/z"], "--set=key=<input>/z"),
        (["--set", f"k=[{OUTDIR}/t/x.tsv,/a/b]"], "--set k=[<outdir>/t/x.tsv,<input>/b]"),
        (["--set", f"k='{OUTDIR}'"], "--set k='<outdir>'"),
        (["--set", f"k={{d: {OUTDIR}/t}}"], "--set k={d: <outdir>/t}"),
        ([f"--set=k={OUTDIR}/x=y/z"], "--set=k=<outdir>/x=y/z"),
    ],
)
def test_redacted_forms(args: list[str], shown: str) -> None:
    got = redact(args)
    assert got == f"gippa run {shown}"
    assert not ABS.search(got), got


@pytest.mark.parametrize(
    "args",
    [
        ["--genome", "genome.fa", "--peptides", "peptides.tsv", "--seed", "7"],
        ["--genome", "data/g.fa"],
        ["--set", "modules=[core,function]"],
        ["--reads", "a.fq,b.fq"],
        ["--outdir", "OUT"],
        ["-n2"],
    ],
)
def test_relative_arguments_are_unchanged(args: list[str]) -> None:
    assert redact(args) == "gippa run " + " ".join(args)


def test_no_command_line_is_empty() -> None:
    from gippa.report.content import redact_command_line

    assert redact_command_line([], OUTDIR) == ""


# ---------------------------------------------------------------- the report page


def _out(tmp: Path, command_line: list[str] | None, *, drop: bool = False) -> Path:
    out = run_out(tmp)
    over: dict[str, Any] = {"outdir": OUTDIR, "profile": PROFILE}
    if command_line is not None:
        over["command_line"] = command_line
    doc = manifest(**over)
    doc["inputs"] = [{**i, "path": GENOME if i.get("role") == "genome" else PEPTIDES}
                     for i in doc.get("inputs", [])]  # fmt: skip
    if drop:
        doc.pop("command_line", None)
    write_manifest(out, doc)
    return out


LINE = ["gippa", "run", "--profile", PROFILE, "--genome", GENOME, "--peptides", PEPTIDES,
        "--outdir", OUTDIR, "--set", f"db={PEPTIDES}"]  # fmt: skip


def _page(out: Path, **kw: Any) -> str:
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    return build_report(RunView.from_outdir(out), {}, {}, cfg(), **kw)


def command_row(page: Page) -> str:
    runs = [t for t in page.root.by_tag("table")
            if t.attrs.get("data-table-provenance") == "run"]  # fmt: skip
    assert len(runs) == 1, "no provenance run table"
    for tr in runs[0].by_tag("tr"):
        th = [c for c in tr.by_tag("th")]
        if th and " ".join(th[0].text().split()) == "command line":
            return " ".join(tr.by_tag("td")[0].text().split())
    raise AssertionError("no command line row")


def ac_13_2_violations(html: str) -> list[str]:
    """AC 13.2: ``/`` followed by a directory component of the fixture's absolute input paths,
    or a profile name from the manifest."""
    bad = [f"/{c}" for c in dir_components(GENOME, PEPTIDES, OUTDIR) if f"/{c}" in html]
    bad += [name for name in PROFILE_NAMES if name in html]
    return bad


def test_ac_13_2_redacted_page(tmp_path: Path) -> None:
    out = _out(tmp_path, LINE)
    html = _page(out, redact=True)
    assert not ac_13_2_violations(html), ac_13_2_violations(html)
    assert command_row(parse(html)) == (
        "gippa run --profile <profile> --genome <input>/nissle_r200k.fna.gz "
        "--peptides <input>/report.pr_matrix.tsv --outdir <outdir> "
        "--set db=<input>/report.pr_matrix.tsv"
    )


def test_default_page_is_unchanged(tmp_path: Path) -> None:
    """Without ``redact`` the row is the M9-02 one (file names, ``OUT``, the profile as run),
    and ``redact=False`` is the default."""
    out = _out(tmp_path, LINE)
    html = _page(out)
    assert html == _page(out, redact=False)
    assert command_row(parse(html)) == (
        f"gippa run --profile {PROFILE} --genome nissle_r200k.fna.gz "
        "--peptides report.pr_matrix.tsv --outdir OUT --set db=report.pr_matrix.tsv"
    )


def test_redact_changes_only_the_command_line(tmp_path: Path) -> None:
    """The fixture's profile names occur nowhere else, so removing them is the command line's
    job alone; every other provenance row is the same."""
    out = _out(tmp_path, LINE)
    plain, red = parse(_page(out)), parse(_page(out, redact=True))
    assert command_row(plain) != command_row(red)

    def rows(page: Page) -> list[str]:
        t = next(t for t in page.root.by_tag("table")
                 if t.attrs.get("data-table-provenance") == "run")  # fmt: skip
        return [" ".join(tr.text().split()) for tr in t.by_tag("tr")
                if "command line" not in tr.text()]  # fmt: skip

    assert rows(plain) == rows(red)


@pytest.mark.parametrize("drop", [True, False])
def test_manifest_without_command_line_has_nothing_to_rewrite(tmp_path: Path, drop: bool) -> None:
    out = _out(tmp_path, [], drop=drop)
    plain, red = parse(_page(out)), parse(_page(out, redact=True))
    assert command_row(red) == command_row(plain)


def _notes(html: str) -> list[str]:
    return [" ".join(n.text().split()) for n in parse(html).root.iter()
            if "note" in n.attrs.get("class", "").split()]  # fmt: skip


def _broken_summary(out: Path) -> Path:
    """A ``detection_summary.tsv`` whose header lacks ``metric``: the reader's error message
    names the file by its absolute path."""
    from helpers.viz_out_tree import write_schema_table

    path = out / "tables" / "detection_summary.tsv"
    write_schema_table(path, "detection_summary",
                       [{"section": "tier", "group_type": "tier", "group": "mayu",
                         "metric": "n_detected", "value": 3}])  # fmt: skip
    text = path.read_text(encoding="utf-8")
    path.write_text(text.replace("metric", "zzmetric", 1), encoding="utf-8")
    return path


def test_unreadable_table_note_is_redacted(tmp_path: Path) -> None:
    """Review r1: ``detection_summary.tsv: could not be read: …`` under ``redact`` holds no
    absolute path (neither the outdir nor any other); without ``redact`` it is unchanged (the
    reader's message, the absolute path included)."""
    out = _out(tmp_path, LINE)
    path = _broken_summary(out)
    plain, red = _notes(_page(out)), _notes(_page(out, redact=True))
    prefix = "detection_summary.tsv: could not be read:"
    assert [n for n in plain if n.startswith(prefix)] == [
        f"{prefix} ValueError: table detection_summary ({path}): schema column 'metric' is missing"
    ]
    mine = [n for n in red if n.startswith(prefix)]
    assert len(mine) == 1, red
    for note in red:
        assert str(tmp_path) not in note and str(out) not in note, note
        assert not re.search(r"(?:^|[\s=,(\[{'\":])/", note), note
    assert [n for n in plain if not n.startswith(prefix)] == [
        n for n in red if not n.startswith(prefix)
    ]


# ---------------------------------------------------------------- the CLI


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIPPA_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIPPA_ROOT", str(tmp_path / "db_empty"))


@pytest.fixture
def fake_registry(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    return swap_registry(monkeypatch)


def _cli(out: Path, *args: str) -> Result:
    from gippa.cli import cli

    fast = ["--set", "report.formats=[svg]", "--set", "report.dpi=50"]
    res = CliRunner().invoke(cli, ["report", str(out), *fast, *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def test_cli_redact(fake_registry: ModuleType, tmp_path: Path) -> None:
    """``gippa report OUT --redact``: AC 13.2 on the written page; a plain rerun restores the
    default row."""
    out = _out(tmp_path, LINE)
    res = _cli(out, "--redact")
    assert res.exit_code == 0, res.output
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert not ac_13_2_violations(html), ac_13_2_violations(html)
    assert "--profile <profile>" in command_row(parse(html))
    res = _cli(out)
    assert res.exit_code == 0, res.output
    html = (out / "report" / "report.html").read_text(encoding="utf-8")
    assert f"--profile {PROFILE}" in command_row(parse(html))


def test_cli_help_lists_redact() -> None:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["report", "--help"])
    assert res.exit_code == 0 and "--redact" in res.output, res.output


# ---------------------------------------------------------------- property


_SEG = st.text(alphabet="abcdefgxyz0123456789_.-", min_size=1, max_size=8).filter(
    lambda s: s not in (".", "..") and not s.startswith("-")
)


#: the ``--set`` YAML forms of review r1 M1: ``{p}`` is the path, ``{r}`` a relative item
_SET_FORMS = {
    "list": "k=[{r},{p}]",
    "list_first": "k=[{p},{r}]",
    "squote": "k='{p}'",
    "dquote": 'k="{p}"',
    "brace": "k={{d: {p}}}",
    "colon": "k=x:{p}",
}


@st.composite
def _abs_arg(draw: st.DrawFn) -> tuple[str, str, str]:
    """``(token, the absolute path inside it, the token redacted)`` in one of the path forms."""
    under = draw(st.booleans())
    segs = draw(st.lists(_SEG, min_size=1, max_size=4))
    rel = draw(_SEG)
    path = "/".join([OUTDIR if under else "/zzhome/zzuser", *segs])
    form = draw(st.sampled_from(["bare", "eq", "set_eq", "comma", "short", "eq_in_path",
                                 *_SET_FORMS]))  # fmt: skip
    if form == "eq_in_path":  # ``--set=key=/home/x=y/z``: the path is ``/home/x=y/z``
        path = f"{path}=y/{rel}"
    shown = ("<outdir>/" + path[len(OUTDIR) + 1 :] if under
             else "<input>/" + PurePosixPath(path).name)  # fmt: skip
    fmt = {"bare": "{p}", "eq": "--genome={p}", "set_eq": "--set=k={p}", "comma": "{r},{p}",
           "short": "-o{p}", "eq_in_path": "--set=key={p}", **_SET_FORMS}[form]  # fmt: skip
    return fmt.format(p=path, r=rel), path, fmt.format(p=shown, r=rel)


@settings(max_examples=150, deadline=None)
@given(parts=st.lists(_abs_arg(), min_size=1, max_size=5),
       profile=st.lists(st.from_regex(r"zz[a-z]{3,8}", fullmatch=True), min_size=1, max_size=3),
       eq=st.booleans())  # fmt: skip
def test_property_nothing_local_survives(
    parts: list[tuple[str, str, str]], profile: list[str], eq: bool
) -> None:
    """No absolute path, no directory component of one, no profile name; each path's basename
    stays; the line is exactly the tokens with each path rewritten in place (review r1 M1)."""
    value = ",".join(profile)
    prof = [f"--profile={value}"] if eq else ["--profile", value]
    got = redact([*prof, *(t for t, _, _ in parts)])
    shown_prof = "--profile=<profile>" if eq else "--profile <profile>"
    assert got == " ".join(["gippa run", shown_prof, *(e for _, _, e in parts)])
    assert not ABS.search(got), got
    # the local prefixes (the drawn segments use none of their letters h, l, m, o, r, s, u, A)
    for comp in dir_components(OUTDIR + "/x", "/zzhome/zzuser/x"):
        assert f"/{comp}" not in got and (comp == "A" or comp not in got), (comp, got)
    for _, path, _ in parts:
        assert PurePosixPath(path).name in got
    assert not [p for p in profile if p in got], got
    assert "<profile>" in got
