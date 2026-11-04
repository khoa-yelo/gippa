"""M12-01 PR 4: ``gippa example`` help and the CLI reference (section 1, section 4 > CLI reference,
AC 16, Q-L4).

``gen_reference.py render_cli`` recurses into groups and renders ``## gippa example fetch`` with its
option table. The option-drift parser of ``tests/unit/test_m11_08_docs.py`` reads one-word
``## gippa <name>`` sections; this file carries its own multi-word parser
(``helpers.example_m12_01.cli_reference``) and checks it on the real ``cli.md``.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import click
import pytest
from helpers.example_m12_01 import CLI_MD, REPO, cli_reference, gippa, option_problems

SCRIPT = REPO / "scripts" / "docs" / "gen_reference.py"


@pytest.fixture(scope="module")
def gr() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gen_reference_m12_01", SCRIPT)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules["gen_reference_m12_01"] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- help


@pytest.mark.parametrize("argv", [["example"], ["example", "--help"]])
def test_group_alone_prints_its_help(argv: list[str]) -> None:
    """``gippa example`` alone prints the group help, exit 0."""
    res = gippa(argv)
    assert res.code == 0, res.err
    assert "fetch" in res.out


def test_fetch_help() -> None:
    res = gippa(["example", "fetch", "--help"])
    assert res.code == 0, res.err
    for opt in ("--outdir", "--version", "--from-file"):
        assert opt in res.out, opt
    assert "bundle version" in " ".join(res.out.split())
    assert "gippa-example" in res.out  # the --outdir default


def test_root_version_is_unchanged() -> None:
    import gippa as pkg

    res = gippa(["--version"])
    assert res.code == 0 and pkg.__version__ in res.out


def test_fetch_version_takes_a_value() -> None:
    """``fetch --version`` is the bundle version (a value), not the shared ``--version`` flag."""
    from gippa.cli import cli

    fetch = cli.commands["example"].commands["fetch"]  # type: ignore[attr-defined]
    opts = {o: p for p in fetch.params for o in getattr(p, "opts", ())}
    assert not opts["--version"].is_flag
    assert opts["--outdir"].default in ("gippa-example", Path("gippa-example"))
    assert not opts["--from-file"].required


# ---------------------------------------------------------------- render_cli recursion


def _toy() -> click.Group:
    @click.group()
    def toy() -> None:
        """Toy tool."""

    @toy.command()
    def alpha() -> None:
        """Run alpha."""

    @toy.group()
    def grp() -> None:
        """A group."""

    @grp.command()
    @click.option("--outdir", default="out", show_default=True, help="Where.")
    @click.option("--version", "bundle", help="Bundle version.")
    def sub(**_: Any) -> None:
        """Run sub."""

    @grp.command("_secret", hidden=True)
    def secret() -> None:
        """Never documented."""

    @toy.command()
    def zeta() -> None:
        """Run zeta."""

    return toy


def _sections(md: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for block in md.split("\n## ")[1:]:
        head, _, body = block.partition("\n")
        out[head.strip()] = body
    return out


def test_render_cli_recurses_into_groups(gr: ModuleType) -> None:
    md = gr.render_cli(_toy(), "toy")
    secs = _sections(md)
    assert "toy grp sub" in secs, list(secs)
    assert "Run sub." in secs["toy grp sub"]
    rows = [ln for ln in secs["toy grp sub"].splitlines() if ln.startswith("| `")]
    assert any(r.startswith("| `--outdir` |") and "`out`" in r for r in rows), rows
    assert any(r.startswith("| `--version` |") and "Bundle version." in r for r in rows), rows
    heads = list(secs)
    assert heads.index("toy alpha") < heads.index("toy grp sub") < heads.index("toy zeta")
    assert not [h for h in heads if "_secret" in h]


def test_render_cli_without_groups_is_unchanged(gr: ModuleType) -> None:
    @click.group()
    def flat() -> None:
        """Flat."""

    @flat.command()
    @click.option("--n", type=int, help="N.")
    def one(n: int) -> None:
        """One."""

    assert list(_sections(gr.render_cli(flat, "flat"))) == ["flat one"]


# ---------------------------------------------------------------- docs/reference/cli.md (AC 16)


def test_cli_md_has_the_fetch_section() -> None:
    ref = cli_reference(CLI_MD.read_text(encoding="utf-8"))
    assert ("example", "fetch") in ref, sorted(ref)
    assert {"--outdir", "--version", "--from-file"} <= ref[("example", "fetch")]


def test_cli_md_is_current(gr: ModuleType) -> None:
    """``gen_reference.py --check``: the committed cli.md is what render_cli writes."""
    from gippa.cli import cli

    assert CLI_MD.read_text(encoding="utf-8") == gr.render_cli(cli, "gippa")


def test_multiword_drift_parser() -> None:
    """The parser this PR needs (Q-L4): the longest ``## gippa <words>`` section wins."""
    ref = cli_reference(CLI_MD.read_text(encoding="utf-8"))
    assert option_problems("gippa example fetch --outdir x --from-file f", ref) == []
    assert option_problems("gippa example fetch --threads 4", ref) == ["example fetch: --threads"]
    assert option_problems("gippa run --genome g --threads 4", ref) == []
