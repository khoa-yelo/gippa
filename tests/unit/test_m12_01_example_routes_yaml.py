"""M12-01 PR 4: the printed route commands (``src/gippa/data/example_routes.yaml``; AC 14, 17b).

Section 2 > Routes: route A converts the DIA-NN report, B runs the protein FASTA plus abundance,
C converts the Cascadia ``.ssl`` files, D converts the generic peptide table; each ends with
``gippa preflight`` and ``gippa run``. Every ``gippa run`` carries the settings Khoa accepted
(2026-11-03): ``--modules core,bgc --set annotation.bakta.db=light --terminators transtermhp
--set annotation.promoter_calculator.enabled=false --threads 4``.
"""

from __future__ import annotations

import re

import pytest
from helpers.example_m12_01 import (
    CLI_MD,
    EXAMPLE_SETTINGS,
    ROUTES_YAML,
    cli_reference,
    load_routes,
    option_problems,
    options,
    subcommand,
)

ROUTE_IDS = ("A", "B", "C", "D", "E")  # route E, 2026-11-06
FORBIDDEN = ("--annotation", "--strict", "--resume", "--executor")
INPUT_OPTS = ("--genome", "--peptides", "--proteins", "--abundance")


@pytest.fixture(scope="module")
def routes() -> dict[str, dict[str, object]]:
    assert ROUTES_YAML.is_file(), "src/gippa/data/example_routes.yaml does not exist"
    return load_routes()


def _cmds(routes: dict[str, dict[str, object]], rid: str, sub: str | None = None) -> list[str]:
    cmds: list[str] = routes[rid]["commands"]  # type: ignore[assignment]
    return [c for c in cmds if sub is None or subcommand(c) == sub]


def _all(routes: dict[str, dict[str, object]]) -> list[str]:
    return [c for rid in routes for c in _cmds(routes, rid)]


def test_five_routes_with_headings(routes: dict[str, dict[str, object]]) -> None:
    assert tuple(sorted(routes)) == ROUTE_IDS
    for rid in ROUTE_IDS:
        assert str(routes[rid]["title"]).strip(), rid


def test_commands_start_with_gippa_and_use_only_the_placeholders(
    routes: dict[str, dict[str, object]],
) -> None:
    for cmd in _all(routes):
        assert cmd.startswith("gippa "), cmd
        assert set(re.findall(r"\{[^}]*\}", cmd)) <= {"{dir}", "{out}"}, cmd
        assert "{" not in re.sub(r"\{(dir|out)\}", "", cmd), cmd


def test_command_sequence_per_route(routes: dict[str, dict[str, object]]) -> None:
    """About three commands: convert (A, C, D), preflight, run, in that order."""
    for rid in ROUTE_IDS:
        subs = [subcommand(c) for c in _cmds(routes, rid)]
        want = ["preflight", "run"] if rid == "B" else ["convert", "preflight", "run"]
        assert subs == want, (rid, subs)


def test_route_inputs(routes: dict[str, dict[str, object]]) -> None:
    conv = {r: options(_cmds(routes, r, "convert")[0]) for r in "ACDE"}
    assert conv["A"]["--from"] == ["diann"] and conv["A"]["--search-db"] == ["same"]
    assert conv["C"]["--from"] == ["cascadia"]
    assert conv["D"]["--from"] == ["generic"]
    assert conv["D"]["--source"] == ["database"] and conv["D"]["--search-db"] == ["same"]
    assert conv["E"]["--from"] == ["diann"] and conv["E"]["--search-db"] == ["same"]
    for r, path in (("A", "{dir}/A_diann/report.parquet"), ("C", "{dir}/C_cascadia/"),
                    ("D", "{dir}/D_generic/peptides.tsv"),
                    ("E", "{dir}/E_sixframe/report.parquet")):  # fmt: skip
        cmd = _cmds(routes, r, "convert")[0]
        assert path in cmd, (r, cmd)
        assert conv[r]["--out"] == [f"{{out}}/{r}/peptides.tsv"], (r, conv[r])
    run_b = options(_cmds(routes, "B", "run")[0])
    assert run_b["--proteins"] == ["{dir}/B_proteins/proteins.faa"]
    assert run_b["--abundance"] == ["{dir}/B_proteins/abundance.tsv"]
    assert run_b["--abundance-type"] == ["maxlfq"]
    assert "--peptides" not in run_b


@pytest.mark.parametrize("rid", ROUTE_IDS)
def test_every_run_carries_the_example_settings(
    routes: dict[str, dict[str, object]], rid: str
) -> None:
    """AC 17b."""
    (cmd,) = _cmds(routes, rid, "run")
    for setting in EXAMPLE_SETTINGS:
        assert setting in cmd, (rid, setting, cmd)
    opts = options(cmd)
    (genome,) = opts["--genome"]
    assert genome.startswith("{dir}/genome/") and genome.endswith(".fna.gz"), genome
    assert opts["--genome-relation"] == ["same"]
    assert opts["--outdir"] == [f"{{out}}/{rid}"]
    assert opts["--threads"] == ["4"] and opts["--modules"] == ["core,bgc"]
    if rid != "B":
        assert opts["--peptides"] == [f"{{out}}/{rid}/peptides.tsv"]


@pytest.mark.parametrize("rid", ROUTE_IDS)
def test_preflight_checks_the_inputs_of_the_run(
    routes: dict[str, dict[str, object]], rid: str
) -> None:
    pre = options(_cmds(routes, rid, "preflight")[0])
    run = options(_cmds(routes, rid, "run")[0])
    for opt in INPUT_OPTS:
        if opt in run:
            assert pre.get(opt) == run[opt], (rid, opt, pre.get(opt), run[opt])


_SETTING_KEYS = ("modules", "annotation.bakta.db", "annotation.terminators.tool",
                 "annotation.promoter_calculator.enabled")  # fmt: skip


def _settings(cmd: str) -> dict[str, str]:
    """The run settings of one command: ``--set k=v``, ``--modules`` and ``--terminators``, with
    ``modules=[a,b]`` read as ``a,b``."""
    opts = options(cmd)
    out: dict[str, str] = {}
    for kv in opts.get("--set", []):
        key, _, value = kv.partition("=")
        out[key] = value.strip("[]") if key == "modules" else value
    for opt, key in (("--modules", "modules"), ("--terminators", "annotation.terminators.tool")):
        for value in opts.get(opt, []):
            assert out.setdefault(key, value) == value, (cmd, key)
    return out


@pytest.mark.parametrize("rid", ROUTE_IDS)
def test_preflight_checks_the_settings_of_the_run(
    routes: dict[str, dict[str, object]], rid: str
) -> None:
    """Review #542: preflight gets the same modules, Bakta DB, terminator tool and Promoter
    Calculator setting as its run (no drift between the two commands)."""
    pre = _settings(_cmds(routes, rid, "preflight")[0])
    run = _settings(_cmds(routes, rid, "run")[0])
    for key in _SETTING_KEYS:
        assert key in run and pre.get(key) == run[key], (rid, key, pre.get(key), run.get(key))


def test_no_forbidden_option(routes: dict[str, dict[str, object]]) -> None:
    for cmd in _all(routes):
        for opt in FORBIDDEN:
            assert not re.search(rf"(?<![\w-]){opt}(?![\w-])", cmd), (opt, cmd)


def test_options_exist_in_the_cli_reference(routes: dict[str, dict[str, object]]) -> None:
    """AC 14: every command names only options of its ``docs/reference/cli.md`` section."""
    ref = cli_reference(CLI_MD.read_text(encoding="utf-8"))
    problems = [p for cmd in _all(routes) for p in option_problems(cmd, ref)]
    assert problems == []
