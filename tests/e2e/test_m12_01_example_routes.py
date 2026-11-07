"""M12-01 slow route test (spec docs/specs/M12-01-example-dataset.md section 5, AC 17); Slurm
only, via scripts/slurm/example_routes.sbatch (its header lists the environment variables).

Fetch mode: ``gippa example fetch --from-file <bundle>``, the printed commands per route, each run
compared with ``expected/<route>.json``. Measure mode (``GIPPA_EXAMPLE_DIR``, release step 6): a
copied build tree, commands from ``gippa.example.route_lines``, no comparison. Both write
``runs/wall_times.tsv`` (never asserted) and ``runs/measured/<route>.json`` (``expected/`` layout).
``$GIPPA_EXAMPLE_RUN_ARGS`` goes on ``gippa run`` only. Route B is included (#536 closed).
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest
from helpers import example_routes as er
from helpers.example_route_e import route_e_misses

pytestmark = [
    pytest.mark.e2e,
    pytest.mark.slow,
    pytest.mark.needs_apptainer,
    pytest.mark.needs_bakta_db,
    pytest.mark.needs_promotech_models,
]

TIERS = {"A": "mayu", "B": "protein_list", "C": "two_peptide", "E": "mayu"}  # AC 17, Q-L10
TIMEOUT_S = 3600  # one command; the measured `gippa run` takes about 5 min
STEP6 = ("before the pin exists (release step 6), measure an unpacked build tree with "
         "GIPPA_EXAMPLE_DIR=<dir> instead")  # fmt: skip
ROOT = Path(os.environ.get("GIPPA_ROOT") or Path(__file__).resolve().parents[2] / "local")
DB = Path(os.environ.get("GIPPA_DB_DIR") or ROOT / "db")


def _skip_or_fail(reason: str) -> None:
    if os.environ.get("GIPPA_EXAMPLE_REQUIRE") == "1":
        pytest.fail(f"GIPPA_EXAMPLE_REQUIRE=1: {reason}", pytrace=False)
    pytest.skip(reason)


def _version() -> str:
    import gippa

    m = re.match(r"\d+\.\d+", os.environ.get("GIPPA_EXAMPLE_VERSION") or gippa.__version__)
    assert m, gippa.__version__
    return m.group(0)


@dataclass
class Bundle:
    directory: Path
    lines: list[str]  # the printed block
    fetched: bool
    gippa: str


@dataclass
class RouteRun:
    exits: list[int] = field(default_factory=list)
    logs: list[str] = field(default_factory=list)
    measured: dict[str, Any] = field(default_factory=dict)


@pytest.fixture(scope="module")
def bundle(tmp_path_factory: pytest.TempPathFactory) -> Bundle:
    from gippa import example

    antismash = DB / "antismash"
    if not (antismash.is_dir() and any(antismash.iterdir())):
        _skip_or_fail(f"no antiSMASH databases at {antismash} (bgc; gippa setup --no-genomad)")
    gippa_bin = shutil.which("gippa")
    if gippa_bin is None:
        _skip_or_fail("no `gippa` on PATH")
    assert gippa_bin is not None
    base = Path(os.environ.get("GIPPA_EXAMPLE_OUT") or tmp_path_factory.mktemp("example"))
    directory = (base / "gippa-example").resolve()
    assert not directory.exists(), f"{directory} exists: give GIPPA_EXAMPLE_OUT a fresh directory"
    base.mkdir(parents=True, exist_ok=True)
    if tree := os.environ.get("GIPPA_EXAMPLE_DIR"):
        assert Path(tree, "genome").is_dir(), f"GIPPA_EXAMPLE_DIR={tree}: not a build tree"
        shutil.copytree(tree, directory, ignore=shutil.ignore_patterns("runs"))
        lines = example.route_lines(_version(), directory, example._routes())
        return Bundle(directory, lines, False, gippa_bin)
    archive = Path(os.environ.get("GIPPA_EXAMPLE_BUNDLE")
                   or ROOT / "example" / f"gippa-example-{_version()}.tar.gz")  # fmt: skip
    if not archive.is_file():
        _skip_or_fail(f"no example bundle at {archive}; {STEP6}")
    cmd = [gippa_bin, "example", "fetch", "--from-file", str(archive), "--outdir", str(directory)]
    if os.environ.get("GIPPA_EXAMPLE_VERSION"):
        cmd += ["--version", _version()]
    p = subprocess.run(cmd, capture_output=True, text=True, check=False, cwd=base)
    if p.returncode != 0:
        msg = f"{shlex.join(cmd)} -> exit {p.returncode}: {p.stderr.strip()}"
        if "no example bundle" in p.stderr:
            _skip_or_fail(f"{msg}; {STEP6}")
        pytest.fail(msg, pytrace=False)
    return Bundle(directory, p.stdout.splitlines(), True, gippa_bin)


def _run(argv: list[str], cwd: Path) -> tuple[int, str]:
    try:
        p = subprocess.run(argv, capture_output=True, text=True, check=False, cwd=cwd,
                           timeout=TIMEOUT_S)  # fmt: skip
        return p.returncode, p.stderr
    except subprocess.TimeoutExpired as e:
        err = e.stderr.decode(errors="replace") if isinstance(e.stderr, bytes) else e.stderr
        return 124, f"timed out after {TIMEOUT_S} s\n{err or ''}"


@pytest.fixture(scope="module")
def runs(bundle: Bundle) -> dict[str, RouteRun]:
    """Every route's printed commands in order (a route stops at its first failure)."""
    from gippa import example

    printed = er.printed_routes(bundle.lines, example._routes())
    assert sorted(printed) == list(er.ROUTES), bundle.lines
    run_args = shlex.split(os.environ.get("GIPPA_EXAMPLE_RUN_ARGS",
                                          "--profile farmshare_local,apptainer"))  # fmt: skip
    out = bundle.directory / "runs"
    (out / "measured").mkdir(parents=True, exist_ok=True)
    result: dict[str, RouteRun] = {}
    with (out / "wall_times.tsv").open("w", encoding="utf-8") as tsv:
        tsv.write("route\tstep\tcommand\texit\twall_s\n")
        for route in er.ROUTES:
            rr = result[route] = RouteRun()
            total, peptides, outdir = 0.0, None, out / route
            for i, line in enumerate(printed[route], 1):
                argv = er.argv(line, bundle.gippa, run_args)
                peptides = er.option(argv, "--out") if argv[1] == "convert" else peptides
                outdir = er.option(argv, "--outdir") if argv[1] == "run" else outdir
                t0 = time.monotonic()
                code, err = _run(argv, out)
                wall = time.monotonic() - t0
                total += wall
                rr.exits.append(code)
                rr.logs.append(f"$ {shlex.join(argv)}\n-> {code}\n{err[-4000:]}")
                tsv.write(f"{route}\t{i}\t{argv[1]}\t{code}\t{wall:.1f}\n")
                tsv.flush()
                if code != 0:
                    break
            tsv.write(f"{route}\ttotal\t\t{max(rr.exits, default=0)}\t{total:.1f}\n")
            tsv.flush()
            assert outdir is not None
            rr.measured = er.measure(route, outdir, peptides) | {"wall_time_s": round(total, 1)}
            (out / "measured" / f"{route}.json").write_text(
                json.dumps(er.as_expected(rr.measured, route), indent=2) + "\n", encoding="utf-8"
            )
    return result


def test_fetch_printed_this_bundle(bundle: Bundle) -> None:
    if not bundle.fetched:
        pytest.skip("measure mode (GIPPA_EXAMPLE_DIR): no fetch")
    assert bundle.lines[0] == f"gippa example {_version()}: {bundle.directory}"


@pytest.mark.parametrize("route", er.ROUTES)
def test_route(runs: dict[str, RouteRun], bundle: Bundle, route: str) -> None:
    rr = runs[route]
    log = "\n".join(rr.logs)
    assert rr.exits and not any(rr.exits), log
    assert rr.measured.get("status") == "completed", (rr.measured, log)
    if route in TIERS:
        assert rr.measured.get("detection_tier") == TIERS[route], rr.measured
    if route == "E":  # AC 22 (amended 2026-11-06 (b)); route_e_misses is unit-tested
        assert not (misses := route_e_misses(rr.measured)), misses
    if bundle.fetched:
        path = bundle.directory / "expected" / f"{route}.json"
        assert path.is_file(), f"{path} missing from the bundle"
        misses = er.compare(rr.measured, json.loads(path.read_text(encoding="utf-8")), route)
        assert not misses, f"route {route} vs {path}:\n" + "\n".join(misses)
