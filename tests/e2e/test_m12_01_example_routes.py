"""M12-01 slow route test: the example bundle's four routes, run for real (spec
docs/specs/M12-01-example-dataset.md section 5 "Slow route test", AC 17).

Submitted with ``sbatch scripts/slurm/example_routes.sbatch <repo-or-worktree>`` (job
``gippa-example-routes``, 4 CPUs, 16 GB, R14); never in CI or the session allocation.

1. ``gippa example fetch --from-file $GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz --outdir DIR``
   (the offline path; ``X.Y`` is the installed GIPPA's ``MAJOR.MINOR``, or
   ``$GIPPA_EXAMPLE_VERSION``, or the file is ``$GIPPA_EXAMPLE_BUNDLE``).
2. Each route's commands are read from the fetch's stdout (the printed block, one ``# <title>``
   heading per route, the titles of ``gippa/data/example_routes.yaml``), never typed here a second
   time. The block must equal ``gippa.example.route_lines`` of that YAML. The commands run in order;
   a route stops at its first non-zero exit.
   The only additions are ``$GIPPA_EXAMPLE_RUN_ARGS`` appended to each ``gippa run`` (default
   ``--profile farmshare_local,apptainer``: the job's own CPUs, as the section 3 measurements;
   ``auto`` inside a Slurm job would submit every task as a further Slurm job).
3. Each run is compared with ``DIR/expected/<route>.json`` by the tolerances stored there.
4. ``DIR/runs/wall_times.tsv`` (next to the outputs) gets one row per command and one ``total``
   row per route. Times are recorded, never asserted. ``DIR/runs/measured/<route>.json`` holds the
   measured keys of the section 2 table (the input for filling ``expected/`` at release prep).

``expected/<route>.json``: per key either ``{"value": v, "tolerance": t}``, or the values at top
level with a ``"tolerance": {key: t}`` map. ``t`` is ``"exact"``, ``"set"`` (order-free exact),
``"informational"`` (not compared), ``{"abs": a, "rel": r}`` meaning
``|got - expected| <= max(a, r * |expected|)``, or the spec's text form ``"max(1, 2 %)"``. A
compared key without a tolerance is a failure (the file must say how to compare it).

Route B (``--proteins`` + ``--abundance``, tier ``protein_list``) is included: #536 (MAPPING ran
only with ``--peptides``) is closed.

Without the bundle file the tests skip with the path in the reason; ``GIPPA_EXAMPLE_REQUIRE=1``
(set by the sbatch) turns that skip into a failure. ``GIPPA_EXAMPLE_OUT`` keeps the outputs (the
sbatch puts them under ``$GIPPA_ROOT/example/routes/<job>``); otherwise a pytest temp dir is used.
"""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

pytestmark = [
    pytest.mark.e2e,
    pytest.mark.slow,
    pytest.mark.needs_apptainer,
    pytest.mark.needs_bakta_db,
    pytest.mark.needs_promotech_models,
]

ROUTES = ("A", "B", "C", "D")
TIERS = {"A": "mayu", "B": "protein_list", "C": "two_peptide"}  # AC 17
CONVERTED = ("A", "C", "D")  # routes whose first command is `gippa convert`
ABUNDANCE = ("A", "B", "D")
DEFAULT_RUN_ARGS = "--profile farmshare_local,apptainer"
TIMEOUT_S = 3600  # one command; the measured `gippa run` is about 5 min


def _root() -> Path:
    return Path(os.environ.get("GIPPA_ROOT") or Path(__file__).resolve().parents[2] / "local")


def _db_dir() -> Path:
    return Path(os.environ.get("GIPPA_DB_DIR") or _root() / "db")


def _bundle_version() -> str:
    if v := os.environ.get("GIPPA_EXAMPLE_VERSION"):
        return v
    import gippa

    m = re.match(r"(\d+)\.(\d+)", gippa.__version__)
    assert m, gippa.__version__
    return f"{m.group(1)}.{m.group(2)}"


def _bundle_file() -> Path:
    if f := os.environ.get("GIPPA_EXAMPLE_BUNDLE"):
        return Path(f)
    return _root() / "example" / f"gippa-example-{_bundle_version()}.tar.gz"


def _skip_or_fail(reason: str) -> None:
    if os.environ.get("GIPPA_EXAMPLE_REQUIRE") == "1":
        pytest.fail(f"GIPPA_EXAMPLE_REQUIRE=1: {reason}", pytrace=False)
    pytest.skip(reason)


# ------------------------------------------------------------------ expected/ and tolerances


def _within(got: float, exp: float, tol: Any) -> bool:
    if isinstance(tol, str):
        m = re.fullmatch(r"\s*max\(\s*([\d.]+)\s*,\s*([\d.]+)\s*%\s*\)\s*", tol)
        assert m, f"unreadable tolerance {tol!r}"
        tol = {"abs": float(m.group(1)), "rel": float(m.group(2)) / 100}
    assert isinstance(tol, Mapping) and set(tol) <= {"abs", "rel"}, f"unknown tolerance {tol!r}"
    bound = max(float(tol.get("abs", 0)), float(tol.get("rel", 0)) * abs(exp))
    return abs(got - exp) <= bound


def _expected(doc: Mapping[str, Any]) -> dict[str, tuple[Any, Any]]:
    """``{key: (value, tolerance)}`` from either layout of ``expected/<route>.json``."""
    out: dict[str, tuple[Any, Any]] = {}
    flat_tol = doc.get("tolerance") if isinstance(doc.get("tolerance"), Mapping) else {}
    for key, val in doc.items():
        if key == "tolerance" and flat_tol:
            continue
        if isinstance(val, Mapping) and "value" in val:
            out[key] = (val["value"], val.get("tolerance"))
        else:
            out[key] = (val, flat_tol.get(key))
    return out


def compare(measured: Mapping[str, Any], doc: Mapping[str, Any]) -> list[str]:
    """Every compared key of ``doc`` against ``measured``; the list of misses (empty: all match)."""
    misses: list[str] = []
    for key, (exp, tol) in sorted(_expected(doc).items()):
        if key in ("route", "version", "gippa_version", "commit", "note", "notes"):
            continue
        if tol in ("informational", "info", "not compared") or key in ("warnings", "wall_time_s"):
            continue
        if key not in measured:
            misses.append(f"{key}: not measured (expected {exp!r})")
            continue
        got = measured[key]
        if tol is None:
            misses.append(f"{key}: no tolerance in expected/ (expected {exp!r}, got {got!r})")
        elif tol == "exact":
            if got != exp:
                misses.append(f"{key}: got {got!r}, expected exactly {exp!r}")
        elif tol == "set":
            if sorted(map(str, got or [])) != sorted(map(str, exp or [])):
                misses.append(f"{key}: got set {sorted(got or [])}, expected {sorted(exp or [])}")
        elif not _within(float(got), float(exp), tol):
            misses.append(f"{key}: got {got}, expected {exp} within {tol}")
    return misses


# ------------------------------------------------------------------ measuring a run


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        rows = (line for line in fh if not line.startswith("#"))
        return list(csv.DictReader(rows, delimiter="\t"))


def _true(v: str | None) -> bool:
    return (v or "").strip().lower() == "true"


def measure(route: str, outdir: Path) -> dict[str, Any]:
    """The section 2 `expected/` keys of one finished route, read from its outputs."""
    got: dict[str, Any] = {}
    man_path = outdir / "run_manifest.json"
    if man_path.is_file():
        man = json.loads(man_path.read_text(encoding="utf-8"))
        got["status"] = man.get("status")
        got["detection_tier"] = man.get("detection_tier")
        got["layers_missing"] = sorted(man.get("layers_missing") or [])
        got["warnings"] = sorted(man.get("warnings") or [])
    if route in CONVERTED and (outdir / "peptides.tsv").is_file():
        got["convert_rows"] = len(_tsv(outdir / "peptides.tsv"))
    genes_path = outdir / "tables" / "genes.tsv"
    if genes_path.is_file():
        genes = _tsv(genes_path)
        got["n_genes"] = len(genes)
        got["n_detected"] = sum(_true(g.get("detected")) for g in genes)
        got["n_detected_2pep"] = sum(_true(g.get("detected_2pep")) for g in genes)
        if route in ABUNDANCE:
            got["n_abundance"] = sum(
                (g.get("abundance_raw") or "NA") not in ("", "NA")
                and not _true(g.get("abundance_imputed"))
                for g in genes
            )
    hits = outdir / "mapping" / "peptide_hits.tsv"
    if hits.is_file():
        got["n_peptides_placed"] = len({h["peptide"] for h in _tsv(hits)})
    return got


# ------------------------------------------------------------------ running the routes


@dataclass
class RouteRun:
    route: str
    commands: list[list[str]]
    exits: list[int] = field(default_factory=list)
    logs: list[str] = field(default_factory=list)
    wall_s: list[float] = field(default_factory=list)


def printed_routes(stdout: str, routes: list[Mapping[str, Any]]) -> dict[str, list[str]]:
    """The fetch's stdout block -> ``{route id: [command lines]}`` by the YAML's titles."""
    by_title = {f"# {r['title']}": str(r["id"]) for r in routes}
    out: dict[str, list[str]] = {}
    current: str | None = None
    for line in stdout.splitlines():
        if line in by_title:
            current = by_title[line]
            out[current] = []
        elif current is not None and line.startswith("gippa "):
            out[current].append(line)
        elif not line.strip():
            current = None
    return out


def _env() -> dict[str, str]:
    env = dict(os.environ)
    env.setdefault("GIPPA_DB_DIR", str(_db_dir()))
    return env


def _argv(line: str, gippa_bin: str, run_args: list[str]) -> list[str]:
    words = shlex.split(line)
    assert words[0] == "gippa", line
    if words[1] == "run":
        words += run_args
    return [gippa_bin, *words[1:]]


@pytest.fixture(scope="module")
def fetched(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, str]:
    """``(DIR, stdout)`` after the offline fetch of the bundle."""
    bundle = _bundle_file()
    if not bundle.is_file():
        _skip_or_fail(f"no example bundle at {bundle} (built at release prep, spec section 3)")
    antismash = _db_dir() / "antismash"
    if not (antismash.is_dir() and any(antismash.iterdir())):
        _skip_or_fail(f"no antiSMASH databases at {antismash} (bgc; gippa setup --no-genomad)")
    gippa_bin = shutil.which("gippa")
    if gippa_bin is None:
        _skip_or_fail("no `gippa` on PATH")
    assert gippa_bin is not None
    base = Path(os.environ.get("GIPPA_EXAMPLE_OUT") or tmp_path_factory.mktemp("example"))
    base.mkdir(parents=True, exist_ok=True)
    directory = (base / "gippa-example").resolve()
    assert not directory.exists(), f"{directory} exists: give GIPPA_EXAMPLE_OUT a fresh directory"
    argv = [gippa_bin, "example", "fetch", "--from-file", str(bundle), "--outdir", str(directory)]
    if os.environ.get("GIPPA_EXAMPLE_VERSION"):
        argv += ["--version", _bundle_version()]
    p = subprocess.run(argv, capture_output=True, text=True, env=_env(), check=False, cwd=base)
    assert p.returncode == 0, f"{shlex.join(argv)} -> {p.returncode}\n{p.stderr}"
    return directory, p.stdout


@pytest.fixture(scope="module")
def runs(fetched: tuple[Path, str]) -> dict[str, RouteRun]:
    """Every route's printed commands, run in order; wall times written to the TSV."""
    from gippa import example

    directory, stdout = fetched
    gippa_bin = shutil.which("gippa")
    assert gippa_bin is not None
    run_args = shlex.split(os.environ.get("GIPPA_EXAMPLE_RUN_ARGS", DEFAULT_RUN_ARGS))
    printed = printed_routes(stdout, example._routes())
    assert sorted(printed) == list(ROUTES), stdout
    out_dir = directory / "runs"
    out_dir.mkdir(exist_ok=True)
    (out_dir / "measured").mkdir(exist_ok=True)
    result: dict[str, RouteRun] = {}
    rows: list[list[str]] = []
    for route in ROUTES:
        rr = RouteRun(route, [_argv(c, gippa_bin, run_args) for c in printed[route]])
        result[route] = rr
        for i, argv in enumerate(rr.commands):
            t0 = time.monotonic()
            p = subprocess.run(argv, capture_output=True, text=True, env=_env(), check=False,
                               cwd=out_dir, timeout=TIMEOUT_S)  # fmt: skip
            rr.wall_s.append(time.monotonic() - t0)
            rr.exits.append(p.returncode)
            rr.logs.append(f"$ {shlex.join(argv)}\n-> {p.returncode}\n{p.stderr[-4000:]}")
            rows.append([route, str(i + 1), argv[1], str(p.returncode), f"{rr.wall_s[-1]:.1f}"])
            if p.returncode != 0:
                break
        rows.append([route, "total", "", str(max(rr.exits, default=0)), f"{sum(rr.wall_s):.1f}"])
        m = measure(route, out_dir / route)
        m["wall_time_s"] = round(sum(rr.wall_s), 1)
        (out_dir / "measured" / f"{route}.json").write_text(
            json.dumps(m, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
    with (out_dir / "wall_times.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["route", "step", "command", "exit", "wall_s"])
        w.writerows(rows)
    return result


# ------------------------------------------------------------------ tests


def test_printed_block_equals_example_routes_yaml(fetched: tuple[Path, str]) -> None:
    """The commands run below are the ones `fetch` printed, which are example_routes.yaml's."""
    from gippa import example

    directory, stdout = fetched
    version = stdout.splitlines()[0].split(":")[0].removeprefix("gippa example ").strip()
    want = example.route_lines(version, directory, example._routes())
    assert stdout.splitlines() == want


@pytest.mark.parametrize("route", ROUTES)
def test_route_commands_exit_0(runs: dict[str, RouteRun], route: str) -> None:
    rr = runs[route]
    assert rr.commands, f"route {route}: no printed commands"
    assert len(rr.exits) == len(rr.commands) and not any(rr.exits), "\n".join(rr.logs)


@pytest.mark.parametrize("route", ROUTES)
def test_route_status_completed(
    runs: dict[str, RouteRun], fetched: tuple[Path, str], route: str
) -> None:
    m = measure(route, fetched[0] / "runs" / route)
    assert m.get("status") == "completed", (m, "\n".join(runs[route].logs))


@pytest.mark.parametrize("route", sorted(TIERS))
def test_route_detection_tier(
    runs: dict[str, RouteRun], fetched: tuple[Path, str], route: str
) -> None:
    m = measure(route, fetched[0] / "runs" / route)
    assert m.get("detection_tier") == TIERS[route], m


@pytest.mark.parametrize("route", ROUTES)
def test_route_matches_expected_within_tolerance(
    runs: dict[str, RouteRun], fetched: tuple[Path, str], route: str
) -> None:
    directory, _ = fetched
    path = directory / "expected" / f"{route}.json"
    assert path.is_file(), f"{path} missing from the bundle"
    doc = json.loads(path.read_text(encoding="utf-8"))
    measured = measure(route, directory / "runs" / route)
    misses = compare(measured, doc)
    assert not misses, f"route {route} vs {path}:\n" + "\n".join(misses)


def test_wall_times_tsv_written(runs: dict[str, RouteRun], fetched: tuple[Path, str]) -> None:
    """Recorded, not asserted: only the file's shape is checked."""
    rows = _tsv(fetched[0] / "runs" / "wall_times.tsv")
    totals = {r["route"]: r for r in rows if r["step"] == "total"}
    assert sorted(totals) == list(ROUTES)
    for r in rows:
        float(r["wall_s"])
