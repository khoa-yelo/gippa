"""``gippa example fetch``: download, verify and unpack the example bundle (M12-01, ADR-0074).

The bundle pins are the ``example_<X>_<Y>`` rows of ``db_sources.yaml`` (``SOURCES``). The
archive is fetched by the ``gippa setup`` downloader (``HttpsFetcher``) into a staging directory
``DIR.partial-*`` next to ``DIR``, checked against the pin (size, then sha256), unpacked by the
setup's safe unpacker, checked file by file against its ``CHECKSUMS.sha256`` and renamed into
place. A failure leaves nothing behind. The printed commands come from ``example_routes.yaml``,
never from the bundle. A marker ``DIR/.gippa-example.json`` (not listed in ``CHECKSUMS.sha256``)
lets a re-fetch of the same bundle be an offline no-op.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from importlib.resources import files
from pathlib import Path, PurePosixPath
from typing import Any

import yaml

import gippa
from gippa.exitcodes import ExitCode
from gippa.setupdb import SOURCES, HttpsFetcher, SourceError
from gippa.setupdb.items import DbItem, file_digest, parse_minor
from gippa.setupdb.run import UNPACK_ERRORS, unpack_archive
from gippa.setupdb.space import DISK_FULL

__all__ = [
    "CHECKSUMS",
    "MARKER",
    "ExampleError",
    "content_problems",
    "fetch",
    "parse_version",
    "pins",
    "route_lines",
    "select_pin",
]

CHECKSUMS = "CHECKSUMS.sha256"
MARKER = ".gippa-example.json"
_PREFIX = "gippa example"
_INSTALLED = re.compile(r"(\d+)\.(\d+)")
SETUP_LINE = (
    "These routes need `gippa setup --no-genomad` (Bakta light and the antiSMASH databases) and a "
    "container runtime (Apptainer or Docker)."
)

Version = tuple[int, int]


class ExampleError(Exception):
    """A refusal or failure with its exit code and message lines (stderr)."""

    def __init__(self, code: ExitCode, *lines: str) -> None:
        super().__init__(lines[0] if lines else "")
        self.code = code
        self.lines = lines


@dataclass(frozen=True, slots=True)
class Selection:
    row: DbItem
    warning: str | None


# ---------------------------------------------------------------- pure parts


def parse_version(text: str) -> Version | None:
    """``X.Y`` or ``X.Y.Z`` (read as ``X.Y``); ``None`` otherwise."""
    return parse_minor(text)


def _vstr(v: Version) -> str:
    return f"{v[0]}.{v[1]}"


def pins(sources: Mapping[str, DbItem]) -> dict[Version, DbItem]:
    """The example bundle rows by version."""
    out: dict[Version, DbItem] = {}
    for name, row in sources.items():
        v = parse_version(row.version) if name.startswith("example_") else None
        if v is not None:
            out[v] = row
    return out


def select_pin(
    rows: Mapping[Version, DbItem], installed: str, requested: Version | None
) -> Selection:
    """Section 1 step 1: the requested pin, else the installed ``MAJOR.MINOR``, else the newest
    older minor of the same major (with the warning); else ``ExampleError``. A requested pin of
    another major is exit 2 (never a fallback, never a warning)."""
    m = _INSTALLED.match(installed)
    if m is None:
        raise ExampleError(ExitCode.ENVIRONMENT, f"cannot read the GIPPA version {installed!r}")
    inst = (int(m.group(1)), int(m.group(2)))
    if requested is not None:
        if requested not in rows:
            listed = ", ".join(_vstr(v) for v in sorted(rows)) or "none"
            raise ExampleError(ExitCode.INVALID_INPUT, f"no example bundle {_vstr(requested)}; "
                               f"pinned versions: {listed}")  # fmt: skip
        if requested[0] != inst[0]:
            raise ExampleError(ExitCode.ENVIRONMENT, f"example bundle {_vstr(requested)} is for "
                               f"GIPPA {requested[0]}.x; this is GIPPA {inst[0]} ({installed}): "
                               f"install GIPPA {requested[0]}.x to use it")  # fmt: skip
        chosen = requested
    else:
        same = [v for v in rows if v[0] == inst[0] and v[1] <= inst[1]]
        if not same:
            raise ExampleError(ExitCode.ENVIRONMENT, f"no example bundle for GIPPA {inst[0]}.x")
        chosen = max(same)
    warning = None
    if chosen < inst:
        c = _vstr(chosen)
        warning = (f"no bundle for {_vstr(inst)}; using {c} (expected numbers were made with "
                   f"GIPPA {c} and may differ, see Upgrade)")  # fmt: skip
    return Selection(rows[chosen], warning)


def _routes() -> list[dict[str, Any]]:
    raw = yaml.safe_load(files("gippa.data").joinpath("example_routes.yaml").read_text("utf-8"))
    routes: list[dict[str, Any]] = raw["routes"]
    return routes


def route_lines(version: str, directory: Path, routes: Sequence[Mapping[str, Any]]) -> list[str]:
    """The stdout block: the bundle line, each route's heading and commands, the setup line."""
    d = str(directory)
    out = [f"{_PREFIX} {version}: {directory}"]
    for route in routes:
        out += ["", f"# {route['title']}"]
        for cmd in route["commands"]:
            words = (w.replace("{out}", "{dir}/runs").replace("{dir}", d)
                     for w in shlex.split(str(cmd)))  # fmt: skip
            out.append(" ".join(shlex.quote(w) for w in words))
    return [*out, "", SETUP_LINE]


def _listed(text: str) -> dict[str, str] | str:
    """``CHECKSUMS.sha256`` lines -> ``{path: hex}``, or the first malformed line."""
    out: dict[str, str] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        m = re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", line)
        name = m.group(2) if m else ""
        p = PurePosixPath(name)
        if m is None or p.is_absolute() or ".." in p.parts or name in out:
            return line
        out[name] = m.group(1)
    return out


def content_problems(
    root: Path, ignore: frozenset[str] = frozenset({CHECKSUMS}), *, unlisted: bool = True
) -> list[str]:
    """Every ``CHECKSUMS.sha256`` line matches its regular file, and (when ``unlisted``) no file
    outside ``ignore`` is unlisted; the problems, each naming its file."""
    try:
        listed = _listed((root / CHECKSUMS).read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError):
        return [f"no readable {CHECKSUMS}"]
    if isinstance(listed, str):
        return [f"malformed {CHECKSUMS} line {listed!r}"]
    problems = []
    for name, want in sorted(listed.items()):
        f = root / name
        if f.is_symlink() or not f.is_file():
            problems.append(f"missing file {name}")
        elif file_digest(f, "sha256") != f"sha256:{want}":
            problems.append(f"checksum mismatch for {name}")
    for f in sorted(root.rglob("*")) if unlisted else ():
        rel = f.relative_to(root).as_posix()
        if (f.is_symlink() or not f.is_dir()) and rel not in listed and rel not in ignore:
            problems.append(f"unlisted file {rel}")
    return problems


# ---------------------------------------------------------------- the command


def _say(*lines: str) -> None:
    for line in lines:
        print(f"{_PREFIX}: {line}", file=sys.stderr, flush=True)


def _marker(row: DbItem, installed: str) -> dict[str, str]:
    return {"archive_sha256": row.checksum, "bundle_version": row.version,
            "gippa_version": installed}  # fmt: skip


def _already_fetched(directory: Path, row: DbItem) -> bool:
    """``DIR`` holds a fetch of ``row``: its marker names the pin and every listed file verifies.
    Unlisted entries (the routes' ``DIR/runs/``) are the user's and ignored."""
    try:
        mark = json.loads((directory / MARKER).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return False
    same = isinstance(mark, dict) and (mark.get("archive_sha256"), mark.get("bundle_version")) == (
        row.checksum, row.version)  # fmt: skip
    return same and not content_problems(directory, frozenset({CHECKSUMS, MARKER}), unlisted=False)


def _first_missing(path: Path) -> Path:
    while not os.path.lexists(path.parent) and path.parent != path:
        path = path.parent
    return path


def _disk_full(directory: Path) -> ExampleError:
    return ExampleError(ExitCode.ENVIRONMENT, f"disk quota exceeded or disk full at "
                        f"{directory.parent}: choose a larger disk with --outdir")  # fmt: skip


def _not_writable(where: Path, exc: OSError | None = None) -> ExampleError:
    why = f" ({exc.strerror or exc})" if exc is not None else ""
    return ExampleError(ExitCode.INVALID_INPUT, f"cannot write to {where}{why}; choose another "
                        "--outdir")  # fmt: skip


def _writable(d: Path) -> bool:
    return os.access(d, os.W_OK | os.X_OK)


def _check_archive(
    archive: Path, row: DbItem, rows: Mapping[Version, DbItem], offline: bool
) -> None:
    """Size, then sha256, against the pin. A wrong ``--from-file`` is exit 1 (and names the other
    pin it matches), a wrong download exit 2."""
    size = archive.stat().st_size
    bad_size = bool(row.size_bytes) and size != row.size_bytes
    got = "" if bad_size and not offline else file_digest(archive, "sha256")
    if not bad_size and got == row.checksum:
        return
    if offline:
        other = next((r for r in rows.values() if r.checksum == got and r is not row), None)
        if other is not None:
            raise ExampleError(ExitCode.INVALID_INPUT, f"{archive} is not example bundle "
                               f"{row.version}: this is example bundle {other.version}; use "
                               f"--version {other.version}")  # fmt: skip
    code = ExitCode.INVALID_INPUT if offline else ExitCode.ENVIRONMENT
    if bad_size:
        raise ExampleError(code, f"size mismatch: {archive.name} has {size} bytes, pinned "
                           f"{row.size_bytes}")  # fmt: skip
    where = archive.name if offline else row.host
    raise ExampleError(code, f"checksum mismatch from {where}: {got}, pinned {row.checksum}")


def _download(row: DbItem, dest: Path, directory: Path) -> None:
    _say(f"downloading {row.url}")
    try:
        HttpsFetcher().fetch(row.url, dest)
    except SourceError as exc:
        lines = [f"cannot download {row.url} ({exc.host}: {exc.reason})",
                 "on a machine with internet access, download that URL, copy the file here and "
                 "run\n  gippa example fetch --from-file "
                 f"{shlex.quote(PurePosixPath(row.path).name)} "
                 f"--outdir {shlex.quote(str(directory))}"]  # fmt: skip
        if "HTTP 404" in exc.reason:
            lines.append(f"the bundle may not be published yet (GIPPA {row.version} release)")
        raise ExampleError(ExitCode.ENVIRONMENT, *lines) from exc


def _unpack(archive: Path, tree: Path, version: str) -> Path:
    """Safe unpack; exactly one top-level directory ``gippa-example-<X.Y>``."""
    tree.mkdir()
    try:
        unsafe = unpack_archive(archive, tree)
    except UNPACK_ERRORS as exc:  # a malformed archive is the source's fault; others are bugs
        raise ExampleError(
            ExitCode.ENVIRONMENT, f"malformed archive {archive.name}: {exc}"
        ) from exc
    except OSError as exc:
        if exc.errno in DISK_FULL:
            raise
        raise ExampleError(ExitCode.ENVIRONMENT, f"cannot unpack {archive.name}: {exc}") from exc
    if unsafe is not None:
        raise ExampleError(ExitCode.ENVIRONMENT, unsafe)
    top = f"gippa-example-{version}"
    entries = sorted(p.name for p in tree.iterdir())
    if entries != [top] or (tree / top).is_symlink() or not (tree / top).is_dir():
        raise ExampleError(ExitCode.ENVIRONMENT, f"malformed archive {archive.name}: expected one "
                           f"top-level directory {top}/, found {entries}")  # fmt: skip
    return tree / top


def _install(
    row: DbItem, rows: Mapping[Version, DbItem], directory: Path, from_file: Path | None,
    installed: str,
) -> None:  # fmt: skip
    """Steps 3-7 inside a staging directory next to ``directory``; removed whatever happens."""
    directory.parent.mkdir(parents=True, exist_ok=True)
    try:
        staging = Path(tempfile.mkdtemp(prefix=f"{directory.name}.partial-", dir=directory.parent))
    except PermissionError as exc:
        raise _not_writable(directory.parent, exc) from exc
    try:
        if from_file is None:
            archive = staging / PurePosixPath(row.path).name
            _download(row, archive, directory)
        else:
            archive = from_file
        _check_archive(archive, row, rows, from_file is not None)
        top = _unpack(archive, staging / "tree", row.version)
        if problems := content_problems(top):
            raise ExampleError(ExitCode.ENVIRONMENT, f"bundle content does not verify against "
                               f"{CHECKSUMS}: " + "; ".join(problems))  # fmt: skip
        (top / MARKER).write_text(json.dumps(_marker(row, installed), sort_keys=True) + "\n",
                                  encoding="utf-8")  # fmt: skip
        try:
            if directory.is_dir():
                for entry in sorted(top.iterdir()):
                    entry.rename(directory / entry.name)
            else:
                top.rename(directory)
        except PermissionError as exc:
            raise _not_writable(directory, exc) from exc
    except OSError as exc:
        if exc.errno in DISK_FULL:
            raise _disk_full(directory) from exc
        raise
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def fetch(outdir: Path, version: Version | None, from_file: Path | None) -> int:
    """``gippa example fetch``; returns the exit code, messages on stderr, commands on stdout."""
    directory = outdir.absolute()
    installed = str(gippa.__version__)
    created: Path | None = None
    existed_empty = False
    try:
        rows = pins(SOURCES)
        sel = select_pin(rows, installed, version)
        if sel.warning:
            _say(sel.warning)
        if os.path.lexists(directory):
            if (
                directory.is_dir()
                and not directory.is_symlink()
                and _already_fetched(directory, sel.row)
            ):
                _say(f"{directory} already fetched (example bundle {sel.row.version})")
                print("\n".join(route_lines(sel.row.version, directory, _routes())), flush=True)
                return ExitCode.SUCCESS
            if not directory.is_dir() or any(directory.iterdir()):
                raise ExampleError(ExitCode.INVALID_INPUT, f"{directory} exists and is not this "
                                   "example; use --outdir or remove it")  # fmt: skip
            if not _writable(directory):
                raise _not_writable(directory)
            existed_empty = True
        else:
            created = _first_missing(directory)
            try:
                directory.parent.mkdir(parents=True, exist_ok=True)
            except OSError as exc:
                raise ExampleError(
                    ExitCode.INVALID_INPUT,
                    f"cannot create {directory} ({exc.strerror or exc}); choose another --outdir",
                ) from exc
            if not _writable(directory.parent):
                raise _not_writable(directory.parent)
        _install(sel.row, rows, directory, from_file, installed)
    except BaseException as exc:
        if created is not None:
            shutil.rmtree(created, ignore_errors=True)
        if existed_empty:  # leave the user's empty directory empty
            for entry in directory.iterdir():
                if entry.is_dir() and not entry.is_symlink():
                    shutil.rmtree(entry)
                else:
                    entry.unlink()
        if isinstance(exc, ExampleError):
            _say(*exc.lines)
            return exc.code
        raise
    print("\n".join(route_lines(sel.row.version, directory, _routes())), flush=True)
    return ExitCode.SUCCESS
