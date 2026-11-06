"""Shared fixtures of the M12-01 ``gippa example fetch`` tests
(docs/specs/M12-01-example-dataset.md).

Tiny synthetic bundles built in ``tmp_path`` (never the real one: R11), pin rows monkeypatched into
``SOURCES``, fake fetchers in place of ``HttpsFetcher``, the installed version through
``gippa.__version__``, and an in-process ``gippa`` call through ``gippa.cli.main`` (so the §5.1 exit
code mapping of ``main`` applies: Click usage errors 1, uncaught exceptions 5).

Seams the implementation must honour (the spec's "the test monkeypatches the selected row"):

- the pins are the ``example_*`` rows of ``gippa.setupdb.items.SOURCES``, read at call time (the
  tests replace ``SOURCES`` in ``gippa.setupdb.items``, ``gippa.setupdb.run``, ``gippa.setupdb``
  and ``gippa.example``);
- the installed version is ``gippa.__version__``, read at call time (also replaced in
  ``gippa.example`` when it holds its own copy);
- the downloader is ``HttpsFetcher()`` from ``gippa.setupdb`` (replaced in ``gippa.setupdb.run``,
  ``gippa.setupdb`` and ``gippa.example``).
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib
import io
import json
import re
import shlex
import tarfile
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
ROUTES_YAML = REPO / "src" / "gippa" / "data" / "example_routes.yaml"
CLI_MD = REPO / "docs" / "reference" / "cli.md"
EXAMPLE_PY = REPO / "src" / "gippa" / "example.py"
ASSET_HOST = "release-assets.githubusercontent.com"
URL_TEMPLATE = (
    "https://github.com/khoa-yelo/gippa/releases/download/v{v}.0/gippa-example-{v}.tar.gz"
)
EXAMPLE_SETTINGS = (
    "--modules core,bgc",
    "--set annotation.bakta.db=light",
    "--terminators transtermhp",
    "--set annotation.promoter_calculator.enabled=false",
    "--threads 4",
)
BOGUS_BUNDLE_COMMAND = "gippa bogus-command-from-the-bundle-readme"
_MODULES = ("gippa.setupdb.items", "gippa.setupdb.run", "gippa.setupdb", "gippa.example")


def url_of(version: str) -> str:
    return URL_TEMPLATE.format(v=version)


def sha256_hex(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


# ---------------------------------------------------------------- bundles


def bundle_files(version: str = "1.1") -> dict[str, bytes]:
    """The layout of spec section 2 with tiny placeholder contents (no row-level data)."""
    meta = {"bundle_version": version, "gippa_version": f"{version}.0", "region": "synthetic"}
    out = {
        "README.md": f"# example {version}\n\n```\n{BOGUS_BUNDLE_COMMAND}\n```\n".encode(),
        "PROVENANCE.md": b"# provenance (test fixture)\n",
        "LICENSE-DATA.md": b"# data licence (test fixture)\n",
        "bundle.json": (json.dumps(meta, sort_keys=True) + "\n").encode(),
        "SELECTION.json": b'{"rule": "test fixture"}\n',
        "genome/region.fna.gz": gzip.compress(b">region\nACGTACGTAC\n", mtime=0),
        "A_diann/report.parquet": b"PAR1 not a real parquet PAR1",
        "B_proteins/proteins.faa": b">p1\nMK\n",
        "B_proteins/abundance.tsv": b"protein_id\tabundance\np1\t1.0\n",
        "C_cascadia/WT_1.ssl": b"file\tscan\n",
        "C_cascadia/WT_2.ssl": b"file\tscan\n",
        "C_cascadia/WT_3.ssl": b"file\tscan\n",
        "D_generic/peptides.tsv": b"peptide\trun\n",
        "E_sixframe/report.parquet": b"PAR1 not a real parquet PAR1",
        "expected/README.md": b"# expected (test fixture)\n",
    }
    for r in "ABCDE":
        out[f"expected/{r}.json"] = (json.dumps({"route": r, "version": version}) + "\n").encode()
    return out


def checksums_text(files: Mapping[str, bytes]) -> bytes:
    """``sha256sum`` format, sorted by path, of every file except ``CHECKSUMS.sha256``."""
    lines = [f"{sha256_hex(data)}  {name}\n" for name, data in sorted(files.items())
             if name != "CHECKSUMS.sha256"]  # fmt: skip
    return "".join(lines).encode()


@dataclass(frozen=True)
class Member:
    """One tar member: a file (``data``), a directory or a symlink (``link``)."""

    name: str
    data: bytes | None = None
    kind: str = "file"  # file | dir | sym
    link: str = ""


def write_tar(dest: Path, members: Iterable[Member]) -> Path:
    """A deterministic tar.gz (sorted members, uid/gid 0, mtime 0, no gzip name or mtime)."""
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for m in sorted(members, key=lambda m: m.name):
            info = tarfile.TarInfo(m.name)
            info.mtime, info.uid, info.gid, info.uname, info.gname = 0, 0, 0, "", ""
            if m.kind == "dir":
                info.type, info.mode = tarfile.DIRTYPE, 0o755
                tar.addfile(info)
            elif m.kind == "sym":
                info.type, info.linkname, info.mode = tarfile.SYMTYPE, m.link, 0o777
                tar.addfile(info)
            else:
                data = m.data or b""
                info.size, info.mode = len(data), 0o644
                tar.addfile(info, io.BytesIO(data))
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("wb") as raw, gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as gz:
        gz.write(buf.getvalue())
    return dest


def tree_members(top: str, files: Mapping[str, bytes]) -> list[Member]:
    dirs = {top} | {f"{top}/{p}" for name in files for p in _parents(name)}
    return [Member(d, kind="dir") for d in sorted(dirs)] + [
        Member(f"{top}/{name}", data) for name, data in files.items()
    ]


def _parents(name: str) -> list[str]:
    parts = name.split("/")[:-1]
    return ["/".join(parts[: i + 1]) for i in range(len(parts))]


@dataclass(frozen=True)
class Bundle:
    version: str
    path: Path
    sha256: str  # "sha256:<hex>"
    size: int
    files: Mapping[str, bytes] = field(default_factory=dict)  # as listed in CHECKSUMS.sha256

    @property
    def url(self) -> str:
        return url_of(self.version)

    @property
    def data(self) -> bytes:
        return self.path.read_bytes()


def _finish(version: str, path: Path, files: Mapping[str, bytes]) -> Bundle:
    data = path.read_bytes()
    return Bundle(version, path, f"sha256:{sha256_hex(data)}", len(data), dict(files))


def build_bundle(
    where: Path,
    version: str = "1.1",
    *,
    alter: Mapping[str, bytes] | None = None,
    extra: Mapping[str, bytes] | None = None,
    drop: Iterable[str] = (),
    members: Iterable[Member] = (),
    top: str | None = None,
) -> Bundle:
    """``gippa-example-<version>.tar.gz`` in ``where``. ``alter`` replaces file contents after
    ``CHECKSUMS.sha256`` was written; ``extra`` adds unlisted files; ``drop`` removes listed files;
    ``members`` adds raw members (unsafe ones, a second top-level entry); ``top`` renames the
    top-level directory."""
    files = bundle_files(version)
    listed = {**files, "CHECKSUMS.sha256": checksums_text(files)}
    shipped = {**listed, **(alter or {}), **(extra or {})}
    for name in drop:
        shipped.pop(name)
    tree = tree_members(top or f"gippa-example-{version}", shipped)
    path = write_tar(where / f"gippa-example-{version}.tar.gz", [*tree, *members])
    return _finish(version, path, listed)


# ---------------------------------------------------------------- pins, version, fetcher


def pin_row(
    bundle: Bundle | None = None, *, version: str | None = None, sha256: str = "", size: int = 0
) -> Any:
    """A ``DbItem`` row ``example_<X>_<Y>`` as spec section 1 > Pin shows it."""
    from gippa.setupdb.items import DbItem

    v = bundle.version if bundle is not None else str(version)
    major, minor = v.split(".")
    return DbItem(
        name=f"example_{major}_{minor}",
        version=v,
        url=url_of(v),
        host="github.com",
        checksum=bundle.sha256 if bundle is not None else sha256 or "sha256:" + "0" * 64,
        path=f"example/gippa-example-{v}.tar.gz",
        tool=None,
        check="sha256",
        size_gib=1,
        redirect_hosts=(ASSET_HOST,),
        size_bytes=bundle.size if bundle is not None else size or 1,  # type: ignore[call-arg]
    )


def _loaded(name: str) -> Any | None:
    try:
        return importlib.import_module(name)
    except ImportError:
        return None


def _patch_all(mp: pytest.MonkeyPatch, attr: str, value: Any, modules: Iterable[str]) -> None:
    for name in modules:
        mod = _loaded(name)
        if mod is not None and (hasattr(mod, attr) or name != "gippa.example"):
            mp.setattr(mod, attr, value, raising=False)


def install_pins(mp: pytest.MonkeyPatch, rows: Iterable[Any]) -> None:
    """``SOURCES`` = the real non-example rows plus ``rows`` (any real ``example_*`` row is
    dropped, so the tests do not depend on the pin PR)."""
    from gippa.setupdb import items

    base = {k: v for k, v in items.SOURCES.items() if not k.startswith("example_")}
    new = MappingProxyType({**base, **{r.name: r for r in rows}})
    _patch_all(mp, "SOURCES", new, _MODULES)


def set_installed(mp: pytest.MonkeyPatch, version: str) -> None:
    import gippa

    mp.setattr(gippa, "__version__", version)
    _patch_all(mp, "__version__", version, ("gippa.example",))


class FetchLog:
    def __init__(self) -> None:
        self.constructed = 0
        self.calls: list[tuple[str, Path]] = []


def install_fetcher(
    mp: pytest.MonkeyPatch, fetch: Callable[[str, Path], None] | None, log: FetchLog | None = None
) -> FetchLog:
    """Replace ``HttpsFetcher`` with a fake whose ``fetch`` calls ``fetch(url, dest)``; with
    ``fetch=None`` constructing the fetcher raises (``--from-file`` and "already fetched")."""
    log = log or FetchLog()

    class FakeFetcher:
        def __init__(self, *_a: Any, **_k: Any) -> None:
            log.constructed += 1
            if fetch is None:
                raise AssertionError("HttpsFetcher constructed: this path must not use the network")

        def fetch(self, url: str, dest: Path) -> None:
            log.calls.append((url, Path(dest)))
            assert fetch is not None
            fetch(url, Path(dest))

    _patch_all(mp, "HttpsFetcher", FakeFetcher, ("gippa.setupdb.run", "gippa.setupdb",
                                                 "gippa.example"))  # fmt: skip
    return log


def serve(payload: bytes) -> Callable[[str, Path], None]:
    def fetch(_url: str, dest: Path) -> None:
        dest.write_bytes(payload)

    return fetch


def raise_(exc: BaseException) -> Callable[[str, Path], None]:
    def fetch(_url: str, dest: Path) -> None:
        raise exc

    return fetch


def block_https(mp: pytest.MonkeyPatch) -> list[str]:
    """``http.client.HTTPSConnection`` raises (spec section 5); returns the hosts asked for."""
    import http.client

    hosts: list[str] = []

    def deny(host: str, *_a: Any, **_k: Any) -> None:
        hosts.append(host)
        raise AssertionError(f"http.client.HTTPSConnection({host!r}) in a no-network test")

    mp.setattr(http.client, "HTTPSConnection", deny)
    return hosts


# ---------------------------------------------------------------- running gippa


@dataclass(frozen=True)
class Result:
    code: int
    out: str
    err: str


def gippa(argv: list[str]) -> Result:
    """``gippa <argv>`` in-process through ``gippa.cli.main`` (its exit code mapping)."""
    from gippa.cli import main

    out, err = io.StringIO(), io.StringIO()
    code: Any = 0
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            main(argv)
        except SystemExit as exc:
            code = exc.code
    if code is None:
        code = 0
    if not isinstance(code, int):  # sys.exit("message")
        err.write(str(code))
        code = 1
    return Result(int(code), out.getvalue(), err.getvalue())


def listing(d: Path) -> list[str]:
    """Every path under ``d`` (relative, sorted): the "nothing left behind" snapshot."""
    return sorted(p.relative_to(d).as_posix() for p in d.rglob("*"))


def verify_tree(d: Path, files: Mapping[str, bytes]) -> None:
    """Every ``CHECKSUMS.sha256`` line verifies, and every listed file has the bundle's bytes."""
    lines = (d / "CHECKSUMS.sha256").read_text(encoding="utf-8").splitlines()
    assert lines, "empty CHECKSUMS.sha256"
    for line in lines:
        digest, _, name = line.partition("  ")
        assert sha256_hex((d / name).read_bytes()) == digest, name
    for name, data in files.items():
        assert (d / name).read_bytes() == data, name


def flat(text: str) -> str:
    """Whitespace collapsed and shell line continuations joined."""
    return " ".join(text.replace("\\\n", " ").split())


# ---------------------------------------------------------------- example_routes.yaml


def load_routes(path: Path = ROUTES_YAML) -> dict[str, dict[str, Any]]:
    """``{route_id: {"title": str, "commands": [str, ...]}}`` from ``example_routes.yaml``.

    The spec does not fix the YAML shape; accepted: a top-level ``routes`` key holding a mapping
    ``id -> route`` or a list of routes with an ``id``; or the ``id -> route`` mapping itself. A
    route holds ``commands`` (a list of strings) and a heading under ``title``, ``heading`` or
    ``name``.
    """
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    node = raw.get("routes", raw) if isinstance(raw, dict) else raw
    if isinstance(node, list):
        node = {str(r["id"]): r for r in node}
    out: dict[str, dict[str, Any]] = {}
    for rid, route in node.items():
        title = next((route[k] for k in ("title", "heading", "name") if k in route), "")
        cmds = route["commands"]
        assert isinstance(cmds, list) and all(isinstance(c, str) for c in cmds), rid
        out[str(rid)] = {"title": str(title), "commands": [flat(c) for c in cmds]}
    return out


def fill(cmd: str, d: str) -> str:
    """A route command with ``{out}`` -> ``{dir}/runs`` and ``{dir}`` -> ``d``."""
    return cmd.replace("{out}", "{dir}/runs").replace("{dir}", d)


def subcommand(cmd: str) -> str:
    return shlex.split(cmd)[1]


def options(cmd: str) -> dict[str, list[str]]:
    """``--option`` -> its values (``""`` for a flag) of one command line."""
    words = shlex.split(cmd)
    out: dict[str, list[str]] = {}
    for i, w in enumerate(words):
        if w.startswith("--"):
            nxt = words[i + 1] if i + 1 < len(words) else ""
            out.setdefault(w, []).append("" if nxt.startswith("--") else nxt)
    return out


# ---------------------------------------------------------------- the CLI reference


def cli_reference(text: str) -> dict[tuple[str, ...], set[str]]:
    """``docs/reference/cli.md`` sections ``## gippa <words>`` -> their documented options; the
    command path may have several words (``## gippa example fetch``; Q-L4)."""
    out: dict[tuple[str, ...], set[str]] = {}
    for block in re.split(r"\n## ", "\n" + text)[1:]:
        head, _, body = block.partition("\n")
        words = tuple(head.split())
        if not words or words[0] != "gippa":
            continue
        rows = [ln.split("|")[1] for ln in body.splitlines() if ln.startswith("| `")]
        out[words[1:]] = {o for r in rows for o in re.findall(r"`(--[a-z0-9-]+)`", r)} | {"--help"}
    return out


def option_problems(cmdline: str, ref: Mapping[tuple[str, ...], set[str]]) -> list[str]:
    """Options of ``gippa <cmd ...> --opt`` that the longest matching section does not list."""
    words = cmdline.split()
    if len(words) < 2 or words[0] != "gippa":
        return []
    path = next((tuple(words[1:n]) for n in range(len(words), 1, -1)
                 if tuple(words[1:n]) in ref), None)  # fmt: skip
    if path is None:
        return [f"no section for {' '.join(words[:3])!r}"]
    opts = re.findall(r"(?<![\w-])--[a-z0-9-]+", cmdline)
    return [f"{' '.join(path)}: {o}" for o in opts if o not in ref[path]]
