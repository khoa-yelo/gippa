"""M12-01 PR 3/PR 7: ``scripts/example/build_example.py`` writes the release bundle.

docs/specs/M12-01-example-dataset.md section 2 (Layout, Subset rules 8, Licence statement) and
section 3 (How the bundle is built), ADR-0074 decisions 1, 4, 5. Pinned here:

1. Defaults: ``--name nissle_r200k``, region ``NZ_CP058217.1:1965001-2165000``
   (``[1965000, 2165000)``), bundle version = the installed GIPPA ``MAJOR.MINOR``, bundle at
   ``$GIPPA_ROOT/example/gippa-example-<X.Y>.tar.gz`` (``--bundle-out DIR`` moves it).
2. Release layout: one top dir ``gippa-example-<X.Y>/``; every ``{dir}`` path of
   ``example_routes.yaml``; ``README.md``, ``PROVENANCE.md``, ``LICENSE-DATA.md``,
   ``bundle.json``, ``SELECTION.json`` and ``CHECKSUMS.sha256`` (every other file, sorted); nothing
   outside the layout (no ``bundles/`` sizing tars); no personal path (R12); and
   ``gippa example fetch --from-file`` accepts it with a pin made from the printout.
3. Determinism: two builds (different work dirs) give byte-identical archives; tar and gzip
   headers carry no time, owner or name.
4. The printout: the archive sha256 and size, a ``db_sources.yaml`` row the loader accepts, and a
   ``validation/example/bundles.tsv`` row.

The pilot is a tiny synthetic stand-in under ``tmp_path`` (``$GIPPA_ROOT/pilot/...``, the paths
``build_example.sources`` reads); the real pilot data is never read. No statistical claim is made
by the build, so there is no calibration test.
"""

from __future__ import annotations

import contextlib
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import random
import re
import tarfile
import textwrap
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from types import ModuleType
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
import yaml
from helpers.example_m12_01 import (
    ASSET_HOST,
    REPO,
    block_https,
    fill,
    flat,
    gippa,
    install_pins,
    load_routes,
    pin_row,
    url_of,
)

SCRIPT = REPO / "scripts" / "example" / "build_example.py"
SOURCES_YAML = REPO / "src" / "gippa" / "data" / "db_sources.yaml"
BUNDLES_TSV = REPO / "validation" / "example" / "bundles.tsv"
CHROM = "NZ_CP058217.1"
START, END = 1_965_000, 2_165_000
NAME = "nissle_r200k"
RUNS = ("WT_1", "WT_2", "WT_3")
PERSONAL = "/home/users/alice/raw"  # an absolute path in the pilot's files (rule 8, R12)
COLUMNS = ["version", "file", "size_bytes", "sha256", "url", "release_tag", "build_commit",
           "build_job", "expected_commit", "asset_digest"]  # fmt: skip
TOP_FILES = {"README.md", "PROVENANCE.md", "LICENSE-DATA.md", "CHECKSUMS.sha256", "bundle.json",
             "SELECTION.json"}  # fmt: skip
HEX64 = re.compile(r"[0-9a-f]{64}")


def _minor() -> str:
    import gippa

    m = re.match(r"(\d+)\.(\d+)", str(gippa.__version__))
    assert m is not None, gippa.__version__
    return f"{m.group(1)}.{m.group(2)}"


def _load_script() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gippa_build_example_m12_01", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- the synthetic pilot

_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODE = {a + b + c: _AA[16 * i + 4 * j + k] for i, a in enumerate("TCAG")
         for j, b in enumerate("TCAG") for k, c in enumerate("TCAG")}  # fmt: skip


def _translate(nt: str) -> str:
    return "".join(_CODE[nt[i : i + 3]] for i in range(0, len(nt) - 2, 3))


def _fasta_gz(path: Path, records: list[tuple[str, str]]) -> None:
    text = "".join(f">{n}\n" + "\n".join(s[i : i + 60] for i in range(0, len(s), 60)) + "\n"
                   for n, s in records)  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode(), mtime=0))


def _stop_free(prot: str, n: int, k: int) -> list[str]:
    """``k`` disjoint stop-free windows of ``n`` residues of ``prot``."""
    out, i = [], 0
    while len(out) < k and i + n <= len(prot):
        w = prot[i : i + n]
        if "*" in w:
            i += w.rindex("*") + 1
        else:
            out.append(w)
            i += n + 50
    assert len(out) == k
    return out


def make_pilot(root: Path) -> None:
    """``$GIPPA_ROOT/pilot/pass3`` and ``A3_cascadia_cpu`` with the files the script reads."""
    rng = random.Random(1917)
    p3, a3 = root / "pilot" / "pass3", root / "pilot" / "A3_cascadia_cpu"
    chrom = "".join(rng.choices("ACGT", k=END + 5_000))
    plasmid = "".join(rng.choices("ACGT", k=3_000))
    _fasta_gz(p3 / "genome" / "nissle1917_wt.fna.gz", [(CHROM, chrom), ("NZ_CP058218.1", plasmid)])

    # region proteins (1-based closed, both strands, CDS and sORF) and two outside the region
    prot = [(f"NIS_{i:04d}", START + 1 + 900 * i, START + 900 * i + 600, "+-"[i % 2],
             "sorf" if i == 5 else "cds") for i in range(1, 7)]  # fmt: skip
    prot += [("NIS_OUT1", 1_000, 1_600, "+", "cds"), ("NIS_EDGE", END - 300, END + 300, "-", "cds")]
    (p3 / "db_a2").mkdir(parents=True, exist_ok=True)
    with (p3 / "db_a2" / "db_id_map.tsv").open("w") as fh:
        fh.write("db_id\tcontig\tstart\tend\tstrand\tkind\n")
        for pid, s, e, strand, kind in prot:
            fh.write(f"{pid}\t{CHROM}\t{s}\t{e}\t{strand}\t{kind}\n")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    _fasta_gz(p3 / "db_a2" / "strain_db.faa.gz",
              [(p[0], "M" + "".join(rng.choices(aa, k=120))) for p in prot])  # fmt: skip

    # DIA-NN 2.3.1 report rows: targets, decoys (same IDs, Decoy = 1), a mixed and a
    # contaminant row (both dropped); File.Name holds an absolute path (rule 8)
    spec = [("NIS_0001", "AAAPEPTLDEK", "AAAPEPTLDEK", 0),
            ("NIS_0002", "LLLMPEPTIDER", "LLLM(UniMod:35)PEPTIDER", 0),
            ("NIS_0003;NIS_0004", "SHAREDPEPK", "SHAREDPEPK", 0),
            ("NIS_0005", "SMALLORFPEPK", "SMALLORFPEPK", 0),
            ("NIS_0001", "KEDLTPEPAAA", "KEDLTPEPAAA", 1),
            ("NIS_0002", "REDITPEPMLLL", "REDITPEPMLLL", 1),
            ("NIS_0001;NIS_OUT1", "MIXEDPEPK", "MIXEDPEPK", 0),
            ("Cont_ALBU_BOVIN", "CONTAMPEPK", "CONTAMPEPK", 0)]  # fmt: skip
    rows: list[dict[str, Any]] = []
    for r, run in enumerate(RUNS):
        for j, (ids, pep, mod, decoy) in enumerate(spec):
            rows.append({
                "File.Name": f"{PERSONAL}/{run}.mzML", "Run": run, "Protein.Group": ids,
                "Protein.Ids": ids, "Genes": "", "Stripped.Sequence": pep,
                "Modified.Sequence": mod, "Precursor.Id": f"{mod}2", "Precursor.Charge": 2,
                "Precursor.Normalised": 1000.0 * (j + 1) + 10 * r, "Q.Value": 0.001 * (j + 1),
                "Global.Q.Value": 0.0005 * (j + 1), "PEP": 0.002 * (j + 1),
                "PG.MaxLFQ": 5000.0 + 100 * j + r, "Decoy": decoy,
            })  # fmt: skip
    (p3 / "A2").mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(pd.DataFrame(rows), preserve_index=False),
                   p3 / "A2" / "A2_report.parquet")  # fmt: skip

    # Cascadia rows: three peptides of the region's frame 0 (kept), one random (dropped)
    peps = [*_stop_free(_translate(chrom[START:END]), 12, 3), "WWWWHHHHWWWWHHHH"]
    for run in RUNS:
        d = a3 / run
        d.mkdir(parents=True, exist_ok=True)
        with (d / f"{run}.ssl").open("w") as fh:
            fh.write("file\tscan\tcharge\tsequence\tscore-type\tscore\n")
            for i, pep in enumerate(peps):
                fh.write(f"{PERSONAL}/{run}.mzML\t{100 + i}\t2\t{pep}\tUNKNOWN\t0.{9 - i}\n")


# ---------------------------------------------------------------- running the build


@dataclass(frozen=True)
class Build:
    root: Path
    archive: Path
    data: bytes
    stdout: str
    tree: Path  # the archive unpacked by the setup's safe unpacker

    @property
    def top(self) -> Path:
        return self.tree / f"gippa-example-{_minor()}"


def _run(root: Path, argv: list[str]) -> tuple[int, str]:
    """``build_example.main(argv)`` with ``GIPPA_ROOT=root``: (exit code, stdout + stderr)."""
    mod = _load_script()
    out = io.StringIO()
    with (
        pytest.MonkeyPatch.context() as mp,
        contextlib.redirect_stdout(out),
        contextlib.redirect_stderr(out),
    ):
        mp.setenv("GIPPA_ROOT", str(root))
        try:
            code = mod.main(argv)
        except SystemExit as exc:  # argparse
            code = exc.code if isinstance(exc.code, int) else 1
    return int(code or 0), out.getvalue()


def _unpack(archive: Path, tree: Path) -> None:
    from gippa.setupdb.run import unpack_archive

    tree.mkdir(parents=True)
    assert unpack_archive(archive, tree) is None, "unsafe member in the bundle"


@pytest.fixture(scope="module")
def pilot(tmp_path_factory: pytest.TempPathFactory) -> Path:
    root = tmp_path_factory.mktemp("gippa_root").resolve()
    make_pilot(root)
    return root


@pytest.fixture(scope="module")
def default_build(pilot: Path, tmp_path_factory: pytest.TempPathFactory) -> Build:
    """``build_example.py`` with no argument at all (every default)."""
    code, out = _run(pilot, [])
    assert code == 0, out
    archive = pilot / "example" / f"gippa-example-{_minor()}.tar.gz"
    assert archive.is_file(), f"no bundle at {archive}"
    tree = tmp_path_factory.mktemp("unpacked") / "tree"
    _unpack(archive, tree)
    return Build(pilot, archive, archive.read_bytes(), out, tree)


def _files(top: Path) -> dict[str, Path]:
    return {p.relative_to(top).as_posix(): p for p in sorted(top.rglob("*")) if p.is_file()}


def _route_paths() -> set[str]:
    """Every ``{dir}/<path>`` that ``example_routes.yaml`` names."""
    out = set()
    for route in load_routes().values():
        for cmd in route["commands"]:
            out |= set(re.findall(r"\{dir\}/(\S+)", cmd))
    assert out, "example_routes.yaml names no {dir} path"
    return out


def _readable(path: Path) -> str:
    """A bundle file as text: parquet cells, gunzipped FASTA, or the bytes."""
    if path.suffix == ".parquet":
        return json.dumps(pq.read_table(path).to_pylist(), default=str)
    raw = path.read_bytes()
    return (gzip.decompress(raw) if path.suffix == ".gz" else raw).decode("utf-8", "replace")


# ---------------------------------------------------------------- 1. defaults


def test_script_needs_no_required_argument(default_build: Build) -> None:
    """``--name``, ``--start`` and ``--end`` default (no argparse error, exit 0)."""
    assert default_build.archive.name == f"gippa-example-{_minor()}.tar.gz"


def test_default_region_and_name(default_build: Build) -> None:
    sel = json.loads((default_build.top / "SELECTION.json").read_text())
    assert sel["name"] == NAME
    reg = sel["region"]
    assert (reg["contig"], reg["start_0based"], reg["end_exclusive"]) == (CHROM, START, END)
    assert reg["one_based_closed"] == "NZ_CP058217.1:1965001-2165000"
    genome = default_build.top / "genome" / f"{NAME}.fna.gz"
    seqs = gzip.decompress(genome.read_bytes()).decode().split(">")[1:]
    assert len(seqs) == 1, "the region is one contig"
    head, *body = seqs[0].splitlines()
    assert len("".join(body)) == END - START
    assert "circular" not in head.lower(), "linear, no topology tag"


def test_bundle_version_is_the_installed_minor(default_build: Build) -> None:
    import gippa

    meta = json.loads((default_build.top / "bundle.json").read_text())
    assert meta["bundle_version"] == _minor()
    assert meta["gippa_version"] == str(gippa.__version__)


def test_bundle_out_moves_the_archive(pilot: Path, tmp_path: Path) -> None:
    out = tmp_path / "rel"
    code, log = _run(pilot, ["--out", str(tmp_path / "work"), "--bundle-out", str(out)])
    assert code == 0, log
    assert sorted(p.name for p in out.iterdir()) == [f"gippa-example-{_minor()}.tar.gz"]


# ---------------------------------------------------------------- 2. release layout


def test_archive_has_one_top_dir(default_build: Build) -> None:
    assert sorted(p.name for p in default_build.tree.iterdir()) == [f"gippa-example-{_minor()}"]
    with tarfile.open(default_build.archive) as tar:
        names = tar.getnames()
    top = f"gippa-example-{_minor()}"
    assert all(n == top or n.startswith(top + "/") for n in names), names


def test_bundle_holds_every_route_path(default_build: Build) -> None:
    files = _files(default_build.top)
    missing = sorted(_route_paths() - set(files))
    assert not missing, f"route inputs missing from the bundle: {missing}"
    assert all(files[p].stat().st_size > 0 for p in _route_paths())


def test_bundle_holds_the_top_files(default_build: Build) -> None:
    assert sorted(TOP_FILES - set(_files(default_build.top))) == []


def test_bundle_holds_only_the_layout(default_build: Build) -> None:
    """Section 2 > Layout: no sizing tars (``bundles/``), no work files."""
    allowed = TOP_FILES | _route_paths()
    extra = [p for p in _files(default_build.top) if p not in allowed
             and not p.startswith("expected/")]  # fmt: skip
    assert extra == []


def test_checksums_list_every_other_file_sorted(default_build: Build) -> None:
    top = default_build.top
    lines = (top / "CHECKSUMS.sha256").read_text().splitlines()
    parsed = [re.fullmatch(r"([0-9a-f]{64}) [ *](.+)", ln) for ln in lines if ln.strip()]
    assert all(parsed), lines
    names = [m.group(2) for m in parsed if m]
    assert names == sorted(names), "CHECKSUMS.sha256 is sorted by path"
    assert set(names) == set(_files(top)) - {"CHECKSUMS.sha256"}
    for m in parsed:
        assert m is not None
        assert hashlib.sha256((top / m.group(2)).read_bytes()).hexdigest() == m.group(1)


def test_bundle_passes_the_fetch_content_check(default_build: Build) -> None:
    from gippa.example import content_problems

    assert content_problems(default_build.top) == []


def test_selection_outputs_match_the_bundle(default_build: Build) -> None:
    top = default_build.top
    sel = json.loads((top / "SELECTION.json").read_text())
    for rel in _route_paths():
        assert rel in sel["outputs"], rel
        got = hashlib.sha256((top / rel).read_bytes()).hexdigest()
        assert sel["outputs"][rel]["sha256"] == got, rel
    for key, src in sel["sources"].items():
        assert HEX64.fullmatch(src["sha256"]), key
        assert not PurePosixPath(str(src["path"])).is_absolute(), (key, src["path"])


def test_bundle_holds_no_personal_path(default_build: Build) -> None:
    """R12 and subset rule 8: no ``$GIPPA_ROOT``, home or build path in any bundle file; the
    ``.ssl`` ``file`` column and DIA-NN's ``File.Name`` are base names."""
    roots = {str(default_build.root), PERSONAL, str(Path.home())}
    for rel, path in _files(default_build.top).items():
        text = _readable(path)
        assert not [r for r in roots if r in text], rel
    report = pq.read_table(default_build.top / "A_diann" / "report.parquet")
    assert set(report.column("File.Name").to_pylist()) == {f"{r}.mzML" for r in RUNS}
    for run in RUNS:
        with (default_build.top / "C_cascadia" / f"{run}.ssl").open() as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert rows and {r["file"] for r in rows} == {f"{run}.mzML"}


def test_licence_statement(default_build: Build) -> None:
    text = flat((default_build.top / "LICENSE-DATA.md").read_text())
    for needle in ("MIT", "PXD078459", "IPX0017169000", "GCF_021559835.1", "INSDC",
                   "grant no licence"):  # fmt: skip
        assert needle in text, needle


def test_readme_has_the_route_commands(default_build: Build) -> None:
    """``example_routes.yaml`` is the single source for the README's commands (section 1 > Pin),
    written with ``{dir}`` as ``gippa-example`` as on the docs page (AC 15)."""
    text = flat((default_build.top / "README.md").read_text())
    for rid, route in load_routes().items():
        assert route["title"] in text, rid
        for cmd in route["commands"]:
            assert fill(cmd, "gippa-example") in text, (rid, cmd)
    for needle in ("PXD078459", "mayu", "format demo", "LICENSE-DATA.md"):
        assert needle in text, needle


def test_provenance(default_build: Build) -> None:
    text = flat((default_build.top / "PROVENANCE.md").read_text())
    for needle in ("PXD078459", "GCF_021559835.1", "NZ_CP058217.1:1965001-2165000", "1965000",
                   "2165000", "SELECTION.json"):  # fmt: skip
        assert needle in text, needle


def test_bundle_json_records_the_region(default_build: Build) -> None:
    text = (default_build.top / "bundle.json").read_text()
    assert "NZ_CP058217.1:1965001-2165000" in text


def test_fetch_accepts_the_bundle(
    default_build: Build, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``gippa example fetch --from-file`` with a pin of the built archive: exit 0, the tree."""
    monkeypatch.chdir(tmp_path)
    block_https(monkeypatch)
    v = _minor()
    sha = hashlib.sha256(default_build.data).hexdigest()
    install_pins(monkeypatch, [pin_row(version=v, sha256="sha256:" + sha,
                                       size=len(default_build.data))])  # fmt: skip
    res = gippa(["example", "fetch", "--from-file", str(default_build.archive), "--outdir", "ex",
                 "--version", v])  # fmt: skip
    assert res.code == 0, res.err
    got = {k: p.read_bytes() for k, p in _files(tmp_path / "ex").items()}
    got.pop(".gippa-example.json", None)
    assert got == {k: p.read_bytes() for k, p in _files(default_build.top).items()}


# ---------------------------------------------------------------- 3. determinism


def test_two_builds_are_byte_identical(default_build: Build, tmp_path: Path) -> None:
    out = tmp_path / "second"
    code, log = _run(default_build.root, ["--out", str(tmp_path / "work2"), "--bundle-out",
                                          str(out)])  # fmt: skip
    assert code == 0, log
    again = (out / default_build.archive.name).read_bytes()
    assert hashlib.sha256(again).hexdigest() == hashlib.sha256(default_build.data).hexdigest()


def test_gzip_header_has_no_name_and_no_time(default_build: Build) -> None:
    head = default_build.data[:10]
    assert head[:3] == b"\x1f\x8b\x08"
    assert head[3] & 0x08 == 0, "gzip FNAME set"
    assert head[4:8] == b"\0\0\0\0", "gzip mtime set"


def _members(archive: Path) -> Iterator[tarfile.TarInfo]:
    with tarfile.open(archive) as tar:
        yield from tar.getmembers()


def test_tar_members_are_normalised(default_build: Build) -> None:
    members = list(_members(default_build.archive))
    names = [m.name for m in members]
    assert names == sorted(names), "members sorted"
    for m in members:
        assert m.isfile() or m.isdir(), m.name
        assert (m.mtime, m.uid, m.gid, m.uname, m.gname) == (0, 0, 0, "", ""), m.name
        assert m.mode == (0o755 if m.isdir() else 0o644), (m.name, oct(m.mode))


# ---------------------------------------------------------------- 4. the printout


def _yaml_row(out: str, key: str) -> dict[str, Any]:
    lines = out.splitlines()
    i = next((n for n, ln in enumerate(lines) if ln.lstrip().startswith(f"{key}:")), None)
    assert i is not None, f"no {key}: row in the printout"
    indent = len(lines[i]) - len(lines[i].lstrip())
    block = [lines[i]]
    for ln in lines[i + 1 :]:
        if ln.strip() and len(ln) - len(ln.lstrip()) <= indent:
            break
        block.append(ln)
    row = yaml.safe_load(textwrap.dedent("\n".join(block)))
    assert isinstance(row, dict) and isinstance(row.get(key), dict), block
    return row[key]  # type: ignore[no-any-return]


def test_printout_has_sha256_and_size(default_build: Build) -> None:
    assert hashlib.sha256(default_build.data).hexdigest() in default_build.stdout
    assert re.search(rf"\b{len(default_build.data)}\b", default_build.stdout)


def test_printout_has_a_db_sources_row(default_build: Build) -> None:
    v = _minor()
    x, y = v.split(".")
    row = _yaml_row(default_build.stdout, f"example_{x}_{y}")
    sha = hashlib.sha256(default_build.data).hexdigest()
    assert str(row["version"]) == v
    assert row["url"] == url_of(v)
    assert row["checksum"] == f"sha256:{sha}"
    assert row["size_bytes"] == len(default_build.data)
    assert row["path"] == f"example/gippa-example-{v}.tar.gz"
    assert (row["tool"], row["check"], row["size_gib"]) == (None, "sha256", 1)
    assert row["redirect_hosts"] == [ASSET_HOST]


def test_printed_db_sources_row_loads(
    default_build: Build, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Appended to ``db_sources.yaml``, the printed row is accepted by the setup loader."""
    from gippa.setupdb import items

    v = _minor()
    key = "example_" + v.replace(".", "_")
    row = _yaml_row(default_build.stdout, key)
    base = yaml.safe_load(SOURCES_YAML.read_text(encoding="utf-8"))
    base = {k: r for k, r in base.items() if not k.startswith("example_")}
    d = tmp_path / "data"
    d.mkdir()
    (d / "db_sources.yaml").write_text(yaml.safe_dump({**base, key: row}), encoding="utf-8")
    monkeypatch.setattr(items, "files", lambda _pkg: d)
    loaded = items._load_sources()[key]
    assert loaded.checksum == f"sha256:{hashlib.sha256(default_build.data).hexdigest()}"
    assert loaded.size_bytes == len(default_build.data)


def test_printout_has_a_bundles_tsv_row(default_build: Build) -> None:
    v = _minor()
    head = BUNDLES_TSV.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert head == COLUMNS
    rows = [ln.split("\t") for ln in default_build.stdout.splitlines()]
    rows = [r for r in rows if len(r) == len(COLUMNS) and r[0] == v]
    assert len(rows) == 1, "one bundles.tsv row in the printout"
    row = dict(zip(COLUMNS, rows[0], strict=True))
    sha = hashlib.sha256(default_build.data).hexdigest()
    assert row["file"] == f"gippa-example-{v}.tar.gz"
    assert int(row["size_bytes"]) == len(default_build.data)
    assert row["sha256"].removeprefix("sha256:") == sha
    assert row["url"] == url_of(v) and row["release_tag"] == f"v{v}.0"
    assert row["asset_digest"] == "", "empty until Khoa uploads the asset"
