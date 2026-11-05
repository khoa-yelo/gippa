"""M12-01 PR 3/PR 7, review round 1 of #546: ``scripts/example/build_example.py`` hardening.

docs/specs/M12-01-example-dataset.md section 2 (Subset rules 1-8, PROVENANCE.md, Licence
statement) and section 3 (reproducible, ``--check``, read-only archive), R12. Pinned here, per
review finding:

- H1: a dirty work tree fails (exit 1, "dirty") unless ``--allow-dirty``; then the commit is
  recorded as ``<sha>-dirty`` in ``bundle.json`` and ``PROVENANCE.md``.
- M1: the selection rules with exact counts in ``SELECTION.json``; rule 7 and the release
  behaviour (refusal keeps the old archive and a matching work tree, ``--force``, mode 0444, no
  ``*.partial`` / ``*.new`` left).
- M2: ``--end`` beyond the chromosome and a missing ``NZ_CP058217.1`` exit 1 with a message.
- M3: PROVENANCE.md labels (DIA-NN image, ``diann_A2.cfg``, ``a2_abundance.tsv.gz``, Cascadia
  weights, Bakta light 6.0).
- M4: ``--check --bundles-tsv TSV`` rebuilds in a temporary dir and compares with the TSV row.
- L1: a personal path inside a cell, in parquet metadata or in an ``--expected`` file fails.
- L2: a pilot symlinked from outside ``$GIPPA_ROOT`` fails (no ``..`` in SELECTION.json).
- L4: route tables keep full float precision. L5: the licence text. L6: the README region length.
- L7: a protein in two groups is one abundance row and one FASTA record; a protein missing from
  the strain database exits 1.

The pilot is a tiny synthetic one (a 50 kb chromosome, region ``[10000, 40000)``) under
``tmp_path``; the real pilot is never read. The build makes no statistical claim, so there is no
calibration test.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import os
import random
import re
import shutil
import subprocess
import tarfile
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from helpers.example_m12_01 import REPO, flat

SCRIPT_DIR = REPO / "scripts" / "example"
SCRIPT = SCRIPT_DIR / "build_example.py"
CHROM = "NZ_CP058217.1"
S, E, L = 10_000, 40_000, 50_000  # region [S, E) of a chromosome of L bp
REGION = ["--start", str(S), "--end", str(E)]
NAME = "nissle_r200k"
RUNS = ("WT_1", "WT_2", "WT_3")
PREC = 0.0006486894965844569  # L4: must round-trip exactly
LEAK = "/home/x/y"
HEX64 = re.compile(r"[0-9a-f]{64}")
COLUMNS = ["version", "file", "size_bytes", "sha256", "url", "release_tag", "build_commit",
           "build_job", "expected_commit", "asset_digest"]  # fmt: skip

# Cascadia peptides planted in the chromosome (route C)
PLUS = "MKWDHRTYPQNG"  # plus strand, inside the region: kept
IL_GENOME, IL_SSL = "GLLDKAVTERSY", "GILDKAVTERSY"  # I/L swapped in the .ssl: kept
MINUS = "WHEKCNFMQRDA"  # minus strand, inside the region: kept
SPAN = "HMWQHMWQHMWQ"  # spans the region end: dropped
RANDOM = "WWWWHHHHWWWWHHHH"  # nowhere: dropped


def _minor() -> str:
    import gippa

    m = re.match(r"(\d+)\.(\d+)", str(gippa.__version__))
    assert m is not None, gippa.__version__
    return f"{m.group(1)}.{m.group(2)}"


# ---------------------------------------------------------------- the synthetic pilot

_AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
_CODONS = [a + b + c for a in "TCAG" for b in "TCAG" for c in "TCAG"]
_ENCODE: dict[str, str] = {}
for _codon, _aa in zip(_CODONS, _AA, strict=True):
    _ENCODE.setdefault(_aa, _codon)


def _encode(prot: str) -> str:
    return "".join(_ENCODE[a] for a in prot)


def _revcomp(nt: str) -> str:
    return nt.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fasta_gz(path: Path, records: list[tuple[str, str]]) -> None:
    text = "".join(f">{n}\n" + "\n".join(s[i : i + 60] for i in range(0, len(s), 60)) + "\n"
                   for n, s in records)  # fmt: skip
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(text.encode(), mtime=0))


# db_id_map rows (1-based closed): the region edges, three group members and an outsider
PROTEINS = [
    ("P_IN_LO", S + 1, S + 600, "+", "cds"),  # start - 1 == S: kept
    ("P_OUT_LO", S, S + 600, "+", "cds"),  # start - 1 == S - 1: dropped
    ("P_IN_HI", E - 600, E, "-", "cds"),  # end == E: kept
    ("P_OUT_HI", E - 600, E + 1, "-", "cds"),  # end == E + 1: dropped
    ("P_MID1", S + 5_000, S + 5_600, "+", "cds"),
    ("P_MID2", S + 7_000, S + 7_300, "-", "sorf"),
    ("P_SHARED", S + 9_000, S + 9_600, "+", "cds"),
    ("P_FAR", 1_000, 1_600, "+", "cds"),  # outside the region
]
REGION_PROTEINS = {"P_IN_LO", "P_IN_HI", "P_MID1", "P_MID2", "P_SHARED"}

# DIA-NN rows (one per run): Protein.Ids, peptide, decoy, PG.MaxLFQ per run
DIANN = [
    ("P_IN_LO", "AAALOEDGEK", 0, (1.0, 2.0, 10.0)),  # median 2, mean 4.33
    ("P_OUT_LO", "OUTLOPEPK", 0, (7.0, 7.0, 7.0)),
    ("P_IN_HI", "HIEDGEPEPK", 0, (PREC, PREC, PREC)),
    ("P_OUT_HI", "OUTHIPEPK", 0, (7.0, 7.0, 7.0)),
    ("P_MID1;P_FAR", "MIXEDPEPK", 0, (7.0, 7.0, 7.0)),
    ("Cont_ALBU_BOVIN", "CONTAMPEPK", 0, (7.0, 7.0, 7.0)),
    ("P_MID1;Cont_ALBU_BOVIN", "MIXCONTK", 0, (7.0, 7.0, 7.0)),
    ("P_IN_LO", "KEGDEOLAAA", 1, (0.0, 0.0, 0.0)),
    ("P_IN_HI", "KPEPEGDEIH", 1, (0.0, 0.0, 0.0)),
    ("P_MID1;P_SHARED", "GROUPONEK", 0, (100.0, 200.0, 300.0)),
    ("P_MID2;P_SHARED", "GROUPTWOK", 0, (400.0, 500.0, 600.0)),
]  # fmt: skip
KEPT_PEPTIDES = {"AAALOEDGEK", "HIEDGEPEPK", "KEGDEOLAAA", "KPEPEGDEIH", "GROUPONEK", "GROUPTWOK"}


def make_pilot(
    root: Path,
    *,
    decoys: bool = True,
    contig: str = CHROM,
    leak_cell: bool = False,
    leak_meta: bool = False,
    drop_from_db: str | None = None,
) -> Path:
    """``root/pilot/pass3`` and ``root/pilot/A3_cascadia_cpu`` with the files the script reads."""
    rng = random.Random(1917)
    p3, a3 = root / "pilot" / "pass3", root / "pilot" / "A3_cascadia_cpu"
    nt = list(rng.choices("ACGT", k=L))

    def plant(pos: int, seq: str) -> None:
        nt[pos : pos + len(seq)] = list(seq)

    plant(S + 1_000, _encode(PLUS))
    plant(S + 2_000, _encode(IL_GENOME))
    plant(S + 3_000, _revcomp(_encode(MINUS)))
    plant(E - 18, _encode(SPAN))  # 6 codons inside, 6 outside
    chrom = "".join(nt)
    plasmid = "".join(rng.choices("ACGT", k=3_000))
    _fasta_gz(p3 / "genome" / "nissle1917_wt.fna.gz", [(contig, chrom), ("NZ_CP058218.1", plasmid)])

    (p3 / "db_a2").mkdir(parents=True, exist_ok=True)
    with (p3 / "db_a2" / "db_id_map.tsv").open("w") as fh:
        fh.write("db_id\tcontig\tstart\tend\tstrand\tkind\n")
        for pid, s, e, strand, kind in PROTEINS:
            fh.write(f"{pid}\t{CHROM}\t{s}\t{e}\t{strand}\t{kind}\n")
    aa = "ACDEFGHIKLMNPQRSTVWY"
    _fasta_gz(p3 / "db_a2" / "strain_db.faa.gz",
              [(p[0], "M" + "".join(rng.choices(aa, k=120))) for p in PROTEINS
               if p[0] != drop_from_db])  # fmt: skip

    rows: list[dict[str, Any]] = []
    for r, run in enumerate(RUNS):
        for j, (ids, pep, decoy, lfq) in enumerate(DIANN):
            if decoy and not decoys:
                continue
            q = PREC if pep == "HIEDGEPEPK" else 0.001 * (j + 1)
            rows.append({
                "File.Name": f"/home/users/alice/raw/{run}.mzML", "Run": run,
                "Protein.Group": ids, "Protein.Ids": ids,
                "Genes": f"a;{LEAK}" if leak_cell else "", "Stripped.Sequence": pep,
                "Modified.Sequence": pep, "Precursor.Id": f"{pep}2", "Precursor.Charge": 2,
                "Precursor.Normalised": 1000.0 * (j + 1) + 10 * r, "Q.Value": q,
                "Global.Q.Value": q, "PEP": q, "PG.MaxLFQ": lfq[r], "Decoy": decoy,
            })  # fmt: skip
    table = pa.Table.from_pandas(pd.DataFrame(rows), preserve_index=False)
    if leak_meta:
        meta = dict(table.schema.metadata or {})
        meta[b"origin"] = f"{LEAK}/search/report.parquet".encode()
        table = table.replace_schema_metadata(meta)
    (p3 / "A2").mkdir(parents=True, exist_ok=True)
    pq.write_table(table, p3 / "A2" / "A2_report.parquet")

    for run in RUNS:
        d = a3 / run
        d.mkdir(parents=True, exist_ok=True)
        with (d / f"{run}.ssl").open("w") as fh:
            fh.write("file\tscan\tcharge\tsequence\tscore-type\tscore\n")
            for i, pep in enumerate((PLUS, IL_SSL, MINUS, SPAN, RANDOM)):
                fh.write(f"/home/users/alice/raw/{run}.mzML\t{100 + i}\t2\t{pep}\tUNKNOWN\t0.9\n")
    return root


# ---------------------------------------------------------------- running the build


def _load(script: Path = SCRIPT) -> ModuleType:
    spec = importlib.util.spec_from_file_location(f"gippa_build_example_{id(script)}", script)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run(root: Path, argv: list[str], script: Path = SCRIPT) -> tuple[int, str]:
    """``build_example.main(argv)`` with ``GIPPA_ROOT=root``: (exit code, stdout + stderr)."""
    mod = _load(script)
    out = io.StringIO()
    with (
        pytest.MonkeyPatch.context() as mp,
        contextlib.redirect_stdout(out),
        contextlib.redirect_stderr(out),
    ):
        mp.setenv("GIPPA_ROOT", str(root))
        mp.delenv("SLURM_JOB_ID", raising=False)
        try:
            code = mod.main(argv)
        except SystemExit as exc:  # argparse
            code = exc.code if isinstance(exc.code, int) else 1
    return int(code or 0), out.getvalue()


@dataclass(frozen=True)
class Dirs:
    root: Path
    work: Path
    rel: Path

    @property
    def archive(self) -> Path:
        return self.rel / f"gippa-example-{_minor()}.tar.gz"

    @property
    def tree(self) -> Path:
        """The work tree (``OUT/<name>/gippa-example-<X.Y>/``, the script's docstring)."""
        return self.work / NAME / f"gippa-example-{_minor()}"

    def argv(self, *extra: str, dirty_ok: bool = True) -> list[str]:
        flags = ["--allow-dirty"] if dirty_ok else []
        return [*REGION, "--out", str(self.work), "--bundle-out", str(self.rel), *flags, *extra]


def dirs(tmp: Path, **pilot: Any) -> Dirs:
    root = tmp / "gippa_root"
    root.mkdir(parents=True, exist_ok=True)
    make_pilot(root.resolve(), **pilot)
    return Dirs(root.resolve(), tmp / "work", tmp / "rel")


def unpack(archive: Path, dest: Path) -> Path:
    from gippa.setupdb.run import unpack_archive

    dest.mkdir(parents=True)
    assert unpack_archive(archive, dest) is None, "unsafe member in the bundle"
    return dest / f"gippa-example-{_minor()}"


def files(top: Path) -> dict[str, bytes]:
    return {p.relative_to(top).as_posix(): p.read_bytes() for p in sorted(top.rglob("*"))
            if p.is_file()}  # fmt: skip


def leftovers(*tops: Path) -> list[str]:
    return [str(p) for t in tops if t.exists() for p in t.rglob("*")
            if p.name.endswith((".partial", ".new")) or ".partial-" in p.name]  # fmt: skip


def tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def read_fasta(path: Path) -> list[str]:
    return [ln[1:].split()[0] for ln in path.read_text().splitlines() if ln.startswith(">")]


@dataclass(frozen=True)
class Built:
    d: Dirs
    out: str
    top: Path

    def sel(self) -> dict[str, Any]:
        return json.loads((self.top / "SELECTION.json").read_text())  # type: ignore[no-any-return]


@pytest.fixture(scope="module")
def base(tmp_path_factory: pytest.TempPathFactory) -> Built:
    """One build of the default synthetic pilot over ``[S, E)``."""
    tmp = tmp_path_factory.mktemp("base")
    d = dirs(tmp)
    code, out = run(d.root, d.argv())
    assert code == 0, out
    return Built(d, out, unpack(d.archive, tmp / "unpacked"))


# ---------------------------------------------------------------- H1: dirty tree


def _git(repo: Path, *args: str) -> str:
    env = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.org",
           "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.org"}  # fmt: skip
    for k in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE"):
        env.pop(k, None)
    res = subprocess.run(["git", "-C", str(repo), "-c", "commit.gpgsign=false", *args],
                         capture_output=True, text=True, check=True, env=env)  # fmt: skip
    return res.stdout.strip()


@pytest.fixture
def repo(tmp_path: Path) -> tuple[Path, str]:
    """A temporary git checkout holding a copy of ``scripts/example/*.py``; (script, HEAD sha)."""
    r = tmp_path / "checkout"
    (r / "scripts" / "example").mkdir(parents=True)
    for p in SCRIPT_DIR.glob("*.py"):
        shutil.copyfile(p, r / "scripts" / "example" / p.name)
    (r / "notes.txt").write_text("clean\n")
    _git(r, "init", "-q")
    _git(r, "add", "-A")
    _git(r, "commit", "-q", "-m", "init")
    return r / "scripts" / "example" / "build_example.py", _git(r, "rev-parse", "HEAD")


def _commit_of(top: Path) -> tuple[str, str]:
    meta = json.loads((top / "bundle.json").read_text())
    return str(meta["build_commit"]), (top / "PROVENANCE.md").read_text()


def test_clean_tree_records_the_plain_sha(tmp_path: Path, repo: tuple[Path, str]) -> None:
    script, sha = repo
    d = dirs(tmp_path)
    code, out = run(d.root, d.argv(dirty_ok=False), script)
    assert code == 0, out
    commit, prov = _commit_of(unpack(d.archive, tmp_path / "u"))
    assert commit == sha
    assert sha in prov and "-dirty" not in prov


def test_dirty_tree_fails_without_allow_dirty(tmp_path: Path, repo: tuple[Path, str]) -> None:
    script, _ = repo
    (script.parents[2] / "notes.txt").write_text("changed\n")
    d = dirs(tmp_path)
    code, out = run(d.root, d.argv(dirty_ok=False), script)
    assert code == 1, out
    assert "dirty" in out.lower(), out
    assert not d.archive.exists()
    assert leftovers(d.work, d.rel) == []


def test_allow_dirty_records_sha_dirty(tmp_path: Path, repo: tuple[Path, str]) -> None:
    script, sha = repo
    (script.parents[2] / "notes.txt").write_text("changed\n")
    d = dirs(tmp_path)
    code, out = run(d.root, d.argv(dirty_ok=True), script)
    assert code == 0, out
    commit, prov = _commit_of(unpack(d.archive, tmp_path / "u"))
    assert commit == f"{sha}-dirty"
    assert f"{sha}-dirty" in prov


# ---------------------------------------------------------------- M1: selection rules


def test_region_edges_select_the_proteins(base: Built) -> None:
    """Rule 2: 1-based start ``S + 1`` and end ``E`` are inside; ``S`` and ``E + 1`` are not."""
    counts = base.sel()["counts"]
    assert counts["region_proteins"] == len(REGION_PROTEINS) == 5
    assert (counts["region_proteins_cds"], counts["region_proteins_sorf"]) == (4, 1)
    ids = pq.read_table(base.top / "A_diann" / "report.parquet").column("Protein.Ids").to_pylist()
    got = {p for x in ids for p in x.split(";")}
    assert {"P_IN_LO", "P_IN_HI"} <= got
    assert not got & {"P_OUT_LO", "P_OUT_HI"}


def test_route_a_counts_are_exact(base: Built) -> None:
    """Rule 3: mixed-ID and ``Cont_`` rows dropped, decoys kept (6 rows per run of 11)."""
    a = base.sel()["counts"]["A"]
    assert (a["rows_total"], a["rows_kept"]) == (33, 18)
    assert (a["target_rows_kept"], a["decoy_rows_kept"]) == (12, 6)
    assert a["rows_per_run"] == {r: 6 for r in RUNS}
    t = pq.read_table(base.top / "A_diann" / "report.parquet").to_pandas()
    assert set(t["Stripped.Sequence"]) == KEPT_PEPTIDES
    assert not t["Protein.Ids"].str.contains("Cont_|P_FAR").any()
    assert set(t.loc[t["Decoy"] == 1, "Protein.Ids"]) == {"P_IN_LO", "P_IN_HI"}


def test_route_c_keeps_minus_strand_and_il_swapped(base: Built) -> None:
    """Rule 5: plus, minus-strand and I/L-swapped peptides kept; spanning the end and random
    dropped (3 of 5 rows per run)."""
    c = base.sel()["counts"]["C"]
    assert c == {r: {"rows_total": 5, "rows_kept": 3} for r in RUNS}
    for r in RUNS:
        kept = set(tsv(base.top / "C_cascadia" / f"{r}.ssl")["sequence"])
        assert kept == {PLUS, IL_SSL, MINUS}, r


def test_route_d_counts_are_exact(base: Built) -> None:
    dd = base.sel()["counts"]["D"]
    assert (dd["rows"], dd["decoy_rows"]) == (18, 6)
    assert len(tsv(base.top / "D_generic" / "peptides.tsv")) == 18


def test_abundance_is_the_median_over_runs(base: Built) -> None:
    """Rule 4: runs 1, 2, 10 give 2 (the median), not 4.33 (the mean)."""
    ab = tsv(base.top / "B_proteins" / "abundance.tsv")
    vals = ab.loc[ab["protein_id"] == "P_IN_LO", "maxlfq"].tolist()
    assert [float(v) for v in vals] == [2.0]


# ---------------------------------------------------------------- M1: release behaviour


def test_no_decoy_fails_with_empty_route(tmp_path: Path) -> None:
    d = dirs(tmp_path, decoys=False)
    code, out = run(d.root, d.argv())
    assert code == 1, out
    assert "empty route" in out
    assert not d.archive.exists()
    assert leftovers(d.work, d.rel, d.root) == []


def _expected(tmp: Path, n: int) -> Path:
    e = tmp / f"expected{n}"
    e.mkdir()
    (e / "A.json").write_text(json.dumps({"n_genes": {"value": n, "tolerance": 0}}) + "\n")
    return e


Refused = tuple[Dirs, bytes, int, str, Path]


@pytest.fixture
def refused(tmp_path: Path) -> Refused:
    """A build, then a second build with another ``--expected``: (dirs, first archive bytes,
    the second build's exit code and log, the second expected dir)."""
    d = dirs(tmp_path)
    code, out = run(d.root, d.argv("--expected", str(_expected(tmp_path, 1))))
    assert code == 0, out
    first = d.archive.read_bytes()
    e2 = _expected(tmp_path, 2)
    code, out = run(d.root, d.argv("--expected", str(e2)))
    return d, first, code, out, e2


def test_refusal_keeps_the_old_archive(refused: Refused) -> None:
    d, first, code, out, _ = refused
    assert code == 1, f"a different archive was not refused:\n{out}"
    assert d.archive.read_bytes() == first
    assert "--force" in out
    assert leftovers(d.work, d.rel) == []


def test_refusal_keeps_the_work_tree_matching_the_archive(refused: Refused, tmp_path: Path) -> None:
    """The build writes into a temporary dir and swaps it in only after the release succeeds."""
    d, _, code, out, _ = refused
    assert code == 1, out
    kept = unpack(d.archive, tmp_path / "kept")
    assert d.tree.is_dir(), d.tree
    assert files(d.tree) == files(kept)


def test_force_replaces_the_archive(refused: Refused, tmp_path: Path) -> None:
    d, first, _, _, e2 = refused
    code, out = run(d.root, d.argv("--expected", str(e2), "--force"))
    assert code == 0, out
    assert d.archive.read_bytes() != first
    top = unpack(d.archive, tmp_path / "forced")
    assert json.loads((top / "expected" / "A.json").read_text())["n_genes"]["value"] == 2
    assert d.archive.stat().st_mode & 0o777 == 0o444


def test_archive_is_read_only(base: Built) -> None:
    assert base.d.archive.stat().st_mode & 0o777 == 0o444


def test_no_partial_file_after_success(base: Built) -> None:
    assert leftovers(base.d.work, base.d.rel) == []


# ---------------------------------------------------------------- M2: region checks


def test_end_beyond_the_chromosome_fails(tmp_path: Path) -> None:
    d = dirs(tmp_path)
    code, out = run(d.root, [*d.argv(), "--end", str(L + 10_000)])
    assert code == 1, out
    assert CHROM in out and re.search(rf"\b{L}\b|\b{L:,}\b", out), out
    assert not d.archive.exists()


def test_missing_chromosome_fails_with_a_message(tmp_path: Path) -> None:
    d = dirs(tmp_path, contig="NZ_OTHER.1")
    code, out = run(d.root, d.argv())  # a KeyError escaping main() fails the test
    assert code == 1, out
    assert CHROM in out and "KeyError" not in out and "Traceback" not in out
    assert not d.archive.exists()


# ---------------------------------------------------------------- M4: --check


def _bundles_tsv(path: Path, sha: str, size: int) -> Path:
    v = _minor()
    row = [v, f"gippa-example-{v}.tar.gz", str(size), f"sha256:{sha}",
           f"https://github.com/khoa-yelo/gippa/releases/download/v{v}.0/gippa-example-{v}.tar.gz",
           f"v{v}.0", "", "", "", ""]  # fmt: skip
    path.write_text("\t".join(COLUMNS) + "\n" + "\t".join(row) + "\n")
    return path


def _check_argv(d: Dirs, tsv_path: Path, rel: Path) -> list[str]:
    return [*REGION, "--allow-dirty", "--check", "--bundles-tsv", str(tsv_path),
            "--out", str(d.work), "--bundle-out", str(rel)]  # fmt: skip


def test_check_matches_the_recorded_row(base: Built, tmp_path: Path) -> None:
    data = base.d.archive.read_bytes()
    t = _bundles_tsv(tmp_path / "bundles.tsv", hashlib.sha256(data).hexdigest(), len(data))
    rel = tmp_path / "rel"
    rel.mkdir()
    (rel / "sentinel").write_text("x")
    code, out = run(base.d.root, _check_argv(Dirs(base.d.root, tmp_path / "w", rel), t, rel))
    assert code == 0, out
    assert sorted(p.name for p in rel.iterdir()) == ["sentinel"], "--check wrote to --bundle-out"


def test_check_mismatch_names_both_shas(base: Built, tmp_path: Path) -> None:
    data = base.d.archive.read_bytes()
    real, fake = hashlib.sha256(data).hexdigest(), "0" * 64
    t = _bundles_tsv(tmp_path / "bundles.tsv", fake, len(data))
    rel = tmp_path / "rel"
    rel.mkdir()
    shutil.copyfile(base.d.archive, rel / base.d.archive.name)
    code, out = run(base.d.root, _check_argv(Dirs(base.d.root, tmp_path / "w", rel), t, rel))
    assert code == 1, out
    assert real in out and fake in out, out
    assert sorted(p.name for p in rel.iterdir()) == [base.d.archive.name]
    assert (rel / base.d.archive.name).read_bytes() == data


# ---------------------------------------------------------------- L1: leak guard


@pytest.mark.parametrize("where", ["cell", "metadata"])
def test_personal_path_in_the_report_fails(tmp_path: Path, where: str) -> None:
    d = dirs(tmp_path, leak_cell=where == "cell", leak_meta=where == "metadata")
    code, out = run(d.root, d.argv())
    assert code == 1, f"{LEAK} ({where}) reached the bundle:\n{out}"
    assert "report.parquet" in out, out
    assert not d.archive.exists()


def test_personal_path_in_expected_fails(tmp_path: Path) -> None:
    d = dirs(tmp_path)
    e = tmp_path / "expected"
    e.mkdir()
    (e / "A.json").write_text(json.dumps({"outdir": f"{LEAK}/runs/A"}) + "\n")
    code, out = run(d.root, d.argv("--expected", str(e)))
    assert code == 1, out
    assert "A.json" in out, out
    assert not d.archive.exists()


# ---------------------------------------------------------------- L2: symlinked pilot


def test_pilot_symlinked_from_outside_the_root_fails(tmp_path: Path) -> None:
    outside = make_pilot((tmp_path / "outside").resolve())
    root = tmp_path / "gippa_root"
    root.mkdir()
    (root / "pilot").symlink_to(outside / "pilot", target_is_directory=True)
    d = Dirs(root.resolve(), tmp_path / "work", tmp_path / "rel")
    code, out = run(d.root, d.argv())
    assert code == 1, out
    assert not d.archive.exists()


def test_selection_source_paths_are_lexical_relative(base: Built) -> None:
    for key, src in base.sel()["sources"].items():
        parts = Path(src["path"]).parts
        assert ".." not in parts and not Path(src["path"]).is_absolute(), (key, src["path"])


# ---------------------------------------------------------------- L4: precision


def test_route_tables_keep_full_precision(base: Built) -> None:
    ab = tsv(base.top / "B_proteins" / "abundance.tsv")
    vals = ab.loc[ab["protein_id"] == "P_IN_HI", "maxlfq"].tolist()
    assert len(vals) == 1 and float(vals[0]) == PREC, vals
    d = tsv(base.top / "D_generic" / "peptides.tsv")
    rows = d[d["peptide"] == "HIEDGEPEPK"]
    assert len(rows) == 3
    assert all(float(v) == PREC for v in rows["q_value"]), rows["q_value"].tolist()
    assert all(float(v) == PREC for v in rows["pep"]), rows["pep"].tolist()


# ---------------------------------------------------------------- L6: README


def test_readme_region_length_follows_the_region(base: Built) -> None:
    text = flat((base.top / "README.md").read_text())
    assert "200 kb" not in text
    assert re.search(r"\b30 kb\b|\b30,?000 bp\b", text), "the README does not say 30 kb"


def test_readme_does_not_call_536_unfixed(base: Built) -> None:
    text = flat((base.top / "README.md").read_text())
    assert "#536" not in text and "issue 536" not in text


# ---------------------------------------------------------------- L7: protein groups


def test_protein_in_two_groups_is_one_row_and_one_record(base: Built) -> None:
    ab = tsv(base.top / "B_proteins" / "abundance.tsv")
    assert ab["protein_id"].tolist().count("P_SHARED") == 1, ab["protein_id"].tolist()
    assert sorted(ab["protein_id"]) == sorted(REGION_PROTEINS)
    names = read_fasta(base.top / "B_proteins" / "proteins.faa")
    assert sorted(names) == sorted(REGION_PROTEINS), names
    b = base.sel()["counts"]["B"]
    assert (b["abundance_rows"], b["fasta_records"]) == (5, 5)


def test_protein_missing_from_the_database_fails(tmp_path: Path) -> None:
    d = dirs(tmp_path, drop_from_db="P_MID2")
    code, out = run(d.root, d.argv())  # a KeyError escaping main() fails the test
    assert code == 1, out
    assert "P_MID2" in out and "KeyError" not in out
    assert not d.archive.exists()


# ---------------------------------------------------------------- L5: licence text


def _sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+", text) if s]


def test_licence_names_the_iprox_record(base: Built) -> None:
    text = flat((base.top / "LICENSE-DATA.md").read_text()).replace("**", "")
    assert "PXD078459 (iProX IPX0017169000)" in text


def test_licence_says_b_and_d_are_derived_tables(base: Built) -> None:
    text = flat((base.top / "LICENSE-DATA.md").read_text()).replace("**", "")
    for needle in ("median", "group", "proforma"):
        assert needle in text.lower(), needle
    both = re.compile(r"B_proteins.*D_generic|routes B and D")
    bd = [s for s in _sentences(text) if both.search(s)]
    assert any("derived" in s.lower() for s in bd), "no sentence says B and D are derived"


def test_licence_does_not_call_b_and_d_unchanged(base: Built) -> None:
    """Any "unchanged" claim names the route A or C files only, never B or D."""
    text = flat((base.top / "LICENSE-DATA.md").read_text()).replace("**", "")
    for s in _sentences(text):
        if "unchanged" not in s.lower():
            continue
        assert re.search(r"A_diann|C_cascadia|routes? A|routes? C", s), s
        assert not re.search(r"B_proteins|D_generic|routes? B|routes? D|and D\b", s), s


# ---------------------------------------------------------------- M3: PROVENANCE labels


def _near_hex(text: str, label: str, span: int = 300) -> bool:
    i = text.find(label)
    return i >= 0 and HEX64.search(text, i, i + len(label) + span) is not None


def test_provenance_lists_tool_and_table_digests(base: Built) -> None:
    text = flat((base.top / "PROVENANCE.md").read_text())
    assert "validation/pilot/diann_A2.cfg" in text
    m = re.search(r"DIA-NN[^.]{0,200}?image", text)
    assert m, "no DIA-NN image label"
    assert HEX64.search(text, m.end(), m.end() + 300), "no DIA-NN image sha256"
    assert _near_hex(text, "a2_abundance.tsv.gz"), "no a2_abundance.tsv.gz sha256"
    m = re.search(r"Cascadia[^.]{0,200}?weights", text)
    assert m, "no Cascadia weights label"
    assert HEX64.search(text, m.end(), m.end() + 300), "no Cascadia weights sha256"
    assert re.search(r"Bakta[^.]{0,100}?light[^.]{0,40}?6\.0", text), "no Bakta light 6.0 DB"


# ---------------------------------------------------------------- the archive itself


def test_archive_members_are_regular(base: Built) -> None:
    with tarfile.open(base.d.archive) as tar:
        assert all(m.isfile() or m.isdir() for m in tar.getmembers())
