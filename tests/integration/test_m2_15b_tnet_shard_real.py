"""M2-15b: the sharded TerminatorNet scan with the real model is byte-identical to upstream
(docs/specs/M2-15b-terminatornet-sharding.md "Tests", ``needs_conda``; examples 3, 5-10), except
on records with non-ACGT bases, where upstream shifts its calls and the sharded scan reports them
in original coordinates (issue #615).

``workflow/bin/gippa_terminatornet_shard`` runs ``plan --core N``, every ``chunk <id> --cpus N``
and ``merge`` in the ``gippa-terminatornet`` prebuilt env (``$GIPPA_ROOT/micromamba/envs``), with
``model.pkl`` symlinked into the cwd as the unsharded module does. Small synthetic FASTAs cut from
the sim100k genome are compared with ``TerminatorNet.py -s -n 2`` (the reference, lead decision 2);
then sim100k at CORE = 10,000 must equal the committed recording
``tests/data/recorded/terminators/terminatornet/sim100k/terminatornet.tsv`` (``-n 17``)
byte for byte.
The genome is ``$GIPPA_ROOT/recordings/terminatornet/sim100k/genome.fa`` (skipped when missing).

``slow`` (minutes of real scanning): run on Slurm (scripts/slurm/test_branch.sbatch / nightly).
"""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "workflow" / "bin" / "gippa_terminatornet_shard"
ROOT = Path(os.environ.get("GIPPA_ROOT", REPO / "local"))
ENV = ROOT / "micromamba" / "envs" / "gippa-terminatornet"
HOME = ENV / "share" / "terminatornet"
GENOME = ROOT / "recordings" / "terminatornet" / "sim100k" / "genome.fa"
RECORDED = (
    REPO
    / "tests"
    / "data"
    / "recorded"
    / "terminators"
    / "terminatornet"
    / "sim100k"
    / "terminatornet.tsv"
)

pytestmark = [
    pytest.mark.slow,
    pytest.mark.needs_conda,
    pytest.mark.skipif(not ENV.is_dir(), reason=f"missing capability: tool env {ENV}"),
    pytest.mark.skipif(not GENOME.is_file(), reason=f"no sim100k recording genome at {GENOME}"),
]


def _env() -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "TERMINATORNET_HOME")}
    env.update(
        TERMINATORNET_HOME=str(HOME), PYTHONNOUSERSITE="1", PATH=f"{ENV / 'bin'}:{env['PATH']}"
    )
    return env


def _py(cwd: Path, *args: str) -> None:
    r = subprocess.run(
        [str(ENV / "bin" / "python"), *args],
        cwd=cwd,
        env=_env(),
        capture_output=True,
        text=True,
        check=False,
    )
    assert r.returncode == 0, f"{args}: {r.stderr[-2000:]}"


def _dir(tmp: Path, name: str, fasta: str) -> Path:
    wd = tmp / name
    wd.mkdir()
    (wd / "genome.fa").write_text(fasta)
    (wd / "model.pkl").symlink_to(HOME / "model.pkl")
    return wd


def _sharded(wd: Path, core: int, cpus: int = 2) -> bytes:
    _py(wd, str(SCRIPT), "plan", "--core", str(core))
    with (wd / "plan.tsv").open(newline="") as fh:
        ids = [r["chunk_id"] for r in csv.DictReader(fh, delimiter="\t")]
    for cid in ids:
        _py(wd, str(SCRIPT), "chunk", cid, "--cpus", str(cpus))
    _py(wd, str(SCRIPT), "merge")
    return (wd / "terminatornet.tsv").read_bytes()


def _unsharded(wd: Path) -> bytes:
    _py(wd, str(HOME / "TerminatorNet.py"), "-f", "genome.fa", "-s", "-n", "2", "-o", "ref.tsv")
    return (wd / "ref.tsv").read_bytes()


def _chr1() -> str:
    seqs: dict[str, list[str]] = {}
    name = ""
    for line in GENOME.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        else:
            seqs[name].append(line.strip())
    return "".join(seqs["chr1"])


def _rc(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _fields(line: str) -> list[str]:
    return line.split("\t")


def _same_except_gapped(got: bytes, ref: bytes, gapped: dict[str, str]) -> None:
    """Issue #615: records without non-ACGT bases equal the unsharded output line for line;
    on a gapped record (unsharded TerminatorNet shifts every call after the first removed base),
    every call is in original coordinates: no removed base inside its span, hairpin (and the
    8-nt tracts when they hold no removed base) equal to the genome there, and the calls ending
    200 nt or more before the first removed base equal the unsharded ones (TERM numbers aside)."""
    g_lines, r_lines = got.decode().splitlines(), ref.decode().splitlines()
    assert [x for x in g_lines if x.split("\t")[0] not in gapped] == [
        x for x in r_lines if x.split("\t")[0] not in gapped
    ]
    checked = 0
    for name, raw in gapped.items():
        up = raw.upper()
        first_gap = next(i for i, c in enumerate(up) if c not in "ACGT")
        calls = [x.split("\t") for x in g_lines if x.split("\t")[0] == name and x.count("\t") == 10]
        checked += len(calls)
        for f in calls:
            s, e = int(f[8]), int(f[9])
            seg = up[s - 1 : e]
            assert set(seg) <= set("ACGT"), f"{f}: span holds a removed base"
            assert (seg if f[10] == "+" else _rc(seg)) == f[3], f"{f}: hairpin not at {s}..{e}"
            wide = up[max(s - 9, 0) : e + 8]
            if s > 8 and set(wide) <= set("ACGT"):
                assert (wide if f[10] == "+" else _rc(wide)) == f[2] + f[3] + f[4], f
        ref_calls = [f for f in map(_fields, r_lines) if f[0] == name and len(f) == 11]
        early_ref = [f[2:] for f in ref_calls if int(f[9]) <= first_gap - 200]
        early_got = [f[2:] for f in calls if int(f[9]) <= first_gap - 200]
        assert early_got == early_ref, f"{name}: calls before the first gap changed"
    if "y" in gapped:  # record y keeps gapped calls; record c's only call spans its N run
        assert checked, "no gapped call was checked"


def test_short_mode_examples_3_5_8_10(tmp_path: Path) -> None:
    s = _chr1()
    dirty = s[1800:1900] + "NNNN" + s[1900:2100].lower()  # 300 cleaned nt, Ns and lower case
    fasta = f">w40\n{s[1850:1890]}\n>e\n>t\nAC\n>c\n{dirty}\n>u\n{s[2800:3250]}\n"
    a, b = _dir(tmp_path, "sharded", fasta), _dir(tmp_path, "ref", fasta)
    got = _sharded(a, core=100)
    _same_except_gapped(got, _unsharded(b), {"c": dirty})
    assert b"t\tNO TERMINATORS FOUND\n" in got and b"e\t" not in got


def test_long_mode_chunk_edges_examples_6_7_8(tmp_path: Path) -> None:
    s = _chr1()
    y = f"{s[30000:31000]}NNacgt{s[31000:33000].lower()}"
    fasta = f">x\n{s[:6000]}\n>y\n{y}\n"
    a, b = _dir(tmp_path, "sharded", fasta), _dir(tmp_path, "ref", fasta)
    _same_except_gapped(_sharded(a, core=1_000), _unsharded(b), {"y": y})


def test_example_9_chunk_bytes_at_1_4_8_cpus(tmp_path: Path) -> None:
    wd = _dir(tmp_path, "wd", f">x\n{_chr1()[:4000]}\n")
    _py(wd, str(SCRIPT), "plan", "--core", "1000")
    blobs = set()
    for cpus in (1, 4, 8):
        _py(wd, str(SCRIPT), "chunk", "c00000p0001", "--cpus", str(cpus))
        blobs.add((wd / "calls_c00000p0001.pkl").read_bytes())
    assert len(blobs) == 1


def test_sim100k_core_10000_equals_the_committed_recording(tmp_path: Path) -> None:
    wd = _dir(tmp_path, "sim100k", GENOME.read_text())
    assert _sharded(wd, core=10_000, cpus=4) == RECORDED.read_bytes()
