"""TerminatorNet calls in original-genome coordinates on gapped assemblies (issue #615, 1.1.1).

Upstream ``clean_sequences`` deletes every non-ACGT character (N runs, IUPAC codes) before the scan
(it keeps U as T and upper-cases; a record with < 75 % ACGT is left uncleaned and < 4 nt is
blanked). ``workflow/bin/gippa_terminatornet_shard`` scans the cleaned sequence, so before the fix
every call after the first removed base is shifted, silently. Lead rule (2026-11-22):

1. ``terminatornet.tsv`` coordinates are in the original record (as upstream ``read_in_sequences``
   returns it) on both strands, whatever the number and place of removed bases (position 0, last
   position, several runs, lower case, IUPAC);
2. a call whose span (hairpin start..stop, the call tuple's first two fields) would contain a
   removed base in the original record is dropped, before upstream ``merge_overlapping``;
3. ``merge`` writes exactly one stderr line with ``WARNING`` and ``non-ACGT`` when the genome lost
   at least one base, carrying the number of removed bases (``<n> non-ACGT``) and of dropped calls
   (``<m> call`` / ``<m> calls``); none on a gap-free genome; ``plan`` and ``chunk`` write none;
4. plan rule (kept, spec M2-15b example 8): ``plan.tsv`` stays in *scan* coordinates, i.e.
   ``seq_length`` is the cleaned length and the cores tile ``nwin(cleaned length)``; the offset map
   is not stored but derived in ``merge`` from the re-read original record by a pure helper;
5. a gap-free genome gives byte-identical output (sha256 recorded on cff0fd68).

Seam pinned here (pure, no upstream import, loaded by path like the M2-15b tests)::

    to_original(calls, raw, cleaned, strand) -> tuple[list[Call], int]

``calls`` are chunk call tuples ``(start, stop, ...)``, 0-based inclusive in the scanned strand
sequence (``cleaned`` for ``+``, its reverse complement for ``-``); ``raw`` is the record before
``clean_sequences`` (mixed case allowed), ``cleaned`` what upstream made of it. It returns the
kept calls in input order with ``start``/``stop`` in the original strand sequence (``raw`` for
``+``, the reverse complement of ``raw`` for ``-``, so ``merge_overlapping(kept, len(raw), '-')``
yields forward original coordinates) and every other element the identical object, and the
number of dropped calls. When ``len(cleaned) == len(raw)`` (gap-free, or a record upstream did not
clean) it returns the calls unchanged and 0.

The CLI tests use a fake upstream of their own (written to tmp_path): the window scorer calls a
12-nt hairpin wherever a window starts with the palindrome ``GGATCC`` (so each site is seen on both
strands), with ``hairpin = SEQUENCE[j:j+12]``, which lets every output line be checked against the
genome. Its ``clean_sequences`` and ``merge_overlapping`` / ``output_terms`` follow upstream
(in-place cleaning, the 75 % rule, the ``-`` strand mapped with ``seq_length`` and sorted by start).
Fillers use only A, C, T so no site appears by chance.
"""

from __future__ import annotations

import hashlib
import importlib.machinery
import importlib.util
import itertools
import os
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

SCRIPT = Path(__file__).resolve().parents[2] / "workflow" / "bin" / "gippa_terminatornet_shard"
MOTIF, SPAN = "GGATCC", 12
_COMP = str.maketrans("ACGT", "TGCA")


def _shard() -> ModuleType:
    name = "gippa_terminatornet_shard"
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_loader(
        name, importlib.machinery.SourceFileLoader(name, str(SCRIPT))
    )
    assert spec is not None and spec.loader is not None, f"cannot load {SCRIPT}"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# ---------------------------------------------------------------- fake upstream

FEATURES_PY = """
import numpy as np
SEQUENCE, GC_PERCENT = "", 0.5
HAIRPIN_LENGTH, TAIL_LENGTH, WINDOW, STEP = 52, 8, 60, 1
FEATURES = ["f_pred"]

def reverse_complement(s):
    s = s.replace('A', '?').replace('T', 'A').replace('?', 'T')
    s = s.replace('C', '?').replace('G', 'C').replace('?', 'G')
    return s[::-1]

def scan_oneprocess_1(j):
    if SEQUENCE[j:j + 6] == "GGATCC":
        return (SEQUENCE[j:j + 52], j, 40.0)
    return "", -1, -1

def scan_oneprocess_2(candidate):
    hrpn, j, tt = candidate
    feats = {"f_pred": 0.5 + (j % 7) / 100, "hairpin_start": 0, "hairpin_end": 11,
             "structure": "..((((....))))..", "GH": -12.5}
    return (feats, SEQUENCE[j:j + 12], "TTTTTTTT", "AAAAAAAA", tt, j)

class FakeScaler:
    def transform(self, x):
        return np.asarray(x, dtype=float)

class FakeModel:
    def predict_proba(self, x):
        p = float(x[0, 0])
        return np.array([[1.0 - p, p]])
"""

TERMINATORNET_PY = """
import sys
from operator import itemgetter
from Features import *
OUT = sys.stdout
BOTH_STRANDS = False

def read_in_sequences(filename):
    seqs, names, s, name = [], [], "", None
    for line in open(filename):
        if line.startswith(">"):
            if s:
                seqs.append(s); names.append(name)
            s, name = "", line[1:].strip()
        else:
            s += line.strip()
    if s:
        seqs.append(s); names.append(name)
    return seqs, names

def clean_sequences(sequences):
    long_seq = False
    for i in range(len(sequences)):
        if len(sequences[i]) < 4:
            sequences[i] = ""
            continue
        s = sequences[i].replace(" ", "").upper().replace("U", "T")
        acgt = sum(s.count(c) for c in "ACGT")
        if acgt < 0.75 * len(s):
            sys.stderr.write("\\n\\nWARNING - provided sequence contains too many "
                             "non-nucleotide characters\\n\\n")
        elif acgt != len(s):
            s = "".join(c for c in s if c in "ACGT")
        sequences[i] = s
        long_seq = long_seq or len(s) >= 501
    return sequences, 0.5, long_seq

def merge_overlapping(terms, seq_length=0, strand="+"):
    out, index = [], 0
    while index < len(terms):
        best = terms[index]
        j = index + 1
        while j < len(terms) and terms[j][0] <= best[1]:
            if terms[j][7] >= best[7]:
                best = terms[j]
            j += 1
        start, stop, a, h, u, st_, gh, pred = best
        st_ = st_[st_.find("("):st_.rfind(")") + 1]
        term = [a, h, u, st_, gh, pred, start + 1, stop + 1]
        index = j
        if BOTH_STRANDS:
            if strand == "-":
                stop_rc, start_rc = seq_length - term[-2], seq_length - term[-1]
                term[-2], term[-1] = start_rc + 1, stop_rc + 1
            term.append(strand)
        out.append(term)
    if BOTH_STRANDS and strand == "-":
        out.sort(key=itemgetter(6))
    return out

def output_terms(terms, name):
    if not terms:
        OUT.write(name + "\\tNO TERMINATORS FOUND\\n")
    for i, t in enumerate(terms):
        OUT.write(name + "\\tTERM_" + str(i + 1) + "\\t" + "\\t".join(map(str, t)) + "\\n")
"""


@pytest.fixture
def home(tmp_path: Path) -> Path:
    h = tmp_path / "tnet_home"
    h.mkdir()
    (h / "Features.py").write_text(FEATURES_PY)
    (h / "TerminatorNet.py").write_text(TERMINATORNET_PY)
    return h


def _env(home: Path, **extra: str) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if k not in ("PYTHONPATH", "TERMINATORNET_HOME")}
    env.update(TERMINATORNET_HOME=str(home), PYTHONNOUSERSITE="1", **extra)
    return env


def _run(home: Path, cwd: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=cwd,
        env=_env(home),
        capture_output=True,
        text=True,
        timeout=120,
        check=False,
    )


def _fasta(records: list[tuple[str, str]]) -> str:
    return "".join(
        f">{n}\n" + "".join(s[i : i + 70] + "\n" for i in range(0, len(s), 70)) for n, s in records
    )


def _pipeline(home: Path, tmp: Path, records: list[tuple[str, str]], core: int = 100) -> dict:
    """plan -> every chunk -> merge in a fresh dir: the plan rows (``rows``), the output (``tsv``)
    and each step's stderr (``plan``, ``chunk``, ``merge``)."""
    wd = tmp / f"wd{len(list(tmp.glob('wd*')))}"
    wd.mkdir()
    (wd / "genome.fa").write_text(_fasta(records))
    code = (
        "import pickle, Features; "
        "pickle.dump((Features.FakeModel(), Features.FakeScaler()), open('model.pkl', 'wb'))"
    )
    subprocess.run(
        [sys.executable, "-c", code], cwd=wd, env=_env(home, PYTHONPATH=str(home)), check=True
    )
    r = _run(home, wd, "plan", "--core", str(core))
    assert r.returncode == 0, r.stderr
    stderr = {"plan": r.stderr, "chunk": ""}
    lines = (wd / "plan.tsv").read_text().splitlines()
    plan = [dict(zip(lines[0].split("\t"), x.split("\t"), strict=True)) for x in lines[1:]]
    for row in plan:
        r = _run(home, wd, "chunk", row["chunk_id"], "--cpus", "1")
        assert r.returncode == 0, r.stderr
        stderr["chunk"] += r.stderr
    r = _run(home, wd, "merge")
    assert r.returncode == 0, r.stderr
    stderr["merge"] = r.stderr
    return {"wd": wd, "rows": plan, "tsv": (wd / "terminatornet.tsv").read_text(), **stderr}


# ---------------------------------------------------------------- independent reference


def _rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _clean(raw: str) -> str:
    if len(raw) < 4:
        return ""
    s = raw.upper().replace("U", "T")
    acgt = sum(s.count(c) for c in "ACGT")
    if acgt < 0.75 * len(s):
        return s
    return "".join(c for c in s if c in "ACGT")


def _kept_index(raw: str, cleaned: str) -> list[int]:
    """Original index of every cleaned position."""
    if len(cleaned) == len(raw):
        return list(range(len(raw)))
    k = [i for i, c in enumerate(raw.upper()) if c in "ACGTU"]
    assert len(k) == len(cleaned)
    return k


def _nwin(n: int) -> int:
    return max(n - 59, 1)


def _expected(records: list[tuple[str, str]]) -> tuple[str, int, int]:
    """(terminatornet.tsv, removed bases, dropped calls) under the #615 rule. Planted sites are at
    least SPAN apart in every layout here, so merge_overlapping keeps every call."""
    blocks: dict[str, list[str]] = {"+": [], "-": []}
    removed = dropped = 0
    for name, raw in records:
        cleaned = _clean(raw)
        lc, lo = len(cleaned), len(raw)
        k = _kept_index(raw, cleaned) if cleaned else []
        if cleaned:
            removed += lo - lc
        for strand in "+-":
            seq = cleaned if strand == "+" else _rc(cleaned)
            terms = []
            for j in range(_nwin(lc)):
                if seq[j : j + 6] != MOTIF:
                    continue
                a, b = (j, j + SPAN - 1) if strand == "+" else (lc - j - SPAN, lc - 1 - j)
                if k[b] - k[a] != b - a:
                    dropped += 1
                    continue
                pred = np.float64(0.5 + (j % 7) / 100)
                terms.append((k[a] + 1, k[b] + 1, seq[j : j + SPAN], pred))
            terms.sort()
            assert all(x[1] < y[0] for x, y in itertools.pairwise(terms)), "overlap"
            if not terms:
                blocks[strand].append(f"{name}\tNO TERMINATORS FOUND\n")
            for i, (s, e, h, p) in enumerate(terms, start=1):
                a_t, u_t = ("AAAAAAAA", "TTTTTTTT")
                blocks[strand].append(
                    f"{name}\tTERM_{i}\t{a_t}\t{h}\t{u_t}\t((((....))))\t-12.5\t{p}\t{s}\t{e}\t"
                    f"{strand}\n"
                )
    return "".join(blocks["+"] + blocks["-"]), removed, dropped


def _filler(n: int, seed: int) -> str:
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACT"), size=n))


def _warnings(stderr: str) -> list[str]:
    return [ln for ln in stderr.splitlines() if "WARNING" in ln and "non-ACGT" in ln]


def _assert_one_warning(stderr: str, removed: int, dropped: int) -> None:
    got = _warnings(stderr)
    assert len(got) == 1, f"want exactly one non-ACGT WARNING line, got {got}\n{stderr}"
    assert re.search(rf"\b{removed} non-ACGT\b", got[0]), (removed, got[0])
    assert re.search(rf"\b{dropped} calls?\b", got[0]), (dropped, got[0])


def _assert_lines_match_genome(tsv: str, records: list[tuple[str, str]]) -> None:
    raw = dict(records)
    for line in tsv.splitlines():
        f = line.split("\t")
        if len(f) != 11:
            continue
        s, e = int(f[8]), int(f[9])
        seg = raw[f[0]][s - 1 : e].upper()
        assert set(seg) <= set("ACGT"), f"{line}: span holds a removed base: {seg}"
        assert (seg if f[10] == "+" else _rc(seg)) == f[3], f"{line}: not at {s}..{e} of {f[0]}"


F = _filler

# Hand example: a 3-N gap between two sites. The + site after the gap sits at original 120..131;
# its - partner covers 114..125 (the 12 nt end at the motif's last base). The first site's -
# window is past the scanned windows of the reverse strand (j >= nwin), so it has no - call.
HAND = F(40, 1) + MOTIF + F(40, 2) + "NNN" + F(30, 3) + MOTIF + F(100, 4)

LAYOUTS: dict[str, list[tuple[str, str]]] = {
    "hand": [("chr", HAND)],
    "gap-at-0": [("chr", "N" + F(70, 5) + MOTIF + F(30, 6) + MOTIF + F(80, 7))],
    "gap-at-last": [("chr", F(70, 8) + MOTIF + F(30, 9) + MOTIF + F(80, 10) + "N")],
    "both-ends-lowercase": [("chr", "nn" + F(70, 11) + MOTIF.lower() + F(90, 12) + "nnn")],
    "several-runs-iupac": [
        (
            "chr",
            F(65, 13)
            + "NNNNN"
            + MOTIF
            + F(20, 14)
            + "R"
            + F(25, 15)
            + MOTIF
            + F(14, 16)
            + "nY"
            + MOTIF
            + F(30, 17)
            + "N" * 9
            + F(20, 18)
            + MOTIF
            + F(90, 19),
        )
    ],
    "crossing-dropped": [
        # site split by an N (exists only in the cleaned sequence): both strands dropped; site whose
        # + span reaches the NN (dropped) but whose - span ends before it (kept)
        ("chr", F(80, 20) + "GGA" + "N" + "TCC" + F(40, 21) + MOTIF + "AC" + "NN" + F(100, 22)),
    ],
    "multi-record": [
        ("gap_free", F(90, 23) + MOTIF + F(80, 24)),
        ("gapped", F(75, 25) + "NN" + MOTIF + F(30, 26) + "N" + MOTIF + F(90, 27)),
        ("mostly_n", F(30, 28) + "N" * 20 + MOTIF + F(30, 29) + "N" * 15 + F(10, 30)),
    ],
}


# ---------------------------------------------------------------- rule 1-3: merge output


def test_hand_example_coordinates_after_a_three_n_gap(home: Path, tmp_path: Path) -> None:
    out = _pipeline(home, tmp_path, LAYOUTS["hand"])
    got = [
        (f[8], f[9], f[10])
        for f in (ln.split("\t") for ln in out["tsv"].splitlines())
        if len(f) == 11
    ]
    assert got == [("41", "52", "+"), ("120", "131", "+"), ("114", "125", "-")], (
        "before the fix the + site after the gap is at 117..128 and its - partner at 111..122"
    )


@pytest.mark.parametrize("layout", sorted(LAYOUTS))
def test_calls_are_in_original_coordinates_on_both_strands(
    home: Path, tmp_path: Path, layout: str
) -> None:
    records = LAYOUTS[layout]
    out = _pipeline(home, tmp_path, records)
    want, removed, dropped = _expected(records)
    _assert_lines_match_genome(out["tsv"], records)
    assert out["tsv"] == want
    assert removed > 0
    _assert_one_warning(out["merge"], removed, dropped)
    assert not _warnings(out["plan"]) and not _warnings(out["chunk"]), (
        "only merge reports removed bases (one WARNING per run)"
    )


def test_crossing_calls_are_dropped_and_counted(home: Path, tmp_path: Path) -> None:
    records = LAYOUTS["crossing-dropped"]
    out = _pipeline(home, tmp_path, records)
    _, removed, dropped = _expected(records)
    assert (removed, dropped) == (3, 3), "reference: split site on both strands + one + span"
    got = [ln.split("\t")[8:11] for ln in out["tsv"].splitlines() if ln.count("\t") == 10]
    site = 80 + 7 + 40  # 0-based start of the intact site
    assert got == [[str(site - 6 + 1), str(site + 6), "-"]], got
    _assert_one_warning(out["merge"], 3, 3)


def test_mostly_n_record_left_uncleaned_by_upstream_is_not_remapped(
    home: Path, tmp_path: Path
) -> None:
    rec = LAYOUTS["multi-record"][2]
    out = _pipeline(home, tmp_path, [rec])
    want, removed, dropped = _expected([rec])
    assert removed == 0 and dropped == 0
    assert out["tsv"] == want
    assert "\t51\t62\t+" in want, "the site sits after 20 N at its raw position"
    assert not _warnings(out["merge"]), "upstream removed nothing: no non-ACGT WARNING"


# ---------------------------------------------------------------- rule 4: plan rule


def test_plan_seq_length_is_the_cleaned_scan_length(home: Path, tmp_path: Path) -> None:
    records = LAYOUTS["multi-record"]
    out = _pipeline(home, tmp_path, records)
    sh = _shard()
    for row in out["rows"]:
        raw = records[int(row["seq_index"])][1]
        assert row["seq_name"] == records[int(row["seq_index"])][0]
        assert int(row["seq_length"]) == len(_clean(raw)), "plan stays in scan coordinates"
    for i, (_n, raw) in enumerate(records):
        for strand in "+-":
            ends = [int(r["core_end"]) for r in out["rows"]
                    if (int(r["seq_index"]), r["strand"]) == (i, strand)]  # fmt: skip
            assert max(ends) == sh.nwin(len(_clean(raw)))


# ---------------------------------------------------------------- rule 5: gap-free regression

GAP_FREE = [
    ("chr", F(100, 31) + MOTIF + F(50, 32) + MOTIF + F(20, 33) + MOTIF + F(600, 34)),
    ("pl", F(70, 35) + MOTIF.lower() + F(120, 36)),
    ("tiny", "AC"),
]
GAP_FREE_SHA256 = "ef9e7929f47e3f1bcb292d5cb4b5014e686f83fc1406e811d4cb74411f30bb97"


def test_gap_free_genome_output_is_byte_identical(home: Path, tmp_path: Path) -> None:
    out = _pipeline(home, tmp_path, GAP_FREE, core=150)
    assert out["tsv"] == _expected(GAP_FREE)[0]
    assert hashlib.sha256(out["tsv"].encode()).hexdigest() == GAP_FREE_SHA256
    assert not _warnings(out["merge"] + out["plan"] + out["chunk"])


# ---------------------------------------------------------------- the pure seam


def _call(s: int, e: int, tag: object) -> tuple:
    return (s, e, "AAAAAAAA", f"h{s}", "TTTTTTTT", "((..))", -10.0, np.float64(0.7))


def test_to_original_hand_example() -> None:
    sh = _shard()
    raw = "ACGTNNNacgtRACGT"  # removed: 4, 5, 6, 11
    cleaned = _clean(raw)
    assert cleaned == "ACGTACGTACGT"
    plus = [_call(0, 3, 0), _call(2, 5, 1), _call(4, 7, 2), _call(8, 11, 3), _call(5, 6, 4)]
    kept, n = sh.to_original(plus, raw, cleaned, "+")
    assert n == 1, "only (2, 5) spans the removed run 4..6"
    assert [c[:2] for c in kept] == [(0, 3), (7, 10), (12, 15), (8, 9)]
    assert [c[3] for c in kept] == ["h0", "h4", "h8", "h5"], "other fields unchanged, input order"
    # - strand: RC(cleaned) index q <-> forward cleaned 11 - q; RC(raw) index <-> forward 15 - x
    minus = [_call(0, 3, 0), _call(3, 4, 1), _call(8, 11, 2)]
    kept, n = sh.to_original(minus, raw, cleaned, "-")
    assert n == 1, "(3, 4) = forward cleaned 7..8 spans the removed R at 11"
    assert [c[:2] for c in kept] == [(0, 3), (12, 15)]


def test_to_original_identity_without_removed_bases() -> None:
    sh = _shard()
    calls = [_call(0, 3, 0), _call(5, 9, 1)]
    for raw in ("ACGTACGTACGT", "acgtacgtACGU"):
        for strand in "+-":
            kept, n = sh.to_original(calls, raw, _clean(raw), strand)
            assert n == 0 and kept == calls
            assert all(
                a is b for x, y in zip(kept, calls, strict=True) for a, b in zip(x, y, strict=True)
            )
    raw = "NNNNNNNNACGT"  # < 75 % ACGT: upstream leaves it as is
    kept, n = sh.to_original(calls, raw, _clean(raw), "+")
    assert n == 0 and kept == calls


_BASES = st.sampled_from(list("ACGTACGTACGTacgtNnRYU"))


@settings(max_examples=200, deadline=None)
@given(
    raw=st.lists(_BASES, min_size=8, max_size=80).map("".join),
    spans=st.lists(st.tuples(st.integers(0, 79), st.integers(0, 11)), max_size=12),
    strand=st.sampled_from("+-"),
)
def test_to_original_property(raw: str, spans: list[tuple[int, int]], strand: str) -> None:
    sh = _shard()
    cleaned = _clean(raw)
    lc, lo = len(cleaned), len(raw)
    if lc == 0:
        return
    k = _kept_index(raw, cleaned)
    calls = [_call(s % lc, min(s % lc + w, lc - 1), i) for i, (s, w) in enumerate(spans)]
    kept, n = sh.to_original(calls, raw, cleaned, strand)
    want = []
    for c in calls:
        s, e = c[0], c[1]
        a, b = (s, e) if strand == "+" else (lc - 1 - e, lc - 1 - s)  # forward cleaned
        if k[b] - k[a] != b - a:
            continue
        oa, ob = k[a], k[b]  # forward original
        want.append((oa, ob, *c[2:]) if strand == "+" else (lo - 1 - ob, lo - 1 - oa, *c[2:]))
    assert kept == want
    assert n == len(calls) - len(want)
    fwd_raw = raw.upper().replace("U", "T")
    seq_raw = fwd_raw if strand == "+" else _rc(fwd_raw)
    seq_cln = cleaned if strand == "+" else _rc(cleaned)
    for c in kept:
        if lc != lo:
            assert set(seq_raw[c[0] : c[1] + 1]) <= set("ACGT"), "kept span holds a removed base"
    for c in kept:
        src = next(x for x in calls if x[3:] == c[3:] and x[1] - x[0] == c[1] - c[0])
        assert seq_raw[c[0] : c[1] + 1] == seq_cln[src[0] : src[1] + 1], "content preserved"
