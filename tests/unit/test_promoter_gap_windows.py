"""Promoter windows over assembly gaps are not sent to Promotech (issue #611, 1.1.1).

Promotech's one-hot encoder drops non-ACGT bases, so a 40-nt window with an ``N`` (or any IUPAC
code) gives a vector shorter than 160 and ANNOTATE:PROMOTECH crashes. Lead rule (2026-11-22):

1. a window with any non-ACGT base (case-insensitive) is not written to the Promotech chunks;
2. ``promoter_layer`` gives it no score (NaN), it is never called a promoter, it is counted and
   the run logs one WARNING with the count;
3. ``promoter_layer`` still raises for a missing score of a window that was sent to Promotech;
4. a genome with only ACGT gives byte-identical outputs.

Placement pinned here (follows the code): ``dedupe_and_chunk`` (shared by ``promoter_windows``
and the novel ``window_chunks``) leaves gap windows out of the chunks but keeps one row per window
in the map, with the same columns; the map must tell ``promoter_layer(window_map, scores)`` which
rows were skipped (signature unchanged), so the step needs no new input. The WARNING count is the
number of skipped window positions (genes here never share positions, so map rows = positions).
"""

from __future__ import annotations

import hashlib
import importlib
import logging
import random
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.annotate.promoters import (
    Window,
    dedupe_and_chunk,
    gene_windows,
    promoter_layer,
    promoters_bed,
)
from gippa.config import GippaConfig
from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene
from gippa.novel.candidates import window_chunks

MAP_COLUMNS = ["header", "contig", "strand", "end3p", "owner", "k"]
_ACGT = set("ACGT")


def _is_gap(seq: str) -> bool:
    return not set(seq.upper()) <= _ACGT


def _win(end3p: int, seq: str, owner: str = "G", strand: str = "+") -> Window:
    assert len(seq) == 40
    return Window("chr1", strand, end3p, seq, owner, end3p)


CLEAN = "ACGTTGCAAC" * 4
WINDOWS = [
    _win(100, CLEAN),
    _win(101, "ACGTACGTACGTACGTNACGTACGTACGTACGTACGTACG"),  # one N
    _win(102, "acgtacgtacgtacgtnacgtacgtacgtacgtacgtacg"),  # lowercase n
    _win(103, "ACGTACGTACGTACGTRACGTACGTACGTACGTACGTACG"),  # IUPAC R
    _win(104, "N" * 40),  # all N
    _win(105, "acgtacgtacgtacgtacgtacgtacgtacgtacgtacgt"),  # soft-masked, all ACGT: kept
    _win(106, "GGGCCCAAAT" * 4),
    _win(107, "N" * 40, owner="H"),  # a second all-N window (same sequence)
]
GAP_END3P = {101, 102, 103, 104, 107}


def _chunk_seqs(chunks: list[list[tuple[str, str]]]) -> dict[str, str]:
    return {h: s for c in chunks for h, s in c}


def _scores_for(chunks: list[list[tuple[str, str]]]) -> dict[str, float]:
    return {h: (i % 9 + 1) / 10 for i, h in enumerate(_chunk_seqs(chunks))}


# ------------------------------------------------------------------ rule 1: chunks are ACGT only


def test_gap_windows_are_not_written_to_promotech_chunks() -> None:
    chunks, _ = dedupe_and_chunk(WINDOWS)
    sent = sorted(_chunk_seqs(chunks).values())
    want = sorted({w.seq for w in WINDOWS if not _is_gap(w.seq)})
    assert sent == want, "only ACGT-only windows (any case) go to Promotech"


@pytest.mark.parametrize(
    "seq",
    [WINDOWS[1].seq, WINDOWS[2].seq, WINDOWS[3].seq, "N" * 40],
    ids=["N", "lowercase-n", "IUPAC-R", "all-N"],
)
def test_each_gap_kind_is_skipped(seq: str) -> None:
    chunks, wmap = dedupe_and_chunk([_win(100, CLEAN), _win(200, seq)])
    assert list(_chunk_seqs(chunks).values()) == [CLEAN]
    assert len(wmap) == 2, "the skipped window keeps its map row (its position is known)"
    layer = promoter_layer(wmap, _scores_for(chunks))
    got = dict(zip(layer["end3p"], layer["score"], strict=True))
    assert np.isnan(got[200])
    assert got[100] == pytest.approx(0.1)


def test_novel_window_chunks_skip_gap_windows_too() -> None:
    chunks, wmap = window_chunks(WINDOWS)
    assert not any(_is_gap(s) for s in _chunk_seqs(chunks).values())
    assert len(wmap) == len(WINDOWS)
    assert list(wmap["orf_id"]) == list(wmap["owner"])


def test_map_keeps_every_window_row_and_columns() -> None:
    _, wmap = dedupe_and_chunk(WINDOWS)
    assert list(wmap.columns) == MAP_COLUMNS
    assert list(wmap["end3p"]) == [w.end3p for w in WINDOWS]
    assert list(wmap["owner"]) == [w.owner for w in WINDOWS]


# --------------------------------------------------------------- rule 2: NaN, never a promoter


def test_layer_gives_nan_to_gap_windows_and_scores_the_rest() -> None:
    chunks, wmap = dedupe_and_chunk(WINDOWS)
    scores = _scores_for(chunks)
    layer = promoter_layer(wmap, scores)
    assert list(layer.columns) == ["contig", "strand", "end3p", "score"]
    assert len(layer) == len(WINDOWS)
    by_pos = dict(zip(layer["end3p"], layer["score"], strict=True))
    assert {p for p, s in by_pos.items() if np.isnan(s)} == GAP_END3P
    header_of = {s: h for h, s in _chunk_seqs(chunks).items()}
    for w in WINDOWS:
        if w.end3p not in GAP_END3P:
            assert by_pos[w.end3p] == scores[header_of[w.seq]]


def test_gap_window_is_never_called_a_promoter() -> None:
    gap = [_win(p, "N" * 40) for p in range(100, 140)] + [_win(200, WINDOWS[3].seq)]
    chunks, wmap = dedupe_and_chunk([*gap, _win(5000, CLEAN)])
    layer = promoter_layer(wmap, dict.fromkeys(_chunk_seqs(chunks), 0.99))
    reps = {"chr1": Replicon("chr1", 10_000, "linear", "chromosome")}
    beds = promoters_bed(layer, reps, 0.5)
    assert [(b.interval.start, b.interval.end) for b in beds] == [(4960, 5000)]


# ------------------------------------------------------------ rule 3: other missing scores raise


def test_missing_score_for_a_sent_window_still_raises() -> None:
    chunks, wmap = dedupe_and_chunk(WINDOWS)
    scores = _scores_for(chunks)
    clean_header = next(h for h, s in _chunk_seqs(chunks).items() if s == CLEAN)
    del scores[clean_header]
    with pytest.raises(ValueError, match="no Promotech score"):
        promoter_layer(wmap, scores)


def test_missing_score_raises_on_a_gap_free_map() -> None:
    chunks, wmap = dedupe_and_chunk([_win(100, CLEAN), _win(101, "GGGCCCAAAT" * 4)])
    scores = _scores_for(chunks)
    scores.pop(next(iter(scores)))
    with pytest.raises(ValueError, match="no Promotech score"):
        promoter_layer(wmap, scores)


# ------------------------------------------------------------------------ property (invariant)


_BASES = st.sampled_from(list("ACGTacgtNnRYKMSWBDHV"))


@settings(max_examples=60, deadline=None)
@given(
    st.lists(
        st.one_of(
            st.text(alphabet="ACGT", min_size=40, max_size=40),
            st.text(alphabet=_BASES, min_size=40, max_size=40),
        ),
        min_size=1,
        max_size=25,
    )
)
def test_property_gap_windows_are_nan_and_clean_ones_scored(seqs: list[str]) -> None:
    ws = [_win(1000 + i, s) for i, s in enumerate(seqs)]
    chunks, wmap = dedupe_and_chunk(ws, chunk_size=4)
    sent = _chunk_seqs(chunks)
    assert not any(_is_gap(s) for s in sent.values())
    assert set(sent.values()) == {s for s in seqs if not _is_gap(s)}
    assert all(len(s) == 40 for s in sent.values())
    layer = promoter_layer(wmap, _scores_for(chunks))
    by_pos = dict(zip(layer["end3p"], layer["score"], strict=True))
    for w in ws:
        assert np.isnan(by_pos[w.end3p]) == _is_gap(w.seq)


# ------------------------------------------------------- step level: promoter_windows → layer


@pytest.fixture
def reg() -> Iterator[ModuleType]:
    from gippa.steps import load_step_modules, registry

    with pytest.MonkeyPatch.context() as mp:
        names = load_step_modules()
        mp.setattr(registry, "_REGISTRY", {})
        for n in names:
            if n != "registry":
                importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
        yield registry


def _gene(tag: str, start: int, end: int, strand: str, contig: str) -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


GENES = [
    _gene("A", 1000, 1300, "+", "c1"),
    _gene("B", 2000, 2300, "-", "c1"),
    _gene("C", 1000, 1300, "+", "p1"),
]
REPS = {
    "c1": Replicon("c1", 3000, "circular", "chromosome"),
    "p1": Replicon("p1", 2000, "linear", "plasmid"),
}


def _genome(gaps: bool) -> dict[str, str]:
    rng = random.Random(611)
    c1 = list("".join(rng.choice("ACGT") for _ in range(3000)))
    p1 = list("".join(rng.choice("ACGT") for _ in range(2000)))
    if gaps:
        c1[900] = "N"  # 40 of A's upstream windows (k = -99 .. -60)
        c1[2400] = "r"  # lowercase IUPAC in 40 of B's windows
        p1[850:910] = list("N" * 30 + "n" * 30)  # a 60-nt gap: 99 of C's windows, 21 all-gap
    return {"c1": "".join(c1), "p1": "".join(p1)}


def _write_inputs(d: Path, genome: dict[str, str]) -> dict[str, Path]:
    d.mkdir(parents=True)
    (d / "genome.fa").write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    pq.write_table(
        pa.table({"contig": ["c1", "p1"], "length": [3000, 2000],
                  "topology": ["circular", "linear"], "kind": ["chromosome", "plasmid"]}),
        d / "replicons.parquet",
    )  # fmt: skip
    cols = {
        "locus_tag": [g.locus_tag for g in GENES],
        "feature_type": [g.feature_type for g in GENES],
        "contig": [g.interval.contig for g in GENES],
        "start": [g.interval.start for g in GENES],
        "end": [g.interval.end for g in GENES],
        "strand": [g.interval.strand for g in GENES],
        "aa": [g.aa for g in GENES],
        "nt": [g.nt for g in GENES],
        "start_codon": [g.start_codon for g in GENES],
    }
    pq.write_table(pa.table(cols), d / "genes.parquet")
    return {"genes": d / "genes.parquet", "replicons": d / "replicons.parquet",
            "genome": d / "genome.fa"}  # fmt: skip


def _fake_promotech(windows_dir: Path, out: Path) -> None:
    """Score every chunk record the way Promotech writes it (tab-separated CHROM/SEQ/PRED)."""
    for fa in sorted(windows_dir.glob("windows_*.fa")):
        recs = [r.split("\n", 1) for r in fa.read_text().split(">")[1:]]
        d = out / f"promotech_{fa.stem.split('_')[1]}"
        d.mkdir(parents=True)
        lines = ["CHROM\tSEQ\tPRED"]
        for h, body in recs:
            seq = body.replace("\n", "")
            assert len(seq) == 40 and set(seq.upper()) <= _ACGT, f"Promotech would crash on {h}"
            lines.append(f"{h}\t{seq}\t{(int(h[1:]) % 97) / 97}")
        (d / "sequences_predictions.csv").write_text("\n".join(lines) + "\n")


def _run(reg: ModuleType, tmp: Path, genome: dict[str, str]) -> dict[str, Path]:
    ins = _write_inputs(tmp / "in", genome)
    out = tmp / "out"
    outs = {"windows": out / "windows", "windows_map": out / "windows_map.parquet"}
    reg.run_step("promoter_windows", GippaConfig(seed=1), ins, outs)
    _fake_promotech(outs["windows"], out / "promotech")
    layer = out / "promoter_windows.parquet"
    reg.run_step(
        "promoter_layer",
        GippaConfig(seed=1),
        {"windows_map": outs["windows_map"], "promotech": out / "promotech"},
        {"promoter_windows": layer},
    )
    return {**outs, "layer": layer}


def _expected_gap_positions(genome: dict[str, str]) -> set[tuple[str, str, int]]:
    return {
        (w.contig, w.strand, w.end3p)
        for g in GENES
        for w in gene_windows(g, REPS[g.interval.contig], genome[g.interval.contig])
        if _is_gap(w.seq)
    }


def _gap_warnings(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [
        r
        for r in caplog.records
        if r.levelno == logging.WARNING and r.name.startswith("gippa.steps.promoter")
    ]


def test_step_skips_gap_windows_and_warns_once_with_the_count(
    tmp_path: Path, reg: ModuleType, caplog: pytest.LogCaptureFixture
) -> None:
    genome = _genome(gaps=True)
    want_gaps = _expected_gap_positions(genome)
    assert len(want_gaps) == 40 + 40 + 99
    with caplog.at_level(logging.INFO, logger="gippa.steps"):
        outs = _run(reg, tmp_path, genome)  # _fake_promotech asserts the chunks are ACGT only

    layer = pd.read_parquet(outs["layer"])
    assert len(layer) == 221 + 221 + 221
    na = {(r.contig, r.strand, r.end3p) for r in layer.itertuples(index=False) if np.isnan(r.score)}
    assert na == want_gaps
    assert layer.loc[~layer["score"].isna(), "score"].between(0.0, 1.0).all()

    warns = _gap_warnings(caplog)
    assert len(warns) == 1, f"one WARNING for the skipped windows, got {[r.message for r in warns]}"
    assert str(len(want_gaps)) in warns[0].getMessage()


def test_step_on_gap_free_genome_logs_no_gap_warning(
    tmp_path: Path, reg: ModuleType, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger="gippa.steps"):
        outs = _run(reg, tmp_path, _genome(gaps=False))
    assert _gap_warnings(caplog) == []
    assert not pd.read_parquet(outs["layer"])["score"].isna().any()


# Recorded on main fd886bfe (before the fix): a gap-free genome must stay byte-identical.
GAP_FREE_SHA256 = {
    "windows_000.fa": "4196df95f6f5d762473ab40ad74e909d4b036b12ca0554af234b07c973b8a351",
    "windows_map.parquet": "81efb3a0d07b2f43aa9d31fa66290fa75ed1b3b0ca66374271ba18dd717c0195",
    "promoter_windows.parquet": "f843be028828c751e654339dc0b25834f9c60bd20e84ffef4c52b46937a45eba",
}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def test_gap_free_genome_outputs_are_byte_identical(tmp_path: Path, reg: ModuleType) -> None:
    outs = _run(reg, tmp_path, _genome(gaps=False))
    got = {
        "windows_000.fa": _sha(outs["windows"] / "windows_000.fa"),
        "windows_map.parquet": _sha(outs["windows_map"]),
        "promoter_windows.parquet": _sha(outs["layer"]),
    }
    assert sorted(p.name for p in outs["windows"].iterdir()) == ["windows_000.fa"]
    assert got == GAP_FREE_SHA256
