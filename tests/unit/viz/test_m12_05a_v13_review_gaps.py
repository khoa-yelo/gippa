"""M12-05a review round 1, test gaps (docs/specs/M12-05a-novel-orf-plot.md § 2 with lead
decision Q1, § 6.1 and AC 5).

- lead decision Q1 (Khoa 2026-11-09, "class q, then peptide q"): among equal ``class_fdr`` the
  lower ``best_peptide_q`` comes first, before "called first" is looked at. A fixed example, so
  that swapping the second and third sort steps fails a plain unit test (until now only a
  Hypothesis property caught it);
- the index title at the cap boundary: exactly ``TOP_N = 20`` candidates is "(20 candidates)",
  21 is "(top 20 of 21; all in tables/novel_orfs.tsv)".
"""

from __future__ import annotations

import importlib
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pytest
from helpers import m12_05a as h
from helpers import viz_locus as vl

import gippa.viz.locus as locus_mod
from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip

LOCUS_COLUMNS = ["track", "id", "contig", "start", "end", "strand", "value", "category", "fill"]


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


@pytest.fixture
def nov() -> ModuleType:
    return importlib.import_module("gippa.viz.novel_orfs")


class _Stub:
    """A cheap ``render_locus``: records the loci, draws nothing."""

    def __init__(self) -> None:
        self.loci: list[locus_mod.Locus] = []

    def __call__(self, lo: locus_mod.Locus, *args: Any, **kw: Any) -> Instance:
        self.loci.append(lo)
        data = pd.DataFrame({c: pd.Series([], dtype=object) for c in LOCUS_COLUMNS})
        return Instance(lo.locus_id, style.new_figure("single", 1.0), data, lo.label, ("start",))


@pytest.fixture
def stub(nov: ModuleType, monkeypatch: pytest.MonkeyPatch) -> Iterator[_Stub]:
    s = _Stub()
    monkeypatch.setattr(locus_mod, "render_locus", s)
    if hasattr(nov, "render_locus"):
        monkeypatch.setattr(nov, "render_locus", s)
    yield s


def _run(tmp: Path, rows: Sequence[dict[str, object]], cfg: GippaConfig) -> list[Instance]:
    importlib.import_module("gippa.viz.novel_orfs")
    view = RunView.from_outdir(vl.build_out(tmp, novel_orfs=rows))
    result = base.registered()["V13"].fn(view, cfg)
    assert not isinstance(result, Skip), result
    return list(result)


def _views(insts: Sequence[Instance]) -> list[str]:
    return [
        str(i.suffix).removeprefix("locus_")
        for i in insts
        if i.suffix is not None and i.suffix.startswith("locus_")
    ]


def _index(insts: Sequence[Instance]) -> Instance:
    found = [i for i in insts if i.suffix == "locus"]
    assert len(found) == 1, [i.suffix for i in insts]
    return found[0]


def _collapse(s: str) -> str:
    return " ".join(s.split())


# --- lead decision Q1: best peptide q before "called first" -------------------------------------

Q1_A = "NOVEL_chr1_5001_p"  # called, rank 1, best peptide q 0.002
Q1_B = "NOVEL_chr1_20001_p"  # not called, rank NA, best peptide q 0.001


def _q1_rows() -> list[dict[str, object]]:
    a = h.orf(Q1_A, cls="antisense", start_longest=5_001, end=5_300, class_fdr="0.0065",
              called=True, rank=1, best_q=0.002)  # fmt: skip
    b = h.orf(Q1_B, cls="antisense", start_longest=20_001, end=20_300, class_fdr="0.0065",
              called=False, rank=None, best_q=0.001)  # fmt: skip
    return [a, b]


def test_q1_fixed_example_matches_the_oracle() -> None:
    rows = _q1_rows()
    assert h.expected_order(rows) == [Q1_B, Q1_A]
    assert h.expected_order(rows[::-1]) == [Q1_B, Q1_A]


@pytest.mark.parametrize("reverse", [False, True], ids=["A_first", "B_first"])
def test_q1_lower_best_peptide_q_beats_called(
    tmp_path: Path, cfg: GippaConfig, stub: _Stub, reverse: bool
) -> None:
    rows = _q1_rows()
    insts = _run(tmp_path, rows[::-1] if reverse else rows, cfg)
    assert _views(insts) == [Q1_B, Q1_A]
    data = _index(insts).data
    assert list(data["orf_id"]) == [Q1_B, Q1_A]
    assert list(data["order"]) == [1, 2]
    assert [bool(v) for v in data["called"]] == [False, True]


@pytest.mark.parametrize("reverse", [False, True], ids=["A_first", "B_first"])
def test_q1_display_order_function(nov: ModuleType, reverse: bool) -> None:
    rows = _q1_rows()
    frame = pd.DataFrame(rows[::-1] if reverse else rows).rename(columns={"class": "cls"})
    frame["rank"] = frame["rank"].astype("Int64")
    assert list(nov.display_order(frame)["orf_id"]) == [Q1_B, Q1_A]


# --- AC 5: the index title at the cap boundary ----------------------------------------------------


@pytest.mark.parametrize(
    ("m", "suffix"),
    [
        (19, "(19 candidates)"),
        (20, "(20 candidates)"),
        (21, "(top 20 of 21; all in tables/novel_orfs.tsv)"),
    ],
)
def test_index_title_helper_at_the_cap(nov: ModuleType, m: int, suffix: str) -> None:
    assert nov._index_title(m) == f"{h.INDEX_TITLE} {suffix}"


def test_index_title_exactly_20_candidates(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    rows = h.e8_rows()[:20]
    insts = _run(tmp_path, rows, cfg)
    ix = _index(insts)
    assert _collapse(ix.title) == f"{h.INDEX_TITLE} (20 candidates)"
    assert len(ix.data) == 20
    assert [bool(v) for v in ix.data["rendered"]] == [True] * 20
    assert _views(insts) == h.expected_order(rows)


def test_index_title_21_candidates(tmp_path: Path, cfg: GippaConfig, stub: _Stub) -> None:
    rows = h.e8_rows()[:21]
    insts = _run(tmp_path, rows, cfg)
    ix = _index(insts)
    want = f"{h.INDEX_TITLE} (top 20 of 21; all in tables/novel_orfs.tsv)"
    assert _collapse(ix.title) == want
    assert len(ix.data) == 21
    assert [bool(v) for v in ix.data["rendered"]] == [True] * 20 + [False]
    assert _views(insts) == h.expected_order(rows)[:20]
