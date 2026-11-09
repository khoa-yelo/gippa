"""Shared fixtures of the M12-05a tests: V13 draws every novel ORF candidate, lowest q first
(docs/specs/M12-05a-novel-orf-plot.md, § 7 examples E1-E10 and the lead decisions of 2026-11-09).

Interface pinned here (PR A test author, beyond the spec text):

- the display order (spec § 2 with lead decision Q1, "class q, then peptide q") is the stable
  sort by ``(class_fdr, best_peptide_q, not called, rank, rank_key)``: ``class_fdr`` read as
  ``_fdr`` (NA, ``unestimated`` and non-finite values sort after every finite value),
  ``best_peptide_q`` ascending with NA last, then called before not called, ``rank`` ascending
  with NA last, then ``gippa.novel.evidence.rank_key``. ``expected_order`` is that oracle;
- ``ORF_KEY_CALLED`` / ``ORF_KEY_UNCALLED`` are the spec § 3 key texts, also importable as
  ``gippa.viz.locus.ORF_KEY_CALLED`` / ``ORF_KEY_UNCALLED``;
- the raw class identifiers that must not appear in any V13 SVG text (lead decision, added
  from the M12-05c review) are the schema enum values that are not plain English words:
  ``nterm_extension``, ``cterm_extension``, ``cds_outframe`` and ``sorf``.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping, Sequence
from typing import Any

from gippa.novel.evidence import rank_key

#: spec § 6.1
INDEX_COLUMNS = [
    "order", "orf_id", "class", "called", "class_fdr", "best_peptide_q", "n_peptides", "rank",
    "rendered",
]  # fmt: skip
TOP_N = 20
ORF_EDGE = "#CC79A7"
ORF_KEY_CALLED = "called novel ORF; tick: OSTIR start; purple pileup: its peptides"
ORF_KEY_UNCALLED = (
    "novel ORF candidate, not called (dashed); tick: OSTIR start; purple pileup: its peptides"
)
INDEX_TITLE = "Novel ORF loci, lowest class FDR first"
#: spec § 5 captions
CAPTION_V13 = (
    "Which new genes or upstream extensions are supported? Panel (a) counts the candidates of "
    "each class by call status, and panel (b) draws up to 20 candidates, called or not, in their "
    "genomic neighbourhood, lowest class FDR first. A filled arrow is a called ORF, a dashed "
    "outline a candidate that was not called."
)
INDEX_CAPTION_V13 = (
    "This index lists the candidates drawn as locus views, lowest class FDR first, with their "
    "call status."
)
QUESTION_V13 = "Which new genes or upstream extensions are supported?"
#: the class enum values that are no plain words (lead decision 2026-11-09)
RAW_CLASSES = ("nterm_extension", "cterm_extension", "cds_outframe", "sorf")
_RAW = re.compile(r"(?<![A-Za-z0-9])(" + "|".join(RAW_CLASSES) + r")(?![A-Za-z0-9])")


def raw_class_tokens(text: str) -> list[str]:
    """The raw class identifiers occurring in ``text`` as whole tokens."""
    return _RAW.findall(text)


def orf(
    orf_id: str,
    *,
    cls: str = "intergenic",
    contig: str = "chr1",
    strand: str = "+",
    start_longest: int | None = 5_001,
    end: int = 5_300,
    start_best_ostir: int | None = None,
    called: bool | None = True,
    rank: int | None = 1,
    class_fdr: str | None = "0.004",
    n_peptides: int = 2,
    peptides: str | None = "PEPTWDEK;SAMPQER",
    best_q: float | None = 0.0001,
    **kw: object,
) -> dict[str, object]:
    """One ``novel_orfs.tsv`` row; coordinates 1-based, in transcription direction (as in
    ``tests/unit/viz/test_v13_novel_orfs.py``)."""
    row: dict[str, object] = {
        "orf_id": orf_id,
        "class": cls,
        "contig": contig,
        "strand": strand,
        "frame": 1 if strand == "+" else -1,
        "start_longest": start_longest,
        "start_best_ostir": start_longest if start_best_ostir is None else start_best_ostir,
        "end": end,
        "length_aa_longest": 99,
        "length_aa_best": 99,
        "start_codon_best": "ATG",
        "n_peptides": n_peptides,
        "peptides": peptides,
        "best_peptide_q": best_q,
        "tir_ostir": 1000.0,
        "promoter_score": 0.9,
        "promoter_supported": True,
        "terminator_downstream": True,
        "terminator_conf": 80.0,
        "terminator_distance_nt": 30,
        "rbs_supported": True,
        "class_fdr": class_fdr,
        "called": called,
        "rank": rank,
        "overlapping_gene": None,
        "flags": None,
    }
    row.update(kw)
    return row


# --- E1: route E of the 1.1.0rc1 example bundle (8 rows, 0 called) -----------------------------

ROUTE_E_PREFIX = "NOVEL_nissle1917_NZ_CP058217_1_1400000_1600000_"


def _e(suffix: str) -> str:
    return ROUTE_E_PREFIX + suffix


#: (suffix, class, strand, start_longest, start_best_ostir, end, start codon, best q, class FDR,
#: TIR, promoter, length aa) from ``runs/E/tables/novel_orfs.tsv``; positions folded onto the
#: 100 kb ``chr1`` of ``helpers.viz_locus`` (mod 100 kb, the ORF IDs keep the run's names)
_ROUTE_E: tuple[tuple[Any, ...], ...] = (
    ("20324_m", "nterm_extension", "-", 22_390, 22_390, 20_324, "TTG", 0.000874751, "1",
     39.3629, 0.603813, 688),
    ("20948_p", "antisense", "+", 20_948, 22_049, 22_573, "ATG", 0.00720727, "0.4", 1325.03,
     0.489257, 174),
    ("25949_p", "intergenic", "+", 25_949, 26_018, 26_059, "ATG", 0.00720728, "1", 77.4018,
     0.667635, 13),
    ("43000_m", "antisense", "-", None, None, 43_000, None, 0.00509304, "0.4", None, None, 21),
    ("58745_p", "antisense", "+", 58_745, 58_745, 59_101, "ATG", 0.0024706, "0.4", 98.247,
     0.541144, 118),
    ("130563_m", "antisense", "-", 30_847, 30_730, 30_563, "TTG", 0.000874751, "0.4", 4.7927,
     0.651464, 55),
    ("149890_p", "cds_outframe", "+", None, None, 49_943, None, 0.00779485, "1", None, None, 17),
    ("195190_m", "antisense", "-", 96_635, 96_407, 95_190, "ATG", 0.00339301, "0.4", 146.571,
     0.770085, 405),
)  # fmt: skip

#: spec E1, the expected display order
E1_ORDER = [
    _e(s)
    for s in ("130563_m", "58745_p", "195190_m", "43000_m", "20948_p", "20324_m", "25949_p",
              "149890_p")
]  # fmt: skip


def route_e_rows() -> list[dict[str, object]]:
    """E1: the 8 route E candidates, none called, ``rank`` NA, one peptide each."""
    rows = []
    for sfx, cls, strand, first, ostir, end, codon, q, fdr, tir, prom, length in _ROUTE_E:
        flags = "no_start_codon" if first is None else None
        rows.append(orf(_e(sfx), cls=cls, strand=strand, start_longest=first,
                        start_best_ostir=ostir, end=end, called=False, rank=None, class_fdr=fdr,
                        n_peptides=1, peptides=None, best_q=q, start_codon_best=codon,
                        tir_ostir=tir, promoter_score=prom, length_aa_best=length,
                        length_aa_longest=length, terminator_conf=None,
                        terminator_downstream=False, terminator_distance_nt=None, flags=flags,
                        promoter_supported=prom is not None and prom > 0.6,
                        rbs_supported=None if tir is None else tir > 1000))  # fmt: skip
    # the table's own row order is by orf_id string, not the display order
    return sorted(rows, key=lambda r: str(r["orf_id"]))


# --- E4-E9 ----------------------------------------------------------------------------------------


def e4_rows() -> tuple[list[dict[str, object]], list[str]]:
    """Tie on q, called vs not: A (called, 2 peptides, rank 1) before B (not called, 5)."""
    a = orf("NOVEL_chr1_9001_m", cls="antisense", strand="-", start_longest=9_300, end=9_001,
            class_fdr="0.0065", called=True, n_peptides=2, rank=1, best_q=0.001)  # fmt: skip
    b = orf("NOVEL_chr1_1001_m", cls="antisense", strand="-", start_longest=1_300, end=1_001,
            class_fdr="0.0065", called=False, n_peptides=5, rank=None, best_q=0.001)  # fmt: skip
    return [b, a], [str(a["orf_id"]), str(b["orf_id"])]


def e5_rows() -> tuple[list[dict[str, object]], list[str]]:
    """Tie on q and call, both uncalled: D (best q 0.0009) before C (0.0025)."""
    c = orf("NOVEL_chr1_1001_p", start_longest=1_001, end=1_300, class_fdr="0.4", called=False,
            rank=None, n_peptides=1, best_q=0.0025)  # fmt: skip
    d = orf("NOVEL_chr1_9001_p", start_longest=9_001, end=9_300, class_fdr="0.4", called=False,
            rank=None, n_peptides=1, best_q=0.0009)  # fmt: skip
    return [c, d], [str(d["orf_id"]), str(c["orf_id"])]


def e6_rows() -> tuple[list[dict[str, object]], list[str]]:
    """Full tie: the string order of ``orf_id`` decides (``1000`` before ``900``)."""
    x = orf("NOVEL_chr1_900_p", start_longest=900, end=1_199, class_fdr="0.4", called=False,
            rank=None, n_peptides=1)  # fmt: skip
    y = orf("NOVEL_chr1_1000_p", start_longest=1_000, end=1_299, class_fdr="0.4", called=False,
            rank=None, n_peptides=1)  # fmt: skip
    return [x, y], ["NOVEL_chr1_1000_p", "NOVEL_chr1_900_p"]


E7_E, E7_F, E7_G = "NOVEL_chr1_5001_p", "NOVEL_chr1_9001_p", "NOVEL_chr1_1001_p"


def e7_rows() -> tuple[list[dict[str, object]], list[str]]:
    """NaN q: F (finite q) first; then among "no q" the called E before the uncalled G."""
    e = orf(E7_E, start_longest=5_001, end=5_300, class_fdr=None, called=True, rank=1)
    f = orf(E7_F, start_longest=9_001, end=9_300, class_fdr="0.9", called=False, rank=None)
    g = orf(E7_G, start_longest=1_001, end=1_300, class_fdr="unestimated", called=False,
            rank=None)  # fmt: skip
    return [g, e, f], [E7_F, E7_E, E7_G]


def e8_rows() -> list[dict[str, object]]:
    """25 rows, none called, class FDR 0.01 + 0.01 i (i = 0 ... 24), listed in reverse."""
    rows = []
    for i in range(25):
        s = 1_001 + 3_000 * i
        rows.append(orf(f"NOVEL_chr1_{s}_p", start_longest=s, end=s + 299,
                        class_fdr=f"{0.01 + 0.01 * i:.2f}", called=False, rank=None))  # fmt: skip
    return rows[::-1]


#: E9: a 10 kb circular contig
C10 = ({"id": "c10", "kind": "chromosome", "length": 10_000, "topology": "circular"},)
E9_ID = "NOVEL_c10_9900_m"


def e9_row(called: bool = False) -> dict[str, object]:
    """One antisense row on ``-`` crossing the origin of ``c10`` (``start_longest = 150``)."""
    return orf(E9_ID, cls="antisense", contig="c10", strand="-", start_longest=150, end=9_900,
               called=called, rank=1 if called else None, class_fdr="0.4", n_peptides=1,
               peptides="PEPTWDEK", best_q=0.002)  # fmt: skip


# --- the order oracle ---------------------------------------------------------------------------


def _finite(v: object) -> float:
    """A finite float, else +inf (NA, ``unestimated``, NaN and infinities)."""
    try:
        f = float(v)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return math.inf
    return f if math.isfinite(f) else math.inf


def is_called(v: object) -> bool:
    return v is not None and not (isinstance(v, float) and math.isnan(v)) and bool(v)


def order_key(row: Mapping[str, object]) -> tuple[object, ...]:
    rank = row.get("rank")
    return (
        _finite(row.get("class_fdr")),
        _finite(row.get("best_peptide_q")),
        not is_called(row.get("called")),
        math.inf if rank is None else float(rank),  # type: ignore[arg-type]
        rank_key(row),
    )


def expected_order(rows: Sequence[Mapping[str, object]]) -> list[str]:
    return [str(r["orf_id"]) for r in sorted(rows, key=order_key)]
