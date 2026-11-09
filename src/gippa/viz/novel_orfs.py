"""V13 ``novel_orfs``: which new genes or upstream extensions are supported? (GIPPA_SPEC.md
§10.2 V13, §7.6, §9.3; M9-12 and its lead decisions 1-10).

Panel (a) counts the candidates per ``class`` (schema enum order), stacked by status: ``called``
is authoritative (decision 5); an uncalled row whose ``class_fdr`` is NA, the literal
``unestimated`` or flagged ``class_fdr_unestimated`` is ``unestimated``; the rest are
``not_called``. Panel (b) draws the first ``TOP_N`` candidates, called or not, in the display
order (M12-05a § 2: class FDR, best peptide q, called first, ``rank``, then
``gippa.novel.evidence.rank_key``; NA last at every step) with the V3 renderer
(``gippa.viz.locus.render_locus``, ``Locus(kind="novel_orf")``), one instance each, plus an index
instance listing every candidate. A called ORF is a filled arrow, an uncalled one a dashed
outline (M12-05a § 3). ``start_longest`` and ``end`` are in transcription direction; ``RunView``
already subtracted 1 from ``start_longest`` / ``start_best_ostir``, which on ``-`` are the high
ends, so 1 is added back there (decision 6). A row crosses the origin when ``end <
start_longest`` on ``+`` or ``start_longest < end`` on ``-``; its span is unwrapped past the
replicon length (decision 7).
"""

from __future__ import annotations

import dataclasses
import math
from typing import TYPE_CHECKING, Any

import numpy as np
import pandas as pd
from matplotlib.ticker import MaxNLocator

from gippa.novel.evidence import rank_key
from gippa.viz import base, locus, style
from gippa.viz.locus import Locus, OrfMark, render_locus
from gippa.viz.words import NOVEL_CLASS_LABEL

if TYPE_CHECKING:
    from gippa.config import GippaConfig
    from gippa.io.run_view import RunView

CLASSES: tuple[str, ...] = (
    "nterm_extension", "cterm_extension", "intergenic", "cds_outframe", "antisense",
    "pseudogene", "sorf",
)  # fmt: skip
STATUSES: tuple[str, ...] = ("called", "not_called", "unestimated")
STATUS_COLOUR = {
    "called": style.DETECTION["detected"],
    "not_called": style.DETECTION["not_detected"],
    "unestimated": style.DETECTION["unestimated"],
}
COUNT_COLUMNS = ("class", "status", "n", "class_fdr_max_called")
#: the index data TSV (M12-05a § 6.1): every candidate, in display order
INDEX_COLUMNS = ("order", "orf_id", "class", "called", "class_fdr", "best_peptide_q", "n_peptides",
                 "rank", "rendered")  # fmt: skip
TOP_N = 20
INDEX_TITLE = "Novel ORF loci, lowest class FDR first"
#: the empty state (M12-04 PR 4, P-14)
EMPTY = ("Nothing to show: no novel ORF candidates in this run. Peptides that map outside "
         "the annotated genes would give candidates.")  # fmt: skip
_OPTIONAL = ("genes_tsv", "peptide_hits", "promoters_bed", "terminators_bed", "regions_tsv",
             "islands")  # fmt: skip


def _fdr(v: object) -> float | None:
    """``class_fdr`` as a finite float; NA and the literal ``unestimated`` are ``None``."""
    try:
        return locus._num(v)
    except ValueError:
        return None


def _is_called(v: object) -> bool:
    """``called``; NA (and NaN) counts as not called."""
    return (
        v is not None
        and v is not pd.NA
        and not (isinstance(v, float) and math.isnan(v))
        and bool(v)
    )


def _label(cls: object) -> str:
    """The class in plain words (no raw enum identifier in figure text)."""
    return NOVEL_CLASS_LABEL.get(str(cls), str(cls).replace("_", " "))


def _status(row: Any) -> str:
    if _is_called(row.called):
        return "called"
    flags = str(row.flags).split(";") if isinstance(row.flags, str) else []
    if _fdr(row.class_fdr) is None or "class_fdr_unestimated" in flags:
        return "unestimated"
    return "not_called"


def _fmt(v: float | None) -> str:
    return "unestimated" if v is None else style.fmt_num(v)


def _counts(orfs: pd.DataFrame) -> pd.DataFrame:
    status = [_status(r) for r in orfs.itertuples(index=False)]
    fdr = [_fdr(v) for v in orfs["class_fdr"]]
    rows = []
    for c in CLASSES:
        for s in STATUSES:
            hit = [i for i, k in enumerate(orfs["cls"]) if k == c and status[i] == s]
            if not hit:
                continue
            vals = [f for i in hit if (f := fdr[i]) is not None] if s == "called" else []
            rows.append([c, s, len(hit), max(vals) if vals else math.nan])
    data = pd.DataFrame(rows, columns=list(COUNT_COLUMNS))
    return data.astype({"n": "int64", "class_fdr_max_called": float})


def _bar_label(d: pd.DataFrame) -> str:
    called = d[d["status"] == "called"]
    if called.empty:
        return "unestimated" if set(d["status"]) == {"unestimated"} else "none called"
    fdr = called["class_fdr_max_called"].dropna()
    return f"FDR ≤ {style.fmt_num(fdr.max())}" if len(fdr) else "unestimated"


def _title(data: pd.DataFrame) -> str:
    n = {s: int(data.loc[data["status"] == s, "n"].sum()) for s in STATUSES}
    if not n["called"]:
        k = n["unestimated"]
        if not k:
            return "No novel ORF called"
        return f"No novel ORF called; {k} {'candidate' if k == 1 else 'candidates'} unestimated"
    fdr = data["class_fdr_max_called"].dropna()
    bound = f"class FDR ≤ {style.fmt_num(fdr.max())}" if len(fdr) else "class FDR unestimated"
    orfs = "ORF" if n["called"] == 1 else "ORFs"
    return f"{n['called']} novel {orfs} called ({bound}), {n['unestimated']} unestimated"


def _panel_a(orfs: pd.DataFrame) -> base.Instance:
    data = _counts(orfs)
    title = _title(data)
    fig = style.new_figure("single", 2.4, layout="constrained")
    ax = fig.add_subplot(label="counts")
    left = np.zeros(len(CLASSES))
    for s in STATUSES:
        sub = data[data["status"] == s]
        n = [int(sub.loc[sub["class"] == c, "n"].sum()) for c in CLASSES]
        ax.barh(range(len(CLASSES)), n, left=left, color=STATUS_COLOUR[s], height=0.7,
                label=s.replace("_", " "))  # fmt: skip
        left += n
    for i, c in enumerate(CLASSES):
        if left[i]:
            label = _bar_label(data[data["class"] == c])
            ax.text(left[i] + 0.02 * left.max(), i, label, va="center", fontsize=6)
    ax.set_yticks(range(len(CLASSES)), [_label(c) for c in CLASSES])
    ax.invert_yaxis()
    ax.set_xlim(0, math.ceil(1.6 * left.max()))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_xlabel("candidates")
    # below the axes: inside, the key met the bar labels (sim100k: "unestimated" twice)
    fig.legend(loc="outside lower center", ncols=len(STATUSES), frameon=False)
    fig.suptitle(title, fontsize=7)
    return base.Instance(None, fig, data, title)


def _int(v: object) -> int | None:
    f = locus._num(v)
    return None if f is None else int(f)


def _span(row: Any, view: RunView, length: int | None) -> tuple[int, int]:
    """The 0-based half-open span (decisions 6, 7); ``start_longest`` NA -> first peptide."""
    minus, end, first = str(row.strand) == "-", int(row.end), _int(row.start_longest)
    if first is None:
        peps = (
            {locus.fold_il(p) for p in row.peptides.split(";")}
            if isinstance(row.peptides, str)
            else set()
        )
        hits = view.table("peptide_hits")
        if hits is not None and len(hits):
            keys = hits["peptide"].astype(str).map(locus.fold_il)
            hits = hits[(hits["contig"] == row.contig) & keys.isin(peps)]
        if hits is None or hits.empty:
            return (end - 1, end + 2) if minus else (end - 3, end)
        return (end - 1, int(hits["end"].max())) if minus else (int(hits["start"].min()), end)
    if minus:
        lo, hi = end - 1, first + 1
    else:
        lo, hi = first, end
    if hi <= lo:  # crosses the origin
        hi += length or 0
    return (lo, hi) if hi > lo else (hi - 1, lo + 1)


def _tick(row: Any, span: tuple[int, int], length: int | None) -> int | None:
    best, first = _int(row.start_best_ostir), _int(row.start_longest)
    if best is None or best == first:
        return None
    x = best + 1 if str(row.strand) == "-" else best
    for k in (0, 1, -1):
        if span[0] <= x + k * (length or 0) <= span[1]:
            return x + k * (length or 0)
    return None


def _view(row: Any, view: RunView, cfg: GippaConfig, scale: tuple[float, float]) -> base.Instance:
    length, _ = locus._replicon(view, str(row.contig))
    span = _span(row, view, length)
    q, fdr = locus._num(row.best_peptide_q), _fdr(row.class_fdr)
    n = _int(row.n_peptides)
    called = _is_called(row.called)
    status = "called" if called else "not called"
    peps_word = "peptide" if n == 1 else "peptides"
    label = (f"{row.orf_id} ({_label(row.cls)}, {status}): {'NA' if n is None else n} "
             f"{peps_word}, best q {style.fmt_q(q)}, class FDR {_fmt(fdr)}")  # fmt: skip
    lo = Locus(str(row.orf_id), "novel_orf", str(row.contig), *span, fdr, label)
    peps = frozenset(str(row.peptides).split(";")) if isinstance(row.peptides, str) else frozenset()
    mark = OrfMark(str(row.strand), _tick(row, span, length), peps, called)
    inst = render_locus(lo, view, cfg, scale, orf=mark)
    return dataclasses.replace(inst, suffix=f"locus_{row.orf_id}")


def _order_key(row: dict[str, Any]) -> tuple[object, ...]:
    """The display order (M12-05a § 2, lead decision Q1): class FDR, best peptide q, called
    first, ``rank``, then ``rank_key``; NA (and ``unestimated``) last at every step."""
    fdr, q, rank = _fdr(row["class_fdr"]), locus._num(row.get("best_peptide_q")), _int(row["rank"])
    return (
        math.inf if fdr is None else fdr,
        math.inf if q is None else q,
        not _is_called(row["called"]),
        math.inf if rank is None else rank,
        rank_key(row),
    )


def display_order(orfs: pd.DataFrame) -> pd.DataFrame:
    """Every candidate, sorted by ``_order_key``; total since ``orf_id`` is unique."""
    rows = orfs.to_dict("records")
    keys = [_order_key({**r, "class": r["cls"]}) for r in rows]
    order = sorted(range(len(rows)), key=lambda i: keys[i])
    return orfs.iloc[order].reset_index(drop=True)


def _index_title(m: int) -> str:
    if m > TOP_N:
        return f"{INDEX_TITLE} (top {TOP_N} of {m}; all in tables/novel_orfs.tsv)"
    return f"{INDEX_TITLE} ({m} candidate{'' if m == 1 else 's'})"


def _index_data(ordered: pd.DataFrame) -> pd.DataFrame:
    m = len(ordered)
    fdr = [_fdr(v) for v in ordered["class_fdr"]]
    data = pd.DataFrame(
        {
            "order": range(1, m + 1),
            "orf_id": [str(v) for v in ordered["orf_id"]],
            "class": [str(v) for v in ordered["cls"]],
            "called": [_is_called(v) for v in ordered["called"]],
            "class_fdr": [math.nan if f is None else f for f in fdr],
            "best_peptide_q": [
                math.nan if (q := locus._num(v)) is None else q for v in ordered["best_peptide_q"]
            ],
            "n_peptides": pd.array([_int(v) for v in ordered["n_peptides"]], dtype="Int64"),
            "rank": pd.array([_int(v) for v in ordered["rank"]], dtype="Int64"),
            "rendered": [i < TOP_N for i in range(m)],
        },
        columns=list(INDEX_COLUMNS),
    )
    return data.astype({"order": "int64", "class_fdr": float, "best_peptide_q": float})


#: the four text columns of an index line, in axes fractions (M12-05a § 5)
_INDEX_X = (0.0, 0.49, 0.69, 0.80)


def _index(ordered: pd.DataFrame) -> base.Instance:
    """The index image (the first ``TOP_N`` lines) and its data (all candidates)."""
    data = _index_data(ordered)
    shown = data.head(TOP_N)
    n = len(shown)
    title = _index_title(len(data))
    head_in, body_in = 0.3, 0.2 + 0.13 * n
    fig = style.new_figure("double", head_in + body_in)
    top = body_in / (head_in + body_in)
    ax = fig.add_axes((0.02, 0.02, 0.96, top - 0.02))
    ax.set_axis_off()
    for i, r in enumerate(shown.to_dict("records")):
        fdr = None if math.isnan(r["class_fdr"]) else float(r["class_fdr"])
        cells = (f"{r['order']}. {r['orf_id']}", _label(r["class"]),
                 "called" if r["called"] else "not called", f"class FDR {_fmt(fdr)}")  # fmt: skip
        y = 1.0 - (i + 0.5) / n
        for x, text in zip(_INDEX_X, cells, strict=True):
            ax.text(x, y, text, fontsize=6, va="center", transform=ax.transAxes)
    fig.suptitle(title, fontsize=7, y=1.0 - 0.05 / (head_in + body_in), va="top")
    return base.Instance("locus", fig, data, title)


@base.plot("V13", "novel_orfs", needs=["novel_orfs"], optional=list(_OPTIONAL))
def novel_orfs(view: RunView, cfg: GippaConfig) -> list[base.Instance] | base.Skip:
    """Panel (a), then one locus view per candidate among the first ``TOP_N`` in display
    order, called or not, then the index of every candidate."""
    orfs = view.table("novel_orfs")
    if orfs is None or orfs.empty:
        return [base.empty_state("V13", "novel_orfs", EMPTY, COUNT_COLUMNS, "single")]
    orfs = orfs.rename(columns={"class": "cls"})  # ``class`` is no attribute name
    out = [_panel_a(orfs)]
    ordered = display_order(orfs)
    scale = locus.run_scale(view.table("genes_tsv"))
    out += [_view(r, view, cfg, scale) for r in ordered.head(TOP_N).itertuples(index=False)]
    return [*out, _index(ordered)]
