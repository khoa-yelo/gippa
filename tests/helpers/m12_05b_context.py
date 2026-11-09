"""M12-05b test helpers: the synthetic gene model of the worked examples, novel_orfs rows and a
brute-force oracle for the host and flanking genes (docs/specs/M12-05b-novel-orf-genes.md).

Independent of ``gippa.novel.gene_context``. The oracle works on position sets (each base of a
circular replicon taken mod its length), so origin-crossing spans and genes need no special case:

* overlap: the ORF span and a gene share at least one base (rule 3, edge case "one base counts");
* a flank candidate is a gene on the ORF's contig that shares no base with the span (rule 4);
* on ``+`` the upstream gap is ``span.start - gene.end`` and the downstream gap
  ``gene.start - span.end``; on a circular replicon both are taken mod the length (the walk wraps
  at most once); on a linear one only genes wholly on that side count; ``-`` swaps the sides;
* ties: the smallest locus tag (rule 4).

Rows follow the ``orf_table`` representation (gippa.novel.outputs docstring): ``start_longest``
and ``end`` are 0-based positions in transcription direction, so a span ``[lo, hi)`` is
``start_longest = lo, end = (hi - 1) mod L`` on ``+`` and ``start_longest = (hi - 1) mod L,
end = lo`` on ``-``, with ``length_aa_longest = (hi - lo - 3) / 3`` (stop codon included).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, GeneModel

#: the new columns, in order (rule 2's table). The spec text says "eight" and "34 columns", but its
#: rule 2 table lists nine (26 + 9 = 35); the table is taken as binding (reported to the lead).
NEW_COLUMNS: tuple[str, ...] = (
    "overlapping_gene_label",
    "overlap_relation",
    "other_overlapping_genes",
    "upstream_gene",
    "upstream_gene_label",
    "upstream_distance_nt",
    "downstream_gene",
    "downstream_gene_label",
    "downstream_distance_nt",
)
DISTANCES = ("upstream_distance_nt", "downstream_distance_nt")
TAG_COLUMNS = ("overlapping_gene", "upstream_gene", "downstream_gene")
#: the 1.1.0rc1 novel_orfs.tsv header (26 columns; tests/data/golden/sim50k_m7)
RC1_COLUMNS: tuple[str, ...] = (
    "orf_id", "class", "contig", "strand", "frame", "start_longest", "start_best_ostir", "end",
    "length_aa_longest", "length_aa_best", "start_codon_best", "n_peptides", "peptides",
    "best_peptide_q", "tir_ostir", "promoter_score", "promoter_supported",
    "terminator_downstream", "terminator_conf", "terminator_distance_nt", "rbs_supported",
    "class_fdr", "called", "rank", "overlapping_gene", "flags",
)  # fmt: skip
_AT = RC1_COLUMNS.index("overlapping_gene") + 1
#: the 1.1 header after M12-05b: the new columns right after ``overlapping_gene``
COLUMNS: tuple[str, ...] = (*RC1_COLUMNS[:_AT], *NEW_COLUMNS, *RC1_COLUMNS[_AT:])


def gene(tag: str, contig: str, start: int, end: int, strand: str, name: str | None) -> Gene:
    return Gene(tag, "cds", Interval(contig, start, end, strand), "M", None, None,  # type: ignore[arg-type]
                gene=name, product=f"product of {tag}")  # fmt: skip


#: the four genes of the worked-examples table (T_0002's name is the literal "NA")
BASE_GENES: tuple[Gene, ...] = (
    gene("T_0001", "c1", 100, 400, "+", "dnaA"),
    gene("T_0002", "c1", 500, 1100, "-", "NA"),
    gene("T_0003", "c1", 1300, 2500, "+", "lacZ"),
    gene("T_0004", "c1", 9500, 9900, "+", None),
)
P_0001 = gene("P_0001", "p1", 200, 1200, "+", "repA")
E10_GENES = (gene("T_0006", "c1", 4800, 5400, "-", None),
             gene("T_0005", "c1", 5000, 5400, "+", None))  # fmt: skip


def model(*, c1: str = "circular", p1: str = "circular", extra: Sequence[Gene] = ()) -> GeneModel:
    """``c1`` (10 000), ``c2`` (5 000, no genes) and the plasmid ``p1`` (3 000, P_0001)."""
    reps = {
        "c1": Replicon("c1", 10_000, c1, "chromosome"),  # type: ignore[arg-type]
        "c2": Replicon("c2", 5_000, "circular", "chromosome"),
        "p1": Replicon("p1", 3_000, p1, "plasmid"),  # type: ignore[arg-type]
    }
    order = {c: i for i, c in enumerate(reps)}
    genes = sorted((*BASE_GENES, *extra, P_0001), key=lambda g: (
        order[g.interval.contig], g.interval.start, g.interval.end, g.interval.strand,
        g.locus_tag))  # fmt: skip
    return GeneModel(reps, tuple(genes), ())


def orf_row(
    contig: str, strand: str, lo: int, hi: int, cls: str, host: str | None, length: int
) -> dict[str, Any]:
    """One novel_orfs row (``orf_table`` representation) with the span ``[lo, hi)``."""
    n = hi - lo
    assert n % 3 == 0 and n >= 6, (lo, hi)
    plus = strand == "+"
    start, end = (lo, (hi - 1) % length) if plus else ((hi - 1) % length, lo)
    n_aa = n // 3 - 1
    base: dict[str, Any] = {
        "orf_id": f"NOVEL_{contig}_{lo + 1}_{'p' if plus else 'm'}", "class": cls,
        "contig": contig, "strand": strand, "frame": (lo % 3 + 1) if plus else -1,
        "start_longest": start, "start_best_ostir": None, "end": end,
        "length_aa_longest": n_aa, "length_aa_best": n_aa, "start_codon_best": "ATG",
        "n_peptides": 1, "peptides": "PEPTIDEK", "best_peptide_q": 0.001, "class_fdr": 0.01,
        "called": True, "rank": 1, "overlapping_gene": host, "flags": None,
    }  # fmt: skip
    return base


def frame(rows: Sequence[Mapping[str, Any]], index: Sequence[int] | None = None) -> pd.DataFrame:
    """A typed novel_orfs frame (``novel_frame``) over the module's current column list; columns
    a row leaves out are NA."""
    from gippa.novel.outputs import NOVEL_ORFS_COLUMNS, novel_frame

    t = novel_frame([{c: r.get(c) for c in NOVEL_ORFS_COLUMNS} for r in rows])
    if index is not None:
        t.index = pd.Index(list(index))
    return t


def is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and v != v)


def cell(v: object) -> object:
    """A frame cell as ``None`` (NA), ``int`` or ``str``."""
    if is_na(v):
        return None
    if isinstance(v, bool):
        return v
    if isinstance(v, int) or (hasattr(v, "dtype") and v.dtype.kind in "iu"):
        return int(v)  # type: ignore[call-overload]
    return v


def label(g: Gene) -> str:
    """Rule 5, independently: the name unless None, blank or the literal NA (any case)."""
    name = g.gene
    if name is None or not name.strip() or name.strip().upper() == "NA":
        return g.locus_tag
    return name


# ---------------------------------------------------------------- oracle


def _bases(start: int, end: int, rep: Replicon) -> set[int]:
    if rep.topology == "circular":
        return {p % rep.length for p in range(start, end)}
    return set(range(start, end))


def span_of(row: Mapping[str, Any], length: int) -> tuple[int, int]:
    """``[lo, hi)`` of a row built by ``orf_row`` (origin-crossing: ``hi > length``)."""
    n = 3 * int(row["length_aa_longest"]) + 3
    lo = int(row["start_longest"]) if row["strand"] == "+" else int(row["end"])
    return lo, lo + n


def oracle(row: Mapping[str, Any], gm: GeneModel) -> dict[str, object]:
    """The expected nine cells of one row (``None`` = NA)."""
    rep = gm.replicons[str(row["contig"])]
    lo, hi = span_of(row, rep.length)
    span = _bases(lo, hi, rep)
    genes = [g for g in gm.genes if g.interval.contig == rep.id]
    over = [g for g in genes if _bases(g.interval.start, g.interval.end, rep) & span]
    host = row["overlapping_gene"]
    host_g = None if is_na(host) else gm.by_locus(str(host))
    others = [g.locus_tag for g in over if g.locus_tag != host]
    out: dict[str, object] = {
        "overlapping_gene_label": None if host_g is None else label(host_g),
        "overlap_relation": "none" if host_g is None
        else "same_strand" if host_g.interval.strand == row["strand"] else "antisense",
        "other_overlapping_genes": ";".join(others) or None,
    }  # fmt: skip
    cands = [g for g in genes if g not in over]
    circ = rep.topology == "circular"

    def left(g: Gene) -> int | None:  # gap to a gene on the lower-coordinate side
        if circ:
            return (lo - g.interval.end) % rep.length
        return lo - g.interval.end if g.interval.end <= lo else None

    def right(g: Gene) -> int | None:
        if circ:
            return (g.interval.start - hi) % rep.length
        return g.interval.start - hi if g.interval.start >= hi else None

    sides = (left, right) if row["strand"] == "+" else (right, left)
    for side, gap in zip(("upstream", "downstream"), sides, strict=True):
        best = min(((d, g.locus_tag, g) for g in cands if (d := gap(g)) is not None),
                   default=None)  # fmt: skip
        out[f"{side}_gene"] = None if best is None else best[1]
        out[f"{side}_gene_label"] = None if best is None else label(best[2])
        out[f"{side}_distance_nt"] = None if best is None else best[0]
    return out
