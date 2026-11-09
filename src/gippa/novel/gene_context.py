"""Host and flanking genes of novel ORFs (M12-05b; GIPPA_SPEC.md §9.3 ``tables/novel_orfs.tsv``).

Pure. For each ``orf_table`` row: the label and strand relation of its host gene
(``overlapping_gene``), the other genes overlapping its span, and the nearest gene on each side
in the ORF's transcription direction with the number of bases strictly between (rules 2-7 of
docs/specs/M12-05b-novel-orf-genes.md). Only the gene model's genes count (no tRNA, rRNA, CRISPR
or other novel ORFs). Circular replicons wrap at most once around; linear ones stop at the end;
a walk never leaves its replicon. Ties go to the smallest locus tag.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import pandas as pd

from gippa.model.coords import Interval, Replicon, split_at_origin
from gippa.model.genes import Gene, GeneModel, gene_label
from gippa.novel.outputs import GENE_CONTEXT_COLUMNS, _opt, _span, is_na

_DISTANCES = ("upstream_distance_nt", "downstream_distance_nt")


def _overlaps(a: Interval, b: Interval, rep: Replicon) -> bool:
    return any(
        p.start < q.end and q.start < p.end
        for p in split_at_origin(a, rep)
        for q in split_at_origin(b, rep)
    )


def _flank(
    span: Interval, cands: list[Gene], rep: Replicon, *, left: bool
) -> tuple[int, str, Gene] | None:
    """The nearest candidate on the lower-coordinate (``left``) or higher side of ``span``."""
    n, circular = rep.length, rep.topology == "circular"
    best: tuple[int, str, Gene] | None = None
    for g in cands:
        iv = g.interval
        d = span.start - iv.end if left else iv.start - span.end
        if circular:
            d %= n
        elif d < 0:
            continue
        key = (d, g.locus_tag, g)
        if best is None or key[:2] < best[:2]:
            best = key
    return best


def _row_context(r: Mapping[str, Any], gm: GeneModel, by_contig: Mapping[str, list[Gene]]) -> (
    dict[str, object]
):  # fmt: skip
    rep = gm.replicons[str(r["contig"])]
    span = _span(r, _opt(r["start_longest"]), int(r["length_aa_longest"]))
    genes = by_contig.get(rep.id, [])
    over = [g for g in genes if _overlaps(span, g.interval, rep)]
    host_tag = r["overlapping_gene"]
    host = None if is_na(host_tag) else gm.by_locus(str(host_tag))
    others = [g.locus_tag for g in over if host is None or g.locus_tag != host.locus_tag]
    out: dict[str, object] = {
        "overlapping_gene_label": None if host is None else gene_label(host),
        "overlap_relation": "none" if host is None
        else "same_strand" if host.interval.strand == r["strand"] else "antisense",
        "other_overlapping_genes": ";".join(others) or None,
    }  # fmt: skip
    cands = [g for g in genes if g not in over]
    plus = r["strand"] == "+"
    for side, left in (("upstream", plus), ("downstream", not plus)):
        best = _flank(span, cands, rep, left=left)
        out[f"{side}_gene"] = None if best is None else best[1]
        out[f"{side}_gene_label"] = None if best is None else gene_label(best[2])
        out[f"{side}_distance_nt"] = None if best is None else best[0]
    return out


def gene_context(table: pd.DataFrame, gm: GeneModel) -> pd.DataFrame:
    """The nine gene-context columns of ``table`` (an ``orf_table`` frame), indexed like it.

    The span is ``start_longest`` (or the first peptide's codon) through the stop codon, as
    ``orf_id`` names it; upstream and downstream follow the ORF's strand."""
    by_contig: dict[str, list[Gene]] = {}
    for g in gm.genes:
        by_contig.setdefault(g.interval.contig, []).append(g)
    rows = [_row_context(r, gm, by_contig) for r in table.to_dict("records")]
    data = {
        c: pd.Series([r[c] for r in rows], dtype="Int64" if c in _DISTANCES else "object",
                     index=table.index)
        for c in GENE_CONTEXT_COLUMNS
    }  # fmt: skip
    return pd.DataFrame(data, columns=list(GENE_CONTEXT_COLUMNS), index=table.index)
