"""Independent oracle of the M12-04 PR 10 ``putative_start`` rule, its attack mutants, and the
scenario checks of the statistician's sign-off (``local/handoff/pr10_statistician.md`` C2, C6;
docs/specs/M12-04-outputs-revision.md "PR 10", Rule 1-5). Written from the rule text; it does
not import the rule under test, only the public chain pieces it is defined on
(``gippa.variants.align`` alignments, ``gippa.variants.flags.start_codons``).

:func:`oracle_flagged` returns the ``variant_id`` of every row the rule flags. ``mutant`` swaps in
one of the attack's wrong rules:

- ``any_first`` (1): some ALT alignment at the codon has the M at peptide position 0;
- ``genomic_order`` (2): positions counted from the lowest coordinate on the minus strand;
- ``hardcoded_starts`` (3): ``{GTG, TTG}`` in place of ``annotation.start_codons`` minus ATG;
- ``no_modulo`` (4): alignment codons compared with ``codon_start`` without ``mod L``;
- ``no_kr_guard`` (5): condition 5 dropped.

Mutants 6 (flag without the cap) and 7 (``putative_start`` exempt in ``id_mixups``) act on the
chain's output, not on the flagged set: see :func:`uncapped_confidence` and
:func:`id_mixups_exempting`.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from typing import Any

import pandas as pd

from gippa.mapping.sixframe import translate_text

FLAG = "putative_start"
RULE_MUTANTS = ("any_first", "genomic_order", "hardcoded_starts", "no_modulo", "no_kr_guard")
CHAIN_MUTANTS = ("no_cap", "exempt_in_id_mixups")
MUTANTS = (*RULE_MUTANTS, *CHAIN_MUTANTS)
_COMP = str.maketrans("ACGTN", "TGCAN")


def _rc(nt: str) -> str:
    return nt.translate(_COMP)[::-1]


def _na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and v != v)


def row_flags(v: object) -> tuple[str, ...]:
    if _na(v):
        return ()
    if isinstance(v, str):
        return tuple(x for x in v.split(";") if x)
    return tuple(str(x) for x in v)  # type: ignore[union-attr]


def _codon_nt(seq: str, length: int, circular: bool, low: int, strand: str) -> str | None:
    if circular:
        nt = "".join(seq[(low + j) % length] for j in range(3))
    elif low >= 0 and low + 3 <= length:
        nt = seq[low : low + 3]
    else:
        return None
    return _rc(nt) if strand == "-" else nt


def preceding_residue(row: Mapping[str, Any], genome: Mapping[str, str], gm: Any,
                      table: int) -> str | None:  # fmt: skip
    """Reference residue of the codon before the site in transcription order; None off a linear
    contig."""
    contig, strand, c = str(row["contig"]), str(row["strand"]), int(row["codon_start"])
    rep = gm.replicons[contig]
    low = c - 3 if strand == "+" else c + 3
    nt = _codon_nt(genome[contig], rep.length, rep.topology == "circular", low, strand)
    return None if nt is None else translate_text(nt, table)


def _positions(row: Mapping[str, Any], aligns: Iterable[Any], gm: Any, mutant: str | None
               ) -> list[int]:  # fmt: skip
    contig, strand, c = str(row["contig"]), str(row["strand"]), int(row["codon_start"])
    rep = gm.replicons[contig]
    circular = rep.topology == "circular"
    out: list[int] = []
    for a in aligns:
        if a.contig != contig or a.strand != strand:
            continue
        for op in a.ops:
            if op.kind != "substitution" or not op.codons:
                continue
            oc = op.codons[0]
            if circular and mutant != "no_modulo":
                oc %= rep.length
            if oc != c:
                continue
            q = op.pep[0]
            if mutant == "genomic_order" and strand == "-":
                q = len(a.peptide) - 1 - q
            out.append(q)
    return out


def oracle_flagged(
    table: pd.DataFrame,
    support: pd.DataFrame,
    aligns: Sequence[Any],
    genome: Mapping[str, str],
    gm: Any,
    cfg: Any,
    mutant: str | None = None,
) -> set[str]:
    """``variant_id`` of the rows the rule (or ``mutant``) flags ``putative_start``.

    ``support``: the ``(variant_id, peptide, ...)`` frame of ``run_variants_qc`` (ALT peptides per
    row); ``aligns``: ``align_hits`` of every peptide hit (as ``run_variants_qc`` builds them)."""
    from gippa.variants.flags import start_codons

    assert mutant is None or mutant in RULE_MUTANTS, mutant
    table_no = cfg.annotation.translation_table
    if mutant == "hardcoded_starts":
        start_set = {"GTG", "TTG"}
    else:
        start_set = {c.upper() for c in cfg.annotation.start_codons} - {"ATG"}
    annotated = start_codons(gm)
    by_pep: dict[str, list[Any]] = {}
    for a in aligns:
        by_pep.setdefault(a.peptide, []).append(a)
    peps: dict[str, set[str]] = {}
    for vid, p in zip(support["variant_id"].astype(str), support["peptide"].astype(str),
                      strict=True):  # fmt: skip
        peps.setdefault(vid, set()).add(p)
    out: set[str] = set()
    for r in table.to_dict("records"):
        if r["alt_type"] != "substitution" or str(r["alt_aa"]) != "M":
            continue
        if str(r["ref_codon"]).upper() not in start_set:
            continue
        if (str(r["contig"]), str(r["strand"]), int(r["codon_start"])) in annotated:
            continue
        vid = str(r["variant_id"])
        qs = _positions(r, [a for p in sorted(peps.get(vid, ())) for a in by_pep.get(p, ())],
                        gm, mutant)  # fmt: skip
        if not qs:
            continue
        if mutant == "any_first":
            if 0 not in qs:
                continue
        elif any(q != 0 for q in qs):
            continue
        if mutant != "no_kr_guard" and not _na(r["gene"]):
            prev = preceding_residue(r, genome, gm, table_no)
            if prev is not None and prev in ("K", "R"):
                continue
        out.add(vid)
    return out


def chain_alignments(peptide_hits: pd.DataFrame, genome: Mapping[str, str], gm: Any,
                     cfg: Any, hits_tolerant: pd.DataFrame | None = None) -> list[Any]:  # fmt: skip
    """``align_hits`` exactly as ``run_variants_qc`` calls it (stop sites, initiator starts)."""
    from gippa.mapping.initiator import annotated_starts, start_keys
    from gippa.mapping.readthrough import annotated_stops, site_keys
    from gippa.variants.align import align_hits

    tt = cfg.annotation.translation_table
    stops = site_keys(annotated_stops(gm, genome, recoded=True, table=tt)[0])
    inits = start_keys(annotated_starts(gm, genome, cfg.annotation.start_codons, table=tt)[0])
    return align_hits(peptide_hits, hits_tolerant, genome, gm.replicons, stop_sites=stops,
                      table=tt, starts=inits)  # fmt: skip


def uncapped_confidence(table: pd.DataFrame, cfg: Any) -> pd.Series:
    """Mutant 6: confidence as M7-02 sets it (``n_alt >= high_min_peptides``), no cap."""
    k = cfg.variants.high_min_peptides
    return pd.Series(["high" if int(n) >= k else "candidate" for n in table["n_alt"]],
                     index=table.index, dtype=object)  # fmt: skip


def id_mixups_exempting(table: pd.DataFrame, support: pd.DataFrame,
                        extra: Iterable[str]) -> tuple[int, int]:  # fmt: skip
    """``id_mixups`` with ``extra`` added to ``EXEMPT_FLAGS`` (mutant 7)."""
    import gippa.variants.qc as qc

    saved = qc.EXEMPT_FLAGS
    try:
        qc.EXEMPT_FLAGS = frozenset(saved | set(extra))
        return qc.id_mixups(table, support)
    finally:
        qc.EXEMPT_FLAGS = saved


# ---------------------------------------------------------------- scenario checks

A_KINDS = ("A1", "A2", "A3")
B_KINDS = ("B1", "B2", "B3", "B_origin")
TRUTH_COLUMNS = ("kind", "contig", "strand", "codon_start", "ref_codon", "alt_aa", "in_gene",
                 "prev_kr", "crosses_origin", "all_first", "peptides")  # fmt: skip


def _key(contig: object, strand: object, codon: object, alt: object) -> tuple[str, str, int, str]:
    return (str(contig), str(strand), int(codon), str(alt))  # type: ignore[call-overload]


def match_truth(table: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    """``truth`` with ``variant_id`` (NA when no ``variants.tsv`` row has its key)."""
    idx = {
        _key(c, s, k, a): v
        for v, c, s, k, a in zip(table["variant_id"], table["contig"], table["strand"],
                                 table["codon_start"], table["alt_aa"], strict=True)
    }  # fmt: skip
    out = truth.copy()
    out["variant_id"] = [
        idx.get(_key(c, s, k, a), pd.NA)
        for c, s, k, a in zip(truth["contig"], truth["strand"], truth["codon_start"],
                              truth["alt_aa"], strict=True)
    ]  # fmt: skip
    return out


def checks(
    table: pd.DataFrame,
    truth: pd.DataFrame,
    flagged: set[str],
    confidence: Mapping[str, str],
    id_mixups: tuple[int, int],
    unguarded: set[str] | None = None,
) -> dict[str, Any]:
    """The sign-off's per-replicate thresholds as counts and verdicts.

    ``flagged``/``confidence``: the (possibly mutated) ``putative_start`` set and confidence per
    ``variant_id``. ``unguarded``: the flagged set of the rule without condition 5 (recorded)."""
    t = match_truth(table, truth)
    has = t["variant_id"].notna()
    a = t[has & t["kind"].isin(A_KINDS)]
    b = t[has & t["kind"].isin(B_KINDS)]
    b2 = t[has & t["kind"].eq("B2")]
    ctrl_km = t[has & t["kind"].eq("control_kM")]
    ctrl_start_rows = int((has & t["kind"].eq("control_start")).sum())
    a_flag = sum(str(v) in flagged for v in a["variant_id"])
    a_high = sum(confidence[str(v)] == "high" for v in a["variant_id"])
    b_flag = sum(str(v) in flagged for v in b["variant_id"])
    km_flag = sum(str(v) in flagged for v in ctrl_km["variant_id"])
    rec = {
        "n_a_planted": int(t["kind"].isin(A_KINDS).sum()),
        "n_a_rows": len(a), "a_flagged": a_flag, "a_high": a_high,
        "n_b_rows": len(b), "b_flagged": b_flag,
        "n_b2_rows": len(b2),
        "b2_unguarded": (None if unguarded is None
                         else sum(str(v) in unguarded for v in b2["variant_id"])),
        "n_control_km_rows": len(ctrl_km), "control_km_flagged": km_flag,
        "control_start_rows": ctrl_start_rows,
        "id_mixups": list(id_mixups),
    }  # fmt: skip
    rec["ok"] = {
        "a_share_1": len(a) > 0 and a_flag == len(a),
        "a_no_high": a_high == 0,
        "b_none": b_flag == 0,
        "controls": km_flag == 0 and ctrl_start_rows == 0,
        "id_mixups": tuple(id_mixups) == (1, 1),
    }
    return rec


def wilson_upper(k: int, n: int, z: float = 1.959963984540054) -> float:
    """Upper end of the Wilson score 95 % interval of ``k / n``."""
    if n == 0:
        return 1.0
    p = k / n
    den = 1 + z * z / n
    centre = p + z * z / (2 * n)
    half = z * ((p * (1 - p) / n + z * z / (4 * n * n)) ** 0.5)
    return (centre + half) / den
