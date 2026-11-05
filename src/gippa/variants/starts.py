"""Putative unannotated starts (docs/specs/M12-04-outputs-revision.md "PR 10", Rule 1-5;
docs/methods/variants.md "Putative starts").

Bacteria read the first codon of a protein as methionine whatever it is, so a peptide that begins
at an unannotated GTG/TTG start maps with one V/L->M difference at its first residue. Such an ALT
row gets the flag ``putative_start`` and its confidence is capped at ``candidate``; its counts,
call and class error rate do not change. Pure functions.
"""

from __future__ import annotations

import math
from collections.abc import Collection, Iterable, Mapping

import pandas as pd

from gippa.mapping.sixframe import translate_text
from gippa.model.coords import Replicon
from gippa.variants.align import HitAlignment

__all__ = ["FLAG", "is_putative_start", "putative_start_codons"]

FLAG = "putative_start"
_StartKey = tuple[str, str, int]
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
_TRYPTIC = frozenset({"K", "R"})


def putative_start_codons(start_codons: Iterable[str]) -> frozenset[str]:
    """Rule 2: the run's ``annotation.start_codons`` (upper case) minus ATG. Pure."""
    return frozenset(c.upper() for c in start_codons) - {"ATG"}


def _preceding_residue(
    seq: str, rep: Replicon, codon_start: int, strand: str, table: int
) -> str | None:
    """Reference residue of the codon before ``codon_start`` in transcription order; None when
    it runs off a linear contig."""
    low = codon_start - 3 if strand == "+" else codon_start + 3
    n = rep.length
    if rep.topology == "circular":
        nt = "".join(seq[(low + j) % n] for j in range(3))
    elif low >= 0 and low + 3 <= n:
        nt = seq[low : low + 3]
    else:
        return None
    nt = nt.upper()
    if strand == "-":
        nt = nt.translate(_COMP)[::-1]
    return translate_text(nt, table)


def _all_first(
    contig: str, strand: str, codon_start: int, rep: Replicon, aligns: Iterable[HitAlignment]
) -> bool:
    """Rule 4: at least one substitution op of ``aligns`` sits at the codon (``mod L`` on a
    circular replicon), and every such op is at peptide (transcription) position 0."""
    circular = rep.topology == "circular"
    seen = False
    for a in aligns:
        if a.contig != contig or a.strand != strand:
            continue
        for op in a.ops:
            if op.kind != "substitution" or not op.codons:
                continue
            oc = op.codons[0] % rep.length if circular else op.codons[0]
            if oc != codon_start:
                continue
            if op.pep[0] != 0:
                return False
            seen = True
    return seen


def _is_na(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def is_putative_start(
    site: Mapping[str, object],
    alt_aligns: Iterable[HitAlignment],
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    start_set: Collection[str],
    annotated: Collection[_StartKey],
    *,
    table: int = 11,
) -> bool:
    """Rules 1-5 for one ``sites`` row; ``alt_aligns`` are the alignments of its ALT peptides,
    ``start_set`` is :func:`putative_start_codons`, ``annotated`` is
    ``gippa.variants.flags.start_codons``. Pure."""
    if site["alt_type"] != "substitution" or str(site["alt_aa"]) != "M":
        return False
    if str(site["ref_codon"]).upper() not in start_set:
        return False
    contig, strand = str(site["contig"]), str(site["strand"])
    c = int(site["codon_start"])  # type: ignore[call-overload]
    if (contig, strand, c) in annotated:
        return False
    rep = replicons[contig]
    if not _all_first(contig, strand, c, rep, alt_aligns):
        return False
    if _is_na(site.get("gene")):
        return True
    return _preceding_residue(genome[contig], rep, c, strand, table) not in _TRYPTIC
