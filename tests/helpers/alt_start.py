"""Scenario ``alt_start_initiation``: peptides that begin at unannotated alternative starts
(M12-04 PR 10; statistician sign-off ``local/handoff/pr10_statistician.md`` C6).

A self-contained world for the ``variants`` chain (``gippa.steps.variants.run_variants_qc``): a
circular chromosome ``chr1`` of ``size`` bp (``100k`` or ``1M``) and a circular plasmid ``p1``
(3 kb), random background, with designed *cassettes* written into it. A cassette is ``W``
transcription codons on one strand (``+``: codon ``i`` at ``base + 3 i``; ``-``: codon ``i`` at
``base - 3 i - 3``; coordinates taken mod L), its site at transcription codon ``S``. Residues
are drawn from a filler alphabet without K, R, M, I, L, V or C, so the only K/R before a site and
the only M in a peptide are the planted ones. Per replicate:

- ``A1`` (20): an in-frame alternative start 5-30 codons upstream of an annotated ATG start (an
  N-terminal extension; the site has no gene). One is the minus-strand site at ``codon_start 3``
  of ``chr1``, whose peptides run through the origin (unwrapped placements,
  ``crosses_origin``);
- ``A2`` (20): initiations in unannotated intergenic ORFs; 5 (25 %) have K/R before the site;
- ``A3`` (10): internal in-gene alternative starts, the residue before not K/R;
- every A site: 2-4 peptides, each starting at the site with the initiator M (position 0), the
  second residue not in GACSTPV; half the A sites on the minus strand;
- ``B1`` (30): true V/L->M SNVs (GTG/TTG/CTG read as M) inside genes, the residue before not K/R,
  read by tryptic peptides (each after a K/R) with the M at position >= 1;
- ``B2`` (30): such SNVs directly after K/R: a peptide starting at the M (cleaved after the K/R)
  and, for most sites, a missed-cleavage peptide starting earlier; the first 4 (all GTG/TTG) have
  only peptides starting at the M (``all_first``);
- ``B3`` (6): out-of-gene alternative-start codons read as M by one peptide at position 0 and one
  at a later position (not ``all_first``);
- ``B_origin`` (1): the B3 pattern on the minus strand of ``p1`` at ``codon_start 30``, its later
  peptide crossing the origin (alignment codon ``L + 30``);
- ``control_start`` (5): annotated non-ATG (GTG/TTG) starts of genes whose protein begins with M,
  read by N-terminal peptides as M (``n_diff 0``: M11-20 credit, no ALT row);
- ``control_kM`` (5): out-of-gene AAG read as M at peptide position 0;
- ``same_db`` (1): an out-of-gene GTG A-like site whose peptides include one ``database`` peptide
  searched against this strain (``search_db = same``): C1's planted ``id_mixups`` case.

Peptide sources: ``related`` arm ``database`` / ``other``; ``same`` arm ``denovo`` / ``same``
(a ``database``/``other`` peptide loses its differences under ``same``, §7.5.1 row 2). A-site
reference codons cycle over ``start_codons`` minus ATG (so the table-4 arm with ``TTA`` plants
TTA sites). Test tooling; deterministic in ``seed``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np
import pandas as pd

from gippa.config import GippaConfig
from gippa.model import Replicon
from gippa.model.coords import Interval
from gippa.model.genes import Gene, GeneModel
from gippasim.codons import code, revcomp

__all__ = ["SIZES", "TRUTH_COLUMNS", "AltStartScenario", "alt_start_initiation"]

SIZES = {"100k": 100_000, "1M": 1_000_000}
PLASMID_BP = 3_000
TRUTH_COLUMNS = (
    "kind", "contig", "strand", "codon_start", "ref_codon", "alt_aa", "in_gene", "prev_kr",
    "crosses_origin", "all_first", "peptides",
)  # fmt: skip
W, S = 100, 40  # cassette codons, site codon
FILL = "AGSTEDQNPFYHW"
SECOND = "EDQNFYHW"  # not in GACSTPV: the initiator Met is kept
CODON = {
    "A": "GCT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT", "K": "AAA",
    "M": "ATG", "N": "AAT", "P": "CCG", "Q": "CAG", "R": "CGT", "S": "TCT", "T": "ACC",
    "W": "TGG", "Y": "TAT",
}  # fmt: skip
COUNTS = {"A1": 20, "A2": 20, "A3": 10, "B1": 30, "B2": 30, "B3": 6, "control_start": 5,
          "control_kM": 5, "same_db": 1}  # fmt: skip
B_CODONS = ("GTG", "TTG", "CTG")


@dataclass(frozen=True)
class AltStartScenario:
    """One replicate: the chain's inputs and the planted truth (one row per site)."""

    genome: dict[str, str]
    gm: GeneModel
    cfg: GippaConfig
    peptide_hits: pd.DataFrame
    peptides: pd.DataFrame
    truth: pd.DataFrame


@dataclass
class _Cassette:
    kind: str
    contig: str
    strand: str
    base: int  # + : genome start of codon 0; - : genome end of codon 0 (unwrapped)
    codons: list[str]
    ref_codon: str
    #: (q, n, source, search_db): the peptide covers codons [S - q, S - q + n)
    peptides: list[tuple[int, int, str, str]] = field(default_factory=list)
    genes: list[tuple[int, int]] = field(default_factory=list)  # (first, end) codons
    in_gene: bool = False
    crosses_origin: bool = False
    initiator: bool = False  # control_start: the M is the credited initiator (n_diff 0)


def _low(c: _Cassette, i: int) -> int:
    return c.base + 3 * i if c.strand == "+" else c.base - 3 * i - 3


def _span(c: _Cassette, a: int, b: int, length: int) -> tuple[int, int]:
    """Genome ``[start, end)`` of codons ``[a, b)``, ``start`` in ``[0, L)`` (unwrapped end)."""
    lo = min(_low(c, a), _low(c, b - 1))
    hi = max(_low(c, a), _low(c, b - 1)) + 3
    shift = (lo // length) * length
    return lo - shift, hi - shift


class _Builder:
    def __init__(self, rng: np.random.Generator, table: int, starts: Sequence[str],
                 relation: str) -> None:  # fmt: skip
        self.rng, self.table, self.relation = rng, table, relation
        self.code = code(table)
        self.alt_codons = sorted({c.upper() for c in starts} - {"ATG"})
        self.src = ("database", "other") if relation == "related" else ("denovo", "same")
        self.n_a = 0

    def filler(self, n: int) -> list[str]:
        return [CODON[str(a)] for a in self.rng.choice(list(FILL), size=n)]

    def cassette(self, kind: str, contig: str, strand: str, base: int, ref: str,
                 prev: str) -> _Cassette:  # fmt: skip
        cods = self.filler(W)
        cods[S] = ref
        cods[S - 1] = CODON[prev]
        cods[S + 1] = CODON[str(self.rng.choice(list(SECOND)))]
        return _Cassette(kind, contig, strand, base, cods, ref)

    def first_peptides(self, k: int, src: tuple[str, str] | None = None
                       ) -> list[tuple[int, int, str, str]]:  # fmt: skip
        """``k`` peptides starting at the site, distinct lengths 7-14."""
        lengths = self.rng.choice(np.arange(7, 15), size=k, replace=False)
        s = src or self.src
        return [(0, int(n), *s) for n in sorted(lengths)]

    def gene(self, c: _Cassette, g0: int, g1: int) -> None:
        c.codons[g0] = "ATG"
        c.codons[g1 - 1] = "TAA"
        c.genes.append((g0, g1))

    def a_ref(self) -> str:
        ref = self.alt_codons[self.n_a % len(self.alt_codons)]
        self.n_a += 1
        return ref

    def prev_not_kr(self) -> str:
        return str(self.rng.choice(list(FILL)))


def _strand(i: int) -> str:
    return "-" if i % 2 else "+"


def _plan(b: _Builder, length: int) -> list[_Cassette]:
    rng = b.rng
    kinds = [k for k, n in COUNTS.items() for _ in range(n)]
    span = 3 * W
    lo, hi = 3 * span, length - 3 * span  # keep the origin of chr1 free
    step = (hi - lo) // len(kinds)
    if step < span + 30:
        raise ValueError(f"chromosome of {length} bp too short for {len(kinds)} cassettes")
    out: list[_Cassette] = []
    # A1 at the origin of chr1: minus strand, codon_start 3, peptides through the origin
    a_i = 0
    origin = b.cassette("A1", "chr1", "-", 3 * S + 6, b.a_ref(), b.prev_not_kr())
    origin.crosses_origin = True
    e = int(rng.integers(5, 31))
    b.gene(origin, S + e, S + e + 25)
    origin.peptides = b.first_peptides(int(rng.integers(2, 5)))
    out.append(origin)
    a_i += 1
    b2_i = 0
    for j, kind in enumerate(kinds):
        pos = lo + j * step + int(rng.integers(0, step - span))
        if kind == "A1" and a_i == COUNTS["A1"]:
            continue  # the origin site is one of the 20
        if kind in ("A1", "A2", "A3"):
            strand = _strand(a_i)
            a_i += 1
        else:
            strand = _strand(j)
        base = pos if strand == "+" else pos + span
        if kind in ("A1", "A2", "A3"):
            n_a2_kr = sum(
                1 for c in out if c.kind == "A2" and c.codons[S - 1] in (CODON["K"], CODON["R"])
            )
            prev = b.prev_not_kr()
            if kind == "A2" and n_a2_kr < 5:
                prev = "K" if n_a2_kr % 2 == 0 else "R"
            c = b.cassette(kind, "chr1", strand, base, b.a_ref(), prev)
            if kind == "A1":
                e = int(rng.integers(5, 31))
                b.gene(c, S + e, S + e + 25)
            elif kind == "A3":
                b.gene(c, 10, 90)
                c.in_gene = True
            c.peptides = b.first_peptides(int(rng.integers(2, 5)))
        elif kind == "B1":
            c = b.cassette(kind, "chr1", strand, base, str(rng.choice(B_CODONS)),
                           b.prev_not_kr())  # fmt: skip
            b.gene(c, 10, 90)
            c.in_gene = True
            qs = sorted(int(q) for q in rng.choice(np.arange(1, 9), size=int(rng.integers(2, 4)),
                                                   replace=False))  # fmt: skip
            for q in qs:
                c.codons[S - q - 1] = CODON[str(rng.choice(["K", "R"]))]
                c.peptides.append((q, q + int(rng.integers(5, 9)), *b.src))
        elif kind == "B2":
            forced = b2_i < 4
            ref = str(rng.choice(B_CODONS[:2] if forced else B_CODONS))
            c = b.cassette(kind, "chr1", strand, base, ref, str(rng.choice(["K", "R"])))
            b.gene(c, 10, 90)
            c.in_gene = True
            c.peptides = b.first_peptides(1)
            if not forced and rng.random() < 0.75:  # missed cleavage at the site's K/R
                q = int(rng.integers(2, 7))
                c.codons[S - q - 1] = CODON[str(rng.choice(["K", "R"]))]
                c.peptides.append((q, q + int(rng.integers(5, 9)), *b.src))
            b2_i += 1
        elif kind == "B3":
            c = b.cassette(kind, "chr1", strand, base, str(rng.choice(b.alt_codons)),
                           b.prev_not_kr())  # fmt: skip
            q = int(rng.integers(1, 6))
            c.peptides = [*b.first_peptides(1), (q, q + int(rng.integers(5, 9)), *b.src)]
        elif kind == "control_start":
            ref = str(rng.choice(["GTG", "TTG"]))
            c = b.cassette(kind, "chr1", strand, base, ref, b.prev_not_kr())
            c.genes.append((S, S + 40))
            c.codons[S + 39] = "TAA"
            c.in_gene, c.initiator = True, True
            c.peptides = b.first_peptides(2)
        elif kind == "control_kM":
            c = b.cassette(kind, "chr1", strand, base, "AAG", b.prev_not_kr())
            c.peptides = b.first_peptides(2)
        else:  # same_db
            c = b.cassette(kind, "chr1", strand, base, "GTG", b.prev_not_kr())
            other = ("database", "other") if b.relation == "related" else ("denovo", "same")
            n = int(rng.integers(7, 14))
            c.peptides = [(0, n, "database", "same"), (0, n + 1, *other)]
        out.append(c)
    # B_origin on p1: minus strand, codon_start 30; the later peptide runs through the origin
    c = b.cassette("B_origin", "p1", "-", 30 + 3 * S + 3, str(rng.choice(b.alt_codons)),
                   b.prev_not_kr())  # fmt: skip
    c.crosses_origin = True
    c.peptides = [(0, 8, *b.src), (4, 20, *b.src)]
    out.append(c)
    return out


def _residue(b: _Builder, c: _Cassette, i: int) -> str:
    return b.code[c.codons[i]]


def _peptide(b: _Builder, c: _Cassette, q: int, n: int) -> str:
    a = S - q
    return "".join("M" if i == S else _residue(b, c, i) for i in range(a, a + n))


def alt_start_initiation(
    seed: int,
    *,
    size: str = "100k",
    genome_relation: str = "related",
    translation_table: int = 11,
    start_codons: Sequence[str] = ("ATG", "GTG", "TTG"),
) -> AltStartScenario:
    """The ``alt_start_initiation`` scenario (module docstring) for replicate ``seed``."""
    if size not in SIZES:
        raise ValueError(f"size {size!r} not in {sorted(SIZES)}")
    if genome_relation not in ("related", "same"):
        raise ValueError(f"genome_relation {genome_relation!r}")
    rng = np.random.default_rng(seed)
    b = _Builder(rng, translation_table, start_codons, genome_relation)
    lengths = {"chr1": SIZES[size], "p1": PLASMID_BP}
    cassettes = _plan(b, lengths["chr1"])
    seqs = {k: list(rng.choice(list("ACGT"), size=n)) for k, n in lengths.items()}
    for c in cassettes:
        L, seq = lengths[c.contig], seqs[c.contig]
        for i, cod in enumerate(c.codons):
            nt = cod if c.strand == "+" else revcomp(cod)
            lo = _low(c, i)
            for j, ch in enumerate(nt):
                seq[(lo + j) % L] = ch
    genome = {k: "".join(v) for k, v in seqs.items()}

    hits: list[dict[str, Any]] = []
    sdb: dict[str, str] = {}
    truth: list[dict[str, Any]] = []
    genes: list[Gene] = []
    n_gene = 0
    for c in cassettes:
        L = lengths[c.contig]
        for g0, g1 in c.genes:
            n_gene += 1
            start, end = _span(c, g0, g1, L)
            aa = "M" + "".join(_residue(b, c, i) for i in range(g0 + 1, g1 - 1))
            genes.append(Gene(f"ALTS_{n_gene:05d}", "cds", Interval(c.contig, start, end, c.strand),
                              aa, None, c.codons[g0]))  # type: ignore[arg-type]  # fmt: skip
        peps = []
        for q, n, source, search_db in c.peptides:
            pep = _peptide(b, c, q, n)
            if pep in sdb:
                raise ValueError(f"peptide {pep} planted twice (seed {seed})")
            start, end = _span(c, S - q, S - q + n, L)
            hits.append({
                "peptide": pep, "peptide_len_aa": n, "source": source, "contig": c.contig,
                "start": start, "end": end, "strand": c.strand,
                "n_diff": 0 if c.initiator else 1,
                "diff_types": pd.NA if c.initiator else "substitution",
                "placement_q": 0.0, "uniqueness": "unique",
            })  # fmt: skip
            sdb[pep] = search_db
            peps.append(pep)
        truth.append({
            "kind": c.kind, "contig": c.contig, "strand": c.strand,
            "codon_start": _low(c, S) % L, "ref_codon": c.ref_codon, "alt_aa": "M",
            "in_gene": c.in_gene, "prev_kr": _residue(b, c, S - 1) in ("K", "R"),
            "crosses_origin": c.crosses_origin,
            "all_first": all(q == 0 for q, *_ in c.peptides), "peptides": ";".join(peps),
        })  # fmt: skip
    reps = {k: Replicon(k, n, "circular") for k, n in lengths.items()}
    order = {k: i for i, k in enumerate(reps)}
    genes.sort(key=lambda g: (order[g.interval.contig], g.interval.start, g.interval.end,
                              g.interval.strand, g.locus_tag))  # fmt: skip
    gm = GeneModel(replicons=reps, genes=tuple(genes), features=())
    cfg = GippaConfig.model_validate({
        "genome_relation": genome_relation, "sample_name": "S1",
        "annotation": {"translation_table": translation_table,
                       "start_codons": list(start_codons)},
    })  # fmt: skip
    return AltStartScenario(
        genome=genome,
        gm=gm,
        cfg=cfg,
        peptide_hits=pd.DataFrame(hits),
        peptides=pd.DataFrame({"peptide": list(sdb), "search_db": list(sdb.values())}),
        truth=pd.DataFrame(truth, columns=list(TRUTH_COLUMNS)),
    )
