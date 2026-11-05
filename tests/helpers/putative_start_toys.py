"""Hand-built worlds for the M12-04 PR 10 ``putative_start`` tests
(docs/specs/M12-04-outputs-revision.md "PR 10"; statistician sign-off
``local/handoff/pr10_statistician.md`` C1-C5). Test tooling only: nothing here imports the rule
under test.

A world is one contig ``c1`` read in **transcription codons** ``0 .. n_codons - 1`` on one strand:

- ``+``: codon ``i`` is genome ``[off + 3 i, off + 3 i + 3)`` (unwrapped; taken mod L on a circular
  contig, so ``off`` places the origin);
- ``-``: codon ``i`` is genome ``[end - 3 i - 3, end - 3 i)``; ``end`` defaults to L.

Every codon holds a filler residue (no K, R, M, I, L, V or C, so the only K/R before a site is the
one a test plants). A :class:`Site` writes its reference codon at codon ``s`` (and the residue
``prev`` at codon ``s - 1``) and places its ALT peptides: peptide ``(q, n)`` covers codons
``[s - q, s - q + n)`` and reads ``M`` at the site, so the M is at peptide position ``q`` in
transcription (peptide) order. Hits are ``peptide_hits`` rows (0-based half-open, unwrapped:
``start`` in ``[0, L)``) with ``n_diff = 1``. Genes are CDSs with an empty ``aa``, so M11-20 does
not credit their starts (``annotated_starts`` skips them) and a site at a gene's first codon is an
ALT row with ``start_codon``.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from gippa.config import GippaConfig
from gippa.mapping.sixframe import translate_text
from gippa.model import Replicon
from helpers.readthrough_toys import CODON, gene, model, rc

FILL = "AGSTEDQNPFYHW"
DEFAULT_PEPTIDES: tuple[tuple[int, int], ...] = ((0, 8), (0, 11))
DB_OTHER = ("database", "other")


@dataclass(frozen=True)
class Site:
    """A planted ALT site at transcription codon ``s``.

    ``peptides``: ``(q, length)`` per ALT peptide (the M at peptide position ``q``). ``kind``
    ``swap``: codon ``s + 1`` is ATG and the peptides read ``MV`` over the genome's ``VM``
    (``ref_codon`` must be GTG). ``sources``: ``(source, search_db)`` per peptide (default all
    ``database`` / ``other``)."""

    s: int
    ref_codon: str = "GTG"
    prev: str | None = "G"
    peptides: tuple[tuple[int, int], ...] = DEFAULT_PEPTIDES
    kind: str = "substitution"
    sources: tuple[tuple[str, str], ...] | None = None
    ref_peptides: tuple[tuple[int, int], ...] = ()  # (q, length): reference reads, n_diff 0


@dataclass(frozen=True)
class World:
    n_codons: int = 120
    strand: str = "+"
    topology: str = "linear"
    off: int | None = None  # + : genome start of codon 0 (default 0); - : genome end of codon 0
    #: (tag, first codon, end codon[, feature_type]) in transcription codons; type default "cds"
    genes: tuple[tuple[Any, ...], ...] = ()
    sites: tuple[Site, ...] = ()
    table: int = 11
    seed: int = 0
    extra_codons: dict[int, str] = field(default_factory=dict)  # transcription codon -> nt

    @property
    def length(self) -> int:
        return 3 * self.n_codons

    def _base(self) -> int:
        if self.off is not None:
            return self.off
        return 0 if self.strand == "+" else self.length

    def low(self, i: int) -> int:
        """Unwrapped lowest genome coordinate of transcription codon ``i``."""
        b = self._base()
        return b + 3 * i if self.strand == "+" else b - 3 * i - 3

    def codon_start(self, i: int) -> int:
        """``variants.tsv`` ``codon_start`` (0-based, wrapped) of transcription codon ``i``."""
        return self.low(i) % self.length

    def span(self, a: int, b: int) -> tuple[int, int]:
        """Genome ``[start, end)`` of transcription codons ``[a, b)``; ``start`` in ``[0, L)``."""
        lo = min(self.low(a), self.low(b - 1))
        hi = max(self.low(a), self.low(b - 1)) + 3
        shift = (lo // self.length) * self.length
        return lo - shift, hi - shift


def _codon_table(w: World) -> list[str]:
    import numpy as np

    rng = np.random.default_rng(w.seed + 17)
    codons = [CODON[str(a)] for a in rng.choice(list(FILL), size=w.n_codons)]
    for i, nt in w.extra_codons.items():
        codons[i] = nt
    for st in w.sites:
        codons[st.s] = st.ref_codon
        if st.prev is not None and st.s >= 1:
            codons[st.s - 1] = CODON[st.prev]
        if st.kind == "swap":
            codons[st.s + 1] = "ATG"
    return codons


def genome_of(w: World) -> str:
    codons = _codon_table(w)
    seq = ["N"] * w.length
    for i, c in enumerate(codons):
        nt = c if w.strand == "+" else rc(c)
        lo = w.low(i)
        for j, ch in enumerate(nt):
            seq[(lo + j) % w.length] = ch
    return "".join(seq)


def peptide_of(w: World, st: Site, q: int, n: int, *, alt: bool = True) -> str:
    """The peptide over codons ``[s - q, s - q + n)``: the ALT read, or (``alt=False``) the
    reference read."""
    codons = _codon_table(w)
    a = st.s - q
    res = []
    for i in range(a, a + n):
        aa = translate_text(codons[i % w.n_codons], w.table)
        if not alt:
            pass
        elif i == st.s:
            aa = "M"
        elif st.kind == "swap" and i == st.s + 1:
            aa = translate_text(st.ref_codon, w.table)
        res.append(aa)
    return "".join(res)


def hits_of(w: World) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    sdb: dict[str, str] = {}
    for st in w.sites:
        srcs = st.sources or (DB_OTHER,) * len(st.peptides)
        for (q, n), (source, search_db) in zip(st.peptides, srcs, strict=True):
            pep = peptide_of(w, st, q, n)
            a = st.s - q
            start, end = w.span(a, a + n)
            rows.append({
                "peptide": pep, "peptide_len_aa": len(pep), "source": source, "contig": "c1",
                "start": start, "end": end, "strand": w.strand, "n_diff": 1,
                "diff_types": "swap" if st.kind == "swap" else "substitution",
                "placement_q": 0.0, "uniqueness": "unique",
            })  # fmt: skip
            sdb[pep] = search_db
        for q, n in st.ref_peptides:
            pep = peptide_of(w, st, q, n, alt=False)
            start, end = w.span(st.s - q, st.s - q + n)
            rows.append({
                "peptide": pep, "peptide_len_aa": len(pep), "source": "database", "contig": "c1",
                "start": start, "end": end, "strand": w.strand, "n_diff": 0,
                "diff_types": pd.NA, "placement_q": 0.0, "uniqueness": "unique",
            })  # fmt: skip
            sdb[pep] = "other"
    hits = pd.DataFrame(rows)
    assert hits["peptide"].is_unique, "toy peptides must be distinct"
    peps = pd.DataFrame({"peptide": list(sdb), "search_db": list(sdb.values())})
    return hits, peps


def gene_model_of(w: World) -> Any:
    reps = {"c1": Replicon("c1", w.length, w.topology)}  # type: ignore[arg-type]
    genes = []
    for g in w.genes:
        tag, g0, g1 = g[:3]
        start, end = w.span(g0, g1)
        genes.append(gene(tag, "c1", start, end, w.strand, g[3] if len(g) > 3 else "cds"))
    return model(reps, genes)


def cfg_of(
    relation: str = "related", *, table: int = 11, start_codons: Sequence[str] | None = None,
    high_min_peptides: int = 2,
) -> GippaConfig:  # fmt: skip
    ann: dict[str, Any] = {"translation_table": table}
    if start_codons is not None:
        ann["start_codons"] = list(start_codons)
    return GippaConfig.model_validate({
        "genome_relation": relation, "sample_name": "S1", "annotation": ann,
        "variants": {"high_min_peptides": high_min_peptides},
    })  # fmt: skip


def run_world(w: World, cfg: GippaConfig | None = None, spectra: pd.DataFrame | None = None
              ) -> tuple[pd.DataFrame, list[Any], pd.DataFrame]:  # fmt: skip
    """``run_variants_qc`` on the world: ``(variants table, VCF records, support frame)``."""
    from gippa.steps.variants import run_variants_qc

    hits, peps = hits_of(w)
    cfg = cfg or cfg_of(table=w.table)
    _pile, table, records, support = run_variants_qc(
        hits, None, peps, gene_model_of(w), {"c1": genome_of(w)}, cfg, spectra
    )
    return table, records, support


def row_at(table: pd.DataFrame, w: World, s: int, alt: str = "M") -> dict[str, Any]:
    """The one ``variants.tsv`` row at transcription codon ``s`` with ALT ``alt``."""
    sel = table[
        (table["codon_start"].astype(int) == w.codon_start(s))
        & (table["strand"] == w.strand)
        & (table["alt_aa"].astype(str) == alt)
    ]
    assert len(sel) == 1, (s, w.codon_start(s), table[["codon_start", "strand", "alt_aa"]])
    return sel.to_dict("records")[0]


def flags(row: dict[str, Any]) -> tuple[str, ...]:
    v = row["flags"]
    return tuple(v) if isinstance(v, (tuple, list)) else tuple(str(v).split(";"))
