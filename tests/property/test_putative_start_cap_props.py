"""M12-04 PR 10 property (sign-off C3, C5): **the cap never promotes and never changes counts**.

On random hand-built worlds (``helpers.putative_start_toys``: either strand, linear or circular
with a random origin, an optional gene, 1-4 planted ALT sites with random reference codon,
preceding residue, peptide positions and sources), the ``variants`` chain
(``gippa.steps.variants.run_variants_qc``) is compared with the same chain without the rule (the
M7-05 flags chain: ``sites`` -> ``flags_for`` / ``left_block_site`` -> ``annotate_sites`` ->
``vcf_records``):

- the same rows (``variant_id``, position, ALT, ``call``, ``n_ref``, ``n_alt``, class error,
  ``nt_changes``, ``alt_codon_vcf``, ``source``);
- flags equal except ``putative_start``, which only substitutions to M at a non-ATG codon of
  ``annotation.start_codons`` carry;
- confidence equal except on ``putative_start`` rows, which are ``candidate``; a ``high`` row is
  ``high`` without the rule (never promotes);
- the flagged rows are exactly those of the independent oracle (``helpers.putative_start_oracle``);
- VCF ``CONF`` equals the TSV ``confidence`` of every represented row (C5).
"""

from __future__ import annotations

from typing import Any

import pandas as pd
from helpers.putative_start_oracle import chain_alignments, oracle_flagged, row_flags
from helpers.putative_start_toys import Site, World, cfg_of, gene_model_of, genome_of, hits_of
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

FLAG = "putative_start"
SLOTS = (15, 45, 75, 100)
SOURCES = st.sampled_from([("database", "other"), ("denovo", "same"), ("database", "same")])
SAME_COLS = ["variant_id", "contig", "codon_start", "strand", "gene", "aa_pos", "ref_codon",
             "ref_aa", "alt_aa", "alt_type", "call", "n_ref", "n_alt", "nt_changes",
             "alt_codon_vcf", "class_error_rate", "class_error_n", "source"]  # fmt: skip


@st.composite
def sites(draw: st.DrawFn, s: int) -> Site:
    pairs = draw(st.lists(st.tuples(st.integers(0, 3), st.integers(7, 12)), min_size=1,
                          max_size=3, unique=True))  # fmt: skip
    return Site(
        s,
        ref_codon=draw(st.sampled_from(["GTG", "TTG", "CTG", "AAG", "GCT", "TTA"])),
        prev=draw(st.sampled_from(["G", "K", "R", "A"])),
        peptides=tuple(pairs),
        sources=tuple(draw(SOURCES) for _ in pairs),
    )


@st.composite
def worlds(draw: st.DrawFn) -> World:
    strand = draw(st.sampled_from(["+", "-"]))
    topology = draw(st.sampled_from(["linear", "circular"]))
    off = None
    if topology == "circular":
        off = draw(st.integers(0, 359))
    slots = draw(st.lists(st.sampled_from(SLOTS), min_size=1, max_size=4, unique=True))
    genes: tuple[tuple[Any, ...], ...] = ()
    if draw(st.booleans()):
        g0 = draw(st.sampled_from([5, 15, 45]))
        g1 = draw(st.sampled_from([60, 90, 118]))
        genes = (("G1", g0, g1),)
    return World(strand=strand, topology=topology, off=off, genes=genes,
                 sites=tuple(draw(sites(s)) for s in sorted(slots)))  # fmt: skip


def _uncapped(w: World, cfg: Any) -> tuple[pd.DataFrame, pd.DataFrame, list[Any]]:
    """The chain without the rule; returns (table, support, alignments)."""
    from gippa.steps.variants import support_of
    from gippa.variants.error_rates import annotate_sites
    from gippa.variants.flags import flags_for, left_block_site, start_codons
    from gippa.variants.pileup import eligible, pileup, sites
    from gippa.variants.qc import site_support
    from gippa.variants.vcf_records import vcf_records

    hits, peps = hits_of(w)
    gm, genome = gene_model_of(w), {"c1": genome_of(w)}
    support = support_of(hits, peps, {g.locus_tag for g in gm.genes})
    aligns = chain_alignments(hits, genome, gm, cfg)
    pile = pileup(eligible(aligns, support, cfg), gm, genome, cfg, support=support)
    sts = sites(pile, cfg)
    alt_peps = sts.attrs["alt_peptides"]
    starts, tol, tt = start_codons(gm), cfg.variants.mass_equivalent_tol_da, w.table
    out = []
    for r, ps in zip(sts.to_dict("records"), alt_peps, strict=True):
        al = [a for a in aligns if a.peptide in ps]
        left = left_block_site(r, al, genome, gm.replicons, tol, table=tt)
        src = {(support[p].source, support[p].search_db) for p in ps}
        out.append(tuple(flags_for(r, starts, cfg, src, left_block=left)))
    sts["flags"] = out
    sts = annotate_sites(sts, [])
    _, table = vcf_records(sts, pile, genome, gm.replicons, table=tt)
    sources = {p: (str(s.source), str(s.search_db)) for p, s in support.items()}
    return table, site_support(table["variant_id"], alt_peps, sources), aligns


@settings(max_examples=80, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    w=worlds(),
    relation=st.sampled_from(["related", "same"]),
    k=st.integers(1, 3),
    ctg=st.booleans(),
)
def test_the_cap_never_promotes_and_never_changes_counts(
    w: World, relation: str, k: int, ctg: bool
) -> None:
    from gippa.steps.variants import run_variants_qc

    starts = ["ATG", "GTG", "TTG", "CTG"] if ctg else None
    cfg = cfg_of(relation, start_codons=starts, high_min_peptides=k)
    hits, peps = hits_of(w)
    gm, genome = gene_model_of(w), {"c1": genome_of(w)}
    _, got, records, _ = run_variants_qc(hits, None, peps, gm, genome, cfg)
    base, support, aligns = _uncapped(w, cfg)

    assert len(got) == len(base)
    pd.testing.assert_frame_equal(
        got[SAME_COLS].reset_index(drop=True), base[SAME_COLS].reset_index(drop=True)
    )
    start_set = {c for c in cfg.annotation.start_codons} - {"ATG"}
    for g, b in zip(got.to_dict("records"), base.to_dict("records"), strict=True):
        gf, bf = row_flags(g["flags"]), row_flags(b["flags"])
        assert set(gf) - {FLAG} == set(bf), (gf, bf)
        if FLAG in gf:
            assert (g["alt_type"], g["alt_aa"]) == ("substitution", "M"), g
            assert str(g["ref_codon"]) in start_set, g
            assert "start_codon" not in gf, g
            assert g["confidence"] == "candidate", g
        else:
            assert g["confidence"] == b["confidence"], (g, b)
        if g["confidence"] == "high":
            assert b["confidence"] == "high", (g, b)

    want = oracle_flagged(base, support, aligns, genome, gm, cfg)
    flagged = {str(v) for v, f in zip(got["variant_id"], got["flags"], strict=True)
               if FLAG in row_flags(f)}  # fmt: skip
    assert flagged == want

    conf = dict(zip(got["variant_id"].astype(str), got["confidence"], strict=True))
    by_pos: dict[tuple[str, int, str], str] = {}
    for g in got.to_dict("records"):
        by_pos[(str(g["contig"]), int(g["codon_start"]), str(g["alt_aa"]))] = str(g["variant_id"])
    for rec in records:
        for alt, c in zip(rec.info["AA_ALT"], rec.info["CONF"], strict=True):  # type: ignore[call-overload]
            vid = by_pos.get((rec.contig, rec.pos0, str(alt)))
            if vid is not None:
                assert c == conf[vid], (rec, conf[vid])
