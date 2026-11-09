"""M12-05b: host and flanking genes of novel ORFs (docs/specs/M12-05b-novel-orf-genes.md AC1, AC2,
AC4, AC5; rules 2-7 and 10; lead decisions 2026-11-09).

Interface pinned here (test-author reading, beyond the spec text):

- ``gippa.novel.gene_context.gene_context(table, gm) -> pd.DataFrame`` returns the new columns of
  rule 2, in that order and only those, indexed like ``table``; the host is ``table``'s own
  ``overlapping_gene`` (rule 1). Distances are an integer dtype (nullable when NA occurs); the
  other cells are strings or NA.
- ``gippa.model.genes.gene_label(gene) -> str`` (rule 5) returns ``Gene.gene`` unchanged when it
  is usable.
- Rule 2's table lists nine new columns, though the spec text says "eight" and "34 columns"; the
  table is binding here (``helpers.m12_05b_context.NEW_COLUMNS``).

Worked examples whose span is not a whole number of codons cannot be a novel_orfs row (``_span``
always spans ``3 * length_aa + 3`` bases), so three are shifted to the nearest codon-whole span and
their numbers recomputed by rule 4 (the rule being tested is unchanged):

- E6: ``[9950, 10060)`` -> ``[9949, 10060)``; upstream 50 -> 49; the downstream 40 the spec
  checks in words is kept.
- E7: ``[9950, 9990)`` -> ``[9950, 9989)``; upstream 50 kept; circular downstream 110 -> 111.
- E10: ``[5500, 5700)`` -> ``[5500, 5701)``; the upstream tie at 100 kept; downstream 3800 -> 3799.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.m12_05b_context import (
    DISTANCES,
    E10_GENES,
    NEW_COLUMNS,
    cell,
    frame,
    gene,
    label,
    model,
    oracle,
    orf_row,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.model.coords import Interval, Replicon
from gippa.model.genes import Gene, GeneModel

L1, LP = 10_000, 3_000


def _gc(table: pd.DataFrame, gm: GeneModel) -> pd.DataFrame:
    from gippa.novel.gene_context import gene_context

    return gene_context(table, gm)


def _cells(out: pd.DataFrame) -> list[tuple[object, ...]]:
    return [tuple(cell(v) for v in r) for r in out[list(NEW_COLUMNS)].itertuples(index=False)]


def _x(ovl: Any, rel: str, other: Any, up: Any, down: Any) -> tuple[object, ...]:
    """Expected cells in the spec's ``label / relation / other ; up ; down`` notation."""
    return (ovl, rel, other, *up, *down)


NA3 = (None, None, None)

# (id, model kwargs, extra genes, row, expected cells after overlapping_gene)
EXAMPLES: list[tuple[str, dict[str, str], tuple[Gene, ...], dict[str, Any], tuple[object, ...]]] = [
    ("E1", {}, (), orf_row("c1", "+", 2600, 2900, "intergenic", None, L1),
     _x(None, "none", None, ("T_0003", "lacZ", 100), ("T_0004", "T_0004", 6600))),
    ("E2", {}, (), orf_row("c1", "-", 2600, 2900, "intergenic", None, L1),
     _x(None, "none", None, ("T_0004", "T_0004", 6600), ("T_0003", "lacZ", 100))),
    ("E3", {}, (), orf_row("c1", "-", 1400, 1700, "antisense", "T_0003", L1),
     _x("lacZ", "antisense", None, ("T_0004", "T_0004", 7800), ("T_0002", "T_0002", 300))),
    ("E4", {}, (), orf_row("c1", "+", 1000, 1450, "cds_outframe", "T_0003", L1),
     _x("lacZ", "same_strand", "T_0002", ("T_0001", "dnaA", 600), ("T_0004", "T_0004", 8050))),
    ("E5", {}, (), orf_row("c1", "-", 500, 1250, "nterm_extension", "T_0002", L1),
     _x("T_0002", "same_strand", None, ("T_0003", "lacZ", 50), ("T_0001", "dnaA", 100))),
    ("E6", {}, (), orf_row("c1", "+", 9949, 10060, "intergenic", None, L1),
     _x(None, "none", None, ("T_0004", "T_0004", 49), ("T_0001", "dnaA", 40))),
    ("E7", {"c1": "linear"}, (), orf_row("c1", "+", 9950, 9989, "intergenic", None, L1),
     _x(None, "none", None, ("T_0004", "T_0004", 50), NA3)),
    ("E7-circular", {}, (), orf_row("c1", "+", 9950, 9989, "intergenic", None, L1),
     _x(None, "none", None, ("T_0004", "T_0004", 50), ("T_0001", "dnaA", 111))),
    ("E8", {}, (), orf_row("p1", "+", 1500, 1800, "intergenic", None, LP),
     _x(None, "none", None, ("P_0001", "repA", 300), ("P_0001", "repA", 1400))),
    ("E8-linear", {"p1": "linear"}, (), orf_row("p1", "+", 1500, 1800, "intergenic", None, LP),
     _x(None, "none", None, ("P_0001", "repA", 300), NA3)),
    ("E9", {}, (), orf_row("c2", "+", 100, 400, "intergenic", None, 5_000),
     _x(None, "none", None, NA3, NA3)),
    ("E10", {}, E10_GENES, orf_row("c1", "+", 5500, 5701, "intergenic", None, L1),
     _x(None, "none", None, ("T_0005", "T_0005", 100), ("T_0004", "T_0004", 3799))),
]  # fmt: skip
IDS = [e[0] for e in EXAMPLES]


# ---------------------------------------------------------------- AC1: the worked examples


@pytest.mark.parametrize(("name", "kw", "extra", "row", "want"), EXAMPLES, ids=IDS)
def test_worked_example(
    name: str, kw: dict[str, str], extra: tuple[Gene, ...], row: dict[str, Any],
    want: tuple[object, ...],
) -> None:  # fmt: skip
    gm = model(extra=extra, **kw)
    out = _gc(frame([row]), gm)
    got = _cells(out)[0]
    labels = dict(zip(NEW_COLUMNS, got, strict=True))
    assert labels == dict(zip(NEW_COLUMNS, want, strict=True)), name


@pytest.mark.parametrize(("name", "kw", "extra", "row", "want"), EXAMPLES, ids=IDS)
def test_worked_example_agrees_with_the_oracle(
    name: str, kw: dict[str, str], extra: tuple[Gene, ...], row: dict[str, Any],
    want: tuple[object, ...],
) -> None:  # fmt: skip
    """The test oracle reproduces the spec's table (guards the oracle of the property below)."""
    gm = model(extra=extra, **kw)
    assert tuple(oracle(row, gm)[c] for c in NEW_COLUMNS) == want, name


def test_all_examples_in_one_table_keep_the_index_and_column_order() -> None:
    """One call over the E1-E6, E8, E9 rows (one model); a non-default index comes back as is."""
    picked = [e for e in EXAMPLES if not e[1] and not e[2]]
    table = frame([e[3] for e in picked], index=[70 + 3 * i for i in range(len(picked))])
    out = _gc(table, model())
    assert list(out.columns) == list(NEW_COLUMNS)
    assert out.index.equals(table.index)
    assert _cells(out) == [e[4] for e in picked]


def test_distances_are_integers_and_nullable() -> None:
    table = frame([e[3] for e in EXAMPLES if not e[1] and not e[2]])
    out = _gc(table, model())
    for c in DISTANCES:
        assert pd.api.types.is_integer_dtype(out[c].dtype), (c, out[c].dtype)
        assert out[c].isna().any()  # E9: an integer column that holds NA


def test_empty_table_gives_the_columns_and_no_rows() -> None:
    out = _gc(frame([]), model())
    assert list(out.columns) == list(NEW_COLUMNS) and len(out) == 0


def test_input_table_is_not_modified() -> None:
    table = frame([EXAMPLES[3][3]])
    before = table.copy()
    _gc(table, model())
    pd.testing.assert_frame_equal(table, before)


def test_relation_comes_from_the_strands_not_the_class() -> None:
    """Rule 2: a class that disagrees with the strands shows up in ``overlap_relation``."""
    row = orf_row("c1", "-", 1400, 1700, "cds_outframe", "T_0003", L1)  # T_0003 is on +
    assert _cells(_gc(frame([row]), model()))[0][1] == "antisense"
    row = orf_row("c1", "+", 1600, 1900, "antisense", "T_0003", L1)
    assert _cells(_gc(frame([row]), model()))[0][1] == "same_strand"


def test_one_base_overlap_is_an_overlap_and_never_a_flank() -> None:
    """Edge case: a gene overlapping the span by one base is an "other" gene, not a flank."""
    # [399, 501) touches T_0001's last base (399) and T_0002's first base (500)
    row = orf_row("c1", "+", 399, 501, "intergenic", None, L1)
    got = dict(zip(NEW_COLUMNS, _cells(_gc(frame([row]), model()))[0], strict=True))
    assert got["other_overlapping_genes"] == "T_0001;T_0002"
    assert got["upstream_gene"] == "T_0004" and got["downstream_gene"] == "T_0003"


def test_adjacent_gene_is_a_flank_at_distance_zero() -> None:
    row = orf_row("c1", "+", 400, 499, "intergenic", None, L1)  # T_0001 ends at 400, T_0002 at 500
    got = dict(zip(NEW_COLUMNS, _cells(_gc(frame([row]), model()))[0], strict=True))
    assert (got["upstream_gene"], got["upstream_distance_nt"]) == ("T_0001", 0)
    assert (got["downstream_gene"], got["downstream_distance_nt"]) == ("T_0002", 1)


def test_partial_orf_open_at_a_linear_end_has_no_flank_on_that_side() -> None:
    """Edge case: a 3'-open ORF ends at the contig end; the downstream flank is NA."""
    gm = model(c1="linear")
    # 3'-open at the linear end: no stop codon, so the span is 3 * length_aa (rule 3, ``_span``)
    row = orf_row("c1", "+", 9940, 10_000, "intergenic", None, L1)
    row |= {"end": 9999, "length_aa_longest": 20, "length_aa_best": 20, "flags": "partial"}
    got = dict(zip(NEW_COLUMNS, _cells(_gc(frame([row]), gm))[0], strict=True))
    assert (got["upstream_gene"], got["upstream_distance_nt"]) == ("T_0004", 40)
    assert got["downstream_gene"] is None and got["downstream_distance_nt"] is None


def test_no_start_codon_takes_the_span_from_the_first_peptide() -> None:
    """Rule 3: ``start_longest`` NA -> the span from the codon the ``orf_id`` names."""
    row = orf_row("c1", "+", 2600, 2900, "intergenic", None, L1) | {
        "start_longest": None,
        "flags": "no_start_codon",
    }
    got = dict(zip(NEW_COLUMNS, _cells(_gc(frame([row]), model()))[0], strict=True))
    assert (got["upstream_gene"], got["upstream_distance_nt"]) == ("T_0003", 100)
    assert (got["downstream_gene"], got["downstream_distance_nt"]) == ("T_0004", 6600)


def test_a_lone_host_is_never_its_own_flank() -> None:
    """Edge case: the host is never a flank; with no other gene the flanks are NA, on a circular
    and on a linear replicon."""
    reps = {"x": Replicon("x", 900, "circular"), "y": Replicon("y", 900, "linear")}
    gm = GeneModel(reps, (gene("X_1", "x", 100, 400, "+", "abc"),
                          gene("Y_1", "y", 100, 400, "+", "abc")), ())  # fmt: skip
    rows = [orf_row("x", "-", 200, 500, "antisense", "X_1", 900),
            orf_row("y", "-", 200, 500, "antisense", "Y_1", 900)]  # fmt: skip
    got = _cells(_gc(frame(rows), gm))
    assert got[0] == ("abc", "antisense", None, *NA3, *NA3)
    assert got[1] == ("abc", "antisense", None, *NA3, *NA3)


def test_noncoding_features_are_not_genes() -> None:
    """Rule 7 (lead decision Q1): a tRNA overlapping or next to the ORF is ignored."""
    from gippa.model.genes import NcFeature

    base = model()
    trna = NcFeature("trna_1", "trna", Interval("c1", 2550, 2700, "+"))
    gm = GeneModel(base.replicons, base.genes, (trna,))
    row = orf_row("c1", "+", 2600, 2900, "intergenic", None, L1)
    assert _cells(_gc(frame([row]), gm)) == [EXAMPLES[0][4]]


# ---------------------------------------------------------------- AC2: gene_label


@pytest.mark.parametrize(
    ("name", "want"),
    [("dnaA", "dnaA"), (None, "T_0009"), ("", "T_0009"), ("  ", "T_0009"), ("NA", "T_0009"),
     ("na", "T_0009"), ("Na", "T_0009"), ("\t", "T_0009")],
)  # fmt: skip
def test_gene_label(name: str | None, want: str) -> None:
    from gippa.model.genes import gene_label

    assert gene_label(gene("T_0009", "c1", 0, 300, "+", name)) == want


_TAG = st.from_regex(r"[A-Z]{1,4}_[0-9]{1,5}", fullmatch=True)
_NAME = st.one_of(
    st.none(),
    st.sampled_from(["", " ", "  ", "\t", "NA", "na", "nA", " NA ", "dnaA", "lacZ_1"]),
    st.text(alphabet="abcdefghijklmnopqrstuvwxyzABCN0123456789_-", min_size=1, max_size=8),
)


@given(_TAG, _NAME)
def test_gene_label_property(tag: str, name: str | None) -> None:
    """Never the literal NA (any case) and never blank; the name itself when usable."""
    from gippa.model.genes import gene_label

    g = gene(tag, "c1", 0, 300, "+", name)
    got = gene_label(g)
    assert got.strip() and got.strip().upper() != "NA"
    assert got == label(g)


# ---------------------------------------------------------------- AC4: write -> read


def _all_examples() -> pd.DataFrame:
    """Every example's row with its expected cells filled in (one frame, mixed models)."""
    rows = []
    for _, _, _, row, want in EXAMPLES:
        rows.append(row | dict(zip(NEW_COLUMNS, want, strict=True)))
    return frame(rows)


def test_round_trip_holds_every_example(tmp_path: Any) -> None:
    from gippa.io.novel_tables import read_novel_orfs, write_novel_orfs
    from gippa.io.tsv import read_tsv

    table = _all_examples()
    assert set(NEW_COLUMNS) <= set(table.columns), "NOVEL_ORFS_COLUMNS lacks the new columns"
    path = tmp_path / "novel_orfs.tsv"
    write_novel_orfs(path, table)
    back = read_novel_orfs(path)
    pd.testing.assert_frame_equal(back, table)
    header, rows = read_tsv(path)
    written = [{c: r[header.index(c)] for c in NEW_COLUMNS} for r in rows]
    for got, (name, *_, want) in zip(written, EXAMPLES, strict=True):
        exp = {c: None if v is None else str(v) for c, v in zip(NEW_COLUMNS, want, strict=True)}
        assert got == exp, name  # distances as computed: not shifted by 1 (rule 10)


def test_writer_adds_no_position_column() -> None:
    from gippa.novel.outputs import POSITION_COLUMNS

    assert POSITION_COLUMNS == ("start_longest", "start_best_ostir", "end")


# ---------------------------------------------------------------- AC5 and the oracle property


def _invariants(out: pd.DataFrame, table: pd.DataFrame) -> None:
    for i in table.index:
        ovl, lab, rel = table.at[i, "overlapping_gene"], out.at[i, "overlapping_gene_label"], (
            out.at[i, "overlap_relation"])  # fmt: skip
        assert (cell(ovl) is None) == (cell(lab) is None), i
        assert (cell(ovl) is None) == (rel == "none"), i
        assert rel in ("same_strand", "antisense", "none"), rel
        for side in ("upstream", "downstream"):
            tag, lab2, d = (cell(out.at[i, f"{side}_{c}"]) for c in ("gene", "gene_label",
                                                                    "distance_nt"))  # fmt: skip
            assert (tag is None) == (lab2 is None) == (d is None), (i, side)
            assert d is None or d >= 0
            assert lab2 != "NA"


def test_ac5_on_the_examples() -> None:
    table = frame([e[3] for e in EXAMPLES if not e[1] and not e[2]])
    _invariants(_gc(table, model()), table)


@st.composite
def cases(draw: st.DrawFn) -> tuple[GeneModel, pd.DataFrame]:
    """A small replicon (circular or linear), up to 7 genes with shuffled tags (so model order and
    tag order differ), a gene-bearing second contig, and 1-4 ORFs with consistent hosts."""
    n = draw(st.integers(60, 240))
    topo = draw(st.sampled_from(["circular", "linear"]))
    reps = {"c1": Replicon("c1", n, topo), "c9": Replicon("c9", 500, "circular")}  # type: ignore[arg-type]
    k = draw(st.integers(0, 7))
    tags = draw(st.permutations([f"G_{i:04d}" for i in range(k)]))
    names = draw(st.lists(st.sampled_from([None, "NA", "", "abc", "xyz"]), min_size=k, max_size=k))
    ivs = []
    for _ in range(k):
        s = draw(st.integers(0, n - 2))
        ln = draw(st.integers(1, (n - s) if topo == "linear" else n - 1))
        ivs.append((s, s + ln, draw(st.sampled_from(["+", "-"]))))
    genes = [gene(t, "c1", s, e, d, nm) for t, (s, e, d), nm in zip(tags, ivs, names, strict=True)]
    rank = {"+": 0, "-": 1}
    genes.sort(key=lambda g: (g.interval.start, g.interval.end, rank[g.interval.strand],
                              g.locus_tag))  # fmt: skip
    genes.append(gene("Z_0001", "c9", 10, 400, "+", "far"))  # never seen from c1
    gm = GeneModel(reps, tuple(genes), ())
    rows = []
    for _ in range(draw(st.integers(1, 4))):
        strand = draw(st.sampled_from(["+", "-"]))
        codons = draw(st.integers(2, (n - 3) // 3))
        width = 3 * codons
        lo = draw(st.integers(0, (n - width) if topo == "linear" else n - 1))
        row = orf_row("c1", strand, lo, lo + width, "intergenic", None, n)
        over = [g for g in genes[:-1] if _hits(g, lo, lo + width, n, topo)]
        if over:
            host = min(over, key=lambda g: g.locus_tag)
            same = host.interval.strand == strand
            row |= {"overlapping_gene": host.locus_tag,
                    "class": "cds_outframe" if same else "antisense"}  # fmt: skip
        rows.append(row)
    index = draw(st.lists(st.integers(-50, 50), min_size=len(rows), max_size=len(rows),
                          unique=True))  # fmt: skip
    return gm, frame(rows, index=index)


def _hits(g: Gene, lo: int, hi: int, n: int, topo: str) -> bool:
    def bases(a: int, b: int) -> set[int]:
        return {p % n for p in range(a, b)} if topo == "circular" else set(range(a, b))

    return bool(bases(g.interval.start, g.interval.end) & bases(lo, hi))


@settings(max_examples=150, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(cases())
def test_gene_context_matches_the_brute_force_oracle(case: tuple[GeneModel, pd.DataFrame]) -> None:
    """Rules 2-5 against position-set geometry: overlap, flanks, wrap, linear ends, ties, labels,
    the minus-strand swap, other replicons never reached; AC5 on every row."""
    gm, table = case
    out = _gc(table, gm)
    assert list(out.columns) == list(NEW_COLUMNS)
    assert out.index.equals(table.index)
    for i, r in zip(table.index, table.to_dict("records"), strict=True):
        want = oracle(r, gm)
        got = {c: cell(out.at[i, c]) for c in NEW_COLUMNS}
        assert got == want, (i, r["strand"], r["start_longest"], r["end"])
    _invariants(out, table)


@settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(cases())
def test_minus_strand_mirror_swaps_the_sides(case: tuple[GeneModel, pd.DataFrame]) -> None:
    """The same span on the other strand swaps upstream and downstream; distances unchanged."""
    gm, table = case
    flipped = table.copy()
    for i in table.index:
        n = gm.replicons[str(table.at[i, "contig"])].length
        r = table.loc[i].to_dict()
        lo = int(r["start_longest"]) if r["strand"] == "+" else int(r["end"])
        hi = lo + 3 * int(r["length_aa_longest"]) + 3
        s = "-" if r["strand"] == "+" else "+"
        start, end = (lo, (hi - 1) % n) if s == "+" else ((hi - 1) % n, lo)
        flipped.at[i, "strand"], flipped.at[i, "start_longest"], flipped.at[i, "end"] = (
            s, start, end)  # fmt: skip
    a, b = _gc(table, gm), _gc(flipped, gm)
    for side, other in (("upstream", "downstream"), ("downstream", "upstream")):
        for c in ("gene", "gene_label", "distance_nt"):
            assert [cell(v) for v in a[f"{side}_{c}"]] == [cell(v) for v in b[f"{other}_{c}"]]
    assert [cell(v) for v in a["other_overlapping_genes"]] == [
        cell(v) for v in b["other_overlapping_genes"]]  # fmt: skip
