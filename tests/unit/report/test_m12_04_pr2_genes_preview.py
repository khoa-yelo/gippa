"""M12-04 PR 2 behaviour 3, AC 2.3: the curated genes preview
(docs/specs/M12-04-outputs-revision.md § PR 2).

The preview (``<table data-table="genes_tsv">`` in the Gene detection section) shows exactly
``PREVIEW_COLUMNS``, the first 20 genes in this order: ``abundance_log2`` descending, NA last,
ties by ``locus_tag``; when the run has no abundance (every ``abundance_log2`` NA), ``detected``
descending, then ``n_peptides_unique`` descending, then ``locus_tag``. The link to the full TSV
stays and ``tables/genes.tsv`` itself is unchanged.

Interface assumed (PR 2 test author): ``locus_tag`` ties break in code-point order, ascending (as
``sorted`` and ``DataFrame.sort_values`` do).
"""

from __future__ import annotations

import math
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.report_fixtures import Page, cfg, parse, table_rows, write_preflight
from helpers.viz_out_tree import gene_row, make_out, write_schema_table
from hypothesis import HealthCheck, example, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import build

PREVIEW_COLUMNS = ["locus_tag", "gene", "product", "detection_status", "n_peptides_unique",
                   "abundance_log2", "region_ids", "flags"]  # fmt: skip
TOP = 20


def expected_order(rows: list[dict[str, Any]]) -> list[str]:
    """The spec's sort, as an oracle."""

    def na(v: Any) -> bool:
        return v is None or (isinstance(v, float) and math.isnan(v))

    def by_abundance(r: dict[str, Any]) -> tuple[Any, ...]:
        a = r["abundance_log2"]
        return (na(a), 0.0 if na(a) else -a, r["locus_tag"])

    def by_detection(r: dict[str, Any]) -> tuple[Any, ...]:
        return (not r["detected"], -(r["n_peptides_unique"] or 0), r["locus_tag"])

    key = by_abundance if any(not na(r["abundance_log2"]) for r in rows) else by_detection
    return [r["locus_tag"] for r in sorted(rows, key=key)][:TOP]


def _preview(page: Page) -> tuple[list[str], list[list[str]]]:
    table = page.table("genes_tsv")
    assert table is not None, "no genes preview"
    det = [s for s in page.sections() if s.attrs.get("id") == "detection"]
    assert det and table in list(det[0].iter()), "the genes preview is not in Gene detection"
    return table_rows(table)


def _report(out: Path) -> Page:
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    return parse(build_report(RunView.from_outdir(out), {}, {}, cfg()))


# ---------------------------------------------------------------- fixtures P and D


@pytest.mark.parametrize("name", ["P", "D"])
def test_preview_has_exactly_the_eight_columns(name: str, tmp_path: Path) -> None:
    header, body = _preview(build(name, tmp_path).page())
    assert header == PREVIEW_COLUMNS
    assert len(body) == TOP
    assert all(len(r) == len(PREVIEW_COLUMNS) for r in body)


def test_preview_sort_with_abundance_fixture_p(tmp_path: Path) -> None:
    """P: abundance 20.0 + 0.01 i on the 99 detected genes, so g00098 comes first."""
    _, body = _preview(build("P", tmp_path).page())
    assert [r[0] for r in body] == [f"g{i:05d}" for i in range(98, 78, -1)]


def test_preview_sort_without_abundance_fixture_d(tmp_path: Path) -> None:
    """D (no abundance): the 4 detected genes first (spec example, run C)."""
    _, body = _preview(build("D", tmp_path).page())
    tags = [r[0] for r in body]
    assert tags[:4] == ["g00000", "g00001", "g00002", "g00003"], tags
    assert tags[4:] == [f"g{i:05d}" for i in range(4, TOP)], tags


def test_full_table_link_stays_and_file_unchanged(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    path = run.out / "tables" / "genes.tsv"
    before = path.read_bytes()
    page = run.page()
    det = next(s for s in page.sections() if s.attrs.get("id") == "detection")
    assert "../tables/genes.tsv" in [a.attrs.get("href") for a in det.by_tag("a")]
    assert path.read_bytes() == before


# ---------------------------------------------------------------- property


_TAG = st.text(alphabet="abcxyzABC019_", min_size=1, max_size=5).map(lambda s: f"g{s}")
_ROW = st.fixed_dictionaries({
    "detected": st.booleans(),
    "n_peptides_unique": st.sampled_from([0, 1, 2, 3]),
    "abundance_log2": st.sampled_from([None, None, 10.0, 12.5, 12.5, 15.25]),
})  # fmt: skip


@st.composite
def gene_tables(draw: st.DrawFn) -> tuple[list[dict[str, Any]], bool]:
    tags = draw(st.lists(_TAG, min_size=1, max_size=30, unique=True))
    rows = [{"locus_tag": t, **draw(_ROW)} for t in tags]
    no_abundance = draw(st.booleans())
    for r in rows:
        if no_abundance:
            r["abundance_log2"] = None
        if not r["detected"]:
            r["n_peptides_unique"] = 0
    return rows, no_abundance


def _write(out: Path, rows: list[dict[str, Any]], no_abundance: bool) -> None:
    body = [gene_row(locus_tag=r["locus_tag"], start=101 + 1000 * i, end=400 + 1000 * i,
                     detected=r["detected"],
                     detection_status="detected" if r["detected"] else "not_detected",
                     detection_tier="two_peptide", n_peptides_unique=r["n_peptides_unique"],
                     abundance_log2=r["abundance_log2"],
                     abundance_type="none" if no_abundance else "ibaq",
                     abundance_scale=None if no_abundance else "molar",
                     product=f"product {i}", gene=f"gen{i}")
            for i, r in enumerate(rows)]  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", body)


@settings(max_examples=40, deadline=None,
          suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large])  # fmt: skip
@given(gene_tables())
@example(([{"locus_tag": "gb", "detected": True, "n_peptides_unique": 2, "abundance_log2": 12.5},
           {"locus_tag": "ga", "detected": True, "n_peptides_unique": 1, "abundance_log2": 12.5},
           {"locus_tag": "gc", "detected": False, "n_peptides_unique": 0,
            "abundance_log2": None}], False))  # fmt: skip
@example(([{"locus_tag": "gb", "detected": True, "n_peptides_unique": 1, "abundance_log2": None},
           {"locus_tag": "gz", "detected": True, "n_peptides_unique": 3, "abundance_log2": None},
           {"locus_tag": "ga", "detected": False, "n_peptides_unique": 0,
            "abundance_log2": None}], True))  # fmt: skip
def test_preview_sort_property(table: tuple[list[dict[str, Any]], bool]) -> None:
    rows, no_abundance = table
    with tempfile.TemporaryDirectory() as tmp:
        out = make_out(Path(tmp))
        write_preflight(out)
        _write(out, rows, no_abundance)
        header, body = _preview(_report(out))
    assert header == PREVIEW_COLUMNS
    assert [r[0] for r in body] == expected_order(rows)
