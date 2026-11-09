"""M12-05b: the novel_orfs columns, schema, old trees, report previews and written runs
(docs/specs/M12-05b-novel-orf-genes.md AC3, AC5-AC9, rules 2, 8, 9; lead decisions 2026-11-09).

Interface pinned here (test-author reading, beyond the spec text):

- Rule 2's table lists nine new columns although the text says "eight" and "34 columns"; the
  table is binding (``helpers.m12_05b_context.NEW_COLUMNS``, 35 columns in all).
- ``gippa.io.migration.ADDED_COLUMNS["novel_orfs"]`` is that tuple, in order; ``RunView`` fills
  every absent one with NA (``overlap_relation`` included) and keeps the schema order.
- The Novel ORFs preview is ``<table data-table="novel_orfs">``; ``gippa.report.build`` holds
  ``NOVEL_PREVIEW``. Its rows are the first ``TOP_ROWS`` rows in table order.
- Genes preview (lead decision "Locus tag in italics"): a ``gene`` cell without a usable name is
  exactly one ``<em>`` whose text is the locus tag; the tooltip is a ``title`` attribute on that
  ``<em>`` or on its ``<td>`` and contains the product. The UniRef50 part is not pinned: no file
  that ``RunView`` reads carries it (``genes.tsv`` has no UniRef column), so its source is the
  implementer's to choose.
- Written runs (AC5, AC6) are the regenerated fixtures ``tests/data/snapshot/*/OUT`` and the
  golden ``tests/data/golden/sim50k_m7/novel_orfs.tsv``, read by their raw header (no fill), so
  a file that was not regenerated fails.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.m12_04_renames import to_1_0
from helpers.m12_05b_context import COLUMNS, DISTANCES, NEW_COLUMNS, RC1_COLUMNS, TAG_COLUMNS
from helpers.report_fixtures import Node, Page, cfg, parse, run_out, table_rows
from helpers.viz_out_tree import gene_row, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.io.tsv import read_tsv

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "src" / "gippa" / "schemas" / "tables" / "novel_orfs.schema.json"
SNAPSHOT = REPO / "tests" / "data" / "snapshot"
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k_m7" / "novel_orfs.tsv"
ROUTE_E = Path("local/example/routes/1945698/gippa-example/runs/E/tables/novel_orfs.tsv")
NOVEL_PREVIEW = ("orf_id", "class", "called", "class_fdr", "strand", "length_aa_best",
                 "overlapping_gene_label", "overlap_relation", "upstream_gene_label",
                 "upstream_distance_nt", "downstream_gene_label", "downstream_distance_nt",
                 "n_peptides", "flags")  # fmt: skip
TYPES = {c: "string" for c in NEW_COLUMNS} | dict.fromkeys(DISTANCES, "integer")


def _fields() -> dict[str, dict[str, Any]]:
    doc = json.loads(SCHEMA.read_text(encoding="utf-8"))
    return {f["name"]: f for f in doc["fields"]}


# ---------------------------------------------------------------- AC3: columns and schema


def test_novel_orfs_columns_insert_the_new_columns_after_overlapping_gene() -> None:
    from gippa.novel.outputs import NOVEL_ORFS_COLUMNS

    assert tuple(NOVEL_ORFS_COLUMNS) == COLUMNS
    assert COLUMNS[-1] == "flags" and len(COLUMNS) == len(RC1_COLUMNS) + len(NEW_COLUMNS)


def test_novel_orfs_types_of_the_new_columns() -> None:
    from gippa.novel.outputs import NOVEL_ORFS_TYPES

    want = {c: "int" if c in DISTANCES else "str" for c in NEW_COLUMNS}
    assert {c: NOVEL_ORFS_TYPES.get(c) for c in NEW_COLUMNS} == want


def test_schema_lists_the_columns_in_order() -> None:
    assert tuple(_fields()) == COLUMNS


@pytest.mark.parametrize("col", NEW_COLUMNS)
def test_schema_type_of_each_new_column(col: str) -> None:
    f = _fields()[col]
    assert f["type"] == TYPES[col]
    assert f.get("description"), col


def test_schema_overlap_relation_enum_and_distance_minimum() -> None:
    f = _fields()
    assert f["overlap_relation"]["constraints"]["enum"] == ["same_strand", "antisense", "none"]
    for c in DISTANCES:
        assert f[c]["constraints"]["minimum"] == 0, c


def test_schema_flank_descriptions_state_the_strand_convention() -> None:
    f = _fields()
    for c in ("upstream_gene", "downstream_gene"):
        assert "strand" in f[c]["description"].lower(), (c, f[c]["description"])


def test_output_schema_version_stays_1_3() -> None:
    """Q3 default: 1.3 is unreleased, so the additions do not bump it."""
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    assert OUTPUT_SCHEMA_VERSION == "1.3"
    pkg = json.loads((SCHEMA.parents[1] / "datapackage.json").read_text(encoding="utf-8"))
    assert pkg["output_schema_version"] == "1.3"


# ---------------------------------------------------------------- AC7: older trees


def test_migration_added_columns() -> None:
    from gippa.io.migration import ADDED_COLUMNS

    assert {t: tuple(c) for t, c in ADDED_COLUMNS.items()} == {"novel_orfs": NEW_COLUMNS}


ROWS = [
    {"orf_id": "NOVEL_chr1_101_p", "class": "intergenic", "contig": "chr1", "strand": "+",
     "frame": 3, "start_longest": 101, "end": 400, "length_aa_longest": 99, "length_aa_best": 99,
     "n_peptides": 2, "class_fdr": "0.004", "called": True, "overlapping_gene": None},
    {"orf_id": "NOVEL_chr1_1401_m", "class": "antisense", "contig": "chr1", "strand": "-",
     "frame": -2, "start_longest": 1700, "end": 1401, "length_aa_longest": 99,
     "length_aa_best": 99, "n_peptides": 1, "class_fdr": "0.02", "called": False,
     "overlapping_gene": "g3"},
]  # fmt: skip


def _old_tree(tmp_path: Path, header: list[str]) -> Path:
    """An OUT/ whose novel_orfs.tsv has ``header`` (old names allowed) and the ``ROWS`` values."""
    from gippa.io.tsv import write_tsv

    out = run_out(tmp_path)
    renamed = dict(zip(to_1_0("novel_orfs", list(RC1_COLUMNS)), RC1_COLUMNS, strict=True))
    write_tsv(out / "tables" / "novel_orfs.tsv", header,
              [[r.get(renamed.get(c, c)) for c in header] for r in ROWS])  # fmt: skip
    return out


def _read(out: Path) -> Any:
    from gippa.io.run_view import RunView

    return RunView.from_outdir(out).table("novel_orfs")


@pytest.mark.parametrize("release", ["1.1.0rc1", "1.0"])
def test_runview_reads_an_old_header_with_the_new_columns_na(tmp_path: Path, release: str) -> None:
    header = list(RC1_COLUMNS) if release == "1.1.0rc1" else to_1_0("novel_orfs", list(RC1_COLUMNS))
    got = _read(_old_tree(tmp_path, header))
    assert list(got.columns) == list(COLUMNS)
    for c in NEW_COLUMNS:
        assert got[c].isna().all(), c  # overlap_relation too: not computed, so NA (rule 8)
    assert list(got["orf_id"]) == [r["orf_id"] for r in ROWS]
    assert list(got["overlapping_gene"].fillna("NA")) == ["NA", "g3"]


def test_runview_still_raises_on_another_missing_column(tmp_path: Path) -> None:
    header = [c for c in RC1_COLUMNS if c != "class"]
    with pytest.raises(ValueError, match="class"):
        _read(_old_tree(tmp_path, header))


@settings(max_examples=30, deadline=None, suppress_health_check=[
    HealthCheck.too_slow, HealthCheck.function_scoped_fixture])  # fmt: skip
@given(st.sets(st.sampled_from(NEW_COLUMNS)))
def test_runview_fills_exactly_the_absent_new_columns(tmp_path: Path, dropped: set[str]) -> None:
    """Some new columns present (values kept), the others NA; the schema order either way."""
    import tempfile

    values = {"overlapping_gene_label": "lacZ", "overlap_relation": "antisense",
              "other_overlapping_genes": "g2", "upstream_gene": "g4",
              "upstream_gene_label": "g4", "upstream_distance_nt": 7,
              "downstream_gene": "g1", "downstream_gene_label": "dnaA",
              "downstream_distance_nt": 0}  # fmt: skip
    with tempfile.TemporaryDirectory(dir=tmp_path) as tmp:
        out = run_out(Path(tmp))
        rows = [ROWS[1] | values]
        write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", rows, drop=dropped)
        got = _read(out)
    assert list(got.columns) == list(COLUMNS)
    for c in NEW_COLUMNS:
        v = got[c].iloc[0]
        if c in dropped:
            assert pd.isna(v), (c, v)
        else:
            assert (int(v) if c in DISTANCES else v) == values[c], (c, v)


def test_read_novel_orfs_keeps_its_exact_header_check(tmp_path: Path) -> None:
    """Rule 8: the internal reader does not fill; an rc1 header is refused."""
    from gippa.io.novel_tables import read_novel_orfs

    out = _old_tree(tmp_path, list(RC1_COLUMNS))
    with pytest.raises(ValueError, match="header"):
        read_novel_orfs(out / "tables" / "novel_orfs.tsv")


@pytest.mark.skipif(not (REPO / ROUTE_E).is_file() and not (REPO.parents[2] / ROUTE_E).is_file(),
                    reason="route E 1.1.0rc1 copy is local only (not in CI)")  # fmt: skip
def test_runview_reads_the_route_e_rc1_table(tmp_path: Path) -> None:
    src = REPO / ROUTE_E if (REPO / ROUTE_E).is_file() else REPO.parents[2] / ROUTE_E
    header, rows = read_tsv(src)
    assert tuple(header) == RC1_COLUMNS
    out = run_out(tmp_path)
    (out / "tables" / "novel_orfs.tsv").write_bytes(src.read_bytes())
    got = _read(out)
    assert len(got) == len(rows) and list(got.columns) == list(COLUMNS)
    assert all(got[c].isna().all() for c in NEW_COLUMNS)


# ---------------------------------------------------------------- AC8: Novel ORFs preview


def _page(out: Path) -> Page:
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    return parse(build_report(RunView.from_outdir(out), {}, {}, cfg()))


PREVIEW_ROWS = [
    ROWS[1] | {"rank": 2, "overlapping_gene_label": "lacZ", "overlap_relation": "antisense",
               "upstream_gene": "g4", "upstream_gene_label": "g4", "upstream_distance_nt": 7800,
               "downstream_gene": "g2", "downstream_gene_label": "g2",
               "downstream_distance_nt": 300, "flags": "rbs_supported"},
    ROWS[0] | {"rank": 1, "overlap_relation": "none", "upstream_gene": "g3",
               "upstream_gene_label": "lacZ", "upstream_distance_nt": 100},
]  # fmt: skip
#: the cells as the preview shows them (NA labels: none / —; everything else as written)
PREVIEW_WANT = [
    ["NOVEL_chr1_1401_m", "antisense", "false", "0.02", "-", "99", "lacZ", "antisense", "g4",
     "7800", "g2", "300", "1", "rbs_supported"],
    ["NOVEL_chr1_101_p", "intergenic", "true", "0.004", "+", "99", "none", "none", "lacZ", "100",
     "—", "NA", "2", "NA"],
]  # fmt: skip


def _novel_preview(tmp_path: Path, rows: list[dict[str, Any]]) -> tuple[list[str], list[list[str]]]:
    out = run_out(tmp_path)
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", rows)
    page = _page(out)
    table = page.table("novel_orfs")
    assert table is not None, "no novel_orfs preview"
    sec = [s for s in page.sections() if s.attrs.get("id") == "novel-orfs"]
    assert sec and table in list(sec[0].iter()), "the preview is not in the Novel ORFs section"
    return table_rows(table)


def test_novel_preview_constant() -> None:
    from gippa.report.build import NOVEL_PREVIEW as got

    assert tuple(got) == NOVEL_PREVIEW


def test_novel_preview_header_and_cells(tmp_path: Path) -> None:
    header, body = _novel_preview(tmp_path, PREVIEW_ROWS)
    assert header == list(NOVEL_PREVIEW)
    assert body == PREVIEW_WANT  # table order (rank 2 first), not re-sorted


def test_novel_preview_keeps_table_order_and_top_rows(tmp_path: Path) -> None:
    from gippa.report.build import TOP_ROWS

    rows = [ROWS[0] | {"orf_id": f"NOVEL_chr1_{1 + 9000 - 300 * i}_p", "rank": i + 1}
            for i in range(TOP_ROWS + 5)]  # fmt: skip
    _, body = _novel_preview(tmp_path, rows)
    assert [r[0] for r in body] == [r["orf_id"] for r in rows[:TOP_ROWS]]


def test_novel_tsv_keeps_na(tmp_path: Path) -> None:
    """Rule 6: the table writes NA; only the preview shows none / —."""
    out = run_out(tmp_path)
    path = write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", PREVIEW_ROWS)
    header, rows = read_tsv(path)
    assert header == list(COLUMNS)
    i, j = header.index("overlapping_gene_label"), header.index("downstream_gene_label")
    assert rows[1][i] is None and rows[1][j] is None


# ---------------------------------------------------------------- genes preview (lead decision)


GENES = [
    {"locus_tag": "g1", "gene": "dnaA", "product": "chromosomal replication initiator protein"},
    {"locus_tag": "g2", "gene": None, "product": "hypothetical protein"},
    {"locus_tag": "g3", "gene": None, "product": 'ABC transporter "permease" & <ATP>-binding'},
    {"locus_tag": "g4", "gene": "  ", "product": "DUF1234 domain-containing protein"},
    {"locus_tag": "g0", "gene": "abcB", "product": "ABC transporter"},
]


def _genes_page(tmp_path: Path) -> tuple[Path, Page]:
    out = run_out(tmp_path)
    body = [gene_row(locus_tag=g["locus_tag"], gene=g["gene"], product=g["product"],
                     start=101 + 1000 * i, end=400 + 1000 * i, n_peptides_unique=5 - i,
                     abundance_log2=20.0 - i, detection_tier="two_peptide",
                     abundance_type="ibaq", abundance_scale="molar")
            for i, g in enumerate(GENES)]  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", body)
    return out, _page(out)


def _gene_cells(page: Page) -> list[tuple[str, Node]]:
    """(locus_tag text, ``gene`` <td>) per preview row."""
    table = page.table("genes_tsv")
    assert table is not None, "no genes preview"
    header, _ = table_rows(table)
    it, ig = header.index("locus_tag"), header.index("gene")
    out = []
    for tr in table.by_tag("tr"):
        tds = [c for c in tr.children if isinstance(c, Node) and c.tag == "td"]
        if tds:
            out.append((tds[it].text().strip(), tds[ig]))
    return out


def _title(td: Node) -> str | None:
    ems = td.by_tag("em")
    return (ems[0].attrs.get("title") if ems else None) or td.attrs.get("title")


def test_genes_preview_gene_cell_never_reads_na(tmp_path: Path) -> None:
    _, page = _genes_page(tmp_path)
    cells = _gene_cells(page)
    assert len(cells) == len(GENES)
    assert all(td.text().strip() not in ("NA", "") for _, td in cells), [
        td.text() for _, td in cells]  # fmt: skip


def test_genes_preview_unnamed_gene_is_the_locus_tag_in_italics(tmp_path: Path) -> None:
    _, page = _genes_page(tmp_path)
    by_tag = dict(_gene_cells(page))
    product = {g["locus_tag"]: g["product"] for g in GENES}
    for tag in ("g2", "g3", "g4"):
        td = by_tag[tag]
        ems = td.by_tag("em")
        assert len(ems) == 1 and ems[0].text() == tag, (tag, td.text())
        assert td.text().strip() == tag
        title = _title(td)
        assert title is not None and product[tag] in title, (tag, title)


def test_genes_preview_named_gene_is_plain(tmp_path: Path) -> None:
    _, page = _genes_page(tmp_path)
    by_tag = dict(_gene_cells(page))
    for tag, name in (("g1", "dnaA"), ("g0", "abcB")):
        assert by_tag[tag].text().strip() == name
        assert not by_tag[tag].by_tag("em"), tag


def test_genes_preview_markup_is_escaped(tmp_path: Path) -> None:
    """The product goes into an attribute: quotes, ``&`` and ``<`` must not break the markup."""
    _, page = _genes_page(tmp_path)
    assert "<ATP>" not in page.html
    by_tag = dict(_gene_cells(page))
    assert _title(by_tag["g3"]) is not None
    assert GENES[2]["product"] in str(_title(by_tag["g3"]))


def test_genes_preview_sorts_by_the_shown_label(tmp_path: Path) -> None:
    """``tables.sort_order`` over the shown text: a label cell is never in the NA bucket."""
    from gippa.report.tables import sort_order

    _, page = _genes_page(tmp_path)
    shown = [td.text().strip() for _, td in _gene_cells(page)]
    assert sort_order(shown) == sorted(range(len(shown)), key=lambda i: shown[i])


def test_genes_tsv_keeps_na(tmp_path: Path) -> None:
    out, _ = _genes_page(tmp_path)
    header, rows = read_tsv(out / "tables" / "genes.tsv")
    i = header.index("gene")
    assert [r[i] for r in rows] == ["dnaA", None, None, "  ", "abcB"]


# ---------------------------------------------------------------- AC5, AC6: written runs


def _written() -> list[Path]:
    return [*sorted(SNAPSHOT.glob("*/OUT/tables/novel_orfs.tsv")), GOLDEN]


@pytest.mark.parametrize("path", _written(), ids=lambda p: str(p.relative_to(REPO)))
def test_written_novel_orfs_has_the_new_header(path: Path) -> None:
    """Edge case "empty table" included (sim_empty): the full header, no rows."""
    header, _ = read_tsv(path)
    assert tuple(header) == COLUMNS


@pytest.mark.parametrize("path", _written(), ids=lambda p: str(p.relative_to(REPO)))
def test_ac5_on_written_runs(path: Path) -> None:
    header, rows = read_tsv(path)
    assert set(NEW_COLUMNS) <= set(header), "not regenerated"
    ix = {c: header.index(c) for c in header}
    for r in rows:
        ovl, lab, rel = r[ix["overlapping_gene"]], r[ix["overlapping_gene_label"]], (
            r[ix["overlap_relation"]])  # fmt: skip
        assert (ovl is None) == (lab is None), r[0]
        assert (ovl is None) == (rel == "none"), r[0]
        assert rel in ("same_strand", "antisense", "none"), (r[0], rel)
        for side in ("upstream", "downstream"):
            cells = [r[ix[f"{side}_{c}"]] for c in ("gene", "gene_label", "distance_nt")]
            assert len({c is None for c in cells}) == 1, (r[0], side, cells)
            assert cells[2] is None or (re.fullmatch(r"[0-9]+", cells[2]) is not None), cells


def test_ac6_tags_join_genes_tsv_in_the_snapshot_run() -> None:
    """Every tag is a genes.tsv locus tag of the same run; each label follows rule 5 there."""
    out = SNAPSHOT / "sim_plots" / "OUT" / "tables"
    gh, grows = read_tsv(out / "genes.tsv")
    names = {r[gh.index("locus_tag")]: r[gh.index("gene")] for r in grows}
    header, rows = read_tsv(out / "novel_orfs.tsv")
    assert set(NEW_COLUMNS) <= set(header), "not regenerated"
    ix = {c: header.index(c) for c in header}
    assert rows, "the snapshot run has novel ORFs"
    seen = 0
    for r in rows:
        tags = [r[ix[c]] for c in TAG_COLUMNS]
        tags += (r[ix["other_overlapping_genes"]] or "").split(";")
        for t in filter(None, tags):
            assert t in names, (r[0], t)
            seen += 1
        for c in TAG_COLUMNS:
            t = r[ix[c]]
            if t is None:
                continue
            name = names[t]
            usable = name is not None and name.strip() and name.strip().upper() != "NA"
            assert r[ix[f"{c}_label"]] == (name if usable else t), (r[0], c)
    assert seen, "no tag checked"


# ---------------------------------------------------------------- AC9: reference page


def test_reference_tables_md_lists_the_new_rows_in_order() -> None:
    text = (REPO / "docs" / "reference" / "tables.md").read_text(encoding="utf-8")
    sec = text.split("## `novel_orfs`", 1)[1].split("\n## ", 1)[0]
    fields = re.findall(r"^\| `([a-z0-9_]+)` \|", sec, flags=re.M)
    assert tuple(fields) == COLUMNS
