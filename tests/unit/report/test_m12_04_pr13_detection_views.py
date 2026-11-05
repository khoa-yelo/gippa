# ruff: noqa: RUF001, RUF003  (the report's wording is Unicode: −, ≤)
"""M12-04 PR 13, T-07 and AC 13.1: the Gene detection section shows two pivoted views of
``tables/detection_summary.tsv`` instead of its raw 20-row preview
(docs/specs/M12-04-outputs-revision.md § PR 13).

- **Tier view:** one row per distinct ``q_threshold`` group of the ``tier`` section, ascending,
  with the columns ``n_detected``, ``n_genes``, ``fdr``, ``fdr_2pep`` and ``fp`` (the spec
  example's header: q, detected, genes, FDR, FDR (2-pep), FP). Each cell is the long file's value
  through ``gippa.viz.style.fmt_num``.
- **Detectability view:** O/E by group: one row per ``(group_type, group)`` of the rows whose
  ``metric`` is ``o_over_e`` (section ``basic``, ``gippa.context.summary.basic_rows``).
- Tiers ``protein_list`` and ``two_peptide`` have no ``q_threshold`` rows: the tier view shows
  the tier's counts only (no FDR or FP column).
- ``detection_summary`` absent: both views are omitted. No new table file is written; the long
  file is unchanged.

Review r1 (lead decisions, binding):
- L1: the tier view drops a column that is NA on every row (``decoy_picked`` writes
  ``n_detected`` and ``fdr`` only: ``q | detected | FDR``);
- L2: ``protein_list`` shows ``tier | detected`` only; ``two_peptide`` adds
  ``detected with ≥ 2 peptides`` (``n_detected_2pep``);
- L3: the section's ``<a href="../tables/detection_summary.tsv">`` resolves from ``OUT/report/``;
- L4: only ``section = tier, group_type = q_threshold`` rows make tier rows and only
  ``section = basic, metric = o_over_e`` rows make O/E rows (distractors in ``other_rows``);
- L5: a non-numeric or non-finite ``q_threshold`` group (``nan``, ``inf``, ``-inf``) sorts after
  every finite one and leaves the finite ones ascending;
- M3: an unreadable long file (garbage bytes, a wrong header, ``value = abc``) builds the page with
  no view and a ``<p class="note">`` naming ``detection_summary.tsv``.

Interface pinned here (PR 13 test author, beyond the spec text):
- the views are ``<table data-view="detection-tier">`` and ``<table data-view="detection-oe">``
  inside ``<section id="detection">``; neither carries ``data-table`` (they are no table file),
  and ``<table data-table="detection_summary">`` (the raw preview) is gone;
- the tier view's first row is a header row of ``<th>`` cells; the first column holds the
  ``q_threshold`` group as written in the file (``0.001``, ``0.005``, ``0.01``; equal to
  ``fmt_alpha``); a body row may use ``<th>`` or ``<td>`` for it;
- counts (``n_detected``, ``n_genes``) are exact: ``fmt_num(int(value))``, so 158 genes reads
  "158", never "160" (``fmt_num`` keeps integers exact; the file's numbers read back as floats);
  every other value is ``fmt_num(value)`` with the default 2 significant digits; NA -> "NA".
"""

from __future__ import annotations

import math
import random
import tempfile
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pytest
from helpers.report_fixtures import Node, Page, cfg, parse, write_preflight
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import Run, build

TIER_VIEW = "detection-tier"
OE_VIEW = "detection-oe"
#: the spec example's header, in order (§ PR 13 Examples)
TIER_HEADER = ["q", "detected", "genes", "FDR", "FDR (2-pep)", "FP"]
TIER_METRICS = ("n_detected", "n_genes", "fdr", "fdr_2pep", "fp")
COUNTS = frozenset({"n_detected", "n_genes"})
NA = float("nan")

#: fixture P's mayu curve: the spec example's first two rows, then made-up rows (n_genes 158
#: so a count of 3 digits is pinned exact; one NA)
P_CURVE: Mapping[str, Mapping[str, float]] = {
    "0.001": {"n_detected": 94, "n_genes": 94, "fdr": 0.00065, "fdr_2pep": 5.9e-07, "fp": 0.061},
    "0.005": {"n_detected": 98, "n_genes": 99, "fdr": 0.0061, "fdr_2pep": 6.1e-05, "fp": 0.60},
    "0.01": {"n_detected": 99, "n_genes": 101, "fdr": 0.0124093, "fdr_2pep": 0.00027,
             "fp": 1.2533},
    "0.02": {"n_detected": 103, "n_genes": 112, "fdr": 0.0213, "fdr_2pep": 0.00081, "fp": 2.3856},
    "0.05": {"n_detected": 118, "n_genes": 131, "fdr": 0.0467, "fdr_2pep": 0.0031, "fp": 6.1177},
    "0.1": {"n_detected": 140, "n_genes": 158, "fdr": 0.0921, "fdr_2pep": NA, "fp": 14.552},
}  # fmt: skip
#: (group_type, group) -> O/E of fixture P's ``basic`` rows (E = 0 gives NA)
P_OE: Mapping[tuple[str, str], float] = {
    ("replicon", "chr1"): 1.21771,
    ("feature_class", "cds"): 1.25006,
    ("feature_class", "ncRNA"): NA,
    ("cog_category", "J"): 1.41232,
}


# ---------------------------------------------------------------- the long file


def _row(section: str, gtype: str, group: str, metric: str, value: float | None,
         **kw: Any) -> dict[str, Any]:  # fmt: skip
    v = None if value is None or (isinstance(value, float) and math.isnan(value)) else value
    return {"section": section, "group_type": gtype, "group": group, "metric": metric,
            "value": v, **kw}  # fmt: skip


def tier_rows(tier: str, n_detected: int, n_2pep: int) -> list[dict[str, Any]]:
    """The ``group_type = tier`` rows every tier writes (``gippa.genetable.detection.detect``)."""
    return [_row("tier", "tier", tier, "n_detected", n_detected),
            _row("tier", "tier", tier, "n_detected_2pep", n_2pep),
            _row("tier", "tier", tier, "n_decoys_paired", 0),
            _row("tier", "tier", tier, "n_decoys_unpaired", 0),
            _row("tier", "tier", tier, "decoy_paired_share", None),
            _row("tier", "tier", tier, "decoy_win_weight", None),
            _row("tier", "tier", tier, "q_floor", None)]  # fmt: skip


def curve_rows(curve: Mapping[str, Mapping[str, float]]) -> list[dict[str, Any]]:
    return [_row("tier", "q_threshold", q, m, v) for q, ms in curve.items() for m, v in ms.items()]


def basic_rows(oe: Mapping[tuple[str, str], float]) -> list[dict[str, Any]]:
    rows = []
    for (gtype, group), ratio in oe.items():
        e = 0.0 if math.isnan(ratio) else 10.0
        o = 0.0 if math.isnan(ratio) else round(ratio * e, 4)
        rows += [_row("basic", gtype, group, "n_genes", 20, n=20),
                 _row("basic", gtype, group, "n_detected", o, n=20),
                 _row("basic", gtype, group, "expected", e, n=20),
                 _row("basic", gtype, group, "o_minus_e", o - e, ci_low=-1.0, ci_high=3.0, n=20),
                 _row("basic", gtype, group, "o_over_e", ratio, n=20),
                 _row("basic", gtype, group, "n_tu", 12, n=20)]  # fmt: skip
    return rows


#: distractor groups (upper case and ``_``: the property's drawn groups are ``[a-z]{3,8}``)
OTHER_SECTION_OE = "X_other_section"
LONELY_GROUP = "X_no_oe"
BASIC_Q = "0.5"


def other_rows() -> list[dict[str, Any]]:
    """``detectability`` and ``qc`` rows and the review r1 L4 distractors: neither view shows
    them (an ``o_over_e`` row outside ``basic``, a ``q_threshold`` row in ``basic``, a ``basic``
    group without ``o_over_e``)."""
    rows = [
        _row("detectability", "bin", str(b), m, v)
        for b in range(3)
        for m, v in (("mean_p", 0.1 * b), ("observed_rate", 0.1 * b + 0.01), ("n", 15))
    ]
    rows += [_row("detectability", "overall", "all", "calibration_slope", 0.97),
             _row("detectability", "overall", "all", "brier", 0.18),
             _row("qc", "run", "run", "unmapped_peptides", 3),
             _row("qc", "unmapped_reason", "no_hit", "n_peptides", 3),
             _row("detectability", "replicon", OTHER_SECTION_OE, "o_over_e", 9.75),
             _row("basic", "q_threshold", BASIC_Q, "n_detected", 777),
             _row("basic", "q_threshold", BASIC_Q, "fdr", 0.333),
             _row("basic", "replicon", LONELY_GROUP, "n_genes", 5, n=5),
             _row("basic", "replicon", LONELY_GROUP, "n_detected", 2, n=5),
             _row("basic", "replicon", LONELY_GROUP, "expected", 3.0, n=5)]  # fmt: skip
    return rows


DISTRACTORS = frozenset({OTHER_SECTION_OE, LONELY_GROUP, BASIC_Q, "777", "9.8", "0.33"})


def write_summary(out: Path, rows: Sequence[Mapping[str, Any]]) -> Path:
    return write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", rows)


def fixture_p(tmp: Path) -> Run:
    """Fixture P (mayu) with a full ``detection_summary``: tier rows, the curve, basic O/E rows,
    detectability and qc rows (in file order tier, basic, detectability, qc)."""
    run = build("P", tmp)
    write_summary(run.out, [*tier_rows("mayu", 99, 91), *curve_rows(P_CURVE),
                            *basic_rows(P_OE), *other_rows()])  # fmt: skip
    return run


# ---------------------------------------------------------------- reading the page


def grid(table: Node) -> tuple[list[str], list[list[str]]]:
    """Header cells and body rows of ``table``; ``<th>`` and ``<td>`` cells both count, in
    order (a body row may head itself with a ``<th>``)."""
    header: list[str] = []
    body: list[list[str]] = []
    for tr in table.by_tag("tr"):
        cells = [c for c in tr.children if isinstance(c, Node) and c.tag in ("th", "td")]
        texts = [" ".join(c.text().split()) for c in cells]
        if not header and not body and cells and all(c.tag == "th" for c in cells):
            header = texts
        elif cells:
            body.append(texts)
    return header, body


def views(page: Page, name: str) -> list[Node]:
    return [t for t in page.root.by_tag("table") if t.attrs.get("data-view") == name]


def view(page: Page, name: str) -> Node:
    found = views(page, name)
    assert len(found) == 1, f'{len(found)} <table data-view="{name}">'
    det = page.section("detection")
    assert found[0] in list(det.iter()), f"{name} is not in the Gene detection section"
    assert "data-table" not in found[0].attrs, "a view is no table file (no data-table)"
    return found[0]


def summary_values(out: Path) -> Any:
    """``detection_summary`` as the report reads it (values as floats after ``%.6g``)."""
    from gippa.io.run_view import RunView

    frame = RunView.from_outdir(out).table("detection_summary")
    assert frame is not None
    return frame


def shown(metric: str, value: Any) -> str:
    """The cell text of ``value`` (the oracle; see the module docstring)."""
    from gippa.viz.style import fmt_num

    v = float(value) if value is not None else NA
    if math.isnan(v):
        return "NA"
    return fmt_num(int(v)) if metric in COUNTS and v.is_integer() else fmt_num(v)


def expected_grid(frame: Any) -> tuple[list[str], list[list[str]]]:
    """The tier view from the long file: ascending ``q_threshold``, ``TIER_METRICS``, a column
    that is NA on every row dropped (review r1 L1)."""
    curve = frame[(frame["section"].astype(str) == "tier")
                  & (frame["group_type"].astype(str) == "q_threshold")]  # fmt: skip
    groups = sorted(set(curve["group"].astype(str)), key=float)
    full = []
    for g in groups:
        sub = curve[curve["group"].astype(str) == g]
        by = dict(zip(sub["metric"].astype(str), sub["value"], strict=True))
        full.append([g, *(shown(m, by.get(m)) for m in TIER_METRICS)])
    keep = [0] + [k for k in range(1, len(TIER_HEADER)) if any(r[k] != "NA" for r in full)]
    return [TIER_HEADER[k] for k in keep], [[r[k] for k in keep] for r in full]


def expected_tier(frame: Any) -> list[list[str]]:
    return expected_grid(frame)[1]


def expected_oe(frame: Any) -> list[tuple[str, str]]:
    """``(group, O/E text)`` of every ``o_over_e`` row of the ``basic`` section."""
    oe = frame[(frame["metric"].astype(str) == "o_over_e")
               & (frame["section"].astype(str) == "basic")]  # fmt: skip
    return [(str(g), shown("o_over_e", v)) for g, v in zip(oe["group"], oe["value"], strict=True)]


def column(header: list[str], rows: list[list[str]], label: str) -> list[str]:
    assert label in header, f"no {label!r} column in {header}"
    k = header.index(label)
    return [r[k] for r in rows]


# ---------------------------------------------------------------- AC 13.1: the tier view


def test_ac_13_1_tier_view_rows_and_values(tmp_path: Path) -> None:
    """On fixture P: one row per distinct ``q_threshold`` of the ``tier`` section, ascending;
    the values equal the long file's through ``fmt_num``."""
    run = fixture_p(tmp_path)
    header, rows = grid(view(run.page(), TIER_VIEW))
    assert header == TIER_HEADER
    assert rows == expected_tier(summary_values(run.out))
    assert [r[0] for r in rows] == list(P_CURVE)  # 0.001, 0.005, 0.01, 0.02, 0.05, 0.1


def test_tier_view_spec_example_rows(tmp_path: Path) -> None:
    """The spec example's first two rows, literally."""
    run = fixture_p(tmp_path)
    _, rows = grid(view(run.page(), TIER_VIEW))
    assert rows[0] == ["0.001", "94", "94", "0.00065", "5.9e−7", "0.061"]
    assert rows[1] == ["0.005", "98", "99", "0.0061", "6.1e−5", "0.60"]


def test_tier_view_counts_are_exact_and_na_is_na(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    header, rows = grid(view(run.page(), TIER_VIEW))
    last = dict(zip(header, rows[-1], strict=True))
    assert last["genes"] == "158" and last["detected"] == "140"
    assert last["FDR (2-pep)"] == "NA"


def test_tier_view_leaves_out_the_tier_rows(tmp_path: Path) -> None:
    """The ``group_type = tier`` rows (``n_detected_2pep``, ``q_floor``, …) are not curve rows."""
    run = fixture_p(tmp_path)
    _, rows = grid(view(run.page(), TIER_VIEW))
    assert len(rows) == len(P_CURVE)
    assert all(r[0] != "mayu" for r in rows)


def test_decoy_picked_tier_view(tmp_path: Path) -> None:
    """``decoy_picked`` writes ``n_detected`` and ``fdr`` per threshold: one row each."""
    run = build("decoy", tmp_path)
    curve = {
        q: {"n_detected": 90 + k, "fdr": 0.0005 * (k + 1)}
        for k, q in enumerate(("0.001", "0.005", "0.01", "0.02", "0.05", "0.1"))
    }
    write_summary(run.out, [*tier_rows("decoy_picked", 92, 80), *curve_rows(curve)])
    header, rows = grid(view(run.page(), TIER_VIEW))
    from gippa.viz.style import fmt_num

    # review r1 L1: genes, FDR (2-pep) and FP are NA on every row, so they are dropped
    assert header == ["q", "detected", "FDR"]
    assert rows == [[q, str(90 + k), fmt_num(0.0005 * (k + 1))] for k, q in enumerate(curve)]


# ---------------------------------------------------------------- the detectability view


def test_oe_view_one_row_per_group(tmp_path: Path) -> None:
    """O/E by group: one row per ``o_over_e`` row; each row shows its group and its O/E through
    ``fmt_num`` (NA when E = 0); detectability and qc rows are not in it."""
    run = fixture_p(tmp_path)
    _, rows = grid(view(run.page(), OE_VIEW))
    want = expected_oe(summary_values(run.out))
    assert len(rows) == len(want) == len(P_OE)
    for group, text in want:
        hits = [r for r in rows if group in r]
        assert len(hits) == 1, f"group {group!r}: {len(hits)} rows in {rows}"
        assert text in hits[0], f"{group}: O/E {text!r} not in {hits[0]}"
    flat = {c for r in rows for c in r}
    assert not {"no_hit", "unmapped_reason", "calibration_slope", "brier"} & flat
    assert not DISTRACTORS & flat, DISTRACTORS & flat


def test_views_ignore_the_distractors(tmp_path: Path) -> None:
    """Review r1 L4: an ``o_over_e`` row outside ``basic``, a ``q_threshold`` row in ``basic``
    and a ``basic`` group without ``o_over_e`` are in no view."""
    page = fixture_p(tmp_path).page()
    for name in (TIER_VIEW, OE_VIEW):
        _, rows = grid(view(page, name))
        flat = {c for r in rows for c in r}
        assert not DISTRACTORS & flat, (name, DISTRACTORS & flat)
    _, tier = grid(view(page, TIER_VIEW))
    assert [r[0] for r in tier] == list(P_CURVE)


def test_full_table_link_resolves(tmp_path: Path) -> None:
    """Review r1 L3: the section links the long file; the link resolves from ``OUT/report/``."""
    run = fixture_p(tmp_path)
    links = [a.attrs["href"] for a in run.page().section("detection").by_tag("a")
             if a.attrs.get("href", "").endswith("detection_summary.tsv")]  # fmt: skip
    assert links == ["../tables/detection_summary.tsv"], links
    target = (run.out / "report" / links[0]).resolve()
    assert target == (run.out / "tables" / "detection_summary.tsv").resolve()
    assert target.is_file()


def test_oe_view_values(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    _, rows = grid(view(run.page(), OE_VIEW))
    by_group = {g: r for r in rows for g in ("chr1", "cds", "ncRNA", "J") if g in r}
    assert "1.2" in by_group["chr1"] and "1.3" in by_group["cds"]
    assert "NA" in by_group["ncRNA"] and "1.4" in by_group["J"]


# ---------------------------------------------------------------- the raw preview goes


def test_raw_detection_summary_preview_is_removed(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    assert page.table("detection_summary") is None, "the raw 20-row preview is still shown"
    heads = [" ".join(h.text().split()) for h in page.section("detection").by_tag("h3")]
    assert not [h for h in heads if h.startswith("detection_summary.tsv (first")], heads
    assert page.table("genes_tsv") is not None, "the genes preview stays"


def test_no_new_table_file_and_long_file_unchanged(tmp_path: Path) -> None:
    """Writing the report adds no file outside ``report/`` (Q8) and leaves the long file as it
    was."""
    from gippa.steps.plots import write_report

    run = fixture_p(tmp_path)
    long_file = run.out / "tables" / "detection_summary.tsv"
    before = {p.relative_to(run.out) for p in run.out.rglob("*")}
    body = long_file.read_bytes()
    write_report(run.view(), None, run.out / "report", run.cfg())
    after = {p.relative_to(run.out) for p in run.out.rglob("*")}
    assert {p for p in after - before if p.parts[0] != "report"} == set()
    assert long_file.read_bytes() == body
    page = parse((run.out / "report" / "report.html").read_text(encoding="utf-8"))
    assert len(views(page, TIER_VIEW)) == 1


# ---------------------------------------------------------------- edge cases


def test_protein_list_view_is_tier_and_detected(tmp_path: Path) -> None:
    """Review r1 L2: ``protein_list`` has no peptide counts: ``tier | detected`` only."""
    run = build("R", tmp_path)
    write_summary(run.out, tier_rows("protein_list", 106, 0))
    header, rows = grid(view(run.page(), TIER_VIEW))
    assert header == ["tier", "detected"]
    assert rows == [["protein_list", "106"]]


def test_two_peptide_view_keeps_the_2pep_count(tmp_path: Path) -> None:
    """Review r1 L2: ``two_peptide`` keeps ``n_detected_2pep`` (values differ on purpose, so a
    swapped column shows)."""
    run = build("D", tmp_path)
    write_summary(run.out, tier_rows("two_peptide", 5, 4))
    header, rows = grid(view(run.page(), TIER_VIEW))
    assert header == ["tier", "detected", "detected with ≥ 2 peptides"]
    assert rows == [["two_peptide", "5", "4"]]


@pytest.mark.parametrize(("name", "tier", "n_detected", "n_2pep"),
                         [("R", "protein_list", 106, 0), ("D", "two_peptide", 4, 4)])  # fmt: skip
def test_count_only_tiers_show_counts_only(
    tmp_path: Path, name: str, tier: str, n_detected: int, n_2pep: int
) -> None:
    """``protein_list`` and ``two_peptide`` have no FDR rows: the tier view shows the counts
    only (no FDR, FDR (2-pep) or FP column) and holds the tier's ``n_detected``."""
    run = build(name, tmp_path)
    write_summary(run.out, tier_rows(tier, n_detected, n_2pep))
    header, rows = grid(view(run.page(), TIER_VIEW))
    assert rows, "the count-only tier view has no row"
    assert not [h for h in header if "FDR" in h or h == "FP"], header
    assert str(n_detected) in {c for r in rows for c in r}, rows


@pytest.mark.parametrize("name", ["R", "P"])
def test_absent_detection_summary_omits_both_views(tmp_path: Path, name: str) -> None:
    run = build(name, tmp_path)
    (run.out / "tables" / "detection_summary.tsv").unlink(missing_ok=True)
    page = run.page()
    assert not views(page, TIER_VIEW) and not views(page, OE_VIEW)
    assert page.section("detection") is not None


def test_header_only_detection_summary_gives_no_tier_row(tmp_path: Path) -> None:
    """A header-only long file: no curve rows, so no tier view body row is invented."""
    run = build("P", tmp_path)
    write_summary(run.out, [])
    page = run.page()  # must not raise
    for t in views(page, TIER_VIEW) + views(page, OE_VIEW):
        _, rows = grid(t)
        assert not [r for r in rows if r and r[0] in P_CURVE], rows


def _break(path: Path, how: str) -> None:
    if how == "garbage":
        path.write_bytes(b"\x00\xff\xfe not a table \x01\x02\n\x9c\x80")
        return
    lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
    if how == "header":
        lines[0] = lines[0].replace("metric", "zzmetric", 1)
    else:  # value = abc on every body row
        k = lines[0].rstrip("\n").split("\t").index("value")
        for i in range(1, len(lines)):
            cells = lines[i].rstrip("\n").split("\t")
            cells[k] = "abc"
            lines[i] = "\t".join(cells) + "\n"
    path.write_text("".join(lines), encoding="utf-8")


@pytest.mark.parametrize("how", ["garbage", "header", "value_abc"])
def test_unreadable_detection_summary_is_a_note(tmp_path: Path, how: str) -> None:
    """Review r1 M3: the page builds, no view is shown and a note names the file."""
    run = fixture_p(tmp_path)
    _break(run.out / "tables" / "detection_summary.tsv", how)
    page = run.page()  # must not raise
    assert not [t for t in page.root.by_tag("table") if "data-view" in t.attrs]
    notes = [" ".join(n.text().split()) for n in page.section("detection").iter()
             if "note" in n.attrs.get("class", "").split()]  # fmt: skip
    assert [n for n in notes if "detection_summary.tsv" in n], notes


@pytest.mark.parametrize("odd", ["nan", "NaN", "inf", "-inf"])
@pytest.mark.parametrize("seed", range(4))
def test_non_finite_threshold_sorts_last(tmp_path: Path, odd: str, seed: int) -> None:
    """Review r1 L5: a ``q_threshold`` group that is no finite number sorts after the finite
    ones, which stay ascending, whatever the file order."""
    run = build("P", tmp_path)
    curve = {q: {"n_detected": 10 + k, "n_genes": 20 + k}
             for k, q in enumerate(("0.001", "0.005", "0.01", "0.05", "0.1"))}  # fmt: skip
    rows = [*curve_rows(curve), _row("tier", "q_threshold", odd, "n_detected", 1),
            _row("tier", "q_threshold", odd, "n_genes", 2)]  # fmt: skip
    random.Random(seed).shuffle(rows)
    write_summary(run.out, rows)
    _, body = grid(view(run.page(), TIER_VIEW))
    assert len(body) == len(curve) + 1, body
    assert [r[0] for r in body[:-1]] == list(curve), body
    assert body[-1][1:3] == ["1", "2"], body


# ---------------------------------------------------------------- property


THRESHOLDS = ("0.001", "0.005", "0.01", "0.02", "0.05", "0.1")


@st.composite
def _summaries(draw: st.DrawFn) -> tuple[list[dict[str, Any]], list[tuple[str, str]]]:
    qs = draw(st.lists(st.sampled_from(THRESHOLDS), min_size=1, max_size=6, unique=True))
    prob = st.one_of(st.none(), st.floats(1e-9, 1.0, allow_nan=False))
    rows: list[dict[str, Any]] = []
    for q in qs:
        n = draw(st.integers(0, 20_000))
        vals = {
            "n_detected": draw(st.integers(0, n)),
            "n_genes": n,
            "fdr": draw(prob),
            "fdr_2pep": draw(prob),
            "fp": draw(st.one_of(st.none(), st.floats(0, 5000))),
        }
        rows += [_row("tier", "q_threshold", q, m, v) for m, v in vals.items()]  # fmt: skip
    groups = draw(st.lists(st.tuples(st.sampled_from(["replicon", "feature_class"]),
                                     st.from_regex(r"[a-z]{3,8}", fullmatch=True)),
                           max_size=5, unique_by=lambda t: t[1]))  # fmt: skip
    oe = {g: draw(st.one_of(st.just(NA), st.floats(0.01, 50))) for g in groups}
    rows += basic_rows(oe) + tier_rows("mayu", 1, 1) + other_rows()
    seed = draw(st.integers(0, 2**16))
    random.Random(seed).shuffle(rows)
    return rows, list(groups)


@settings(
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(case=_summaries())
def test_property_views_pivot_the_long_file(
    case: tuple[list[dict[str, Any]], list[tuple[str, str]]],
) -> None:
    """For any long file (rows in any order, distractor sections included): the tier view is
    one row per distinct threshold, ascending, values through ``fmt_num``, all-NA columns
    dropped; the O/E view has one row per ``o_over_e`` group of ``basic``."""
    from gippa.io.run_view import RunView
    from gippa.report.build import build_report

    rows, groups = case
    with tempfile.TemporaryDirectory() as d:
        out = make_out(Path(d))
        write_preflight(out)
        write_summary(out, rows)
        page = parse(build_report(RunView.from_outdir(out), {}, {}, cfg()))
        frame = summary_values(out)
        header, body = grid(view(page, TIER_VIEW))
        assert (header, body) == expected_grid(frame)
        if groups:
            _, oe_rows = grid(view(page, OE_VIEW))
            assert len(oe_rows) == len(groups)
            for group, text in expected_oe(frame):
                hit = [r for r in oe_rows if group in r]
                assert len(hit) == 1 and text in hit[0], (group, text, oe_rows)
