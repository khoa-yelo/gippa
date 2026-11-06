# ruff: noqa: RUF001, RUF003  (figure text is Unicode: −, …)
"""M12-04 PR 6 review round 1 (#564): V17 reasons, the V07 twin, prose from ``words``, V18 raw
tokens, the V07 title by kind, the V18 valid-p threshold, mutant killers, ``shorten_middle`` and
small points (docs/specs/M12-04-outputs-revision.md § PR 6, § 1 rule 2; lead decisions of the
review round).

Interface pinned here (test author; reported to the lead):

- ``gippa.viz.qc.FLAG_COLUMNS``: run flag -> the panel-(e) columns it leaves uncomputed; every
  key has a ``gippa.viz.words.NOT_COMPUTED`` sentence starting "not computed: ". At least
  ``cai_reference_missing`` -> cai, fop; ``layer_missing:codon`` -> the codon-usage group;
  ``layer_missing:ostir`` -> tir_ostir; ``layer_missing:promoters`` -> promoter_score,
  promoter_score_tu; ``layer_missing:terminators`` -> terminator_conf_downstream;
  ``layer_missing:tx_rates`` -> tx_rate_max_tu, tx_rate_sum_tu.
- ``gippa.viz.words`` holds ``V07_MORE_ROWS`` (a ``{k}`` template), ``V07_HEADER`` (block ->
  header), ``V18_SMALL_NOTE``, ``V18_ALL_SMALL``, ``V18_FAMILY_NAMES`` and ``RAW_ABUNDANCE``;
  ``detection.MORE_ROWS``/``V07_HEADER`` and ``diagnostics.SMALL_NOTE``/``ALL_SMALL`` are those
  objects, ``diagnostics`` reads ``V18_FAMILY_NAMES``, and none of the prose literals is written in
  ``detection.py``/``diagnostics.py``/``qc.py``.
- ``style.shorten_middle(text, limit)``: ``limit < 1`` raises ``ValueError``; otherwise the result
  has ``min(len(text), limit)`` characters, ``text`` itself when it fits, else head + "…" + tail
  with the tail taking the odd character (limit 1 -> "…", limit 2 -> "…" + last character).
"""

from __future__ import annotations

import inspect
import itertools
import random
import re
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import m12_04_pr5 as h5
from helpers import m12_04_pr6 as h
from helpers import support_analysis_data as sa
from helpers import twins as tw
from helpers.viz_artists import markers
from helpers.viz_out_tree import write_schema_table
from hypothesis import given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.io.tsv import write_tsv
from gippa.viz import base, style, words
from gippa.viz.base import Instance, Skip
from gippa.viz.words import FEATURE_LABEL

XLABEL_PP = "observed − expected detected (percentage points)"
RAW_TOKEN = re.compile(r"\w+_\w+")
ADR = re.compile(r"ADR-\d{4}")
#: the minimum columns each run flag must cover (lead, review round 1)
FLAG_MIN: dict[str, set[str]] = {
    "cai_reference_missing": {"cai", "fop"},
    "layer_missing:codon": {"cai", "tai", "enc", "fop", "gc3", "cpb", "at_2_10",
                            "rare_codons_2_10"},
    "layer_missing:ostir": {"tir_ostir"},
    "layer_missing:promoters": {"promoter_score", "promoter_score_tu"},
    "layer_missing:terminators": {"terminator_conf_downstream"},
    "layer_missing:tx_rates": {"tx_rate_max_tu", "tx_rate_sum_tu"},
}  # fmt: skip


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _render(pid: str, out: Path, cfg: GippaConfig) -> Instance:
    from gippa.steps.plots import load_plot_modules

    load_plot_modules()
    result = base.registered()[pid].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result[0]


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in inst.figure.axes]
    return found[0]


def _hist(inst: Instance) -> dict[str, Axes]:
    return {ax.get_label()[2:]: ax for ax in inst.figure.axes if ax.get_label().startswith("b:")}


# --- 1. V17 reasons from run flags ---------------------------------------------------------------


def _labelled(texts: list[Text], label: str) -> list[Text]:
    pat = re.compile(rf"^{re.escape(label)}(?![ ]?[A-Za-z0-9])")
    return [t for t in texts if pat.match(t.get_text())]


def _reason_of_row(inst: Instance, col: str) -> str | None:
    """The "not computed: …" text on the panel-(e) row of ``col`` (its label or a text at the
    row's height), or ``None``."""
    ax = _ax(inst, "e")
    drawn = h.ax_texts(inst.figure, ax)
    rows = _labelled([t for t in drawn if t in ax.get_yticklabels()], FEATURE_LABEL[col])
    assert len(rows) == 1, (col, [t.get_text() for t in drawn])
    row = rows[0]
    if "not computed: " in row.get_text():
        return row.get_text()[row.get_text().index("not computed: ") :]
    ticks = sorted(ax.transData.transform([(0, y) for y in ax.get_yticks()])[:, 1])
    pitch = min((b - a for a, b in itertools.pairwise(ticks)), default=10.0)
    e = row.get_window_extent()
    yc = (e.y0 + e.y1) / 2
    for t in drawn:
        if t.get_text().startswith("not computed: "):
            te = t.get_window_extent()
            if abs((te.y0 + te.y1) / 2 - yc) <= pitch / 2:
                return t.get_text()
    return None


def test_flag_columns_cover_the_lead_mapping() -> None:
    from gippa.viz import qc

    every = {c for cols in qc.FEATURE_GROUPS.values() for c in cols}
    for flag, cols in FLAG_MIN.items():
        assert flag in qc.FLAG_COLUMNS, sorted(qc.FLAG_COLUMNS)
        assert cols <= set(qc.FLAG_COLUMNS[flag]), (flag, qc.FLAG_COLUMNS[flag])
    for flag, cols in qc.FLAG_COLUMNS.items():
        assert set(cols) <= every, (flag, cols)
        assert words.NOT_COMPUTED[flag].startswith("not computed: "), flag


def test_cai_reference_missing_names_cai_and_fop(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = [h.v17_gene(f"g{i}", detected=True, measured=True, missing=("cai", "fop"))
             for i in range(6)]  # fmt: skip
    inst = _render("V17", h.v17_out(tmp_path, genes=genes, flags=["cai_reference_missing"]), cfg)
    for col in ("cai", "fop"):
        assert _reason_of_row(inst, col) == "not computed: no CAI reference", col


@pytest.mark.parametrize("flag", sorted(f for f in FLAG_MIN if f.startswith("layer_missing:")))
def test_layer_flag_gives_a_reason_on_every_column(
    tmp_path: Path, cfg: GippaConfig, flag: str
) -> None:
    cols = sorted(FLAG_MIN[flag])
    genes = [h.v17_gene(f"g{i}", detected=True, measured=True, missing=cols) for i in range(6)]
    inst = _render("V17", h.v17_out(tmp_path, genes=genes, flags=[flag]), cfg)
    want = words.NOT_COMPUTED[flag]
    for col in cols:
        assert _reason_of_row(inst, col) == want, (flag, col)


# --- 2. the V07 twin follows the static plot ------------------------------------------------------


def _v7_points(s: dict[str, Any]) -> list[dict[str, Any]]:
    out = []
    for t in tw.traces(s, "scatter", "scattergl"):
        xs, ys = tw.values(t.get("x")), tw.values(t.get("y"))
        ex = t.get("error_x") or {}
        plus, minus = tw.values(ex.get("array")), tw.values(ex.get("arrayminus"))
        for k, (x, y) in enumerate(zip(xs, ys, strict=True)):
            out.append({"x": x, "y": str(y), "plus": plus[k] if plus else None,
                        "minus": minus[k] if minus else None})  # fmt: skip
    return out


def test_twin_plots_pp_when_present(tmp_path: Path, cfg: GippaConfig) -> None:
    from gippa.viz import detection

    inst = _render("V07", h.v07_out(tmp_path), cfg)
    data = tw.instance_frame(inst.data)
    s = tw.spec("V07", data)
    assert tw.title_of(s["layout"]["xaxis"]) == detection.V07_XLABEL_PP
    pts = [p for p in _v7_points(s) if p["x"] is not None]
    for r in data.to_dict("records"):
        if pd.isna(r["o_minus_e_pp"]):
            continue
        hit = [p for p in pts if p["y"].startswith(str(r["group"]))]
        assert len(hit) == 1, (r["group"], pts)
        assert hit[0]["x"] == pytest.approx(r["o_minus_e_pp"])
        if not (pd.isna(r["ci_low_pp"]) or pd.isna(r["ci_high_pp"])) and hit[0]["plus"]:
            assert hit[0]["plus"] == pytest.approx(r["ci_high_pp"] - r["o_minus_e_pp"])


def test_twin_of_a_1_0_tsv_still_plots_genes() -> None:
    from gippa.viz import detection

    data = tw.v07()  # the 1.0 columns, no _pp
    s = tw.spec("V07", data)
    assert tw.title_of(s["layout"]["xaxis"]) == detection.V07_XLABEL
    pts = [p for p in _v7_points(s) if p["x"] is not None]
    assert sorted(p["x"] for p in pts) == pytest.approx(sorted(data["o_minus_e"].dropna()))


# --- 3. no inline prose ---------------------------------------------------------------------------


def test_prose_constants_live_in_words() -> None:
    from gippa.viz import detection, diagnostics

    assert "{k}" in words.V07_MORE_ROWS
    assert words.V07_MORE_ROWS.format(k="3") == "3 more rows in the data TSV"
    assert detection.MORE_ROWS is words.V07_MORE_ROWS
    assert dict(words.V07_HEADER) == {
        "replicon": "Replicon",
        "feature_class": "Feature class",
        "cog_category": "COG category",
        "region": "Region",
    }
    assert detection.V07_HEADER is words.V07_HEADER
    assert words.V18_SMALL_NOTE == "Not shown (fewer than 10 tests): "
    assert diagnostics.SMALL_NOTE is words.V18_SMALL_NOTE
    assert words.V18_ALL_SMALL == (
        "Nothing to show: no test family has 10 or more tests in this run."
    )
    assert diagnostics.ALL_SMALL is words.V18_ALL_SMALL
    assert isinstance(words.V18_FAMILY_NAMES, Mapping)
    assert words.RAW_ABUNDANCE == "raw abundance"


def _source(name: str) -> str:
    import importlib

    return inspect.getsource(importlib.import_module(f"gippa.viz.{name}"))


@pytest.mark.parametrize(
    ("module", "literals"),
    [
        ("detection", ("more rows in the data TSV", "\"Replicon\"", "\"Feature class\"")),
        ("diagnostics", ("Not shown (fewer", "Nothing to show", "Region tests", "Join count",
                         "Moran's I", "gene-label", "TU-block")),
        ("qc", ("raw abundance", "not computed")),
    ],
)  # fmt: skip
def test_no_prose_literal_in_the_plot_modules(module: str, literals: tuple[str, ...]) -> None:
    src = _source(module)
    assert [lit for lit in literals if lit in src] == []


def test_diagnostics_reads_the_family_names_from_words() -> None:
    assert "V18_FAMILY_NAMES" in _source("diagnostics")


# --- 4. V18: no raw tokens, no ADR in any figure text ---------------------------------------------


def _add_determinant_nulls(out: Path) -> None:
    rng = random.Random(11)
    rows = [{"model": "elastic_net", "metric": "cv_r2_mean", "value": 0.21, "n": 500,
             "flags": "null_uncalibrated"},
            {"model": "elastic_net", "metric": "null_p", "value": 0.002, "n_perm": 999,
             "flags": "null_uncalibrated"}]  # fmt: skip
    write_schema_table(out / "tables" / "determinants_model.tsv", "determinants_model", rows)
    path = out / "plots" / "support" / "null_determinants.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, list(sa.COLUMNS["null_determinants"]),
              [["elastic_net", b, repr(round(rng.gauss(0, 0.05), 4)), repr(0.21)]
               for b in range(1, 200)])  # fmt: skip


def _region_ids(out: Path) -> set[str]:
    """The run's own region and contig IDs (``annotation/regions.tsv``,
    ``tables/regions_stats.tsv``, ``tables/islands.tsv``): identifiers from the data, which
    V18 text may show (lead, review rounds 1 and 2)."""
    ids: set[str] = set()
    for rel, col in (("annotation/regions.tsv", "region_id"),
                     ("tables/regions_stats.tsv", "region_id"),
                     ("tables/islands.tsv", "contig")):  # fmt: skip
        if (out / rel).is_file():
            ids |= {r[col] for r in h5.read_table(out, rel) if r.get(col)}
    return ids


def _bad_tokens(inst: Instance, out: Path) -> list[str]:
    """Texts with an ADR number or a ``\\w+_\\w+`` token that is not one of the run's region IDs."""
    allowed = _region_ids(out)
    found = [inst.title, *h.texts(inst.figure)]

    def ok(m: str) -> bool:  # a data ID, or a piece of one shortened with an ellipsis
        return any(m in a for a in allowed)

    raw = [t for t in found if not all(ok(m) for m in RAW_TOKEN.findall(t))]
    return raw + [t for t in found if ADR.search(t) and t not in raw]


def test_v18_hand_fixture_has_no_raw_token_or_adr(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_regions=12, n_join_tu=12, n_moran_gl=3,
                    n_pairs=12)  # fmt: skip
    h.add_island_nulls(out)
    _add_determinant_nulls(out)
    assert _bad_tokens(_render("V18", out, cfg), out) == []


def test_v18_sim_plots_has_no_raw_token_or_adr(tmp_path: Path) -> None:
    out = h5.sim_copy(tmp_path)
    inst = h5.only(h5.render(out, ["V18"]), "V18")
    assert "prophage_1" in _region_ids(out)
    assert _bad_tokens(inst, out) == []


def test_uncalibrated_note_has_no_adr() -> None:
    from gippa.viz import null_panels

    assert not ADR.search(null_panels.UNCALIBRATED_NOTE), null_panels.UNCALIBRATED_NOTE


# --- 5. V07 title by kind -------------------------------------------------------------------------


def _basic_rows(gtype: str, names: list[str]) -> list[dict[str, object]]:
    return [r for g in names
            for r in h.basic(gtype, g, n_genes=20, n_detected=10, expected=10.0,
                             ci=(-2.0, 2.0))]  # fmt: skip


def test_v07_title_counts_each_kind(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = [
        *_basic_rows("replicon", ["chr1", "plasmid_1"]),
        *_basic_rows("feature_class", ["cds", "hypothetical", "pseudogene", "sorf"]),
        *_basic_rows("cog_category", list("CEGJKLMSTV")),
    ]
    inst = _render("V07", h.v07_out(tmp_path, summary=summary), cfg)
    assert inst.title == ("No group departs from detectability at q ≤ 0.05 (2 regions, "
                          "2 replicons, 4 feature classes, 10 COG categories tested)")  # fmt: skip


def test_v07_title_singular_kinds(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render("V07", h.v07_out(tmp_path, regions=h.A_REGIONS[:1]), cfg)
    assert inst.title == ("No group departs from detectability at q ≤ 0.05 (1 region, "
                          "1 replicon, 1 feature class, 1 COG category tested)")  # fmt: skip


def test_v07_title_leaves_out_kinds_with_zero(tmp_path: Path, cfg: GippaConfig) -> None:
    summary = [*_basic_rows("replicon", ["chr1"]), *_basic_rows("feature_class", ["cds", "sorf"])]
    inst = _render("V07", h.v07_out(tmp_path, summary=summary), cfg)
    assert inst.title == ("No group departs from detectability at q ≤ 0.05 (2 regions, "
                          "1 replicon, 2 feature classes tested)")  # fmt: skip
    only_regions = _render("V07", h.v07_out(tmp_path / "r", summary=[]), cfg)
    assert (
        only_regions.title == "No group departs from detectability at q ≤ 0.05 (2 regions tested)"
    )
    assert "classes tested" not in only_regions.title


# --- 6. V18: the threshold counts only p in [0, 1] ------------------------------------------------


def _with_outside(out: Path, valid: int, outside: list[float]) -> None:
    rows = [h.region_p(p, i) for i, p in enumerate(h.spread(valid))]
    rows += [h.region_p(p, 100 + i) for i, p in enumerate(outside)]
    write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats", rows)


def test_eight_valid_and_three_outside_is_not_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_join_tu=12)
    _with_outside(out, 8, [1.5, 2.0, -0.2])
    inst = _render("V18", out, cfg)
    assert "regions_stats / o_minus_e" not in _hist(inst)
    notes = [t for t in h.texts(inst.figure) if t.startswith("Not shown")]
    assert len(notes) == 1 and "region O − E (8)" in " ".join(notes[0].split()), notes


def test_ten_valid_and_three_outside_is_drawn(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_join_tu=12)
    _with_outside(out, 10, [1.5, 2.0, -0.2])
    inst = _render("V18", out, cfg)
    assert "regions_stats / o_minus_e" in _hist(inst)


# --- 7. surviving mutants -------------------------------------------------------------------------


def test_m8_all_small_with_null_draws_is_not_empty(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=True, n_regions=2, n_join_tu=4)
    h.add_island_nulls(out)
    inst = _render("V18", out, cfg)
    assert not inst.empty
    assert any(t.startswith("Not shown (fewer than 10 tests): ") for t in h.texts(inst.figure))
    assert any(ax.get_label().startswith("c:islands") for ax in inst.figure.axes)


def test_m12_ties_at_the_row_cap_go_to_table_order(tmp_path: Path, cfg: GippaConfig) -> None:
    """60 regions: 25 at |pp| 30, 30 at 20, 5 at 10 (interleaved, signs mixed). 40 are drawn:
    the 25 largest and the first 15 of the 20-pp tie in table order."""
    regions = []
    for i in range(60):
        level = 30.0 if i % 12 < 5 else 20.0 if i % 12 < 11 else 10.0
        sign = 1 if i % 2 else -1
        effect = sign * level / 4  # n = 25: pp = 4 · effect
        regions.append(h.region(f"prophage_{i:02d}", "prophage", n_eligible=25, observed=5.0,
                                expected=5.0 - effect, effect=effect,
                                ci=(effect - 0.5, effect + 0.5), q=0.5))  # fmt: skip
    inst = _render("V07", h.v07_out(tmp_path, summary=[], regions=regions), cfg)
    recs = inst.data.to_dict("records")
    order = sorted(range(len(recs)), key=lambda i: (-abs(h.num(recs[i]["o_minus_e_pp"])), i))
    want = sorted(order[:40])
    labels = [t.get_text().split(" (n = ")[0] for t in _forest_axes(inst).get_yticklabels()]
    assert sorted(labels) == sorted(str(recs[i]["group"]) for i in want)


def _forest_axes(inst: Instance) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_xlabel() == XLABEL_PP]
    assert len(found) == 1
    return found[0]


def test_m15_cog_row_label_is_exact(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render("V07", h.v07_out(tmp_path), cfg)
    labels = [t.get_text() for t in _forest_axes(inst).get_yticklabels()]
    assert "J Translation (n = 12)" in labels, labels


# --- 8. style.shorten_middle ----------------------------------------------------------------------


def test_shorten_middle_exact_length_is_unchanged() -> None:
    assert style.shorten_middle("abcdef", 6) == "abcdef"
    assert style.shorten_middle("abcdefg", 6) == "ab…efg"


def test_shorten_middle_limit_one_and_two() -> None:
    assert style.shorten_middle("abcdef", 1) == "…"
    assert style.shorten_middle("abcdef", 2) == "…f"
    assert style.shorten_middle("a", 1) == "a"


@pytest.mark.parametrize("limit", [0, -1, -24])
def test_shorten_middle_limit_below_one_raises(limit: int) -> None:
    with pytest.raises(ValueError):
        style.shorten_middle("abcdef", limit)


@settings(max_examples=60, deadline=None)
@given(
    text=st.text(st.characters(blacklist_characters="…"), max_size=60),
    limit=st.integers(1, 40),
)
def test_property_shorten_middle(text: str, limit: int) -> None:
    out = style.shorten_middle(text, limit)
    assert len(out) == min(len(text), limit)
    if len(text) <= limit:
        assert out == text
    else:
        head, _, tail = out.partition("…")
        assert out.count("…") == 1
        assert text.startswith(head) and text.endswith(tail)
        assert len(tail) - len(head) in (0, 1)


# --- 9. small points ------------------------------------------------------------------------------


def test_v17_title_without_a_detected_column(tmp_path: Path, cfg: GippaConfig) -> None:
    """A gene table without ``detected`` (RunView rejects such a genes.tsv, so the view is
    wrapped): the title states the measured share only, never "0 % detected"."""
    import gippa.viz.qc  # noqa: F401  (registers V17)

    genes = [h.v17_gene(f"g{i}", detected=True, measured=i < 7) for i in range(10)]
    view = RunView.from_outdir(h.v17_out(tmp_path, genes=genes))
    table = view.table

    def no_detected(key: str) -> pd.DataFrame | None:
        frame = table(key)
        return frame.drop(columns=["detected"]) if key == "genes_tsv" else frame

    view.table = no_detected  # type: ignore[method-assign]
    result = base.registered()["V17"].fn(view, cfg)
    assert not isinstance(result, Skip)
    title = result[0].title
    assert title.endswith("; 70 % with measured abundance"), title
    assert "detected" not in title


def test_v17_panel_e_rows_sorted_by_human_label(tmp_path: Path, cfg: GippaConfig) -> None:
    from gippa.viz import qc

    genes = [h.v17_gene(f"g{i}", detected=True, measured=True) for i in range(6)]
    inst = _render("V17", h.v17_out(tmp_path, genes=genes), cfg)
    raw = getattr(words, "RAW_ABUNDANCE", "raw abundance")  # pinned in the words test
    group_of = {(raw if c == "abundance_raw" else words.feature_label(c)): g
                for g, cols in qc.FEATURE_GROUPS.items() for c in cols}  # fmt: skip
    ax = _ax(inst, "e")
    h.draw(inst.figure)
    ticks = sorted(ax.get_yticklabels(), key=lambda t: -t.get_window_extent().y0)
    labels = [t.get_text() for t in ticks if t.get_text()]
    assert set(labels) == set(group_of), labels
    runs: list[tuple[str, list[str]]] = []
    for lab in labels:
        if not runs or runs[-1][0] != group_of[lab]:
            runs.append((group_of[lab], []))
        runs[-1][1].append(lab)
    assert len({g for g, _ in runs}) == len(runs), runs  # each group contiguous
    for g, labs in runs:
        assert labs == sorted(labs, key=str.casefold), (g, labs)


def test_v07_x_limits_come_from_the_drawn_rows(tmp_path: Path, cfg: GippaConfig) -> None:
    """45 regions with |pp| 8..52 and ±1 pp CIs, plus one at 0 pp with a ±90 pp CI: that row is
    cut (smallest |pp|) and must not stretch the axis."""
    regions = []
    for i in range(45):
        effect = (8 + i) / 4 * (1 if i % 2 else -1)  # n = 25
        regions.append(h.region(f"prophage_{i:02d}", "prophage", n_eligible=25, observed=5.0,
                                expected=5.0 - effect, effect=effect,
                                ci=(effect - 0.25, effect + 0.25), q=0.5))  # fmt: skip
    regions.append(h.region("prophage_wide", "prophage", n_eligible=25, observed=5.0,
                            expected=5.0, effect=0.0, ci=(-22.5, 22.5), q=0.5))  # fmt: skip
    inst = _render("V07", h.v07_out(tmp_path, summary=[], regions=regions), cfg)
    ax = _forest_axes(inst)
    assert all(abs(m.x) > 0 for m in markers(ax))  # the wide row is not drawn
    lo, hi = ax.get_xlim()
    assert lo > -70 and hi < 70, (lo, hi)


def test_v18_small_note_has_a_right_margin(tmp_path: Path, cfg: GippaConfig) -> None:
    out = h.v18_out(tmp_path, constant=False, n_regions=3, n_hl=4, n_join_tu=12, n_moran_gl=5,
                    n_pairs=6)  # fmt: skip
    inst = _render("V18", out, cfg)
    fig = inst.figure
    notes = [t for t in h.drawn_texts(fig) if t.get_text().startswith("Not shown")]
    assert len(notes) == 1
    assert notes[0].get_window_extent().x1 <= fig.bbox.x1 - 8, notes[0].get_window_extent()


@pytest.mark.parametrize("seed", [1, 2])
def test_v18_note_margin_property(tmp_path: Path, cfg: GippaConfig, seed: int) -> None:
    rng = random.Random(seed)
    sizes = {k: rng.randint(1, 9) for k in ("n_regions", "n_hl", "n_moran_gl", "n_pairs")}
    with tempfile.TemporaryDirectory() as tmp:
        inst = _render("V18", h.v18_out(Path(tmp), constant=False, n_join_tu=12, **sizes), cfg)
    notes = [t for t in h.drawn_texts(inst.figure) if t.get_text().startswith("Not shown")]
    assert notes and all(n.get_window_extent().x1 <= inst.figure.bbox.x1 - 8 for n in notes)
