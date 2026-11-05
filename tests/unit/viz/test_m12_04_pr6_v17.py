# ruff: noqa: RUF003  (figure text is Unicode: −, …)
"""M12-04 PR 6, V17 labels (P-16; docs/specs/M12-04-outputs-revision.md § PR 6, AC 6.4, 6.5).

Interface assumed where the spec leaves it open (test author; reported to the lead):

- Panel (e) (axes label ``"e"``) names each feature by ``gippa.viz.words.FEATURE_LABEL``: a
  drawn text of the axes equals the label, or starts with it followed by a separator (not a
  further word), e.g. "highest transcription rate in the transcription unit — not computed: …".
- A 100 %-missing feature whose cause is a run flag carries the reason on its row (in its label
  or in a text of axes ``"e"`` at the row's height): ``layer_missing:tx_rates`` -> "not
  computed: transcription-rate layer off" (``tx_rate_max_tu``, ``tx_rate_sum_tu``),
  ``cai_reference_missing`` -> "not computed: no CAI reference" (``cai``).
- The title has correct plurals ("1 preflight warning", "1 error") and states "<X> with measured
  abundance (<Y> detected)", X and Y = ``fmt_pct`` over the non-partial genes.
- Panel (b): a contig label that does not fit is drawn inside the axes ``"b"`` or shortened
  with one middle ellipsis to at most 24 characters.
"""

from __future__ import annotations

import itertools
import re
import tempfile
from collections.abc import Sequence
from pathlib import Path

import pytest
from helpers import m12_04_pr5 as h5
from helpers import m12_04_pr6 as h
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from matplotlib.axes import Axes
from matplotlib.text import Text

from gippa.config import GippaConfig, resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance, Skip, plural
from gippa.viz.words import FEATURE_LABEL

TX_REASON = "not computed: transcription-rate layer off"
CAI_REASON = "not computed: no CAI reference"
#: the panel-(e) columns of the M9 V17 (``gippa.viz.qc.FEATURE_GROUPS``) that have a label
LABELLED = sorted(c for c in h.FEATURE_VALUES if c in FEATURE_LABEL)
RAW_TOKENS = sorted(c for c in FEATURE_LABEL if "_" in c)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@pytest.fixture
def cfg() -> GippaConfig:
    return resolve_config()


def _render(out: Path, cfg: GippaConfig) -> Instance:
    import gippa.viz.qc  # noqa: F401  (registers V17)

    result = base.registered()["V17"].fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    assert len(result) == 1
    return result[0]


def _ax(inst: Instance, label: str) -> Axes:
    found = [ax for ax in inst.figure.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in inst.figure.axes]
    return found[0]


def _genes(n: int = 10, missing: Sequence[str] = ()) -> list[dict[str, object]]:
    return [h.v17_gene(f"g{i}", detected=i < 6, measured=i < 7, missing=missing)
            for i in range(n)]  # fmt: skip


def _labelled(texts: list[Text], label: str) -> list[Text]:
    pat = re.compile(rf"^{re.escape(label)}(?![ ]?[A-Za-z0-9])")
    return [t for t in texts if pat.match(t.get_text())]


# --- panel (e): human names and reasons ----------------------------------------------------------


def test_panel_e_uses_feature_labels(tmp_path: Path, cfg: GippaConfig) -> None:
    inst = _render(h.v17_out(tmp_path, genes=_genes()), cfg)
    ax = _ax(inst, "e")
    drawn = h.ax_texts(inst.figure, ax)
    names = [t.get_text() for t in drawn]
    for col in LABELLED:
        assert col not in names, f"raw column {col!r} in panel e"
        assert _labelled(drawn, FEATURE_LABEL[col]), f"no label for {col}: {names}"


def _reason_on_row(inst: Instance, col: str, reason: str) -> bool:
    ax = _ax(inst, "e")
    drawn = h.ax_texts(inst.figure, ax)
    rows = _labelled([t for t in drawn if t in ax.get_yticklabels()], FEATURE_LABEL[col])
    assert len(rows) == 1, (col, [t.get_text() for t in drawn])
    row = rows[0]
    if reason in row.get_text():
        return True
    ticks = sorted(ax.transData.transform([(0, y) for y in ax.get_yticks()])[:, 1])
    pitch = min((b - a for a, b in itertools.pairwise(ticks)), default=10.0)
    e = row.get_window_extent()
    yc = (e.y0 + e.y1) / 2
    for t in drawn:
        if reason in t.get_text():
            te = t.get_window_extent()
            if abs((te.y0 + te.y1) / 2 - yc) <= pitch / 2:
                return True
    return False


def test_tx_rate_rows_say_the_layer_was_off(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = _genes(missing=("tx_rate_max_tu", "tx_rate_sum_tu"))
    inst = _render(h.v17_out(tmp_path, genes=genes, flags=["layer_missing:tx_rates"]), cfg)
    assert _reason_on_row(inst, "tx_rate_max_tu", TX_REASON)
    assert _reason_on_row(inst, "tx_rate_sum_tu", TX_REASON)


def test_cai_row_says_there_is_no_reference(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = _genes(missing=("cai",))
    inst = _render(h.v17_out(tmp_path, genes=genes, flags=["cai_reference_missing"]), cfg)
    assert _reason_on_row(inst, "cai", CAI_REASON)


def test_no_reason_without_the_flag_or_when_partly_measured(
    tmp_path: Path, cfg: GippaConfig
) -> None:
    """A reason is shown only for a 100 %-missing feature with its run flag."""
    genes = _genes(missing=("tx_rate_max_tu",))
    genes[0]["cai"] = None  # cai 10 % missing, flag set: no reason
    out = h.v17_out(tmp_path, genes=genes, flags=["cai_reference_missing"])
    found = h.texts(_render(out, cfg).figure)
    assert not any(TX_REASON in t or CAI_REASON in t for t in found), found


# --- AC 6.5: no raw column names, plurals ---------------------------------------------------------


def test_no_raw_column_name_in_any_text(tmp_path: Path, cfg: GippaConfig) -> None:
    genes = _genes(missing=("tx_rate_max_tu", "tx_rate_sum_tu", "cai", "at_2_10"))
    out = h.v17_out(tmp_path, genes=genes, flags=["layer_missing:tx_rates"])
    inst = _render(out, cfg)
    found = [inst.title, *h.texts(inst.figure)]
    for col in RAW_TOKENS:
        pat = re.compile(rf"(?<![\w]){re.escape(col)}(?![\w])")
        assert not any(pat.search(t) for t in found), (col, [t for t in found if col in t])


def test_one_warning_is_singular(tmp_path: Path, cfg: GippaConfig) -> None:
    """Example A: "1 preflight warnings"."""
    inst = _render(h.v17_out(tmp_path, genes=_genes()), cfg)
    found = [inst.title, *h.texts(inst.figure)]
    assert not any("1 preflight warnings" in t for t in found)
    assert inst.title.startswith("1 preflight warning,") or inst.title.startswith(
        "1 preflight warning;"
    ), inst.title


_LEVELS = st.lists(st.sampled_from(["ok", "info", "warn", "error"]), min_size=0, max_size=6)
_GENES = st.lists(st.tuples(st.booleans(), st.booleans()), min_size=1, max_size=30)


@settings(max_examples=20, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(levels=_LEVELS, genes=_GENES)
def test_property_title_plurals_and_both_shares(
    levels: list[str], genes: list[tuple[bool, bool]]
) -> None:
    """The title: plural-correct counts and "<X> with measured abundance (<Y> detected)"."""
    style.apply()
    checks = [h.check(f"rule.{i}", lvl) for i, lvl in enumerate(levels)]
    rows = [h.v17_gene(f"g{i}", detected=d, measured=m) for i, (d, m) in enumerate(genes)]
    with tempfile.TemporaryDirectory() as tmp:
        title = _render(h.v17_out(Path(tmp), checks=checks, genes=rows), resolve_config()).title
    n_warn = levels.count("warn")
    assert re.match(rf"^{re.escape(plural(n_warn, 'preflight warning'))}\b", title), title
    assert not re.search(r"(?<![\d,])1 errors\b", title), title
    measured = sum(m for _, m in genes)
    detected = sum(d for d, _ in genes)
    if measured:
        want = (f"{style.fmt_pct(measured, len(genes))} with measured abundance "
                f"({style.fmt_pct(detected, len(genes))} detected)")  # fmt: skip
        assert want in title, title


def test_example_a_title_states_both_quantities(tmp_path: Path, cfg: GippaConfig) -> None:
    """7 of 10 genes measured, 6 detected: "70 % with measured abundance (60 % detected)"."""
    title = _render(h.v17_out(tmp_path, genes=_genes()), cfg).title
    assert "70 % with measured abundance (60 % detected)" in title, title


# --- panel (b): long contig labels; AC 6.4 --------------------------------------------------------


def test_long_contig_label_goes_inside_or_is_shortened(tmp_path: Path, cfg: GippaConfig) -> None:
    replicons = [{"id": h.LONG_CONTIG, "kind": "chromosome", "length": 200_000,
                  "topology": "linear"}, h.TWO_REPLICONS[1]]  # fmt: skip
    inst = _render(h.v17_out(tmp_path, replicons=replicons, genes=_genes()), cfg)
    ax = _ax(inst, "b")
    box = ax.get_window_extent()
    found = []
    for t in h.ax_texts(inst.figure, ax):
        text = t.get_text()
        if text == h.LONG_CONTIG:
            e = t.get_window_extent()
            found.append(e.x0 >= box.x0 - 1 and e.x1 <= box.x1 + 1)
        elif h.ELLIPSIS in text and h.middle_shortened(text, h.LONG_CONTIG):
            found.append(True)
    assert found and all(found), [t.get_text() for t in h.ax_texts(inst.figure, ax)]


def test_texts_stay_in_the_figure_with_a_long_contig(tmp_path: Path, cfg: GippaConfig) -> None:
    replicons = [{"id": h.LONG_CONTIG, "kind": "chromosome", "length": 200_000,
                  "topology": "linear"}]  # fmt: skip
    genes = _genes(missing=("tx_rate_max_tu", "tx_rate_sum_tu"))
    out = h.v17_out(tmp_path, replicons=replicons, genes=genes, flags=["layer_missing:tx_rates"])
    inst = _render(out, cfg)
    assert h.out_of_bounds(inst.figure) == []
    assert any(TX_REASON in t for t in h.texts(inst.figure))


def test_sim_plots_v17_has_human_labels_in_bounds(tmp_path: Path) -> None:
    inst = h5.only(h5.render(h5.sim_copy(tmp_path), ["V17"]), "V17")
    assert h.out_of_bounds(inst.figure) == []
    found = [inst.title, *h.texts(inst.figure)]
    assert not any("1 preflight warnings" in t for t in found)
    names = [t.get_text() for t in h.ax_texts(inst.figure, _ax(inst, "e"))]
    assert not set(names) & set(FEATURE_LABEL), names
