"""M12-04 PR 10 rendering (P-17): V12 (b) writes residue J as "I/L", shows the top 15 sites and
states an unestimated class error once (docs/specs/M12-04-outputs-revision.md "PR 10" Rendering,
AC 10.2-10.4).

Figure text is every ``Text`` artist of the V12 figure (suptitle, axes titles, panel (b) lines,
legend and tick labels), plus ``Instance.title``. A panel (b) site block starts with its header
line, the one ``Text`` holding " n_alt " (``gippa.viz.variants._header``). ``variants.tsv`` keeps
``J`` (its schema says so); only figure text changes.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
from helpers import viz_variants as vv
from matplotlib.text import Text

from gippa.config import resolve_config
from gippa.io.run_view import RunView
from gippa.viz import base, style
from gippa.viz.base import Instance

ARROW = "→"


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


def _render(out: Path) -> Instance:
    import gippa.viz.variants  # noqa: F401  (registers V12)

    got = base.registered()["V12"].fn(RunView.from_outdir(out), resolve_config())
    assert isinstance(got, list) and len(got) == 1, got
    return got[0]


def _texts(inst: Instance) -> list[str]:
    return [t.get_text() for t in inst.figure.findobj(Text) if t.get_text()]


def _headers(inst: Instance) -> list[str]:
    (ax,) = [a for a in inst.figure.axes if a.get_label() == "b"]
    return [t.get_text() for t in ax.texts if " n_alt " in t.get_text()]


def _sites(n: int, *, rate: str | None = None, alt: str = "K", ref: str = "Q"
           ) -> list[dict[str, object]]:  # fmt: skip
    """``rate=None``: an NA class error, flagged ``class_error_unestimated`` (M12-04 PR 11c)."""
    flags = "snv" if rate is not None else f"snv;{vv.CE_FLAG}"
    return [
        vv.variant(f"VAR{i:04d}", aa_pos=5 + 3 * i, ref_aa=ref, alt_aa=alt, n_alt=40 - i,
                   class_error_rate=rate, flags=flags)
        for i in range(n)
    ]  # fmt: skip


# ---------------------------------------------------------------- AC 10.2: J as I/L


def test_ac_10_2_no_j_arrow_in_figure_text_and_i_l_is_used(tmp_path: Path) -> None:
    sites = [
        vv.variant("VAR0001", aa_pos=10, ref_aa="F", alt_aa="J", n_alt=5),
        vv.variant("VAR0002", aa_pos=40, ref_aa="M", alt_aa="J", n_alt=4),
        vv.variant("VAR0003", aa_pos=70, ref_aa="Q", alt_aa="K", n_alt=3),
    ]
    inst = _render(vv.build_out(tmp_path, sites, []))
    texts = [*_texts(inst), inst.title]
    joined = "\n".join(texts)
    assert f"{ARROW}J" not in joined and f"J{ARROW}" not in joined, joined
    assert f"F{ARROW}I/L" in joined, joined
    assert f"M{ARROW}I/L" in joined, joined
    assert f"Q{ARROW}K" in joined, joined


def test_the_plotted_data_keeps_the_table_letter(tmp_path: Path) -> None:
    """Only figure text is reworded; the V12 data keeps ``alt_aa = J`` as ``variants.tsv``."""
    sites = [vv.variant("VAR0001", aa_pos=10, ref_aa="F", alt_aa="J", n_alt=5)]
    inst = _render(vv.build_out(tmp_path, sites, []))
    b = inst.data[inst.data["panel"] == "b"]
    assert list(b["alt_aa"]) == ["J"]
    assert any(f"F{ARROW}I/L" in t for t in _headers(inst))


# ---------------------------------------------------------------- AC 10.3: top 15


@pytest.mark.parametrize("n", [16, 30, 40])
def test_ac_10_3_at_most_15_site_blocks_and_top_15_of_in_the_title(tmp_path: Path, n: int) -> None:
    inst = _render(vv.build_out(tmp_path, _sites(n, rate="0.01"), []))
    b = inst.data[inst.data["panel"] == "b"]
    assert len(b) == 15, len(b)
    assert list(b["variant_id"]) == [f"VAR{i:04d}" for i in range(15)]  # top-site order kept
    assert len(_headers(inst)) == 15, _headers(inst)
    assert "top 15 of" in inst.title, inst.title
    assert any("top 15 of" in t for t in _texts(inst)), inst.title
    m = re.search(r"top 15 of (\d+)", inst.title)
    assert m is not None and int(m.group(1)) == n, inst.title


# ---------------------------------------------------------------- AC 10.4: unestimated once


@pytest.mark.parametrize("n", [3, 15, 20])
def test_ac_10_4_unestimated_for_every_shown_site_is_stated_once(tmp_path: Path, n: int) -> None:
    from gippa.viz.variants import UNESTIMATED

    inst = _render(vv.build_out(tmp_path, _sites(n), []))
    joined = "\n".join(_texts(inst))
    assert joined.count(UNESTIMATED) == 1, joined
    assert inst.title.count(UNESTIMATED) == 1, inst.title
    for line in _headers(inst):
        assert UNESTIMATED not in line, line


# ---------------------------------------------------------------- review round 1 (#566)


@pytest.mark.parametrize("n", [3, 15])
def test_no_top_clause_when_every_site_is_shown(tmp_path: Path, n: int) -> None:
    """n <= 15: every ALT site is drawn, so the title has no "top 15 of" clause."""
    inst = _render(vv.build_out(tmp_path, _sites(n, rate="0.01"), []))
    assert len(_headers(inst)) == n
    assert "top " not in inst.title, inst.title
    assert not any("top 15 of" in t for t in _texts(inst))


def test_unestimated_stated_once_when_it_holds_for_every_shown_site(tmp_path: Path) -> None:
    """Lead decision (review round 1): "stated once" is judged over the shown sites. The 15
    shown sites have no class error rate; a 16th, unshown ALT site has one. The phrase
    (``UNESTIMATED``) occurs exactly once in the figure text, on no site line."""
    from gippa.viz.variants import UNESTIMATED

    sites = _sites(15)
    sites.append(vv.variant("VAR0099", aa_pos=280, n_alt=1, class_error_rate="0.01"))
    inst = _render(vv.build_out(tmp_path, sites, []))
    heads = _headers(inst)
    assert len(heads) == 15 and not any("aa 280" in h for h in heads), heads
    joined = "\n".join(_texts(inst))
    assert joined.count(UNESTIMATED) == 1, joined
    for line in heads:
        assert UNESTIMATED not in line, line
