"""M12-04 PR 4, the plot half of P-14: V12-V15 empty states read "Nothing to show: <reason>.
<what would change it>" (AC 4.8; docs/specs/M12-04-outputs-revision.md § PR 4).

AC 4.8: every empty-state title of V12-V15 starts with "Nothing to show:" and holds one more
sentence. The spec's example (A, V15) is pinned verbatim: "Nothing to show: the determinants
model needs ≥ 100 genes with abundance in ≥ 25 TUs (this run: 99 genes, 34 TUs)." V14 at function
coverage 0 is PR 8's message (it starts with the same prefix); V15 gating is PR 7b's, so only the
V15 empty-state title is tested here.

Interface assumed (PR 4 test author): the minimum in the V15 message comes from
``gippa.determinants.models.MIN_ROWS`` / ``MIN_TUS``, and the run's counts are the model rows of
``determinants_features.tsv`` (detected, measured, not imputed, not ribosomal) and their TUs,
written with ``style.fmt_count``.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from helpers import m12_04_pr4 as h
from helpers import snapshot_fixture as sf
from helpers.m12_04_pr5 import render
from helpers.report_content import write_flags
from helpers.viz_out_tree import make_out, write_schema_table
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import write_not_fitted

from gippa.determinants.models import MIN_ROWS, MIN_TUS
from gippa.viz import base, style

PIDS = ("V12", "V13", "V14", "V15")
V15_EXAMPLE = (
    "Nothing to show: the determinants model needs ≥ 100 genes with abundance in "
    "≥ 25 TUs (this run: 99 genes, 34 TUs)."
)


def check(title: str) -> None:
    """AC 4.8 and § 1 rule 2: the prefix, a further sentence, no raw token."""
    assert title.startswith(h.NOTHING + " "), title
    rest = title[len(h.NOTHING) :].strip()
    assert len(rest) >= 10 and rest.endswith("."), title
    assert rest[0].isalnum(), title
    assert h.raw_tokens([title]) == [] and not h.SNAKE.search(title), title


def _empties(result: Any) -> list[base.Instance]:
    assert not isinstance(result, base.Skip), result
    return [i for i in result if isinstance(i, base.Instance) and i.empty]


def _call(pid: str, out: Path, *overrides: str) -> Any:
    from gippa.config import resolve_config
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules

    load_plot_modules()
    cfg = resolve_config(overrides=list(overrides))
    view = RunView.from_outdir(out)
    return sf.styled(lambda: base.registered()[pid].fn(view, cfg))


# ---------------------------------------------------------------- V12, V13, V14


@pytest.mark.parametrize("relation", ["related", "same"])
def test_v12_no_alt_site(relation: str, tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "variants.tsv", "variants", [])
    found = _empties(_call("V12", out, f"genome_relation={relation}"))
    assert len(found) == 1
    check(found[0].title)


def test_v13_no_candidates(tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", [])
    found = _empties(_call("V13", out))
    assert len(found) == 1
    check(found[0].title)


@pytest.mark.parametrize("resources", [True, False])
def test_v14_no_enrichment_rows(resources: bool, tmp_path: Path) -> None:
    out = make_out(tmp_path)
    write_schema_table(out / "tables" / "enrichment_gsea.tsv", "enrichment_gsea", [])
    write_schema_table(out / "tables" / "enrichment_detection.tsv", "enrichment_detection", [])
    if not resources:
        write_flags(
            out, "function", ["layer_missing:go", "layer_missing:kegg", "layer_missing:cog"]
        )
    found = _empties(_call("V14", out))
    assert len(found) == 1
    check(found[0].title)


# ---------------------------------------------------------------- V15


def _features(out: Path, n_rows: int, n_tus: int) -> None:
    """``n_rows`` model rows over ``n_tus`` TUs, plus rows the model leaves out (undetected,
    imputed, ribosomal)."""
    rows = [
        {
            "locus_tag": f"g{i:05d}",
            "tu_id": f"TU{i % n_tus:04d}",
            "detected": True,
            "abundance_log2": 20.0 + 0.01 * i,
            "abundance_imputed": False,
            "ribosomal_protein": False,
        }
        for i in range(n_rows)
    ]
    rows += [
        {
            "locus_tag": "x1",
            "tu_id": "TUX1",
            "detected": False,
            "abundance_log2": None,
            "abundance_imputed": False,
            "ribosomal_protein": False,
        },
        {
            "locus_tag": "x2",
            "tu_id": "TUX2",
            "detected": True,
            "abundance_log2": 19.0,
            "abundance_imputed": True,
            "ribosomal_protein": False,
        },
        {
            "locus_tag": "x3",
            "tu_id": "TUX3",
            "detected": True,
            "abundance_log2": 25.0,
            "abundance_imputed": False,
            "ribosomal_protein": True,
        },
    ]
    write_schema_table(out / "tables" / "determinants_features.tsv", "determinants_features", rows)


def _v15(out: Path) -> base.Instance:
    found = _empties(render(out, ["V15"])["V15"])
    assert len(found) == 1
    return found[0]


def test_v15_spec_example_verbatim(tmp_path: Path) -> None:
    """A: 99 genes with abundance in 34 TUs, ``too_few_rows``."""
    out = make_out(tmp_path)
    write_not_fitted(out, 99)
    _features(out, 99, 34)
    inst = _v15(out)
    assert inst.title == V15_EXAMPLE
    check(inst.title)


@settings(
    max_examples=8,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.function_scoped_fixture],
)
@given(n=st.integers(min_value=2, max_value=400), tus=st.integers(min_value=2, max_value=60))
def test_property_v15_states_the_minimum_and_this_run(
    tmp_path_factory: pytest.TempPathFactory, n: int, tus: int
) -> None:
    """Plural counts (n, TUs ≥ 2) below either minimum, as ``too_few_rows`` implies."""
    tus = min(tus, n)
    assume(n < MIN_ROWS or tus < MIN_TUS)
    out = make_out(tmp_path_factory.mktemp("v15"))
    write_not_fitted(out, n)
    _features(out, n, tus)
    inst = _v15(out)
    want = (
        f"{h.NOTHING} the determinants model needs ≥ {style.fmt_count(MIN_ROWS)} genes with "
        f"abundance in ≥ {style.fmt_count(MIN_TUS)} TUs (this run: {style.fmt_count(n)} "
        f"genes, {style.fmt_count(tus)} TUs)."
    )
    assert inst.title == want


# ---------------------------------------------------------------- the empty snapshot fixture


@pytest.mark.parametrize("pid", PIDS)
def test_every_empty_state_on_the_empty_fixture(pid: str) -> None:
    why = sf.missing(sf.EMPTY)
    if why is not None:
        pytest.fail(why)
    result = render(sf.fixture_out(sf.EMPTY), [pid])[pid]
    found = _empties(result)
    assert found, f"{pid} draws no empty state on {sf.EMPTY}"
    for inst in found:
        check(inst.title)
