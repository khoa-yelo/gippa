# ruff: noqa: RUF001, RUF002, RUF003  (the report's wording is Unicode: ≤, ≥, ², −, –, ±)
"""M12-04 PR 1, review round 1: the lead's decisions M1-M6 and L1-L4 on the run box, the cards,
the CSS and the demo banner (docs/specs/M12-04-outputs-revision.md § PR 1; PR #553).

Interface assumed (PR 1 review round 1 test author):
- M1: an instance drawn by ``base.empty_state`` (``Instance.empty``) and written by
  ``write_instance`` gives no key finding, like a skipped plot. The tests go through
  ``read_plot_files`` + ``head.findings(plots, skipped)`` and through the report; how the empty
  state is persisted (an SVG attribute, a support marker) is the implementer's choice. On
  ``sim_empty`` the box has no "0 of 0" and no "NA" either (V01's title there is
  "0 of 0 genes detected (NA); ...", an instance that is not ``empty``);
- M2: the R² card prints the SD to the mean's number of decimals, more when needed to keep one
  significant digit of the SD ("0.21 ± 0.004"); a true minus;
- M3: decoy-picked estimate "picked target–decoy gene-level FDR {fdr} at q ≤ {q}";
- M4: the protein card is exactly "N of M proteins placed" (no ``.estimate`` element); a fitted
  R² without a null p has no "unestimated" and no dangling ";";
- M5: the genes card takes its tier from ``head.tier_of`` (manifest first); an empty gene table
  gives no "0 of 0" (a "no genes" card or none);
- M6: novel-ORF and ALT card numbers go through ``fmt_q`` / ``fmt_num``; the de novo suffix is
  "; N of M de novo-only, no database or RBS support";
- L1: the input roles come from the same record as the inputs (manifest ``inputs``, else
  preflight ``inputs``);
- L2: an unreadable ``genes.tsv`` gives "genes: could not be read";
- L3: ``report.css`` has ``@media print`` with ``nav#toc`` ``position: static``, and the anchors
  have a ``scroll-margin-top``;
- L4: the demo banner comes directly after the ``<h1>``, before ``div#this-run``.
"""

from __future__ import annotations

import json
import math
import re
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.report_content import words, write_variants
from helpers.report_fixtures import Node, Page, figure, frame
from helpers.viz_out_tree import write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import (
    FIXTURES,
    GENOME_STEM,
    STEMS,
    TITLE_V01,
    TITLE_V02,
    TITLE_V07,
    Run,
    build,
    fixture_d,
    fixture_p,
    fixture_r,
    text,
    write_fitted,
    write_genes,
    write_run_manifest,
)

from gippa.viz import style

MINUS = "−"
NUMBER = re.compile(r"(?<![\w.])[−-]?\d[\d,  ]*(?:\.\d+)?")
#: the tier phrase of each tier on the genes card
CARD_TIER = {"mayu": "Mayu", "decoy_picked": "target–decoy", "two_peptide": "two-peptide rule",
             "protein_list": "protein list"}  # fmt: skip


def _box(page: Page) -> Node:
    found = page.root.by_attr("id", "this-run")
    assert len(found) == 1, f"#this-run: {len(found)} found"
    return found[0]


def _box_text(run: Run) -> str:
    return text(_box(run.page()))


def _cards(run: Run) -> dict[str, str]:
    return {k: text(v) for k, v in run.page().cards().items()}


def _card(run: Run, name: str) -> str:
    found = _cards(run)
    assert name in found, f"{run.name}: no {name} card; cards: {sorted(found)}"
    return found[name]


def _box_tier(box: str) -> str | None:
    from gippa.viz.words import TIER_TEXT

    found = [t for t, s in TIER_TEXT.items() if s.split("{")[0].strip() in box]
    return found[0] if len(found) == 1 else None


def _card_tiers(card: str) -> set[str]:
    return {t for t, phrase in CARD_TIER.items() if phrase in card}


def sig_digits(token: str) -> int:
    digits = re.sub(r"[^0-9.]", "", token)
    if "." in digits:
        return len(digits.replace(".", "").lstrip("0")) or 1
    return len(digits.rstrip("0")) or 1


# ---------------------------------------------------------------- M1: empty plots, no finding


def _write_plots(plots_dir: Path, empty: set[str]) -> dict[str, Any]:
    """V01, V07 and V02 through ``write_instance``: ``base.empty_state`` for the ids in
    ``empty`` (its message the plot's ``TITLES`` entry), a drawn instance otherwise."""
    from unit.report.test_m12_04_pr1_fixtures import TITLES

    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.viz import base

    for pid in ("V01", "V07", "V02"):
        name = STEMS[pid][4:]
        if pid in empty:
            inst = base.empty_state(pid, name, TITLES[pid], ["x"], "single")
        else:
            inst = base.Instance(suffix=None, figure=figure(), data=frame(), title=TITLES[pid])
        write_instance(inst, pid, name, plots_dir, ["svg"], 50)
    return read_plot_files(plots_dir)


def test_m1_findings_skip_an_empty_plot(tmp_path: Path) -> None:
    from gippa.report import head

    plots = _write_plots(tmp_path / "plots", empty={"V07"})
    assert head.findings(plots, {}) == [TITLE_V01, TITLE_V02]


def test_m1_findings_all_empty_gives_none(tmp_path: Path) -> None:
    from gippa.report import head

    plots = _write_plots(tmp_path / "plots", empty={"V01", "V07", "V02"})
    assert head.findings(plots, {}) == []


def test_m1_findings_keep_real_plots(tmp_path: Path) -> None:
    """The control: drawn instances written by ``write_instance`` all give their title."""
    from gippa.report import head

    plots = _write_plots(tmp_path / "plots", empty=set())
    assert head.findings(plots, {}) == [TITLE_V01, TITLE_V07, TITLE_V02]


def test_m1_box_has_no_empty_plot_title(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    run.plots = _write_plots(tmp_path / "plots", empty={"V02"})
    box = _box_text(run)
    assert TITLE_V01 in box and TITLE_V07 in box, box
    assert TITLE_V02 not in box, box


@pytest.fixture(scope="module")
def sim_empty_page(tmp_path_factory: pytest.TempPathFactory) -> Page:
    from helpers.report_fixtures import parse
    from helpers.snapshot_fixture import build_report, missing

    why = missing("sim_empty")
    if why is not None:
        pytest.fail(why)
    page = build_report("sim_empty", tmp_path_factory.mktemp("sim_empty"))
    return parse(page.read_text(encoding="utf-8"))


def test_m1_sim_empty_box_has_no_empty_state_finding(sim_empty_page: Page) -> None:
    box = text(_box(sim_empty_page))
    assert "No genes" not in box, box
    assert "No detection summary" not in box, box
    assert not re.search(r"\b0 of 0\b", box), box
    assert not re.search(r"\bNA\b", box), box


def test_m5_sim_empty_genes_card_has_no_zero_of_zero(sim_empty_page: Page) -> None:
    cards = {k: text(v) for k, v in sim_empty_page.cards().items()}
    assert not any(re.search(r"\b0 of 0\b", t) for t in cards.values()), cards
    card = cards.get("genes_detected")
    if card is not None:
        tier = _box_tier(text(_box(sim_empty_page)))
        assert _card_tiers(card) <= ({tier} if tier else set()), (tier, card)


# ---------------------------------------------------------------- M2: the R² card decimals


@pytest.mark.parametrize(
    ("r2", "sd", "shown"),
    [
        (0.21, 0.03, "0.21 ± 0.03"),
        (0.13, 0.04, "0.13 ± 0.04"),
        (-0.0174037, 0.0656274, f"{MINUS}0.017 ± 0.066"),
        (0.21, 0.004, "0.21 ± 0.004"),
        (0.21, 0.0004, "0.21 ± 0.0004"),
        (0.21, 0.123, "0.21 ± 0.12"),
    ],
)
def test_m2_sd_to_the_means_decimals(r2: float, sd: float, shown: str, tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    write_fitted(run.out, r2, sd, 0.3)
    card = _card(run, "model_r2")
    assert f"Model R²: {shown} (CV)" in card, card
    assert not re.search(r"(?<![\w])-\d", card), "a true minus, not a hyphen"


def _decimals(token: str) -> int:
    return len(token.split(".", 1)[1]) if "." in token else 0


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(r2=st.integers(-900, 900).filter(lambda x: abs(x) >= 1), sd=st.integers(1, 9000))
def test_m2_property_sd_never_rounds_to_zero(r2: int, sd: int) -> None:
    """``m ± s``: ``s`` has the decimals of ``m``, or more to keep 1 significant digit; it is
    never all zeros and is ``sd`` rounded to its decimals (values in thousandths and
    ten-thousandths, so ``%.6g`` round-trips them)."""
    m, s = r2 / 1000, sd / 100_000
    with tempfile.TemporaryDirectory() as tmp:
        run = fixture_r(Path(tmp))
        write_fitted(run.out, m, s, 0.3)
        card = _card(run, "model_r2")
    found = re.search(r"Model R²: ([−]?[\d.]+) ± ([\d.]+) \(CV\)", card)
    assert found, card
    mean, spread = found.groups()
    first_sig = -math.floor(math.log10(s))
    want = max(_decimals(mean), first_sig)
    assert _decimals(spread) == want, (m, s, card)
    assert float(spread) > 0, (m, s, card)
    assert abs(float(spread) - s) <= 0.5 * 10**-want + 1e-12, (m, s, card)


# ---------------------------------------------------------------- M3: decoy-picked wording


def test_m3_decoy_picked_says_gene_level_once(tmp_path: Path) -> None:
    card = _card(build("decoy", tmp_path), "genes_detected")
    fdr = style.fmt_q(0.0098)
    assert f"picked target–decoy gene-level FDR {fdr} at q ≤ 0.01" in card, card
    assert card.count("gene-level") == 1, card


# ---------------------------------------------------------------- M4: no "unestimated"


def test_m4_protein_card_is_exactly_n_of_m_placed(tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    node = run.page().cards()["proteins_placed"]
    assert text(node) == "106 of 106 proteins placed", text(node)
    assert not node.by_class("estimate"), "an empty estimate is not rendered"


def test_m4_protein_card_in_a_both_run(tmp_path: Path) -> None:
    card = _card(build("both", tmp_path), "proteins_placed")
    assert card == "12 of 15 proteins placed", card


@pytest.mark.parametrize(("r2", "sd"), [(0.21, 0.03), (0.05, 0.06)])
def test_m4_fitted_r2_without_null_p(r2: float, sd: float, tmp_path: Path) -> None:
    from helpers.report_content import det_rows, write_det_model

    run = fixture_r(tmp_path)
    rows = [r for r in det_rows("elastic_net", r2, sd, None) if r["metric"] != "null_p"]
    write_det_model(run.out, rows)
    card = _card(run, "model_r2")
    assert card.startswith("Model R²:") and "(CV)" in card, card
    assert "unestimated" not in card, card
    assert not card.rstrip().endswith(";") and ";" not in card, card
    assert "null p" not in card, card


# ---------------------------------------------------------------- M5: the genes card tier


def test_m5_genes_card_tier_comes_from_the_manifest(tmp_path: Path) -> None:
    """The manifest says mayu, the gene table decoy_picked: card and box both say mayu."""
    run = fixture_p(tmp_path)
    write_genes(run.out, n_detected=99, n_2pep=91, tier="decoy_picked",
                abundance=("ibaq", "molar"))  # fmt: skip
    page = run.page()
    card = text(page.cards()["genes_detected"])
    assert _box_tier(text(_box(page))) == "mayu", text(_box(page))
    assert _card_tiers(card) == {"mayu"}, card
    assert "Mayu gene-level FDR" in card, card


def test_m5_genes_card_tier_from_manifest_when_the_table_is_mixed(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    path = run.out / "tables" / "genes.tsv"
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    col = header.index("detection_tier")
    body = []
    for i, line in enumerate(lines[1:]):
        cells = line.split("\t")
        cells[col] = "two_peptide" if i % 2 else "mayu"
        body.append("\t".join(cells))
    path.write_text("\n".join([lines[0], *body]) + "\n", encoding="utf-8")
    card = _card(run, "genes_detected")
    assert _card_tiers(card) == {"mayu"}, card
    assert "unestimated" not in card, card


@pytest.mark.parametrize("name", FIXTURES)
def test_m5_genes_card_agrees_with_the_box(name: str, tmp_path: Path) -> None:
    page = build(name, tmp_path).page()
    tier = _box_tier(text(_box(page)))
    card = text(page.cards()["genes_detected"])
    assert tier is not None and _card_tiers(card) == {tier}, (tier, card)


def test_m5_empty_gene_table_gives_no_zero_of_zero(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    write_schema_table(run.out / "tables" / "genes.tsv", "genes", [])
    page = run.page()
    cards = {k: text(v) for k, v in page.cards().items()}
    assert not any(re.search(r"\b0 of 0\b", t) for t in cards.values()), cards
    card = cards.get("genes_detected")
    if card is not None:
        assert "no genes" in card.lower(), card
        assert _card_tiers(card) <= {_box_tier(text(_box(page)))}, card


# ---------------------------------------------------------------- M6: novel ORF and ALT numbers

DENOVO = "denovo_only_unsupported"


def _novel_380(out: Path) -> None:
    """380 called rows (116 flagged de novo-only), max class FDR 0.872576, and 20 not called."""
    rows = []
    for i in range(400):
        called = i < 380
        fdr = "0.872576" if i == 0 else "0.00383386" if called else "0.95"
        flags = DENOVO if called and i < 116 else "rbs_supported"
        rows.append({"orf_id": f"N{i}", "class": "intergenic", "called": called,
                     "class_fdr": fdr, "flags": flags})  # fmt: skip
    write_schema_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", rows)


def test_m6_novel_card_numbers_are_formatted(tmp_path: Path) -> None:
    run = fixture_d(tmp_path)
    _novel_380(run.out)
    page = run.page()
    card = text(page.cards()["novel_orfs"])
    assert "0.872576" not in card and "0.00383386" not in card, card
    assert any(f"class FDR ≤ {v}" in card
               for v in (style.fmt_q(0.872576), style.fmt_num(0.872576, 2),
                         style.fmt_num(0.872576, 3))), card  # fmt: skip
    assert all(sig_digits(t) <= 3 for t in NUMBER.findall(card)), card
    assert len(card) <= 160, (len(card), card)


def test_m6_novel_card_de_novo_suffix(tmp_path: Path) -> None:
    run = fixture_d(tmp_path)
    _novel_380(run.out)
    est = words(run.page().cards()["novel_orfs"].by_class("estimate")[0].text())
    assert est.endswith("; 116 of 380 de novo-only, no database or RBS support"), est
    assert "calls de novo-only without" not in est, est


def test_m6_alt_card_numbers_are_formatted(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    write_variants(run.out, [("V1", "high", "0.00383386"), ("V2", "candidate", "0.872576")])
    card = _card(run, "alt_calls")
    assert "0.00383386" not in card and "0.872576" not in card, card
    lo = {style.fmt_q(0.00383386), style.fmt_num(0.00383386, 2), style.fmt_num(0.00383386, 3)}
    hi = {style.fmt_q(0.872576), style.fmt_num(0.872576, 2), style.fmt_num(0.872576, 3)}
    assert any(v in card for v in lo) and any(v in card for v in hi), card
    assert all(sig_digits(t) <= 3 for t in NUMBER.findall(card)), card
    assert len(card) <= 160, (len(card), card)


# ---------------------------------------------------------------- L1: input roles


def test_l1_no_manifest_takes_inputs_from_preflight(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    (run.out / "run_manifest.json").unlink()
    box = _box_text(run)
    assert "input: not recorded" not in box, box
    assert "peptides (774)" in box, box
    assert GENOME_STEM in box, box


def test_l1_manifest_without_inputs_takes_them_from_preflight(tmp_path: Path) -> None:
    """A 1.0 manifest without ``inputs``, a preflight with them: genome and inputs from it."""
    run = fixture_p(tmp_path)
    write_run_manifest(run.out, roles=None, tier="mayu", modules=run.modules)
    box = _box_text(run)
    assert f"genome: {GENOME_STEM}" in box, box
    assert "genome: not recorded" not in box, box
    assert "peptides (774)" in box and "input: not recorded" not in box, box


def test_l1_no_manifest_proteins_run_keeps_its_card(tmp_path: Path) -> None:
    run = fixture_r(tmp_path)
    (run.out / "run_manifest.json").unlink()
    box = _box_text(run)
    assert "proteins (106)" in box, box
    assert "106 of 106 proteins placed" in _card(run, "proteins_placed")


# ---------------------------------------------------------------- L2: unreadable genes.tsv


def test_l2_unreadable_genes_table(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    (run.out / "tables" / "genes.tsv").write_text("not\ta genes table\n\x00\x01\n", "utf-8")
    box = _box_text(run)
    assert "genes: could not be read" in box, box
    assert "genes: not recorded" not in box, box


# ---------------------------------------------------------------- L3: CSS


def _css(page: Page) -> str:
    return " ".join(n.text() for n in page.root.by_tag("style"))


def _block(css: str, at: str) -> str | None:
    """The body of the first ``at { ... }`` block of ``css`` (nested braces counted)."""
    m = re.search(re.escape(at) + r"[^{]*\{", css)
    if m is None:
        return None
    depth, i = 1, m.end()
    while i < len(css) and depth:
        depth += {"{": 1, "}": -1}.get(css[i], 0)
        i += 1
    return css[m.end() : i - 1]


def _rules(css: str) -> list[tuple[str, str]]:
    return [(s.strip(), b) for s, b in re.findall(r"([^{}]+)\{([^{}]*)\}", css)]


def test_l3_print_css_unsticks_the_toc(tmp_path: Path) -> None:
    css = _css(fixture_p(tmp_path).page())
    body = _block(css, "@media print")
    assert body is not None, "no @media print block"
    hits = [d for s, d in _rules(body) if "#toc" in s]
    assert any(re.search(r"position:\s*static", d) for d in hits), body


def test_l3_anchors_have_scroll_margin_top(tmp_path: Path) -> None:
    page = fixture_p(tmp_path).page()
    css = _css(page)
    anchor = page.root.by_attr("id", "sec-summary")[0]
    classes = (anchor.attrs.get("class") or "").split()
    sel = re.compile(r"\[id\^=['\"]?sec-|:target|\bsection\b"
                     + "".join(rf"|\.{re.escape(c)}\b" for c in classes))  # fmt: skip
    found = [(s, d) for s, d in _rules(css) if sel.search(s) and "scroll-margin-top" in d]
    assert found, "no scroll-margin-top on the anchor targets"


# ---------------------------------------------------------------- L4: the demo banner


def test_l4_demo_banner_directly_after_h1(tmp_path: Path) -> None:
    run = fixture_p(tmp_path)
    path = run.out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    doc["demo"] = True
    path.write_text(json.dumps(doc, sort_keys=True), encoding="utf-8")
    body = run.page().root.by_tag("body")[0]
    kids = [c for c in body.children if isinstance(c, Node)]
    h1 = next(i for i, n in enumerate(kids) if n.tag == "h1")
    assert kids[h1 + 1].attrs.get("data-demo") == "true", [(n.tag, n.attrs) for n in kids[:4]]
    assert kids[h1 + 2].attrs.get("id") == "this-run", [(n.tag, n.attrs) for n in kids[:4]]


def test_l4_no_demo_box_still_first_after_h1(tmp_path: Path) -> None:
    body = fixture_p(tmp_path).page().root.by_tag("body")[0]
    kids = [c for c in body.children if isinstance(c, Node)]
    h1 = next(i for i, n in enumerate(kids) if n.tag == "h1")
    assert kids[h1 + 1].attrs.get("id") == "this-run"
    assert not body.by_attr("data-demo"), "no banner without demo: true"


# ---------------------------------------------------------------- fixtures still sane


@pytest.mark.parametrize("name", ["P", "R", "D"])
def test_box_tier_is_readable_on_every_main_fixture(name: str, tmp_path: Path) -> None:
    """The control for ``_box_tier``: the box of each fixture names its own tier."""
    run = build(name, tmp_path)
    assert _box_tier(_box_text(run)) == run.tier
