"""M12-04 PR 3, review round 1: the methods notes state the method correctly (H1-H4, L1-L3) and
the Methods section shows a note only where it applies (M1-M3, mutation gaps).

Interface pinned here (PR 3 review 1 test author, from the lead's decisions):
- H1 ``denovo-precision``: length classes whose reference match rate is within 0.05 of the decoy
  rate are left out; the block names the reference rate and the 0.05 margin and no longer says
  "target and decoy match rates are too close";
- H2 ``abundance``: every sentence that calls abundance "molar" names iBAQ or Top3; the sentence
  on a supplied protein table says "relative";
- H3 ``novel-orfs``: seven classes, four with a class FDR; the two-distinct-peptide rule sentence
  covers the other classes and any class whose FDR cannot be estimated (not only "without
  decoys"); no "in four classes" total;
- H4 ``determinants``: ≥ 100 genes with abundance in ≥ 25 transcription units;
- L1 ``detection-tiers`` decoy_picked: "divided by", not "scaled up by"; L2 ``placement-fdr``:
  names the strata "0, 1 and 2 differences" (digits or words); L3 ``function``: names GSEA
  (gene-set enrichment);
- M1: a topic is not shown when its report section (``TOPIC_SECTION`` below) is listed in
  ``#sec-not-applicable``; ``operons`` and ``determinants`` also need a finite
  ``genes.abundance_log2``, like ``abundance``; the tier note is always shown;
- M2: without ``run_manifest.json[steps]`` a table counts only with at least one row;
- M3: every step of ``METHOD_TOPICS`` is registered, every fallback table is a ``PATHS`` key;
- a ``failed`` step adds no note, ``steps: []`` leaves the tier note alone, the tier comes from
  the manifest or (without it) from ``genes.tsv``, and an all-``inf`` ``abundance_log2`` is no
  abundance.
"""

from __future__ import annotations

import json
import re
import tempfile
from pathlib import Path

import pandas as pd
import pytest
from helpers.m12_04_pr2 import not_applicable_entries, sim_report
from helpers.report_content import write_variants
from helpers.viz_out_tree import write_schema_table
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from unit.report.test_m12_04_pr1_fixtures import Run, build
from unit.report.test_m12_04_pr3_methods import (
    STEPS_A,
    STEPS_CORE,
    snippet_nodes,
    tier_topics,
    topics_of,
    with_steps,
)
from unit.report.test_m12_04_pr3_snippets import MAX_WORDS, blocks

#: methods topic -> the id of the report section it belongs to (gippa.report.build.SECTIONS)
TOPIC_SECTION = {
    "detection-tiers": "detection",
    "placement-fdr": "mapping",
    "denovo-precision": "mapping",
    "abundance": "detection",
    "operons": "operons",
    "islands": "genome-context",
    "variants": "variants",
    "novel-orfs": "novel-orfs",
    "function": "function",
    "determinants": "determinants",
}
NEEDS_ABUNDANCE = ("abundance", "operons", "determinants")


def block(page: str, key: str | None = None) -> str:
    """The single report block of ``page`` (``key`` for a keyed one), whitespace collapsed."""
    found = [b for k, b in blocks(page) if k == key]
    assert len(found) == 1, f"{page}.{key}: {len(found)} blocks"
    text = " ".join(found[0].split())
    assert 0 < len(text.split()) <= MAX_WORDS, f"{page}.{key}: {len(text.split())} words"
    return text


def sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?])\s+", text) if s]


def section_ids() -> dict[str, str]:
    """Section heading -> section id."""
    from gippa.report.build import SECTIONS

    return {s.heading: s.id for s in SECTIONS}


def off_sections(page: object) -> set[str]:
    ids = section_ids()
    return {ids[name] for name, _ in not_applicable_entries(page) if name in ids}  # type: ignore[arg-type]


def shown_in_off_sections(page: object) -> list[tuple[str, str]]:
    off = off_sections(page)
    topics = snippet_nodes(page)  # type: ignore[arg-type]
    return [(t, TOPIC_SECTION[t.split(".", 1)[0]]) for t in topics
            if TOPIC_SECTION.get(t.split(".", 1)[0]) in off]  # fmt: skip


def genes_path(run: Run) -> Path:
    from gippa.io.run_view import PATHS

    return run.out / PATHS["genes_tsv"]


def header_only(run: Run, key: str) -> None:
    from gippa.io.run_view import PATHS

    write_schema_table(run.out / PATHS[key], key, [])


def edit_manifest(run: Run, **changes: object) -> None:
    """Set (or, with ``None``, drop) keys of ``run_manifest.json``."""
    path = run.out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    for k, v in changes.items():
        if v is None:
            doc.pop(k, None)
        else:
            doc[k] = v
    path.write_text(json.dumps(doc, sort_keys=True, indent=1), encoding="utf-8")


# ================================================================ H1-H4 and L1-L3: the notes


def test_h1_denovo_precision_names_the_reference_rate_margin() -> None:
    text = block("denovo-precision")
    assert "target and decoy match rates are too close" not in text, text
    hits = [s for s in sentences(text) if "0.05" in s]
    assert hits, f"no 0.05 margin: {text}"
    assert any(re.search(r"reference", s, re.I) and re.search(r"decoy", s, re.I) for s in hits), (
        f"the 0.05 sentence must compare the reference rate with the decoy rate: {hits}"
    )
    assert re.search(r"\bclass", text, re.I), f"length classes are not named: {text}"
    assert re.search(r"no estimate|not estimated|is not reported|has no", text, re.I), text


def test_h2_abundance_protein_table_is_relative_unless_ibaq_or_top3() -> None:
    text = block("abundance")
    for s in sentences(text):
        if re.search(r"\bmolar\b", s, re.I):
            assert re.search(r"iBAQ|Top3", s), f"'molar' without the iBAQ/Top3 qualifier: {s!r}"
    table = [s for s in sentences(text) if re.search(r"protein table", s, re.I)]
    assert table, f"no sentence on a supplied protein table: {text}"
    assert any(re.search(r"\brelative\b", s, re.I) for s in table), table


def test_h3_novel_orfs_seven_classes_four_with_fdr() -> None:
    text = block("novel-orfs")
    assert not re.search(r"\bin four classes\b", text, re.I), text
    assert re.search(r"\b(seven|7)\b", text, re.I), f"the seven classes are not named: {text}"
    assert re.search(r"\b(four|4)\b", text, re.I), f"the four FDR classes are not named: {text}"
    assert "Without decoys, at least two distinct peptides" not in text, text
    rule = [s for s in sentences(text) if re.search(r"\b(two|2) distinct peptides", s, re.I)]
    assert rule, f"no two-distinct-peptide rule: {text}"
    assert any(re.search(r"class", s, re.I) and re.search(r"estimat", s, re.I) for s in rule), (
        f"the two-peptide rule must cover the other classes and unestimable FDRs: {rule}"
    )


def test_h4_determinants_needs_100_genes_in_25_units() -> None:
    text = block("determinants")
    assert re.search(r"\b100\b", text), text
    assert re.search(r"\b25\b[^.]*transcription units", text), text
    assert re.search(r"abundance", text, re.I), text


def test_l1_decoy_picked_divides_by_the_paired_share() -> None:
    text = block("detection-tiers", "decoy_picked")
    assert "divided by" in text, text
    assert "scaled up by" not in text, text


def test_l2_placement_fdr_names_three_strata() -> None:
    text = block("placement-fdr")
    assert re.search(r"\b(0|zero), (1|one),? and (2|two) differences", text, re.I), text


def test_l3_function_mentions_gsea() -> None:
    text = block("function")
    assert "GSEA" in text, text
    assert re.search(r"gene[- ]set enrichment", text, re.I), text


# ================================================================ M1: not-applicable sections


def test_m1_d_without_abundance_shows_no_operons_or_determinants(tmp_path: Path) -> None:
    """D has no abundance; with its pair table present the Operons section applies, so only the
    abundance rule can keep the operons note out."""
    run = with_steps(build("D", tmp_path), STEPS_CORE)
    header_only(run, "pairs")
    topics = topics_of(run)
    assert "operons" not in topics, topics
    assert "determinants" not in topics, topics
    assert "detection-tiers.two_peptide" in topics, topics


def test_m1_all_inf_abundance_is_no_abundance(tmp_path: Path) -> None:
    run = with_steps(build("P", tmp_path), STEPS_A)
    header_only(run, "pairs")
    path = genes_path(run)
    genes = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    genes["abundance_log2"] = "inf"
    genes.to_csv(path, sep="\t", index=False)
    topics = topics_of(run)
    for t in NEEDS_ABUNDANCE:
        assert t not in topics, (t, topics)


@pytest.mark.parametrize("name", ["P", "R", "D", "decoy"])
@pytest.mark.parametrize("steps", ["A", "none"])
def test_m1_no_note_in_a_not_applicable_section(name: str, steps: str, tmp_path: Path) -> None:
    run = with_steps(build(name, tmp_path), STEPS_A if steps == "A" else None)
    page = run.page()
    assert not shown_in_off_sections(page), (shown_in_off_sections(page), off_sections(page))
    assert len(tier_topics(snippet_nodes(page))) == 1


def test_m1_p_novel_orfs_note_hidden_when_its_section_is_not_applicable(tmp_path: Path) -> None:
    """P lists the novel_orfs step but writes no novel ORF table: the Novel ORFs section is not
    applicable, so its note goes too."""
    run = with_steps(build("P", tmp_path), STEPS_A)
    page = run.page()
    assert "novel-orfs" in off_sections(page)
    assert "novel-orfs" not in snippet_nodes(page)


def test_m1_sim_empty_notes_and_not_applicable_are_disjoint() -> None:
    page, _ = sim_report("sim_empty")
    topics = set(snippet_nodes(page))
    assert not shown_in_off_sections(page), shown_in_off_sections(page)
    assert len(tier_topics(topics)) == 1, topics
    if "abundance" not in topics:  # no finite abundance in this run
        assert not {"operons", "determinants"} & topics, sorted(topics)


# ================================================================ M2: fallback counts rows


def test_m2_r_without_steps_has_no_placement_fdr(tmp_path: Path) -> None:
    """R's ``mapping_fdr`` is header-only (no peptides): no placement note."""
    run = with_steps(build("R", tmp_path), None)
    topics = topics_of(run)
    assert "placement-fdr" not in topics, topics
    assert "detection-tiers.protein_list" in topics, topics


def test_m2_header_only_determinants_model_is_no_determinants(tmp_path: Path) -> None:
    run = with_steps(build("R", tmp_path), None)
    assert "determinants" in topics_of(run)  # fitted, with abundance
    header_only(run, "determinants_model")
    assert "determinants" not in topics_of(run)


def test_m2_header_only_variants_is_no_variants(tmp_path: Path) -> None:
    run = with_steps(build("P", tmp_path), None)
    write_variants(run.out, [("v1", "high", "0.01")])  # P's own table is header-only
    assert "variants" in topics_of(run)
    header_only(run, "variants")
    assert "variants" not in topics_of(run)


# ================================================================ M3: names exist


def test_m3_method_topic_steps_are_registered() -> None:
    from gippa.report.content import METHOD_TOPICS
    from gippa.steps import load_step_modules
    from gippa.steps.registry import all_steps

    load_step_modules()
    names = {s.name for s in all_steps()}
    bad = {t: [s for s in steps if s not in names] for t, (_, steps, _) in METHOD_TOPICS.items()}
    assert not any(bad.values()), bad


def test_m3_method_topic_fallback_tables_are_run_view_paths() -> None:
    from gippa.io.run_view import PATHS
    from gippa.report.content import METHOD_TOPICS

    bad = {t: [k for k in tables if k not in PATHS] for t, (_, _, tables) in METHOD_TOPICS.items()}
    assert not any(bad.values()), bad


def test_m3_topic_section_map_covers_method_topics() -> None:
    from gippa.report.content import METHOD_TOPICS

    assert set(METHOD_TOPICS) <= set(TOPIC_SECTION)
    assert set(TOPIC_SECTION.values()) <= set(section_ids().values())


# ================================================================ mutation gaps


def test_failed_step_contributes_no_note(tmp_path: Path) -> None:
    run = build("P", tmp_path)
    edit_manifest(run, steps=[{"name": "variants", "status": "failed"},
                              {"name": "determinants", "status": "completed"}])  # fmt: skip
    topics = topics_of(run)
    assert "variants" not in topics, topics
    assert "determinants" in topics, topics


def test_empty_steps_list_gives_only_the_tier_note(tmp_path: Path) -> None:
    run = with_steps(build("P", tmp_path), [])
    assert topics_of(run) == {"detection-tiers.mayu"}


def test_tier_from_the_manifest_alone(tmp_path: Path) -> None:
    run = with_steps(build("D", tmp_path), STEPS_CORE)
    path = genes_path(run)
    genes = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    genes.drop(columns=["detection_tier"]).to_csv(path, sep="\t", index=False)
    assert tier_topics(topics_of(run)) == {"detection-tiers.two_peptide"}


def test_tier_from_genes_tsv_alone(tmp_path: Path) -> None:
    run = with_steps(build("D", tmp_path), STEPS_CORE)
    edit_manifest(run, detection_tier=None)
    assert tier_topics(topics_of(run)) == {"detection-tiers.two_peptide"}


def test_manifest_tier_wins_over_genes_tsv(tmp_path: Path) -> None:
    run = with_steps(build("D", tmp_path), STEPS_CORE)
    edit_manifest(run, detection_tier="mayu")
    assert tier_topics(topics_of(run)) == {"detection-tiers.mayu"}


# ================================================================ property: steps -> notes

_P_DIR = Path(tempfile.mkdtemp(prefix="gippa-m1204-pr3-r1-"))
_P_RUN: list[Run] = []


def _p_run() -> Run:
    """One P run (pair and novel ORF tables header-only, so every section applies)."""
    if not _P_RUN:
        run = build("P", _P_DIR)
        header_only(run, "pairs")
        header_only(run, "novel_orfs")
        _P_RUN.append(run)
    return _P_RUN[0]


@settings(max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    st.dictionaries(
        st.sampled_from(sorted(STEPS_A)),
        st.sampled_from(["completed", "cached", "failed"]),
        max_size=12,
    )
)
def test_property_notes_follow_non_failed_steps(statuses: dict[str, str]) -> None:
    """For P (database peptides, finite abundance, every section applicable) the notes are the
    tier note plus each topic with a listed, non-failed step, minus de novo precision."""
    from gippa.report.content import METHOD_TOPICS

    run = _p_run()
    edit_manifest(run, steps=[{"name": n, "status": s, "n_tasks": 1, "n_retries": 0}
                              for n, s in sorted(statuses.items())])  # fmt: skip
    ran = {n for n, s in statuses.items() if s != "failed"}
    want = {t for t, (_, steps, _) in METHOD_TOPICS.items() if ran & set(steps)}
    want.discard("denovo-precision")
    page = run.page()
    assert set(snippet_nodes(page)) == want | {"detection-tiers.mayu"}
    assert not shown_in_off_sections(page)
