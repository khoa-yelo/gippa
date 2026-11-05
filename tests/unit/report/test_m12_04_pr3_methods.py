"""M12-04 PR 3, behaviour 1 (R-07): the report's Methods section holds the snippets of the steps
that ran and of the run's detection tier (docs/specs/M12-04-outputs-revision.md § PR 3, AC 3.2).

Fixtures P (A-like: peptides, tier ``mayu``, iBAQ, modules ``core,bgc,function``), D (C-like: de
novo peptides, tier ``two_peptide``, no abundance, ``core,bgc``), R (``protein_list``) and
``decoy`` (``decoy_picked``) come from tests/unit/report/test_m12_04_pr1_fixtures.py. Their
manifests carry the two baseline steps of ``helpers.manifest_data``; ``with_steps`` rewrites
``run_manifest.json[steps]`` with the step list of an evidence-like run (or drops the key, an old
run).

Interface pinned here (PR 3 test author, beyond the spec text):
- every snippet the Methods section shows is one element with class ``snippet`` and
  ``data-topic`` equal to the snippet's file stem: ``<topic>`` for an unkeyed block and
  ``<topic>.<key>`` for a keyed one (``detection-tiers.two_peptide``), as the template does now;
- inside that element the snippet body (the package-data HTML) appears in full and the element's
  text ends with ``Full note: <a>`` whose ``href`` is ``gippa.report.content.methods_url(topic)``
  (``https://khoa-yelo.github.io/gippa/methods/<topic>/``; the page of the key-less topic);
- the tier is the manifest's ``detection_tier``; exactly one ``detection-tiers.*`` snippet shows;
- a step is "run" when it is listed in ``run_manifest.json[steps]``; a step that maps to no
  snippet adds nothing, so a run whose steps map to none shows only its tier snippet;
- without ``steps`` (an old run) the snippets follow the tables present: ``mapping_fdr`` →
  ``placement-fdr``, ``variants`` → ``variants``, ``determinants_model`` → ``determinants``,
  ``islands`` → ``islands``; ``abundance`` needs a finite ``genes.abundance_log2`` (spec example C
  has no abundance snippet although ``genes.tsv`` is present).
"""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

import pytest
from helpers.m12_04_pr2 import METHODS_PREFIX
from helpers.report_content import words, write_variants
from helpers.report_fixtures import Node, Page, island_rows, parse
from helpers.viz_out_tree import write_schema_table
from unit.report.test_m12_04_pr1_fixtures import Run, build

ROOT = Path(__file__).resolve().parents[3]
SNIPPETS = ROOT / "src" / "gippa" / "report" / "snippets"

#: the steps an A-like run executes (core, bgc, function; peptides, no proteins)
STEPS_CORE = (
    "codon_usage", "context_layers", "context_stats", "detectability", "determinants",
    "gene_table", "loci", "manifest", "map_exact", "map_tolerant", "novel_candidates",
    "novel_orfs", "novel_promoter_layer", "novel_rbs_layer", "operon_pairs", "operons", "ostir",
    "ostir_input", "ostir_novel", "parse_bakta", "placement", "plots", "promotech",
    "promotech_novel", "promoter_layer", "promoter_windows", "rbs_layer", "region_layer",
    "report", "sixframe_index", "terminator_layer", "transterm_input", "transtermhp", "tx_layer",
    "uniop", "uniop_input", "variants",
)  # fmt: skip
STEPS_A = (*STEPS_CORE, "function")
#: steps that hold no method a snippet describes (inputs, layers, bookkeeping)
STEPS_NO_SNIPPET = ("codon_usage", "ostir", "ostir_input", "manifest", "plots", "report")

TIERS = ("decoy_picked", "mayu", "two_peptide", "protein_list")


def with_steps(run: Run, names: Iterable[str] | None) -> Run:
    """Rewrite ``run_manifest.json[steps]`` (``None`` drops the key: a pre-1.1 manifest)."""
    path = run.out / "run_manifest.json"
    doc = json.loads(path.read_text(encoding="utf-8"))
    if names is None:
        doc.pop("steps", None)
    else:
        doc["steps"] = [{"name": n, "status": "completed", "n_tasks": 1, "n_retries": 0}
                        for n in names]  # fmt: skip
    path.write_text(json.dumps(doc, sort_keys=True, indent=1), encoding="utf-8")
    return run


def snippet_nodes(page: Page) -> dict[str, Node]:
    nodes = page.section("methods").by_class("snippet")
    topics = [n.attrs.get("data-topic", "") for n in nodes]
    assert len(topics) == len(set(topics)), f"a snippet shown twice: {topics}"
    return dict(zip(topics, nodes, strict=True))


def topics_of(run: Run) -> set[str]:
    return set(snippet_nodes(run.page()))


def tier_topics(topics: Iterable[str]) -> set[str]:
    return {t for t in topics if t.startswith("detection-tiers")}


@pytest.fixture
def p_run(tmp_path: Path) -> Run:
    """P as an A-like run: its pair and novel ORF tables exist (header-only), so the Operons and
    Novel ORFs sections apply and their notes can show (review 1, M1)."""
    run = with_steps(build("P", tmp_path), STEPS_A)
    for key in ("pairs", "novel_orfs"):
        write_schema_table(run.out / "tables" / f"{key}.tsv", key, [])
    return run


@pytest.fixture
def d_run(tmp_path: Path) -> Run:
    return with_steps(build("D", tmp_path), STEPS_CORE)


# ---------------------------------------------------------------- AC 3.2


def test_ac32_fixture_d_tier_and_denovo_but_no_abundance(d_run: Run) -> None:
    topics = topics_of(d_run)
    assert "detection-tiers.two_peptide" in topics, topics
    assert "denovo-precision" in topics, topics
    assert "abundance" not in topics, topics


def test_ac32_fixture_p_mayu_not_two_peptide(p_run: Run) -> None:
    topics = topics_of(p_run)
    assert "detection-tiers.mayu" in topics, topics
    assert "detection-tiers.two_peptide" not in topics, topics


# ---------------------------------------------------------------- spec examples A and C


def test_example_a_holds_tier_placement_abundance_operons_function_islands(p_run: Run) -> None:
    """A: "the tier snippet (mayu), placement FDR, abundance, operons, function and islands"."""
    topics = topics_of(p_run)
    want = {"detection-tiers.mayu", "placement-fdr", "abundance", "operons", "function",
            "islands"}  # fmt: skip
    assert want <= topics, sorted(want - topics)


def test_example_a_has_no_denovo_precision(p_run: Run) -> None:
    """A has database peptides only: no de novo precision note (example A's list)."""
    assert "denovo-precision" not in topics_of(p_run)


def test_example_c_holds_tier_placement_denovo_islands(d_run: Run) -> None:
    topics = topics_of(d_run)
    want = {"detection-tiers.two_peptide", "placement-fdr", "denovo-precision", "islands"}
    assert want <= topics, sorted(want - topics)


def test_function_snippet_only_when_the_function_step_ran(d_run: Run) -> None:
    """C runs ``core,bgc``: no function step, so no function snippet."""
    assert "function" not in topics_of(d_run)


# ---------------------------------------------------------------- tier


@pytest.mark.parametrize(
    "name, tier", [("P", "mayu"), ("D", "two_peptide"), ("R", "protein_list"),
                   ("decoy", "decoy_picked")],
)  # fmt: skip
def test_exactly_the_runs_tier_snippet(name: str, tier: str, tmp_path: Path) -> None:
    run = with_steps(build(name, tmp_path), STEPS_CORE)
    assert tier_topics(topics_of(run)) == {f"detection-tiers.{tier}"}


def test_decoy_picked_has_its_own_block(tmp_path: Path) -> None:
    """Edge case: tier ``decoy_picked`` has its own block (not the mayu one)."""
    run = with_steps(build("decoy", tmp_path), STEPS_A)
    body = words(parse((SNIPPETS / "detection-tiers.decoy_picked.html").read_text("utf-8"))
                 .root.text())  # fmt: skip
    node = snippet_nodes(run.page())["detection-tiers.decoy_picked"]
    assert body and body in words(node.text())
    mayu = words(parse((SNIPPETS / "detection-tiers.mayu.html").read_text("utf-8")).root.text())
    assert mayu != body


# ---------------------------------------------------------------- steps


def test_steps_without_a_snippet_contribute_nothing(tmp_path: Path) -> None:
    """Edge case: a step with no snippet contributes nothing; only the tier snippet is left."""
    run = with_steps(build("P", tmp_path), STEPS_NO_SNIPPET)
    assert topics_of(run) == {"detection-tiers.mayu"}


def test_snippets_follow_the_steps_not_every_file(tmp_path: Path) -> None:
    """Before PR 3 every snippet was rendered; now a run whose steps skip islands, operons,
    variants, novel ORFs, determinants and function shows none of those."""
    run = with_steps(build("P", tmp_path), ("parse_bakta", "map_exact", "loci", "placement",
                                            "gene_table", "detectability"))  # fmt: skip
    topics = topics_of(run)
    for absent in ("islands", "operons", "variants", "novel-orfs", "determinants", "function"):
        assert absent not in topics, (absent, topics)
    assert "placement-fdr" in topics


def test_variants_novel_orfs_determinants_follow_their_steps(p_run: Run) -> None:
    topics = topics_of(p_run)
    for step_topic in ("variants", "novel-orfs", "determinants"):
        assert step_topic in topics, (step_topic, topics)


# ---------------------------------------------------------------- old manifest (no steps)


def test_fallback_without_steps_follows_the_tables_d(tmp_path: Path) -> None:
    run = with_steps(build("D", tmp_path), None)
    topics = topics_of(run)
    assert "detection-tiers.two_peptide" in topics
    assert "denovo-precision" in topics
    assert "placement-fdr" in topics
    assert "abundance" not in topics


def test_fallback_without_steps_follows_the_tables_p(tmp_path: Path) -> None:
    """P's tables: genes (with abundance), detection_summary, peptide_hits, mapping_fdr,
    determinants_model, variants (given one row here). No islands, operons, novel ORFs or
    function tables."""
    run = with_steps(build("P", tmp_path), None)
    write_variants(run.out, [("v1", "high", "0.01")])  # a header-only table counts as absent
    topics = topics_of(run)
    for present in ("detection-tiers.mayu", "placement-fdr", "abundance", "variants",
                    "determinants"):  # fmt: skip
        assert present in topics, (present, topics)
    for absent in ("islands", "operons", "novel-orfs", "function", "denovo-precision"):
        assert absent not in topics, (absent, topics)


def test_fallback_islands_appear_with_the_islands_table(tmp_path: Path) -> None:
    run = with_steps(build("P", tmp_path), None)
    assert "islands" not in topics_of(run)
    write_schema_table(run.out / "tables" / "islands.tsv", "islands", island_rows(1))
    assert "islands" in topics_of(run)


# ---------------------------------------------------------------- shape of each snippet


@pytest.mark.parametrize("name", ["P", "D", "R", "decoy"])
def test_each_snippet_shows_its_body_and_ends_with_full_note(name: str, tmp_path: Path) -> None:
    from gippa.report.content import methods_url

    run = with_steps(build(name, tmp_path), STEPS_A)
    nodes = snippet_nodes(run.page())
    assert nodes, "no snippet in the Methods section"
    for stem, node in nodes.items():
        topic = stem.split(".", 1)[0]
        src = SNIPPETS / f"{stem}.html"
        assert src.is_file(), f"shown snippet {stem} has no file {src.name}"
        body = words(parse(src.read_text(encoding="utf-8")).root.text())
        text = words(node.text())
        assert body and body in text, stem
        links = node.by_tag("a")
        assert links, f"{stem}: no Full note link"
        last = links[-1]
        assert last.attrs.get("href") == methods_url(topic) == f"{METHODS_PREFIX}{topic}/", stem
        assert text.endswith(f"Full note: {words(last.text())}"), (stem, text[-120:])
