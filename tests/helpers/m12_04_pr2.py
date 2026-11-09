# ruff: noqa: RUF003  (the report's wording is Unicode: −, ≤, π₀, ℓ̂)
"""Shared pieces of the M12-04 PR 2 tests: captions, plain wording, the curated genes preview,
conditional sections, empty-state notices and the glossary
(docs/specs/M12-04-outputs-revision.md § PR 2).

Builds on the PR 1 fixtures (tests/unit/report/test_m12_04_pr1_fixtures.py: P, R, D, P_nokegg,
P_core, decoy, both) and the M9-25 snapshot fixtures (``sim_plots``, ``sim_empty``).

Interface pinned here (PR 2 test author, beyond the spec text):
- ``QUESTIONS``: the 16 §10.2 "Question" cells, as a literal (the code keeps none; lead decision
  2). ``gippa.report.content.GUIDES`` (M12-05c, replacing ``CAPTIONS``) maps each registered
  plot ID to its structured guide, whose ``question`` is the Question;
- the AC 2.1 flag vocabulary is every key of ``gippa.viz.words.FLAG_TEXT``, ``SKIP_TEXT`` and
  ``SCOPED_TEXT`` that holds an underscore or a colon (lead decision 1). A ``prefix:`` key
  matches the prefix followed by a parameter (``layer_missing:kegg``);
- prose is the text of ``p, li, figcaption, .card, .notice, h2, h3, dd`` (AC 2.1);
- ``dress(run, tmp)`` gives a PR 1 fixture what ``gippa.steps.plots`` would pass the report: the
  skip reasons of ``skipped_plots(view)`` and, for every other plot of a section, one small
  drawn instance written by ``write_instance`` (so every slot has a ``<figure>``).
"""

from __future__ import annotations

import functools
import re
import tempfile
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

from helpers.report_fixtures import NAMES, Node, Page, figure, frame, parse

REPO = Path(__file__).resolve().parents[2]
METHODS_DIR = REPO / "docs" / "methods"

#: GIPPA_SPEC.md §10.2, the "Question" column, verbatim (IDs zero-padded as registered)
QUESTIONS: Mapping[str, str] = {
    "V01": "Where in the genome is the proteome?",
    "V02": "Is abundance clustered along the chromosome?",
    "V03": "What happens at this prophage, BGC, island, flagged operon or novel ORF?",
    "V05": "Where did my sequences go?",
    "V06": "Can I trust tolerant hits and de novo peptides?",
    "V07": "Which regions and categories are detected less than detectability predicts?",
    "V08": "Which functions dominate abundance?",
    "V10": "Is there structure beyond operons?",
    "V11": "Do operon partners look alike, and which other neighbours do too?",
    "V12": "Where do peptides support ALT residues?",
    "V13": "Which new genes or upstream extensions are supported?",
    "V14": "Which functions are high or low?",
    "V15": "What explains abundance?",
    "V16": "How do search strategies differ?",
    "V17": "Are the inputs sane?",
    "V18": "Are the tests calibrated on this dataset?",
}

#: the report's section plots (gippa.report.build.SECTIONS, M9-02), for ``dress``
REPORT_PLOTS = ("V17", "V05", "V06", "V07", "V18", "V01", "V02", "V10", "V03", "V11", "V12",
                "V13", "V08", "V14", "V15")  # fmt: skip

#: AC 2.1: internal references never allowed in prose
INTERNAL_REFS = (re.compile(r"ADR-\d{4}"), re.compile(r"§\s?\d"), re.compile("GIPPA_SPEC"),
                 re.compile("docs/methods"))  # fmt: skip

#: behaviour 2: an absolute docs-site methods link (Q12)
METHODS_URL = re.compile(r"^https://khoa-yelo\.github\.io/gippa/methods/([A-Za-z0-9_-]+)/$")
METHODS_PREFIX = "https://khoa-yelo.github.io/gippa/methods/"

PROSE_TAGS = frozenset({"p", "li", "figcaption", "h2", "h3", "dd"})
PROSE_CLASSES = ("card", "notice")


def words(text: str) -> str:
    return " ".join(text.split())


#: skip reasons outside ``gippa.viz.words`` that prose must not show either (review r1 L3:
#: ``gippa.steps.plots.skipped_plots`` gives ``not_registered``)
EXTRA_TOKENS = ("not_registered",)


@functools.cache
def flag_vocabulary() -> tuple[str, ...]:
    """Flag and skip tokens with an underscore or a colon (lead decision 1), plus
    ``EXTRA_TOKENS``."""
    from gippa.viz.words import FLAG_TEXT, SCOPED_TEXT, SKIP_TEXT

    keys = set(FLAG_TEXT) | set(SKIP_TEXT) | {k for m in SCOPED_TEXT.values() for k in m}
    keys |= set(EXTRA_TOKENS)
    return tuple(sorted(k for k in keys if "_" in k or ":" in k))


@functools.cache
def _flag_patterns() -> tuple[tuple[str, re.Pattern[str]], ...]:
    out = []
    for tok in flag_vocabulary():
        if tok.endswith(":"):  # a prefix template: the prefix with a parameter
            pat = re.compile(r"(?<![\w:])" + re.escape(tok) + r"(?=[\w./-])")
        else:
            pat = re.compile(r"(?<![\w:])" + re.escape(tok) + r"(?![\w:])")
        out.append((tok, pat))
    return tuple(out)


def flag_tokens(text: str) -> list[str]:
    """The flag-vocabulary tokens that occur in ``text`` as whole words."""
    return [tok for tok, pat in _flag_patterns() if pat.search(text)]


def internal_refs(text: str) -> list[str]:
    return [m.group(0) for pat in INTERNAL_REFS for m in pat.finditer(text)]


def prose_nodes(root: Node) -> list[Node]:
    """``p, li, figcaption, .card, .notice, h2, h3, dd`` under ``root`` (AC 2.1)."""
    return [n for n in root.iter()
            if n.tag in PROSE_TAGS or any(n.has_class(c) for c in PROSE_CLASSES)]  # fmt: skip


def prose_violations(page: Page) -> list[str]:
    """``"<tag>: <hit> in <text>"`` for every internal reference or flag token in prose."""
    bad = []
    for n in prose_nodes(page.root):
        text = words(n.text())
        for hit in internal_refs(text) + flag_tokens(text):
            bad.append(f"<{n.tag} class={n.attrs.get('class', '')!r}> {hit!r} in: {text[:160]}")
    return sorted(set(bad))


def methods_topics() -> set[str]:
    """The pages under docs/methods (file stems)."""
    return {p.stem for p in METHODS_DIR.glob("*.md")}


def methods_links(root: Node) -> list[str]:
    return [n.attrs["href"] for n in root.by_tag("a")
            if n.attrs.get("href", "").startswith(METHODS_PREFIX)]  # fmt: skip


# ---------------------------------------------------------------- fixtures


def write_fake_plots(plots_dir: Path, pids: Iterable[str],
                     empty: Mapping[str, str] | None = None) -> dict[str, Any]:  # fmt: skip
    """One instance per plot ID through ``write_instance`` (SVG only): ``base.empty_state`` with
    the message ``empty[pid]`` for the IDs in ``empty``, a small drawn figure otherwise."""
    from gippa.io.plot_files import read_plot_files, write_instance
    from gippa.viz import base

    empty = empty or {}
    for pid in pids:
        name = NAMES[pid]
        if pid in empty:
            inst = base.empty_state(pid, name, empty[pid], ["x"], "single")
        else:
            inst = base.Instance(suffix=None, figure=figure(), data=frame(),
                                 title=f"{pid}: a drawn finding")  # fmt: skip
        write_instance(inst, pid, name, plots_dir, ["svg"], 50)
    return dict(read_plot_files(plots_dir))


def skipped(view: Any) -> dict[str, str]:
    from gippa.steps.plots import load_plot_modules, skipped_plots

    load_plot_modules()
    return dict(skipped_plots(view))


def dress(run: Any, tmp: Path, empty: Mapping[str, str] | None = None) -> Any:
    """``run`` with the real skip reasons and a drawn (or empty-state) instance for every other
    report plot; returns ``run``."""
    reasons = skipped(run.view())
    drawn = [p for p in REPORT_PLOTS if p not in reasons]
    run.skipped = reasons
    run.plots = write_fake_plots(tmp / f"plots-{run.name}", drawn, empty)
    return run


@functools.cache
def sim_report(fixture: str) -> tuple[Page, Path]:
    """The ``report.html`` of a snapshot fixture (plots rendered for real; once per process),
    with the directory holding its ``plots/``."""
    from helpers.snapshot_fixture import build_report, missing

    why = missing(fixture)
    if why is not None:
        raise RuntimeError(why)
    tmp = Path(tempfile.mkdtemp(prefix=f"gippa-m1204-pr2-{fixture}-"))
    page = build_report(fixture, tmp)
    return parse(page.read_text(encoding="utf-8")), tmp


def section(page: Page, sid: str) -> Node | None:
    found = [s for s in page.sections() if s.attrs.get("id") == sid]
    return found[0] if found else None


def not_applicable(page: Page) -> Node | None:
    found = [s for s in page.sections() if s.attrs.get("id") == "sec-not-applicable"]
    assert len(found) <= 1, f"{len(found)} section#sec-not-applicable"
    return found[0] if found else None


def not_applicable_entries(page: Page) -> list[tuple[str, str]]:
    """``(name, reason)`` of every ``li`` of ``#sec-not-applicable`` ("<Section>: <reason>")."""
    sec = not_applicable(page)
    if sec is None:
        return []
    out = []
    for li in sec.by_tag("li"):
        text = words(li.text())
        name, sep, reason = text.partition(": ")
        assert sep, f"not-applicable entry is not '<Section>: <reason>': {text!r}"
        out.append((name, reason))
    return out
