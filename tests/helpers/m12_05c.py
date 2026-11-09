"""Shared pieces of the M12-05c tests: the structured plot guides of the HTML report
(docs/specs/M12-05c-plot-guides.md, lead decisions 2026-11-09: the five labels, always shown, no
extra heading, no ``<details>``).

Interface pinned here (M12-05c PR 1 test author, beyond the spec text):
- ``gippa.report.content.GUIDES: Mapping[str, Guide]``; ``Guide`` a frozen dataclass with the
  fields ``question, panels, read, caveats`` in that order (decision 7); ``content.question(pid)``
  returns ``GUIDES[pid].question``;
- a figure's guide is the one ``<dl class="guide">`` that is a direct child of its
  ``<figcaption>``; its fields are the direct ``dt``/``dd`` children, in order; the Panels ``dd``
  holds one ``<ul>`` with one ``<li>`` per ``Guide.panels`` item, in order;
- the guide's text nodes are separated by whitespace (spec example 10: "Question Which ...
  predicts? Panels One row per group"), so ``words(figcaption.text())`` reads the fields and the
  Panels items as separate words;
- the Data ``dd`` reads ``plots/data/<stem>.tsv; figure files: SVG PNG PDF`` (one label per image
  file, in the order ``_plot_slot`` builds them) or, without a data TSV, ``no data table for this
  figure; figure files: SVG PNG PDF`` (example 9);
- a short guide's How to read starts "As for the first" and links ``#guide-<pid>``; the V13 index
  instance (``V13_novel_orfs_locus``) has the short field set with its own index sentence and is
  not required to link the full guide (example 4);
- ``build_report`` reads ``content.GUIDES`` at call time (a monkeypatched ``GUIDES`` without an
  entry makes it raise ``KeyError``, § 6).
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from typing import Any

from helpers.m12_04_pr2 import words
from helpers.report_fixtures import Node, Page

FULL_LABELS = ("Question", "Panels", "How to read", "Caveats", "Data")
FULL_LABELS_NO_CAVEATS = ("Question", "Panels", "How to read", "Data")
SHORT_LABELS = ("Question", "How to read", "Data")
LETTERS = "abcdefghij"

#: the V13 index instance (M12-04 review r1 L2); its How to read is the index sentence
INDEX_STEMS: Mapping[str, str] = {"V13": "V13_novel_orfs_locus"}

#: § 4: the number of Panels items of each guide (lettered when more than one)
PANEL_COUNTS: Mapping[str, int] = {
    "V01": 1, "V02": 1, "V03": 1, "V05": 1, "V06": 3, "V07": 1, "V08": 1, "V10": 2, "V11": 3,
    "V12": 2, "V13": 2, "V14": 3, "V15": 4, "V16": 2, "V17": 5, "V18": 3,
}  # fmt: skip
#: § 4: the only guide without Caveats (example 2)
NO_CAVEATS = frozenset({"V17"})

#: § 4 V02 and V07, verbatim (AC 11): (question, panels, read, caveats)
V02_GUIDE: tuple[str, tuple[str, ...], str, str | None] = (
    "Is abundance clustered along the chromosome?",
    (
        "One panel per replicon: x is the position (kb), y the log2 abundance of each detected "
        "gene with a measured value (one point per gene); undetected genes are ticks along the "
        "bottom; the line is a LOESS trend; shaded bands are abundance islands; vertical lines "
        "mark the origin and terminus when they are known. With no measured abundance, each "
        "panel is one strip of ticks coloured by detection status.",
    ),
    "A shaded band is a run of neighbouring genes more (or less) abundant than the rest of the "
    "replicon, beyond chance at the family-wise error level in the title. The title gives the "
    "number of islands and the lag-1 Moran's I with its q. A wavy trend line without bands is not "
    "evidence of clustering.",
    "The trend line is a visual guide and is not tested. The island test shuffles whole "
    "transcription units, so operon partners looking alike do not by themselves make an island.",
)
V07_GUIDE: tuple[str, tuple[str, ...], str, str | None] = (
    "Which regions and categories are detected less than detectability predicts?",
    (
        "One row per group, in blocks: replicon, feature class, COG category, region; at most 40 "
        "rows, and a note names the rest. The point is observed minus expected detected genes in "
        "percentage points of the group's genes, the bar its 95 % confidence interval; the group "
        "size n is in the row label, and region rows also carry their q at the right.",
    ),
    "A row left of zero has fewer detected genes than its genes' detectability predicts, a row "
    "right of zero more. Region rows carry a q: only those at or below the level in the title are "
    "findings, and the title names the strongest. Other rows have a confidence interval only: "
    "one that excludes 0 is a departure.",
    "The expected counts come from the detectability model, so check its reliability in V18 (a) "
    "first: where the model is off, every group can drift.",
)

#: a ``<dl class="guide ...">`` element in the raw page (``dl`` never nests in a guide)
GUIDE_HTML = re.compile(r"<dl\b[^>]*\bclass=\"guide[\" ][^>]*>.*?</dl>", re.S)


def guides() -> dict[str, Any]:
    from gippa.report.content import GUIDES

    return dict(GUIDES)


def kids(node: Node) -> list[Node]:
    return [c for c in node.children if isinstance(c, Node)]


def guide_dls(figcaption: Node) -> list[Node]:
    """The ``dl.guide`` direct children of ``figcaption``."""
    return [c for c in kids(figcaption) if c.tag == "dl" and c.has_class("guide")]


def fields(dl: Node) -> list[tuple[str, Node]]:
    """``(dt text, dd)`` of a guide, in order; every ``dt`` is followed by its ``dd``."""
    out: list[tuple[str, Node]] = []
    label: str | None = None
    for c in kids(dl):
        if c.tag == "dt":
            assert label is None, f"two dt in a row: {label!r}, {words(c.text())!r}"
            label = words(c.text())
        elif c.tag == "dd":
            assert label is not None, f"a dd without its dt: {words(c.text())!r}"
            out.append((label, c))
            label = None
        else:
            raise AssertionError(f"<{c.tag}> directly in a guide")
    assert label is None, f"dt {label!r} without its dd"
    return out


def field(dl: Node, label: str) -> Node:
    found = [dd for lab, dd in fields(dl) if lab == label]
    assert len(found) == 1, f"{len(found)} dd for {label!r}"
    return found[0]


def plot_of(node: Node) -> str | None:
    cur: Node | None = node
    while cur is not None:
        if "data-plot" in cur.attrs and cur.tag != "dl":
            return cur.attrs["data-plot"]
        cur = cur.parent
    return None


def figcaption_of(fig: Node) -> Node:
    found = fig.by_tag("figcaption")
    assert len(found) == 1, f"{len(found)} figcaptions in a figure"
    return found[0]


def slot_figures(page: Page, pid: str) -> list[Node]:
    return [f for f in page.root.by_tag("figure") if plot_of(f) == pid]


def guide_of(fig: Node) -> Node:
    cap = figcaption_of(fig)
    dls = guide_dls(cap)
    assert len(dls) == 1, f"{cap.attrs.get('data-stem')}: {len(dls)} figcaption > dl.guide"
    return dls[0]


def labels(dl: Node) -> tuple[str, ...]:
    return tuple(lab for lab, _ in fields(dl))


def full_labels(guide: Any) -> tuple[str, ...]:
    return FULL_LABELS if guide.caveats is not None else FULL_LABELS_NO_CAVEATS


def panel_items(dl: Node) -> list[str]:
    dd = field(dl, "Panels")
    uls = [c for c in kids(dd) if c.tag == "ul"]
    assert len(uls) == 1, f"Panels dd holds {len(uls)} ul"
    return [words(li.text()) for li in kids(uls[0]) if li.tag == "li"]


def links(node: Node) -> list[tuple[str, str]]:
    return [(words(a.text()), a.attrs.get("href", "")) for a in node.by_tag("a")]


def data_text(stem: str, inst: Any) -> str:
    """The Data ``dd`` text of an ``InstanceFiles`` (AC 5, example 9)."""
    labs = " ".join(rel.rsplit(".", 1)[1].upper() for rel in inst.images)
    head = f"plots/data/{stem}.tsv" if inst.data is not None else (
        "no data table for this figure")  # fmt: skip
    return f"{head}; figure files: {labs}"


def data_links(inst: Any) -> list[tuple[str, str]]:
    out = []
    if inst.data is not None:
        out.append((f"plots/data/{inst.stem}.tsv", f"../plots/data/{inst.stem}.tsv"))
    out += [(rel.rsplit(".", 1)[1].upper(), f"../plots/{rel}") for rel in inst.images]
    return out


def guide_html(html: str) -> list[str]:
    """The outer HTML of every ``dl.guide`` of a page, in order (AC 6)."""
    return GUIDE_HTML.findall(html)
