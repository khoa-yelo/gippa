# ruff: noqa: RUF001  (the report's wording is Unicode)
"""The static ``report/report.html`` (GIPPA_SPEC.md §10.3 bullets 1, 2, 4; M9-02 behaviour 2-6).

Pure: a ``RunView``, the plot files read by ``gippa.io.plot_files.read_plot_files`` and the skip
reasons in, one HTML string out. The 13 ``SECTIONS`` always exist; each shows its plots (SVG text
inlined, ids prefixed with the instance stem so no two SVGs share an id), then the first
``TOP_ROWS`` rows of its tables as written, with a link to the full TSV. The template inlines its
CSS and ``tables.js``; ``twins`` (instance stem -> ``(div id, figure JSON)``, M9-14) adds one Plotly
div per instance beside its static SVG and inlines plotly.js once in ``<head>`` (only when a twin
exists). Every link is relative (``../plots/...``, ``../tables/...``) and names a file that exists,
apart from the absolute links to the docs-site methods pages (links, never loads; M12-04 Q12).

M12-04 PR 2: a section that does not apply to the run (``not_applicable``) is listed once in
``#sec-not-applicable`` at the end instead of being rendered; an empty-state instance is a
``.notice``; a 0-row preview is "none" with its link; the genes preview is ``GENES_PREVIEW``.

M12-04 PR 13: Gene detection shows the pivoted views of ``detection_summary``
(``gippa.report.detection``) instead of its raw preview; ``redact`` rewrites the provenance
command line (``content.redact_command_line``); a page whose fully inlined size reaches
``REPORT_LINK_THRESHOLD`` is built again with the V13 locus views linked (``plots/`` SVGs).
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import jinja2
import pandas as pd
from markupsafe import Markup

import gippa
from gippa.config import GippaConfig
from gippa.determinants.models import MIN_ROWS, MIN_TUS
from gippa.io.plot_files import PlotFiles, svg_title
from gippa.io.run_view import PATHS, RunView
from gippa.io.tsv import format_value
from gippa.report import content, detection, head
from gippa.viz import style
from gippa.viz.words import plain

TOP_ROWS = 20
#: bytes: at or over this fully inlined page size the V13 locus views are linked, not inlined;
#: 80 % of ``gippa.steps.plots.REPORT_BUDGET`` (25 MB, §10.3; M12-04 PR 13). Read at call time.
REPORT_LINK_THRESHOLD = 20_000_000
#: the instances linked when the page is too large
LINKED_PREFIX = "V13_novel_orfs_locus_"
UNESTIMATED = "unestimated"
#: a Plotly twin: ``(div id, figure JSON)``
Twin = tuple[str, str]
_BANNER = re.compile(r"plotly\.js v(\S+)")


@dataclass(frozen=True)
class Section:
    id: str
    heading: str
    plots: tuple[str, ...] = ()
    tables: tuple[str, ...] = ()


SECTIONS: tuple[Section, ...] = (
    Section("summary", "Summary"),
    Section("inputs-qc", "Inputs and quality control", ("V17",)),
    Section("mapping", "Mapping", ("V05", "V06")),
    Section("detection", "Gene detection", ("V07", "V18"), ("genes_tsv",)),
    Section(
        "genome-context",
        "Genome context",
        ("V01", "V02", "V10", "V03"),
        ("islands", "regions_stats", "autocorrelation"),
    ),
    Section("operons", "Operons", ("V11",), ("pairs", "operon_flags")),
    Section("variants", "Variants", ("V12",)),
    Section("novel-orfs", "Novel ORFs", ("V13",), ("novel_orfs",)),
    Section("function", "Function", ("V08", "V14")),
    Section("determinants", "Determinants of abundance", ("V15",)),
    Section("sensitivity", "Sensitivity", (), ("sensitivity",)),
    Section("methods", "Methods"),
    Section("provenance", "Provenance"),
)
#: the curated genes preview, in order (R-06); ``tables/genes.tsv`` keeps every column
GENES_PREVIEW = ("locus_tag", "gene", "product", "detection_status", "n_peptides_unique",
                 "abundance_log2", "region_ids", "flags")  # fmt: skip
NOT_APPLICABLE = "Not applicable in this run"
SENSITIVITY_NONE = "No flagged results, so there is nothing to compare across sensitivity reruns."
#: section id -> the tables whose absence makes it inapplicable (all of them absent)
_NEEDS: Mapping[str, tuple[str, ...]] = {
    "operons": ("pairs", "operon_flags"),
    "variants": ("variants",),
    "novel-orfs": ("novel_orfs",),
}
#: section id -> the calling that needs peptides (named in its not-applicable reason)
_NEEDS_PEPTIDES: Mapping[str, str] = {"variants": "variant calling",
                                      "novel-orfs": "novel ORF calling"}  # fmt: skip
#: section id -> the module it needs
_MODULE: Mapping[str, str] = {"function": "function"}

_ENV = jinja2.Environment(
    loader=jinja2.PackageLoader("gippa.report", "templates"),
    autoescape=jinja2.select_autoescape(["html", "j2"]),
    keep_trailing_newline=True,
    trim_blocks=True,
    lstrip_blocks=True,
)
_ID = re.compile(r'(\sid=")([^"]*)"')
_HREF = re.compile(r'(href=")#([^"]*)"')
_URL = re.compile(r"url\(\s*#([^)\s]+)\s*\)")


def embed_svg(text: str, prefix: str) -> str:
    """``text`` from its ``<svg`` on (XML prolog and doctype dropped), every id and local
    reference (``href="#x"``, ``url(#x)``) prefixed with ``<prefix>-``."""
    start = text.find("<svg")
    body = text[start:] if start >= 0 else ""
    body = _ID.sub(lambda m: f'{m.group(1)}{prefix}-{m.group(2)}"', body)
    body = _HREF.sub(lambda m: f'{m.group(1)}#{prefix}-{m.group(2)}"', body)
    return _URL.sub(lambda m: f"url(#{prefix}-{m.group(1)})", body)


def _script_json(text: str) -> str:
    """JSON that cannot end its ``<script>`` element (``<`` written as ``\\u003c``)."""
    return text.replace("<", "\\u003c")


def not_applicable(
    view: RunView, plots: Mapping[str, PlotFiles] | None = None
) -> list[tuple[str, str, str]]:
    """``(section id, name, plain reason)`` of each section that does not apply to the run (its
    module off and none of its plots drawn, or every table it needs absent), then the BGC content
    when ``bgc`` is off; the mobile module is a Genome context note instead (behaviour 4, lead
    decisions 1-3; review r1 M1: a drawn plot keeps a module's section)."""
    off = view.modules_off()
    plots = plots or {}
    out = []
    for s in SECTIONS:
        mod = _MODULE.get(s.id)
        needs = _NEEDS.get(s.id)
        drawn = any(not i.empty for p in s.plots if p in plots for i in plots[p].instances)
        if mod is not None and mod in off and not drawn:
            out.append((s.id, s.heading, plain(f"module_off:{mod}")))
        elif needs and not any(view.has(k) for k in needs):
            reason = plain(f"input_absent:{needs[0]}")
            if s.id in _NEEDS_PEPTIDES and not _has_rows(view, "peptide_hits"):
                reason += f" ({_NEEDS_PEPTIDES[s.id]} needs peptides)"
            out.append((s.id, s.heading, reason))
    if "bgc" in off:
        out.append(("", "Biosynthetic gene clusters (BGC)", "not searched, so no BGC regions "
                    "or BGC locus views (" + plain("module_off:bgc") + ")"))  # fmt: skip
    return out


def _has_rows(view: RunView, key: str) -> bool:
    try:
        frame = view.table(key) if view.has(key) else None
    except Exception:
        return False
    return frame is not None and len(frame) > 0


def _plot_slot(
    pid: str,
    files: PlotFiles | None,
    skipped: Mapping[str, str],
    twins: Mapping[str, Twin],
    notices: bool = False,
    linked: bool = False,
) -> dict[str, Any]:
    """One plot's slot; with ``notices`` an empty-state instance is a ``notice`` (its title, no
    image, no twin; M12-04 PR 2) and each drawn instance carries its ``caption``; with ``linked``
    a ``LINKED_PREFIX`` instance with an SVG file is an ``img`` of it, not inlined (PR 13)."""
    if files is None or not files.instances:
        note = f"not shown: {plain(skipped[pid])}" if pid in skipped else "plot file absent"
        return {"pid": pid, "note": note, "instances": [], "notices": []}
    insts: list[dict[str, Any]] = []
    empty: list[dict[str, Any]] = []
    for inst in files.instances:
        links = [(rel.rsplit(".", 1)[1].upper(), f"../plots/{rel}") for rel in inst.images]
        if inst.data is not None:
            links.append(("data TSV", f"../plots/{inst.data}"))
        if notices and inst.empty:  # the title and its file links, no image (review r1 M2)
            title = None if inst.svg is None else svg_title(inst.svg)
            empty.append({"title": " ".join((title or "Nothing to show in this run.").split()),
                          "links": links})  # fmt: skip
            continue
        img = next((f"../plots/{r}" for r in inst.images if r.endswith(".svg")), None)
        if not (linked and inst.stem.startswith(LINKED_PREFIX)) or inst.svg is None:
            img = None
        svg = None if inst.svg is None or img else Markup(embed_svg(inst.svg, inst.stem))
        twin = twins.get(inst.stem)
        spec = None if twin is None else {"id": twin[0], "spec": Markup(_script_json(twin[1]))}
        insts.append({"stem": inst.stem, "svg": svg, "img": img, "links": links, "twin": spec,
                      "caption": content.caption(pid, inst.stem)})  # fmt: skip
    drawn = empty or any(i["svg"] is not None or i["img"] for i in insts)
    return {"pid": pid, "note": None if drawn else "plot file absent", "instances": insts,
            "notices": empty}  # fmt: skip


def _genes_preview(view: RunView) -> tuple[list[str], list[list[str]]]:
    """The first ``TOP_ROWS`` genes in ``GENES_PREVIEW`` columns (those present): by
    ``abundance_log2`` descending, NA last, ties by ``locus_tag``; with no abundance at all, by
    ``detected``, then ``n_peptides_unique`` (both descending), then ``locus_tag`` (R-06)."""
    genes = view.table("genes_tsv")
    if genes is None:
        return [], []
    cols = [c for c in GENES_PREVIEW if c in genes.columns]

    def num(col: str) -> pd.Series:
        if col not in genes.columns:
            return pd.Series(float("nan"), index=genes.index)
        return pd.to_numeric(genes[col], errors="coerce").astype(float)

    ab = num("abundance_log2")
    tag = (
        genes["locus_tag"].astype(str) if "locus_tag" in genes.columns else genes.index.astype(str)
    )
    if bool(ab.notna().any()):
        keys = list(zip(ab.isna(), -ab.fillna(0.0), tag, strict=True))
    else:
        det, pep = num("detected").fillna(0.0), num("n_peptides_unique").fillna(0.0)
        keys = list(zip(det <= 0, -pep, tag, strict=True))
    order = sorted(range(len(genes)), key=lambda i: keys[i])[:TOP_ROWS]
    recs = genes[cols].iloc[order].to_dict("records")
    return cols, [[format_value(r[c]) for c in cols] for r in recs]


def _table(view: RunView, key: str) -> dict[str, Any]:
    name = PATHS[key].rsplit("/", 1)[1]
    try:
        head = view.head(key, TOP_ROWS)
        whole = view.table(key) if key == "sensitivity" and head is not None else None
        if key == "genes_tsv" and head is not None and head[1]:
            head = _genes_preview(view)
    except Exception as exc:  # an unreadable file is a note, never a crash (review B5)
        return {"key": key, "name": name, "error": _unreadable(key, exc)}
    if head is None:
        return {"key": key, "name": name, "absent": True}
    header, rows = head
    found: dict[str, Any] = {
        "key": key,
        "name": name,
        "header": header,
        "rows": rows,
        "href": f"../{PATHS[key]}",
    }
    if whole is not None and len(whole):  # rendered whole, by family (M9-16 behaviour 5)
        if _all_zero(whole, ("n_primary", "n_variant")):  # never on NA (review r1 L4)
            return found | {"sentence": SENSITIVITY_NONE}
        found |= {"header": list(whole.columns), "groups": content.sensitivity_groups(whole)}
    elif not rows:  # a 0-row preview is hidden (AC 2.6)
        found["none"] = True
    return found


def _detection(view: RunView) -> dict[str, Any] | None:
    """The Gene detection views and the link to the long file; ``None`` when it is absent; an
    unreadable file is a note, never a crash (review B5)."""
    if not view.has("detection_summary"):
        return None
    try:
        frame = view.table("detection_summary")
        found: dict[str, Any] = {} if frame is None else dict(detection.views(frame))
    except Exception as exc:
        return {"error": _unreadable("detection_summary", exc)}
    return found | {"href": f"../{PATHS['detection_summary']}"}


def _all_zero(frame: pd.DataFrame, cols: tuple[str, ...]) -> bool:
    """Whether every ``cols`` cell of ``frame`` is a number equal to 0 (an NA or a missing column
    is not)."""
    if not set(cols) <= set(frame.columns):
        return False
    counts = frame[list(cols)].apply(pd.to_numeric, errors="coerce")
    return bool((counts.notna() & (counts == 0)).all(axis=None))


def _unreadable(key: str, exc: BaseException) -> str:
    return f"{PATHS[key].rsplit('/', 1)[1]}: could not be read: {type(exc).__name__}: {exc}"


def _count(frame: pd.DataFrame, column: str, alpha: float | None = None) -> int:
    col = frame[column]
    if alpha is None:
        return int(col.fillna(False).astype(bool).sum())
    return int((pd.to_numeric(col, errors="coerce") <= alpha).sum())


def _tier_fdr(ds: pd.DataFrame | None, alpha: float) -> str | None:
    if ds is None:
        return None
    hit = ds[(ds["section"] == "tier") & (ds["group_type"] == "q_threshold")
             & (ds["group"] == format_value(alpha)) & (ds["metric"] == "fdr")]  # fmt: skip
    vals = pd.to_numeric(hit["value"], errors="coerce").dropna()
    return None if vals.empty else style.fmt_q(float(vals.iloc[0]))


def _novel_fdr(novel: pd.DataFrame) -> str | None:
    """The max class FDR over CALLED rows; ``None`` (unestimated) if any called row is (B1).
    M7-17: a suffix counts the called rows flagged ``denovo_only_unsupported`` (ADR-0050)."""
    called = novel[novel["called"].fillna(False).astype(bool)]
    fdrs = pd.to_numeric(called["class_fdr"], errors="coerce")
    flags = called["flags"].fillna("").astype(str).str.split(";")
    k = int(flags.map(lambda f: "denovo_only_unsupported" in f).sum())
    suffix = (f"; {style.fmt_count(k)} of {style.fmt_count(len(called))} de novo-only, no database"
              " or RBS support") if k else ""  # fmt: skip
    if called.empty or bool(fdrs.isna().any()):
        return UNESTIMATED + suffix if k else None
    return f"class FDR {_at_most(float(fdrs.max()))}{suffix}"


def _at_most(q: float) -> str:
    """``"≤ 0.87"``; ``"< 0.001"`` for a q below 0.001 (``fmt_q``)."""
    text = style.fmt_q(q)
    return text if text.startswith("<") else f"≤ {text}"


#: tier -> the genes card estimate; ``{fdr}`` and ``{q}`` are filled in (M12-04 PR 1)
_TIER_FDR = {
    "decoy_picked": "picked target–decoy gene-level FDR {fdr} at q ≤ {q}",
    "mayu": "Mayu gene-level FDR {fdr}, from peptides at q ≤ {q}",
}
_NO_FDR = {
    "two_peptide": "two-peptide rule: no FDR estimated",
    "protein_list": "protein list: no FDR estimated",
}
#: tiers where the two-peptide count is a separate statement
_TWO_PEP_LINE = ("decoy_picked", "mayu")
NO_POWER = "no predictive power"
#: why a card shows no ``.estimate`` (``checks.unlabelled_numbers`` accepts these)
EXACT_COUNT = "an exact count"
SD_IN_VALUE = "the CV SD is shown with the value"


def _alt_estimate(alt: pd.DataFrame) -> str | None:
    """The ``class_error_rate`` range over the numeric rows, "(N unestimated)" appended (lead
    decision 5); ``None`` when no row has a number."""
    rates = pd.to_numeric(alt["class_error_rate"], errors="coerce")
    known, n_un = rates.dropna(), int(rates.isna().sum())
    if known.empty:
        return None
    lo, hi = style.fmt_q(float(known.min())), style.fmt_q(float(known.max()))
    tail = f" ({style.fmt_count(n_un)} unestimated)" if n_un else ""
    return f"class error rate {lo if lo == hi else f'{lo}–{hi}'} (pooled per class){tail}"


def _metric(en: pd.DataFrame, metric: str) -> tuple[float | None, list[str]]:
    rows = en[en["metric"] == metric]
    vals = pd.to_numeric(rows["value"], errors="coerce").dropna()
    flags = [f for v in rows["flags"] if isinstance(v, str) for f in v.split(";")]
    return (None if vals.empty else float(vals.iloc[0])), flags


def _not_fitted(view: RunView, genes: pd.DataFrame | None) -> str:
    """Why the determinants model was not fitted: no abundance, or the minimum sizes with this
    run's counts."""
    counts = head.abundance_counts(genes)
    if "layer_missing:abundance" in view.run_flags() or (counts is not None and counts[0] == 0):
        return "no abundance in this run"
    need = f"needs ≥ {MIN_ROWS} genes with abundance in ≥ {MIN_TUS} TUs"
    if counts is None:
        return need
    return (
        f"{need} (this run: {style.fmt_count(counts[0])} genes, {style.fmt_count(counts[1])} TUs)"
    )


def cards(view: RunView, cfg: GippaConfig) -> tuple[list[dict[str, Any]], list[str]]:
    """The summary cards (M9-02 behaviour 5, M9-16 behaviour 1): a number and its error estimate
    or ``unestimated``, and under ``sources`` the ``<path>:<column>`` entries read (decision 3);
    a card whose table is absent is left out. Also the notes naming each table that could not be
    read (its card left out, the other cards kept; review B5)."""
    a = cfg.stats.alpha
    out: list[dict[str, Any]] = []
    notes: dict[str, str] = {}

    def table(key: str) -> pd.DataFrame | None:
        try:
            return view.table(key)
        except Exception as exc:
            notes.setdefault(key, _unreadable(key, exc))
            return None

    def src(key: str, *cols: str) -> list[str]:
        return [f"{PATHS[key]}:{c}" for c in cols] if view.has(key) else []

    def add(name: str, label: str, value: int | str, estimate: str | None, sources: list[str],
            *, no_estimate: str | None = None, **extra: str) -> None:  # fmt: skip
        """A card; a missing ``estimate`` reads ``unestimated``, or is left out when
        ``no_estimate`` says why none applies (an exact count, an SD in the value)."""
        card: dict[str, Any] = {"name": name, "label": label, "value": str(value),
                                "estimate": estimate or ("" if no_estimate else UNESTIMATED),
                                "sources": sources, **extra}  # fmt: skip
        if no_estimate and not estimate:
            card["no_estimate"] = no_estimate
        out.append(card)

    def n_unmapped(kind: str) -> int:
        unmapped = table("unmapped")
        return 0 if unmapped is None else int((unmapped["kind"] == kind).sum())

    roles = head.manifest_roles(view)  # gate on the manifest record only
    count = style.fmt_count

    def genes(t: pd.DataFrame) -> None:
        if t.empty:
            return  # no genes: no card, never "0 of 0" (M12-04 PR 1 review r1)
        tier = head.tier_of(view, t)  # the box's tier: manifest first
        sources = src("genes_tsv", "detected", "detection_tier", "detected_2pep")
        label = f"of {count(len(t))} genes detected"
        n = _count(t, "detected")
        est = _NO_FDR.get(str(tier))
        fdr = None if tier not in _TIER_FDR else _tier_fdr(table("detection_summary"),
                                                           a.identification)  # fmt: skip
        if fdr is not None:
            sources += src("detection_summary", "value")
            est = _TIER_FDR[str(tier)].format(fdr=fdr, q=style.fmt_alpha(a.identification))
        extra = {}
        if tier in _TWO_PEP_LINE:
            extra["detail"] = f"{count(_count(t, 'detected_2pep'))} with ≥ 2 distinct peptides"
        add("genes_detected", label, count(n), est, sources, **extra)

    def peptides(hits: pd.DataFrame) -> None:
        n_in = hits["peptide"].nunique() + n_unmapped("peptide")
        if (roles is not None and "peptides" not in roles) or n_in == 0:
            return  # no peptide input: no card, never "0 of 0" (M12-04 PR 1)
        placed = hits[pd.to_numeric(hits["placement_q"], errors="coerce") <= a.placement]
        sources = src("peptide_hits", "peptide", "placement_q") + src("unmapped", "kind")
        add("peptides_placed", f"of {count(n_in)} peptides placed",
            count(placed["peptide"].nunique()),
            f"placement q ≤ {style.fmt_alpha(a.placement)}", sources)  # fmt: skip

    def proteins(hits: pd.DataFrame) -> None:
        n = hits["protein_id"].nunique()
        n_in = n + n_unmapped("protein")
        if (roles is not None and "proteins" not in roles) or n_in == 0:
            return
        add("proteins_placed", f"of {count(n_in)} proteins placed", count(n), None,
            src("protein_hits", "protein_id") + src("unmapped", "kind"),
            no_estimate=EXACT_COUNT)  # fmt: skip

    def islands(t: pd.DataFrame) -> None:
        add("islands", "abundance islands", _count(t, "p_value", a.context),
            f"FWER ≤ {format_value(a.context)}", src("islands", "p_value"))  # fmt: skip

    def operons(t: pd.DataFrame) -> None:
        q = format_value(a.operon_flags)
        add(
            "operon_flags",
            f"operon pairs flagged at q ≤ {q} (discordant or co-abundant)",
            _count(t, "q_value", a.operon_flags),
            f"FDR ≤ {q}",
            src("operon_flags", "q_value"),
        )

    def variants(t: pd.DataFrame) -> None:
        alt = t[t["call"] == "ALT"]  # ``mixed`` is not an ALT call (lead decision 5)
        counts = alt["confidence"].astype(str).value_counts()
        known = [c for c in ("high", "candidate") if c in counts.index]
        order = known + sorted(c for c in counts.index if c not in known)
        by = ", ".join(f"{count(int(counts[c]))} {c}" for c in order)
        add("alt_calls", f"ALT calls ({by or 'none'})", count(len(alt)), _alt_estimate(alt),
            src("variants", "call", "confidence", "class_error_rate"))  # fmt: skip

    def novel(t: pd.DataFrame) -> None:
        add("novel_orfs", "novel ORFs called", count(_count(t, "called")), _novel_fdr(t),
            src("novel_orfs", "called", "class_fdr"))  # fmt: skip

    def model(t: pd.DataFrame) -> None:
        """ADR-0033: the elastic-net whole-model rows only; ADR-0034: "(uncalibrated)"; label
        first, "no predictive power" when the mean or mean minus SD is ≤ 0 (M12-04 PR 1)."""
        en = t[(t["model"] == "elastic_net") & (t["feature_or_block"] == "all")]
        (r2, _), (sd, _) = _metric(en, "cv_r2_mean"), _metric(en, "cv_r2_sd")
        p, flags = _metric(en, "null_p")
        sources = src("determinants_model", "model", "feature_or_block", "metric", "value",
                      "flags")  # fmt: skip
        if r2 is None:
            add("model_r2", "", "not fitted:", _not_fitted(view, table("genes_tsv")), sources,
                prefix="Model R²:")  # fmt: skip
            return
        value = style.fmt_num(r2) if sd is None else head.mean_sd(r2, sd)
        no_power = r2 <= 0 or (sd is not None and r2 - sd <= 0)
        label = f"(CV), {NO_POWER}" if no_power else "(CV)"
        unc = " (uncalibrated)" if "null_uncalibrated" in flags else ""
        est = None if p is None else f"null p {style.fmt_q(p)}{unc}"
        add("model_r2", label if est is None else f"{label};", value, est, sources,
            no_estimate=None if sd is None else SD_IN_VALUE, prefix="Model R²:")  # fmt: skip

    specs = [("genes_tsv", genes), ("peptide_hits", peptides), ("protein_hits", proteins),
             ("islands", islands), ("operon_flags", operons), ("variants", variants),
             ("novel_orfs", novel), ("determinants_model", model)]  # fmt: skip
    for key, card in specs:
        frame = table(key)
        if frame is None:
            continue
        try:
            card(frame)
        except Exception as exc:
            notes.setdefault(key, _unreadable(key, exc))
    return out, list(notes.values())


def build_report(
    view: RunView,
    plots: Mapping[str, PlotFiles],
    skipped: Mapping[str, str],
    cfg: GippaConfig,
    twins: Mapping[str, Twin] | None = None,
    redact: bool = False,
) -> str:
    """The report page (behaviour 2-6); deterministic, offline, no absolute path. ``redact``
    rewrites the provenance command line (R-09). At or over ``REPORT_LINK_THRESHOLD`` UTF-8 bytes
    fully inlined, the page is built again with the V13 locus views linked (PR 13)."""
    page = _page(view, plots, skipped, cfg, twins or {}, redact, linked=False)
    if len(page.encode("utf-8")) < REPORT_LINK_THRESHOLD:
        return page
    return _page(view, plots, skipped, cfg, twins or {}, redact, linked=True)


def _page(
    view: RunView,
    plots: Mapping[str, PlotFiles],
    skipped: Mapping[str, str],
    cfg: GippaConfig,
    twins: Mapping[str, Twin],
    redact: bool,
    linked: bool,
) -> str:
    absent = not_applicable(view, plots)
    off = {sid for sid, _, _ in absent}
    shown_sections = [s for s in SECTIONS if s.id not in off]
    sections: list[dict[str, Any]] = [
        {"id": s.id, "heading": s.heading,
         "plots": [_plot_slot(p, plots.get(p), skipped, twins, True, linked)
                   for p in s.plots],
         "tables": [_table(view, k) for k in s.tables]}
        for s in shown_sections
    ]  # fmt: skip
    summary, notes = cards(view, cfg)
    doc, note = content.read_manifest(view)
    det = _detection(view)
    if redact:  # the unreadable-table notes, paths rewritten (review r1)
        outdir = str((doc or {}).get("outdir") or "")
        notes = [content.redact_note(n, outdir) for n in notes]
        for t in (t for s in sections for t in s["tables"]):
            if t.get("error"):
                t["error"] = content.redact_note(t["error"], outdir)
        if det and det.get("error"):
            det["error"] = content.redact_note(det["error"], outdir)
    shown = any(i.stem in twins and not i.empty for s in shown_sections for p in s.plots
                if p in plots for i in plots[p].instances)  # fmt: skip
    entries = [(s.id, s.heading) for s in shown_sections]
    if absent:
        entries.append(("not-applicable", NOT_APPLICABLE))
    plotlyjs = _plotlyjs() if shown else None
    version = None if plotlyjs is None else _BANNER.search(plotlyjs)
    table = cfg.annotation.translation_table
    return _ENV.get_template("report.html.j2").render(
        sections=sections, cards=summary, notes=notes, version=gippa.__version__,
        box=head.run_box(view, plots, skipped, cfg),
        toc=head.toc(view, entries), not_applicable=[(n, r) for _, n, r in absent],
        not_applicable_heading=NOT_APPLICABLE, glossary=content.GLOSSARY,
        detection=det,
        provenance=content.provenance(view, doc, note, redact), extras=content.extras(view, doc),
        warnings=content.summary_warnings(view), snippets=content.snippets(),
        code_note=content.genetic_code_note(table),
        code_methods=content.genetic_code_methods(table),
        licence_notes=content.licence_notes(doc),
        demo_banner=content.demo_banner(doc),
        plotlyjs=plotlyjs and Markup(plotlyjs),
        plotly_version=version and version.group(1),
    )  # fmt: skip


def _plotlyjs() -> str:
    """plotly.js as bundled in the plotly package (no network)."""
    import plotly.offline

    text: str = plotly.offline.get_plotlyjs()
    return text
