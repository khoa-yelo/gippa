# ruff: noqa: RUF001, E501, E731, B023
"""Decision prototype of a simplified report.html (Khoa's request; targets 1.0.1). NOT shipped code.

    PYTHONPATH=src python scripts/dev/report_prototype.py OUT out.html [--titles]

Reads a finished run through the same readers and twin builders as ``gippa report``, then lays the
page out differently: key points first, one finding per section, grouped figures split into one
large plot per panel, everything else in collapsed detail. It computes no statistic: every number
is read from the run's tables, cards or finding titles. Long decimals are rounded for display only;
the full value stays in the cell's tooltip. Every section and main plot also carries a static "How
to read this" and a "Numbers" block computed from the figure's data TSV; a number that cannot be
read prints "not available" and is listed on stderr. Output links are absolute paths into OUT, so the page
is a local artefact: do not commit one built from pilot data.
"""

from __future__ import annotations

import html
import json
import re
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import plotly.offline
import yaml

from gippa import __version__
from gippa.config import GippaConfig
from gippa.io.plot_files import read_plot_files, svg_title
from gippa.io.run_view import PATHS, RunView
from gippa.report import content
from gippa.report.build import SECTIONS, cards, embed_svg
from gippa.steps.plots import load_plot_modules, twin_meta, twins

TICK, AXIS, LEGEND, TITLE = 12, 13, 12, 14

#: how each grouped figure is split: panels matching ``show`` stay visible, the rest are collapsed.
#: "first" = the first panel (or the first replicon's panels); ``None`` = every panel.
PANELS: dict[str, dict[str, Any]] = {
    "V17": {"show": ["abundance", "peptide lengths"], "h": 400, "more_cols": 3, "more_h": 420},
    "V05": {"show": None, "cols": 1, "h": 500},
    "V06": {"show": None, "h": 430},
    "V07": {"show": None, "cols": 1, "h": 1040},
    "V18": {"show": ["reliability"], "h": 430, "more_cols": 3, "more_h": 270},
    "V02": {"show": "first", "cols": 1, "h": 430, "more_cols": 2},
    "V10": {"show": "first_replicon", "h": 420, "more_cols": 2},
    "V11": {"show": ["pair class", "partners"], "h": 440},
    "V08": {"show": None, "cols": 1, "h": 520, "top": 34, "legend": False},
    "V14": {"show": ["NES", "running"], "h": 520, "more_cols": 1, "more_h": 560},
    "V15": {"show": ["block", "observed"], "h": 450, "more_cols": 2, "more_h": 360},
}
#: table -> (sort column, ascending, key columns). Other columns stay behind "all columns".
TABLES: dict[str, tuple[str | None, bool, list[str]]] = {
    "genes_tsv": ("abundance_log2", False, ["locus_tag", "gene", "product", "detected", "detection_tier", "gene_q_value", "n_peptides_unique", "coverage_frac", "abundance_log2", "p_detect", "cog_category"]),
    "detection_summary": (None, True, ["section", "group_type", "group", "metric", "value", "ci_low", "ci_high", "n"]),
    "regions_stats": ("q_value", True, ["region_id", "region_type", "stat", "n_genes", "observed", "expected", "effect", "ci_low", "ci_high", "q_value"]),
    "islands": ("q_value", True, ["island_id", "direction", "contig", "start", "end", "n_genes", "mean_abundance_log2", "detected_frac", "q_value"]),
    "autocorrelation": ("q_value", True, ["contig", "statistic", "lag", "null_model", "value", "null_lo", "null_hi", "q_value"]),
    "pairs": ("q_value", True, ["gene_a", "gene_b", "class", "test", "intergenic_bp", "delta_log2", "q_value", "lfdr"]),
    "operon_flags": ("q_value", True, ["gene_a", "gene_b", "uniop_class", "abundance_pattern", "flag", "q_value", "lfdr"]),
    "novel_orfs": ("rank", True, ["orf_id", "class", "contig", "strand", "length_aa_best", "n_peptides", "best_peptide_q", "ostir_tir", "class_fdr", "called", "overlapping_gene"]),
    "variants": ("n_alt", False, ["gene", "aa_pos", "ref_aa", "alt_aa", "alt_type", "call", "confidence", "n_ref", "n_alt", "class_error_rate"]),
}
LABELS = {
    "o_minus_e": "observed − expected", "q_value": "q", "abundance_log2": "log2 abundance",
    "p_detect": "expected P(detect)", "gene_q_value": "gene q", "coverage_frac": "coverage",
    "n_peptides_unique": "unique peptides", "cog_category": "COG", "delta_log2": "Δ log2",
    "n_genes": "genes", "ci_low": "CI low", "ci_high": "CI high", "lfdr": "local fdr",
    "n_ref": "ref peptides", "n_alt": "ALT peptides", "class_error_rate": "class error rate",
    "class_fdr": "class FDR", "best_peptide_q": "best peptide q", "ostir_tir": "OSTIR TIR",
    "length_aa_best": "length (aa)", "n_peptides": "peptides", "intergenic_bp": "gap (bp)",
    "null_lo": "null low", "null_hi": "null high", "mean_abundance_log2": "mean log2 abundance",
    "detected_frac": "detected fraction",
}
#: the sentence that opens each section: the finding title of this plot (generated from tables)
LEAD = {"inputs-qc": "V17", "mapping": "V05", "detection": "V07", "genome-context": "V01",
        "operons": "V11", "variants": "V12", "novel-orfs": "V13", "function": "V14",
        "determinants": "V15"}  # fmt: skip
NAV = {"summary": "Key points", "inputs-qc": "QC", "mapping": "Mapping", "detection": "Detection",
       "genome-context": "Genome", "operons": "Operons", "variants": "Variants",
       "novel-orfs": "Novel ORFs", "function": "Function", "determinants": "Determinants",
       "sensitivity": "Sensitivity", "methods": "Methods", "provenance": "Provenance"}  # fmt: skip
_LONG = re.compile(r"(?<![\w.])(\d+\.\d{4,}(?:[eE][+-]?\d+)?)(?![\w.])")
_LONGNUM = re.compile(r"^[+-]?\d*\.\d{4,}([eE][+-]?\d+)?$|^[+-]?\d+(\.\d+)?[eE][+-]?\d+$")
E = html.escape


def short(text: str) -> str:
    """Long decimals in a sentence rounded to 3 significant figures (display only)."""
    return _LONG.sub(lambda m: f"{float(m.group(1)):.3g}", text)


def cell(value: Any) -> str:
    raw = "NA" if value is None or (isinstance(value, float) and pd.isna(value)) else str(value)
    if re.fullmatch(r"-?\d+\.0", raw):
        raw = raw[:-2]
    if _LONGNUM.match(raw):
        return f'<td title="{E(raw)}">{float(raw):.3g}</td>'
    if len(raw) > 48:
        return f'<td title="{E(raw)}">{E(raw[:46])}…</td>'
    return f"<td>{E(raw)}</td>"


# --- figures: one Plotly figure per panel -----------------------------------------------------


def _suffix(ref: str) -> str:
    return ref[1:].split(" ")[0]


def explode(figure: dict[str, Any]) -> list[dict[str, Any]]:
    """One ``{title, data, layout}`` per subplot of a twin; blank panels are dropped."""
    lay, panels = figure["layout"], []
    keep = {k: lay[k] for k in KEEP_LAYOUT if k in lay}
    subtitles = [a for a in lay.get("annotations", []) if a.get("xref") == "paper" and a.get("yref") == "paper" and a.get("yanchor") == "bottom"]
    keys = sorted((k for k in lay if re.fullmatch(r"xaxis\d*", k)), key=lambda k: int(k[5:] or 1))
    for key in keys:
        xid = "x" + key[5:]
        x = dict(lay[key])
        yid = str(x.get("anchor", "y"))
        y = dict(lay.get("yaxis" + yid[1:], {}))
        traces = [dict(t) for t in figure["data"] if t.get("xaxis", "x") == xid and "domain" not in t]
        if not traces:
            continue
        centre, top = sum(x["domain"]) / 2 if "domain" in x else 0.5, (y.get("domain") or [0, 1])[1]
        title = ""
        for a in subtitles:
            if abs(a["x"] - centre) < 0.03 and abs(a["y"] - top) < 0.05:
                title = a["text"]
        for t in traces:
            t["xaxis"], t["yaxis"] = "x", "y"
        def local(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
            out = []
            for it in items:
                xr, yr = str(it.get("xref", "x")), str(it.get("yref", "y"))
                if xr.split(" ")[0] == xid or yr.split(" ")[0] == yid:
                    it = dict(it)
                    if xr.split(" ")[0] == xid:
                        it["xref"] = xr.replace(xid, "x", 1)
                    elif xr != "paper":
                        continue
                    if yr.split(" ")[0] == yid:
                        it["yref"] = yr.replace(yid, "y", 1)
                    elif yr != "paper":
                        continue
                    out.append(it)
            return out
        for ax in (x, y):
            for gone in ("domain", "anchor", "matches", "overlaying", "scaleanchor"):
                ax.pop(gone, None)
        panels.append({"title": title, "data": traces, "x": x, "y": y,
                       "annotations": local(lay.get("annotations", [])), "shapes": local(lay.get("shapes", [])), "keep": keep})  # fmt: skip
    for t in figure["data"]:
        if "domain" in t or t.get("type") in ("sankey", "table"):
            title = next((a["text"] for a in subtitles if t.get("type") == "table" and "warning" in a["text"]), "")
            panels.append({"title": title, "data": [{**t, "domain": {"x": [0, 1], "y": [0, 1]}}], "x": None, "y": None, "annotations": [], "shapes": [], "keep": {}})
    return panels


KEEP_LAYOUT = ("barmode", "barnorm", "bargap", "bargroupgap", "violinmode", "violingap", "boxmode", "hovermode")


def panel_figure(panel: dict[str, Any], template: Any, height: int, top: int = 12, legend: bool = True) -> str:
    """The panel as a figure at reading sizes (the twins use the print sizes, 6-7, as CSS px)."""
    named = any(t.get("showlegend", True) and t.get("name") for t in panel["data"])
    many = sum(1 for t in panel["data"] if t.get("showlegend", True) and t.get("name")) > 1
    for t in panel["data"]:
        if t.get("type") == "sankey":
            t["textfont"] = {"size": TICK}
    layout: dict[str, Any] = {
        "template": template, "height": height, "autosize": True,
        "font": {"size": TICK, "family": "DejaVu Sans, sans-serif"},
        "margin": {"l": 70, "r": 24, "t": top, "b": 60},
        "showlegend": bool(named and many and legend),
        **panel.get("keep", {}),
        "legend": {"orientation": "h", "y": -0.16, "yanchor": "top", "font": {"size": LEGEND}},
        "shapes": panel["shapes"],
        "annotations": [{**a, "font": {**a.get("font", {}), "size": max(a.get("font", {}).get("size", 0), 11)}} for a in panel["annotations"]],
    }
    if layout["showlegend"]:
        layout["margin"]["b"] = 110
    for name, ax in (("xaxis", panel["x"]), ("yaxis", panel["y"])):
        if ax is None:
            continue
        ax = {**ax, "automargin": True, "tickfont": {"size": TICK}}
        title = ax.get("title")
        if isinstance(title, dict):
            ax["title"] = {**title, "font": {"size": AXIS}}
        layout[name] = ax
    return json.dumps({"data": panel["data"], "layout": layout}, separators=(",", ":")).replace("<", "\\u003c")


class Page:
    def __init__(self) -> None:
        self.n = 0

    def plot(self, panel: dict[str, Any], template: Any, height: int, **opts: Any) -> str:
        self.n += 1
        spec = panel_figure(panel, template, height, **opts)
        title = re.sub(r"<br\s*/?>", " · ", panel["title"] or "")
        head = f'<div class="ptitle">{title}</div>' if title else ""
        return (f'<div class="panel">{head}<div class="pl" style="height:{height}px" data-spec="s{self.n}"></div>'
                f'<script type="application/json" id="s{self.n}">{spec}</script></div>')  # fmt: skip

    def grid(self, panels: list[dict[str, Any]], template: Any, cols: int, height: int, **opts: Any) -> str:
        body = "".join(self.plot(p, template, height, **opts) for p in panels)
        return f'<div class="grid c{cols}">{body}</div>' if body else ""


def links(inst: Any, out: Path) -> str:
    items = [(r.rsplit(".", 1)[1].upper(), r) for r in inst.images]
    if inst.data is not None:
        items.append(("data", inst.data))
    return "".join(f'<a href="file://{out}/plots/{r}">{label}</a>' for label, r in items)


def more(summary: str, body: str, open_: bool = False) -> str:
    if not body:
        return ""
    return f'<details class="more"{" open" if open_ else ""}><summary>{E(summary)}</summary>{body}</details>'


def split(pid: str, panels: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    show = PANELS[pid]["show"]
    if show is None:
        return panels, []
    if show == "first":
        return panels[:1], panels[1:]
    if show == "first_replicon":
        first = (panels[0]["title"] or "").split(" ")[-1]
        keep = [p for p in panels if first and (p["title"] or "").endswith(first)]
        return keep, [p for p in panels if p not in keep]
    keep = [p for p in panels if any(s in (p["title"] or "") for s in show)]
    return keep, [p for p in panels if p not in keep]


def twin_block(pid: str, inst: Any, spec: str, page: Page, out: Path) -> tuple[str, str]:
    conf = PANELS[pid]
    figure = json.loads(spec)
    template = figure["layout"].get("template")
    shown, rest = split(pid, explode(figure))
    cols, height = conf.get("cols", 2), conf["h"]
    opts = {"top": conf.get("top", 12), "legend": conf.get("legend", True)}
    visible = page.grid(shown, template, 1 if len(shown) == 1 else cols, height, **opts)
    hidden = ""
    if rest:
        tall = [p for p in rest if "SHAP" in (p["title"] or "")]
        hidden = page.grid(tall, template, 1, 560) + page.grid([p for p in rest if p not in tall], template, conf.get("more_cols", 2), conf.get("more_h", 340))
    label = f"{len(rest)} more panel{'s' if len(rest) != 1 else ''}"
    # figure-level notes (paper coordinates) are caveats, not panel furniture: keep them as text
    anns = figure["layout"].get("annotations", [])
    texts = [a["text"] for a in anns if a.get("xref") == "paper" and a.get("yref") == "paper" and a.get("yanchor") != "bottom"]
    dropped = [sh for sh in figure["layout"].get("shapes", []) if sh.get("xref") == "paper" and sh.get("yref") == "paper"]
    if dropped:
        print(f"{pid}: {len(dropped)} paper-level shapes not carried over", file=sys.stderr)
    notes = "".join(f'<p class="fnote">{t}</p>' for t in texts)
    tools = f'{notes}<div class="flinks">{pid} {links(inst, out)}</div>'
    return visible + (more(label, hidden) if hidden else "") + tools, inst.stem


def static_block(inst: Any, width: int, out: Path) -> str:
    svg = embed_svg(inst.svg, inst.stem)
    native = re.search(r'<svg[^>]*\swidth="([\d.]+)pt"', svg)
    if native:  # 1.8 x the print size: a single-column figure must not fill the page
        width = min(width, round(float(native.group(1)) * 4 / 3 * 1.8))
    return f'<figure class="static" style="max-width:{width}px">{svg}<div class="flinks">{E(inst.stem)} {links(inst, out)}</div></figure>'


def loci_by_rank(out: Path, instances: list[Any]) -> tuple[list[Any], list[Any]]:
    """Loci in the index's rank order (the current report shows them alphabetically); the index
    instance returned separately."""
    index = [i for i in instances if i.stem in ("V03_locus", "V13_novel_orfs_locus")]
    loci = [i for i in instances if i not in index]
    if index and index[0].data:
        tsv = pd.read_csv(out / "plots" / index[0].data, sep="\t", dtype=str, keep_default_na=False)
        if "locus_id" in tsv:
            rank = {f"{index[0].stem}_{v}": n for n, v in enumerate(tsv["locus_id"])}
            loci.sort(key=lambda i: rank.get(i.stem, len(rank)))
    return loci, index


# --- tables -----------------------------------------------------------------------------------


def table_block(view: RunView, key: str, top: int = 20) -> str:
    name = PATHS[key].rsplit("/", 1)[1]
    try:
        frame = view.table(key)
    except Exception as exc:
        return f'<p class="note">{E(name)}: could not be read: {E(type(exc).__name__)}</p>'
    if frame is None:
        return f'<p class="note">{E(name)}: not produced in this run.</p>'
    sort, asc, keep = TABLES[key]
    if not len(frame):
        return f'<p class="note">{E(name)}: no rows in this run.</p>'
    total = len(frame)
    if sort and sort in frame:
        frame = frame.assign(_s=pd.to_numeric(frame[sort], errors="coerce")).sort_values("_s", ascending=asc, na_position="last", kind="stable").drop(columns="_s")
    frame = frame.head(top)
    cols = [c for c in keep if c in frame] + [c for c in frame if c not in keep]
    head = "".join(f'<th class="{"" if c in keep else "x"}" title="{E(c)}">{E(LABELS.get(c, c.replace("_", " ")))}</th>' for c in cols)
    rows = []
    for rec in frame[cols].to_dict("records"):
        tds = [cell(rec[c]) for c in cols]
        rows.append("<tr>" + "".join(td if c in keep else td.replace("<td", '<td class="x"', 1) for td, c in zip(tds, cols, strict=True)) + "</tr>")
    order = f", sorted by {LABELS.get(sort, sort)}" if sort else ""
    n_x = len(cols) - len([c for c in keep if c in frame])
    toggle = f'<label class="xt"><input type="checkbox" data-allcols> all {len(cols)} columns</label>' if n_x else ""
    return (f'<div class="tbl"><div class="thead"><b>{E(name)}</b> <span>top {len(frame)} of {total} rows{order}</span>'
            f'<input type="search" data-filter="{key}" placeholder="filter rows" aria-label="filter {E(name)}">{toggle}'
            f'<a href="file://{view_root(view)}/{PATHS[key]}">full table</a></div>'
            f'<div class="table-wrap"><table data-table="{key}"><thead><tr>{head}</tr></thead><tbody>{"".join(rows)}</tbody></table></div></div>')  # fmt: skip


_OUT: list[Path] = []


def view_root(_: RunView) -> Path:
    return _OUT[0]


def extras_block(items: list[dict[str, Any]]) -> tuple[str, int]:
    body, n = "", 0
    for x in items:
        if "rows" in x:
            head = "".join(f"<th>{E(str(h))}</th>" for h in x["header"])
            rows = "".join("<tr>" + "".join(cell(v) for v in r) + "</tr>" for r in x["rows"])
            body += f'<p>{E(x["caption"])}:</p><div class="table-wrap"><table data-table="{E(x["key"])}"><thead><tr>{head}</tr></thead><tbody>{rows}</tbody></table></div>'
        else:
            body += f'<p class="caveat">{E(x["text"])}</p>'
            n += 1
    return body, n


# --- key points -------------------------------------------------------------------------------


def tile(c: dict[str, Any]) -> str:
    label = (c.get("prefix", "") + " " + c["label"]).strip().rstrip(";")
    value, est = short(c["value"]), c["estimate"]
    empty = c["name"] == "proteins_placed" and c["value"] == "0"
    brief = est.split("; descriptive")[0] if c["name"] == "peptides_placed" else est
    detail = f'<div class="td">{E(short(c["detail"]))}</div>' if c.get("detail") else ""
    return (f'<div class="tile{" dim" if empty else ""}" data-card="{E(c["name"])}" title="{E(est)}"><div class="tv">{E(value)}</div>'
            f'<div class="tl">{E(label)}</div><div class="te">{E(short(brief))}</div>{detail}</div>')  # fmt: skip


def titles_of(files: dict[str, Any]) -> dict[str, str]:
    out = {}
    for pid, pf in files.items():
        for inst in pf.instances:
            if inst.svg and pid not in out:
                t = svg_title(inst.svg)
                if t:
                    out[pid] = t
    return out


def locus_titles(pf: Any) -> dict[str, str]:
    return {i.stem: svg_title(i.svg) or "" for i in pf.instances if i.svg}


def key_points(view: RunView, cfg: GippaConfig, files: dict[str, Any], titles: dict[str, str], summary: list[dict[str, Any]], notes: list[str], doc: Any) -> str:
    by = {c["name"]: c for c in summary}
    pre = view.preflight() or {}
    reps = pre.get("replicons", []) if isinstance(pre, dict) else []
    size = sum(int(r["length"]) for r in reps)
    who = f"{len(reps)} replicon{'s' if len(reps) != 1 else ''}, {size / 1e6:.2f} Mb" if reps else ""
    tier = pre.get("evidence_tier") if isinstance(pre, dict) else None
    meta = " · ".join(x for x in (who, f"evidence: {tier}" if isinstance(tier, str) else "", f"GIPPA {(doc or {}).get('gippa_version', __version__)} (run)") if x)

    finds: list[tuple[str, str, str]] = []
    def find(sec: str, label: str, pid: str) -> None:
        if pid in titles:
            finds.append((sec, label, titles[pid]))
    find("genome-context", "Overall", "V01")
    find("mapping", "Mapping", "V05")
    find("detection", "Detected vs expected", "V07")
    if "V03" in files:
        loci, _ = loci_by_rank(_OUT[0], files["V03"].instances)
        tops = [svg_title(i.svg) for i in loci[:2] if i.svg]
        if tops:
            finds.append(("genome-context", "Top-ranked loci", "; ".join(t for t in tops if t)))
    find("function", "Function", "V14")
    find("determinants", "Determinants", "V15")
    find("detection", "Calibration", "V18")

    care: list[str] = []
    checks = pre.get("checks", []) if isinstance(pre, dict) else []
    # the report must not carry absolute paths (M9-02 behaviour 6): a shipped version redacts properly
    care += [f'Preflight {c.get("level")}: ' + re.sub(r"/(?:[\w.-]+/)+([\w.-]+)", r"\1", str(c.get("message") or c.get("rule_id"))) for c in checks if isinstance(c, dict) and c.get("level") in ("warn", "error")]
    for c in summary:
        if c["name"] in ("alt_calls", "novel_orfs") and c["value"] != "0" and "unestimated" in c["estimate"]:
            care.append(f'{c["value"]} {c["label"]}: error rate unestimated.'.replace("1 novel ORFs", "1 novel ORF"))
    if "V06" in titles and ("No stratum" in titles["V06"] or "unestimated" in titles["V06"]):
        care.append(f'Tolerant mapping: {titles["V06"]}.')
    r2 = by.get("model_r2")
    if r2 and "uncalibrated" in r2["estimate"]:
        care.append("Model R²: the permutation null p is uncalibrated, not a discovery claim (ADR-0034).")
    care += [x["text"] for x in content.extras(view, doc).get("summary", []) if "text" in x]
    care += notes
    code = content.genetic_code_note(cfg.annotation.translation_table)
    if code:
        care.append(code["text"])
    li = "".join(f'<li><a href="#{s}">{E(lab)}</a> {E(short(t))}</li>' for s, lab, t in finds)
    warn = "".join(f"<li>{E(short(t))}</li>" for t in care) or "<li>No warnings from this run.</li>"
    order = ["genes_detected", "peptides_placed", "model_r2", "alt_calls", "novel_orfs", "islands", "operon_flags", "proteins_placed"]
    tiles = "".join(tile(by[n]) for n in order if n in by)
    return (f'<p class="meta">{E(meta)}</p><div class="tiles">{tiles}</div>'
            f'<div class="kp"><div><h3>What this run shows</h3><ul>{li}</ul></div>'
            f'<div class="care"><h3>Read with care</h3><ul>{warn}</ul></div></div>')  # fmt: skip


# --- sections ---------------------------------------------------------------------------------


# --- "How to read this" and "Numbers in this plot" ----------------------------------------------
#
# The text is static and generic to any run; it paraphrases GIPPA_SPEC.md and docs/methods/ (the
# source of each entry is named beside it). Every number is read at build time from the figure's
# data TSV or from the run table the figure draws; none is typed here.

#: per-section text. The summary: §10.3 "Wording"; the others as named.
SECTION_READ: dict[str, str] = {
    "summary": "Each tile is one headline number with its error estimate below it, or “unestimated” when the run has none. “Uncalibrated” and “unestimated” mean the value describes this run and is not a validated error rate. Each line links to the section with its plot and caveats.",  # §10.3
    "inputs-qc": "These plots ask whether the inputs look sensible before any result is read. They describe the inputs and do not show that the results are correct.",  # §10.2 V17
    "mapping": "This section follows every input sequence to where it ended up, then shows how reliable the placements are. A placement says where a peptide matches the genome; gene detection is decided later, in the next section.",  # §10.2 V5, V6; §7.3
    "detection": "This section compares how many genes were detected with how many their physical properties predict. “Not detected” is not “not expressed”: an undetected gene may lie below the detection limit.",  # §10.3 wording; docs/methods/islands.md (left-censored); detectability.md
    "genome-context": "This section asks where along the genome detected proteins lie and whether abundance clusters beyond operons. Failing to find a cluster here is not evidence that none exists.",  # §10.2 V1, V2, V10; docs/methods/islands.md
    "operons": "This section compares the abundance of neighbouring genes. Both flag families rarely fire: an empty list at the q cut-off is the expected output, and is evidence neither of concordance nor against missed operons.",  # docs/methods/operons.md, Discordance and Power
    "variants": "A site is ALT when only ALT peptides cover it and mixed when both REF and ALT peptides do; confidence is high with at least 2 distinct ALT peptides and candidate with 1. Calls come from peptide counts, so an ALT site is not a validated genetic variant, and an “unestimated” class error means the run has no error rate for it.",  # §7.5 rules 3-4; docs/methods/variants.md
    "novel-orfs": "Candidates are ORFs supported by accepted peptides that do not fit an annotated gene in frame. With a numeric class FDR a call needs that q at the cut-off; with “unestimated” a call only meets the class rule and the minimum peptide count, so no error rate backs it.",  # docs/methods/novel-orfs.md, Call rule
    "function": "This section asks which functions have high abundance or are detected more or less than expected. Annotation coverage leans toward well-studied genes, and “no significant set” never means “no enrichment”.",  # §7.7.1; docs/methods/function-gsea.md
    "determinants": "This section asks how much of the variation in abundance sequence-derived features predict in held-out genes. All features are predictions and the data come from one condition, so these are associations, and native genes add regulation and turnover the features do not capture.",  # §7.8 caveats
    "sensitivity": "Each row reruns one analysis family under a changed analysis choice, such as a missing-data policy, and compares the flagged set with the primary run; Jaccard 1 means identical sets. No Jaccard cut-off is defined and it has no interval, so read it with the counts, because a single flip moves a small set a long way. Some rows cannot change by construction: missing=<p> keeps detection, so O − E flags are unchanged, and determinants rows under it are exact copies that say nothing about robustness.",  # docs/methods/sensitivity.md, Interpretation and "What the variants can change"
    "methods": "These notes say how each statistic was computed and where its claim holds; the full notes are in docs/methods/. q is a false-discovery-adjusted value, FWER is the chance of any false call in the family, and a TU is a transcription unit. The numbers below are the cut-offs this run used: settings, not measures of quality.",  # §8.1; config
}

#: per-plot text: (short name, text). Sources: GIPPA_SPEC.md §10.2 and the module docstrings in
#: gippa.viz, plus the methods note named beside the entry.
PLOT_READ: dict[str, tuple[str, str]] = {
    "V17": ("V17 input checks", "Peptide lengths: distinct peptides (y) per length in aa (x), placed or unmapped. Abundance: genes (y) per log2 abundance bin (x), measured or imputed; the collapsed panels add preflight counts, contig lengths and the share of genes missing each feature. Preflight errors stop a run; no cut-off is defined for the other panels."),  # §7.0
    "V05": ("V05 mapping flow", "Read left to right: input class, then route, then outcome; band width is the number of distinct sequences. Contaminant, rejected-by-FDR and unexplained bands are sequences not placed; no cut-off is defined for their size. A placed sequence is not by itself a detected gene."),
    "V06": ("V06 placement error", "(a) Targets and decoys per length class and number of differences, on a log axis with 0 at the floor; decoys estimate chance hits, so fewer relative to targets is better. (b) Stratum FDR per number of differences against the dashed placement cut-off; a point above the line means that stratum's estimated error exceeds it. These values are descriptive, not q-values: peptides are accepted on placement q per family, and a stratum with no target has no estimate."),  # docs/methods/placement-fdr.md; build.py peptides card
    "V07": ("V07 observed vs expected", "Each row is one group: the dot is observed minus expected detected genes and the bar its 95 % CI, so left of 0 means fewer detected than detectability predicts; no cut-off is defined for O − E. Only individual regions have a q, from a separate studentised test, so CI and q can disagree. Flag ci_few_units means do not read the CI as 95 %, and a group covering nearly all genes sits near 0 by construction."),  # docs/methods/detection-summary.md; ADR-0018
    "V18": ("V18 statistical diagnostics", "Reliability: mean predicted detection (x) against observed rate (y) in 10 equal-count bins; points on the dashed diagonal mean predictions match observations, and slope, intercept and Brier are descriptive, with no test or per-run cut-off. Each small panel is the p-value histogram of one test family: flat is what a family with no signal looks like, a left spike is signal or a too-liberal test (this plot cannot tell which), and no cut-off is defined. These panels do not validate the tests, and detection ORA rows flagged ora_uncalibrated stay uncalibrated."),  # docs/methods/detectability.md (§8.1 has no detectability claim); docs/methods/function.md Pass-1 limit
    "V01": ("V01 genome overview", "One sector per replicon; rings show genes by strand and detection, abundance, peptide density, transcription units, promoter and terminator ticks, regions, islands and GC skew, with markers for ALT sites and novel ORFs. Use it for orientation: it draws no test and has no cut-off."),
    "V02": ("V02 abundance along the genome", "Each point is a measured detected gene: position (x) against log2 abundance (y); undetected genes form the rug. The trend line is for display and nothing is claimed from it. Shading marks islands at the FWER cut-off, and no shading is not evidence that no island exists."),  # docs/methods/islands.md
    "V10": ("V10 autocorrelation", "Lag in genes (x) against Moran's I, or the join statistic for detection (y); bands show what shuffled gene labels and shuffled TU blocks produce, and a filled marker has TU-block q at the cut-off. Above the TU-block band suggests structure beyond operons. The TU-block null assumes exact TU boundaries and is anti-conservative at short lags when TUs are split, so read a short-lag signal with that in mind."),  # ADR-0061
    "V03": ("V03 loci", "Gene arrows are coloured by abundance (grey is undetected), with a peptide pileup, promoter and terminator marks and shaded regions. These views show the evidence at a locus and test nothing themselves; the q beside each locus is its region's O − E test. A missing promoter or terminator call never proves a missing element."),  # §7.4 quoted caveat SIGMA70
    "V11": ("V11 operon pairs", "(a) How far apart partners' abundances are (|Δ|, x) as a density per pair class: OP are predicted operon pairs, NO same-strand non-operon pairs, and CV and SW are decoys that cannot share an operon, so OP left of the others is what partners looking alike gives. (c) Upstream against downstream abundance; a pair counts as flagged only at q at the cut-off. (b, collapsed) Local fdr ranks NO pairs and does not make calls; no cut-off is defined for |Δ|."),  # §7.4; docs/methods/operons.md
    "V12": ("V12 variants", "(a) Bars count sites by call and confidence, split by flag; a site with several flags counts once in each, so segments can exceed the total. (b, collapsed) Peptide-against-genome alignments at the top ALT sites. A mass_equivalent flag means an equal-mass reading can explain the ALT, and strain_difference is applied to any ALT when the genome may be another isolate."),  # §7.5
    "V13": ("V13 novel ORFs", "Bars count candidates per class by status: called, not called, or unestimated (not called, with no class FDR). The locus views show the peptides at each called ORF. Without a numeric class FDR, a called ORF has no error rate."),  # docs/methods/novel-orfs.md
    "V08": ("V08 abundance by COG category", "Log2 abundance (y) of measured genes per COG category (x), sorted by median. This is descriptive: it tests nothing, and imputed values are not drawn."),
    "V14": ("V14 enrichment", "(a) Running score along TUs ranked by abundance: a peak near the start means the set ranks high, a trough near the end low. (b) NES (x; positive is high abundance, no interval) against −log10 q (y); hollow markers are not calibrated, redundant sets are small and grey, and GSEA is conditional on detection. (c, collapsed) Detected TUs observed over expected: above 1 is more than expected; no row has a CI, and rows flagged ora_uncalibrated stay uncalibrated."),  # ADR-0025, 0026, 0027; docs/methods/function.md
    "V15": ("V15 determinants", "(a) R² lost when a feature block is dropped; elastic net is primary, LightGBM a secondary check, and the intervals are uncalibrated. (d) Observed against out-of-fold predicted abundance: closer to the diagonal is better predicted, and marked genes (over 2 SD of residuals) are candidates for regulation or stability, not findings. No R² cut-off is defined, CAI, tAI and GC3 are collinear so read blocks and not single features, and the permutation p is uncalibrated, not a discovery claim."),  # §7.8; ADR-0031, 0033, 0034
}


class Raw(str):
    """Already rendered HTML."""


NOT_AVAILABLE: list[str] = []


def missing(pid: str, what: str) -> Raw:
    NOT_AVAILABLE.append(f"{pid}: {what}")
    return Raw('<span class="na">not available</span>')


def isna_(v: Any) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and pd.isna(v))


def num(v: Any) -> Raw:
    """A number for display: integers whole, above 100 one decimal, else 3 significant figures;
    the value as read stays in the tooltip."""
    x = float(v)
    if x == int(x) and abs(x) < 1e15:
        return Raw(str(int(x)))
    text = f"{x:.1f}" if abs(x) >= 100 else f"{x:.3g}"
    return Raw(f'<span title="{E(str(v))}">{text}</span>')


def pct(v: Any) -> Raw:
    return Raw(f'<span title="{E(str(v))}">{100 * float(v):.3g} %</span>')


def show(pid: str, what: str, v: Any) -> str:
    if isna_(v):
        return missing(pid, what)
    if isinstance(v, Raw):
        return v
    if isinstance(v, bool | str):
        return E(str(v))
    return num(v)


def tbl(pid: str, head: list[str], rows: list[list[Any]], ident: int = 1) -> str:
    """A compact table; a missing cell prints "not available" and is listed, named by the first
    ``ident`` cells of its row and its column."""
    body = ""
    for r in rows:
        who = " ".join(str(x) for x in r[:ident])
        body += "<tr>" + "".join(f"<td>{show(pid, f'{who} {h}', v)}</td>" for h, v in zip(head, r, strict=True)) + "</tr>"
    return f'<div class="sc"><table><thead><tr>{"".join(f"<th>{E(h)}</th>" for h in head)}</tr></thead><tbody>{body}</tbody></table></div>'


def kv(pid: str, items: list[tuple[str, Any]]) -> str:
    return '<p class="kv">' + "".join(f"<span><i>{E(k)}</i> {show(pid, k, v)}</span>" for k, v in items) + "</p>"


def sub(summary: str, body: str) -> str:
    return f'<details class="numsub"><summary>{E(summary)}</summary>{body}</details>'


def how(text: str, name: str = "") -> str:
    label = f"How to read this: {E(name)}" if name else "How to read this"
    return f'<div class="how"><b>{label}</b>{E(text)}</div>'


def numbers(body: str, label: str = "Numbers in this plot") -> str:
    return f'<div class="nums"><b>{E(label)}</b>{body}</div>'


class Ctx:
    def __init__(self, out: Path, view: RunView, cfg: GippaConfig, files: dict[str, Any], live: dict[str, Any], cards_: list[dict[str, Any]]) -> None:
        self.out, self.view, self.cfg, self.files, self.live, self.cards = out, view, cfg, files, live, cards_
        self.a = cfg.stats.alpha

    def data(self, inst: Any) -> pd.DataFrame:
        return pd.read_csv(self.out / "plots" / inst.data, sep="\t", keep_default_na=False, na_values=["NA"], low_memory=False)

    def first(self, pid: str) -> pd.DataFrame:
        return self.data(next(i for i in self.files[pid].instances if i.data))

    def stem(self, pid: str, stem: str) -> pd.DataFrame | None:
        inst = next((i for i in self.files[pid].instances if i.stem == stem and i.data), None)
        return None if inst is None else self.data(inst)

    def table(self, key: str) -> pd.DataFrame | None:
        return self.view.table(key) if self.view.has(key) else None


def nv(s: Any) -> Any:
    return pd.to_numeric(s, errors="coerce")


def alt_counts(ctx: Ctx) -> str:
    """The two ALT counts on the page, each with what it counts."""
    var = ctx.table("variants")
    if var is None:
        return kv("ALT count", [("ALT sites", None)])
    alt = var[var["call"] == "ALT"]
    n_mixed = int((var["call"] == "mixed").sum())
    high, cand = int((alt["confidence"] == "high").sum()), int((alt["confidence"] == "candidate").sum())
    marks = None
    if "V01" in ctx.files:
        v01 = ctx.first("V01")
        marks = int((v01["track"] == "alt_site").sum())
    if marks is None:
        tail = "The V01 marker count could not be read, so what the V01 title counts cannot be told from this run."
    elif marks == len(alt) + n_mixed:
        tail = f"{marks} is the V01 title and its markers: call = ALT plus call = mixed ({n_mixed})."
    else:
        tail = f"The V01 title and markers count {marks}, which matches neither sum, so what they count cannot be told from the data."
    return f'<p class="alt"><b>ALT counts on this page:</b> {len(alt)} is the summary tile and the V12 title: sites with call = ALT ({high} high, {cand} candidate). {tail}</p>'


def nums_summary(ctx: Ctx) -> str:
    by = {c["name"]: c for c in ctx.cards}
    items: list[tuple[str, Any]] = []
    genes = ctx.table("genes_tsv")
    if genes is not None:
        det = int(genes["detected"].astype(str).str.lower().eq("true").sum())
        items += [("genes in genes.tsv", len(genes)), ("detected", det), ("share detected", pct(det / len(genes)))]
    if "V05" in ctx.files:
        v05 = ctx.first("V05")
        pep = v05[v05["input_class"] == "peptide"]
        placed, total = int(pep.loc[pep["outcome"] == "placed", "n"].sum()), int(pep["n"].sum())
        items += [("peptides in", total), ("placed", placed), ("share placed", pct(placed / total) if total else None)]
    if "model_r2" in by:
        items += [("elastic-net CV R²", Raw(E(short(by["model_r2"]["value"]))))]
    return numbers(kv("Summary", items) + alt_counts(ctx), "Numbers behind the tiles")


def nums_v17(ctx: Ctx) -> str:
    d = ctx.first("V17")
    a = d[(d["panel"] == "a") & (d["group"] == "level")].set_index("key")["value"]
    out = kv("V17", [(f"preflight {k}", a.get(k)) for k in ("ok", "info", "warn", "error")])
    c, h = d[d["panel"] == "c"], d[d["panel"] == "d"]
    pl, un = c[c["group"] == "placed"], c[c["group"] == "unmapped"]
    meas, imp = h[h["group"] == "measured"], h[h["group"] == "imputed"]
    items: list[tuple[str, Any]] = [("peptides placed", nv(pl["value"]).sum()), ("unmapped", nv(un["value"]).sum())]
    if len(c):
        items += [("length range (aa)", Raw(f"{nv(c['key']).min():g}–{nv(c['key']).max():g}"))]
    items += [("genes measured", nv(meas["value"]).sum()), ("imputed", nv(imp["value"]).sum())]
    if len(meas):
        peak = meas.loc[nv(meas["value"]).idxmax()]
        items += [("log2 range (bin centres)", Raw(f"{nv(meas['key']).min():.1f}–{nv(meas['key']).max():.1f}")), ("tallest bin", Raw(f"{float(peak['key']):.1f}"))]
    out += kv("V17", items)
    rows: list[list[Any]] = [["contig length", str(r["key"]), str(r["group"]), r["value"]] for _, r in d[d["panel"] == "b"].iterrows()]
    e = d[d["panel"] == "e"].assign(v=lambda f: nv(f["value"]))
    for g, sub_ in e.groupby("group", sort=True):
        top = sub_.loc[sub_["v"].idxmax()]
        rows.append(["highest missing share", str(g), str(top["key"]), pct(top["v"])])
    return out + tbl("V17", ["panel", "item", "detail", "value"], rows, 2)


def nums_v05(ctx: Ctx) -> str:
    d = ctx.first("V05")
    rows: list[list[Any]] = []
    for cls in ("peptide", "protein"):
        g = d[d["input_class"] == cls]
        total = int(g["n"].sum())
        if not total:
            rows.append([cls, "all", "no input", 0, Raw("–")])
            continue
        rows += [[cls, str(r["route"]), str(r["outcome"]), r["n"], pct(r["n"] / total)] for _, r in g.iterrows() if r["n"]]
        rows.append([cls, "all", "input", total, pct(1.0)])
    conflicts = d[d["outcome"] == "conflicts"]["n"]
    tail = kv("V05", [("placed and unmapped (conflicts)", conflicts.iloc[0] if len(conflicts) else None)])
    return tbl("V05", ["input", "route", "outcome", "sequences", "share of input"], rows, 3) + tail


def nums_v06(ctx: Ctx) -> str:
    d, alpha = ctx.first("V06"), ctx.a.placement
    b = d[d["panel"] == "b"]
    a = d[d["panel"] == "a"]
    strata = sorted({int(m) for m in nv(a["n_diff"]).dropna()} | {int(m) for m in nv(b["n_diff"]).dropna()})
    rows = []
    for m in strata:
        sb = b[nv(b["n_diff"]) == m].set_index("metric")["value"]
        f = sb.get("fdr_stratum")
        side = "no estimate" if isna_(f) else ("above cut-off" if float(f) > alpha else "at or below")
        rows.append([f"{m} difference{'' if m == 1 else 's'}", f, sb.get("n_accepted"), side])
    out = kv("V06", [("placement cut-off (dashed)", alpha)]) + tbl("V06", ["stratum", "stratum FDR", "accepted", f"against {alpha:g}"], rows)
    wide = []
    for lc in dict.fromkeys(a["length_class"].astype(str)):
        row: list[Any] = [lc]
        for m in strata:
            g = a[(a["length_class"].astype(str) == lc) & (nv(a["n_diff"]) == m)].set_index("metric")["value"]
            row.append(Raw(f"{show('V06', f'{lc} m={m} targets', g.get('targets'))} / {show('V06', f'{lc} m={m} decoys', g.get('decoys'))}"))
        wide.append(row)
    out += tbl("V06", ["length class (aa)", *[f"{m} diff: targets / decoys" for m in strata]], wide)
    c = d[d["panel"] == "c"]
    if len(c):
        out += tbl("V06", ["score bin", "de novo precision", "CI low", "CI high"], [[r["score_bin"], r["value"], r["ci_low"], r["ci_high"]] for _, r in c.iterrows()])
    return out


def nums_v07(ctx: Ctx) -> str:
    d, alpha = ctx.first("V07"), ctx.a.context
    rows: list[list[Any]] = []
    for _, r in d.iterrows():
        n = float(r["n"]) if pd.notna(r["n"]) else 0.0
        reg = r["block"] == "region"
        q: Any = r["q_value"] if pd.notna(r["q_value"]) else ("untested" if reg else "no test")
        share = lambda col: pct(float(r[col]) / n) if n and pd.notna(r[col]) else None
        flags = "" if pd.isna(r["flags"]) else str(r["flags"])
        rows.append([str(r["block"]).replace("_", " "), str(r["group"]), r["n"], r["observed"], r["expected"], share("observed"), share("expected"), r["o_minus_e"], r["ci_low"], r["ci_high"], q, flags])
    tested = d[pd.notna(nv(d["q_value"]))]
    head = kv("V07", [("regions tested", len(tested)), (f"regions at q ≤ {alpha:g}", int((nv(tested["q_value"]) <= alpha).sum()))])
    return head + tbl("V07", ["block", "group", "genes", "observed", "expected", "% obs", "% exp", "O − E", "CI low", "CI high", "q", "flags"], rows, 2)


def nums_v18(ctx: Ctx, drawn: list[tuple[str, str]]) -> str:
    d = ctx.first("V18")
    a = d[d["panel"] == "a"]
    m = a[a["metric"].notna()].set_index("metric")["value"]
    pts = a[a["metric"].isna()].sort_values("bin", key=nv)
    out = kv("V18", [("slope", m.get("calibration_slope")), ("intercept", m.get("calibration_intercept")), ("Brier", m.get("brier")), ("eligible genes", m.get("n_eligible")), ("detected", m.get("n_detected"))])
    head = ["bin", *[str(int(b)) for b in nv(pts["bin"])]]
    out += tbl("V18", head, [["predicted", *nv(pts["x"])], ["observed", *nv(pts["y"])], ["genes", *nv(pts["n"])]])
    b = d[(d["panel"] == "b") & d["metric"].isna()]
    rows = []
    for table, fam in drawn:
        g = b[(b["table"] == table) & (b["family"] == fam)].assign(lo=lambda f: nv(f["bin"]), k=lambda f: nv(f["n"]))
        total = g["k"].sum()
        rows.append([table.replace("_", " "), fam, total, g.loc[g["lo"] < 0.05 - 1e-9, "k"].sum(), total / 20, g.loc[g["lo"] >= 0.95 - 1e-9, "k"].sum(), g["k"].max()])
    out += tbl("V18", ["table", "family", "p-values", "p < 0.05", "flat bin height", "p ≥ 0.95", "tallest bin"], rows, 2)
    fams = b.assign(k=lambda f: nv(f["n"])).groupby(["table", "family"], sort=False)["k"].sum()
    rest = [f"{t} {f} ({int(n)})" for (t, f), n in fams.items() if (t, f) not in drawn]
    if rest:
        out += f'<p class="kv"><span><i>Families not drawn (n)</i> {E("; ".join(rest))}</span></p>'
    c = d[(d["panel"] == "c") & (d["metric"] == "observed")]
    if len(c):
        out += sub(f"Null panels, SVG only: observed statistic ({len(c)})", tbl("V18", ["table", "family", "observed"], [[str(r["table"]), str(r["family"]), r["value"]] for _, r in c.iterrows()], 2))
    return out


def nums_v01(ctx: Ctx) -> str:
    d = ctx.first("V01")
    count = lambda t: int((d["track"] == t).sum())
    rows: list[list[Any]] = []
    for _, r in d[d["track"] == "sector"].iterrows():
        rows.append(["replicon length (bp)", str(r["contig"]), r["value"]])
    for t in ("gene_plus", "gene_minus", "region", "alt_site"):
        for cat, n in d.loc[d["track"] == t, "category"].value_counts().sort_index().items():
            rows.append([t.replace("_", " "), str(cat), int(n)])
    rows.append(["abundance heat track", "genes", count("abundance")])
    pd_ = nv(d.loc[d["track"] == "peptide_density", "value"])
    rows.append(["peptide density", "bins", len(pd_)])
    rows.append(["peptide density", "highest bin", pd_.max() if len(pd_) else None])
    for t in ("tu_start", "promoter", "terminator", "island", "novel_orf"):
        rows.append([t.replace("_", " "), "marks", count(t)])
    for t in ("oric", "ter"):
        g = d[d["track"] == t]
        for _, r in g.iterrows():
            rows.append([t, str(r["contig"]), r["start"]])
    return tbl("V01", ["track", "item", "value"], rows, 2)


def nums_v02(ctx: Ctx) -> str:
    d = ctx.first("V02")
    rows = []
    for contig, g in d.groupby("contig", sort=False):
        pts = g[g["kind"] == "point"]
        val = nv(pts["value"])
        rug = g[g["kind"] == "rug"]["category"]
        rows.append([str(contig), len(pts), int((rug == "not_detected").sum()), int((rug == "partial_excluded").sum()), val.min() if len(val) else Raw("–"), val.max() if len(val) else Raw("–"), int((g["kind"] == "island").sum())])
    tail = kv("V02", [(k, r["position_bp"]) for k, sel in (("oriC (bp)", "oric"), ("terminus (bp)", "ter")) for _, r in d[d["kind"] == sel].iterrows()])
    return tbl("V02", ["replicon", "points (measured)", "rug: not detected", "rug: partial", "lowest log2", "highest log2", "islands shaded"], rows) + tail


def nums_v10(ctx: Ctx) -> str:
    d, alpha = ctx.first("V10"), ctx.a.context
    rows = []
    for (contig, stat, null), g in d.groupby(["contig", "statistic", "null_model"], sort=False):
        flags = g["flags"].fillna("").astype(str)
        na = flags.str.contains("not_estimable|few_values") | nv(g["value"]).isna()
        ok = g[~na]
        one = ok[nv(ok["lag"]) == 1]
        if len(one):
            r = one.iloc[0]
            obs = float(r["effect"]) + float(r["null_mean"]) if stat == "join_count_bb" else r["value"]
            vals: list[Any] = [obs, r["null_lo"], r["null_hi"], r["q_value"]]
        else:
            vals = ["not estimable"] * 4
        rows.append([str(contig), "Moran's I" if stat == "morans_i" else "BB join", "TU-block" if null == "tu_block" else "gene-label", *vals, int((nv(ok["q_value"]) <= alpha).sum()), len(ok)])
    return kv("V10", [("q cut-off", alpha)]) + tbl("V10", ["replicon", "statistic", "null", "lag 1", "band low", "band high", "q", "lags at cut-off", "lags estimable"], rows, 3)


def nums_v03(ctx: Ctx) -> str:
    index = ctx.stem("V03", "V03_locus")
    if index is None:
        return kv("V03", [("locus index", None)])
    rows = []
    for _, r in index.iterrows():
        loc = ctx.stem("V03", f"V03_locus_{r['locus_id']}")
        def n(track: str, cat: str | None = None) -> Any:
            if loc is None:
                return None
            g = loc[loc["track"] == track]
            return len(g if cat is None else g[g["category"] == cat])
        q = r["error"] if pd.notna(r["error"]) else "no test"
        rows.append([r["rank"], str(r["locus_id"]), str(r["kind"]), str(r["contig"]), r["start"], r["end"], q, n("gene", "detected"), n("gene", "not_detected"), n("promoter"), n("terminator"), str(r["rendered"])])
    return tbl("V03", ["rank", "locus", "kind", "contig", "start", "end", "q", "genes detected", "genes not", "promoter marks", "terminator marks", "drawn"], rows, 2)


def nums_v11(ctx: Ctx) -> str:
    d, alpha = ctx.first("V11"), ctx.a.operon_flags
    a = d[(d["panel"] == "a") & d["n"].notna()].set_index("class")["n"]
    c = d[d["panel"] == "c"]
    rows = []
    for cls in ("OP", "NO", "CV", "SW", "DV"):
        if cls in a.index:
            med = nv(c.loc[c["class"] == cls, "abs_delta_log2"])
            rows.append([cls, a[cls], med.median() if len(med) else Raw("–"), "yes" if cls != "DV" else "no: table only"])
    out = tbl("V11", ["class", "pairs in (a)", "median |Δ| (log2)", "density drawn"], rows)
    flags = c[c["flag"].notna()]
    q = nv(flags["q_value"])
    lf = nv(d.loc[d["panel"] == "b", "lfdr"])
    out += kv("V11", [("flag cut-off q", alpha), ("flagged at cut-off", int((q <= alpha).sum())), ("candidate rows", len(flags)), ("smallest q", q.min() if len(q) else Raw("–")), ("lfdr values (NO)", len(lf)), (f"lfdr ≤ {alpha:g}", int((lf <= alpha).sum())), ("smallest lfdr", lf.min() if len(lf) else Raw("–"))])
    for f, n in flags["flag"].value_counts().items():
        out += f'<p class="kv"><span><i>{E(str(f))}</i> {int(n)}</span></p>'
    return out


def nums_v12(ctx: Ctx) -> str:
    d = ctx.first("V12")
    a = d[d["panel"] == "a"]
    b = d[d["panel"] == "b"]
    flags = list(dict.fromkeys(a["flag"].astype(str)))
    flags = [f for f in flags if f != "all"] + ["all"]
    rows = []
    for (call, conf), g in a.groupby(["call", "confidence"], sort=False):
        s = g.set_index("flag")["n"]
        rows.append([str(call), str(conf), *[int(s[f]) if f in s.index else 0 for f in flags]])
    out = alt_counts(ctx) + tbl("V12", ["call", "confidence", *[f.replace("_", " ") for f in flags]], rows, 2)
    if len(b):
        top = [[r["rank"], str(r["variant_id"]), str(r["gene"]), r["aa_pos"], f"{r['ref_aa']}→{r['alt_aa']}", str(r["alt_type"]), r["n_alt"], r["n_ref"], str(r["confidence"]), str(r["class_error_rate"]) if pd.notna(r["class_error_rate"]) else None] for _, r in b.iterrows()]
        if (b["alt_aa"].astype(str) == "J").any() or (b["ref_aa"].astype(str) == "J").any():
            out += '<p>J in the change column means I or L, which mass spectrometry cannot tell apart.</p>'
        out += sub(f"Panel (b): the {len(b)} sites aligned", tbl("V12", ["rank", "site", "gene", "aa", "change", "type", "ALT peptides", "REF peptides", "confidence", "class error"], top, 2))
    return out


def nums_v13(ctx: Ctx) -> str:
    d = ctx.stem("V13", "V13_novel_orfs")
    if d is None:
        return kv("V13", [("class counts", None)])
    rows = [[str(r["class"]), str(r["status"]), r["n"], r["class_fdr_max_called"] if pd.notna(r["class_fdr_max_called"]) else "unestimated"] for _, r in d.iterrows()]
    out = tbl("V13", ["class", "status", "ORFs", "largest class FDR of called"], rows, 2)
    orfs = ctx.table("novel_orfs")
    if orfs is not None:
        called = orfs[orfs["called"].astype(str).str.lower() == "true"]
        cols = [c for c in ("orf_id", "class", "strand", "length_aa_best", "n_peptides", "best_peptide_q", "class_fdr") if c in called]
        if len(called):
            out += tbl("V13", [c.replace("_", " ") for c in cols], [[r[c] if not isinstance(r[c], str) else str(r[c]) for c in cols] for _, r in called.iterrows()])
    return out


def nums_v08(ctx: Ctx) -> str:
    d = ctx.first("V08")
    s = d[d["locus_tag"].isna()]
    rows = []
    for _, r in s.iterrows():
        med = r["median_log2"] if pd.notna(r["median_log2"]) else "none measured"
        frac = pct(float(r["n_detected"]) / float(r["n_genes"])) if float(r["n_genes"]) else None
        rows.append([str(r["category"]), str(r["super_category"]), med, r["n_genes"], r["n_detected"], frac, r["n_imputed"]])
    return tbl("V08", ["COG", "group", "median log2", "genes", "detected", "detected share", "imputed"], rows)


def nums_v14(ctx: Ctx) -> str:
    d, alpha = ctx.first("V14"), ctx.a.enrichment
    a = d[d["panel"] == "a"].drop_duplicates(["collection", "term"])
    b, c = d[d["panel"] == "b"], d[d["panel"] == "c"]
    flag_text = lambda v: "" if pd.isna(v) else str(v).replace("ci_unestimated;", "").replace("ci_unestimated", "")
    out = tbl("V14", ["(a) set", "name", "collection", "NES", "q", "flags"], [[str(r["term"]), str(r["name"]), str(r["collection"]), r["nes"], r["q_value"], flag_text(r["flags"])] for _, r in a.iterrows()])
    rows = []
    red = b["redundant"].astype(str).str.lower() == "true"
    bad = b["flags"].fillna("").astype(str).str.contains("genomic_cluster|detection_confounded")
    q = nv(b["q_value"])
    for coll in dict.fromkeys(b["collection"].astype(str)):
        m = b["collection"].astype(str) == coll
        nes = nv(b.loc[m, "nes"])
        rows.append([coll, int(m.sum()), int((m & (q <= alpha)).sum()), int((m & (q <= alpha) & ~red).sum()), int((m & (q <= alpha) & ~red & bad).sum()), nes.min(), nes.max()])
    rows.append(["all", len(b), int((q <= alpha).sum()), int(((q <= alpha) & ~red).sum()), int(((q <= alpha) & ~red & bad).sum()), nv(b["nes"]).min(), nv(b["nes"]).max()])
    out += tbl("V14", ["(b) collection", "sets drawn", f"q ≤ {alpha:g}", "of those non-redundant", "of those hollow", "lowest NES", "highest NES"], rows)
    if len(c):
        qc = nv(c["q_value"])
        oe = nv(c["o_over_e"])
        n_unc = int(c["flags"].fillna("").astype(str).str.contains("ora_uncalibrated").sum())
        body = kv("V14", [("rows drawn", len(c)), (f"q ≤ {alpha:g}", int((qc <= alpha).sum())), ("over", int((c["direction"] == "over").sum())), ("under", int((c["direction"] == "under").sum())), ("lowest O/E", oe.min()), ("highest O/E", oe.max()), ("flagged ora_uncalibrated", n_unc)])
        order = c.assign(_q=qc).sort_values("_q", kind="stable")
        body += sub(f"Panel (c): all {len(c)} rows by q", tbl("V14", ["category", "name", "collection", "direction", "O/E", "q", "coverage"], [[str(r["term"]), str(r["name"]), str(r["collection"]), str(r["direction"]), r["o_over_e"], r["q_value"], r["coverage"]] for _, r in order.iterrows()], 3))
        out += body
    return out


def nums_v15(ctx: Ctx) -> str:
    d = ctx.first("V15")
    a, dd = d[d["panel"] == "a"], d[d["panel"] == "d"]
    items: list[tuple[str, Any]] = []
    model = ctx.table("determinants_model")
    if model is not None:
        for name, label in (("elastic_net", "elastic-net"), ("lightgbm", "LightGBM")):
            m = model[(model["model"] == name) & (model["feature_or_block"] == "all")].set_index("metric")["value"]
            items += [(f"{label} CV R² mean", nv(m.get("cv_r2_mean"))), ("SD", nv(m.get("cv_r2_sd"))), ("null p (uncalibrated)", nv(m.get("null_p")))]
    out = kv("V15", items)
    rows = []
    for block in dict.fromkeys(a["block"].astype(str)):
        row: list[Any] = [block]
        for name in ("elastic_net", "lightgbm"):
            g = a[(a["block"] == block) & (a["model"] == name)]
            r = g.iloc[0] if len(g) else None
            row += [r["x"], r["ci_low"], r["ci_high"]] if r is not None else [None, None, None]
        rows.append(row)
    out += tbl("V15", ["(a) block", "elastic-net ΔR²", "low", "high", "LightGBM ΔR²", "low", "high"], rows)
    cats = dd["category"].value_counts()
    out += kv("V15", [("(d) model genes", int(cats.get("model", 0)) + int(cats.get("unexplained", 0))), ("of those unexplained", int(cats.get("unexplained", 0))), ("ribosomal reference (not in R²)", int(cats.get("ribosomal_reference", 0)))])
    shap = ctx.view.support("determinants_shap")
    if shap is not None and len(shap):
        m = nv(shap["shap"]).abs().groupby(shap["feature"].astype(str)).mean().sort_values(ascending=False, kind="stable").head(20)
        out += sub("Panel (c): top 20 features, mean |SHAP| (log2)", tbl("V15", ["feature", "mean |SHAP|"], [[k, float(v)] for k, v in m.items()]))
    else:
        out += kv("V15", [("(c) SHAP values", None)])
    b = d[(d["panel"] == "b") & (d["category"] == "pd_mean")]
    if len(b):
        pdr = []
        for f, g in b.groupby("feature", sort=False):
            x, y = nv(g["x"]), nv(g["y"])
            rug = d[(d["panel"] == "b") & (d["category"] == "rug") & (d["feature"] == f)]
            pdr.append([str(f), x.min(), x.max(), y.min(), y.max(), len(rug)])
        out += sub("Panel (b): partial dependence ranges", tbl("V15", ["feature", "lowest x", "highest x", "lowest prediction", "highest prediction", "genes in rug"], pdr))
    else:
        out += kv("V15", [("(b) partial dependence", None)])
    return out


def nums_sensitivity(ctx: Ctx) -> str:
    """Each family's primary set against the count another section shows: the definitions differ."""
    whole = ctx.table("sensitivity")
    if whole is None or not len(whole):
        return numbers(kv("Sensitivity", [("sensitivity.tsv", None)]), "Numbers in this section")
    prim = whole.drop_duplicates("family").set_index("family")["n_primary"]
    a = ctx.a
    rows: list[list[Any]] = []
    if "regions" in prim:
        r = ctx.table("regions_stats")
        n = None
        if r is not None:
            keep = (r["null_model"].astype(str) == "circular_shift") & (nv(r["q_value"]) <= a.context)
            n = int(keep.sum())
        rows.append(["regions", prim["regions"], n, f"regions_stats rows with q ≤ {a.context:g}; both count (region, statistic)"])
    if "islands" in prim:
        i = ctx.table("islands")
        rows.append(["islands", prim["islands"], None if i is None else int(nv(i["n_genes"]).sum()), "genes in reported islands; the rerun is a separate Monte Carlo draw"])
    if "operon_flags" in prim:
        o = ctx.table("operon_flags")
        n = None if o is None else int((nv(o["q_value"]) <= a.operon_flags).sum())
        rows.append(["operon_flags", prim["operon_flags"], n, f"operon_flags.tsv rows with q ≤ {a.operon_flags:g} (tile); the sensitivity count takes every row of that table{'' if o is None else f', {len(o)} here'}, whatever its q"])
    if "enriched_categories" in prim:
        g, ora = ctx.table("enrichment_gsea"), ctx.table("enrichment_detection")
        n = None if g is None or ora is None else int((nv(g["q_value"]) <= a.enrichment).sum() + (nv(ora["q_value"]) <= a.enrichment).sum())
        rows.append(["enriched_categories", prim["enriched_categories"], n, f"GSEA and ORA rows with q ≤ {a.enrichment:g}, redundant included; the rerun is a separate draw at fewer permutations. The V14 title leaves redundant sets out"])
    if "determinants" in prim:
        att = ctx.table("determinants_attribution")
        n = None
        if att is not None and "V15" in ctx.files:
            d = ctx.first("V15")
            n = int(((d["panel"] == "d") & (d["category"] == "unexplained")).sum())
        rows.append(["determinants", prim["determinants"], n, "genes marked unexplained in V15 (d); the rerun uses its own fold draw, so the sets need not match"])
    return numbers(tbl("Sensitivity", ["family", "n primary (this table)", "same family elsewhere", "what the second count is"], rows), "Numbers in this section")


def nums_methods(ctx: Ctx) -> str:
    a, s = ctx.a, ctx.cfg.stats
    items: list[tuple[str, Any]] = [(f"alpha.{k}", getattr(a, k)) for k in ("identification", "placement", "context", "operon_flags", "variants", "novel_orfs", "enrichment", "determinants")]
    items += [("stats.permutations", s.permutations), ("stats.bootstrap", s.bootstrap), ("genetic code", ctx.cfg.annotation.translation_table)]
    return numbers(kv("Methods", items), "Numbers in this section")


NUMS = {"V17": nums_v17, "V05": nums_v05, "V06": nums_v06, "V07": nums_v07, "V01": nums_v01, "V02": nums_v02,
        "V10": nums_v10, "V03": nums_v03, "V11": nums_v11, "V12": nums_v12, "V13": nums_v13, "V08": nums_v08,
        "V14": nums_v14, "V15": nums_v15}  # fmt: skip


def read_block(pid: str, ctx: Ctx) -> str:
    """The plot's reading guide beside its numbers. A failed read prints "not available": listed."""
    name, text = PLOT_READ[pid]
    try:
        if pid == "V18":
            inst = next(i for i in ctx.files[pid].instances if i.stem in ctx.live)
            found = (re.match(r"(.+?)<br>(.+), n = ", p["title"] or "") for p in explode(json.loads(ctx.live[inst.stem][1])))
            body = nums_v18(ctx, [(m.group(1), m.group(2)) for m in found if m])
        else:
            body = NUMS[pid](ctx)
    except Exception as exc:
        body = str(missing(pid, f"numbers could not be read ({type(exc).__name__}: {exc})"))
    return f'<div class="rn">{how(text, name)}{numbers(body)}</div>'


def section_body(plots: tuple[str, ...], files: dict[str, Any], live: dict[str, Any], page: Page, out: Path, ctx: Ctx) -> tuple[str, str, str]:
    """(visible figures, collapsed figures, count summary)."""
    vis, hid = "", ""
    for pid in plots:
        pf = files.get(pid)
        if pf is None or not pf.instances:
            vis += f'<p class="note">{pid}: not produced in this run.</p>'
            continue
        vis += read_block(pid, ctx)
        if pid == "V03":
            loci, index = loci_by_rank(out, pf.instances)
            vis += "".join(static_block(i, 1240, out) for i in loci[:2] if i.svg)
            rest = "".join(static_block(i, 1240, out) for i in [*index, *loci[2:]] if i.svg)
            hid += more(f"{len(loci[2:])} more loci and the ranked index", rest)
            continue
        if pid == "V12":
            hid += more("Peptide-vs-genome alignments at the top ALT sites (V12)", "".join(static_block(i, 1100, out) for i in pf.instances if i.svg))
            continue
        if pid == "V13":
            first, *others = pf.instances
            if first.svg:
                vis += static_block(first, 760, out)
            loci, index = loci_by_rank(out, others)
            hid += more(f"Locus views of {len(loci)} novel ORFs", "".join(static_block(i, 1000, out) for i in [*index, *loci] if i.svg))
            continue
        for inst in pf.instances:
            if inst.stem in live and pid in PANELS:
                block, _ = twin_block(pid, inst, live[inst.stem][1], page, out)
                vis += block
            elif inst.svg:
                vis += static_block(inst, 1000 if pid == "V01" else 1100, out)
    return vis, hid, ""


CSS = """
body{font-family:"DejaVu Sans",Arial,sans-serif;margin:0;color:#222;font-size:15px;line-height:1.45;background:#fff}
main{max-width:1280px;margin:0 auto;padding:0 24px 60px}
nav{position:sticky;top:0;z-index:5;background:#123;color:#fff;padding:6px 24px;display:flex;flex-wrap:wrap;gap:4px 14px;font-size:13px}
nav a{color:#fff;text-decoration:none} nav a:hover{text-decoration:underline} nav b{margin-right:10px}
h1{font-size:1.5em;margin:18px 0 0} section{padding:10px 0 18px;border-top:1px solid #ddd;scroll-margin-top:40px}
#summary{border-top:0} h2{font-size:1.25em;margin:10px 0 4px} h3{font-size:1em;margin:0 0 6px}
.meta{color:#555;margin:2px 0 12px}
.tiles{display:grid;grid-template-columns:repeat(4,1fr);gap:10px}
.tile{border:1px solid #ccd;border-radius:6px;padding:10px 14px;background:#fafbfd}
.tile.dim{opacity:.55}.tv{font-size:1.7em;font-weight:700;line-height:1.2}.tl{font-size:.92em}.te,.td{font-size:.8em;color:#556}
.kp{display:grid;grid-template-columns:3fr 2fr;gap:18px;margin-top:14px}
.kp ul{margin:0;padding-left:1.1em}.kp li{margin-bottom:5px}.kp li a{font-weight:600;margin-right:3px}
.care{border-left:4px solid #8a5a00;background:#fff8e6;padding:10px 14px;border-radius:0 6px 6px 0}
.lead{font-size:1.05em;margin:0 0 10px;font-weight:600;color:#123}
.grid{display:grid;gap:14px 18px;margin:8px 0}.grid.c1{grid-template-columns:1fr}.grid.c2{grid-template-columns:1fr 1fr}.grid.c3{grid-template-columns:1fr 1fr 1fr}
.panel{min-width:0}.ptitle{font-weight:600;font-size:14px;margin-bottom:2px}
.static{margin:10px 0}.static svg{width:100%;height:auto;display:block}
.flinks{font-size:.78em;color:#777;margin-top:2px}.flinks a{margin-left:8px;color:#369}
details.more{margin:10px 0;border:1px solid #dde;border-radius:6px;padding:0 14px}
details.more>summary{cursor:pointer;padding:9px 0;font-weight:600;color:#245}
details.more[open]{padding-bottom:10px}
.tbl{margin:8px 0 16px}.thead{display:flex;gap:12px;align-items:center;flex-wrap:wrap;font-size:.85em;margin-bottom:4px}.thead span{color:#667}
.thead input[type=search]{font-size:.95em;padding:2px 6px}.xt{color:#667}
.table-wrap{overflow-x:auto}table{border-collapse:collapse;font-size:.82em}
th,td{border-bottom:1px solid #e6e6ea;padding:3px 9px;text-align:left;white-space:nowrap}th{background:#f3f4f7;cursor:pointer}
.tbl:not(.allcols) .x{display:none}
.fnote{font-size:.85em;color:#555;margin:.2em 0}.note{color:#666;font-style:italic}.caveat{font-size:.9em;margin:.35em 0}
.prov td,.prov th{white-space:normal;overflow-wrap:anywhere}.prov{width:100%;table-layout:auto}
.rn{display:grid;grid-template-columns:1fr 2fr;gap:12px;margin:10px 0 6px;align-items:start}
.how{background:#f1f6fb;border-left:3px solid #369;padding:6px 10px;font-size:.82em;line-height:1.34;border-radius:0 4px 4px 0;margin:6px 0}
.rn .how{margin:0}.how b,.nums>b{display:block;font-size:.95em;color:#123;margin-bottom:2px}
.nums{border:1px solid #dde;border-radius:4px;padding:5px 9px;font-size:.8em;max-height:164px;overflow:auto;min-width:0}
section>.nums{max-height:none;margin:6px 0}.nums p{margin:2px 0 4px}
.nums table{font-size:1em;margin:3px 0 5px}.nums th{cursor:default;position:sticky;top:-5px;white-space:normal;vertical-align:bottom;font-weight:600}.nums td,.nums th{padding:1px 7px}
.kv{margin:2px 0 4px;display:flex;flex-wrap:wrap;gap:1px 12px}.kv i{font-style:normal;color:#556}
.na{color:#a33;font-style:italic}.alt{margin:3px 0 5px}.numsub>summary{cursor:pointer;color:#245;margin:2px 0}
@media(max-width:900px){.tiles{grid-template-columns:1fr 1fr}.kp,.grid.c2,.grid.c3,.rn{grid-template-columns:1fr}}
"""

JS = """
(function(){
var done=new WeakSet();
var config={responsive:true,displayModeBar:false};
function draw(d){if(done.has(d))return;done.add(d);var s=JSON.parse(document.getElementById(d.dataset.spec).textContent);Plotly.newPlot(d,s.data,s.layout,config);}
var io=new IntersectionObserver(function(es){es.forEach(function(e){if(e.isIntersecting){draw(e.target);io.unobserve(e.target);}});},{rootMargin:'900px'});
document.querySelectorAll('.pl').forEach(function(d){io.observe(d);});
window.__drawAll=function(){document.querySelectorAll('.pl').forEach(function(d){if(d.offsetParent!==null)draw(d);});};
document.querySelectorAll('[data-allcols]').forEach(function(c){c.addEventListener('change',function(){c.closest('.tbl').classList.toggle('allcols',c.checked);});});
})();
"""


def build(out: Path) -> str:
    _OUT[:] = [out]
    cfg = GippaConfig.model_validate(yaml.safe_load((out / "resolved_config.yaml").read_text()))
    load_plot_modules()
    view = RunView.from_outdir(out).with_modules(cfg.modules)
    files = read_plot_files(out / "plots")
    live = twins(out / "plots", files, twin_meta(view, cfg))
    summary, notes = cards(view, cfg)
    doc, note = content.read_manifest(view)
    titles = titles_of(files)
    extras = content.extras(view, doc)
    page = Page()
    NOT_AVAILABLE.clear()
    ctx = Ctx(out, view, cfg, files, live, summary)
    tables_js = (Path(content.__file__).parent / "templates" / "tables.js").read_text()
    parts = []
    for s in SECTIONS:
        vis, hid, _ = section_body(s.plots, files, live, page, out, ctx)
        lead_pid = LEAD.get(s.id)
        lead = f'<p class="lead">{E(short(titles[lead_pid]))}</p>' if lead_pid in titles else ""
        body = ""
        if s.id == "summary":
            demo = content.demo_banner(doc)
            body = (f'<p class="note">{E(demo)}</p>' if demo else "") + key_points(view, cfg, files, titles, summary, notes, doc)
        tkeys = [k for k in s.tables if k != "sensitivity"] + (["variants"] if s.id == "variants" else [])
        tbl = "".join(table_block(view, k) for k in tkeys)
        cav, n_cav = ("", 0)
        if s.id != "summary":
            cav, n_cav = extras_block(extras.get(s.id, []))
        foot = hid
        if tbl:
            foot += more(f"Tables ({len(tkeys)}): top rows, key columns", tbl, open_=s.id == "novel-orfs")
        if cav:
            foot += more(f"Caveats and assumptions ({n_cav})" if n_cav else "Diagnostics", cav)
        if s.id == "sensitivity":
            whole = view.table("sensitivity") if view.has("sensitivity") else None
            if whole is not None and len(whole):
                j = pd.to_numeric(whole["jaccard"], errors="coerce").dropna()
                if len(j):
                    lead = f'<p class="lead">Flagged sets under the missing-data variants: Jaccard {j.min():.2f} to {j.max():.2f} against the primary run.</p>'
                rows = [r for _, g in content.sensitivity_groups(whole) for r in g]
                head = "".join(f"<th>{E(c)}</th>" for c in whole.columns)
                tr = "".join("<tr>" + "".join(f"<td>{E(str(v))}</td>" for v in r) + "</tr>" for r in rows)
                body = f'<div class="table-wrap"><table data-table="sensitivity"><thead><tr>{head}</tr></thead><tbody>{tr}</tbody></table></div>'
        if s.id == "methods":
            snippets = "".join(f'<h3>{E(t)}</h3>{b}' for t, b in content.snippets())
            lines = f"<p>{E(content.genetic_code_methods(cfg.annotation.translation_table))}</p>"
            lic = content.licence_notes(doc)
            if lic:
                lines += f"<p>Licence notes: {E('; '.join(lic))}</p>"
            foot += more("Methods notes", lines + snippets + cav)
        if s.id == "provenance":
            prov = content.provenance(view, doc, note)
            rows = dict(prov.get("lines", []))
            lead = f'<p class="lead">GIPPA {E(rows.get("GIPPA version", __version__))}, commit {E(rows.get("commit", "not recorded")[:10])}, seed {E(rows.get("seed", "not recorded"))}.</p>'
            tab = lambda head, rs: f'<table class="prov"><thead><tr>{"".join(f"<th>{E(h)}</th>" for h in head)}</tr></thead><tbody>' + "".join("<tr>" + "".join(f"<td>{E(str(v))}</td>" for v in r) + "</tr>" for r in rs) + "</tbody></table>"
            det = tab(["item", "value"], prov.get("lines", [])) + f'<p>Resolved configuration SHA-256: {E(str(prov.get("config")))}</p>' + tab(["file", "role", "SHA-256"], prov["inputs"])
            if "tools" in prov:
                det += tab(["tool", "version", "digest"], prov["tools"]) + tab(["database", "version", "checksum", "retrieved"], prov.get("databases", []))
            foot += more("Full provenance: inputs, tools, databases, digests", det)
        heading = "Summary: key points" if s.id == "summary" else s.heading
        guide = how(SECTION_READ[s.id]) if s.id in SECTION_READ else ""
        extra = {"summary": nums_summary, "sensitivity": nums_sensitivity, "methods": nums_methods}.get(s.id)
        if extra:
            try:
                guide += extra(ctx)
            except Exception as exc:
                guide += numbers(str(missing(s.id, f"numbers could not be read ({type(exc).__name__}: {exc})")), "Numbers in this section")
        parts.append(f'<section id="{s.id}"><h2>{E(heading)}</h2>{lead}{guide}{body}{vis}{foot}</section>')
    nav = "".join(f'<a href="#{s.id}">{NAV[s.id]}</a>' for s in SECTIONS)
    return (f'<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>GIPPA report (prototype)</title><style>{CSS}</style>'
            f"<script>{plotly.offline.get_plotlyjs()}</script></head><body><nav><b>GIPPA report</b>{nav}</nav><main>"
            f'<h1>GIPPA report</h1>{"".join(parts)}</main><script>{tables_js}</script><script>{JS}</script></body></html>')  # fmt: skip


if __name__ == "__main__":
    out_dir, target = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(build(out_dir), encoding="utf-8")
    print(target, target.stat().st_size)
    print(f"not available: {len(NOT_AVAILABLE)}", file=sys.stderr)
    for item in NOT_AVAILABLE:
        print(f"  {item}", file=sys.stderr)
