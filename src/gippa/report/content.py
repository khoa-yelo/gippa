# ruff: noqa: RUF001  (the report's wording is Unicode)
"""Report content beyond the M9-02 previews (M9-16; GIPPA_SPEC.md §10.3 sections 11-13): the
mandatory caveats and small tables of each section, the sensitivity tables by family, the
provenance block and the methods snippets; the plot guides (M12-05c) and the glossary (M12-04
PR 2).

Pure apart from the reads a ``RunView`` does and the package-data snippets: a view in, template
data out. A caveat that an ADR or the spec mandates keeps its wording, with the internal
reference replaced by a link to the docs-site methods page (``methods_url``).
"""

from __future__ import annotations

import dataclasses
import math
import re
from collections.abc import Collection, Mapping, Sequence
from importlib import resources
from pathlib import PurePosixPath
from typing import Any

import pandas as pd
from markupsafe import Markup

import gippa
import gippa.determinants
from gippa.context.diagnostics import render_l_hat
from gippa.data import calibration
from gippa.io.run_view import RunView
from gippa.io.tsv import format_value
from gippa.model.genetic_code import codon_table, table_name
from gippa.viz import function_terms
from gippa.viz.style import fmt_num, fmt_q
from gippa.viz.words import plain

#: GIPPA_SPEC.md §7.4, verbatim
SIGMA70 = (
    "Promotech and Promoter Calculator see σ70 promoters only, and TerminatorNet and TransTermHP "
    "see intrinsic terminators only. A missing call therefore never proves a missing element."
)
#: the docs-site methods pages the report links (Q12): a link, never a load
METHODS_BASE = "https://khoa-yelo.github.io/gippa/methods/"


def methods_url(topic: str) -> str:
    """The absolute docs-site URL of the methods page ``docs/methods/<topic>.md``."""
    return f"{METHODS_BASE}{topic}/"


#: ADR-0032, ADR-0033, ADR-0034 (the references are the link to the determinants methods page)
DETERMINANTS = (
    "The univariate Spearman ρ and AUROC CIs (TU-block bootstrap) are uncalibrated.",
    "The elastic-net ΔR² is the primary block attribution; LightGBM ΔR² and SHAP shares are "
    "secondary, nonlinear checks.",
    "The permutation null p is uncalibrated: it is valid only for exchangeable TUs and is not a "
    "discovery claim.",
)
#: ADR-0026, ADR-0027
GSEA = (
    "GSEA results are conditional on detection, with no adjustment for detection probability.",
    "NES is shown without an interval.",
)
#: ADR-0025
GSEA_CLUSTER = (
    "P-values are valid under TU exchangeability, and there is no level claim for sets whose "
    "genes cluster along the genome."
)
#: docs/methods/islands.md, "Validity domain (ADR-0021 Amendment 1)"
ISLANDS = "The FWER is claimed for a dependence regime, not for a replicon size."
KEGG_TERMS = "KEGG data: academic use only (rest.kegg.jp terms)."
#: the function section and provenance when the run had no KEGG (M12-03, ``gippa setup --no-kegg``)
KEGG_SKIPPED = "KEGG not installed: pathway enrichment skipped."
KEGG_ABSENT = "KEGG not installed"
NOT_RECORDED_RUN = "not recorded in this run"
#: ADR-0061 decision 1 (rows flagged ``tu_source_distance_model``), verbatim
MORAN_DISTANCE = (
    "Moran's I under the TU-block null is anti-conservative at short lags when TUs come from the "
    "distance model; rejection at α = 0.05 was 0.34 in simulation"
)
DISCORDANCE_DISTANCE = (
    "conservative under the distance-model fallback (merged TU boundaries widen the "
    "within-operon fit); fewer discordant pairs are found"
)
#: ADR-0061 decision 4 (rows flagged ``tu_block_tu_error``, every run), verbatim
MORAN_TU_ERROR = (
    "Moran's I under the TU-block null assumes exact TU boundaries and is anti-conservative at "
    "short lags when TUs are split; rejection at α = 0.05 was 0.67 in simulation with UniOP's "
    "K-12 error rates"
)
#: the row clause of an autocorrelation caveat (review r1 L1)
ROWS_CLAUSE = "(these rows carry the flag ‘{}’ in the table)"
DISTANCE_MODEL = "distance_model"
NOT_WRITTEN = "run_manifest.json not written by this GIPPA version."
NOT_RECORDED = "not recorded"
MORAN_LABEL = "lag-1 Moran's I of the detrended TU scores (used to choose the block length)"
_DIAG_HEADER = ["replicon", "ℓ̂ (l_hat)", "ℓ′ (block_len)", MORAN_LABEL, "block rule"]
_BIN_COLUMNS = ("bin", "lower_bp", "upper_bp", "N_NO", "N_cv", "N_sw")
#: the decoy counts by their plain names (behaviour 2)
_BIN_HEADER = ["bin", "lower_bp", "upper_bp", "non-operon pairs", "convergent pairs",
               "sandwich pairs"]  # fmt: skip
#: the Genome context note of a run with the mobile module off (lead decision 1)
MOBILE_OFF = "Mobile elements: prophage and plasmid regions not searched (the mobile module was "\
    "not run)."  # fmt: skip


def _p(text: str, topic: str | None = None, cls: str = "caveat") -> dict[str, Any]:
    """A paragraph; ``topic`` adds a link to its docs-site methods page."""
    return {"text": text, "cls": cls, "link": None if topic is None else methods_url(topic)}


def _bound(v: Any) -> str:
    if isinstance(v, int | float) and math.isinf(v):
        return "∞" if v > 0 else "−∞"
    return format_value(v)


def read_manifest(view: RunView) -> tuple[Mapping[str, Any] | None, str | None]:
    """The run manifest and, when it exists but cannot be read, a note naming it."""
    try:
        return view.manifest(), None
    except Exception as exc:
        return None, f"run_manifest.json could not be read: {type(exc).__name__}"


#: the banner of a ``gippa demo`` report (M11-15 item 2)
DEMO_BANNER = "Synthetic demo data (SIM100K, gippasim seed 20261014)"


def demo_banner(doc: Mapping[str, Any] | None) -> str | None:
    """``DEMO_BANNER`` when the run manifest says ``demo: true``, else None."""
    return DEMO_BANNER if isinstance(doc, Mapping) and doc.get("demo") is True else None


#: preflight rules whose ``warn`` opens the summary section (M11-14w; Q3: this one only in 1.0)
SUMMARY_WARNINGS = ("genome.genetic_code",)


def summary_warnings(view: RunView) -> list[dict[str, str]]:
    """The ``rule_id`` and ``message`` of each ``SUMMARY_WARNINGS`` preflight check at level
    ``warn`` in the run's ``preflight.json`` (M11-14w); empty when it is absent or unreadable."""
    try:
        doc = view.preflight()
    except Exception:
        return []
    checks = (doc or {}).get("checks") if isinstance(doc, Mapping) else None
    return [{"rule_id": str(c["rule_id"]), "message": str(c.get("message") or "")}
            for c in _list(checks)
            if isinstance(c, Mapping) and c.get("rule_id") in SUMMARY_WARNINGS
            and c.get("level") == "warn"]  # fmt: skip


def genetic_code_methods(table: int) -> str:
    """The Methods sentence of every run (M11-16 item 7)."""
    return (
        f"Genes were called and all six frames translated with NCBI genetic code {table} "
        f"({table_name(table)})."
    )


def genetic_code_note(table: int) -> dict[str, str] | None:
    """The summary note of a run at table 4 or 25 (``None`` at 11; M11-16 item 7, D8)."""
    if table == 11:
        return None
    tga = codon_table(table)["TGA"]
    text = f"Genetic code {table} ({table_name(table)}): TGA is read as {tga}."
    return {"table": str(table), "text": text}


def _list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _databases(doc: Mapping[str, Any] | None) -> list[Mapping[str, Any]]:
    return [d for d in _list((doc or {}).get("databases")) if isinstance(d, Mapping)]


def kegg_date(doc: Mapping[str, Any] | None) -> str | None:
    """The ``retrieved`` dates of the manifest's ``kegg*`` databases (lead decisions 2, 8)."""
    kegg = [d for d in _databases(doc) if str(d.get("name", "")).startswith("kegg")]
    dates = sorted({str(d["retrieved"]) for d in kegg if d.get("retrieved")})
    return ", ".join(dates) or None


def _safe(view: RunView, key: str) -> pd.DataFrame | None:
    try:
        return view.table(key)
    except Exception:  # the table preview of the section already carries the note
        return None


#: the second sentence of the no-terms notice (M12-04 PR 8)
NO_TERMS_LIGHT = ("With the Bakta light database this is expected: rerun with the full database "
                  "to get function results.")  # fmt: skip
NO_TERMS_CHECK = "Check the annotation source."


def function_no_terms(view: RunView) -> bool:
    """Whether the Function section is the no-terms notice alone (M12-04 PR 8)."""
    return function_terms.no_terms(function_terms.collections(view))


def _function(view: RunView, doc: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    flags = view.run_flags()
    colls = function_terms.collections(view)
    if colls is not None and function_terms.no_terms(colls):
        second = NO_TERMS_LIGHT if "bakta_db_light" in flags else NO_TERMS_CHECK
        return [_p(f"{function_terms.no_terms_sentence(colls)} {second}", cls="notice")]
    gsea = _safe(view, "enrichment_gsea") if view.has("enrichment_gsea") else None
    if view.has("enrichment_gsea"):
        out += [_p(t, "function-gsea") for t in GSEA]
        if gsea is not None and gsea["flags"].astype(str).str.contains("genomic_cluster").any():
            out.append(_p(GSEA_CLUSTER, "function-gsea"))
    no_kegg = "layer_missing:kegg" in flags
    cats = _safe(view, "category_summary") if view.has("category_summary") else None
    if not view.has("category_summary"):
        return [*out, _p(KEGG_SKIPPED)] if out and no_kegg else out
    cov = "" if cats is None else ", ".join(
        f"{c} {format_value(v)}"
        for c, v in cats.drop_duplicates("collection")[["collection", "coverage"]].itertuples(False)
    )  # fmt: skip
    if colls is not None:  # the event first, then category_summary (M12-04 PR 8)
        cov = function_terms.coverage_note(colls)
    out.append(_p(f"Annotation coverage leans toward well-studied genes; coverage per "
                  f"collection: {cov or NOT_RECORDED}.", "function"))  # fmt: skip
    if "bakta_db_light" in flags:
        out.append(_p("Annotation used the Bakta light DB, which assigns fewer functions than "
                      "the full DB."))  # fmt: skip
    if no_kegg:
        out.append(_p(KEGG_SKIPPED))
    elif colls is None or function_terms.kegg_terms(colls):
        date = kegg_date(doc) or NOT_RECORDED_RUN
        out.append(_p(f"{KEGG_TERMS} KEGG retrieval date: {date}."))
    return out


def _diagnostics(view: RunView) -> list[dict[str, Any]]:
    if not view.has("island_diagnostics"):
        if not view.has("islands"):
            return []
        return [_p("island_diagnostics.parquet: not produced in this run, so the block lengths "
                   "and the lag-1 Moran's I per replicon are not shown.")]  # fmt: skip
    diag = _safe(view, "island_diagnostics")
    if diag is None:
        return [_p("island_diagnostics.parquet: could not be read.")]
    rows = [[str(r["contig"]), render_l_hat(r["l_hat"], str(r["block_rule"])),
             format_value(r["block_len"]), fmt_num(r["morans_i_lag1"]), str(r["block_rule"])]
            for r in diag.to_dict("records")]  # fmt: skip
    return [{"key": "island_diagnostics", "header": _DIAG_HEADER, "rows": rows,
             "caption": "Island diagnostics per replicon",
             "link": methods_url("islands")}]  # fmt: skip


def _missed(view: RunView) -> list[dict[str, Any]]:
    try:
        found = view.events("operons.missed_operons")
    except Exception:
        return [_p("operons.events.jsonl: could not be read.")]
    if not found:
        return []
    data = found[-1]
    bins = [b for b in data.get("bin_counts", []) if isinstance(b, Mapping)]
    rows = [[(_bound if c.endswith("_bp") else format_value)(b.get(c)) for c in _BIN_COLUMNS]
            for b in bins]  # fmt: skip
    caption = f"Missed-operon decoy counts per distance bin (π₀ = {format_value(data.get('pi0'))})"
    return [{"key": "missed_operons_bins", "header": list(_BIN_HEADER), "rows": rows,
             "caption": caption, "link": methods_url("operons")}]  # fmt: skip


def _num(value: Any) -> str:
    return f"{value:.3f}" if isinstance(value, int | float) else NOT_RECORDED


def operon_source_sentence(reason: str | None) -> str:
    """The M11-14d summary and Methods sentence of a distance-model run; the operating point
    (14s2) is read from ``calibration.yaml`` at render time."""
    op = calibration.load().get("operon_distance") or {}

    def ci(key: str) -> str:
        return f"{_num(op.get(key))} [{_num(op.get(key + '_ci_low'))}, " \
               f"{_num(op.get(key + '_ci_high'))}]"  # fmt: skip

    return (
        f"Operons/TUs: intergenic-distance model ({reason or NOT_RECORDED}); "
        f"precision {ci('precision_at_pair_threshold')}, recall {ci('recall_at_pair_threshold')}"
        " on RegulonDB K-12 at the pair threshold; TU-based tests use these units."
    )


def _tu_block_morans(view: RunView) -> bool:
    if not view.has("autocorrelation"):
        return False
    ac = _safe(view, "autocorrelation")
    if ac is None or not {"statistic", "null_model"} <= set(ac.columns):
        return False
    return bool(((ac["statistic"] == "morans_i") & (ac["null_model"] == "tu_block")).any())


def extras(view: RunView, doc: Mapping[str, Any] | None) -> dict[str, list[dict[str, Any]]]:
    """Section id -> its caveat paragraphs (``text``) and small tables (``key``, ``header``,
    ``rows``, ``caption``), each present when its condition holds (M9-16 behaviour 3, 4), plus
    the operon-source sentence and the ADR-0061 caveats (M11-14d)."""
    out: dict[str, list[dict[str, Any]]] = {}
    distance = (doc or {}).get("operon_source") == DISTANCE_MODEL
    if distance:
        sentence = _p(operon_source_sentence((doc or {}).get("operon_source_reason")))
        out["summary"] = [sentence]
        out["methods"] = [sentence]
    if view.has("pairs") or view.has("operon_flags"):
        out["operons"] = [_p(SIGMA70, "operons"), *_missed(view)]
        if distance and view.has("pairs"):
            text = "Discordance tests of pairs with both genes detected are"
            out["operons"].append(_p(f"{text} {DISCORDANCE_DISTANCE}.", "operons"))
    if model_fitted(_safe(view, "determinants_model") if view.has("determinants_model") else None):
        out["determinants"] = [_p(t) for t in gippa.determinants.CAVEATS]
        out["determinants"] += [_p(t, "determinants") for t in DETERMINANTS]
    out["function"] = _function(view, doc)
    ctx = [_p(ISLANDS, "islands")] if view.has("islands") else []
    if _tu_block_morans(view):
        rows = ROWS_CLAUSE.format(plain("tu_block_tu_error"))
        ctx.append(_p(f"{MORAN_TU_ERROR} {rows}.", "autocorrelation"))
        if distance:
            rows = ROWS_CLAUSE.format(plain("tu_source_distance_model"))
            ctx.append(_p(f"{MORAN_DISTANCE} {rows}.", "autocorrelation"))
    if "mobile" in view.modules_off():
        ctx.append(_p(MOBILE_OFF, cls="note"))
    out["genome-context"] = ctx + _diagnostics(view)
    return out


def sensitivity_groups(frame: pd.DataFrame) -> list[tuple[str, list[list[str]]]]:
    """``sensitivity.tsv`` rows by ``family`` in first-seen order, cells formatted, ``jaccard``
    as ``%.2f``, ``NA`` or "no flagged set in either run" (M9-16 behaviour 5)."""
    groups: dict[str, list[list[str]]] = {}
    for rec in frame.to_dict("records"):
        cells = []
        for col in frame.columns:
            v = rec[col]
            if col == "jaccard":
                both0 = _is_zero(rec.get("n_primary")) and _is_zero(rec.get("n_variant"))
                missing = v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))
                v = "no flagged set in either run" if both0 else "NA" if missing else f"{v:.2f}"
            cells.append(format_value(v))
        groups.setdefault(str(rec["family"]), []).append(cells)
    return list(groups.items())


def _is_zero(value: Any) -> bool:
    """Whether ``value`` is a number equal to 0 (``None``, NA and text are not)."""
    try:
        return float(value) == 0.0
    except (TypeError, ValueError):
        return False


def _tools(doc: Mapping[str, Any]) -> list[list[str]]:
    tools = [t for t in _list(doc.get("tools")) if isinstance(t, Mapping)]
    return [[str(t.get("name", "")), str(t.get("version") or NOT_RECORDED),
             str(t.get("container_digest") or t.get("conda_lock_hash") or NOT_RECORDED)]
            for t in tools]  # fmt: skip


#: a path attached to a short option (``-o/x``)
_SHORT_PATH = re.compile(r"-[A-Za-z]/")


def _local_path(arg: str, outdir: str, redact: bool = False) -> str:
    """``arg`` without an absolute local path (M9-02 behaviour 6, §6.6; #323 B1): a path under
    the run's ``outdir`` as ``OUT/<relative>``, any other absolute path as its file name (with
    ``redact``: ``<outdir>/<relative>`` and ``<input>/<file name>``, M12-04 PR 13). The path may
    be the whole token, follow its last ``=``, be an item of a comma list, or be attached to a
    short option (``-o/x``); relative parts are unchanged."""
    if arg.startswith("/"):
        head, value = "", arg
    elif _SHORT_PATH.match(arg):
        head, value = arg[:2], arg[2:]
    else:
        cut = arg.rfind("=") + 1
        head, value = arg[:cut], arg[cut:]
    return head + ",".join(_shown_path(item, outdir, redact) for item in value.split(","))


def _shown_path(value: str, outdir: str, redact: bool = False) -> str:
    """One path item: ``OUT/<relative>`` under ``outdir``, else the file name if absolute;
    with ``redact`` ``<outdir>/<relative>`` and ``<input>/<file name>``."""
    if not value.startswith("/"):
        return value
    path, out = PurePosixPath(value), PurePosixPath(outdir) if outdir.startswith("/") else None
    if out is not None and (path == out or out in path.parents):
        return str(PurePosixPath("<outdir>" if redact else "OUT", path.relative_to(out)))
    if redact:
        return str(PurePosixPath("<input>", path.name))
    return path.name or value


#: a URL (left unchanged) or a path: a ``/`` at the start or after ``=``, ``,``, ``[``, ``{``,
#: whitespace, a quote or ``:``, running to the next ``,``, ``]``, ``}``, quote or whitespace
_STOP = r",\]}'\"\s"
_URL = rf"[A-Za-z][A-Za-z0-9+.\-]*://[^{_STOP}]*"
_ARG_PATH = re.compile(rf"(?P<url>{_URL})|(?:^|(?<=[=,\[{{\s'\":]))/[^{_STOP}]*")
#: the same in prose (an unreadable-table note): also after ``(``, also ending at ``)``
_NOTE_PATH = re.compile(rf"(?P<url>{_URL})|(?:^|(?<=[=,\[{{(\s'\":]))/[^{_STOP})]*")


def _redact_paths(text: str, outdir: str, pattern: re.Pattern[str] = _ARG_PATH) -> str:
    """``text`` with every path of ``pattern`` redacted (``_shown_path``); URLs unchanged."""
    return pattern.sub(lambda m: m[0] if m["url"] else _shown_path(m[0], outdir, True), text)


def redact_note(text: str, outdir: str) -> str:
    """An unreadable-table note for ``gippa report --redact``: its absolute paths rewritten as
    in ``redact_command_line`` (a path may also follow ``(`` and end at ``)``)."""
    return _redact_paths(text, outdir, _NOTE_PATH)


def redact_command_line(args: Sequence[str], outdir: str) -> str:
    """``args`` joined by single spaces for a shareable report (``gippa report --redact``,
    R-09): every absolute path as ``<input>/<file name>``, a path under ``outdir`` as
    ``<outdir>/<relative>`` (``<outdir>`` itself), the ``--profile`` value as ``<profile>``
    (next token or after ``--profile=``); relative arguments and URLs are unchanged. A path
    starts a token, follows ``=`` (``--set=k=/a=b/c`` is the path ``/a=b/c``), ``,``, ``[``,
    ``{``, a space, a quote or ``:``, or is attached to a short option (``-o/x``); it runs to
    the next ``,``, ``]``, ``}``, quote or whitespace (review r1 M1)."""
    shown, prev = [], ""
    for arg in map(str, args):
        if prev == "--profile":
            shown.append("<profile>")
        elif arg.startswith("--profile="):
            shown.append("--profile=<profile>")
        else:
            cut = 2 if _SHORT_PATH.match(arg) else 0
            shown.append(arg[:cut] + _redact_paths(arg[cut:], outdir))
        prev = arg
    return " ".join(shown)


def _command_line(doc: Mapping[str, Any], redact: bool = False) -> str:
    """The manifest's ``command_line`` joined by single spaces, absolute paths made local
    (``redact``: ``redact_command_line``)."""
    outdir = str(doc.get("outdir") or "")
    if redact:
        return redact_command_line([str(a) for a in doc.get("command_line") or []], outdir)
    return " ".join(_local_path(str(a), outdir) for a in doc.get("command_line") or [])


def licence_notes(doc: Mapping[str, Any] | None) -> list[str]:
    """The run manifest's ``licence_notes`` (M11-15c2), for the Methods section."""
    return [str(n) for n in _list((doc or {}).get("licence_notes")) if isinstance(n, str)]


def provenance(
    view: RunView, doc: Mapping[str, Any] | None, note: str | None, redact: bool = False
) -> dict[str, Any]:
    """Provenance (M9-16 behaviour 6; lead decisions 1, 7, 8): from ``run_manifest.json`` when it
    is readable, else the pass-1 subset (preflight SHA-256s) with a note; ``redact`` rewrites
    the command line only (M12-04 PR 13)."""
    rows = doc.get("inputs") if doc is not None else (view.preflight() or {}).get("inputs")
    inputs = [[str(i.get("file", "")).replace("\\", "/").rsplit("/", 1)[-1],
               str(i.get("role", "")), str(i.get("sha256", ""))]
              for i in _list(rows) if isinstance(i, Mapping)]  # fmt: skip
    base = {
        "version": gippa.__version__,
        "inputs": inputs,
        "config": view.digest("resolved_config"),
    }
    if doc is None:
        return {**base, "note": note or NOT_WRITTEN}
    seeds = doc.get("seeds")
    dbs = [[str(d.get("name", "")), str(d.get("version", "")), str(d.get("checksum", "")),
            str(d.get("retrieved", ""))] for d in _databases(doc)]  # fmt: skip
    go = [str(d.get("version")) for d in _databases(doc) if d.get("name") == "go_basic"]
    dirty = " (uncommitted changes)" if doc.get("git_dirty") else ""
    # M12-03: the flag wins over a lone KEGG item's date (a partial KEGG install is not used)
    kegg_line = (KEGG_ABSENT if "layer_missing:kegg" in view.run_flags()
                 else kegg_date(doc) or NOT_RECORDED_RUN)  # fmt: skip
    lines = [
        ("GIPPA version", str(doc.get("gippa_version", NOT_RECORDED))),
        ("commit", f"{doc.get('git_commit') or NOT_RECORDED}{dirty}"),
        ("command line", _command_line(doc, redact) or NOT_RECORDED),
        ("seed", str(seeds.get("master")) if isinstance(seeds, Mapping) else NOT_RECORDED),
        ("KEGG retrieval date", kegg_line),
        ("GO release", ", ".join(go) or NOT_RECORDED),
        ("manifest digest (SHA-256 of run_manifest.json)", view.digest("run_manifest") or ""),
    ]
    if doc.get("provenance_digest"):
        lines.append(("provenance digest", str(doc["provenance_digest"])))
    return {**base, "lines": lines, "tools": _tools(doc), "databases": dbs}


#: methods topic -> (its heading, the steps whose run shows it, the tables that show it when the
#: manifest lists no steps); the report shows them in this order, after the run's tier (R-07)
METHOD_TOPICS: Mapping[str, tuple[str, tuple[str, ...], tuple[str, ...]]] = {
    "placement-fdr": ("Placement FDR", ("placement",), ("mapping_fdr",)),
    "denovo-precision": ("De novo precision", ("placement",), ("mapping_fdr",)),
    "abundance": ("Abundance", ("gene_table",), ("genes_tsv",)),
    "operons": ("Operons", ("operon_pairs", "operons"), ("pairs", "operons", "operon_flags")),
    "islands": ("Islands", ("context_stats",), ("islands",)),
    "variants": ("Variants", ("variants",), ("variants",)),
    "novel-orfs": ("Novel ORFs", ("novel_orfs",), ("novel_orfs",)),
    "function": ("Function enrichment", ("function",),
                 ("enrichment_detection", "enrichment_gsea", "category_summary")),
    "determinants": ("Determinants", ("determinants",), ("determinants_model",)),
}  # fmt: skip
#: methods topic -> the report section it belongs to (``gippa.report.build.SECTIONS``): a note
#: whose section is not applicable to the run is not shown (M12-04 PR 3 review r1 M1)
TOPIC_SECTION: Mapping[str, str] = {
    "placement-fdr": "mapping", "denovo-precision": "mapping", "abundance": "detection",
    "operons": "operons", "islands": "genome-context", "variants": "variants",
    "novel-orfs": "novel-orfs", "function": "function", "determinants": "determinants",
}  # fmt: skip
#: topics that need a finite ``abundance_log2`` (review r1 M1)
NEEDS_ABUNDANCE = frozenset({"abundance", "operons", "determinants"})
#: detection tier -> the heading of its methods snippet
TIER_TITLES: Mapping[str, str] = {
    "decoy_picked": "Detection tier: picked decoys",
    "mayu": "Detection tier: MAYU",
    "two_peptide": "Detection tier: two peptides",
    "protein_list": "Detection tier: protein list",
}


def _snippet(stem: str) -> Markup | None:
    """The package-data snippet ``<stem>.html`` (generated from ``docs/methods`` by
    scripts/dev/sync_report_snippets.py), or ``None``."""
    path = resources.files("gippa.report").joinpath("snippets", f"{stem}.html")
    return Markup(path.read_text(encoding="utf-8")) if path.is_file() else None


def _steps_run(doc: Mapping[str, Any] | None) -> frozenset[str] | None:
    """The steps the manifest lists as run (not failed); ``None`` without a ``steps`` list."""
    steps = doc.get("steps") if isinstance(doc, Mapping) else None
    if not isinstance(steps, list):
        return None
    return frozenset(str(s.get("name")) for s in steps
                     if isinstance(s, Mapping) and s.get("status") != "failed")  # fmt: skip


def _run_tier(view: RunView, doc: Mapping[str, Any] | None) -> str | None:
    """The manifest's ``detection_tier``, else the first one ``genes.tsv`` carries."""
    tier = doc.get("detection_tier") if isinstance(doc, Mapping) else None
    if tier is None:
        genes = _safe(view, "genes_tsv")
        if genes is not None and "detection_tier" in genes.columns:
            tier = next(iter(genes["detection_tier"].dropna()), None)
    return None if tier is None else str(tier)


def _has_denovo(view: RunView) -> bool:
    """Whether the run placed de novo peptides (or wrote de novo precision rows)."""
    hits, fdr = _safe(view, "peptide_hits"), _safe(view, "mapping_fdr")
    return bool((hits is not None and "source" in hits.columns
                 and (hits["source"] == "denovo").any())
                or (fdr is not None and "section" in fdr.columns
                    and (fdr["section"] == "denovo_precision").any()))  # fmt: skip


def _has_abundance(view: RunView) -> bool:
    genes = _safe(view, "genes_tsv")
    if genes is None or "abundance_log2" not in genes.columns:
        return False
    vals = pd.to_numeric(genes["abundance_log2"], errors="coerce")
    return bool(vals.map(lambda v: math.isfinite(float(v))).any())


def has_rows(view: RunView, key: str) -> bool:
    """Whether the run wrote table ``key`` with at least one row (an unreadable one: no)."""
    frame = _safe(view, key) if view.has(key) else None
    return frame is not None and len(frame) > 0


def methods_snippets(
    view: RunView, doc: Mapping[str, Any] | None, off: Collection[str] = ()
) -> list[dict[str, Any]]:
    """The Methods section's snippets (R-07): the run's tier (always), then each topic of
    ``METHOD_TOPICS`` whose step ran (the manifest's ``steps``, not failed; a run without them:
    one of whose tables has a row) and whose section (``TOPIC_SECTION``) is not in ``off``, the
    not-applicable section ids of ``gippa.report.build.not_applicable``. De novo precision also
    needs de novo peptides, and the ``NEEDS_ABUNDANCE`` topics a finite ``abundance_log2``. Each:
    ``topic`` (the snippet's file stem), ``title``, ``body`` and ``link`` (its methods page)."""
    ran = _steps_run(doc)
    tier = _run_tier(view, doc) or ""
    picked = [(f"detection-tiers.{tier}", TIER_TITLES[tier])] if tier in TIER_TITLES else []
    abundance = _has_abundance(view)
    for topic, (title, steps, tables) in METHOD_TOPICS.items():
        shown = (any(s in ran for s in steps) if ran is not None
                 else any(has_rows(view, t) for t in tables))  # fmt: skip
        if (not shown or TOPIC_SECTION.get(topic) in off
                or (topic in NEEDS_ABUNDANCE and not abundance)
                or (topic == "denovo-precision" and not _has_denovo(view))):  # fmt: skip
            continue
        picked.append((topic, title))
    out = []
    for stem, title in picked:
        body = _snippet(stem)
        if body is not None:
            topic = stem.split(".", 1)[0]
            out.append({"topic": stem, "title": title, "body": body, "link": methods_url(topic)})
    return out


def model_fitted(model: pd.DataFrame | None) -> bool:
    """Whether ``determinants_model`` has a finite ``cv_r2_mean`` (any model): the Determinants
    caveats print only then (behaviour 4)."""
    if model is None or not {"metric", "value"} <= set(model.columns):
        return False
    vals = pd.to_numeric(model.loc[model["metric"] == "cv_r2_mean", "value"], errors="coerce")
    return bool(vals.map(lambda v: math.isfinite(float(v))).any())


@dataclasses.dataclass(frozen=True)
class Guide:
    """How to read one plot (M12-05c): its GIPPA_SPEC.md §10.2 Question, one item per panel
    (lettered "(a) ", "(b) ", ... as the plot prints them), what a pattern means and, optionally,
    when not to trust it. The report renders it as a labelled ``<dl>`` in the figure caption."""

    question: str
    panels: tuple[str, ...]
    read: str
    caveats: str | None = None


#: plot ID -> its guide (M12-05c § 4); static text, the same in every report
GUIDES: Mapping[str, Guide] = {
    "V01": Guide(
        "Where in the genome is the proteome?",
        ("One circle; each replicon is a sector with ticks in kb. Numbered rings, outer to "
         "inner, as listed in the ring key: regions (prophage, plasmid, BGC, CRISPR and AMR "
         "arcs); genes on the + strand, then the − strand, coloured by detection status; log2 "
         "abundance as a heat ring (white: no value or imputed); peptide density per 10-kb bin; "
         "transcription-unit starts; promoter and terminator ticks; abundance islands; "
         "cumulative GC skew; variant sites and novel-ORF candidates. Radial lines mark the "
         "origin and terminus when they are known.",),
        "Arcs where the gene rings are mostly \"detected\" and the abundance and peptide rings "
        "are dense are the parts of the genome the proteome covers. A region arc over undetected "
        "genes and no peptides is a region the sample did not detect; V07 tests whether that is "
        "more than detectability predicts.",
        "A ring with nothing to draw is left out, and the title names the missing layers. The "
        "circle is descriptive; the tests are in V02, V07 and the Genome context tables.",
    ),
    "V02": Guide(
        "Is abundance clustered along the chromosome?",
        ("One panel per replicon: x is the position (kb), y the log2 abundance of each detected "
         "gene with a measured value (one point per gene); undetected genes are ticks along the "
         "bottom; the line is a LOESS trend; shaded bands are abundance islands; vertical lines "
         "mark the origin and terminus when they are known. With no measured abundance, each "
         "panel is one strip of ticks coloured by detection status.",),
        "A shaded band is a run of neighbouring genes more (or less) abundant than the rest of "
        "the replicon, beyond chance at a family-wise error of 0.05. The title gives the number "
        "of islands and the lag-1 Moran's I with its q. A wavy trend line without bands is not "
        "evidence of clustering.",
        "The trend line is a visual guide and is not tested. The island test shuffles whole "
        "transcription units, so operon partners looking alike do not by themselves make an "
        "island.",
    ),
    "V03": Guide(
        "What happens at this prophage, BGC, island, flagged operon or novel ORF?",
        ("One locus per figure, with flanks; x is the position (kb). Gene arrows point in the "
         "direction of transcription, coloured by log2 abundance on one scale for the whole run; "
         "undetected genes are white with a grey outline and hatching. Mapped peptides are "
         "stacked above the arrows. Promoters are bent arrows, terminators lollipops; the region "
         "or island is shaded. A dashed line marks the origin when the window crosses it.",),
        "The title gives the locus, its detected share against the share detectability "
        "predicts, and the q. Look at which genes carry peptides and how abundant they are "
        "against the flanks. The views follow the \"Locus views, in display order\" table: the "
        "best locus of each kind in turn, then the second best, and so on, at most 50.",
        "Novel-ORF candidates are drawn under Novel ORFs (V13), not here. A gene supported by a "
        "single peptide is marked, and its abundance is not estimated.",
    ),
    "V05": Guide(
        "Where did my sequences go?",
        ("One flow diagram, left to right: input class (peptide or protein), then the matching "
         "route (exact, tolerant, miniprot, tblastn or unmapped), then the outcome (placed, "
         "contaminant, rejected by FDR or unexplained). Band widths are proportional to the "
         "number of distinct sequences, and the labels give the counts. When every sequence "
         "takes one route to one outcome, the diagram is a single bar.",),
        "Most of the flow should end in \"placed\". A large unmapped, unexplained or rejected "
        "share means the genome, the search database and the input do not match well; V06 checks "
        "the tolerant route, and the unmapped table lists the sequences.",
        "Each sequence counts once: one placed by one route and unmapped by another counts as "
        "placed.",
    ),
    "V06": Guide(
        "Can I trust tolerant hits and de novo peptides?",
        ("(a) targets vs. decoys: target and decoy placements per peptide length class and "
         "number of differences, on a log scale (a zero count sits at the floor).",
         "(b) stratum FDR: the estimated false discovery rate per number of differences, with the "
         "acceptance level dashed; a stratum with no targets is NA.",
         "(c) de novo precision: the share of correct de novo peptides per score bin, with its "
         "confidence interval; drawn only for de novo input."),
        "Where decoys come close to targets in (a), those placements are mostly chance; (b) "
        "turns this into an FDR, and only strata at or below the dashed line are accepted. Bins "
        "with low precision in (c) add errors to every later count. When every peptide matched "
        "exactly, the plot is one line saying so.",
        "Without reference peptides, de novo precision is not measured and an assumed value is "
        "used; panel (c) names it.",
    ),
    "V07": Guide(
        "Which regions and categories are detected less than detectability predicts?",
        ("One row per group, in blocks: replicon, feature class, COG category, region. The point "
         "is observed minus expected detected genes in percentage points of the group's genes, "
         "the bar its 95 % confidence interval; the group size n is in the row label and the q "
         "at the right.",),
        "A row left of zero has fewer detected genes than its genes' detectability predicts, a "
        "row right of zero more. Only rows whose q is at or below the level in the title are "
        "findings; the title counts them.",
        "The expected counts come from the detectability model, so check its reliability in "
        "V18 (a) first: where the model is off, every group can drift.",
    ),
    "V08": Guide(
        "Which functions dominate abundance?",
        ("One column per COG category, highest median abundance first and \"none\" (no COG "
         "annotation) last: each measured gene's log2 abundance is a point, or a violin with its "
         "median above 50 genes; the label above each column is detected / total genes. With no "
         "measured abundance, each column shows the detected fraction instead.",),
        "Categories to the left hold the most abundant proteins. Compare the medians, not the "
        "number of points.",
        "Imputed values are not drawn. The plot is descriptive; the tested comparison of "
        "functions is V14.",
    ),
    "V10": Guide(
        "Is there structure beyond operons?",
        ("(a) Moran's I of abundance against the lag (genes apart, 1 to 20) per replicon, with "
         "two shaded 95 % null envelopes: gene shuffles (grey) and transcription-unit block "
         "shuffles (blue).",
         "(b) the same for detection: the join statistic of neighbouring pairs in which both "
         "genes are detected, as a Pearson residual, against the same two envelopes."),
        "A filled marker lies outside the blue envelope with q at or below the level: genes that "
        "far apart look alike beyond what operons explain. Open markers inside the blue envelope "
        "at every lag mean operons explain the similarity. The title gives the lag-1 result.",
        "The grey envelope ignores operons and is too narrow; judge by the blue one. Lags with "
        "too few gene pairs are not drawn.",
    ),
    "V11": Guide(
        "Do operon partners look alike, and which other neighbours do too?",
        ("(a) |Δ| by pair class: histograms of the absolute log2 abundance difference of "
         "adjacent genes, for same-operon (OP), non-operon (NO), convergent (CV) and sandwiched "
         "(SW) pairs, with n in the legend.",
         "(b) the local FDR of each pair against |Δ|, with the flag level dashed.",
         "(c) partners: upstream against downstream log2 abundance of OP and NO pairs, coloured "
         "by class; flagged pairs are larger outlined markers, one shape per flag."),
        "OP pairs at smaller |Δ| than NO pairs in (a) mean operon partners are made at similar "
        "levels, as expected. Pairs below the dashed line in (b) are flagged, for example an "
        "operon pair whose partners differ or a non-operon pair that behaves like one. In (c), "
        "points near the diagonal look alike.",
        "Pairs with an imputed abundance are left out of (a) and (c). A flag is a lead to check "
        "in the locus view (V03), not a confirmed promoter or terminator.",
    ),
    "V12": Guide(
        "Where do peptides support ALT residues?",
        ("(a) variant sites by call (ALT or mixed), one bar per confidence level, stacked by "
         "flag (a site with several flags counts in each); the total is printed above each bar.",
         "(b) alignments of the top ALT sites (at most 15): the reference protein window, then "
         "the covering peptide with the ALT residue in bold colour, and the number of covering "
         "peptides."),
        "ALT means the peptides show only the alternative residue; mixed means both residues are "
        "seen. High-confidence ALT sites without flags are the strongest. In (b), the peptide "
        "should match the reference everywhere except the highlighted residue.",
        "The class error rate is pooled over every site of a class, not per site; a rate of "
        "exactly 1 (marked *) means it was not measured. I and L have the same mass and are "
        "written I/L.",
    ),
    # V13 follows M12-05a: every candidate, class q first, then the best peptide q
    "V13": Guide(
        "Which new genes or upstream extensions are supported?",
        ("(a) candidates per class (for example N-terminal extension, intergenic, small ORF), "
         "stacked by status: called, not called or unestimated.",
         "(b) every candidate, called or not, smallest q first (the class FDR; ties go to the "
         "called candidate, then the smallest best peptide q; unestimated last), each in its "
         "genomic neighbourhood: the candidate (filled when called, a dashed outline when not), "
         "the nearby gene arrows, its peptides and the promoter and terminator marks, with its "
         "class, status and q in the title."),
        "A called candidate passed its class FDR; its q is the lowest FDR at which it would be "
        "called. Read (b) from the top: the first candidates are the best supported, and a "
        "candidate with several peptides, a promoter upstream and a start codon with a "
        "ribosome-binding site is the most convincing.",
        "When a class has too few decoys its q is unestimated and nothing in it is called; the "
        "candidates are still drawn. At most 20 candidates are drawn; the index title gives the "
        "total. A single peptide is weak evidence on its own.",
    ),
    "V14": Guide(
        "Which functions are high or low?",
        ("(a) GSEA running sums for the top six gene sets: the running enrichment score along "
         "the genes ranked by abundance, the set's genes as ticks and the leading edge shaded.",
         "(b) each set's normalised enrichment score (NES) against −log10 q, coloured by "
         "collection; flagged sets are hollow, redundant sets small and grey, and the ten most "
         "significant are labelled.",
         "(c) detection over-representation: observed over expected detected genes per category "
         "on a log2 axis, coloured by direction, with −log10 q in a strip beside it and the level "
         "line."),
        "In (b), sets to the right are more abundant than average and sets to the left less; "
        "higher means a smaller q. A running sum in (a) that peaks early is driven by the most "
        "abundant genes. In (c), categories to the right of 1 are detected more often than "
        "expected.",
        "A set marked \"not calibrated\" (dashed outline in (a), hollow in (b)) is clustered in "
        "the genome or tied to detection, so its q is too small. Panel (c) has no confidence "
        "intervals.",
    ),
    "V15": Guide(
        "What explains abundance?",
        ("(a) block ΔR²: the drop in cross-validated R² when each feature block (transcription, "
         "translation, codon, protein, context) is left out, for the elastic net (the main "
         "model) and LightGBM (a check), with intervals.",
         "(b) partial dependence of abundance on the codon adaptation index, the tRNA adaptation "
         "index, the OSTIR initiation rate and the promoter score.",
         "(c) SHAP values of the top 20 features: one point per gene, coloured by the feature's "
         "value.",
         "(d) observed against out-of-fold predicted log2 abundance; genes more than 2 SD off "
         "are marked as unexplained."),
        "A large ΔR² in (a) means that block carries information the others lack. (b) and (c) "
        "show the direction, for example higher codon adaptation going with higher abundance. In "
        "(d), points near the diagonal are well predicted; the title gives the cross-validated "
        "R².",
        "When the model does not predict abundance (R² at or below 0, or its interval includes "
        "0), only (a) and (d) are drawn, under a banner. Feature blocks are correlated, and the "
        "intervals are uncalibrated because the cross-validation repeats share data.",
    ),
    "V16": Guide(
        "How do search strategies differ?",
        ("(a) genes detected per exact set of runs: each bar counts the genes detected by "
         "exactly the runs marked by the dots below it; each run's total is at the left; the 30 "
         "largest sets are drawn and the rest summed as \"other\".",
         "(b) pairwise concordance: below the diagonal, Spearman ρ of the shared measured log2 "
         "abundances with its 95 % CI; above it, the shared values of each pair against the "
         "identity line."),
        "A tall bar under every dot means the runs agree on those genes; a bar under a single "
        "dot counts genes only that run finds. In (b), ρ near 1 means the runs rank abundance "
        "alike; high ρ with points off the identity line is a difference of scale only.",
        "Panel (b) uses only genes measured in both runs. Runs with different versions may not "
        "be comparable; the Runs table says when.",
    ),
    "V17": Guide(
        "Are the inputs sane?",
        ("(a) preflight checks: counts by level (ok, info, warn, error).",
         "(b) contig lengths.",
         "(c) distinct peptide lengths (aa), placed and unmapped.",
         "(d) log2 abundance, measured and imputed.",
         "(e) missing values (%) per gene feature over complete genes, with the reason when a "
         "feature was not computed."),
        "An error in (a) needs action and a warning a look; the list under the panels names "
        "them. Unmapped peptides much longer or shorter than the placed ones point to a wrong "
        "digest or database. A large imputed share in (d), or a feature near 100 % missing in "
        "(e), weakens every later plot that uses it.",
    ),
    "V18": Guide(
        "Are the tests calibrated on this dataset?",
        ("(a) detectability reliability: the observed detection rate against the mean predicted "
         "probability per bin, with the slope, intercept and Brier score.",
         "(b) p-value histograms, one per test family, with n; families with fewer than 10 tests "
         "are listed, not drawn.",
         "(c) null distributions for regions, islands and the determinants model, with the "
         "observed statistic marked and its p."),
        "In (a), points on the diagonal (slope near 1, intercept near 0) mean the detectability "
        "model predicts well, so V07 can be trusted. In (b), a calibrated test is flat away from "
        "zero; a spike near zero is signal, a pile-up near 1 a conservative test. In (c), an "
        "observed value far in the tail is unlikely under the null.",
        "Gene-shuffle nulls ignore operons and are too narrow; judge by the transcription-unit "
        "block null.",
    ),
}  # fmt: skip

#: plot ID -> the stem of its index instance and that instance's How to read (review r1 L2;
#: M12-05c § 4 V13); the V03 index is a table, not an image (M12-04 PR 4, P-12)
INDEX_STEMS: Mapping[str, str] = {"V13": "V13_novel_orfs_locus"}
INDEX_READS: Mapping[str, str] = {
    "V13": "This index lists every candidate drawn in (b), in the same order.",
}
#: plot ID -> what a short guide calls its first figure ("As for the first <name> ...")
GUIDE_NAMES: Mapping[str, str] = {"V03": "locus view", "V13": "novel ORF figure"}


def short_read(pid: str, stem: str) -> str | None:
    """The How to read of a later instance's short guide: the index instance's own sentence, or
    ``None`` for the "As for the first ..." pointer to the full guide (M12-05c decision 3)."""
    return INDEX_READS[pid] if INDEX_STEMS.get(pid) == stem else None


def guide_name(pid: str) -> str:
    """The name a short guide of ``pid`` gives its first figure."""
    return GUIDE_NAMES.get(pid, "figure of this plot")


#: plot ID -> (the stem of its image-free index table, the table's caption)
INDEX_TABLES: Mapping[str, tuple[str, str]] = {
    "V03": ("V03_locus", "Locus views, in display order"),
}
LOCUS_INDEX_HEADER = ("rank", "locus", "family", "estimate", "q")


def locus_index(header: Sequence[str], rows: Sequence[Sequence[str | None]]) -> list[list[str]]:
    """The report rows of ``data/V03_locus.tsv``: rank, locus, family (in words), the name of the
    family's estimate and its value with 2 significant digits (``fmt_q``), in TSV order."""
    col = {name: i for i, name in enumerate(header)}

    def cell(row: Sequence[str | None], name: str) -> str | None:
        return row[col[name]] if name in col else None

    out = []
    for r in rows:
        family = cell(r, "family") or cell(r, "kind") or ""
        out.append([cell(r, "rank") or "", cell(r, "locus_id") or "",
                     family.replace("_", " "), cell(r, "estimate") or "q",
                     fmt_q(cell(r, "error"))])  # fmt: skip
    return out


def question(pid: str) -> str:
    """The GIPPA_SPEC.md §10.2 Question of plot ``pid`` (``GUIDES[pid].question``)."""
    return GUIDES[pid].question


#: the Methods glossary (behaviour 2): term -> plain definition
GLOSSARY: tuple[tuple[str, str], ...] = (
    ("FDR", "False discovery rate: the expected share of false findings among the findings "
     "reported."),
    ("FWER", "Family-wise error rate: the probability of at least one false finding in a whole "
     "family of tests."),
    ("q-value", "The smallest FDR at which a finding would be reported."),
    ("local FDR (lfdr)", "The probability that one particular finding is false, given its "
     "value."),
    ("O − E", "Observed minus expected: how many more (or fewer) genes were detected than the "
     "detectability model predicts."),
    ("π₀", "The estimated share of tested items with no true effect."),
    ("block length (ℓ′)", "The number of consecutive transcription units resampled together in "
     "the block bootstrap, chosen from the estimated correlation length ℓ̂ as min(20, 2ℓ̂ + 1, "
     "number of TUs), counted per replicon."),
    ("BB join count", "The number of neighbouring gene pairs where both genes are detected, "
     "compared with its value under random placement."),
)  # fmt: skip
