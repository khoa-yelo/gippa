# ruff: noqa: RUF001  (the report's wording is Unicode)
"""Report content beyond the M9-02 previews (M9-16; GIPPA_SPEC.md §10.3 sections 11-13): the
mandatory caveats and small tables of each section, the sensitivity tables by family, the
provenance block and the methods snippets; the figure captions and the glossary (M12-04 PR 2).

Pure apart from the reads a ``RunView`` does and the package-data snippets: a view in, template
data out. A caveat that an ADR or the spec mandates keeps its wording, with the internal
reference replaced by a link to the docs-site methods page (``methods_url``).
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
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
from gippa.viz.style import fmt_num
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


def _function(view: RunView, doc: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    flags = view.run_flags()
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
    out.append(_p(f"Annotation coverage leans toward well-studied genes; coverage per "
                  f"collection: {cov or NOT_RECORDED}.", "function"))  # fmt: skip
    if "bakta_db_light" in flags:
        out.append(_p("Annotation used the Bakta light DB, which assigns fewer functions than "
                      "the full DB."))  # fmt: skip
    if no_kegg:
        out.append(_p(KEGG_SKIPPED))
    else:
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


def _local_path(arg: str, outdir: str) -> str:
    """``arg`` without an absolute local path (M9-02 behaviour 6, §6.6; #323 B1): a path under
    the run's ``outdir`` as ``OUT/<relative>``, any other absolute path as its file name. The
    path may be the whole token, follow its last ``=``, be an item of a comma list, or be
    attached to a short option (``-o/x``); relative parts are unchanged."""
    if arg.startswith("/"):
        head, value = "", arg
    elif _SHORT_PATH.match(arg):
        head, value = arg[:2], arg[2:]
    else:
        cut = arg.rfind("=") + 1
        head, value = arg[:cut], arg[cut:]
    return head + ",".join(_shown_path(item, outdir) for item in value.split(","))


def _shown_path(value: str, outdir: str) -> str:
    """One path item: ``OUT/<relative>`` under ``outdir``, else the file name if absolute."""
    if not value.startswith("/"):
        return value
    path, out = PurePosixPath(value), PurePosixPath(outdir) if outdir.startswith("/") else None
    if out is not None and (path == out or out in path.parents):
        return str(PurePosixPath("OUT", path.relative_to(out)))
    return path.name or value


def _command_line(doc: Mapping[str, Any]) -> str:
    """The manifest's ``command_line`` joined by single spaces, absolute paths made local."""
    outdir = str(doc.get("outdir") or "")
    return " ".join(_local_path(str(a), outdir) for a in doc.get("command_line") or [])


def licence_notes(doc: Mapping[str, Any] | None) -> list[str]:
    """The run manifest's ``licence_notes`` (M11-15c2), for the Methods section."""
    return [str(n) for n in _list((doc or {}).get("licence_notes")) if isinstance(n, str)]


def provenance(view: RunView, doc: Mapping[str, Any] | None, note: str | None) -> dict[str, Any]:
    """Provenance (M9-16 behaviour 6; lead decisions 1, 7, 8): from ``run_manifest.json`` when it
    is readable, else the pass-1 subset (preflight SHA-256s) with a note."""
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
        ("command line", _command_line(doc) or NOT_RECORDED),
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


def methods_snippets(view: RunView, doc: Mapping[str, Any] | None) -> list[dict[str, Any]]:
    """The Methods section's snippets (R-07): the run's tier, then each topic of
    ``METHOD_TOPICS`` whose step ran (the manifest's ``steps``; a run without them: whose table
    is present). De novo precision also needs de novo peptides and abundance a finite
    ``abundance_log2``. Each: ``topic`` (the snippet's file stem), ``title``, ``body`` and
    ``link`` (its methods page)."""
    ran = _steps_run(doc)
    tier = _run_tier(view, doc) or ""
    picked = [(f"detection-tiers.{tier}", TIER_TITLES[tier])] if tier in TIER_TITLES else []
    checks = {"denovo-precision": _has_denovo, "abundance": _has_abundance}
    for topic, (title, steps, tables) in METHOD_TOPICS.items():
        shown = (any(s in ran for s in steps) if ran is not None
                 else any(view.has(t) for t in tables))  # fmt: skip
        if shown and checks.get(topic, lambda _v: True)(view):
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


#: plot ID -> its caption: the GIPPA_SPEC.md §10.2 Question, then how to read the figure (R-04)
CAPTIONS: Mapping[str, str] = {
    "V01": "Where in the genome is the proteome? From the outside in, the rings show features, "
    "genes on each strand coloured by detection status, abundance, peptide density, "
    "transcription units with promoters and terminators, abundance islands, GC skew and the "
    "variant and novel-ORF sites.",
    "V02": "Is abundance clustered along the chromosome? Each point is a detected gene's log2 "
    "abundance at its position; the line is a LOESS trend; shaded bands are abundance islands.",
    "V03": "What happens at this prophage, BGC, island, flagged operon or novel ORF? Arrows are "
    "the genes of the window coloured by log2 abundance, with the mapped peptides stacked above "
    "them and promoters, terminators and regions marked along the track.",
    "V05": "Where did my sequences go? Bands follow the distinct input peptides or proteins from "
    "their input class through the matching route to the outcome, with widths proportional to "
    "the number of sequences.",
    "V06": "Can I trust tolerant hits and de novo peptides? Panel (a) counts target and decoy "
    "matches by peptide length and number of differences, panel (b) shows the estimated FDR per "
    "number of differences against the dashed acceptance level, and panel (c), when present, "
    "shows de novo precision per score bin.",
    "V07": "Which regions and categories are detected less than detectability predicts? Each "
    "row is a group's observed minus expected number of detected genes with its confidence "
    "interval; rows left of zero are detected less often than predicted.",
    "V08": "Which functions dominate abundance? Each column is a COG category, showing the log2 "
    "abundances of its measured genes as points or a violin with the median, or its detected "
    "fraction when the run has no abundance.",
    "V10": "Is there structure beyond operons? Panel (a) shows Moran's I and panel (b) the "
    "join-count statistic at each lag along the replicon, against shaded null envelopes; filled "
    "markers are significant under the TU-block null.",
    "V11": "Do operon partners look alike, and which other neighbours do too? Panel (a) compares "
    "the abundance differences of adjacent gene pairs by pair class, panel (b) plots the local "
    "FDR against that difference, and panel (c) plots upstream against downstream abundance with "
    "flagged pairs outlined.",
    "V12": "Where do peptides support ALT residues? Panel (a) counts variant sites by call and "
    "confidence, and panel (b) aligns the covering peptide of each top site to the reference "
    "with the ALT residue highlighted.",
    "V13": "Which new genes or upstream extensions are supported? Panel (a) counts the "
    "candidates of each class by call status, and panel (b) draws the top called candidates in "
    "their genomic neighbourhood.",
    "V14": "Which functions are high or low? Panel (a) shows the GSEA running sums of the top "
    "sets, panel (b) plots each set's normalised enrichment score against −log10 q, and panel "
    "(c) shows detection over-representation as observed over expected.",
    "V15": "What explains abundance? Panel (a) shows how much the cross-validated R² drops when "
    "each feature block is left out, panel (b) the partial dependence of key features, panel (c) "
    "the SHAP values of the top features, and panel (d) observed against predicted abundance.",
    "V16": "How do search strategies differ? Panel (a) counts the genes detected by each "
    "combination of runs, and panel (b) shows the abundance correlation and the shared log2 "
    "abundances of each pair of runs.",
    "V17": "Are the inputs sane? Panels show the preflight summary, contig lengths, peptide "
    "lengths, abundances and missing values per feature.",
    "V18": "Are the tests calibrated on this dataset? Panel (a) compares predicted with observed "
    "detection, panel (b) shows the p-value histogram of each test family, which a calibrated "
    "test keeps flat away from zero, and panel (c) places each observed statistic in its null "
    "distribution.",
}

#: plot ID -> the stem of its index instance and that instance's caption (review r1 L2)
INDEX_STEMS: Mapping[str, str] = {"V03": "V03_locus", "V13": "V13_novel_orfs_locus"}
INDEX_CAPTIONS: Mapping[str, str] = {
    "V03": "This index lists every selected locus in display order, with its kind and the error "
    "estimate of its test family, and marks the loci not drawn as locus views.",
    "V13": "This index lists the top called candidates drawn as locus views, by rank.",
}


def question(pid: str) -> str:
    """The GIPPA_SPEC.md §10.2 Question that opens ``CAPTIONS[pid]``."""
    return CAPTIONS[pid].split("? ", 1)[0] + "?"


def caption(pid: str, stem: str) -> str:
    """The caption of instance ``stem`` of plot ``pid``: the index instance has its own (the
    Question, then what the index lists), any other the plot's ``CAPTIONS`` entry."""
    if INDEX_STEMS.get(pid) == stem:
        return f"{question(pid)} {INDEX_CAPTIONS[pid]}"
    return CAPTIONS.get(pid, stem)


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
