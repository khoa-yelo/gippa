"""The head of ``report.html`` (M12-04 PR 1): the "This run" box, the version line and the
table of contents.

Pure apart from the reads a ``RunView`` does: a view, the plot files, the skip reasons and the
config in, template data out. Numbers go through the ``gippa.viz.style`` formatters.
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from typing import Any

import pandas as pd

import gippa
from gippa.config import GippaConfig
from gippa.io.plot_files import PlotFiles, svg_title
from gippa.io.run_view import RunView
from gippa.viz import style
from gippa.viz.style import fmt_alpha, fmt_count
from gippa.viz.words import TIER_NAME, TIER_TEXT

# ---------------------------------------------------------------- inputs

NOT_RECORDED = "not recorded"
UNREADABLE = "could not be read"
GENOME_SUFFIXES = (".gz", ".fna", ".fasta", ".fa")
FINDING_PLOTS = ("V01", "V07", "V02")
#: section id -> the table whose row count the TOC shows
TOC_COUNTS: Mapping[str, str] = {
    "detection": "genes_tsv",
    "variants": "variants",
    "novel-orfs": "novel_orfs",
}


#: what reading a present but malformed table raises (``ValueError`` covers decoding, parquet
#: and schema errors; ``KeyError`` a missing column; ``EOFError`` a truncated ``.gz``)
READ_ERRORS = (OSError, ValueError, KeyError, EOFError)


def _safe_table(view: RunView, key: str) -> pd.DataFrame | None:
    try:
        return view.table(key) if view.has(key) else None
    except READ_ERRORS:  # the section preview carries the note
        return None


def _unreadable(view: RunView, key: str) -> bool:
    """Whether ``key`` is present but cannot be read."""
    if not view.has(key):
        return False
    try:
        view.table(key)
    except READ_ERRORS:
        return True
    return False


def _safe_doc(read: Any) -> Mapping[str, Any] | None:
    try:
        doc = read()
    except Exception:
        return None
    return doc if isinstance(doc, Mapping) else None


def recorded_inputs(view: RunView) -> list[Mapping[str, Any]] | None:
    """The input records of the run: the manifest's ``inputs``, else (no manifest, or a 1.0
    manifest without them) the preflight's; ``None`` when neither records any."""
    for read in (view.manifest, view.preflight):
        rows = (_safe_doc(read) or {}).get("inputs")
        found = [r for r in rows if isinstance(r, Mapping)] if isinstance(rows, list) else []
        if found:
            return found
    return None


def manifest_roles(view: RunView) -> frozenset[str] | None:
    """The input roles of the run manifest's own ``inputs``; ``None`` without them. The cards
    gate on these only: a run without a manifest record falls back to its tables."""
    rows = (_safe_doc(view.manifest) or {}).get("inputs")
    if not isinstance(rows, list) or not rows:
        return None
    return frozenset(str(r.get("role", "")) for r in rows if isinstance(r, Mapping))


def input_roles(view: RunView) -> frozenset[str] | None:
    """The input roles of ``recorded_inputs``, or ``None`` when no inputs are recorded (an old
    run: the cards fall back to the tables, the box says "not recorded")."""
    rows = recorded_inputs(view)
    return None if rows is None else frozenset(str(r.get("role", "")) for r in rows)


def genome_stem(name: str) -> str:
    """``nissle.fna.gz`` -> ``nissle``: the file name without its directory, ``.gz`` and FASTA
    suffix."""
    stem = name.replace("\\", "/").rsplit("/", 1)[-1]
    for suffix in GENOME_SUFFIXES:
        stem = stem.removesuffix(suffix)
    return stem


def n_unmapped(view: RunView, kind: str) -> int:
    unmapped = _safe_table(view, "unmapped")
    return 0 if unmapped is None else int((unmapped["kind"] == kind).sum())


def n_input(view: RunView, kind: str) -> int | None:
    """Distinct input sequences of ``kind`` (``peptide`` or ``protein``): placed plus unmapped;
    ``None`` when the hits table is absent."""
    key, col = ("peptide_hits", "peptide") if kind == "peptide" else ("protein_hits", "protein_id")
    hits = _safe_table(view, key)
    if hits is None:
        return None
    return int(hits[col].nunique()) + n_unmapped(view, kind)


def _input_line(view: RunView, roles: frozenset[str] | None) -> str:
    if roles is None:
        return f"input: {NOT_RECORDED}"
    parts = []
    if "peptides" in roles:
        hits = _safe_table(view, "peptide_hits")
        denovo = (hits is not None and len(hits) > 0 and "source" in hits
                  and set(hits["source"].dropna().astype(str)) == {"denovo"})  # fmt: skip
        n = n_input(view, "peptide")
        name = "de novo peptides" if denovo else "peptides"
        parts.append(name if n is None else f"{name} ({fmt_count(n)})")
    if "proteins" in roles:
        n = n_input(view, "protein")
        parts.append("proteins" if n is None else f"proteins ({fmt_count(n)})")
    return f"input: {', '.join(parts) or NOT_RECORDED}"


def _abundance_line(genes: pd.DataFrame | None) -> str:
    label = None
    if genes is not None and "abundance_type" in genes:
        types = genes["abundance_type"].dropna()
        if not types.empty:
            scales = genes.get("abundance_scale")
            scale = None
            if scales is not None:
                found = scales[genes["abundance_type"] == types.iloc[0]].dropna()
                scale = None if found.empty else found.iloc[0]
            label = style.abundance_label(types.iloc[0], scale)
    return f"abundance: {label or 'none'}"


def _replicon_line(view: RunView) -> str:
    reps = (_safe_doc(view.preflight) or {}).get("replicons")
    if not isinstance(reps, list) or not reps:
        return f"replicons: {NOT_RECORDED}"
    raw = [r.get("length") for r in reps if isinstance(r, Mapping)]
    lengths = [int(x) for x in raw if isinstance(x, int | float)]
    n = len(reps)
    word = "replicon" if n == 1 else "replicons"
    if not lengths or len(lengths) != len(reps):
        return f"{n} {word}"
    return f"{n} {word}, {fmt_count(sum(lengths))} bp"


def tier_of(view: RunView, genes: pd.DataFrame | None) -> str | None:
    """The run's detection tier: the manifest's ``detection_tier``, else the single tier of the
    gene table; ``None`` when neither gives one."""
    doc = _safe_doc(view.manifest)
    tier = doc.get("detection_tier") if doc is not None else None
    if isinstance(tier, str) and tier:
        return tier
    if genes is not None and "detection_tier" in genes:
        tiers = sorted({str(x) for x in genes["detection_tier"].dropna()})
        if len(tiers) == 1:
            return tiers[0]
    return None


def tier_sentence(tier: str | None, alpha: float) -> str:
    """``TIER_TEXT`` of ``tier`` with the run's alpha, after the tier's plain name."""
    if tier is None or tier not in TIER_TEXT:
        return f"detection tier: {NOT_RECORDED}"
    return f"{TIER_NAME.get(tier, tier)}: {TIER_TEXT[tier].format(alpha=fmt_alpha(alpha))}"


def version_line(doc: Mapping[str, Any] | None) -> str:
    line = f"Report built by GIPPA {gippa.__version__}"
    made = doc.get("gippa_version") if doc is not None else None
    if isinstance(made, str) and made and made != gippa.__version__:
        line += f" from a run made with GIPPA {made}"
    return line


def findings(plots: Mapping[str, PlotFiles], skipped: Mapping[str, str]) -> list[str]:
    """The finding titles of V01, V07 and V02 (in that order), from the SVG ``<title>`` of each
    plot's first instance; a skipped plot, an empty-state plot or one without a title gives
    none."""
    out = []
    for pid in FINDING_PLOTS:
        files = plots.get(pid)
        if pid in skipped or files is None:
            continue
        inst = next((i for i in files.instances if i.svg is not None), None)
        if inst is None or inst.empty:
            continue
        title = None if inst.svg is None else svg_title(inst.svg)
        if title and title.strip():
            out.append(" ".join(title.split()))
    return out


def run_box(
    view: RunView, plots: Mapping[str, PlotFiles], skipped: Mapping[str, str], cfg: GippaConfig
) -> dict[str, Any]:
    """The "This run" box: its lines, the key findings and the version line (behaviour 1, 2).
    A run without genes has no key findings: its plots say only that there is nothing to show."""
    rows = recorded_inputs(view)
    genome = next((r for r in rows or [] if r.get("role") == "genome"), None)
    name = None if genome is None else genome.get("file") or genome.get("path")
    genes = _safe_table(view, "genes_tsv")
    if genes is not None:
        gene_line = f"{fmt_count(len(genes))} genes"
    else:
        gene_line = f"genes: {UNREADABLE if _unreadable(view, 'genes_tsv') else NOT_RECORDED}"
    lines = [
        f"genome: {genome_stem(str(name)) if name else NOT_RECORDED}",
        _replicon_line(view),
        gene_line,
        _input_line(view, input_roles(view)),
        _abundance_line(genes),
        f"modules: {', '.join(str(m) for m in cfg.modules)}",
    ]
    return {
        "lines": lines,
        "tier": tier_sentence(tier_of(view, genes), cfg.stats.alpha.identification),
        "findings": [] if genes is not None and genes.empty else findings(plots, skipped),
        "version": version_line(_safe_doc(view.manifest)),
    }


def toc(view: RunView, sections: Sequence[tuple[str, str]]) -> list[dict[str, str]]:
    """One TOC entry per ``(id, heading)``: ``href`` ``#sec-<id>`` and the heading, with the row
    count of the section's main table in brackets when that table exists (behaviour 4)."""
    out = []
    for sid, heading in sections:
        text = heading
        key = TOC_COUNTS.get(sid)
        frame = None if key is None else _safe_table(view, key)
        if frame is not None:
            text = f"{heading} ({fmt_count(len(frame))})"
        out.append({"href": f"#sec-{sid}", "text": text})
    return out


def mean_sd(mean: float, sd: float) -> str:
    """``"m ± s"``: the mean through ``fmt_num``, the SD to the mean's decimals, more when needed
    so it keeps one significant digit (``0.21 ± 0.004``); a true minus."""
    m = style.fmt_num(mean)
    if "e" in m or not math.isfinite(sd):
        return f"{m} ± {style.fmt_num(sd)}"
    decimals = len(m.split(".", 1)[1]) if "." in m else 0
    s = abs(sd)
    if s > 0:
        decimals = max(decimals, -math.floor(math.log10(s)))
    return f"{m} ± {s:.{decimals}f}"


def abundance_counts(genes: pd.DataFrame | None) -> tuple[int, int] | None:
    """Genes with a finite ``abundance_log2`` and their distinct ``tu_id``; ``None`` without a
    gene table."""
    if genes is None or "abundance_log2" not in genes:
        return None
    vals = pd.to_numeric(genes["abundance_log2"], errors="coerce").astype(float)
    ok = genes[vals.map(math.isfinite)]
    tus = ok["tu_id"].dropna().nunique() if "tu_id" in ok else 0
    return len(ok), int(tus)
