# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −, ≤, …)
"""Shared builders and artist checks for the M12-04 PR 6 tests: V07 layout, V18 small-n, V17
labels (docs/specs/M12-04-outputs-revision.md § PR 6; P-08, P-10, P-16).

- ``drawn_texts(fig)``: every text artist matplotlib actually draws, after an Agg draw (the
  writer's renderer): visible, non-empty, and for tick labels only the ticks inside the axis
  view interval (matplotlib keeps label objects for ticks it does not draw).
- ``out_of_bounds(fig)``: AC 6.4, the drawn texts whose window extent leaves the figure bbox
  (1 px tolerance), the same check M9 uses for clipping.
- ``v07_out``/``v17_out``/``v18_out``: hand-built ``OUT/`` trees (R11: nothing row-level from
  the example bundle), in the formats the producers write (as the M9 V07/V17/V18 tests build
  them).
"""

from __future__ import annotations

import json
import math
import random
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
from matplotlib.axes import Axes
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.figure import Figure
from matplotlib.text import Text

from helpers import support_analysis_data as sa
from helpers.viz_out_tree import gene_row, make_out, write_events, write_schema_table

MINUS = "−"
ELLIPSIS = "…"

# ---------------------------------------------------------------- drawn text


def draw(fig: Figure) -> FigureCanvasAgg:
    canvas = FigureCanvasAgg(fig)
    canvas.draw()
    return canvas


def _undrawn_tick_labels(fig: Figure) -> set[int]:
    """``id`` of every tick label whose tick lies outside its axis view interval."""
    hidden: set[int] = set()
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            span = max(abs(hi - lo), 1e-12)
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tick.get_loc()
                if loc is None or not (lo - 1e-9 * span <= loc <= hi + 1e-9 * span):
                    hidden.update({id(tick.label1), id(tick.label2)})
    return hidden


def drawn_texts(fig: Figure) -> list[Text]:
    draw(fig)
    hidden = _undrawn_tick_labels(fig)
    return [t for t in fig.findobj(Text)
            if t.get_visible() and t.get_text().strip() and id(t) not in hidden]  # fmt: skip


def texts(fig: Figure) -> list[str]:
    return [t.get_text() for t in drawn_texts(fig)]


def ax_texts(fig: Figure, ax: Axes) -> list[Text]:
    """The drawn texts of one axes (tick labels, axis labels, titles, annotations)."""
    mine = {id(t) for t in ax.findobj(Text)}
    return [t for t in drawn_texts(fig) if id(t) in mine]


def out_of_bounds(fig: Figure) -> list[str]:
    box = fig.bbox
    bad = []
    for t in drawn_texts(fig):
        e = t.get_window_extent()
        if e.x0 < box.x0 - 1 or e.x1 > box.x1 + 1 or e.y0 < box.y0 - 1 or e.y1 > box.y1 + 1:
            bad.append(t.get_text())
    return bad


def is_bold(t: Text) -> bool:
    w = t.get_fontweight()
    return w in ("bold", "heavy", "black", "extra bold", "semibold", "demibold", "demi") or (
        isinstance(w, int | float) and w >= 600
    )


def middle_shortened(token: str, full: str, limit: int = 24) -> bool:
    """``token`` is ``full`` shortened with one middle ellipsis to at most ``limit`` chars."""
    if len(token) > limit or token.count(ELLIPSIS) != 1:
        return False
    head, tail = token.split(ELLIPSIS)
    return bool(head) and bool(tail) and full.startswith(head) and full.endswith(tail)


def isna(v: object) -> bool:
    return v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v))


def num(v: object) -> float:
    return math.nan if isna(v) else float(str(v))


# ---------------------------------------------------------------- V07


def basic(
    gtype: str,
    group: str,
    *,
    n_genes: int,
    n_detected: int,
    expected: float,
    ome: float | None = None,
    ci: tuple[float | None, float | None] = (None, None),
    flags: str | None = None,
) -> list[dict[str, object]]:
    """The seven ``basic`` rows of one group, as ``gippa.context.summary.basic_rows`` writes."""
    value = n_detected - expected if ome is None else ome

    def row(metric: str, v: object, lo: object = None, hi: object = None, f: object = None) -> dict:
        return {"section": "basic", "group_type": gtype, "group": group, "metric": metric,
                "value": v, "ci_low": lo, "ci_high": hi, "n": n_genes, "flags": f}  # fmt: skip

    return [
        row("n_genes", float(n_genes)),
        row("n_detected", float(n_detected)),
        row("pct_detected", 100.0 * n_detected / n_genes if n_genes else None),
        row("expected", expected),
        row("o_minus_e", value, ci[0], ci[1], flags),
        row("o_over_e", n_detected / expected if expected else None),
        row("n_tu", float(n_genes)),
    ]


def region(
    region_id: str,
    region_type: str,
    *,
    n_eligible: int = 25,
    observed: float | None = 10.0,
    expected: float | None = 10.0,
    effect: float | None = 0.0,
    ci: tuple[float | None, float | None] = (-1.0, 1.0),
    q: float | None = 0.5,
    flags: str | None = None,
    label: str | None = None,
) -> dict[str, object]:
    return {
        "region_id": region_id, "region_type": region_type, "subset": "all", "contig": "chr1",
        "start": 1001, "end": 9000, "n_genes": n_eligible, "n_eligible": n_eligible,
        "observed": observed, "expected": expected, "o_minus_e": effect, "o_over_e": None,
        "hl_shift_log2": None, "stat": "o_minus_e", "effect": effect, "ci_low": ci[0],
        "ci_high": ci[1], "p_value": q, "q_value": q,
        "method": "replicon_contrast" if region_type == "plasmid" else "circular_shift_o_minus_e",
        "null_model": "none" if region_type == "plasmid" else "circular_shift",
        "n_perm": 0 if region_type == "plasmid" else 999, "label": label, "flags": flags,
    }  # fmt: skip


#: an A-like V07 input: three basic blocks (one COG row) and two tested regions, none passing
A_BASIC = [
    *basic("cog_category", "J", n_genes=12, n_detected=3, expected=6.4, ci=(-5.0, -1.5)),
    *basic("feature_class", "cds", n_genes=40, n_detected=24, expected=19.0, ci=(-1.0, 9.0)),
    *basic("replicon", "chr1", n_genes=40, n_detected=20, expected=17.9, ome=2.1, ci=(-1.0, 5.2)),
]
A_REGIONS = [
    region("prophage_1", "prophage", n_eligible=25, observed=8.0, expected=10.0, effect=-2.0,
           ci=(-5.0, 1.0), q=0.2),
    region("bgc_1", "bgc", n_eligible=10, observed=6.0, expected=5.0, effect=1.0,
           ci=(-1.0, 3.0), q=0.4),
    region("plasmid_1", "plasmid", n_eligible=6, effect=0.8, ci=(-0.5, 2.1), q=None,
           label="dosage_plus_regulation", flags="dosage_plus_regulation"),
    region("crispr_2", "crispr", n_eligible=0, observed=None, expected=None, effect=None,
           ci=(None, None), q=None, flags="small_region;few_shifts"),
]  # fmt: skip


def v07_out(
    tmp: Path,
    *,
    summary: Sequence[Mapping[str, object]] | None = A_BASIC,
    regions: Sequence[Mapping[str, object]] | None = A_REGIONS,
    flags: Sequence[str] = (),
) -> Path:
    out = make_out(tmp)
    if summary is not None:
        write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", summary)
    if regions is not None:
        write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats", regions)
    if flags:
        write_events(out / "logs" / "events.jsonl",
                     [{"event": "layer_missing", "data": {"flag": f}} for f in flags])  # fmt: skip
    return out


def many_regions(n: int, seed: int = 0) -> list[dict[str, object]]:
    """``n`` prophage regions with distinct |O − E| / n (no ties at the display cut)."""
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        n_el = rng.randint(5, 60)
        effect = round((i + 1) * 0.01 * (1 if i % 2 else -1) * n_el / 10, 6)
        rows.append(region(f"prophage_{i + 1}", "prophage", n_eligible=n_el, observed=5.0,
                           expected=5.0 - effect, effect=effect,
                           ci=(effect - 1.0, effect + 1.0), q=0.5))  # fmt: skip
    return rows


# ---------------------------------------------------------------- V17


def check(rule_id: str, level: str, message: str = "") -> dict[str, Any]:
    return {"rule_id": rule_id, "input": "genome", "level": level,
            "rule": f"rule text of {rule_id}", "value": None, "message": message,
            "exit_code": 1 if level == "error" else 0}  # fmt: skip


TWO_REPLICONS = [
    {"id": "chr1", "kind": "chromosome", "length": 100_000, "topology": "circular"},
    {"id": "plasmid_1", "kind": "plasmid", "length": 8_000, "topology": "circular"},
]


def preflight(
    checks: Sequence[Mapping[str, Any]], replicons: Sequence[Mapping[str, Any]] = TWO_REPLICONS
) -> dict[str, Any]:
    levels = ("ok", "info", "warn", "error")
    return {
        "schema_version": "1.0", "gippa_version": "0.0.1", "inputs": [], "checks": list(checks),
        "summary": {f"n_{lvl}": sum(c["level"] == lvl for c in checks) for lvl in levels},
        "evidence_tier": "decoy_picked",
        "routing": {"n_proteins_long": 0, "n_proteins_peptide_path": 0, "n_proteins_dropped": 0,
                    "n_peptides": 0, "abundance_used": True, "abundance_scale": "log2"},
        "replicons": [dict(r) for r in replicons],
        "contig_map": [{"id": r["id"], "original_id": r["id"]} for r in replicons],
    }  # fmt: skip


#: every panel-(e) feature column with a value (the 100 %-missing ones are set per test)
FEATURE_VALUES: dict[str, object] = {
    "abundance_raw": 1000.0, "cai": 0.5, "tai": 0.3, "enc": 45.0, "fop": 0.4, "gc3": 0.5,
    "cpb": 0.01, "at_2_10": 0.6, "rare_codons_2_10": 1, "promoter_score": 2.0,
    "promoter_score_tu": 2.5, "tx_rate_max_tu": 3.0, "tx_rate_sum_tu": 4.0,
    "terminator_conf_downstream": 80.0, "mw_da": 30000.0, "pi": 6.5, "gravy": -0.3,
    "tm_segments": 0, "n_observable_peptides": 7, "tir_ostir": 1000.0, "sd_dg": -5.0,
    "mrna_dg": -10.0, "spacer_bp": 7, "start_codon": "ATG",
}  # fmt: skip


def v17_gene(tag: str, *, detected: bool, measured: bool, imputed: bool = False,
             missing: Sequence[str] = (), **kw: object) -> dict[str, object]:  # fmt: skip
    log2 = 20.0 if measured or imputed else None
    values = {k: (None if k in missing else v) for k, v in FEATURE_VALUES.items()}
    if not measured:
        values["abundance_raw"] = None
    row = gene_row(locus_tag=tag, feature_type="CDS", detected=detected,
                   detection_status="detected" if detected else "not_detected",
                   abundance_log2=log2, abundance_imputed=imputed,
                   abundance_type="ibaq" if log2 is not None else "none", **values)  # fmt: skip
    row.update(kw)
    return row


def v17_out(
    tmp: Path,
    *,
    checks: Sequence[Mapping[str, Any]] = (check("genome.parse", "ok"),
                                           check("genome.contiguity", "warn", "small")),
    replicons: Sequence[Mapping[str, Any]] = TWO_REPLICONS,
    genes: Sequence[Mapping[str, object]] | None = None,
    flags: Sequence[str] = (),
) -> Path:  # fmt: skip
    out = make_out(tmp)
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(preflight(checks, replicons), sort_keys=True), encoding="utf-8")
    if genes is not None:
        write_schema_table(out / "tables" / "genes.tsv", "genes", genes)
    if flags:
        write_events(out / "logs" / "events.jsonl",
                     [{"event": "flag", "data": {"flag": f}} for f in flags])  # fmt: skip
    return out


# ---------------------------------------------------------------- V18


def ds(group_type: str, group: object, metric: str, value: float | None) -> dict[str, object]:
    return {"section": "detectability", "group_type": group_type, "group": str(group),
            "metric": metric, "value": value}  # fmt: skip


def overall(constant: bool) -> list[dict[str, object]]:
    rows = [
        ds("overall", "all", "brier", 0.11),
        ds("overall", "all", "calibration_intercept", None if constant else 0.02),
        ds("overall", "all", "calibration_slope", None if constant else 0.97),
        ds("overall", "all", "n_detected", 40.0),
        ds("overall", "all", "n_eligible", 100.0),
    ]
    if constant:
        rows.append(ds("overall", "all", "detectability_constant", 1.0))
    return rows  # fmt: skip


def reliability(constant: bool) -> list[dict[str, object]]:
    """A fitted detectability model (two bins) or a constant one (one bin, no curve)."""
    pts = [(0.4, 0.4, 100)] if constant else [(0.1, 0.12, 50), (0.6, 0.55, 50)]
    rows = []
    for i, (mp, rate, n) in enumerate(pts, start=1):
        rows += [ds("bin", i, "mean_p", mp), ds("bin", i, "observed_rate", rate),
                 ds("bin", i, "n", n)]  # fmt: skip
    return rows + overall(constant)


def spread(n: int, seed: int = 0) -> list[float]:
    """``n`` p-values in (0, 1), away from the bin edges."""
    rng = random.Random(seed)
    return [round(0.0125 + 0.05 * rng.randrange(20), 4) for _ in range(n)]


def region_p(p: float, i: int, stat: str = "o_minus_e") -> dict[str, object]:
    return {"region_id": f"R{i:04d}", "region_type": "prophage", "subset": "all",
            "contig": "chr1", "start": 101 + 1000 * i, "end": 900 + 1000 * i, "n_genes": 5,
            "stat": stat, "p_value": p, "method": "x", "null_model": "circular_shift"}  # fmt: skip


def autocorr_p(p: float, statistic: str, null_model: str, lag: int) -> dict[str, object]:
    return {"contig": "chr1", "statistic": statistic, "lag": lag, "null_model": null_model,
            "value": 0.1, "p_value": p}  # fmt: skip


def pair_p(p: float, test: str, i: int) -> dict[str, object]:
    return {"pair_id": f"P{i:04d}", "class": "OP", "gene_a": f"a{i}", "gene_b": f"b{i}",
            "strand": "+", "intergenic_bp": 20, "test": test, "p_value": p}  # fmt: skip


def v18_out(
    tmp: Path,
    *,
    constant: bool = True,
    n_regions: int = 0,
    n_hl: int = 0,
    n_join_tu: int = 0,
    n_moran_gl: int = 0,
    n_pairs: int = 0,
    seed: int = 0,
) -> Path:
    """Detectability (fitted or constant) and p-value families of the given sizes: regions
    ``o_minus_e`` and ``hl_shift``, autocorrelation ``join_count_bb / tu_block`` and
    ``morans_i / gene_label``, pairs ``both_detected_t``."""
    out = make_out(tmp)
    write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary",
                       reliability(constant))  # fmt: skip
    regions = [region_p(p, i) for i, p in enumerate(spread(n_regions, seed))]
    regions += [region_p(p, 500 + i, "hl_shift") for i, p in enumerate(spread(n_hl, seed + 1))]
    if regions:
        write_schema_table(out / "tables" / "regions_stats.tsv", "regions_stats", regions)
    ac = [
        autocorr_p(p, "join_count_bb", "tu_block", lag)
        for lag, p in enumerate(spread(n_join_tu, seed + 2), start=1)
    ]
    ac += [autocorr_p(p, "morans_i", "gene_label", lag)
           for lag, p in enumerate(spread(n_moran_gl, seed + 3), start=1)]  # fmt: skip
    if ac:
        write_schema_table(out / "tables" / "autocorrelation.tsv", "autocorrelation", ac)
    if n_pairs:
        pairs = [pair_p(p, "both_detected_t", i) for i, p in enumerate(spread(n_pairs, seed + 4))]
        write_schema_table(out / "tables" / "pairs.tsv", "pairs", pairs)
    return out


#: a long contig name like the one the spec quotes (40 characters)
LONG_CONTIG = "nissle1917_NZ_CP058217_1_1965000_2165000"


def island_row(contig: str, p: float) -> dict[str, object]:
    return {"island_id": "ISL1", "direction": "high", "contig": contig, "start": 1001,
            "end": 5400, "tu_first": "TU0001", "tu_last": "TU0003", "n_tus": 3, "n_genes": 5,
            "genes": "g1;g2;g3;g4;g5", "mean_abundance_log2": 21.5, "detected_frac": 0.8,
            "scan_s": 4.25, "block_len": 3, "stat": "scan_max", "effect": 1.5, "p_value": p,
            "method": "scan_max_stat", "null_model": "circular_block_bootstrap",
            "n_perm": 999}  # fmt: skip


def add_island_nulls(out: Path, contig: str = LONG_CONTIG, n: int = 50) -> None:
    """``null_islands`` for one scanned contig (and its ``islands`` row)."""
    from gippa.io.tsv import write_tsv

    rng = random.Random(7)
    write_schema_table(out / "tables" / "islands.tsv", "islands", [island_row(contig, 0.03)])
    path = out / "plots" / "support" / "null_islands.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = [[contig, "abs", b, repr(round(rng.gauss(0.0, 1.5), 4)), repr(4.2)]
            for b in range(1, n + 1)]  # fmt: skip
    write_tsv(path, list(sa.COLUMNS["null_islands"]), rows)
