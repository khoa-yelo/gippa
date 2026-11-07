"""Pure helpers of the M12-01 slow route test (tests/e2e/test_m12_01_example_routes.py): the
`expected/` keys and tolerances of spec section 2, measuring a finished route, comparing it, and
splitting the printed command block. ``expected/<route>.json`` is ``{key: {"value": v,
"tolerance": t}}`` with ``t`` one of ``"exact"``, ``"set"``, ``"informational"`` or
``{"abs": a, "rel": r}`` (``|got - expected| <= max(a, r * |expected|)``)."""

from __future__ import annotations

import csv
import json
import math
import shlex
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

ROUTES = ("A", "B", "C", "D", "E")
CONVERTED = frozenset("ACDE")  # first command `gippa convert`
ABUNDANCE = frozenset("ABDE")
NOVEL = frozenset("E")  # the six-frame route: the novel-ORF keys (route E, 2026-11-06)
NOVEL_KEYS = ("n_novel_candidates", "n_novel_called", "n_novel_unpaired", "novel_classes")
NOVEL_CLASS_FIELDS = (
    "regime",
    "class",
    "estimated",
    "n_targets",
    "n_target_wins",
    "n_decoys",
    "q_floor",
    "floor_above_alpha",
)
_COUNT_2 = {"abs": 1, "rel": 0.02}
TOLERANCES: dict[str, Any] = {  # spec section 2, `expected/` (key numbers)
    "status": "exact",
    "detection_tier": "exact",
    "layers_missing": "set",
    "convert_rows": "exact",
    "n_genes": {"abs": 1, "rel": 0.01},
    "n_detected": _COUNT_2,
    "n_detected_2pep": _COUNT_2,
    "n_peptides_placed": _COUNT_2,
    "n_abundance": _COUNT_2,
    "n_novel_candidates": _COUNT_2,
    "n_novel_called": "exact",
    "n_novel_unpaired": _COUNT_2,
    "novel_classes": "informational",
    "warnings": "informational",
    "wall_time_s": "informational",
}


def required_keys(route: str) -> set[str]:
    keys = set(TOLERANCES)
    if route not in CONVERTED:
        keys.discard("convert_rows")
    if route not in ABUNDANCE:
        keys.discard("n_abundance")
    if route not in NOVEL:
        keys -= set(NOVEL_KEYS)
    return keys


def within(got: float, exp: float, tol: Mapping[str, float]) -> bool:
    return abs(got - exp) <= max(float(tol.get("abs", 0)), float(tol.get("rel", 0)) * abs(exp))


def _miss(key: str, got: Any, exp: Any, tol: Any) -> str | None:
    if tol == "exact":
        return None if got == exp else f"{key}: got {got!r}, expected exactly {exp!r}"
    if tol == "set":
        ok = isinstance(got, list) and isinstance(exp, list) and set(got) == set(exp)
        return None if ok else f"{key}: got {got!r}, expected the set {exp!r}"
    if not (isinstance(tol, Mapping) and tol and set(tol) <= {"abs", "rel"}):
        return f"{key}: unknown tolerance {tol!r}"
    try:
        g, e = float(got), float(exp)
    except (TypeError, ValueError):
        return f"{key}: not numeric (got {got!r}, expected {exp!r})"
    if not (math.isfinite(g) and math.isfinite(e)) or not within(g, e, tol):
        return f"{key}: got {got!r}, expected {exp!r} within {dict(tol)}"
    return None


def compare(measured: Mapping[str, Any], doc: Mapping[str, Any], route: str) -> list[str]:
    """Misses of ``measured`` against ``expected/<route>.json`` (empty: all within tolerance)."""
    misses = [f"{k}: absent from expected/" for k in sorted(required_keys(route) - set(doc))]
    for key in sorted(doc):
        entry = doc[key]
        if not (isinstance(entry, Mapping) and {"value", "tolerance"} <= set(entry)):
            misses.append(f"{key}: not {{'value': ..., 'tolerance': ...}} ({entry!r})")
        elif entry["tolerance"] == "informational":
            continue
        elif key not in measured:
            misses.append(f"{key}: not measured (expected {entry['value']!r})")
        elif (m := _miss(key, measured[key], entry["value"], entry["tolerance"])) is not None:
            misses.append(m)
    return misses


def as_expected(measured: Mapping[str, Any], route: str) -> dict[str, dict[str, Any]]:
    """``measured/<route>.json``: the nested layout with the spec tolerances, ready for
    ``build_example.py --expected``."""
    return {k: {"value": measured.get(k), "tolerance": TOLERANCES[k]}
            for k in sorted(required_keys(route))}  # fmt: skip


def _tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))


def _true(v: str | None) -> bool:
    return (v or "").strip().lower() == "true"


def measure(route: str, outdir: Path, peptides: Path | None = None) -> dict[str, Any]:
    """The section 2 keys of one route, read from its run directory and converted table."""
    got: dict[str, Any] = {}
    if (outdir / "run_manifest.json").is_file():
        man = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
        got["status"] = man.get("status")
        got["detection_tier"] = man.get("detection_tier")
        got["layers_missing"] = sorted(man.get("layers_missing") or [])
        got["warnings"] = sorted(man.get("warnings") or [])
    if route in CONVERTED and peptides is not None and peptides.is_file():
        got["convert_rows"] = len(_tsv(peptides))
    if (outdir / "tables" / "genes.tsv").is_file():
        genes = _tsv(outdir / "tables" / "genes.tsv")
        got["n_genes"] = len(genes)
        got["n_detected"] = sum(_true(g.get("detected")) for g in genes)
        got["n_detected_2pep"] = sum(_true(g.get("detected_2pep")) for g in genes)
        if route in ABUNDANCE:
            got["n_abundance"] = sum((g.get("abundance_raw") or "NA") not in ("", "NA")
                                     and not _true(g.get("abundance_imputed"))
                                     for g in genes)  # fmt: skip
    if (hits := outdir / "mapping" / "peptide_hits.tsv").is_file():
        got["n_peptides_placed"] = len({h["peptide"] for h in _tsv(hits)})
    if route in NOVEL:
        got |= _novel(outdir)
    return got


def _novel(outdir: Path) -> dict[str, Any]:
    """Route E's novel-ORF keys (spec section 2, `expected/`). ``n_novel_unpaired`` counts the
    ``no_makedb_pair`` rows of ``intermediate/novel/novel_class_fdr.parquet``: 0 when it has no
    data rows or is absent from a finished run."""
    got: dict[str, Any] = {}
    if (orfs := outdir / "tables" / "novel_orfs.tsv").is_file():
        rows = _tsv(orfs)
        got["n_novel_candidates"] = len(rows)
        got["n_novel_called"] = sum(_true(r.get("called")) for r in rows)
    fdr = next(  # the runner writes it under intermediate/ (runner/local.py)
        (p for p in (outdir / "intermediate" / "novel" / "novel_class_fdr.parquet",
                     outdir / "novel" / "novel_class_fdr.parquet") if p.is_file()),
        None,
    )  # fmt: skip
    if fdr is not None:
        import pandas as pd

        table = pd.read_parquet(fdr)
        reason = table["reason"] if "reason" in table.columns else pd.Series([], dtype=object)
        got["n_novel_unpaired"] = int((reason == "no_makedb_pair").sum())  # 0 with no data rows
    elif orfs.is_file() or (outdir / "run_manifest.json").is_file():
        # not written because the run has no candidate: 0, never null (amendment 2026-11-06 (b))
        got["n_novel_unpaired"] = 0
    if (man := outdir / "run_manifest.json").is_file():
        block = json.loads(man.read_text(encoding="utf-8")).get("novel_orfs") or {}
        got["novel_classes"] = [{k: c.get(k) for k in NOVEL_CLASS_FIELDS}
                                for c in block.get("classes") or []]  # fmt: skip
    return got


def printed_routes(lines: Iterable[str], routes: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    """The printed block -> ``{route id: [command lines]}``, headed by the YAML titles."""
    by_title = {f"# {r['title']}": str(r["id"]) for r in routes}
    out: dict[str, list[str]] = {}
    current: str | None = None
    for line in lines:
        if line in by_title:
            current = by_title[line]
            out[current] = []
        elif current is not None and line.startswith("gippa "):
            out[current].append(line)
        else:
            current = None
    return out


def argv(line: str, gippa_bin: str, run_args: list[str]) -> list[str]:
    """A printed command as argv; ``run_args`` go on ``gippa run`` only."""
    words = shlex.split(line)
    assert words[0] == "gippa", line
    return [gippa_bin, *words[1:], *(run_args if words[1] == "run" else [])]


def option(words: list[str], name: str) -> Path | None:
    """The value of ``name`` in an argv (``--out``, ``--outdir``), else ``None``."""
    return Path(words[words.index(name) + 1]) if name in words[:-1] else None
