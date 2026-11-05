# ruff: noqa: RUF001, RUF002, RUF003  (figure text is Unicode: −, ≤, Δ)
"""Shared fixtures and artist lookups for the M12-04 PR 5 tests (palette, labels and number
style in plots; docs/specs/M12-04-outputs-revision.md § PR 5).

- ``render(out, pids)``: the registered plots ``pids`` on an ``OUT/`` tree under the snapshot
  style (``helpers.snapshot_fixture.styled``), with the tree's own ``resolved_config.yaml`` when
  it has one; V16 (compare input, not a run) is never drawn here.
- ``sim_copy(tmp)``: a writable copy of the committed ``sim_plots`` snapshot ``OUT/`` (the
  richest gippasim run: peptides, iBAQ (molar), every table); ``rewrite_table`` edits one of its
  tables row by row (gzip or plain in, plain out; ``RunView`` reads both).
- ``k_like(tmp)``: fixture "K" for AC 5.2. The spec's K is a large K-12-like run; no such fixture
  exists, so this is ``sim_plots`` with the numbers K produced that broke the style rules:
  tiny q values (``4.9e-07`` and ``1.234e-05``), q values with 3+ significant digits (``0.432``,
  ``0.0197``) and a negative lag-1 Moran's I (``−0.147806``, A's value).
- ``figure_texts(inst)``: every visible text of an instance after a draw (titles, suptitles,
  axis labels, tick labels, legend texts, annotations, colour-bar labels) plus ``inst.title``.
- ``gene_patches(ax)``: the gene arrows of a locus view: closed polygons (patches or members of a
  patch collection) whose vertices all lie in ``|y| ≤ 0.3`` of the gene axes (region and island
  shading spans the axes height, the novel ORF band sits at y = 1.7).
"""

from __future__ import annotations

import csv
import gzip
import io
import re
import shutil
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib.colors as mcolors
import numpy as np
from matplotlib.axes import Axes
from matplotlib.collections import Collection
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Polygon
from matplotlib.text import Text

from helpers import snapshot_fixture as sf

SIM_OUT = sf.fixture_out("sim_plots")

#: a q value written after "q", "q =" or "q:" (a true or ASCII minus allowed, an exponent too)
Q_NUMBER = re.compile(r"\bq(?:\s*[=:]\s*|\s+)([−\-]?\d[\d.]*(?:[eE][−\-+]?\d+)?)")
#: an ASCII hyphen-minus written as a sign: before a digit and not joining two word characters
ASCII_MINUS = re.compile(r"(?<![\w.])-\d")


# ---------------------------------------------------------------- rendering


def render(out: Path, pids: Iterable[str] | None = None) -> dict[str, Any]:
    """``pid -> list[Instance] | Skip`` of the registered plots ``pids`` (all but V16 when
    ``None``) on ``out``, drawn under the snapshot style."""
    from gippa.config import resolve_config
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base

    wanted = None if pids is None else set(pids)

    def run() -> dict[str, Any]:
        load_plot_modules()
        cfg_path = out / "resolved_config.yaml"
        cfg = resolve_config(cfg_path if cfg_path.is_file() else None)
        view = RunView.from_outdir(out)
        found: dict[str, Any] = {}
        for pid, spec in base.registered().items():
            if pid == "V16" or (wanted is not None and pid not in wanted):
                continue
            found[pid] = spec.fn(view, cfg)
        return found

    result: dict[str, Any] = sf.styled(run)
    return result


def instances(results: dict[str, Any]) -> list[tuple[str, Any]]:
    """``(pid, Instance)`` of every drawn instance (skipped plots dropped)."""
    from gippa.viz import base

    return [(pid, inst) for pid, r in results.items() if not isinstance(r, base.Skip)
            for inst in r]  # fmt: skip


def only(results: dict[str, Any], pid: str, suffix: object = None) -> Any:
    """The one instance of ``pid`` with ``suffix``."""
    found = [i for p, i in instances(results) if p == pid and i.suffix == suffix]
    assert len(found) == 1, (pid, suffix, [(p, i.suffix) for p, i in instances(results)])
    return found[0]


# ---------------------------------------------------------------- fixture copies


def _read(path: Path) -> tuple[list[str], list[str], list[dict[str, str]]]:
    raw = path.read_bytes()
    text = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode("utf-8")
    lines = text.splitlines(keepends=True)
    comments = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if not ln.startswith("#")]
    reader = csv.DictReader(io.StringIO("".join(body)), delimiter="\t")
    rows = list(reader)
    return comments, list(reader.fieldnames or []), rows


def rewrite_table(
    out: Path, rel: str, fn: Callable[[dict[str, str]], dict[str, str] | None]
) -> None:
    """Apply ``fn`` to every row of ``out/rel`` (``None`` drops the row); written plain."""
    path = out / rel
    comments, columns, rows = _read(path)
    kept = [r2 for r in rows if (r2 := fn(dict(r))) is not None]
    buf = io.StringIO()
    buf.write("".join(comments))
    writer = csv.DictWriter(buf, columns, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(kept)
    path.write_text(buf.getvalue(), encoding="utf-8")


def read_table(out: Path, rel: str) -> list[dict[str, str]]:
    return _read(out / rel)[2]


def sim_copy(tmp: Path, name: str = "OUT") -> Path:
    """A writable copy of the ``sim_plots`` snapshot ``OUT/``."""
    dst = tmp / name
    shutil.copytree(SIM_OUT, dst)
    return dst


def set_abundance(out: Path, atype: str, scale: str) -> None:
    """Every gene of ``out``'s ``genes.tsv`` gets ``abundance_type``/``abundance_scale``."""

    def fn(r: dict[str, str]) -> dict[str, str]:
        r["abundance_type"], r["abundance_scale"] = atype, scale
        return r

    rewrite_table(out, "tables/genes.tsv", fn)


#: fixture K's numbers (see the module docstring)
K_TINY_Q = 4.9e-07
K_SMALL_Q = 1.234e-05
K_LONG_Q = 0.432
K_MORAN = -0.147806


def k_like(tmp: Path) -> Path:
    out = sim_copy(tmp, "K")

    def regions(r: dict[str, str]) -> dict[str, str]:
        if r["stat"] == "o_minus_e" and r["region_type"] == "prophage":
            r["p_value"] = r["q_value"] = repr(K_TINY_Q)
        elif r["stat"] == "o_minus_e":
            r["p_value"] = r["q_value"] = repr(K_LONG_Q)
        return r

    def orfs(r: dict[str, str]) -> dict[str, str]:
        r["best_peptide_q"] = repr(K_TINY_Q)
        return r

    def ac(r: dict[str, str]) -> dict[str, str]:
        if r["statistic"] == "morans_i" and r["lag"] == "1":
            r["value"] = repr(K_MORAN)
            r["q_value"] = repr(K_SMALL_Q if r["null_model"] == "tu_block" else K_LONG_Q)
            r["p_value"] = r["q_value"]
        return r

    def flags(r: dict[str, str]) -> dict[str, str]:
        r["q_value"] = repr(K_SMALL_Q)
        return r

    rewrite_table(out, "tables/regions_stats.tsv", regions)
    rewrite_table(out, "tables/novel_orfs.tsv", orfs)
    rewrite_table(out, "tables/autocorrelation.tsv", ac)
    rewrite_table(out, "tables/operon_flags.tsv", flags)
    return out


# ---------------------------------------------------------------- text


def figure_texts(inst: Any) -> list[str]:
    """``inst.title`` and every visible, non-empty text artist of its figure after a draw."""
    fig: Figure = inst.figure
    fig.canvas.draw()
    found = [inst.title] if inst.title else []
    found += [t.get_text() for t in fig.findobj(Text) if t.get_visible() and t.get_text()]
    return found


def sig_digits(number: str) -> int:
    """Significant digits of a written number (``"0.050"`` -> 2, ``"4.9e-07"`` -> 2)."""
    mant = re.split(r"[eE]", number.replace("−", "-").lstrip("-+"))[0]
    digits = mant.replace(".", "").lstrip("0")
    return max(len(digits), 1)


def style_violations(text: str) -> list[str]:
    """What breaks AC 5.2 in one text: a q with more than 2 significant digits, an ASCII ``-``
    before a digit, or ``e-0``."""
    bad = [f"q {m.group(1)!r} has {sig_digits(m.group(1))} significant digits"
           for m in Q_NUMBER.finditer(text) if sig_digits(m.group(1)) > 2]  # fmt: skip
    if ASCII_MINUS.search(text):
        bad.append("ASCII '-' before a digit")
    if "e-0" in text:
        bad.append("'e-0'")
    return bad


# ---------------------------------------------------------------- artists


def axes(fig: Figure, label: str) -> Axes:
    found = [ax for ax in fig.axes if ax.get_label() == label]
    assert len(found) == 1, [ax.get_label() for ax in fig.axes]
    return found[0]


def colorbar_labels(fig: Figure) -> list[str]:
    return [ax.get_xlabel() or ax.get_ylabel() for ax in fig.axes
            if ax.get_label() == "<colorbar>"]  # fmt: skip


def hexc(colour: Any) -> str:
    return str(mcolors.to_hex(colour, keep_alpha=False)).upper()


@dataclass(frozen=True)
class GenePatch:
    x0: float
    x1: float
    face: str
    edge: str
    hatch: str | None
    hatch_colour: str | None


def _hatch_colour(art: Any) -> str | None:
    getter = getattr(art, "get_hatchcolor", None)
    if getter is None:
        return None
    try:
        value = getter()
    except Exception:
        return None
    try:
        return hexc(value[0] if np.ndim(value) == 2 else value)
    except (ValueError, TypeError, IndexError):
        return None


def gene_patches(ax: Axes) -> list[GenePatch]:
    out: list[GenePatch] = []
    for p in ax.patches:
        if not isinstance(p, Polygon):
            continue
        xy = np.asarray(p.get_xy(), dtype=float)
        if len(xy) and np.all(np.abs(xy[:, 1]) <= 0.3):
            out.append(GenePatch(float(xy[:, 0].min()), float(xy[:, 0].max()),
                                 hexc(p.get_facecolor()), hexc(p.get_edgecolor()), p.get_hatch(),
                                 _hatch_colour(p)))  # fmt: skip
    for c in ax.collections:
        if not isinstance(c, Collection):
            continue
        faces, edges = c.get_facecolors(), c.get_edgecolors()
        for i, path in enumerate(c.get_paths()):
            v = np.asarray(path.vertices, dtype=float)
            if not len(v) or not np.all(np.abs(v[:, 1]) <= 0.3) or len(v) < 4:
                continue
            face = hexc(faces[i % len(faces)]) if len(faces) else "NONE"
            edge = hexc(edges[i % len(edges)]) if len(edges) else "NONE"
            out.append(GenePatch(float(v[:, 0].min()), float(v[:, 0].max()), face, edge,
                                 c.get_hatch(), _hatch_colour(c)))  # fmt: skip
    return out


def patches_of(ax: Axes, start: float, end: float) -> list[GenePatch]:
    """The gene arrows drawn inside ``[start, end]`` (window-frame bp)."""
    return [g for g in gene_patches(ax) if g.x0 >= start - 1e-6 and g.x1 <= end + 1e-6]


def legend_entries(fig: Figure) -> list[tuple[str, Any]]:
    """``(label, handle)`` of every figure and axes legend."""
    legends = list(fig.legends) + [ax.get_legend() for ax in fig.axes if ax.get_legend()]
    out: list[tuple[str, Any]] = []
    for lg in legends:
        handles = getattr(lg, "legend_handles", None) or getattr(lg, "legendHandles", [])
        out += [(t.get_text(), h) for t, h in zip(lg.get_texts(), handles, strict=False)]
    return out


def patch_style(handle: Any) -> tuple[str, str, str | None]:
    """``(face, edge, hatch)`` of a legend patch handle."""
    assert isinstance(handle, Patch), type(handle)
    return hexc(handle.get_facecolor()), hexc(handle.get_edgecolor()), handle.get_hatch()


def labels_of(fig: Figure, names: Sequence[str]) -> dict[str, tuple[str, str]]:
    """``axes label -> (xlabel, ylabel)`` for the axes named ``names``."""
    return {ax.get_label(): (ax.get_xlabel(), ax.get_ylabel()) for ax in fig.axes
            if ax.get_label() in names}  # fmt: skip
