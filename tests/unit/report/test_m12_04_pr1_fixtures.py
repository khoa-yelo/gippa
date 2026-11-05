# ruff: noqa: RUF002, RUF003  (the report's wording is Unicode: ≤, ≥, ², −)
"""Fixtures P, R, D (and their variants) for the M12-04 PR 1 tests: the run box, the
input-aware summary cards and the TOC (docs/specs/M12-04-outputs-revision.md § PR 1).

The fixtures are hand-built, A-, B- and C-like ``OUT/`` trees (R11: nothing row-level from the
example bundle). The counts mirror the evidence runs the spec quotes, so the spec's examples
read directly as expectations:

- ``P``: peptides (DIA-NN-like, ``source = database``), tier ``mayu``, iBAQ (molar), 158 genes,
  99 detected (91 with ≥ 2 distinct peptides), 774 distinct peptides on 783 hit rows, Mayu
  gene-level FDR 0.0124093 at q = 0.01, empty ``protein_hits.tsv``, and a determinants model
  that was not fitted (``too_few_rows``; 99 genes with abundance in 34 TUs);
- ``R``: proteins only, tier ``protein_list``, MaxLFQ (relative), 158 genes, 106 detected from
  106 proteins, an empty ``peptide_hits.tsv``, and an elastic-net R² of −0.0174037 ± 0.0656274
  with null p 0.414585 (``null_uncalibrated``);
- ``D``: de novo peptides (``source = denovo``), tier ``two_peptide``, no abundance
  (``abundance_type = none``, run flag ``layer_missing:abundance``), 158 genes, 4 detected,
  21 peptides all placed, 16 descriptive placement strata in ``mapping_fdr.tsv`` (the ≈ 300
  character string the spec moves off the card), and univariate-only determinants rows;
- variants: ``P_nokegg`` (P with ``layer_missing:kegg``), ``P_core`` (P with ``--modules core``),
  ``decoy`` (P re-tiered ``decoy_picked``), ``both`` (peptides and proteins).

Interface assumed beyond the spec text (PR 1 test author):
- the run box reads the genome file name and the input roles from ``run_manifest.json ›
  inputs`` (``helpers.report_content.manifest`` builds a schema-valid one); the preflight JSON
  carries the same inputs and the replicons; ``preflight.json › replicons`` gives the replicon
  count and total length;
- modules come from the step's ``cfg.modules``; every fixture writes the same list to
  ``resolved_config.yaml`` and the manifest ``config.modules``, and the report is built from
  ``RunView.from_outdir(out).with_modules(cfg.modules)`` as ``gippa.steps.plots`` does;
- the abundance line uses ``abundance_label(genes.abundance_type, genes.abundance_scale)``;
- the key findings are the finding titles the SVGs carry as ``<title>`` (``write_instance``;
  ``gippa.io.plot_files.svg_title``), passed to ``build_report`` as ``PlotFiles``.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest
from helpers.report_content import det_rows, manifest, words, write_det_model, write_flags
from helpers.report_fixtures import Page, cfg, parse
from helpers.viz_out_tree import (
    gene_row,
    make_out,
    write_resolved_config,
    write_schema_table,
)

GENOME_FILE = "nissle_r200k.fna.gz"
GENOME_STEM = "nissle_r200k"
REPLICON_BP = 200_000
N_GENES = 158
CONTIG = "chr1"

#: the module lists of the evidence runs (A: core,bgc,function; B, C: core,bgc)
MODULES_A = ("core", "bgc", "function")
MODULES_BC = ("core", "bgc")

#: finding titles of the fake plots (the spec's run-A V01 title and two made-up ones)
TITLE_V01 = "99 of 158 genes detected (63 %); 0 islands at FWER ≤ 0.05; 0 ALT sites, 0 novel ORFs"
TITLE_V07 = "Observed vs expected detection: 99 detected, 81 expected"
TITLE_V02 = "Abundance along the chromosome: 0 islands"
TITLES = {"V01": TITLE_V01, "V07": TITLE_V07, "V02": TITLE_V02}
STEMS = {"V01": "V01_circos", "V07": "V07_obs_vs_exp", "V02": "V02_abundance_landscape"}

FIXTURES = ("P", "R", "D", "P_nokegg", "P_core", "decoy", "both")


@dataclass
class Run:
    """One built fixture: its ``OUT/``, the build config and the plot inputs of the report."""

    name: str
    out: Path
    tier: str
    modules: tuple[str, ...]
    overrides: tuple[str, ...] = ()
    plots: dict[str, Any] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)

    def cfg(self) -> Any:
        return cfg(f"modules=[{', '.join(self.modules)}]", *self.overrides)

    def view(self) -> Any:
        from gippa.io.run_view import RunView

        return RunView.from_outdir(self.out).with_modules(list(self.modules))

    def html(self) -> str:
        from gippa.report.build import build_report

        return build_report(self.view(), self.plots, self.skipped, self.cfg())

    def page(self) -> Page:
        return parse(self.html())

    def cards(self) -> list[dict[str, Any]]:
        from gippa.report.build import cards

        found, _ = cards(self.view(), self.cfg())
        return found


# ---------------------------------------------------------------- writers


def write_preflight(out: Path, inputs: Sequence[Mapping[str, Any]], n_replicons: int = 1) -> None:
    path = out / "preflight" / "preflight.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    size = REPLICON_BP // n_replicons
    body = {
        "schema_version": "1",
        "gippa_version": "0.0.0",
        "inputs": [dict(i) for i in inputs],
        "replicons": [{"id": f"chr{k + 1}", "length": size, "topology": "circular"}
                      for k in range(n_replicons)],
    }  # fmt: skip
    path.write_text(json.dumps(body, sort_keys=True, indent=1), encoding="utf-8")


def inputs_for(roles: Sequence[str]) -> list[dict[str, Any]]:
    files = {"genome": GENOME_FILE, "peptides": "peptides.tsv", "proteins": "proteins.faa",
             "abundance": "abundance.tsv"}  # fmt: skip
    return [{"role": r, "file": files[r], "path": files[r], "sha256": "c" * 64, "size": 100}
            for r in roles]  # fmt: skip


def write_run_manifest(
    out: Path, *, roles: Sequence[str] | None, tier: str, modules: Sequence[str],
    version: str | None = None,
) -> None:  # fmt: skip
    """``run_manifest.json``: schema-valid, with ``inputs`` of ``roles`` (dropped when ``None``,
    an old run), ``detection_tier`` and ``config.modules``."""
    import gippa
    from gippa.manifest import canonical_json

    doc = manifest()
    doc["gippa_version"] = version or gippa.__version__
    doc["detection_tier"] = tier
    if isinstance(doc.get("config"), dict):
        doc["config"]["modules"] = list(modules)
    if roles is None:
        doc.pop("inputs", None)
    else:
        doc["inputs"] = inputs_for(roles)
    (out / "run_manifest.json").write_text(canonical_json(doc), encoding="utf-8")


def write_genes(
    out: Path, *, n_detected: int, n_2pep: int, tier: str, abundance: tuple[str, str] | None,
    n_tus: int = 34,
) -> None:  # fmt: skip
    """``N_GENES`` genes on ``CONTIG``; the first ``n_detected`` detected (the first ``n_2pep``
    with ``detected_2pep``), with a log2 abundance when ``abundance`` is ``(type, scale)``; TU
    ``i % n_tus`` so the detected genes span ``min(n_detected, n_tus)`` TUs."""
    a_type, a_scale = abundance if abundance is not None else ("none", None)
    rows = []
    for i in range(N_GENES):
        det = i < n_detected
        start = 101 + 1000 * i
        rows.append(gene_row(
            locus_tag=f"g{i:05d}", contig=CONTIG, start=start, end=start + 299,
            feature_type="cds", tu_id=f"TU{i % n_tus + 1:04d}", detected=det,
            detected_2pep=i < n_2pep, detection_status="detected" if det else "not_detected",
            detection_tier=tier, n_peptides_unique=(3 if i < n_2pep else 1) if det else 0,
            abundance_log2=20.0 + 0.01 * i if det and abundance is not None else None,
            abundance_type=a_type, abundance_scale=a_scale,
        ))  # fmt: skip
    write_schema_table(out / "tables" / "genes.tsv", "genes", rows)


def write_tier_fdr(out: Path, fdr: float) -> None:
    rows = [{"section": "tier", "group_type": "q_threshold", "group": "0.01", "metric": "fdr",
             "value": fdr},
            {"section": "tier", "group_type": "q_threshold", "group": "0.05", "metric": "fdr",
             "value": 0.0411}]  # fmt: skip
    write_schema_table(out / "tables" / "detection_summary.tsv", "detection_summary", rows)


def write_peptides(out: Path, n: int, *, source: str, dup: int = 0) -> None:
    """``n`` distinct peptides, the first ``dup`` of them on two rows (multi-mapping), all at
    ``placement_q = 0``; ``n = 0`` writes the header only (as the evidence runs do)."""
    rows = []
    for i in range(n):
        for copy in range(2 if i < dup else 1):
            start = 101 + 30 * i + 10_000 * copy
            rows.append({"peptide": f"PEPTIDE{i:04d}K", "peptide_len": 12, "source": source,
                         "contig": CONTIG, "start": start, "end": start + 36, "strand": "+",
                         "route": "exact", "n_diff": 0, "placement_q": 0.0,
                         "length_class": "12-14", "uniqueness": "unique"})  # fmt: skip
    write_schema_table(out / "mapping" / "peptide_hits.tsv", "peptide_hits", rows)


def write_proteins(out: Path, n: int) -> None:
    rows = [{"protein_id": f"JNOIEE_{i:05d}", "protein_len_aa": 100, "route": "miniprot",
             "contig": CONTIG, "start": 101 + 1000 * i, "end": 400 + 1000 * i, "strand": "+",
             "rank": 1, "class": "identical"} for i in range(n)]  # fmt: skip
    write_schema_table(out / "mapping" / "protein_hits.tsv", "protein_hits", rows)


def write_unmapped(out: Path, n_peptides: int = 0, n_proteins: int = 0) -> None:
    rows = [{"sequence_id": f"u{i}", "sequence": "PEPTIDEK", "kind": "peptide",
             "label": "unmapped", "reason": "no_hit"} for i in range(n_peptides)]  # fmt: skip
    rows += [{"sequence_id": f"up{i}", "sequence": "MPEPTIDE", "kind": "protein",
              "label": "unmapped", "reason": "below_identity"}
             for i in range(n_proteins)]  # fmt: skip
    write_schema_table(out / "mapping" / "unmapped.tsv", "unmapped", rows)


def write_strata(out: Path, n_diff: Sequence[int] = (0,), lengths: Sequence[str] = ("7",)) -> None:
    """Descriptive placement strata of ``mapping_fdr.tsv``, both decoy types."""
    rows = [{"section": "placement", "decoy_type": d, "n_diff": n, "length_class": lc,
             "targets": 10, "decoys": 0, "fdr_stratum": 0.0, "n_accepted": 10}
            for d in ("reverse", "shuffle") for n in n_diff for lc in lengths]  # fmt: skip
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", rows)


def write_variants(out: Path) -> None:
    """A header-only ``variants.tsv`` (the TOC count "Variants (0)")."""
    write_schema_table(out / "tables" / "variants.tsv", "variants", [])


def write_not_fitted(out: Path, n_genes: int) -> None:
    """A's determinants rows: both whole models ``too_few_rows`` (value and n NA), univariate
    rows with ``n = n_genes``."""
    rows: list[dict[str, Any]] = [
        {"model": m, "feature_or_block": "all", "metric": "cv_r2_mean", "value": None, "n": None,
         "n_perm": 0, "flags": "too_few_rows;ci_unestimated"}
        for m in ("lightgbm", "elastic_net")
    ]  # fmt: skip
    rows += [{"model": "univariate", "feature_or_block": f, "metric": "spearman_rho",
              "value": 0.1, "ci_low": -0.1, "ci_high": 0.3, "n": n_genes, "n_perm": 0,
              "flags": "ci_unestimated"} for f in ("cpb", "gc3")]  # fmt: skip
    write_det_model(out, rows)


def write_univariate_auroc(out: Path) -> None:
    """C's determinants rows: univariate AUROC only (no abundance, so no regression model)."""
    rows = [{"model": "univariate", "feature_or_block": f, "metric": "auroc", "value": 0.6,
             "ci_low": 0.4, "ci_high": 0.8, "n": N_GENES, "n_perm": 0,
             "flags": "ci_unestimated"} for f in ("cpb", "gc3")]  # fmt: skip
    write_det_model(out, rows)


def write_fitted(
    out: Path, r2: float, sd: float, p: float, flags: str = "null_uncalibrated"
) -> None:
    write_det_model(out, det_rows("elastic_net", r2, sd, p, flags))


# ---------------------------------------------------------------- fixtures


def _base(tmp: Path, name: str, *, roles: Sequence[str] | None, tier: str,
          modules: Sequence[str], version: str | None = None,
          preflight_inputs: bool = True) -> Path:  # fmt: skip
    out = make_out(tmp, name)
    write_preflight(out, inputs_for(roles or ()) if preflight_inputs else [])
    write_run_manifest(out, roles=roles, tier=tier, modules=modules, version=version)
    write_resolved_config(out / "resolved_config.yaml", modules)
    return out


def fixture_p(tmp: Path, name: str = "P", *, modules: Sequence[str] = MODULES_A,
              tier: str = "mayu", **kw: Any) -> Run:  # fmt: skip
    out = _base(tmp, name, roles=kw.pop("roles", ("genome", "peptides")), tier=tier,
                modules=modules, **kw)  # fmt: skip
    write_genes(out, n_detected=99, n_2pep=91, tier=tier, abundance=("ibaq", "molar"))
    write_tier_fdr(out, 0.0124093)
    write_peptides(out, 774, source="database", dup=9)
    write_proteins(out, 0)
    write_unmapped(out)
    write_strata(out)
    write_not_fitted(out, 99)
    write_variants(out)
    return Run(name, out, tier, tuple(modules))


def fixture_r(tmp: Path, name: str = "R", **kw: Any) -> Run:
    out = _base(tmp, name, roles=("genome", "proteins", "abundance"), tier="protein_list",
                modules=MODULES_BC, **kw)  # fmt: skip
    write_genes(out, n_detected=106, n_2pep=0, tier="protein_list",
                abundance=("maxlfq", "relative"))  # fmt: skip
    write_peptides(out, 0, source="database")
    write_proteins(out, 106)
    write_unmapped(out)
    write_schema_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", [])
    write_fitted(out, -0.0174037, 0.0656274, 0.414585)
    return Run(name, out, "protein_list", MODULES_BC)


def fixture_d(tmp: Path, name: str = "D", **kw: Any) -> Run:
    out = _base(tmp, name, roles=("genome", "peptides"), tier="two_peptide", modules=MODULES_BC,
                **kw)  # fmt: skip
    write_genes(out, n_detected=4, n_2pep=4, tier="two_peptide", abundance=None)
    write_peptides(out, 21, source="denovo")
    write_proteins(out, 0)
    write_unmapped(out)
    write_strata(out, n_diff=(0,), lengths=("7", "8", "9", "10-11", "12-14", "15-18", "19-24",
                                            "25+"))  # fmt: skip
    write_univariate_auroc(out)
    write_flags(out, "abundance", ["layer_missing:abundance", "c_assumed"])
    return Run(name, out, "two_peptide", MODULES_BC)


def fixture_both(tmp: Path) -> Run:
    """Peptides and proteins both given (edge case: both cards shown)."""
    run = fixture_p(tmp, "both", roles=("genome", "peptides", "proteins"))
    write_proteins(run.out, 12)
    write_unmapped(run.out, n_proteins=3)
    return run


def build(name: str, tmp: Path) -> Run:
    if name == "P":
        return fixture_p(tmp)
    if name == "R":
        return fixture_r(tmp)
    if name == "D":
        return fixture_d(tmp)
    if name == "P_nokegg":
        run = fixture_p(tmp, name)
        write_flags(run.out, "function", ["layer_missing:kegg"])
        return run
    if name == "P_core":
        return fixture_p(tmp, name, modules=("core",))
    if name == "decoy":
        run = fixture_p(tmp, name, tier="decoy_picked")
        write_tier_fdr(run.out, 0.0098)
        return run
    if name == "both":
        return fixture_both(tmp)
    raise KeyError(name)


def fake_plots(pids: Sequence[str]) -> dict[str, Any]:
    """``PlotFiles`` whose SVG carries the finding title of ``TITLES`` as its ``<title>``."""
    from gippa.io.plot_files import InstanceFiles, PlotFiles

    out = {}
    for pid in pids:
        svg = ('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10">'
               f"<title>{TITLES[pid]}</title><g id=\"figure_1\"></g></svg>")  # fmt: skip
        stem = STEMS[pid]
        inst = InstanceFiles(stem=stem, svg=svg, images=(f"{stem}.svg",), data=None)
        out[pid] = PlotFiles(pid, (inst,))
    return out


def text(node: Any) -> str:
    return words(node.text())


# ---------------------------------------------------------------- the fixtures are what they say


@pytest.mark.parametrize("name", FIXTURES)
def test_fixture_tables_read_back(name: str, tmp_path: Path) -> None:
    """Every fixture is a readable run: genes, the tier, the manifest and the modules."""
    from gippa.io.run_view import RunView

    run = build(name, tmp_path)
    view = RunView.from_outdir(run.out)
    genes = view.table("genes_tsv")
    assert genes is not None and len(genes) == N_GENES
    assert set(genes["detection_tier"].dropna()) == {run.tier}
    doc = view.manifest()
    assert doc is not None and doc["detection_tier"] == run.tier
    assert view.preflight() is not None
    assert run.cfg().modules == list(run.modules)


def test_fixture_counts_mirror_the_evidence_runs(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    p = RunView.from_outdir(fixture_p(tmp_path).out)
    hits = p.table("peptide_hits")
    assert hits is not None and hits["peptide"].nunique() == 774 and len(hits) == 783
    genes = p.table("genes_tsv")
    assert genes is not None and int(genes["detected"].sum()) == 99
    ab = genes[genes["abundance_log2"].notna()]
    assert len(ab) == 99 and ab["tu_id"].nunique() == 34
    r = RunView.from_outdir(fixture_r(tmp_path).out)
    prot, pep = r.table("protein_hits"), r.table("peptide_hits")
    assert prot is not None and prot["protein_id"].nunique() == 106
    assert pep is not None and len(pep) == 0
    d = RunView.from_outdir(fixture_d(tmp_path).out)
    dh = d.table("peptide_hits")
    assert dh is not None and len(dh) == 21 and set(dh["source"]) == {"denovo"}
    assert "layer_missing:abundance" in d.run_flags()
