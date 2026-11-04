"""Issue #536: a proteins-only run (``--proteins``, optional ``--abundance``, no ``--peptides``)
produces the gene table through the protein path (Khoa 2026-11-03, option 1).

Contract pinned here (test author), from docs/user-guide/inputs.md (proteins are required
"unless peptides are given"; records of 40 aa or fewer go to the peptide path, longer ones to
miniprot, then tblastn, ``mapping/protein_hits.tsv``) and docs/user-guide/interpreting-results.md
#detection-tiers (``protein_list``: "proteins only, no peptides"; detected means "a mapped input
protein"; gene error ``unestimated``):

1. The step that writes the ``mapping`` directory (``loci``) runs without a ``peptides`` key: with
   only ``proteins_long`` and the miniprot / tblastn outputs (recorded fixtures, injected; no tool
   runs) it writes the seven §7.2.9 files, ``protein_hits.tsv`` holding every protein and
   ``peptide_hits.tsv`` header-only.
2. ``gene_table`` on a loci directory with protein hits only and no C9 peptides: every gene's
   ``detection_tier`` is ``protein_list``; a non-partial gene is ``detected`` exactly when a rank-1
   protein hit names it; ``gene_q_value`` is NA; no peptide counts.
3. ``--abundance`` is joined by the input protein ID through the rank-1 protein hits (C12), not only
   by locus tag.
4. ``gippa run --proteins`` (no ``--peptides``) passes ``proteins_long`` and ``abundance`` to
   Nextflow and never ``peptides`` (PR #537 review, lead decision, fix option a): the short records
   (≤ 40 aa, within ``input.peptide_len``) are not used, the tier stays ``protein_list``, and
   preflight warns ``proteins.short_records_unused`` with their count. With ``--peptides`` the short
   records join the peptide table as before.
5. PR #537 review M2: ``loci`` given ``peptides`` but not every peptide-chain output raises, naming
   the missing key; on a proteins-only run, proteins no tool aligns are ``kind = protein`` rows of
   ``unmapped.tsv``.

The Nextflow side (main.nf runs MAPPING / GENETABLE / CONTEXT on proteins only) is pinned in
workflow/tests/main_proteins_only.nf.test. Small examples only; no statistical claim is made, so no
calibration test (§8.4): ``protein_list`` detection is a deterministic rule.
"""

from __future__ import annotations

import json
import math
import os
import random
import re
import shutil
import stat
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE.parent / "helpers"))
sys.path.insert(0, str(_HERE / "genetable"))
import m4_contract_fixture as mf  # noqa: E402

N_GENES = 20


def _tag(i: int) -> str:
    return mf.tag(i)


# ---------------------------------------------------------------- hand-made contract (M4-06)


def _protein_hit(pid: str, gene_index: int, rank: int, fr: dict[str, Any]) -> dict[str, Any]:
    g = fr["genes"].iloc[gene_index]
    return {
        "protein_id": pid, "protein_len_aa": 99, "route": "miniprot", "contig": str(g["contig"]),
        "start": int(g["start"]), "end": int(g["end"]), "strand": str(g["strand"]), "rank": rank,
        "identity": 1.0, "query_coverage": 1.0, "reciprocal_overlap": 1.0, "class": "identical",
        "gene": _tag(gene_index), "flags": [],
    }  # fmt: skip


def proteins_only_frames(
    hits: list[tuple[str, int, int]],
    abundance: list[tuple[str, float]] | None = None,
    partial: tuple[int, ...] = (),
) -> dict[str, Any]:
    """The M4-06 hand fixture with no C9 peptides, a header-only C11, C12 = ``hits``
    (``(protein_id, gene index, rank)``) and C10 = ``abundance`` (``(protein_id, value)``)."""
    fr = mf.build_frames(N_GENES)
    for i in partial:
        fr["genes"].loc[i, "partial"] = True
        fr["layers"].at[i, "flags"] = ["partial"]
    fr["peptides"] = None
    fr["peptide_hits"] = mf.empty_like(fr["peptide_hits"])
    fr["unmapped"] = mf.empty_like(fr["unmapped"])
    rows = [_protein_hit(p, g, r, fr) for p, g, r in hits]
    fr["protein_hits"] = mf._typed_frame(rows, mf.PROTEIN_HITS_COLUMNS, mf.PROTEIN_HITS_TYPES)
    if abundance is None:
        fr["abundance"] = None
    else:
        fr["abundance"] = pd.DataFrame(
            {
                "protein_id": [p for p, _ in abundance],
                "value": pd.Series([v for _, v in abundance], dtype="float64"),
                "raw_id": [p for p, _ in abundance],
            }
        )
    return fr


def write_proteins_only_work(root: Path, fr: dict[str, Any]) -> dict[str, Path]:
    """Contract files in the ``gippa run`` layout (M3 TSVs under ``mapping/``, C14 under
    ``annotation/``); no ``peptides`` file, ``abundance`` only when given."""
    from gippa.io.work_tables import default_paths

    flat = default_paths(root) | {"preflight": root / "preflight.json"}
    skip = ["peptides"] + ([] if fr["abundance"] is not None else ["abundance"])
    mf.write_work(fr, flat, skip=skip)
    for name, sub in (("peptide_hits.tsv", "mapping"), ("protein_hits.tsv", "mapping"),
                      ("unmapped.tsv", "mapping"), ("id_map.tsv", "annotation")):  # fmt: skip
        (root / sub).mkdir(exist_ok=True)
        shutil.move(root / name, root / sub / name)
    keys = {k: p for k, p in flat.items() if k not in skip}
    keys = {k: p for k, p in keys.items() if k not in ("peptide_hits", "protein_hits", "unmapped")}
    return keys | {"mapping": root / "mapping", "id_map": root / "annotation" / "id_map.tsv"}


def run_gene_table(root: Path, fr: dict[str, Any]) -> pd.DataFrame:
    from gippa.config import GippaConfig
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    spec = reg.get("gene_table")
    keys = write_proteins_only_work(root / "work", fr)
    allowed = set(spec.inputs) | set(spec.optional_inputs)
    assert "peptides" not in keys
    inputs = {k: p for k, p in keys.items() if k in allowed}
    out = root / "out"
    outputs = {k: out / k for k in spec.outputs}
    outputs["genes_pre"] = out / "genetable" / "genes_pre.parquet"
    reg.run_step("gene_table", GippaConfig(seed=1), inputs, outputs)
    return pd.read_parquet(outputs["genes_pre"]).set_index("locus_tag")


def _na(x: Any) -> bool:
    return x is None or x is pd.NA or (isinstance(x, float) and math.isnan(x))


# ---------------------------------------------------------------- 2. protein_list tier


HITS = [("PROT_A", 2, 1), ("PROT_B", 5, 1), ("PROT_C", 9, 1), ("PROT_D", 12, 2)]


def test_gene_table_from_protein_hits_only(tmp_path: Path) -> None:
    """Protein hits only: tier ``protein_list``; detected = rank-1 protein hit on a non-partial
    gene; a rank-2 hit alone does not detect; a partial gene stays ``partial_excluded``."""
    df = run_gene_table(tmp_path, proteins_only_frames(HITS, partial=(9,)))
    assert len(df) == N_GENES
    assert set(df["detection_tier"]) == {"protein_list"}
    detected = {t for t, d in df["detected"].items() if bool(d)}
    assert detected == {_tag(2), _tag(5)}
    assert df.loc[_tag(9), "detection_status"] == "partial_excluded"
    assert df.loc[_tag(12), "detection_status"] == "not_detected"  # rank 2 only
    others = [t for t in df.index if t not in {_tag(2), _tag(5), _tag(9)}]
    assert set(df.loc[others, "detection_status"]) == {"not_detected"}
    assert df["gene_q_value"].isna().all()  # unestimated in this tier
    assert (df["n_peptides_unique"].fillna(0).astype(int) == 0).all()
    assert not df["detected_2pep"].fillna(False).astype(bool).any()


def test_gene_table_proteins_only_without_any_hit(tmp_path: Path) -> None:
    """No protein maps: the table still has every gene, nothing detected, tier ``protein_list``."""
    df = run_gene_table(tmp_path, proteins_only_frames([]))
    assert len(df) == N_GENES
    assert set(df["detection_tier"]) == {"protein_list"}
    assert not df["detected"].astype(bool).any()


@settings(max_examples=40, deadline=None, suppress_health_check=[HealthCheck.too_slow])
@given(
    rank1=st.sets(st.integers(0, N_GENES - 1), max_size=N_GENES),
    rank2=st.sets(st.integers(0, N_GENES - 1), max_size=6),
    partial=st.sets(st.integers(0, N_GENES - 1), max_size=4),
)
def test_property_protein_list_detects_exactly_rank1_hits(
    rank1: set[int], rank2: set[int], partial: set[int]
) -> None:
    """Invariant of the ``protein_list`` tier (pure ``detect``): detected genes are the non-partial
    genes named by a rank-1 C12 hit, whatever rank-2 hits exist; partial genes are excluded."""
    from gippa.config import GippaConfig
    from gippa.genetable.covariates import covariate_table
    from gippa.genetable.detection import detect
    from gippa.genetable.evidence import accepted_peptides, gene_evidence

    hits = [(f"P{g}", g, 1) for g in sorted(rank1)]
    hits += [(f"Q{g}", g, 2) for g in sorted(rank2)]
    fr = proteins_only_frames(hits, partial=tuple(sorted(partial)))
    inp = mf.make_inputs(fr)
    cfg = GippaConfig(seed=1)
    cov = covariate_table(inp.genes, cfg)
    acc = accepted_peptides(inp.peptide_hits, inp.peptides, cfg)
    ev = gene_evidence(inp.genes, inp.peptide_hits, acc, inp.replicons)
    det = detect(inp.genes, ev, acc, inp.peptides, inp.peptide_hits, inp.protein_hits, cov, cfg)
    assert det.tier == "protein_list"
    t = det.table.set_index("locus_tag")
    want = {_tag(g) for g in rank1 - partial}
    assert {x for x, d in t["detected"].items() if bool(d)} == want
    assert {x for x, s in t["detection_status"].items() if s == "partial_excluded"} == {
        _tag(g) for g in partial
    }


# ---------------------------------------------------------------- 3. abundance by protein ID


def test_abundance_is_joined_by_input_protein_id(tmp_path: Path) -> None:
    """C10 IDs are the input protein IDs (not locus tags): each value reaches the gene of that
    protein's rank-1 hit; a gene without a mapped protein value has no observed abundance."""
    ab = [("PROT_A", 21.5), ("PROT_B", 18.25), ("PROT_X", 30.0)]  # PROT_X maps nowhere
    df = run_gene_table(tmp_path, proteins_only_frames(HITS, abundance=ab))
    assert set(df["detection_tier"]) == {"protein_list"}
    # fixture preflight.json: routing.abundance_scale = log2, so log2 is the input value
    assert df.loc[_tag(2), "abundance_log2"] == pytest.approx(21.5)
    assert df.loc[_tag(5), "abundance_log2"] == pytest.approx(18.25)
    assert not bool(df.loc[_tag(2), "abundance_imputed"])
    assert _na(df.loc[_tag(0), "abundance_raw"])  # no protein of gene 0 in C10
    assert set(df["abundance_type"].dropna()) <= {"ibaq", "top3", "maxlfq", "intensity",
                                                  "unknown"}  # fmt: skip
    assert "none" not in set(df.loc[[_tag(2), _tag(5)], "abundance_type"])


# ---------------------------------------------------------------- 1. loci without peptides


@pytest.fixture(scope="module")
def sim100k(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """SIM100K (the genome of the recorded miniprot / tblastn fixtures) after ``gippa preflight``
    with its peptides, and the resolved configuration."""
    import dataclasses
    import hashlib

    from helpers.mapping_steps import preflight

    from gippa.config import resolve_config
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SIM100K, SIM100K_SEED, simulate

    tmp = tmp_path_factory.mktemp("sim100k")
    res = simulate(dataclasses.replace(SIM100K, proteomics=ProteomicsConfig()), SIM100K_SEED,
                   tmp / "sim")  # fmt: skip
    sha = hashlib.sha256(res.genome_fasta.read_bytes()).hexdigest()
    assert sha.startswith("dd10216e"), "SIM100K genome differs from the recorded fixtures"
    run = preflight(res, tmp / "pf")
    return {"res": res, "run": run, "cfg": resolve_config(run.config_yaml, ["seed=3"])}


def _protein_only_drive(
    sim: dict[str, Any],
    tmp: Path,
    extra: Sequence[tuple[str, str]] = (),
) -> dict[str, Any]:
    """The ``mapping`` producer with roots genome + bakta_json + proteins_long only (recorded
    miniprot / tblastn outputs injected; half the proteins go to tblastn). ``extra`` records are
    appended to the protein FASTA and align nowhere (neither tool output names them)."""
    from helpers.mapping_steps import (
        RECORDED_MINIPROT,
        RECORDED_TBLASTN,
        drive,
        fresh_steps,
        proteins_key,
    )

    tmp.mkdir(parents=True, exist_ok=True)
    src = RECORDED_MINIPROT / "proteins.faa"
    ids = [ln[1:].split()[0] for ln in src.read_text().splitlines() if ln.startswith(">")]
    faa = tmp / "proteins.long.faa"
    text = src.read_text().rstrip("\n") + "\n"
    faa.write_text(text + "".join(f">{i}\n{seq}\n" for i, seq in extra))
    drop = set(ids[::2])
    gff = _gff_without(RECORDED_MINIPROT / "miniprot.gff", drop, tmp / "miniprot.part.gff")
    tbl = _tsv_only(RECORDED_TBLASTN / "tblastn.tsv", drop, tmp / "tblastn.part.tsv")
    run, res = sim["run"], sim["res"]
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        producers = [s.name for s in reg.all_steps() if "mapping" in s.outputs]
        assert producers == ["loci"], producers
        roots = {"genome": run.genome, "bakta_json": res.bakta_json, proteins_key(reg): faa}
        (mp_key,) = reg.get("miniprot").outputs
        (tb_key,) = reg.get("tblastn").outputs
        avail = drive(reg, sim["cfg"], tmp / "work", roots,
                      inject={"miniprot": {mp_key: gff}, "tblastn": {tb_key: tbl}},
                      targets=tuple(producers))  # fmt: skip
    return {"ids": ids, "drop": drop, "work": tmp / "work", "avail": avail}


@pytest.fixture(scope="module")
def protein_only_loci(
    sim100k: dict[str, Any], tmp_path_factory: pytest.TempPathFactory
) -> dict[str, Any]:
    return _protein_only_drive(sim100k, tmp_path_factory.mktemp("proteins_only"))


#: proteins no tool output names (random 60-aa sequences)
UNALIGNED = [(f"NOHIT_{k}", "M" + "".join(random.Random(k).choice("ACDEFGHIKLMNPQRSTVWY")
                                         for _ in range(59))) for k in (1, 2)]  # fmt: skip


@pytest.fixture(scope="module")
def protein_only_loci_unaligned(
    sim100k: dict[str, Any], tmp_path_factory: pytest.TempPathFactory
) -> dict[str, Any]:
    return _protein_only_drive(sim100k, tmp_path_factory.mktemp("unaligned"), extra=UNALIGNED)


@pytest.fixture(scope="module")
def peptide_chain(sim100k: dict[str, Any], tmp_path_factory: pytest.TempPathFactory) -> Any:
    """The peptide chain up to ``placement`` on SIM100K's peptides (real work files for loci)."""
    from helpers.mapping_steps import drive, fresh_steps, roots_of

    tmp = tmp_path_factory.mktemp("peptide_chain")
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, sim100k["cfg"], tmp / "work", roots_of(sim100k["run"]),
                      targets=("placement",))  # fmt: skip
    return {"avail": avail, "tmp": tmp}


def _gff_without(src: Path, drop: set[str], dst: Path) -> Path:
    out: list[str] = []
    keep = True
    for ln in src.read_text().splitlines():
        if ln.startswith("##PAF"):
            keep = ln.split("\t")[1] not in drop
        elif ln.startswith("##"):
            keep = True
        if keep:
            out.append(ln)
    dst.write_text("\n".join(out) + "\n")
    return dst


def _tsv_only(src: Path, ids: set[str], dst: Path) -> Path:
    rows = [ln for ln in src.read_text().splitlines() if ln.split("\t")[0] in ids]
    dst.write_text("\n".join(rows) + "\n")
    return dst


def test_loci_runs_without_peptides(protein_only_loci: dict[str, Any]) -> None:
    from helpers.mapping_steps import MODULE2_FILES, module2_files, read_tsv_df

    avail = protein_only_loci["avail"]
    files = module2_files(protein_only_loci["work"], avail)
    assert set(files) == set(MODULE2_FILES)
    assert read_tsv_df(files["peptide_hits.tsv"]).empty  # header-only: there are no peptides
    hits = read_tsv_df(files["protein_hits.tsv"])
    assert set(hits["protein_id"]) == set(protein_only_loci["ids"])
    routes = hits.groupby("protein_id")["route"].agg(set)
    drop = protein_only_loci["drop"]
    for pid, r in routes.items():
        assert r == ({"tblastn"} if pid in drop else {"miniprot"}), (pid, r)
    first = hits[hits["rank"].astype(int) == 1].set_index("protein_id")
    assert (first["gene"] == first.index).all()  # recorded IDs are the locus tags
    un = read_tsv_df(files["unmapped.tsv"])
    assert un[un["kind"] == "protein"].empty


def test_proteins_only_loci_output_loads_as_the_m3_contract(
    protein_only_loci: dict[str, Any],
) -> None:
    """The loci directory is valid C11-C13 input for gene_table (no contract error)."""
    from helpers.mapping_steps import module2_files

    from gippa.io.work_tables import load_m4_inputs

    avail = protein_only_loci["avail"]
    files = module2_files(protein_only_loci["work"], avail)
    d = files["protein_hits.tsv"].parent
    _m4, issues = load_m4_inputs({"peptide_hits": d / "peptide_hits.tsv",
                                  "protein_hits": d / "protein_hits.tsv",
                                  "unmapped": d / "unmapped.tsv",
                                  "genes": avail["genes"],
                                  "replicons": avail["replicons"]})  # fmt: skip
    bad = [i for i in issues if i.severity == "error" and i.code.split(".")[1] in
           ("peptide_hits", "protein_hits", "unmapped")]  # fmt: skip
    assert not bad, bad


PEPTIDE_WORK = ("placement", "mapping_fdr", "hits_exact", "hits_tolerant")


@pytest.mark.parametrize("missing", PEPTIDE_WORK)
def test_loci_with_peptides_but_a_missing_chain_output_raises(
    peptide_chain: dict[str, Any], tmp_path: Path, missing: str
) -> None:
    """PR #537 review M2: ``peptides`` given means the peptide chain ran, so loci without one of
    its outputs fails (naming the key) instead of silently writing a proteins-only table."""
    from helpers.mapping_steps import fresh_steps

    from gippa.config import GippaConfig

    avail = peptide_chain["avail"]
    keys = ["genome", "genes", "features", "replicons", "peptides", "sixframe", *PEPTIDE_WORK]
    inputs = {k: avail[k] for k in keys if k != missing}
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        with pytest.raises(Exception) as err:  # StepError through run_step
            reg.run_step("loci", GippaConfig(seed=3), inputs, {"mapping": tmp_path / "mapping"})
    text = f"{err.value} {err.value.__cause__}"
    assert missing in text, f"loci failed without naming {missing}: {text}"


def test_proteins_only_loci_lists_unaligned_proteins_as_unmapped(
    protein_only_loci_unaligned: dict[str, Any],
) -> None:
    """PR #537 review M2: proteins neither miniprot nor tblastn aligns are ``kind = protein``
    rows of ``unmapped.tsv``; aligned proteins are not."""
    from helpers.mapping_steps import module2_files, read_tsv_df

    run = protein_only_loci_unaligned
    files = module2_files(run["work"], run["avail"])
    un = read_tsv_df(files["unmapped.tsv"])
    prot = un[un["kind"] == "protein"]
    assert set(prot["sequence_id"]) == {i for i, _ in UNALIGNED}, prot
    seqs = dict(UNALIGNED)
    for r in prot.itertuples(index=False):
        assert r.sequence == seqs[r.sequence_id]
    hits = read_tsv_df(files["protein_hits.tsv"])
    assert set(hits["protein_id"]) == set(run["ids"])


# ---------------------------------------------------------------- 4. gippa run params.json


LONG_PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
SHORT_PROTEIN = "MKVLAAGIVGRSTEQLLK"  # 18 aa: the peptide path (input.peptide_len 7..40)


@pytest.fixture
def run_env(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    from helpers.fake_env import fake_gippa_root, write_setup_manifest

    root = fake_gippa_root(tmp_path / "gippa_root")
    monkeypatch.setenv("GIPPA_ROOT", str(root))
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text("#!/bin/sh\nexit 0\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    db = tmp_path / "db"
    write_setup_manifest(db)
    return db


def _write(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def _gippa_run(
    tmp: Path, db: Path, proteins: str, abundance: str | None, peptides: str | None = None
) -> dict[str, str]:
    from click.testing import CliRunner

    from gippa.cli import cli
    from gippa.exitcodes import ExitCode

    rng = random.Random(11)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    genome = _write(tmp / "in" / "genome.fa", f">chr1 [topology=circular]\n{seq}\n")
    faa = _write(tmp / "in" / "proteins.faa", proteins)
    out = tmp / "OUT"
    argv = ["run", "--genome", str(genome), "--proteins", str(faa)]
    if abundance is not None:
        argv += ["--abundance", str(_write(tmp / "in" / "abundance.tsv", abundance))]
    if peptides is not None:
        argv += ["--peptides", str(_write(tmp / "in" / "peptides.tsv", peptides))]
    argv += ["--outdir", str(out), "--db-dir", str(db), "--modules", "core"]
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == ExitCode.SUCCESS, res.output
    params = json.loads((out / "work" / "params.json").read_text(encoding="utf-8"))
    return {k: v for k, v in params["inputs"].items() if v}


def _preflight_json(tmp: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((tmp / "OUT" / "preflight" / "preflight.json").read_text())
    return doc


def test_run_proteins_and_abundance_reach_nextflow(tmp_path: Path, run_env: Path) -> None:
    inputs = _gippa_run(
        tmp_path, run_env, f">P1\n{LONG_PROTEIN}\n", "protein_id\tvalue\nP1\t21.5\n"
    )
    work = (tmp_path / "OUT" / "work" / "preflight").absolute()
    assert Path(inputs["proteins_long"]) == work / "proteins.long.faa", inputs
    assert Path(inputs["abundance"]) == work / "abundance.norm.parquet", inputs
    assert "peptides" not in inputs, inputs  # no record of 40 aa or fewer
    # the join key: abundance protein IDs are the IDs of proteins.long.faa
    faa_ids = {ln[1:].split()[0] for ln in (work / "proteins.long.faa").read_text().splitlines()
               if ln.startswith(">")}  # fmt: skip
    assert set(pd.read_parquet(inputs["abundance"])["protein_id"]) <= faa_ids


# PR #537 review, lead decision (fix option a): without --peptides the short records
# (input.peptide_len, 7..40 aa) are not used (staging them flips the tier to two_peptide and
# undetects every long protein); preflight warns, naming the count. With --peptides they keep
# today's peptide path.
SHORT_RULE = "proteins.short_records_unused"
SHORT_2 = "MKVLAAGIVGRSTEQLLKAW"  # 20 aa
TOO_SHORT = "MKVLA"  # 5 aa: dropped (proteins.length_classes), not counted here
PROTEINS_WITH_SHORT = (
    f">P1\n{LONG_PROTEIN}\n>S1\n{SHORT_PROTEIN}\n>S2\n{SHORT_2}\n>T1\n{TOO_SHORT}\n"
)
USER_PEPTIDES = "peptide\tprotein_ids\trun\tintensity\nAAGLIEVQAPILSR\tP1\tr1\t50\n"


def _short_checks(tmp: Path) -> list[dict[str, Any]]:
    return [c for c in _preflight_json(tmp)["checks"] if c["rule_id"] == SHORT_RULE]


def test_run_proteins_only_does_not_stage_short_records(tmp_path: Path, run_env: Path) -> None:
    """Item 1: no ``peptides`` key in params.json; the run's tier stays ``protein_list``."""
    inputs = _gippa_run(tmp_path, run_env, PROTEINS_WITH_SHORT, None)
    assert "proteins_long" in inputs, inputs
    assert "peptides" not in inputs, f"short --proteins records staged as peptides: {inputs}"
    assert _preflight_json(tmp_path)["evidence_tier"] == "protein_list"


def test_run_proteins_only_warns_about_unused_short_records(tmp_path: Path, run_env: Path) -> None:
    """Item 2: one ``warn`` check on input ``proteins`` naming the count (2: S1, S2; the 5-aa T1 is
    dropped, not counted), saying the records are not used and how to pass them (--peptides)."""
    _gippa_run(tmp_path, run_env, PROTEINS_WITH_SHORT, None)
    checks = _short_checks(tmp_path)
    assert len(checks) == 1, _preflight_json(tmp_path)["checks"]
    (c,) = checks
    assert (c["level"], c["input"]) == ("warn", "proteins"), c
    assert re.search(r"\b2\b", c["message"]), c
    assert "--peptides" in c["message"], c
    assert "not used" in c["message"], c
    v = c["value"]
    assert v == 2 or (isinstance(v, dict) and 2 in v.values()), c


def test_run_proteins_only_long_records_no_short_warning(tmp_path: Path, run_env: Path) -> None:
    _gippa_run(tmp_path, run_env, f">P1\n{LONG_PROTEIN}\n", None)
    assert all(c["level"] != "warn" for c in _short_checks(tmp_path)), _short_checks(tmp_path)


def test_run_with_peptides_keeps_short_records_on_the_peptide_path(
    tmp_path: Path, run_env: Path
) -> None:
    """Item 3: with --peptides the short protein records join the normalised peptide table, as
    before, and no unused-records warning fires."""
    inputs = _gippa_run(tmp_path, run_env, PROTEINS_WITH_SHORT, None, peptides=USER_PEPTIDES)
    work = (tmp_path / "OUT" / "work" / "preflight").absolute()
    assert Path(inputs["peptides"]) == work / "peptides.norm.parquet", inputs
    got = set(pd.read_parquet(inputs["peptides"])["peptide"])
    assert got == {"AAGLIEVQAPILSR", SHORT_PROTEIN, SHORT_2}, got
    assert all(c["level"] != "warn" for c in _short_checks(tmp_path)), _short_checks(tmp_path)
