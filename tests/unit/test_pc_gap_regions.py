"""TU upstream regions over assembly gaps are not sent to Promoter Calculator (issue #614, 1.1.1).

Promoter Calculator's reverse complement knows only A, C, G, T, U, so a TU upstream region with an
``N`` (or any IUPAC code) raises ``KeyError`` in ANNOTATE:PROMOTER_CALCULATOR. Lead rule
(2026-11-22, the #611 / #612 pattern):

1. a TU whose upstream region (after strand-aware extraction, origin wrap included) holds any base
   outside ACGT (case-insensitive) gets no ``tu_upstream/<tu_id>.fa``; every other file is
   byte-identical to today's; ``write_tu_regions`` returns only the written paths;
2. ``operon_pairs`` logs exactly one WARNING with the number of skipped TUs (none without gaps);
3. a skipped TU has ``tx_rate_max`` / ``tx_rate_sum`` NA (not 0, which means "scored, no hit")
   and ``tss_distance_bp`` NA, through ``layers.parquet`` (``tx_rate_max_tu`` /
   ``tx_rate_sum_tu``) and the determinants matrix; a scored TU without a hit keeps 0;
4. a gap-free genome gives byte-identical outputs.

Seam (test author): the chain ``operon_pairs`` → fake Promoter Calculator (one CSV per written
FASTA, as the Nextflow module) → ``tx_layer`` → ``context_layers`` is run through ``run_step``,
each step getting every declared (optional) input that an earlier step produced. So the skipped
set may reach ``tx_layer`` either as "a TU of ``tus.parquet`` without a CSV" (no new key) or as a
new ``operon_pairs`` output the step declares, or ``context_layers`` may re-apply the ACGT rule to
the regions itself (it has ``tus`` and ``genome``). This file pins none of these, nor how (or
whether) ``tx_rates.parquet`` marks a skipped TU; it pins only that the scored rows (non-null
``tx_rate``) are exactly today's hits and that unscored rows add no between-gene evidence.
(``test_tx_layer_is_a_registered_py_step`` still pins the ``tx_layer`` inputs; a new key there
needs the lead.)
"""

from __future__ import annotations

import hashlib
import importlib
import logging
import random
import re
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import pytest
from helpers.promoter_calculator_fixtures import header_only, write_pc_csv
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.annotate.operons import Tu, region_sequence, write_tu_regions
from gippa.annotate.promoter_calculator import region_fasta_text
from gippa.config import GippaConfig
from gippa.model.coords import Interval

_ACGT = set("ACGT")


def _is_gap(seq: str) -> bool:
    return not set(seq.upper()) <= _ACGT


def _tu(tu_id: str, contig: str, strand: str, up: tuple[int, int]) -> Tu:
    # the gene interval is irrelevant to write_tu_regions: reuse the region
    return Tu(
        tu_id=tu_id,
        contig=contig,
        interval=Interval(contig, up[0], up[1], strand),  # type: ignore[arg-type]
        strand=strand,
        genes=(f"g_{tu_id}",),
        upstream=Interval(contig, up[0], up[1], strand),  # type: ignore[arg-type]
    )


# ------------------------------------------------------------- rule 1: write_tu_regions (unit)


def _unit_genome() -> dict[str, str]:
    rng = random.Random(614)
    c1 = list("".join(rng.choice("ACGT") for _ in range(2000)))
    c1[150] = "N"  # TU_N (+)
    c1[560] = "n"  # TU_n (+), lowercase
    c1[950] = "R"  # TU_R (-), IUPAC; complemented strand-aware
    c1[1350] = "y"  # TU_y (-), lowercase IUPAC
    c1[1450:1480] = list("acgt" * 7 + "ac")  # TU_soft (+): soft-masked, ACGT only -> kept
    c1[10] = "N"  # TU_wrap: the part of the region past the origin
    return {"c1": "".join(c1), "p1": "".join(rng.choice("ACGT") for _ in range(600))}


UNIT_TUS = [
    _tu("TU_1", "c1", "+", (60, 100)),
    _tu("TU_N", "c1", "+", (100, 200)),
    _tu("TU_2", "c1", "-", (200, 400)),
    _tu("TU_n", "c1", "+", (500, 600)),
    _tu("TU_R", "c1", "-", (900, 1000)),
    _tu("TU_3", "c1", "-", (1100, 1200)),
    _tu("TU_y", "c1", "-", (1300, 1400)),
    _tu("TU_soft", "c1", "+", (1440, 1500)),
    _tu("TU_wrap", "c1", "+", (1950, 2050)),  # unwrapped interval across the origin
    _tu("TU_wrap_ok", "c1", "-", (1990, 2005)),  # wraps, but 1990..2004 holds no gap
    _tu("TU_p", "p1", "+", (100, 400)),
]
UNIT_SKIPPED = {"TU_N", "TU_n", "TU_R", "TU_y", "TU_wrap"}


def test_unit_fixture_regions_are_what_the_test_claims() -> None:
    g = _unit_genome()
    got = {t.tu_id for t in UNIT_TUS if _is_gap(region_sequence(t.upstream, g[t.contig]))}
    assert got == UNIT_SKIPPED


def test_gap_regions_get_no_fasta_and_the_rest_are_unchanged(tmp_path: Path) -> None:
    g = _unit_genome()
    paths = write_tu_regions(tmp_path / "tu_upstream", UNIT_TUS, g)
    kept = [t for t in UNIT_TUS if t.tu_id not in UNIT_SKIPPED]
    assert paths == [tmp_path / "tu_upstream" / f"{t.tu_id}.fa" for t in kept], (
        "only the written paths, in TU order"
    )
    assert sorted(p.name for p in (tmp_path / "tu_upstream").iterdir()) == sorted(
        f"{t.tu_id}.fa" for t in kept
    )
    for t in kept:
        # today's bytes: one record, one sequence line (= region_fasta_text)
        assert (tmp_path / "tu_upstream" / f"{t.tu_id}.fa").read_text() == region_fasta_text(t, g)


@pytest.mark.parametrize(
    ("base", "strand"),
    [("N", "+"), ("n", "+"), ("R", "+"), ("r", "-"), ("Y", "-"), ("K", "+"), ("-", "+")],
)
def test_each_non_acgt_base_skips_the_tu(tmp_path: Path, base: str, strand: str) -> None:
    seq = list("ACGT" * 50)
    seq[120] = base
    g = {"c1": "".join(seq)}
    tus = [_tu("TU_a", "c1", strand, (100, 150)), _tu("TU_b", "c1", strand, (0, 60))]
    paths = write_tu_regions(tmp_path / "tu", tus, g)
    assert [p.name for p in paths] == ["TU_b.fa"]
    assert not (tmp_path / "tu" / "TU_a.fa").exists()


def test_all_acgt_regions_any_case_are_all_written(tmp_path: Path) -> None:
    g = {"c1": "acgtACGT" * 50}
    tus = [_tu("TU_a", "c1", "+", (0, 100)), _tu("TU_b", "c1", "-", (100, 300))]
    paths = write_tu_regions(tmp_path / "tu", tus, g)
    assert [p.name for p in paths] == ["TU_a.fa", "TU_b.fa"]
    for t in tus:
        assert (tmp_path / "tu" / f"{t.tu_id}.fa").read_text() == region_fasta_text(t, g)


_BASES = "ACGTacgtNnRYKMSWBDHV"


@settings(
    max_examples=60, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    genome=st.text(alphabet=_BASES, min_size=200, max_size=200)
    | st.text(alphabet="ACGTacgt", min_size=200, max_size=200),
    regions=st.lists(
        st.tuples(st.integers(0, 199), st.integers(1, 120), st.sampled_from(["+", "-"])),
        min_size=1,
        max_size=12,
    ),
)
def test_property_written_set_is_exactly_the_acgt_regions(
    tmp_path_factory: pytest.TempPathFactory, genome: str, regions: list[tuple[int, int, str]]
) -> None:
    g = {"c1": genome}
    tus = [_tu(f"TU_{i}", "c1", s, (a, a + n)) for i, (a, n, s) in enumerate(regions)]
    out = tmp_path_factory.mktemp("prop") / "tu"
    paths = write_tu_regions(out, tus, g)
    clean = [t for t in tus if not _is_gap(region_sequence(t.upstream, genome))]
    assert paths == [out / f"{t.tu_id}.fa" for t in clean]
    assert sorted(p.name for p in out.iterdir()) == sorted(p.name for p in paths)
    for t in clean:
        text = (out / f"{t.tu_id}.fa").read_text()
        assert text == region_fasta_text(t, g)
        assert set(text.split("\n")[1].upper()) <= _ACGT


# ----------------------------------------------------------------- step chain (rules 2, 3, 4)


@pytest.fixture
def reg() -> Iterator[ModuleType]:
    from gippa.steps import load_step_modules, registry

    with pytest.MonkeyPatch.context() as mp:
        names = load_step_modules()
        mp.setattr(registry, "_REGISTRY", {})
        for n in names:
            if n != "registry":
                importlib.reload(importlib.import_module(f"gippa.steps.{n}"))
        yield registry


#: (locus_tag, contig, start, end, strand); far apart, so the distance model makes singletons
GENES = [
    ("A", "c1", 100, 400, "+"),  # region [5800, 6150) wraps the origin; N at c1[50] -> skipped
    ("B", "c1", 1500, 1800, "+"),  # clean, hits
    ("C", "c1", 2600, 2900, "-"),  # region [2850, 3200); R at c1[3000] -> skipped
    ("D", "c1", 3800, 4100, "+"),  # N at c1[4000] inside the gene, outside the region; no hit
    ("E", "c1", 4800, 5100, "-"),  # clean, hits
    ("F", "p1", 1000, 1300, "+"),  # 30 N at p1[800:830] -> skipped
    ("G", "p1", 2000, 2300, "+"),  # clean, no hit
]
SKIPPED_GENES = {"A", "C", "F"}
HIT_GENES = {"B", "E"}
NO_HIT_GENES = {"D", "G"}
LENGTHS = {"c1": 6000, "p1": 3000}


def _genome(gaps: bool) -> dict[str, str]:
    rng = random.Random(614)
    c1 = list("".join(rng.choice("ACGT") for _ in range(LENGTHS["c1"])))
    p1 = list("".join(rng.choice("ACGT") for _ in range(LENGTHS["p1"])))
    if gaps:
        c1[50] = "N"
        c1[3000] = "R"
        c1[4000] = "N"
        p1[800:830] = list("N" * 30)
    return {"c1": "".join(c1), "p1": "".join(p1)}


def _write_inputs(d: Path, genome: dict[str, str]) -> dict[str, Path]:
    d.mkdir(parents=True)
    (d / "genome.fa").write_text("".join(f">{c}\n{s}\n" for c, s in genome.items()))
    pq.write_table(
        pa.table({"contig": ["c1", "p1"], "length": [LENGTHS["c1"], LENGTHS["p1"]],
                  "topology": ["circular", "linear"], "kind": ["chromosome", "plasmid"]}),
        d / "replicons.parquet",
    )  # fmt: skip
    n = len(GENES)
    cols = {
        "locus_tag": [g[0] for g in GENES],
        "feature_type": ["cds"] * n,
        "contig": [g[1] for g in GENES],
        "start": [g[2] for g in GENES],
        "end": [g[3] for g in GENES],
        "strand": [g[4] for g in GENES],
        "aa": ["M" * 100] * n,
        "nt": [None] * n,
        "start_codon": ["ATG"] * n,
    }
    pq.write_table(pa.table(cols), d / "genes.parquet")
    return {"genes": d / "genes.parquet", "replicons": d / "replicons.parquet",
            "genome": d / "genome.fa"}  # fmt: skip


def _outs(reg: ModuleType, name: str, d: Path, files: dict[str, str]) -> dict[str, Path]:
    """Every declared output: the known work-directory names, ``d / key`` for a new key."""
    return {k: d / files.get(k, k) for k in reg.get(name).outputs}


def _ins(reg: ModuleType, name: str, pool: dict[str, Path]) -> dict[str, Path]:
    spec = reg.get(name)
    missing = [k for k in spec.inputs if k not in pool]
    assert not missing, f"{name}: required inputs {missing} not produced by the chain"
    return {k: pool[k] for k in (*spec.inputs, *spec.optional_inputs) if k in pool}


_FILES = {
    "uniop_cds": "uniop_cds.gff", "uniop_map": "uniop_map.tsv", "pairs_raw": "pairs_raw.parquet",
    "tus": "tus.parquet", "tu_upstream": "tu_upstream", "operon_meta": "operon_meta.json",
    "tx_rates": "tx_rates.parquet", "annotation": "annotation", "layers": "layers.parquet",
}  # fmt: skip
#: 1-based region TSS of the fake hits (region length 350): both kept, on the TU's strand
HITS = [(251, "+", 1234.5), (281, "+", 30.0), (100, "-", 999.0)]


def _fake_pc(tu_upstream: Path, pc: Path, hit_tus: set[str]) -> list[Path]:
    """One CSV per written FASTA (the PROMOTER_CALCULATOR module); asserts the tool would run."""
    out = []
    for fa in sorted(tu_upstream.glob("*.fa")):
        seq = "".join(ln for ln in fa.read_text().splitlines() if not ln.startswith(">"))
        assert set(seq) <= set("ACGTU"), f"Promoter Calculator would raise KeyError on {fa.name}"
        csv = pc / f"{fa.stem}.csv"
        out.append(write_pc_csv(csv, HITS) if fa.stem in hit_tus else header_only(csv))
    return out


def _tu_of(tus: pd.DataFrame) -> dict[str, str]:
    return {str(g): str(t) for t, gs in zip(tus["tu_id"], tus["genes"], strict=True) for g in gs}


def _chain(reg: ModuleType, tmp: Path, genome: dict[str, str]) -> dict[str, Path]:
    cfg = GippaConfig(seed=1)
    pool: dict[str, Path] = dict(_write_inputs(tmp / "in", genome))
    out = tmp / "out"
    out.mkdir()
    for name in ("uniop_input", "operon_pairs"):
        outs = _outs(reg, name, out, _FILES)
        reg.run_step(name, cfg, _ins(reg, name, pool), outs)
        pool |= outs
    tu_of = _tu_of(pd.read_parquet(pool["tus"]))
    assert len(set(tu_of.values())) == len(GENES), "fixture: every gene is its own TU"
    csvs = _fake_pc(pool["tu_upstream"], out / "pc", {tu_of[g] for g in HIT_GENES})
    pool["pc"] = csvs  # type: ignore[assignment]  # a list, as Nextflow passes it
    for name in ("tx_layer", "context_layers"):
        outs = _outs(reg, name, out, _FILES)
        reg.run_step(name, cfg, _ins(reg, name, pool), outs)
        pool |= outs
    return pool


def _skip_warnings(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    """WARNINGs of the operon steps other than the distance-model fallback ones."""
    fallback = {"operon_pairs.run_flag", "operon_pairs.layer_missing"}
    return [
        r
        for r in caplog.records
        if r.levelno == logging.WARNING
        and r.name in ("gippa.steps.operon_pairs", "gippa.steps.uniop_input")
        and getattr(r, "event", None) not in fallback
    ]


def test_chain_skips_gap_tus_and_warns_once_with_the_count(
    tmp_path: Path, reg: ModuleType, caplog: pytest.LogCaptureFixture
) -> None:
    genome = _genome(gaps=True)
    with caplog.at_level(logging.INFO, logger="gippa.steps"):
        pool = _chain(reg, tmp_path, genome)  # _fake_pc asserts every FASTA is ACGT only
    tus = pd.read_parquet(pool["tus"])
    tu_of = _tu_of(tus)
    written = {p.stem for p in pool["tu_upstream"].glob("*.fa")}
    assert written == {tu_of[g] for g, *_ in GENES} - {tu_of[g] for g in SKIPPED_GENES}

    warns = _skip_warnings(caplog)
    assert len(warns) == 1, f"one WARNING for the skipped TUs, got {[r.message for r in warns]}"
    assert re.search(rf"\b{len(SKIPPED_GENES)}\b", warns[0].getMessage())


def test_chain_written_regions_are_todays_bytes(tmp_path: Path, reg: ModuleType) -> None:
    genome = _genome(gaps=True)
    pool = _chain(reg, tmp_path, genome)
    tus = pd.read_parquet(pool["tus"])
    g = {c: s.upper() for c, s in genome.items()}
    for r in tus.itertuples(index=False):
        fa = pool["tu_upstream"] / f"{r.tu_id}.fa"
        if not fa.exists():
            continue
        tu = Tu(str(r.tu_id), str(r.contig), Interval(r.contig, r.start, r.end, r.strand),
                str(r.strand), tuple(r.genes),
                Interval(r.contig, r.upstream_start, r.upstream_end, r.strand))  # fmt: skip
        assert fa.read_text() == region_fasta_text(tu, g)


def test_skipped_tu_is_na_and_scored_no_hit_tu_is_zero_in_layers(
    tmp_path: Path, reg: ModuleType
) -> None:
    pool = _chain(reg, tmp_path, _genome(gaps=True))
    lay = pd.read_parquet(pool["layers"]).set_index("locus_tag")
    for g in SKIPPED_GENES:
        assert pd.isna(lay.loc[g, "tx_rate_max_tu"]), f"{g}: skipped TU, not 'scored, no hit'"
        assert pd.isna(lay.loc[g, "tx_rate_sum_tu"]), g
        assert pd.isna(lay.loc[g, "tss_distance_bp"]), g
    for g in NO_HIT_GENES:
        assert lay.loc[g, "tx_rate_max_tu"] == 0.0, g
        assert lay.loc[g, "tx_rate_sum_tu"] == 0.0, g
        assert pd.isna(lay.loc[g, "tss_distance_bp"]), g
    for g in HIT_GENES:
        assert lay.loc[g, "tx_rate_max_tu"] == pytest.approx(1234.5), g
        assert lay.loc[g, "tx_rate_sum_tu"] == pytest.approx(1264.5), g
        assert not pd.isna(lay.loc[g, "tss_distance_bp"]), g


def test_tx_rates_scored_rows_are_exactly_the_hits(tmp_path: Path, reg: ModuleType) -> None:
    pool = _chain(reg, tmp_path, _genome(gaps=True))
    tx = pd.read_parquet(pool["tx_rates"])
    assert list(tx.columns) == ["contig", "strand", "tss", "tx_rate", "tu_id"]
    tu_of = _tu_of(pd.read_parquet(pool["tus"]))
    scored = tx[tx["tx_rate"].notna()]
    assert set(scored["tu_id"]) == {tu_of[g] for g in HIT_GENES}
    assert len(scored) == 2 * len(HIT_GENES)
    assert scored["tss"].notna().all()
    assert tx["tx_rate"].notna().all(), "a skipped TU adds no marker row"


def test_skipped_tu_does_not_break_pair_evidence_or_determinants(
    tmp_path: Path, reg: ModuleType
) -> None:
    from gippa.annotate.layers import genes_of, replicons_of, tus_of
    from gippa.determinants.features import feature_frame, model_matrix
    from gippa.operons.pairs import Evidence, pair_frame
    from gippa.steps.operons import read_layer

    pool = _chain(reg, tmp_path, _genome(gaps=True))
    cfg = GippaConfig(seed=1)
    rep_df = pd.read_parquet(pool["replicons"])
    reps = replicons_of(rep_df)
    tus = pd.read_parquet(pool["tus"])
    genes = genes_of(pd.read_parquet(pool["genes"]))
    pairs_raw = pd.read_parquet(pool["pairs_raw"])
    n = len(GENES)

    # operons step: the tx_rates layer as it reads it; unscored rows add no between-gene evidence
    tx = read_layer(cfg, "tx_rates", pool["tx_rates"], reps, tus)
    gtbl = pd.DataFrame({"locus_tag": [g[0] for g in GENES], "detected": [True] * n,
                         "p_detect": [0.9] * n, "abundance_log2": [20.0] * n})  # fmt: skip
    rna = pd.DataFrame({"contig": [], "start": [], "end": [], "type": []})

    def frame(t: pd.DataFrame) -> pd.DataFrame:
        return pair_frame(gtbl, genes, reps, pairs_raw, Evidence(None, None, t, rna), cfg)

    full = frame(tx)
    hits_only = frame(tx[tx["tx_rate"].notna()].astype({"tss": "int64"}))
    assert full["tx_rate_between"].notna().any(), "fixture: some pair has tx evidence"
    pd.testing.assert_series_equal(full["tx_rate_between"], hits_only["tx_rate_between"])

    # determinants: the gene table copies the layers; NA stays NA, 0 stays 0, nothing crashes
    lay = pd.read_parquet(pool["layers"])
    g_df = pd.read_parquet(pool["genes"])
    tbl = g_df[["locus_tag", "contig", "start", "end", "strand"]].merge(lay, on="locus_tag")
    tbl = tbl.assign(partial=False, p_detect=0.9, detected=True, abundance_log2=20.0,
                     abundance_imputed=False, replicon_type="chromosome")  # fmt: skip
    feats = feature_frame(tbl, genes, tus_of(tus), lay, None, None, reps).set_index("locus_tag")
    mm, _ = model_matrix(feats.reset_index(drop=True))
    mm.index = feats.index
    for col in ("tx_rate_max_tu", "tx_rate_sum_tu"):
        assert col in mm.columns, "a partly-NA feature is kept"
        assert mm.loc[sorted(SKIPPED_GENES), col].isna().all(), col
        assert (mm.loc[sorted(NO_HIT_GENES), col] == 0.0).all(), col
        assert np.isfinite(mm.loc[sorted(HIT_GENES), col]).all(), col


def test_gap_free_chain_logs_no_skip_warning(
    tmp_path: Path, reg: ModuleType, caplog: pytest.LogCaptureFixture
) -> None:
    with caplog.at_level(logging.INFO, logger="gippa.steps"):
        pool = _chain(reg, tmp_path, _genome(gaps=False))
    assert _skip_warnings(caplog) == []
    assert len(list(pool["tu_upstream"].glob("*.fa"))) == len(GENES)
    lay = pd.read_parquet(pool["layers"])
    assert not lay["tx_rate_max_tu"].isna().any()


# Recorded on main cff0fd68 (before the fix): a gap-free genome must stay byte-identical.
GAP_FREE_SHA256: dict[str, str] = {
    "tu_upstream/TU_1.fa": "1772f699f37bb94b66e70df7a9aedb3de225fb7987b299dff6c6bcbbadd0706f",
    "tu_upstream/TU_2.fa": "b86fca52fe082a6ee90b47109a701fbd571a2fd10cb2b5e7e5c417caee5dcc36",
    "tu_upstream/TU_3.fa": "05ca885d07ded3d5e2b95f7ee997216259f8d6fa190f6984d99f4c4b062422fb",
    "tu_upstream/TU_4.fa": "bf2bb619a950e34092296dfbd8c465cfb10e67bec91b65e5efe8318e2037e734",
    "tu_upstream/TU_5.fa": "78c9a2f4c8199a0cb5950f5d424bcebd7db51e7304ab473998eb6cd606a10269",
    "tu_upstream/TU_6.fa": "3590fdc3f252f30c86c9b47dc5e737fa4f53f5b5c934db5f8a6d7c420b5bee23",
    "tu_upstream/TU_7.fa": "21a3cfa1ad87fd9943173b53cb3fc1a048886ee04c2b1ce9e3ab71ec3958314d",
    "tus.parquet": "36b5362d1217dc181541d602191773c97634e135f11e107747d2284c68052988",
    "pairs_raw.parquet": "74e24b9f977e2a0beb66c63c05831a66c35d38277a3cc26dcc0d296f5dd005c7",
    "tx_rates.parquet": "a377c24f21dd4d3afd733669b8e0e4814ac26283f420c05021cb2ad444ad8736",
    "layers.parquet": "d4f6a9e52069cb9be091272c113ef49c88cb0b3dcaf35cbd7fa040844ac61acf",
}


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _gap_free_hashes(pool: dict[str, Path]) -> dict[str, str]:
    got = {f"tu_upstream/{p.name}": _sha(p) for p in sorted(pool["tu_upstream"].glob("*"))}
    for key in ("tus", "pairs_raw", "tx_rates", "layers"):
        got[_FILES[key]] = _sha(pool[key])
    return got


def test_gap_free_genome_outputs_are_byte_identical(tmp_path: Path, reg: ModuleType) -> None:
    pool = _chain(reg, tmp_path, _genome(gaps=False))
    assert _gap_free_hashes(pool) == GAP_FREE_SHA256
