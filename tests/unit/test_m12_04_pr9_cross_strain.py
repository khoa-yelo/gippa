"""M12-04 PR 9, T-02: ``cross_strain_protein`` only on non-identical hits with a foreign ID; the
run flag ``foreign_protein_ids`` (docs/specs/M12-04-outputs-revision.md > PR 9, AC 9.2-9.4, § 5.4;
ADR-0075).

Rule (spec > PR 9 Behaviour, from the review's correction: the flag comes from
``gippa.mapping.loci.annotate_protein_hits``, "``protein_id`` is not a locus tag of the model"):
a protein hit gets ``cross_strain_protein`` exactly when its ``protein_id`` is not a locus tag of
the model **and** its ``class`` is not ``identical``. When at least one protein ID is not a locus
tag, the run flag ``foreign_protein_ids`` goes into ``run_manifest.json > flags`` once, whatever the
class. The gene flag in ``genes.tsv`` follows the hit flag, as today. Route B of the example (all
106 hits identity 1, coverage 1, class identical, ``JNOIEE_*`` IDs against ``EOIENK_*`` genes)
then carries no row flag.

Pinned here (test author):

- Fixture R (AC 9.2) is the recorded miniprot 0.18 run on SIM100K (tests/data/recorded/miniprot,
  the issue #536 proteins-only fixture) with every ``GSIM_`` protein ID renamed ``JNOIEE_`` in
  ``proteins.faa`` and ``miniprot.gff``: same sequences and alignments, foreign IDs. ``loci`` runs
  through ``gippa.steps.registry.run_step`` (tests/helpers/mapping_steps.py ``drive``; tblastn
  finds nothing). The control is the same run with the original locus-tag IDs.
- The run flag reaches the manifest the way every step run flag does: a log record of the
  ``loci`` step (logger ``gippa.steps.loci``) with an ``event`` and a ``data`` mapping holding
  ``"flag": "foreign_protein_ids"`` or a ``"flags"`` list with it; ``gippa.manifest`` collects
  ``data.flag``/``data.flags`` of the step events. The event name is the implementer's choice.
- "No gene has ``cross_strain_protein``" is checked through ``gippa.genetable.abundance
  .protein_abundance`` (the only place the gene flag comes from, M4-14 decision 2) fed with the
  written ``protein_hits.tsv``.
- AC 9.4: §9.4's Run row lives in code as ``gippa.manifest.RUN_FLAGS`` (GIPPA_SPEC.md is
  read-only; ADR-0075 records the token); the manifest schema's ``flags`` pattern must accept it;
  ``gippa.viz.words.FLAG_TEXT`` (M12-04 PR 0) holds the spec's sentence verbatim.

No statistical claim is made (a deterministic flag rule), so there is no calibration test (§8.4).
"""

from __future__ import annotations

import dataclasses
import hashlib
import logging
import sys
from collections.abc import Iterator, Mapping, Sequence
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))

CROSS = "cross_strain_protein"
TWO = "two_proteins_one_gene"
FOREIGN = "foreign_protein_ids"
FOREIGN_TEXT = "proteins come from a related annotation; identical matches are not flagged per gene"
CLASSES = ("identical", "near", "divergent")
TAGS = ("EOIENK_00005", "EOIENK_00010", "EOIENK_00015", "EOIENK_00020")
SETTINGS = settings(max_examples=100, deadline=None, suppress_health_check=[HealthCheck.too_slow])


def flags_of(x: object) -> set[str]:
    if x is None or (not isinstance(x, str | list | tuple) and pd.isna(x)):
        return set()
    if isinstance(x, list | tuple):
        return {str(f) for f in x}
    return {f for f in str(x).split(";") if f and f != "NA"}


def model(tags: Sequence[str] = TAGS) -> Any:
    """The part of a ``GeneModel`` ``annotate_protein_hits`` reads: the genes' locus tags."""
    return SimpleNamespace(genes=tuple(SimpleNamespace(locus_tag=t) for t in tags))


def hit(pid: str, gene: str | None, cls: str, rank: int = 1, flags: str = "") -> dict[str, Any]:
    ident = {"identical": 1.0, "near": 0.97, "divergent": 0.8}[cls]
    return {
        "protein_id": pid, "protein_len_aa": 99, "route": "miniprot", "contig": "chr1",
        "start": 100 * rank, "end": 100 * rank + 300, "strand": "+", "rank": rank,
        "identity": ident, "query_coverage": 1.0, "reciprocal_overlap": 1.0, "class": cls,
        "gene": gene, "flags": flags,
    }  # fmt: skip


def annotate(rows: Sequence[Mapping[str, Any]], tags: Sequence[str] = TAGS) -> list[set[str]]:
    from gippa.mapping.loci import annotate_protein_hits

    df = annotate_protein_hits(pd.DataFrame(list(rows)), model(tags))
    assert list(df["protein_id"]) == [r["protein_id"] for r in rows]
    return [flags_of(f) for f in df["flags"]]


# ---------------------------------------------------------------- the hit rule (unit)


def test_route_b_foreign_identical_hits_are_not_flagged() -> None:
    """Spec example: B's 106 hits are all identity 1, coverage 1, class identical, with IDs of a
    related annotation; after PR 9 no row carries ``cross_strain_protein``."""
    rows = [hit(f"JNOIEE_{i:05d}", t, "identical") for i, t in enumerate(TAGS)]
    assert annotate(rows) == [set()] * len(TAGS)


def test_foreign_near_is_flagged() -> None:
    """Edge case and AC 9.3: a foreign ID with ``class = near`` keeps the flag."""
    assert annotate([hit("WP_000001", TAGS[0], "near")]) == [{CROSS}]


def test_foreign_divergent_is_flagged() -> None:
    assert annotate([hit("WP_000002", TAGS[1], "divergent")]) == [{CROSS}]


def test_locus_tag_divergent_is_not_flagged() -> None:
    """Edge case: the same annotation, so no cross-strain claim."""
    assert annotate([hit(TAGS[1], TAGS[1], "divergent"), hit(TAGS[2], TAGS[2], "near")]) == [
        set(),
        set(),
    ]


def test_two_proteins_one_gene_is_unchanged_on_identical_foreign_hits() -> None:
    """Only ``cross_strain_protein`` changes: two input proteins on one gene keep
    ``two_proteins_one_gene`` (spec example 15 of M3-14, now with an identical foreign hit)."""
    out = annotate([hit("P0A7V8", TAGS[0], "identical"), hit(TAGS[0], TAGS[0], "identical")])
    assert out == [{TWO}, {TWO}]


def test_input_flags_are_kept() -> None:
    out = annotate([hit("WP_9", TAGS[0], "identical", flags="fs;st"),
                    hit("WP_8", TAGS[1], "near", flags="fs")])  # fmt: skip
    assert out == [{"fs", "st"}, {"fs", CROSS}]


@SETTINGS
@given(
    rows=st.lists(
        st.tuples(
            st.one_of(st.sampled_from(TAGS), st.sampled_from(("JNOIEE_1", "JNOIEE_2", "WP_3"))),
            st.one_of(st.none(), st.sampled_from(TAGS)),
            st.sampled_from(CLASSES),
        ),
        min_size=1,
        max_size=8,
    )
)
def test_flag_iff_foreign_and_not_identical(rows: list[tuple[str, str | None, str]]) -> None:
    """Property: ``cross_strain_protein`` ⇔ (ID not a locus tag) ∧ (class ≠ identical), on every
    row, whatever the other rows hold."""
    recs = [hit(p, g, c, rank=i + 1) for i, (p, g, c) in enumerate(rows)]
    for (p, _, c), fl in zip(rows, annotate(recs), strict=True):
        assert (CROSS in fl) == (p not in TAGS and c != "identical"), (p, c, fl)


def has_foreign(rows: Sequence[Mapping[str, Any]], tags: Sequence[str] = TAGS) -> bool:
    from gippa.mapping.loci import has_foreign_protein_ids

    return bool(has_foreign_protein_ids(pd.DataFrame(list(rows)), model(tags)))


def test_run_flag_on_mixed_ids() -> None:
    """One foreign ID among locus-tag IDs raises the run flag, whatever its class (it is any, not
    all, of the protein IDs)."""
    assert has_foreign([hit(TAGS[0], TAGS[0], "identical"),
                        hit("WP_1", TAGS[1], "identical", rank=2)])  # fmt: skip
    assert has_foreign([hit("WP_1", TAGS[1], "identical"), hit(TAGS[0], TAGS[0], "near", rank=2),
                        hit(TAGS[2], TAGS[2], "divergent", rank=3)])  # fmt: skip


def test_run_flag_off_on_locus_tag_ids() -> None:
    assert not has_foreign([hit(t, t, c, rank=i + 1) for i, (t, c) in
                            enumerate(zip(TAGS, ("identical", "near", "divergent", "near"),
                                          strict=True))])  # fmt: skip


@SETTINGS
@given(
    rows=st.lists(
        st.tuples(
            st.one_of(st.sampled_from(TAGS), st.sampled_from(("JNOIEE_1", "JNOIEE_2", "WP_3"))),
            st.sampled_from(CLASSES),
        ),
        min_size=1,
        max_size=8,
    )
)
def test_run_flag_iff_any_foreign_id(rows: list[tuple[str, str]]) -> None:
    """Property: the run flag ⇔ at least one protein ID is not a locus tag of the model."""
    recs = [hit(p, TAGS[0], c, rank=i + 1) for i, (p, c) in enumerate(rows)]
    assert has_foreign(recs) == any(p not in TAGS for p, _ in rows), rows


# ---------------------------------------------------------------- fixture R (AC 9.2, 9.3)


@pytest.fixture(scope="module")
def sim100k(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from helpers.mapping_steps import preflight

    from gippa.config import resolve_config
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SIM100K, SIM100K_SEED, simulate

    tmp = tmp_path_factory.mktemp("sim100k_pr9")
    res = simulate(dataclasses.replace(SIM100K, proteomics=ProteomicsConfig()), SIM100K_SEED,
                   tmp / "sim")  # fmt: skip
    sha = hashlib.sha256(res.genome_fasta.read_bytes()).hexdigest()
    assert sha.startswith("dd10216e"), "SIM100K genome differs from the recorded miniprot fixture"
    run = preflight(res, tmp / "pf")
    return {"res": res, "run": run, "cfg": resolve_config(run.config_yaml, ["seed=3"])}


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _renamed(src: Path, dst: Path, old: str, new: str) -> Path:
    dst.write_text(src.read_text(encoding="utf-8").replace(old, new), encoding="utf-8")
    return dst


def _drive_loci(sim: dict[str, Any], tmp: Path, prefix: str | None) -> dict[str, Any]:
    """``loci`` on proteins only; ``prefix`` replaces ``GSIM_`` in every protein ID (None keeps
    the locus tags). Returns the written ``protein_hits.tsv`` frame and the step's log records."""
    from helpers.mapping_steps import RECORDED_MINIPROT, drive, fresh_steps, proteins_key

    tmp.mkdir(parents=True, exist_ok=True)
    old, new = "GSIM_", prefix or "GSIM_"
    faa = _renamed(RECORDED_MINIPROT / "proteins.faa", tmp / "proteins.long.faa", old, new)
    gff = _renamed(RECORDED_MINIPROT / "miniprot.gff", tmp / "miniprot.gff", old, new)
    tbl = tmp / "tblastn.tsv"
    tbl.write_text("", encoding="utf-8")
    rec = _Recorder()
    log = logging.getLogger("gippa.steps")
    log.addHandler(rec)
    old_level = log.level
    log.setLevel(logging.DEBUG)
    try:
        with pytest.MonkeyPatch.context() as mp:
            reg = fresh_steps(mp)
            run, res = sim["run"], sim["res"]
            roots = {"genome": run.genome, "bakta_json": res.bakta_json, proteins_key(reg): faa}
            (mp_key,) = reg.get("miniprot").outputs
            (tb_key,) = reg.get("tblastn").outputs
            avail = drive(reg, sim["cfg"], tmp / "work", roots,
                          inject={"miniprot": {mp_key: gff}, "tblastn": {tb_key: tbl}},
                          targets=("loci",))  # fmt: skip
    finally:
        log.removeHandler(rec)
        log.setLevel(old_level)
    path = Path(avail["mapping"]) / "protein_hits.tsv"
    hits = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, na_values=["NA"])
    loci_records = [r for r in rec.records if r.name == "gippa.steps.loci"]
    return {"hits": hits, "records": loci_records, "tags": _model_tags(res.bakta_json)}


def _model_tags(bakta_json: Path) -> frozenset[str]:
    import json

    doc = json.loads(Path(bakta_json).read_text(encoding="utf-8"))
    return frozenset(str(f["locus"]) for f in doc["features"] if f.get("locus"))


@pytest.fixture(scope="module")
def fixture_r(sim100k: dict[str, Any], tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _drive_loci(sim100k, tmp_path_factory.mktemp("fixture_r"), "JNOIEE_")


@pytest.fixture(scope="module")
def control(sim100k: dict[str, Any], tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _drive_loci(sim100k, tmp_path_factory.mktemp("control"), None)


def _event_flags(records: Sequence[logging.LogRecord]) -> Iterator[str]:
    for r in records:
        data = getattr(r, "data", None)
        if getattr(r, "event", None) is None or not isinstance(data, Mapping):
            continue
        if isinstance(data.get("flag"), str):
            yield data["flag"]
        if isinstance(data.get("flags"), list):
            yield from (f for f in data["flags"] if isinstance(f, str))


def test_fixture_r_has_foreign_ids_and_identical_hits(fixture_r: dict[str, Any]) -> None:
    hits = fixture_r["hits"]
    assert len(hits) > 0
    assert hits["protein_id"].str.startswith("JNOIEE_").all()
    assert not set(hits["protein_id"]) & fixture_r["tags"]
    # every recorded alignment is identity 1 / coverage 1: the route B case (AC 9.3 is unit-tested)
    assert (hits["class"] == "identical").all(), hits["class"].value_counts().to_dict()


def test_ac_9_2_no_identical_hit_carries_the_flag(fixture_r: dict[str, Any]) -> None:
    hits = fixture_r["hits"]
    ident = hits[hits["class"] == "identical"]
    bad = [(p, f) for p, f in zip(ident["protein_id"], ident["flags"], strict=True)
           if CROSS in flags_of(f)]  # fmt: skip
    assert not bad, f"{len(bad)} identical foreign hits carry {CROSS}, e.g. {bad[:3]}"


def test_ac_9_2_no_gene_carries_the_flag(fixture_r: dict[str, Any]) -> None:
    """The gene flag follows the hit flag: genes whose rank-1 protein hits are all identical get
    no ``cross_strain_protein`` from ``protein_abundance``."""
    from gippa.genetable.abundance import protein_abundance

    hits = fixture_r["hits"].copy()
    hits["rank"] = hits["rank"].astype(int)
    hits["flags"] = [sorted(flags_of(f)) for f in hits["flags"]]
    top = hits[(hits["rank"] == 1) & hits["gene"].notna()]
    pmap = pd.DataFrame({"protein_id": top["protein_id"].astype(object),
                         "gene": top["gene"].astype(object),
                         "method": ["protein_hit"] * len(top)})  # fmt: skip
    ab = pd.DataFrame({"protein_id": top["protein_id"].tolist(),
                       "value": [20.0 + 0.01 * i for i in range(len(top))]})  # fmt: skip
    genes = protein_abundance(ab, pmap, "maxlfq", "log2", protein_hits=hits)
    assert len(genes) > 0
    by_gene = hits.groupby("gene")["class"].agg(lambda c: set(c))
    identical_only = {str(g) for g, cl in by_gene.items() if cl == {"identical"}}
    bad = [t for t, f in zip(genes["locus_tag"], genes["flags"], strict=True)
           if t in identical_only and CROSS in flags_of(f)]  # fmt: skip
    assert not bad, f"{len(bad)} genes with only identical hits carry {CROSS}: {bad[:5]}"


def test_ac_9_2_loci_emits_the_run_flag(fixture_r: dict[str, Any]) -> None:
    found = list(_event_flags(fixture_r["records"]))
    assert FOREIGN in found, f"loci logged no {FOREIGN} run flag; event flags {found}"


def test_ac_9_2_manifest_flags_hold_the_run_flag_once(fixture_r: dict[str, Any]) -> None:
    from helpers.manifest_data import build, event, validate

    evs = [event("loci", str(r.event), dict(r.data), level=r.levelname)  # type: ignore[attr-defined]
           for r in fixture_r["records"]
           if getattr(r, "event", None) is not None and isinstance(getattr(r, "data", None),
                                                                     Mapping)]  # fmt: skip
    m = build(events=evs + evs)  # a resumed step logs its events again
    assert m["flags"].count(FOREIGN) == 1, m["flags"]
    validate(m)


def test_locus_tag_ids_raise_no_run_flag(control: dict[str, Any]) -> None:
    hits = control["hits"]
    assert len(hits) > 0 and set(hits["protein_id"]) <= control["tags"]
    assert FOREIGN not in set(_event_flags(control["records"]))
    assert not any(CROSS in flags_of(f) for f in hits["flags"])


# ---------------------------------------------------------------- AC 9.4 vocabulary


def test_ac_9_4_run_vocabulary_has_the_token() -> None:
    from gippa.manifest import RUN_FLAGS

    assert FOREIGN in RUN_FLAGS


def test_ac_9_4_manifest_schema_accepts_the_token() -> None:
    from helpers.manifest_data import build, event, schema_errors

    m = build(events=[event("loci", "loci.run_flag", {"flag": FOREIGN}, level="WARNING")])
    assert FOREIGN in m["flags"]
    assert not schema_errors(m)


def test_ac_9_4_flag_text_has_the_sentence() -> None:
    from gippa.viz.words import FLAG_TEXT

    assert FLAG_TEXT[FOREIGN] == FOREIGN_TEXT
