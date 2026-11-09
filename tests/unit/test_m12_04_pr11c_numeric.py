"""M12-04 PR 11c: numeric ``variants.class_error_rate`` and ``novel_orfs.class_fdr`` (T-05).

docs/specs/M12-04-outputs-revision.md > PR 11c (behaviour, examples, edge cases, AC 11c.1-11c.3),
§ 5.3 (type changes) and § 5.4 (flag changes); Q3 (1.0 output folders are mapped on read).

Interface pinned here (test-author reading, beyond the spec text):

- Schemas: ``class_error_rate`` (variants) and ``class_fdr`` (novel_orfs) have ``type: number``
  and are not ``required`` (NA is a legal value); the variants ``flags`` pattern accepts
  ``class_error_unestimated``.
- Producers: ``gippa.variants.error_rates.annotate_sites`` writes NaN (not the word) for a class
  without observations; ``gippa.novel.class_call.class_fdr_cell`` returns NaN for an unestimated
  ORF (still 1.0 for ``lost_pick``); ``gippa.novel.candidates.minimal_call`` writes NaN.
- ``gippa.variants.vcf_records.vcf_records`` builds the ``variants.tsv`` frame: its
  ``class_error_rate`` cells are floats (NaN/NA when unestimated), and a row's ``flags`` tuple
  holds ``class_error_unestimated`` if and only if the rate is NA (added when missing, dropped
  when stale), after the other Variants flags (the K VAR_25 example). ``CLASS_ERR`` is ``.`` for
  an NA allele.
- ``gippa.novel.outputs.orf_table`` gives a ``float64`` ``class_fdr`` column; NA iff the row
  has ``class_fdr_unestimated`` (the existing rule, now keyed on NA).
- ``write_variants``/``write_novel_orfs`` write NA as ``NA``; ``read_variants``/
  ``read_novel_orfs`` give a float column back.
- ``RunView.table`` on a 1.0 tree whose cell holds ``unestimated`` returns NaN and adds the
  matching flag when it is absent (no duplicate when present).
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Sequence
from functools import cache
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.m12_04_renames import to_1_0
from helpers.novel_outputs_fixture import (
    novel_io,
    orfs_of,
    scored,
    table,
    two_contigs,
)
from helpers.variant_tables_fixture import (
    GENOME,
    GM_,
    GP,
    REPS,
    body_lines,
    info_dict,
    peptides_at,
    records,
    run,
    sub,
    vr,
    vt,
    write,
)
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
TABLES = REPO / "src" / "gippa" / "schemas" / "tables"
CE, CE_FLAG = "class_error_rate", "class_error_unestimated"
FDR, FDR_FLAG = "class_fdr", "class_fdr_unestimated"
#: K VAR_25 (spec PR 11c example): the 1.0 flags, then the 1.1 flags
K_VAR25_FLAGS = ("mass_equivalent", "snv", "strain_difference", "database_derived", "start_codon")
K_VAR25_1_1 = "mass_equivalent;snv;strain_difference;database_derived;start_codon;" + CE_FLAG

#: every fixture and golden holding one of the two tables (AC 11c.1, 11c.2)
FIXTURES: tuple[tuple[str, str], ...] = (
    ("tests/data/golden/sim50k_m7/variants.tsv", "variants"),
    ("tests/data/golden/sim50k_m7/novel_orfs.tsv", "novel_orfs"),
    ("tests/data/snapshot/sim_plots/OUT/tables/variants.tsv", "variants"),
    ("tests/data/snapshot/sim_plots/OUT/tables/novel_orfs.tsv", "novel_orfs"),
    ("tests/data/snapshot/sim_empty/OUT/tables/variants.tsv", "variants"),
    ("tests/data/snapshot/sim_empty/OUT/tables/novel_orfs.tsv", "novel_orfs"),
)
#: the goldens hold numbers only, so 11c leaves them byte-identical (spec example)
GOLDEN_SHA256 = {
    # M12-05b appends nine gene-context columns; the 26 earlier columns are byte-identical
    "tests/data/golden/sim50k_m7/novel_orfs.tsv":
        "08a29207c88dcfdc56ea76b41b69995011e169ebdbdba94701ab2f914b3f5065",
    # PR 10 (putative_start) caps the TTG->M row, so this golden's content changes there
    "tests/data/golden/sim50k_m7/variants.tsv":
        "22c3d32ce87cdc899ad594a63ec3835aa79adec5c30226328bb91370f00426f2",
}  # fmt: skip
COLUMN = {"variants": (CE, CE_FLAG), "novel_orfs": (FDR, FDR_FLAG)}


def _schema(stem: str) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((TABLES / f"{stem}.schema.json").read_text(encoding="utf-8"))
    return doc


def _field(stem: str, name: str) -> dict[str, Any]:
    (f,) = [f for f in _schema(stem)["fields"] if f["name"] == name]
    return dict(f)


def _na(v: object) -> bool:
    if v is None or v is pd.NA:
        return True
    return isinstance(v, float) and math.isnan(v)


def _flag_set(v: object) -> set[str]:
    if isinstance(v, list | tuple | set | frozenset):
        return {str(x) for x in v}
    if _na(v) or v == "":
        return set()
    return set(str(v).split(";"))


def _is_number_or_na(v: object) -> bool:
    return _na(v) or (isinstance(v, int | float) and not isinstance(v, bool))


# ---------------------------------------------------------------- schemas (AC 11c.1, 11c.3)


@pytest.mark.parametrize(("stem", "col"), [("variants", CE), ("novel_orfs", FDR)])
def test_column_has_schema_type_number_and_admits_na(stem: str, col: str) -> None:
    f = _field(stem, col)
    assert f["type"] == "number", f
    assert not (f.get("constraints") or {}).get("required"), f"{stem}.{col} must admit NA: {f}"
    assert "NA" in _schema(stem)["missingValues"]
    assert "unestimated" not in json.dumps(f.get("constraints") or {}), f


def _flags_pattern(stem: str) -> re.Pattern[str]:
    return re.compile(_field(stem, "flags")["constraints"]["pattern"])


def test_variants_flags_pattern_has_class_error_unestimated() -> None:
    pat = _flags_pattern("variants")
    assert pat.fullmatch(CE_FLAG)
    assert pat.fullmatch(K_VAR25_1_1)
    assert pat.fullmatch(f"vcf_unrepresentable;{CE_FLAG}")
    assert not pat.fullmatch(FDR_FLAG)  # the novel ORF flag stays out of the Variants vocabulary


def test_novel_orfs_flags_pattern_keeps_class_fdr_unestimated() -> None:
    assert _flags_pattern("novel_orfs").fullmatch(f"{FDR_FLAG};no_start_codon")


def test_class_error_unestimated_is_in_the_variants_vocabulary_constants() -> None:
    assert CE_FLAG in vr().VARIANT_FLAGS
    from gippa.variants import flags as vflags

    assert CE_FLAG in vflags.FLAG_ORDER


def test_flag_text_has_class_error_unestimated_in_plain_words() -> None:
    from gippa.viz.words import FLAG_TEXT, plain

    assert CE_FLAG in FLAG_TEXT
    text = plain(CE_FLAG)
    assert text and text != CE_FLAG and "_" not in text, text


def test_flags_reference_page_lists_the_new_flag_under_variants() -> None:
    """The §9.4 docs table is ``docs/reference/flags.md`` (``gen_reference.py render_flags``)."""
    text = (REPO / "docs" / "reference" / "flags.md").read_text(encoding="utf-8")
    m = re.search(r"^## `variants`\n(.*?)(?=^## )", text, re.MULTILINE | re.DOTALL)
    assert m, "no variants section in docs/reference/flags.md"
    assert f"- `{CE_FLAG}`" in m.group(1), m.group(1)


# ---------------------------------------------------------------- producers


def test_annotate_sites_writes_nan_for_a_class_without_observations() -> None:
    from gippa.variants.error_rates import ClassRate, annotate_sites

    sites = pd.DataFrame({
        "variant_id": ["VAR_1", "VAR_2", "VAR_3"],
        "alt_type": ["substitution", "swap", "substitution"],
        "flags": [("snv",), ("mass_equivalent",), ("mass_equivalent",)],
    })  # fmt: skip
    rates = [ClassRate("substitution", False, 4, 0, 0.0), ClassRate("swap", True, 2, 2, 1.0)]
    out = annotate_sites(sites, rates)
    assert pd.api.types.is_float_dtype(out[CE]), out[CE].dtype
    got = out[CE].tolist()
    assert got[0] == 0.0 and got[1] == 1.0, got  # 0 and 1 stay numbers
    assert math.isnan(got[2]), got
    assert out["class_error_n"].isna().tolist() == [False, False, True]


def _orf_fdr(q: float, reason: str | None) -> Any:
    from gippa.novel.class_call import OrfClassFdr

    return OrfClassFdr("NOVEL_c1_1_p", "picked_decoy", "intergenic", q, reason, False)


def test_class_fdr_cell_is_a_number_or_nan() -> None:
    from gippa.novel.class_call import class_fdr_cell

    assert class_fdr_cell(_orf_fdr(0.0, None)) == 0.0
    assert class_fdr_cell(_orf_fdr(0.02, None)) == 0.02
    assert class_fdr_cell(_orf_fdr(math.nan, "lost_pick")) == 1.0  # Q1 of M7-11b, unchanged
    cell = class_fdr_cell(_orf_fdr(math.nan, "no_decoys"))
    assert isinstance(cell, float) and math.isnan(cell), cell


def test_minimal_call_class_fdr_is_nan() -> None:
    from helpers.novel_outputs_fixture import cfg

    from gippa.novel.candidates import minimal_call

    df = minimal_call(orfs_of(two_contigs()), cfg())
    assert len(df) > 0
    assert pd.api.types.is_float_dtype(df[FDR]), df[FDR].dtype
    assert df[FDR].isna().all()


# ---------------------------------------------------------------- variants writer + VCF


@cache
def _bag() -> tuple[Any, pd.DataFrame]:
    """``(pile, sites)`` with four sites on both strands (two ALT substitutions, a block)."""
    aligns = (
        peptides_at(GP, 44, sub(GP, 44, "S"), 2)
        + peptides_at(GP, 50, ("block", (GP.codon(50),), "N", "GG"), 1)
        + peptides_at(GM_, 10, sub(GM_, 10, "S"), 1)
        + peptides_at(GP, 80, sub(GP, 80, "F"), 1)
    )
    return run(aligns)


def _vcf_table(rates: Sequence[object], flags: Sequence[tuple[str, ...]] | None = None
               ) -> pd.DataFrame:  # fmt: skip
    pile, sites = _bag()
    st_ = sites.copy()
    st_[CE] = pd.Series(list(rates), index=st_.index, dtype=object)
    if flags is not None:
        st_["flags"] = list(flags)
    _, frame = vr().vcf_records(st_, pile, GENOME, REPS)
    return frame


def _tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def test_k_var25_example_flags_and_na(tmp_path: Path) -> None:
    """Spec example: an unestimated K site gets NA and the flag after its other flags."""
    _, frame, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3),
                             {"flags": [K_VAR25_FLAGS], CE: [math.nan]})  # fmt: skip
    (v,) = frame[CE].tolist()
    assert _na(v), v
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, frame)
    (row,) = _tsv(path)
    assert row[CE] == "NA", row
    assert row["flags"] == K_VAR25_1_1, row


def test_missing_class_error_column_means_na_and_flag(tmp_path: Path) -> None:
    _, frame, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3))
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, frame)
    (row,) = _tsv(path)
    assert row[CE] == "NA" and row["flags"].split(";")[-1] == CE_FLAG, row


@pytest.mark.parametrize(("rate", "cell"), [(0.0, "0"), (1.0, "1"), (0.0123456789, "0.0123457")])
def test_zero_and_one_stay_numbers_without_the_flag(tmp_path: Path, rate: float, cell: str) -> None:
    _, frame, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 3),
                             {"flags": [("snv", CE_FLAG)], CE: [rate]})  # fmt: skip
    (v,) = frame[CE].tolist()
    assert isinstance(v, float) and v == rate, v
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, frame)
    (row,) = _tsv(path)
    assert row[CE] == cell and row["flags"] == "snv", row  # a stale flag is dropped


_RATE = st.one_of(st.none(), st.just(math.nan), st.just(0.0), st.just(1.0),
                  st.floats(0, 1, allow_nan=False))  # fmt: skip


@settings(max_examples=40, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(data=st.data())
def test_variants_na_iff_flag_through_the_writer(tmp_path_factory: Any, data: Any) -> None:
    """AC 11c.2 on random rows: frame and file NA exactly where the flag is, numbers kept."""
    _, sites = _bag()
    n = len(sites)
    rates = data.draw(st.lists(_RATE, min_size=n, max_size=n))
    stale = data.draw(st.lists(st.booleans(), min_size=n, max_size=n))
    flags = [("snv", CE_FLAG) if s else ("snv",) for s in stale]
    frame = _vcf_table(rates, flags)
    for v, fl in zip(frame[CE], frame["flags"], strict=True):
        assert _is_number_or_na(v), v
        assert _na(v) == (CE_FLAG in _flag_set(fl)), (v, fl)
    path = Path(tmp_path_factory.mktemp("v")) / "variants.tsv"
    vt().write_variants(path, frame)
    rows = _tsv(path)
    for r, want in zip(rows, rates, strict=True):
        assert (r[CE] == "NA") == (CE_FLAG in r["flags"].split(";")), r
        if _na(want):
            assert r[CE] == "NA", r
        else:
            assert float(r[CE]) == float(f"{want:.6g}"), (r[CE], want)
    back = vt().read_variants(path)
    assert pd.api.types.is_float_dtype(back[CE]), back[CE].dtype
    assert back[CE].isna().tolist() == [r[CE] == "NA" for r in rows]


def test_written_variants_with_na_validate_against_the_schema(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    frame = _vcf_table([math.nan, 0.0, 1.0, None])
    path = tmp_path / "variants.tsv"
    vt().write_variants(path, frame)
    d = {k: v for k, v in _schema("variants").items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(tmp_path), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_vcf_class_err_is_dot_for_na_and_a_number_otherwise(tmp_path: Path) -> None:
    """The VCF carries ``CLASS_ERR`` (Number=A); an NA allele reads ``.``, never a word."""
    assert any(f.id == "CLASS_ERR" for f in vr().INFO_FIELDS)
    recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), {CE: [math.nan]})
    write(tmp_path / "na.vcf", recs)
    assert info_dict(body_lines(tmp_path / "na.vcf")[0][7])["CLASS_ERR"] == "."
    for rate, cell in ((0.0, "0"), (1.0, "1"), (0.0123456789, "0.0123457")):
        recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), {CE: [rate]})
        write(tmp_path / "v.vcf", recs)
        assert info_dict(body_lines(tmp_path / "v.vcf")[0][7])["CLASS_ERR"] == cell


def test_vcf_never_writes_the_word_unestimated(tmp_path: Path) -> None:
    recs, *_ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1))  # no rate column at all
    write(tmp_path / "v.vcf", recs)
    body = "\n".join("\t".join(r) for r in body_lines(tmp_path / "v.vcf"))
    assert "unestimated" not in body.replace(CE_FLAG, ""), body


# ---------------------------------------------------------------- novel_orfs writer


_FDR = st.one_of(st.none(), st.just(math.nan), st.just(0.0), st.just(1.0),
                 st.floats(0, 1, allow_nan=False))  # fmt: skip


@settings(max_examples=30, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(data=st.data())
def test_novel_orfs_na_iff_flag_through_the_writer(tmp_path_factory: Any, data: Any) -> None:
    case = two_contigs()
    orfs = orfs_of(case)
    ids = [o.orf_id for o in orfs]
    fdrs = data.draw(st.lists(_FDR, min_size=len(ids), max_size=len(ids)))
    t = table(case, orfs, scored(orfs, class_fdr=dict(zip(ids, fdrs, strict=True))))
    assert pd.api.types.is_float_dtype(t[FDR]), t[FDR].dtype
    for v, fl in zip(t[FDR], t["flags"], strict=True):
        assert _na(v) == (FDR_FLAG in _flag_set(fl)), (v, fl)
    path = Path(tmp_path_factory.mktemp("n")) / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, t)
    rows = _tsv(path)
    for r in rows:
        assert (r[FDR] == "NA") == (FDR_FLAG in r["flags"].split(";")), r
        assert r[FDR] == "NA" or re.fullmatch(r"[0-9.eE+-]+", r[FDR]), r
    back = novel_io().read_novel_orfs(path)
    assert pd.api.types.is_float_dtype(back[FDR]), back[FDR].dtype
    again = path.with_name("again.tsv")
    novel_io().write_novel_orfs(again, back)
    assert again.read_bytes() == path.read_bytes()


def test_novel_orfs_zero_stays_zero(tmp_path: Path) -> None:
    """Spec example C: ``0`` -> ``0``, no flag."""
    case = two_contigs()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs, class_fdr={o.orf_id: 0.0 for o in orfs}))
    path = tmp_path / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, t)
    for r in _tsv(path):
        assert r[FDR] == "0" and FDR_FLAG not in r["flags"].split(";"), r


def test_written_novel_orfs_with_na_validate_against_the_schema(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    case = two_contigs()
    orfs = orfs_of(case)
    t = table(case, orfs, scored(orfs))  # minimal_call only: every class FDR unestimated
    path = tmp_path / "novel_orfs.tsv"
    novel_io().write_novel_orfs(path, t)
    assert {r[FDR] for r in _tsv(path)} == {"NA"}
    d = {k: v for k, v in _schema("novel_orfs").items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(tmp_path), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


# ---------------------------------------------------------------- fixtures and goldens


def _read_fixture(rel: str, stem: str) -> pd.DataFrame:
    from gippa.io.run_view import read_schema_tsv

    return read_schema_tsv(REPO / rel, stem, _schema(stem))


@pytest.mark.parametrize(("rel", "stem"), FIXTURES, ids=[r for r, _ in FIXTURES])
def test_every_fixture_reads_as_float_and_na_iff_flag(rel: str, stem: str) -> None:
    """AC 11c.1 (pandas reads a float column through the schema) and AC 11c.2."""
    frame = _read_fixture(rel, stem)
    col, flag = COLUMN[stem]
    assert pd.api.types.is_float_dtype(frame[col]), (rel, frame[col].dtype)
    for v, fl in zip(frame[col], frame["flags"], strict=True):
        assert _na(v) == (flag in _flag_set(fl)), (rel, v, fl)


@pytest.mark.parametrize("rel", sorted(GOLDEN_SHA256))
def test_numeric_goldens_are_unchanged(rel: str) -> None:
    got = hashlib.sha256((REPO / rel).read_bytes()).hexdigest()
    assert got == GOLDEN_SHA256[rel], f"{rel} changed; it holds numbers only (spec PR 11c example)"


# ---------------------------------------------------------------- 1.0 output folders (Q3)


def _write_1_0(out: Path, stem: str, rows: list[dict[str, str]]) -> None:
    """A 1.0 table: the 1.0 header (11b renames put back), cells as given, the rest by type."""
    fill = {"integer": "1", "number": "0.5", "boolean": "true"}
    names = [f["name"] for f in _schema(stem)["fields"]]
    kinds = {f["name"]: str(f.get("type")) for f in _schema(stem)["fields"]}
    full = [{c: r.get(c, fill.get(kinds[c], "x")) for c in names} for r in rows]
    path = out / "tables" / f"{stem}.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["\t".join(to_1_0(stem, names))] + ["\t".join(r[c] for c in names) for r in full]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _variant_row(i: int, rate: str, flags: str) -> dict[str, str]:
    return {"variant_id": f"VAR_{i}", "contig": "chr1", "codon_start": str(100 + 3 * i),
            "strand": "+", "gene": "g1", "aa_pos": str(10 + i), "ref_codon": "GCT",
            "ref_aa": "A", "alt_aa": "V", "alt_type": "substitution", "call": "ALT",
            "confidence": "high", "n_ref": "0", "n_alt": "3", "nt_changes": "1",
            "alt_codon_vcf": "GTT", CE: rate, "class_error_n": "NA", "source": "denovo",
            "flags": flags}  # fmt: skip


def test_runview_maps_a_1_0_variants_table(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    out = tmp_path / "OUT"
    _write_1_0(out, "variants", [
        _variant_row(1, "unestimated", ";".join(K_VAR25_FLAGS)),
        _variant_row(2, "unestimated", "NA"),
        _variant_row(3, "0", "snv"),
        _variant_row(4, "1", "NA"),
        _variant_row(5, "0.0123", "snv"),
    ])  # fmt: skip
    frame = RunView.from_outdir(out).table("variants")
    assert frame is not None
    assert pd.api.types.is_float_dtype(frame[CE]), frame[CE].dtype
    rates = frame[CE].tolist()
    assert math.isnan(rates[0]) and math.isnan(rates[1]), rates
    assert rates[2:] == [0.0, 1.0, 0.0123], rates
    flags = frame["flags"].tolist()
    assert flags[0] == K_VAR25_1_1, flags
    assert flags[1] == CE_FLAG, flags
    assert flags[2] == "snv" and _na(flags[3]) and flags[4] == "snv", flags


def _orf_row(i: int, fdr: str, flags: str) -> dict[str, str]:
    return {"orf_id": f"NOVEL_c1_{i}_p", "class": "intergenic", "contig": "c1", "strand": "+",
            "called": "false", FDR: fdr, "flags": flags, "start_codon_best": "ATG",
            "peptides": "PEPTIDE", "overlapping_gene": "NA"}  # fmt: skip


def test_runview_maps_a_1_0_novel_orfs_table(tmp_path: Path) -> None:
    from gippa.io.run_view import RunView

    out = tmp_path / "OUT"
    _write_1_0(out, "novel_orfs", [
        _orf_row(1, "unestimated", f"{FDR_FLAG};no_start_codon"),
        _orf_row(2, "unestimated", "NA"),
        _orf_row(3, "0", "rbs_supported"),
        _orf_row(4, "1", "NA"),
    ])  # fmt: skip
    frame = RunView.from_outdir(out).table("novel_orfs")
    assert frame is not None
    assert pd.api.types.is_float_dtype(frame[FDR]), frame[FDR].dtype
    fdrs = frame[FDR].tolist()
    assert math.isnan(fdrs[0]) and math.isnan(fdrs[1]) and fdrs[2:] == [0.0, 1.0], fdrs
    flags = frame["flags"].tolist()
    assert flags[0].split(";").count(FDR_FLAG) == 1, flags  # present: not added twice
    assert flags[1] == FDR_FLAG, flags
    assert flags[2] == "rbs_supported" and _na(flags[3]), flags


# ---------------------------------------------------------------- migration note (§ 5.3, § 5.4)


def test_changelog_removal_fragment_states_the_type_change() -> None:
    from helpers.removal_notes import removal_texts

    text = "\n".join(removal_texts())
    for token in (f"`{CE}`", f"`{FDR}`", f"`{CE_FLAG}`", f"`{FDR_FLAG}`"):
        assert token in text, f"{token} not in a changelog.d/*.removal.md fragment"


def test_outputs_changes_in_1_1_lists_the_type_and_flag_changes() -> None:
    text = (REPO / "docs" / "user-guide" / "outputs.md").read_text(encoding="utf-8")
    sec = text[text.index("## Changes in 1.1") :]
    for col, flag in ((CE, CE_FLAG), (FDR, FDR_FLAG)):
        assert any(f"`{col}`" in ln and f"`{flag}`" in ln for ln in sec.splitlines()), col
    assert any(f"`{CE_FLAG}`" in ln and "new" in ln.lower() for ln in sec.splitlines())


def test_schema_exceptions_document_the_retyped_columns() -> None:
    doc = json.loads((REPO / "scripts" / "release" / "schema_exceptions.json").read_text("utf-8"))
    retyped = doc["v1.0.0"]["retyped"]
    assert CE in retyped.get("variants", {}), retyped
    assert FDR in retyped.get("novel_orfs", {}), retyped
