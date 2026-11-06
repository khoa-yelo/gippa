"""M12-04 PR 11c, review round 1: the VCF ``CLASS_ERR`` header, a strict ``class_error_rate``
cell, and the docs that describe the NA cells (lead decisions M1, L1-L4 on #558).

Interface pinned here (test-author reading of the lead decisions):

- M1: ``gippa.variants.vcf_records.INFO_FIELDS`` declares ``CLASS_ERR`` as ``Number=A``,
  ``Type=Float`` with the description ``Class-specific error rate of each allele; . when not
  estimated``; the Nextflow stub header in ``workflow/subworkflows/local/variants.nf`` writes the
  same ``##INFO``/``##FORMAT`` lines as the Python writer; the golden
  ``tests/data/golden/sim50k_m7/variants.vcf`` carries the new line. htslib (pysam) reads a
  number as a float and ``.`` as missing.
- L1: ``vcf_records`` raises ``ValueError`` for a ``class_error_rate`` cell that is neither a
  number nor NA (the word ``unestimated``, ``abc``), as ``gippa.novel.outputs.orf_table`` does.
- L2: the methods/user-guide pages never describe these cells as holding the word
  ``unestimated`` (the flag names stay allowed).
- L3: the changelog fragment and outputs.md > Changes in 1.1 name the VCF ``FLAGS`` value
  ``class_error_unestimated`` and ``CLASS_ERR`` ``Type=Float`` with ``.``.
- L4: the generated ``docs/reference/tables.md`` intro says "renamed or retyped".
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest
from helpers.variant_tables_fixture import (
    GP,
    body_lines,
    info_dict,
    peptides_at,
    records,
    sub,
    write,
)

REPO = Path(__file__).resolve().parents[2]
CE, CE_FLAG, FDR_FLAG = "class_error_rate", "class_error_unestimated", "class_fdr_unestimated"
CLASS_ERR_LINE = (
    "##INFO=<ID=CLASS_ERR,Number=A,Type=Float,"
    'Description="Class-specific error rate of each allele; . when not estimated">'
)
STUB = REPO / "workflow" / "subworkflows" / "local" / "variants.nf"
GOLDEN_VCF = REPO / "tests" / "data" / "golden" / "sim50k_m7" / "variants.vcf"


def _site(rate: object) -> list[object]:
    return list(records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), {CE: [rate]})[0])


def _meta_lines(path: Path) -> list[str]:
    return [ln for ln in path.read_text(encoding="utf-8").splitlines()
            if ln.startswith(("##INFO=", "##FORMAT="))]  # fmt: skip


def _python_meta(tmp_path: Path) -> list[str]:
    """The ``##INFO`` and ``##FORMAT`` lines the Python writer emits, in order."""
    path = tmp_path / "empty.vcf"
    write(path, [])
    return _meta_lines(path)


def _stub_meta() -> list[str]:
    """The stub's ``##INFO``/``##FORMAT`` lines, built as its ``vcfMeta`` builds them."""
    text = STUB.read_text(encoding="utf-8")
    out = []
    for kind, var in (("INFO", "info"), ("FORMAT", "fmt")):
        m = re.search(rf"def {var} = (\[\[.*?\]\])\n", text, re.S)
        assert m, f"variants.nf stub has no `def {var}` list"
        for f in re.findall(r"\['(\w+)', '([^']*)', '(\w+)', '([^']*)'\]", m.group(1)):
            out.append(f'##{kind}=<ID={f[0]},Number={f[1]},Type={f[2]},Description="{f[3]}">')
    return out


# ---------------------------------------------------------------- M1: the CLASS_ERR header


def test_python_class_err_header_line_is_exact(tmp_path: Path) -> None:
    lines = [ln for ln in _python_meta(tmp_path) if ln.startswith("##INFO=<ID=CLASS_ERR,")]
    assert lines == [CLASS_ERR_LINE]


def test_python_info_field_class_err_is_float() -> None:
    from gippa.variants.vcf_records import INFO_FIELDS

    (f,) = [f for f in INFO_FIELDS if f.id == "CLASS_ERR"]
    assert (f.number, f.type) == ("A", "Float"), f


def test_nf_stub_class_err_line_is_exact() -> None:
    lines = [ln for ln in _stub_meta() if ln.startswith("##INFO=<ID=CLASS_ERR,")]
    assert lines == [CLASS_ERR_LINE]


def test_nf_stub_meta_lines_equal_the_python_writer(tmp_path: Path) -> None:
    """Full ``##INFO``/``##FORMAT`` lines, not just the IDs: a stub run publishes the header a
    real run writes."""
    stub = _stub_meta()
    assert stub == _python_meta(tmp_path)
    assert sum(ln.startswith("##INFO=") for ln in stub) >= 10, stub


def test_golden_vcf_has_the_float_class_err_header() -> None:
    meta = _meta_lines(GOLDEN_VCF)
    assert [ln for ln in meta if ln.startswith("##INFO=<ID=CLASS_ERR,")] == [CLASS_ERR_LINE]
    for row in body_lines(GOLDEN_VCF):
        for v in info_dict(row[7]).get("CLASS_ERR", ".").split(","):
            assert v == "." or math.isfinite(float(v)), row


def test_htslib_reads_class_err_as_float_and_dot_as_missing(
    tmp_path: Path, capfd: pytest.CaptureFixture[str]
) -> None:
    pysam = pytest.importorskip("pysam")
    got = {}
    for name, rate in (("num", 0.0123456789), ("zero", 0.0), ("na", math.nan)):
        path = tmp_path / f"{name}.vcf"
        write(path, _site(rate))
        capfd.readouterr()
        with pysam.VariantFile(str(path)) as vf:
            assert vf.header.info["CLASS_ERR"].type == "Float"
            assert vf.header.info["CLASS_ERR"].number == "A"
            (rec,) = list(vf)
            got[name] = rec.info["CLASS_ERR"]
        err = capfd.readouterr().err
        assert err == "", err
    (num,) = got["num"]
    assert isinstance(num, float) and num == pytest.approx(0.0123457, rel=1e-6), got
    (zero,) = got["zero"]
    assert isinstance(zero, float) and zero == 0.0, got
    assert got["na"] == (None,), got


# ---------------------------------------------------------------- L1: strict cells


@pytest.mark.parametrize("cell", ["unestimated", "abc"])
def test_vcf_records_rejects_a_word_in_class_error_rate(cell: str) -> None:
    with pytest.raises(ValueError, match=re.escape(cell)):
        records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), {CE: [cell]})


@pytest.mark.parametrize(("cell", "want"), [("0.25", 0.25), ("0", 0.0), (0.5, 0.5)])
def test_vcf_records_reads_a_numeric_cell(cell: object, want: float) -> None:
    _, frame, _, _ = records(peptides_at(GP, 44, sub(GP, 44, "S"), 1), {CE: [cell]})
    (v,) = frame[CE].tolist()
    assert isinstance(v, float) and v == want, v


# ---------------------------------------------------------------- L3: changelog and outputs.md


def _sentences(text: str) -> list[str]:
    flat = " ".join(text.split())
    return re.split(r"(?<=[.;])\s+(?=[A-Z`(])", flat)


def _assert_vcf_changes(text: str, where: str) -> None:
    sents = _sentences(text)
    assert any("`FLAGS`" in s and f"`{CE_FLAG}`" in s for s in sents), (
        f"{where}: no sentence says the VCF `FLAGS` field carries `{CE_FLAG}`"
    )
    assert any("`CLASS_ERR`" in s and "Type=Float" in s for s in sents), (
        f"{where}: no sentence says `CLASS_ERR` is now Type=Float"
    )
    assert any("`CLASS_ERR`" in s and "`.`" in s for s in sents), (
        f"{where}: no sentence says an unestimated `CLASS_ERR` is `.`"
    )


def test_changelog_fragment_names_the_vcf_changes() -> None:
    text = "\n".join(p.read_text(encoding="utf-8")
                     for p in sorted((REPO / "changelog.d").glob("*.removal.md"))
                     if CE in p.read_text(encoding="utf-8"))  # fmt: skip
    assert text, "no changelog.d/*.removal.md fragment for class_error_rate"
    _assert_vcf_changes(text, "changelog fragment")


def test_outputs_changes_in_1_1_names_the_vcf_changes() -> None:
    text = (REPO / "docs" / "user-guide" / "outputs.md").read_text(encoding="utf-8")
    sec = text[text.index("## Changes in 1.1") :]
    nxt = re.search(r"\n## ", sec[3:])
    _assert_vcf_changes(sec[: nxt.start() + 3] if nxt else sec, "outputs.md > Changes in 1.1")


# ---------------------------------------------------------------- L2: no word in these cells

L2_FILES = (
    "docs/methods/novel-orfs.md",
    "docs/methods/variants.md",
    "docs/user-guide/interpreting-results.md",
    "docs/get-started/from-spectra.md",
)
_COLUMN = re.compile(
    r"class_fdr(?!_)|class_error_rate|CLASS_ERR|class FDR|class error rate", re.IGNORECASE
)
_WORD = re.compile(r"`unestimated`")  # the bare word; `class_fdr_unestimated` does not match
_BAD = (
    re.compile(r"(class_fdr|class_error_rate|CLASS_ERR)`?\s*(=|==|is)\s*`?unestimated(?![\w-])"),
    re.compile(r"\bgives?\s+`?unestimated(?![\w-])"),
)


def _chunks(text: str) -> list[str]:
    """Sentences of each paragraph, and each table row on its own."""
    out = []
    for block in re.split(r"\n\s*\n", text):
        parts = block.splitlines() if block.lstrip().startswith("|") else [block]
        for p in parts:
            out += _sentences(p)
    return out


@pytest.mark.parametrize("rel", L2_FILES)
def test_docs_do_not_put_the_word_unestimated_in_these_cells(rel: str) -> None:
    text = (REPO / rel).read_text(encoding="utf-8")
    bad = [
        s
        for s in _chunks(text)
        if any(b.search(s) for b in _BAD) or (_COLUMN.search(s) and _WORD.search(s))
    ]
    assert not bad, f"{rel} describes class_fdr/class_error_rate as the word unestimated: {bad}"


def test_l2_check_allows_the_flag_names() -> None:
    ok = (
        f"An unestimated class FDR is NA and the row carries `{FDR_FLAG}`. "
        f"The class error rate is NA with the flag `{CE_FLAG}`."
    )
    assert not any(
        any(b.search(s) for b in _BAD) or (_COLUMN.search(s) and _WORD.search(s))
        for s in _chunks(ok)
    )
    for s in ("and `class_fdr = unestimated`.", "a database-only input give `unestimated`.",
              "the class FDR is `unestimated`, and"):  # fmt: skip
        assert any(b.search(s) for b in _BAD) or (_COLUMN.search(s) and _WORD.search(s)), s


# ---------------------------------------------------------------- L4: tables.md intro


def test_tables_reference_intro_says_renamed_or_retyped() -> None:
    page = (REPO / "docs" / "reference" / "tables.md").read_text(encoding="utf-8")
    intro = page[: page.index("\n## ")]
    assert "Columns renamed or retyped in 1.1 are listed in" in intro, intro
    gen = (REPO / "scripts" / "docs" / "gen_reference.py").read_text(encoding="utf-8")
    assert "Columns renamed or retyped in 1.1" in gen
