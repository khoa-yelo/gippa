"""M12-04 PR 11b: one name per quantity (T-04): schemas, writers, readers, goldens, docs.

docs/specs/M12-04-outputs-revision.md > PR 11b, § 5.1 and AC 11b.1, 11b.4; ADR-0075. Lead
defaults: Q1 ``pair_class`` in both tables, Q2 ``operon_flags.flag`` -> ``flags``, Q15 included.
No aliases: an old name is never written next to its new one.

Test-author readings (beyond the spec text):

- AC 11b.1 lists ``class`` among the forbidden field names, but ``novel_orfs.class`` and
  ``protein_hits.class`` are not in the § 5.1 rename table. The schema check is per table: the
  old names of § 5.1 are gone from *their* tables (``class`` from ``pairs``, ``flag`` from
  ``operon_flags``, ``mw`` from ``determinants_features``); the other tables are untouched.
- AC 11b.4's grep covers the distinctive old names (``GREP_OLD``) as quoted column keys or
  attribute access in ``src/gippa/**/*.py`` and as words in ``workflow/**/*.nf`` headers. The
  config key ``input.peptide_len`` (attribute access ``cfg.input.peptide_len``, the
  ``peptide_len:`` pydantic field) and the functions ``promoter_between`` and
  ``terminator_between`` stay; the allow-list is the migration module ``gippa/io/migration.py``.
  The generic old names (``class``, ``flag``, ``mw``) are covered through the writers' column
  constants and frames instead.
- ``OUTPUT_SCHEMA_VERSION`` goes to at least ``1.3`` with major 1 (Q14); a later 1.1 PR may bump
  it again.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd
import pytest
from helpers.m12_04_renames import GREP_OLD, RENAMES, TRIPLES, VALUE_RENAMES

REPO = Path(__file__).resolve().parents[2]
SCHEMAS = REPO / "src" / "gippa" / "schemas"
TABLES = SCHEMAS / "tables"
DATA = REPO / "tests" / "data"
MIGRATION_MODULE = REPO / "src" / "gippa" / "io" / "migration.py"


def _fields(stem: str) -> list[str]:
    doc = json.loads((TABLES / f"{stem}.schema.json").read_text(encoding="utf-8"))
    return [f["name"] for f in doc["fields"]]


# ---------------------------------------------------------------- AC 11b.1: schemas


@pytest.mark.parametrize(("table", "old", "new"), TRIPLES, ids=[f"{t}.{o}" for t, o, _ in TRIPLES])
def test_schema_has_the_new_name_and_not_the_old(table: str, old: str, new: str) -> None:
    names = _fields(table)
    assert new in names, f"{table}.schema.json lacks {new!r}: {names}"
    assert old not in names, f"{table}.schema.json still has the 1.0 name {old!r} (no aliases)"


@pytest.mark.parametrize("table", sorted(RENAMES))
def test_schema_columns_are_unique_and_keep_their_positions(table: str) -> None:
    """A rename keeps the column where it was: the 1.1 field list is the 1.0 list renamed."""
    names = _fields(table)
    assert len(names) == len(set(names)), names
    back = {new: old for old, new in RENAMES[table].items()}
    restored = [back.get(n, n) for n in names]
    assert not set(RENAMES[table].values()) & set(restored), restored


def test_operon_flags_and_pairs_share_pair_class_and_flags() -> None:
    """Q1, Q2: one name per pair quantity in both tables (a join on pair_id needs suffixes)."""
    for table in ("pairs", "operon_flags"):
        names = _fields(table)
        assert "pair_class" in names and "flags" in names, (table, names)


def test_unchanged_on_purpose() -> None:
    """§ 5.1: the yes/no calls and pairs.tx_rate_between keep their names."""
    assert {"promoter_between_yes", "terminator_between_yes"} <= set(_fields("operon_flags"))
    assert "tx_rate_between" in _fields("pairs")


def test_pair_class_keeps_the_class_enum() -> None:
    def enum(table: str, col: str) -> list[str]:
        doc = json.loads((TABLES / f"{table}.schema.json").read_text(encoding="utf-8"))
        (f,) = [f for f in doc["fields"] if f["name"] == col]
        return list(f.get("constraints", {}).get("enum", []))

    assert set(enum("pairs", "pair_class")) == {"OP", "NO", "CV", "DV", "SW"}
    assert set(enum("operon_flags", "pair_class")) == {"OP", "NO", "CV", "DV", "SW"}


def test_determinants_pd_feature_enum_uses_the_new_value() -> None:
    doc = json.loads((SCHEMAS / "plot_support" / "determinants_pd.json").read_text("utf-8"))
    (feat,) = [f for f in doc["fields"] if f["name"] == "feature"]
    values = feat["constraints"]["enum"]
    assert "promoter_score" in values and "promotech_max" not in values, values


def test_output_schema_version_is_at_least_1_3_with_major_1() -> None:
    """Q14: 1.2 -> 1.3; the major stays 1 under ADR-0075."""
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    major, minor = (int(x) for x in OUTPUT_SCHEMA_VERSION.split("."))
    assert major == 1 and minor >= 3, OUTPUT_SCHEMA_VERSION
    package = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    assert package["output_schema_version"] == OUTPUT_SCHEMA_VERSION


# ---------------------------------------------------------------- writers


def test_writer_column_constants_equal_the_schemas() -> None:
    from gippa.determinants.features import FEATURE_COLUMNS
    from gippa.mapping.outputs import PEPTIDE_HITS_COLUMNS
    from gippa.novel.outputs import NOVEL_ORFS_COLUMNS
    from gippa.operons.flags import FLAG_COLUMNS
    from gippa.operons.pairs import PAIR_COLUMNS

    for table, cols in (
        ("novel_orfs", NOVEL_ORFS_COLUMNS),
        ("determinants_features", FEATURE_COLUMNS),
        ("peptide_hits", PEPTIDE_HITS_COLUMNS),
        ("pairs", PAIR_COLUMNS),
        ("operon_flags", FLAG_COLUMNS),
    ):
        assert list(cols) == _fields(table), table
        assert not set(RENAMES[table]) & set(cols), (table, cols)


def test_genetable_contract_peptide_hits_uses_peptide_len_aa() -> None:
    from gippa.genetable.contract import CONTRACT_COLUMNS

    cols = CONTRACT_COLUMNS["peptide_hits"]
    assert "peptide_len_aa" in cols and "peptide_len" not in cols, cols


def test_mapping_stub_header_uses_peptide_len_aa() -> None:
    """The loci stub of workflow/subworkflows/local/mapping.nf writes the §9.3 header."""
    text = (REPO / "workflow" / "subworkflows" / "local" / "mapping.nf").read_text("utf-8")
    m = re.search(r"peptide_hits:\s*'([^']*)'", text)
    assert m is not None
    assert m.group(1).split() == _fields("peptide_hits")


def test_pair_frame_writes_pair_class() -> None:
    from helpers.pairs_data import BASE, frame, genes_of, ring

    probs = {("G1", "G2"): 0.8, ("G3", "G4"): 0.9, ("G5", "G1"): 0.1}  # test_pairs P_BASE
    df = frame(genes_of(BASE), ring(), probs)
    assert list(df.columns[: len(_fields("pairs"))]) == _fields("pairs")
    assert "class" not in df.columns
    assert set(df["pair_class"]) <= {"OP", "NO", "CV", "DV", "SW"}
    assert "OP" in set(df["pair_class"])


def test_flag_pairs_reads_and_writes_the_new_names() -> None:
    """``flag_pairs`` takes the 1.1 pair frame (``pair_class``, ``promoter_score_between``,
    ``terminator_score_between``) and writes ``pair_class`` and ``flags``; spec example 1."""
    from helpers.flags_data import cfg, frame, row

    from gippa.operons.flags import flag_pairs

    pairs = frame(row("OP", prom=0.81, term=0.10, delta=1.4, q=0.02)).rename(
        columns=dict(RENAMES["pairs"])
    )
    out = flag_pairs(pairs, cfg())
    assert list(out.columns) == _fields("operon_flags")
    assert list(out["pair_class"]) == ["OP"]
    assert list(out["flags"]) == ["internal_promoter_supported"]


def test_write_operon_tables_headers(tmp_path: Path) -> None:
    from gippa.io.operon_tables import write_operon_flags_tsv, write_pairs_tsv
    from gippa.operons.flags import FLAG_COLUMNS
    from gippa.operons.pairs import PAIR_COLUMNS

    write_pairs_tsv(tmp_path / "pairs.tsv", pd.DataFrame(columns=list(PAIR_COLUMNS)))
    write_operon_flags_tsv(tmp_path / "flags.tsv", pd.DataFrame(columns=list(FLAG_COLUMNS)))
    for name, table in (("pairs.tsv", "pairs"), ("flags.tsv", "operon_flags")):
        header = (tmp_path / name).read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert header == _fields(table), name


def test_determinants_features_take_the_genes_values() -> None:
    """``promoter_score``/``promoter_score_tu``/``mw_da`` hold exactly the genes.tsv columns of
    the same name (the 1.0 table already copied them)."""
    from gippa.determinants.features import BLOCKS, FEATURE_BLOCKS, GENE_TABLE_SOURCES

    for name in VALUE_RENAMES.values():
        assert GENE_TABLE_SOURCES[name] == name
    for old in VALUE_RENAMES:
        assert old not in GENE_TABLE_SOURCES
        assert all(old not in fs for fs in FEATURE_BLOCKS.values())
        assert all(old not in cs for cs in BLOCKS.values())  # feature_or_block values
    assert FEATURE_BLOCKS["transcription"][:2] == ("promoter_score", "promoter_score_tu")
    assert "mw_da" in FEATURE_BLOCKS["protein"]


def test_pd_features_and_v15_label() -> None:
    """V15(b): "promoter score (promotech_max)" becomes "promoter score"."""
    from gippa.determinants import support
    from gippa.viz import determinants as v15

    want = ("cai", "tai", "tir_ostir", "promoter_score")
    assert want == support.PD_FEATURES
    assert want == v15.PD_FEATURES
    assert v15.PD_XLABEL["promoter_score"] == "promoter score"
    assert v15.PD_LAYER["promoter_score"] == "promoters"
    for mapping in (v15.PD_XLABEL, v15.PD_LAYER):
        assert "promotech_max" not in mapping
    assert all("promotech" not in str(v) for v in v15.PD_XLABEL.values())


# ---------------------------------------------------------------- AC 11b.4: source grep

_QUOTED = re.compile(r"""["']({})["']""".format("|".join(GREP_OLD)))
_ATTR = re.compile(r"\.({})\b".format("|".join(n for n in GREP_OLD if n != "peptide_len")))
_NF_WORD = re.compile(r"(?<![\w-])({})(?![\w-])".format("|".join(GREP_OLD)))


def _py_hits() -> list[str]:
    hits = []
    for path in sorted((REPO / "src" / "gippa").rglob("*.py")):
        if path == MIGRATION_MODULE:
            continue
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if _QUOTED.search(line) or _ATTR.search(line):
                hits.append(f"{path.relative_to(REPO)}:{i}: {line.strip()}")
    return hits


def test_no_src_reference_to_an_old_column_name() -> None:
    hits = _py_hits()
    assert hits == [], "1.0 column names still used as keys:\n" + "\n".join(hits)


def test_no_workflow_header_uses_an_old_column_name() -> None:
    hits = []
    for path in sorted((REPO / "workflow").rglob("*.nf")):
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if "'" in line and _NF_WORD.search(line):
                hits.append(f"{path.relative_to(REPO)}:{i}")
    assert hits == [], hits


def test_the_migration_module_exists_for_the_allow_list() -> None:
    """The grep's allow-list must point at a real file (else the test above guards nothing)."""
    assert MIGRATION_MODULE.is_file(), MIGRATION_MODULE


# ---------------------------------------------------------------- goldens and fixtures

_FIXTURES = [
    "golden/sim50k_m7/novel_orfs.tsv",
    "golden/sim50k/mapping/peptide_hits.tsv",
    *(
        f"snapshot/{s}/OUT/{rel}"
        for s in ("sim_plots", "sim_empty")
        for rel in (
            "tables/novel_orfs.tsv",
            "mapping/peptide_hits.tsv",
            "tables/pairs.tsv",
            "tables/operon_flags.tsv",
            "tables/determinants_features.tsv",
        )
    ),
]
_STEM = {
    "novel_orfs.tsv": "novel_orfs",
    "peptide_hits.tsv": "peptide_hits",
    "pairs.tsv": "pairs",
    "operon_flags.tsv": "operon_flags",
    "determinants_features.tsv": "determinants_features",
}


def _gzipped(path: Path) -> bool:
    """Some snapshot fixtures are gzip under a ``.tsv`` name (the M9 fixture packer)."""
    with path.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def _header(path: Path) -> list[str]:
    if _gzipped(path):
        with gzip.open(path, "rt", encoding="utf-8") as fh:
            return fh.readline().rstrip("\r\n").split("\t")
    with path.open(encoding="utf-8") as fh:
        return fh.readline().rstrip("\r\n").split("\t")


@pytest.mark.parametrize("rel", _FIXTURES)
def test_golden_and_fixture_headers_use_the_new_names(rel: str) -> None:
    path = DATA / rel
    assert path.is_file(), path
    table = _STEM[path.name]
    header = _header(path)
    assert not set(RENAMES[table]) & set(header), (rel, header)
    assert set(RENAMES[table].values()) <= set(header), (rel, header)
    assert header == _fields(table), rel


#: sim_empty holds none of the renamed features (its model has too few rows), so only sim_plots
_VALUE_FILES = [
    ("snapshot/sim_plots/OUT/tables/determinants_model.tsv", "feature_or_block"),
    ("snapshot/sim_plots/OUT/plots/support/determinants_pd.tsv", "feature"),
    ("snapshot/sim_plots/OUT/plots/support/determinants_shap.tsv.gz", "feature"),
]


@pytest.mark.parametrize(("rel", "col"), _VALUE_FILES)
def test_fixture_feature_values_use_the_new_names(rel: str, col: str) -> None:
    path = DATA / rel
    assert path.is_file(), path
    comp = "gzip" if _gzipped(path) else None
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, compression=comp)
    values = set(frame[col])
    assert not values & set(VALUE_RENAMES), (rel, sorted(values & set(VALUE_RENAMES)))


def test_the_sim_plots_fixture_still_has_the_renamed_features() -> None:
    """The value check above is not vacuous: the model fixture names the renamed features."""
    path = DATA / "snapshot/sim_plots/OUT/tables/determinants_model.tsv"
    values = set(pd.read_csv(path, sep="\t", dtype=str)["feature_or_block"])
    assert set(VALUE_RENAMES.values()) <= values, sorted(values)


# ---------------------------------------------------------------- docs


def _doc_sections(text: str) -> dict[str, list[str]]:
    """``## `table` `` section -> the field names of its column table."""
    out: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        m = re.match(r"^## `([^`]+)`", line)
        if m:
            current = m.group(1)
            out[current] = []
        elif line.startswith("## "):
            current = None
        elif current is not None:
            cell = re.match(r"^\| `([^`]+)` \|", line)
            if cell:
                out[current].append(cell.group(1))
    return out


def test_tables_reference_uses_the_new_names() -> None:
    sections = _doc_sections((REPO / "docs" / "reference" / "tables.md").read_text("utf-8"))
    for table, renames in RENAMES.items():
        names = sections[table]
        assert names == _fields(table), table
        assert not set(renames) & set(names), table


def test_value_renames_follow_the_column_renames() -> None:
    """§ 5.1: the feature values change exactly as the determinants_features columns do."""
    assert dict(VALUE_RENAMES) == dict(RENAMES["determinants_features"])
