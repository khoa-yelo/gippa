"""M12-04 PR 11a, T-03: ``tables/genes.tsv`` in the groups identity -> detection -> abundance ->
function -> context -> features -> flags (docs/specs/M12-04-outputs-revision.md § PR 11a, § 5.2;
ADR-0075 decision 2).

AC covered here:

- 11a.1: the header the writer emits equals the § 5.2 list exactly, and so does the schema
  ``fields`` order (and ``GENES_COLUMNS``, which the writer, the reader and the Parquet files use);
- the column set, names and types do not change (§ PR 11a "Behaviour");
- 11a.5: a written ``genes.tsv`` validates with frictionless against the schema, whose descriptor
  order is the file order;
- ``docs/reference/tables.md`` lists the genes fields in the § 5.2 order (it is generated from the
  schema, so ``gen_reference.py --check`` keeps it in step);
- Q14: ``OUTPUT_SCHEMA_VERSION`` and ``datapackage.json`` ``output_schema_version`` leave 1.2; the
  major stays 1 (the spec states the bump to 1.3 under PR 11b; whichever of 11a/11b lands first
  bumps it, so the test asks for >= 1.3 with major 1).

``ORDER_1_1`` is § 5.2 verbatim: (column, group, 1.0 position). ``ORDER_1_0`` is derived from the
1.0 positions and is the order on main before this PR.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCHEMAS = REPO / "src" / "gippa" / "schemas"
GENES_SCHEMA = SCHEMAS / "tables" / "genes.schema.json"
TABLES_MD = REPO / "docs" / "reference" / "tables.md"

sys.path.insert(0, str(REPO / "tests" / "helpers"))

#: § 5.2: (column, group, 1.0 position), in the 1.1 order.
ORDER_1_1_TABLE: tuple[tuple[str, str, int], ...] = (
    ("locus_tag", "identity", 1),
    ("source_id", "identity", 2),
    ("gene", "identity", 3),
    ("product", "identity", 4),
    ("feature_type", "identity", 5),
    ("hypothetical", "identity", 6),
    ("contig", "identity", 7),
    ("start", "identity", 8),
    ("end", "identity", 9),
    ("strand", "identity", 10),
    ("length_aa", "identity", 11),
    ("partial", "identity", 12),
    ("crosses_origin", "identity", 13),
    ("replicon_type", "identity", 14),
    ("detected", "detection", 48),
    ("detection_status", "detection", 50),
    ("detection_tier", "detection", 49),
    ("gene_q_value", "detection", 51),
    ("detected_2pep", "detection", 52),
    ("n_peptides_unique", "detection", 53),
    ("n_peptides_multi_identical", "detection", 54),
    ("n_peptides_shared", "detection", 55),
    ("n_peptides_tolerant", "detection", 56),
    ("coverage_frac", "detection", 57),
    ("best_peptide_q", "detection", 58),
    ("p_detect", "detection", 59),
    ("abundance_log2", "abundance", 61),
    ("abundance_raw", "abundance", 60),
    ("top3_log2", "abundance", 66),
    ("abundance_imputed", "abundance", 62),
    ("abundance_type", "abundance", 64),
    ("abundance_scale", "abundance", 65),
    ("imputation_policy", "abundance", 63),
    ("cog_id", "function", 67),
    ("cog_category", "function", 68),
    ("go_terms", "function", 69),
    ("kegg_ko", "function", 70),
    ("ec", "function", 71),
    ("annotation_source", "function", 72),
    ("tu_id", "context", 15),
    ("tu_position", "context", 16),
    ("tu_size", "context", 17),
    ("uniop_p_prev", "context", 18),
    ("uniop_p_next", "context", 19),
    ("region_ids", "context", 20),
    ("island_id", "context", 21),
    ("replichore", "context", 22),
    ("replichore_position", "context", 23),
    ("gc_local_dev", "context", 24),
    ("promoter_score", "features", 25),
    ("promoter_score_tu", "features", 26),
    ("tx_rate_max_tu", "features", 27),
    ("tx_rate_sum_tu", "features", 28),
    ("terminator_conf_downstream", "features", 29),
    ("tir_ostir", "features", 30),
    ("sd_dg", "features", 31),
    ("mrna_dg", "features", 32),
    ("spacer_bp", "features", 33),
    ("start_codon", "features", 34),
    ("cai", "features", 35),
    ("tai", "features", 36),
    ("enc", "features", 37),
    ("fop", "features", 38),
    ("gc3", "features", 39),
    ("cpb", "features", 40),
    ("at_2_10", "features", 41),
    ("rare_codons_2_10", "features", 42),
    ("mw_da", "features", 43),
    ("pi", "features", 44),
    ("gravy", "features", 45),
    ("tm_segments", "features", 46),
    ("n_observable_peptides", "features", 47),
    ("flags", "flags", 73),
)
ORDER_1_1: tuple[str, ...] = tuple(c for c, _, _ in ORDER_1_1_TABLE)
ORDER_1_0: tuple[str, ...] = tuple(c for c, _, _ in sorted(ORDER_1_1_TABLE, key=lambda r: r[2]))
#: § 5.2 group ranges (1-based, inclusive)
GROUPS: dict[str, tuple[int, int]] = {
    "identity": (1, 14), "detection": (15, 26), "abundance": (27, 33), "function": (34, 39),
    "context": (40, 49), "features": (50, 72), "flags": (73, 73),
}  # fmt: skip
#: § PR 11a "Examples" (run A): 1-based positions before -> after
MOVES = {"detected": (48, 15), "detection_status": (50, 16), "abundance_log2": (61, 27),
         "cog_category": (68, 35)}  # fmt: skip

_SCHEMA_TYPE = {"str": "string", "int": "integer", "float": "number", "bool": "boolean",
                "list": "string"}  # fmt: skip


def _schema() -> dict:
    return json.loads(GENES_SCHEMA.read_text(encoding="utf-8"))


def _schema_names() -> tuple[str, ...]:
    return tuple(f["name"] for f in _schema()["fields"])


def _frame(n: int = 3) -> pd.DataFrame:
    """A small gene table in the internal representation (``gene_table_fixture``)."""
    import gene_table_fixture as gt

    t = gt.build_world(frames=gt.frames_with()).assemble()
    return t.head(n).reset_index(drop=True)


# ----------------------------------------------------------------------------- spec self-checks


def test_spec_table_is_a_permutation_of_73_columns() -> None:
    assert len(ORDER_1_1) == 73 and len(set(ORDER_1_1)) == 73
    assert sorted(p for _, _, p in ORDER_1_1_TABLE) == list(range(1, 74))


def test_spec_groups_are_contiguous_and_in_order() -> None:
    for i, (_, group, _) in enumerate(ORDER_1_1_TABLE, start=1):
        lo, hi = GROUPS[group]
        assert lo <= i <= hi, (i, group)
    assert list(dict.fromkeys(g for _, g, _ in ORDER_1_1_TABLE)) == list(GROUPS)


# ----------------------------------------------------------------------------- AC 11a.1


def test_schema_fields_follow_the_section_5_2_order() -> None:
    assert _schema_names() == ORDER_1_1


def test_genes_columns_follow_the_section_5_2_order() -> None:
    from gippa.genetable.assemble import GENES_COLUMNS

    assert GENES_COLUMNS == ORDER_1_1


def test_schema_order_matches_the_writer_order() -> None:
    from gippa.genetable.assemble import GENES_COLUMNS

    assert _schema_names() == GENES_COLUMNS


@pytest.mark.parametrize(("column", "move"), sorted(MOVES.items()))
def test_example_a_positions(column: str, move: tuple[int, int]) -> None:
    before, after = move
    assert ORDER_1_0.index(column) + 1 == before  # the spec's corrected 1.0 positions
    assert _schema_names().index(column) + 1 == after


def test_written_header_equals_section_5_2(tmp_path: Path) -> None:
    from gippa.io.genes_tsv import write_genes_tsv

    path = tmp_path / "genes.tsv"
    write_genes_tsv(path, _frame())
    header = path.read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
    assert tuple(header) == ORDER_1_1


@settings(
    max_examples=25, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(perm=st.permutations(ORDER_1_1))
def test_property_header_is_section_5_2_whatever_the_frame_order(
    tmp_path: Path, perm: list[str]
) -> None:
    """The writer emits schema order for any input column order, and the reader inverts it."""
    from gippa.io.genes_tsv import read_genes_tsv, write_genes_tsv

    frame = _frame()
    path = tmp_path / "genes_perm.tsv"
    write_genes_tsv(path, frame[list(perm)])
    header = path.read_text(encoding="utf-8").split("\n", 1)[0].split("\t")
    assert tuple(header) == ORDER_1_1
    back = read_genes_tsv(path)
    assert tuple(back.columns) == ORDER_1_1
    assert back["locus_tag"].tolist() == frame["locus_tag"].tolist()


def test_parquet_intermediate_uses_the_same_order(tmp_path: Path) -> None:
    import pyarrow.parquet as pq

    from gippa.io.genes_tsv import write_genes_parquet

    path = tmp_path / "genes.parquet"
    write_genes_parquet(path, _frame())
    assert tuple(pq.read_schema(path).names) == ORDER_1_1


# ---------------------------------------------------------------------- set and types unchanged


def test_column_set_and_names_unchanged() -> None:
    assert set(_schema_names()) == set(ORDER_1_0)


def test_schema_types_agree_with_the_writer_types() -> None:
    from gippa.genetable.assemble import GENES_TYPES

    for f in _schema()["fields"]:
        assert f.get("type", "string") == _SCHEMA_TYPE[GENES_TYPES[f["name"]]], f["name"]


# ----------------------------------------------------------------------------- AC 11a.5


def test_frictionless_validates_a_written_genes_tsv(tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    from gippa.io.genes_tsv import write_genes_tsv

    write_genes_tsv(tmp_path / "genes.tsv", _frame(5))
    descriptor = {k: v for k, v in _schema().items() if k != "coordinates"}
    schema = frictionless.Schema.from_descriptor(descriptor)
    assert [f.name for f in schema.fields] == list(ORDER_1_1)
    resource = frictionless.Resource(path="genes.tsv", basepath=str(tmp_path), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_frictionless_rejects_a_1_0_ordered_file_against_the_1_1_descriptor(
    tmp_path: Path,
) -> None:
    """The descriptor order is the file order: position-matched validation sees the change."""
    frictionless = pytest.importorskip("frictionless")
    from gippa.io.tsv import write_tsv

    write_tsv(tmp_path / "old.tsv", list(ORDER_1_0), [])
    descriptor = {k: v for k, v in _schema().items() if k != "coordinates"}
    resource = frictionless.Resource(
        path="old.tsv", basepath=str(tmp_path),
        schema=frictionless.Schema.from_descriptor(descriptor),
    )  # fmt: skip
    assert not resource.validate().valid


# ----------------------------------------------------------------------------- tables.md


def _genes_section_fields(text: str) -> list[str]:
    m = re.search(r"^## `genes`\n(.*?)(?=^## |\Z)", text, re.S | re.M)
    assert m, "no genes section in docs/reference/tables.md"
    return re.findall(r"^\| `([a-z0-9_]+)` \|", m.group(1), re.M)


def test_tables_md_lists_genes_fields_in_section_5_2_order() -> None:
    assert tuple(_genes_section_fields(TABLES_MD.read_text(encoding="utf-8"))) == ORDER_1_1


def test_tables_md_is_what_the_generator_renders_for_genes() -> None:
    sys.path.insert(0, str(REPO / "scripts" / "docs"))
    import gen_reference

    rendered = gen_reference.render_tables(SCHEMAS)
    committed = TABLES_MD.read_text(encoding="utf-8")
    assert _genes_section_fields(rendered) == _genes_section_fields(committed)


# ----------------------------------------------------------------------------- Q14


def _version(text: str) -> tuple[int, int]:
    major, minor = text.split(".")
    return int(major), int(minor)


def test_output_schema_version_bumped_minor_under_adr_0075() -> None:
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    major, minor = _version(OUTPUT_SCHEMA_VERSION)
    assert major == 1 and minor >= 3, OUTPUT_SCHEMA_VERSION


def test_datapackage_output_schema_version_matches() -> None:
    from gippa.schemas import OUTPUT_SCHEMA_VERSION

    package = json.loads((SCHEMAS / "datapackage.json").read_text(encoding="utf-8"))
    assert package["output_schema_version"] == OUTPUT_SCHEMA_VERSION
    assert _version(package["output_schema_version"]) >= (1, 3)
