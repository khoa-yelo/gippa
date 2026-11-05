"""M12-04 PR 11b: reading the 1.0 determinants feature *values* (lead decision, VALUE_RENAMES).

docs/specs/M12-04-outputs-revision.md > PR 11b; ADR-0075. The ``determinants_features`` columns
``promotech_max``, ``promotech_tu`` and ``mw`` became ``promoter_score``, ``promoter_score_tu``
and ``mw_da``. A 1.0 run also holds those names as *values*: ``determinants_model.
feature_or_block`` and the ``feature`` column of the ``determinants_pd`` and
``determinants_shap`` support files. ``RunView`` renames them on read (``table()`` and
``support()``), so ``gippa report`` and V15 work on a 1.0 tree.

Interface pinned here: ``gippa.io.migration.VALUE_RENAMES`` (1.0 feature name -> 1.1 name, the
same pairs as ``COLUMN_RENAMES["determinants_features"]``), ``VALUE_COLUMNS`` (table or support
stem -> its feature column) and ``migrate_values(table, frame) -> DataFrame``.

The 1.0 tree is the ``sim_plots`` snapshot run with the 1.1 names put back to 1.0 (headers of
``determinants_features`` and the three value columns), so the 1.0 and 1.1 reads can be compared
cell for cell.
"""

from __future__ import annotations

import gzip
import shutil
from collections.abc import Callable
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import pytest
from helpers.m12_04_renames import RENAMES, VALUE_RENAMES
from hypothesis import given
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SIM_PLOTS = REPO / "tests" / "data" / "snapshot" / "sim_plots" / "OUT"
BACK = {new: old for old, new in VALUE_RENAMES.items()}
#: the files holding feature names as values: (path under OUT, RunView reader, key, column)
VALUE_FILES = (
    ("tables/determinants_model.tsv", "table", "determinants_model", "feature_or_block"),
    ("plots/support/determinants_pd.tsv", "support", "determinants_pd", "feature"),
    ("plots/support/determinants_shap.tsv.gz", "support", "determinants_shap", "feature"),
)


# ---------------------------------------------------------------- the map


def test_value_renames_equal_the_determinants_features_column_renames() -> None:
    from gippa.io.migration import COLUMN_RENAMES
    from gippa.io.migration import VALUE_RENAMES as CODE

    assert dict(CODE) == dict(VALUE_RENAMES) == dict(RENAMES["determinants_features"])
    assert dict(COLUMN_RENAMES["determinants_features"]) == dict(VALUE_RENAMES)


def test_value_columns_name_the_three_feature_columns() -> None:
    from gippa.io.migration import VALUE_COLUMNS

    assert dict(VALUE_COLUMNS) == {key: col for _, _, key, col in VALUE_FILES}


# ---------------------------------------------------------------- migrate_values


@pytest.mark.parametrize(("table", "col"), [(k, c) for _, _, k, c in VALUE_FILES])
def test_migrate_values_renames_the_1_0_names_and_nothing_else(table: str, col: str) -> None:
    from gippa.io.migration import migrate_values

    frame = pd.DataFrame(
        {col: ["promotech_max", "cai", "mw", "promotech_tu", "all", "promoter_score"],
         "value": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]}
    )  # fmt: skip
    before = frame.copy()
    got = migrate_values(table, frame)
    assert list(got[col]) == [
        "promoter_score", "cai", "mw_da", "promoter_score_tu", "all", "promoter_score",
    ]  # fmt: skip
    assert list(got["value"]) == [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    pd.testing.assert_frame_equal(frame, before)  # the input is not changed in place


@pytest.mark.parametrize("table", ["determinants_features", "genes", "pairs", "determinants_xyz"])
def test_migrate_values_leaves_other_tables_alone(table: str) -> None:
    from gippa.io.migration import migrate_values

    frame = pd.DataFrame({"feature": ["mw", "promotech_max"], "feature_or_block": ["mw", "x"]})
    pd.testing.assert_frame_equal(migrate_values(table, frame), frame)


def test_migrate_values_without_the_feature_column_is_unchanged() -> None:
    from gippa.io.migration import migrate_values

    frame = pd.DataFrame({"locus_tag": ["g1"], "shap": [0.5]})
    pd.testing.assert_frame_equal(migrate_values("determinants_shap", frame), frame)
    empty = pd.DataFrame({"feature": pd.Series([], dtype=object)})
    pd.testing.assert_frame_equal(migrate_values("determinants_pd", empty), empty)


_NAMES = st.sampled_from(
    [*VALUE_RENAMES, *VALUE_RENAMES.values(), "cai", "tai", "tir_ostir", "all", "codon", "mw_x"]
)


@given(st.lists(_NAMES, max_size=12), st.sampled_from([k for _, _, k, _ in VALUE_FILES]))
def test_migrate_values_is_idempotent_and_maps_each_cell(names: list[str], table: str) -> None:
    from gippa.io.migration import VALUE_COLUMNS, migrate_values

    col = VALUE_COLUMNS[table]
    frame = pd.DataFrame({col: pd.Series(names, dtype=object)})
    once = migrate_values(table, frame)
    assert list(once[col]) == [VALUE_RENAMES.get(n, n) for n in names]
    assert not set(once[col]) & set(VALUE_RENAMES)
    pd.testing.assert_frame_equal(migrate_values(table, once), once)


# ---------------------------------------------------------------- a 1.0 tree through RunView


def _map_tsv(path: Path, fn: Callable[[list[str], list[list[str]]], None]) -> None:
    """Rewrite a TSV (plain or gzipped, by its magic bytes: some snapshot ``.tsv`` files are
    gzipped) by ``fn(header, rows)`` editing the cells in place."""
    raw = path.read_bytes()
    gz = raw[:2] == b"\x1f\x8b"
    raw = gzip.decompress(raw) if gz else raw
    lines = raw.decode("utf-8").splitlines()
    header, rows = lines[0].split("\t"), [ln.split("\t") for ln in lines[1:]]
    fn(header, rows)
    text = "\n".join(["\t".join(header), *("\t".join(r) for r in rows)]) + "\n"
    path.write_bytes(gzip.compress(text.encode("utf-8"), mtime=0) if gz else text.encode())


def _values_to_1_0(col: str) -> Callable[[list[str], list[list[str]]], None]:
    def fn(header: list[str], rows: list[list[str]]) -> None:
        i = header.index(col)
        for r in rows:
            r[i] = BACK.get(r[i], r[i])

    return fn


def _header_to_1_0(header: list[str], rows: list[list[str]]) -> None:
    back = {new: old for old, new in RENAMES["determinants_features"].items()}
    header[:] = [back.get(c, c) for c in header]


@pytest.fixture(scope="module")
def trees(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, Path]:
    """(the 1.1 sim_plots run, the same run written with the 1.0 names)."""
    root = tmp_path_factory.mktemp("pr11b_values")
    new = shutil.copytree(SIM_PLOTS, root / "OUT_1_1")
    old = shutil.copytree(SIM_PLOTS, root / "OUT_1_0")
    for rel, _, _, col in VALUE_FILES:
        _map_tsv(old / rel, _values_to_1_0(col))
    _map_tsv(old / "tables/determinants_features.tsv", _header_to_1_0)
    return new, old


def _cells(path: Path, col: str) -> set[str]:
    comp = "gzip" if path.read_bytes()[:2] == b"\x1f\x8b" else None
    frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, compression=comp)
    return set(frame[col])


@pytest.mark.parametrize(("rel", "how", "key", "col"), VALUE_FILES, ids=[k for *_, k, _ in
                                                                         VALUE_FILES])  # fmt: skip
def test_a_1_0_tree_reads_with_the_1_1_feature_names(
    trees: tuple[Path, Path], rel: str, how: str, key: str, col: str
) -> None:
    from gippa.io.run_view import RunView

    new, old = trees
    present = _cells(old / rel, col) & set(VALUE_RENAMES)
    assert present, f"precondition: the 1.0 {rel} holds a 1.0 feature name"
    want = getattr(RunView.from_outdir(new), how)(key)
    got = getattr(RunView.from_outdir(old), how)(key)
    assert want is not None and got is not None
    assert not set(got[col].astype(str)) & set(VALUE_RENAMES)
    assert {VALUE_RENAMES[n] for n in present} <= set(got[col].astype(str))
    pd.testing.assert_frame_equal(got, want)


def test_a_1_0_features_table_reads_with_the_1_1_columns(trees: tuple[Path, Path]) -> None:
    from gippa.io.run_view import RunView

    new, old = trees
    path = old / "tables/determinants_features.tsv"
    raw = path.read_bytes()
    head = (gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw).decode().split("\n")[0]
    assert set(VALUE_RENAMES) <= set(head.split("\t"))  # precondition: a 1.0 header
    want = RunView.from_outdir(new).table("determinants_features")
    got = RunView.from_outdir(old).table("determinants_features")
    assert want is not None and got is not None
    pd.testing.assert_frame_equal(got, want)


def test_v15_draws_the_same_panels_on_a_1_0_run(trees: tuple[Path, Path]) -> None:
    from gippa.config import resolve_config
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base, style
    from gippa.viz.base import Instance, Skip

    style.apply()
    load_plot_modules()
    spec = base.registered()["V15"]
    new, old = trees
    cfg = resolve_config()
    want = spec.fn(RunView.from_outdir(new), cfg)
    got = spec.fn(RunView.from_outdir(old), cfg)
    assert not isinstance(want, Skip), want
    assert not isinstance(got, Skip), got
    assert len(got) == len(want) == 1
    g, w = got[0], want[0]
    assert isinstance(g, Instance) and isinstance(w, Instance)
    assert g.title == w.title
    pd.testing.assert_frame_equal(g.data, w.data)
    assert "promoter_score" in set(g.data["feature"].astype(str))
    assert not set(g.data["feature"].astype(str)) & set(VALUE_RENAMES)
    labels = [ax.get_label() for ax in g.figure.axes]
    assert labels == [ax.get_label() for ax in w.figure.axes]
    assert "b_promoter_score" in labels
    xlabels = [ax.get_xlabel() for ax in g.figure.axes]
    assert xlabels == [ax.get_xlabel() for ax in w.figure.axes]
    plt.close(g.figure)
    plt.close(w.figure)
