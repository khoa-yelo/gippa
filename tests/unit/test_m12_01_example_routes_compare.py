"""Fast checks of the M12-01 slow route test's helpers (tests/helpers/example_routes.py): the
`expected/` comparison, the measurement and the printed-block split (spec section 2, AC 17)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from helpers import example_routes as er
from hypothesis import given
from hypothesis import strategies as st


def _doc(route: str, **values: object) -> dict[str, dict[str, object]]:
    base = {"status": "completed", "detection_tier": "mayu", "layers_missing": ["a", "b"],
            "convert_rows": 5, "n_genes": 100, "n_detected": 100, "n_detected_2pep": 0,
            "n_peptides_placed": 0, "n_abundance": 0, "warnings": [],
            "wall_time_s": 1.0}  # fmt: skip
    return er.as_expected(base | values, route)


def _measured(**values: object) -> dict[str, object]:
    return {k: v["value"] for k, v in _doc("A").items()} | values


@pytest.mark.parametrize(("got", "exp", "ok"), [(102, 100, True), (103, 100, False),
                                                (1, 0, True), (2, 0, False)])  # fmt: skip
def test_count_tolerance_max_1_or_2_percent(got: int, exp: int, ok: bool) -> None:
    doc = _doc("A", n_detected=exp)
    assert (er.compare(_measured(n_detected=got), doc, "A") == []) is ok


def test_exact_set_and_informational() -> None:
    doc = _doc("A")
    doc["n_detected_2pep"]["tolerance"] = "informational"
    m = _measured(layers_missing=["b", "a"], warnings=["x"], wall_time_s=999, n_detected_2pep=50)
    assert er.compare(m, doc, "A") == []
    bad = _measured(status="failed", layers_missing=["a"])
    assert [x.split(":")[0] for x in er.compare(bad, doc, "A")] == ["layers_missing", "status"]


def test_bad_values_and_missing_keys_are_misses() -> None:
    doc = _doc("A")
    assert er.compare(_measured(n_genes="NA"), doc, "A") == [
        "n_genes: not numeric (got 'NA', expected 100)"]  # fmt: skip
    del doc["convert_rows"]
    assert er.compare(_measured(), doc, "A") == ["convert_rows: absent from expected/"]
    assert er.compare(_measured(), _doc("B"), "B") == []  # B has no convert_rows
    doc = _doc("C")
    assert "n_abundance" not in doc
    doc["n_genes"] = {"value": 100, "tolerance": "max(1, 1 %)"}
    assert er.compare(_measured(), doc, "C") == ["n_genes: unknown tolerance 'max(1, 1 %)'"]
    m = _measured()
    del m["n_genes"]
    assert er.compare(m, _doc("C"), "C") == ["n_genes: not measured (expected 100)"]


@given(exp=st.integers(0, 10_000), delta=st.integers(-500, 500))
def test_abs_rel_bound(exp: int, delta: int) -> None:
    assert er.within(exp + delta, exp, {"abs": 1, "rel": 0.02}) == (
        abs(delta) <= max(1, 0.02 * exp)
    )


def test_measure_round_trips_through_as_expected(tmp_path: Path) -> None:
    (tmp_path / "tables").mkdir()
    (tmp_path / "mapping").mkdir()
    (tmp_path / "run_manifest.json").write_text(json.dumps({
        "status": "completed", "detection_tier": "mayu", "layers_missing": ["b", "a"],
        "warnings": ["s:w"]}))  # fmt: skip
    (tmp_path / "tables" / "genes.tsv").write_text(
        "locus_tag\tdetected\tdetected_2pep\tabundance_raw\tabundance_imputed\n"
        "g1\ttrue\ttrue\t10.5\tfalse\ng2\ttrue\tfalse\t3\ttrue\n"
        "g3\tfalse\tfalse\tNA\tfalse\n")  # fmt: skip
    (tmp_path / "mapping" / "peptide_hits.tsv").write_text("peptide\nAAK\nAAK\nCCR\n")
    (tmp_path / "pep.tsv").write_text("peptide\trun\nAAK\tWT_1\nCCR\tWT_1\n")
    got = er.measure("A", tmp_path, tmp_path / "pep.tsv") | {"wall_time_s": 3.0}
    assert got == {"status": "completed", "detection_tier": "mayu", "layers_missing": ["a", "b"],
                   "warnings": ["s:w"], "convert_rows": 2, "n_genes": 3, "n_detected": 2,
                   "n_detected_2pep": 1, "n_abundance": 1, "n_peptides_placed": 2,
                   "wall_time_s": 3.0}  # fmt: skip
    assert er.compare(got, er.as_expected(got, "A"), "A") == []
    assert "convert_rows" not in er.measure("B", tmp_path, tmp_path / "pep.tsv")
    assert "n_abundance" not in er.measure("C", tmp_path)


def test_printed_routes_split_by_yaml_titles() -> None:
    from gippa import example

    routes = example._routes()
    lines = example.route_lines("1.1", Path("/x/gippa-example"), routes)
    got = er.printed_routes(lines, routes)
    assert sorted(got) == list(er.ROUTES)
    for r in routes:
        assert len(got[str(r["id"])]) == len(r["commands"])
    title = f"# {routes[0]['title']}"
    tail = er.printed_routes([title, "gippa convert x", "not a command", "gippa run y"], routes)
    assert tail == {str(routes[0]["id"]): ["gippa convert x"]}


def test_argv_and_option() -> None:
    extra = ["--profile", "farmshare_local,apptainer"]
    conv = er.argv("gippa convert --from diann a --out /o/A/p.tsv", "/bin/g", extra)
    assert conv == ["/bin/g", "convert", "--from", "diann", "a", "--out", "/o/A/p.tsv"]
    assert er.option(conv, "--out") == Path("/o/A/p.tsv")
    assert er.option(conv, "--outdir") is None
    run = er.argv("gippa run --genome g --outdir o", "/bin/g", extra)
    assert run == ["/bin/g", "run", "--genome", "g", "--outdir", "o", *extra]
