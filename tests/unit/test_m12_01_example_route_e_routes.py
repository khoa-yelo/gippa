"""M12-01 route E (amendment 2026-11-06), work-split row 11a: route E in
``src/gippa/data/example_routes.yaml`` and the ``gippa example fetch`` printout (AC 1, 14), the
slow route test's helpers (``tests/helpers/example_routes.py``: ``ROUTES``, ``CONVERTED``,
``ABUNDANCE``, ``NOVEL``, the novel keys; AC 23), ``TIERS["E"]`` and the AC 22 check.

Spec docs/specs/M12-01-example-dataset.md section 2 (Routes, Route E's commands, ``expected/``)
and section 5 (Slow route test, route E). No network, no run.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import example_routes as er
from helpers.example_m12_01 import REPO, fill, flat, load_routes, options, subcommand
from helpers.example_route_e import route_e_misses

TITLE_E = "Route E: DIA-NN report, six-frame database"
COMMANDS_E = [  # section 2 > Route E's commands, verbatim (route E region, 2026-11-06: `_e` genome)
    "gippa convert --from diann {dir}/E_sixframe/report.parquet --out {out}/E/peptides.tsv "
    "--search-db same",
    "gippa preflight --genome {dir}/genome/nissle_r200k_e.fna.gz --peptides {out}/E/peptides.tsv "
    "--set 'modules=[core,bgc]' --set annotation.bakta.db=light "
    "--set annotation.terminators.tool=transtermhp "
    "--set annotation.promoter_calculator.enabled=false",
    "gippa run --genome {dir}/genome/nissle_r200k_e.fna.gz --peptides {out}/E/peptides.tsv "
    "--genome-relation same --outdir {out}/E --modules core,bgc "
    "--set annotation.bakta.db=light --terminators transtermhp "
    "--set annotation.promoter_calculator.enabled=false --threads 4",
]
NOVEL_KEYS = {"n_novel_candidates", "n_novel_called", "n_novel_unpaired", "novel_classes"}
BASE_KEYS = {"status", "detection_tier", "layers_missing", "convert_rows", "n_genes",
             "n_detected", "n_detected_2pep", "n_peptides_placed", "n_abundance", "warnings",
             "wall_time_s"}  # fmt: skip
CLASS_FIELDS = ("regime", "class", "estimated", "n_targets", "n_target_wins", "n_decoys",
                "q_floor", "floor_above_alpha")  # fmt: skip


@pytest.fixture(scope="module")
def routes() -> dict[str, dict[str, Any]]:
    return load_routes()


# ---------------------------------------------------------------- example_routes.yaml (AC 14)


def test_five_routes(routes: dict[str, dict[str, Any]]) -> None:
    assert tuple(sorted(routes)) == ("A", "B", "C", "D", "E")
    assert list(routes)[-1] == "E", "route E comes after D"


def test_route_e_title_and_commands(routes: dict[str, dict[str, Any]]) -> None:
    assert routes["E"]["title"] == TITLE_E
    assert routes["E"]["commands"] == [flat(c) for c in COMMANDS_E]


def test_route_e_is_route_a_on_the_six_frame_report(routes: dict[str, dict[str, Any]]) -> None:
    """They are route A's with ``A_diann`` -> ``E_sixframe``, ``{out}/A`` -> ``{out}/E`` and
    (route E region, 2026-11-06) ``nissle_r200k.fna.gz`` -> ``nissle_r200k_e.fna.gz``."""
    want = [c.replace("A_diann", "E_sixframe").replace("{out}/A", "{out}/E")
            .replace("nissle_r200k.fna.gz", "nissle_r200k_e.fna.gz")
            for c in routes["A"]["commands"]]  # fmt: skip
    assert routes["E"]["commands"] == want


def test_route_e_inputs(routes: dict[str, dict[str, Any]]) -> None:
    """AC 14 for E: ``--from diann``, ``--search-db same``, the six-frame report, ``--out
    {out}/E/peptides.tsv``; no extra flag on the run (no ``--annotation``)."""
    cmds = routes["E"]["commands"]
    assert [subcommand(c) for c in cmds] == ["convert", "preflight", "run"]
    conv = options(cmds[0])
    assert conv["--from"] == ["diann"] and conv["--search-db"] == ["same"]
    assert conv["--out"] == ["{out}/E/peptides.tsv"]
    assert "{dir}/E_sixframe/report.parquet" in cmds[0]
    run = options(cmds[2])
    assert run["--peptides"] == ["{out}/E/peptides.tsv"] and run["--outdir"] == ["{out}/E"]
    assert "--annotation" not in run and "--proteins" not in run


def test_fetch_printout_has_route_e() -> None:
    """AC 1 (five route headings): ``gippa example fetch`` prints route E's heading and its three
    commands with the absolute bundle directory."""
    from gippa import example

    d = Path("/x/gippa-example")
    lines = example.route_lines("1.1", d, example._routes())
    i = lines.index(f"# {TITLE_E}")
    assert lines[i + 1 : i + 4] == [fill(flat(c), str(d)) for c in COMMANDS_E]
    assert sum(ln.startswith("# Route ") for ln in lines) == 5


def test_build_bundle_paths_include_e_sixframe() -> None:
    """The build's layout check reads the YAML's ``{dir}`` paths."""
    paths = {p for r in load_routes().values() for c in r["commands"]
             for p in c.split() if p.startswith("{dir}/")}  # fmt: skip
    assert "{dir}/E_sixframe/report.parquet" in paths


# ---------------------------------------------------------------- helpers (AC 23)


def test_helper_route_sets() -> None:
    assert er.ROUTES == ("A", "B", "C", "D", "E")
    assert sorted(er.CONVERTED) == list("ACDE")
    assert sorted(er.ABUNDANCE) == list("ABDE")
    assert sorted(er.NOVEL) == ["E"]


def test_required_keys_of_e_add_the_novel_keys() -> None:
    assert er.required_keys("E") == BASE_KEYS | NOVEL_KEYS


@pytest.mark.parametrize(("route", "drop"), [("A", set()), ("B", {"convert_rows"}),
                                             ("C", {"n_abundance"}), ("D", set())])  # fmt: skip
def test_required_keys_of_a_to_d_unchanged(route: str, drop: set[str]) -> None:
    assert er.required_keys(route) == BASE_KEYS - drop


def test_novel_tolerances() -> None:
    """Section 2 > ``expected/``: candidates and unpaired max(1, 2 %), called exact, classes
    informational."""
    t = er.TOLERANCES
    assert t["n_novel_candidates"] == {"abs": 1, "rel": 0.02}
    assert t["n_novel_unpaired"] == {"abs": 1, "rel": 0.02}
    assert t["n_novel_called"] == "exact"
    assert t["novel_classes"] == "informational"


def _classes() -> list[dict[str, Any]]:
    return [{"regime": "picked_decoy", "class": "antisense", "estimated": True, "n_targets": 3,
             "n_target_wins": 2, "n_decoys": 1, "q_floor": 0.5, "floor_above_alpha": True,
             "extra": "dropped"}]  # fmt: skip


def _run_dir(tmp: Path, *, novel: bool = True) -> Path:
    (tmp / "tables").mkdir(parents=True)
    (tmp / "novel").mkdir()
    (tmp / "mapping").mkdir()
    (tmp / "mapping" / "peptide_hits.tsv").write_text("peptide\nAAK\n")
    man = {"status": "completed", "detection_tier": "mayu", "layers_missing": [], "warnings": [],
           "novel_orfs": {"classes": _classes()} if novel else None}  # fmt: skip
    (tmp / "run_manifest.json").write_text(json.dumps(man))
    (tmp / "tables" / "genes.tsv").write_text(
        "locus_tag\tdetected\tdetected_2pep\tabundance_raw\tabundance_imputed\n"
        "g1\ttrue\ttrue\t10.5\tfalse\n")  # fmt: skip
    (tmp / "tables" / "novel_orfs.tsv").write_text(
        "# gippa novel_orfs\n# schema 1\norf_id\tcalled\tclass\n"
        "o1\tfalse\tantisense\no2\tfalse\tintergenic\no3\ttrue\tantisense\n")  # fmt: skip
    pd.DataFrame({"orf_id": ["o1", "o2", "o3"],
                  "reason": ["no_makedb_pair", None, "paired_share_low"]}).to_parquet(
        tmp / "novel" / "novel_class_fdr.parquet")  # fmt: skip
    (tmp / "pep.tsv").write_text("peptide\trun\nAAK\tWT_1\n")
    return tmp


def test_measure_reads_the_novel_keys_for_e(tmp_path: Path) -> None:
    got = er.measure("E", _run_dir(tmp_path), tmp_path / "pep.tsv")
    assert got["n_novel_candidates"] == 3  # `#` lines skipped
    assert got["n_novel_called"] == 1
    assert got["n_novel_unpaired"] == 1
    assert got["novel_classes"] == [{k: _classes()[0][k] for k in CLASS_FIELDS}]
    assert got["convert_rows"] == 1 and got["n_abundance"] == 1


def test_measure_reads_no_novel_key_for_a_to_d(tmp_path: Path) -> None:
    d = _run_dir(tmp_path)
    for route in "ABCD":
        assert not NOVEL_KEYS & set(er.measure(route, d, d / "pep.tsv")), route


def test_measure_without_novel_classes(tmp_path: Path) -> None:
    got = er.measure("E", _run_dir(tmp_path, novel=False), tmp_path / "pep.tsv")
    assert got["novel_classes"] == []


def test_expected_e_round_trip(tmp_path: Path) -> None:
    got = er.measure("E", _run_dir(tmp_path), tmp_path / "pep.tsv") | {"wall_time_s": 1.0}
    doc = er.as_expected(got, "E")
    assert set(doc) == er.required_keys("E")
    assert er.compare(got, doc, "E") == []
    assert er.compare(got | {"n_novel_called": 2}, doc, "E")[0].startswith("n_novel_called")
    del doc["n_novel_unpaired"]
    assert er.compare(got, doc, "E") == ["n_novel_unpaired: absent from expected/"]


# ---------------------------------------------------------------- the slow test (AC 22)


def test_e2e_tier_of_route_e() -> None:
    path = REPO / "tests" / "e2e" / "test_m12_01_example_routes.py"
    spec = importlib.util.spec_from_file_location("gippa_e2e_example_routes", path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # its dataclasses look the module up
    try:
        spec.loader.exec_module(mod)
    finally:
        del sys.modules[spec.name]
    assert mod.TIERS["E"] == "decoy_picked"


def test_e2e_applies_the_ac22_check() -> None:
    text = (REPO / "tests" / "e2e" / "test_m12_01_example_routes.py").read_text()
    assert "route_e_misses" in text


_OK = {"detection_tier": "decoy_picked", "n_novel_called": 0, "n_novel_candidates": 4,
       "n_novel_unpaired": 1, "novel_classes": [{"regime": "picked_decoy"}]}  # fmt: skip


@pytest.mark.parametrize(
    ("change", "field"),
    [({}, None), ({"n_novel_candidates": 0, "novel_classes": []}, "n_novel_candidates"),
     ({"detection_tier": "mayu"}, "detection_tier"),
     ({"n_novel_called": 1}, "n_novel_called"),
     ({"n_novel_unpaired": 4}, "n_novel_unpaired"),
     ({"novel_classes": [{"regime": "denovo_subgroup"}]}, "novel_classes"),
     ({"n_novel_candidates": None}, "n_novel_candidates")],
)  # fmt: skip
def test_route_e_misses(change: dict[str, Any], field: str | None) -> None:
    got = route_e_misses(_OK | change)
    assert [m.split(":")[0] for m in got] == ([field] if field else [])
