"""M12-01 route E (amendment 2026-11-06), work-split row 11a: ``scripts/example/build_example.py``
writes ``E_sixframe/report.parquet`` from the six-frame arm A4.

Spec docs/specs/M12-01-example-dataset.md: subset rules 3E, 7 (route E), 8 and 9 with the eight
worked examples, section 2 > PROVENANCE.md and the licence statement, section 3 (``sources()``
gains ``a4_report``, ``a4_db_id_map``, ``a4_strain_db``); AC 20, 21, 21b, 21c. ADR-0074 amendment
2026-11-06.

Seams pinned here (the spec leaves the names to the implementer; these are the ones the tests
read, mirroring route A's):

- ``SELECTION.json > counts > E``: ``rows_total``, ``rows_kept``, ``target_rows_kept``,
  ``decoy_rows_kept``, ``rows_with_orf_ids``, ``rows_per_run`` (as ``counts > A``);
- ``SELECTION.json > orf_id_rewrite``: ``columns`` (the columns rule 9 rewrote) and
  ``distinct_ids`` (the number of distinct ORF IDs rewritten);
- ``SELECTION.json > sources``: ``a4_report``, ``a4_db_id_map``, ``a4_strain_db``;
- ``SELECTION.json > rewritten_columns > E_sixframe/report.parquet`` (rule 8, as for A);
- ``bundle.json > counts > E`` equals ``SELECTION.json``'s; ``bundles/sizes.json`` has
  ``E_sixframe``;
- a failed rule-7 or rule-9 check exits 1 with ``build_example: ...`` naming route E or the ID.

The A4 stand-in is synthetic (tests/helpers/example_a4.py); no A4 row is ever read (R11). The
build makes no statistical claim, so there is no calibration test; rule 3E and 9 hold for every
placement, which the Hypothesis property checks.
"""

from __future__ import annotations

import contextlib
import importlib.util
import io
import json
import re
import tarfile
import tempfile
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow.parquet as pq
import pytest
from helpers.example_a4 import (
    CHROM,
    FAA_REL,
    KEPT,
    MAP_REL,
    ORF_COLUMNS,
    REPORT_REL,
    REWRITTEN_KEYS,
    ROWS,
    RUNS,
    Entry,
    Row,
    bad_class,
    entries,
    make_pilot_e,
    region_contig,
    set_route_e,
    with_rows,
)
from helpers.example_m12_01 import REPO, flat
from hypothesis import HealthCheck, Phase, given, settings
from hypothesis import strategies as st

SCRIPT = REPO / "scripts" / "example" / "build_example.py"
START, END = 1_965_000, 2_165_000  # the default region
C = "nissle1917_NZ_CP058217_1_1965000_2165000"
E_REPORT = "E_sixframe/report.parquet"
S2, E2, L2 = 10_000, 40_000, 50_000  # the small region of the variant builds


def _minor() -> str:
    import gippa

    m = re.match(r"(\d+)\.(\d+)", str(gippa.__version__))
    assert m is not None
    return f"{m.group(1)}.{m.group(2)}"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gippa_build_example_route_e", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@dataclass(frozen=True)
class Built:
    root: Path
    code: int
    out: str
    work: Path
    rel: Path
    tree: Path | None

    def sel(self) -> dict[str, Any]:
        assert self.tree is not None
        return json.loads((self.tree / "SELECTION.json").read_text())  # type: ignore[no-any-return]

    def report(self) -> list[dict[str, Any]]:
        assert self.tree is not None
        path = self.tree / E_REPORT
        assert path.is_file(), f"{E_REPORT} missing from the bundle"
        return pq.read_table(path).to_pylist()  # type: ignore[no-any-return]

    @property
    def archive(self) -> Path:
        return self.rel / f"gippa-example-{_minor()}.tar.gz"


def build(root: Path, tmp: Path, region: tuple[int, int] | None = None) -> Built:
    """``build_example.main`` with ``GIPPA_ROOT=root``; the archive unpacked when it exists."""
    work, rel = tmp / "work", tmp / "rel"
    argv = ["--out", str(work), "--bundle-out", str(rel), "--allow-dirty"]
    if region is not None:
        argv += ["--start", str(region[0]), "--end", str(region[1])]
    mod, buf = _load(), io.StringIO()
    set_route_e(mod, *(region or (START, END)))  # A4 sits in the shared region (amendment b)
    with (
        pytest.MonkeyPatch.context() as mp,
        contextlib.redirect_stdout(buf),
        contextlib.redirect_stderr(buf),
    ):
        mp.setenv("GIPPA_ROOT", str(root))
        mp.delenv("SLURM_JOB_ID", raising=False)
        try:
            code = mod.main(argv)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
    archive = rel / f"gippa-example-{_minor()}.tar.gz"
    tree = None
    if archive.is_file():
        with tarfile.open(archive) as tar:
            tar.extractall(tmp / "unpacked", filter="data")
        tree = tmp / "unpacked" / f"gippa-example-{_minor()}"
    return Built(root, int(code or 0), buf.getvalue(), work, rel, tree)


@pytest.fixture(scope="module")
def base(tmp_path_factory: pytest.TempPathFactory) -> Built:
    """The default region with the spec's rule-9 examples (one build)."""
    tmp = tmp_path_factory.mktemp("route_e")
    root = (tmp / "gippa_root").resolve()
    make_pilot_e(root, START, END, END + 5_000)
    b = build(root, tmp)
    assert b.code == 0, b.out
    return b


def _small(tmp: Path, **a4: Any) -> Built:
    root = (tmp / "gippa_root").resolve()
    make_pilot_e(root, S2, E2, L2, **a4)
    return build(root, tmp, (S2, E2))


def _by_peptide(rows: list[dict[str, Any]]) -> dict[str, list[dict[str, Any]]]:
    out: dict[str, list[dict[str, Any]]] = {}
    for r in rows:
        out.setdefault(str(r["Stripped.Sequence"]), []).append(r)
    return out


ENT = entries(START, END)


# ---------------------------------------------------------------- AC 20: rule 3E


def test_route_e_report_is_in_the_bundle(base: Built) -> None:
    assert base.tree is not None and (base.tree / E_REPORT).is_file()
    out = base.sel()["outputs"]
    assert E_REPORT in out, sorted(out)


def test_kept_rows_are_exactly_the_region_rows(base: Built) -> None:
    """AC 20: (a) region CDS + region ORF and (b) a region ORF only are kept, target and decoy
    alike; ORFs crossing either region end, an ORF or CDS outside and a ``Cont_`` ID drop the
    row."""
    got = _by_peptide(base.report())
    assert set(got) == set(KEPT), sorted(set(got) ^ set(KEPT))
    for pep, rows in got.items():
        assert sorted(r["Run"] for r in rows) == sorted(RUNS), pep


@pytest.mark.parametrize("row", ROWS, ids=lambda r: r.peptide)
def test_rule_3e_per_row(base: Built, row: Row) -> None:
    got = _by_peptide(base.report())
    assert (row.peptide in got) is row.kept, (row.ids, row.decoy)
    if row.kept:
        assert {int(r["Decoy"]) for r in got[row.peptide]} == {row.decoy}


def test_minus_strand_region_orf_is_kept_as_target_and_decoy(base: Built) -> None:
    got = _by_peptide(base.report())
    assert ENT["O_AS"].strand == "-"
    assert "ASTARGETK" in got and "ASDECOYK" in got


def test_route_e_counts(base: Built) -> None:
    e = base.sel()["counts"]["E"]
    per_run = len(KEPT)
    assert e["rows_total"] == len(ROWS) * len(RUNS)
    assert e["rows_kept"] == per_run * len(RUNS)
    n_dec = sum(r.decoy for r in ROWS if r.kept)
    assert (e["target_rows_kept"], e["decoy_rows_kept"]) == (
        (per_run - n_dec) * len(RUNS), n_dec * len(RUNS))  # fmt: skip
    orf = {k for k, v in ENT.items() if v.kind == "orf"}
    n_orf = sum(1 for r in ROWS if r.kept and set(r.ids) & orf)
    assert e["rows_with_orf_ids"] == n_orf * len(RUNS)
    assert e["rows_per_run"] == {r: per_run for r in RUNS}


def test_all_columns_are_kept(base: Built) -> None:
    src = pq.read_table(base.root / REPORT_REL).column_names
    assert base.tree is not None
    assert pq.read_table(base.tree / E_REPORT).column_names == src


def test_other_cells_are_unchanged(base: Built) -> None:
    """Rules 8 and 9 change only ``File.Name`` and the ORF IDs; every other cell is the source's."""
    src = [r for r in pq.read_table(base.root / REPORT_REL).to_pylist()
           if r["Stripped.Sequence"] in KEPT]  # fmt: skip
    skip = {"File.Name", *ORF_COLUMNS}

    def key(r: dict[str, Any]) -> tuple[str, str]:
        return str(r["Run"]), str(r["Stripped.Sequence"])

    got = sorted(base.report(), key=key)
    want = sorted(src, key=key)
    assert len(got) == len(want)
    for g, w in zip(got, want, strict=True):
        assert {k: v for k, v in g.items() if k not in skip} == {
            k: v for k, v in w.items() if k not in skip}  # fmt: skip


def test_route_a_is_unchanged_by_route_e(base: Built) -> None:
    """Rule 3E: route A still reads A2 (no ORF ID in ``A_diann``)."""
    assert base.tree is not None
    ids = pq.read_table(base.tree / "A_diann" / "report.parquet").column("Protein.Ids").to_pylist()
    assert ids and not any("GIPPAORF_" in str(x) for x in ids)


# ---------------------------------------------------------------- AC 21: rule 9


# The spec's eight examples, literally (S = 1965000, C as above): (input Protein.Ids, the row's
# peptide in the fixture, expected Protein.Ids or None when dropped).
SPEC_EXAMPLES = [
    ("CDS_X", "CDSTARGETK", "CDS_X"),
    ("GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m", "ASDECOYK",
     f"GIPPAORF_antisense_{C}_101_400_m"),
    ("CDS_X;GIPPAORF_cds_outframe_NZ_CP058217.1_2000002_2000301_p", "MIXTARGETK",
     f"CDS_X;GIPPAORF_cds_outframe_{C}_35002_35301_p"),
    ("GIPPAORF_intergenic_NZ_CP058217.1_2164900_2165120_p", "XENDK", None),
    ("CDS_X;GIPPAORF_antisense_NZ_CP058217.1_10_309_m", "FARORFK", None),
    ("Cont_P1;GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m", "CONTORFK", None),
    ("GIPPAORF_intergenic_NZ_CP058217.1_1965001_1965300_p", "FIRSTBASEK",
     f"GIPPAORF_intergenic_{C}_1_300_p"),
    ("GIPPAORF_intergenic_NZ_CP058217.1_2164701_2165000_p", "LASTBASEK",
     f"GIPPAORF_intergenic_{C}_199701_200000_p"),
]  # fmt: skip


def test_fixture_contig_is_the_spec_contig() -> None:
    assert region_contig(START, END) == C


@pytest.mark.parametrize(("ids", "peptide", "want"), SPEC_EXAMPLES,
                         ids=[x[1] for x in SPEC_EXAMPLES])  # fmt: skip
def test_spec_rule_9_examples(base: Built, ids: str, peptide: str, want: str | None) -> None:
    src = {r["Protein.Ids"] for r in pq.read_table(base.root / REPORT_REL).to_pylist()
           if r["Stripped.Sequence"] == peptide}  # fmt: skip
    assert src == {ids}, "fixture drift: the row does not carry the spec's input"
    got = _by_peptide(base.report()).get(peptide)
    if want is None:
        assert got is None, f"{ids} must be dropped"
    else:
        assert got is not None, f"{ids} must be kept"
        assert {r["Protein.Ids"] for r in got} == {want}


def test_spec_parse_example() -> None:
    """The spec's second example: ``parse_orf_id`` gives ``(antisense, C, 100, 400, -)``."""
    from gippa.novel.fdr import parse_orf_id

    assert parse_orf_id(f"GIPPAORF_antisense_{C}_101_400_m") == ("antisense", C, 100, 400, "-")
    assert ENT["O_AS"].parsed(START, END) == ("antisense", C, 100, 400, "-")


@pytest.mark.parametrize("column", ORF_COLUMNS)
def test_every_orf_id_round_trips_through_parse_orf_id(base: Built, column: str) -> None:
    """AC 21: with the real parser, every ORF ID of ``Protein.Ids``, ``Protein.Group`` and any
    other string column holding one gives the region contig, ``start1 - S - 1``, ``end1 - S``,
    class and strand unchanged."""
    from gippa.novel.fdr import parse_orf_id

    want = {ENT[k].rewritten(START, END): ENT[k].parsed(START, END) for k in REWRITTEN_KEYS}
    seen = set()
    for r in base.report():
        for pid in str(r[column]).split(";"):
            if pid.startswith("GIPPAORF_"):
                assert pid in want, (column, pid)
                assert parse_orf_id(pid) == want[pid], pid
                seen.add(pid)
    if column == "Protein.Ids":
        assert seen == set(want)


def test_ids_rewritten_in_order_with_cds_ids_unchanged(base: Built) -> None:
    by = _by_peptide(base.report())
    for row in ROWS:
        if not row.kept:
            continue
        want = ";".join(ENT[k].rewritten(START, END) for k in row.ids)
        group = ";".join(ENT[k].rewritten(START, END) for k in (row.group or row.ids))
        for r in by[row.peptide]:
            assert r["Protein.Ids"] == want, row.peptide
            assert r["Protein.Group"] == group, row.peptide
            assert r["Protein.Names"] == want, row.peptide


def test_minus_strand_offset_and_suffix(base: Built) -> None:
    """Minus strand: both ends shift by ``S`` and the ``_m`` suffix stays; plus keeps ``_p``."""
    got = {r["Protein.Ids"] for r in _by_peptide(base.report())["NTERMK"]}
    e = ENT["O_NT"]
    assert got == {f"GIPPAORF_nterm_extension_{C}_{e.s1 - START}_{e.e1 - START}_m"}
    last = {r["Protein.Ids"] for r in _by_peptide(base.report())["LASTBASEK"]}
    assert last == {f"GIPPAORF_intergenic_{C}_199701_200000_p"}


def test_no_chromosome_orf_id_is_left(base: Built) -> None:
    for r in base.report():
        for col, v in r.items():
            if isinstance(v, str):
                assert not re.search(rf"GIPPAORF_\w+?_{re.escape(CHROM)}_", v), (col, v)


def test_nterm_extension_is_not_translation_checked(base: Built) -> None:
    """Rule 9 check 3 skips ``nterm_extension`` (its record is upstream part + junction)."""
    assert "NTERMK" in _by_peptide(base.report())


def test_selection_records_the_rewrite(base: Built) -> None:
    sel = base.sel()
    rw = sel["orf_id_rewrite"]
    assert set(rw["columns"]) == set(ORF_COLUMNS)
    assert rw["distinct_ids"] == len(REWRITTEN_KEYS) == 5
    assert "Route E" in sel["rule"] and "region coordinates" in flat(sel["rule"])


# ---------------------------------------------------------------- AC 21b: rule 9 checks


def _fails_naming(b: Built, *needles: str) -> None:
    assert b.code == 1, b.out
    assert any(n in b.out for n in needles), (needles, b.out)
    assert not b.archive.exists()


def test_small_region_builds(tmp_path: Path) -> None:
    """The variant fixture without a fault builds (so the failures below are the checks')."""
    b = _small(tmp_path)
    assert b.code == 0, b.out
    assert set(_by_peptide(b.report())) == set(KEPT)


def test_map_coordinates_disagreeing_with_the_id_fail(tmp_path: Path) -> None:
    e = entries(S2, E2)["O_AS"]
    _fails_naming(_small(tmp_path, map_shift="O_AS"), e.db_id, e.rewritten(S2, E2))


def test_translation_differing_from_the_record_fails(tmp_path: Path) -> None:
    e = entries(S2, E2)["O_OUT"]
    _fails_naming(_small(tmp_path, bad_translation="O_OUT"), e.db_id, e.rewritten(S2, E2))


def test_rewritten_id_that_does_not_parse_fails(tmp_path: Path) -> None:
    """A rewritten ID that ``parse_orf_id`` (the real one) cannot read back fails the build."""
    from gippa.novel.fdr import parse_orf_id

    ents = entries(S2, E2)
    ents["O_OUT"] = bad_class(ents["O_OUT"])
    assert parse_orf_id(ents["O_OUT"].rewritten(S2, E2)) is None
    _fails_naming(_small(tmp_path, ents=ents), ents["O_OUT"].db_id, ents["O_OUT"].rewritten(S2, E2))


# ---------------------------------------------------------------- AC 21c: rule 7 for route E


@pytest.mark.parametrize("variant", [{"decoys": False}, {"targets": False}, {"orf_rows": False}],
                         ids=["no_decoy", "no_target", "no_orf_row"])  # fmt: skip
def test_empty_route_e_fails(tmp_path: Path, variant: dict[str, bool]) -> None:
    rows = with_rows(ROWS, **variant)
    assert rows, variant
    b = _small(tmp_path, rows=rows)
    assert b.code == 1, b.out
    line = next((ln for ln in b.out.splitlines() if "empty route" in ln), "")
    assert re.search(r"\bE\b", line.split("empty route", 1)[-1]), b.out
    assert not b.archive.exists()


# ---------------------------------------------------------------- rule 8, sources, counts


def test_route_e_file_names_are_base_names(base: Built) -> None:
    assert {r["File.Name"] for r in base.report()} == {f"{r}.mzML" for r in RUNS}
    assert "File.Name" in base.sel()["rewritten_columns"][E_REPORT]


def test_a4_sources_are_recorded(base: Built) -> None:
    src = base.sel()["sources"]
    for key, rel in (("a4_report", REPORT_REL), ("a4_db_id_map", MAP_REL),
                     ("a4_strain_db", FAA_REL)):  # fmt: skip
        assert src[key]["path"] == rel.as_posix(), key
        assert re.fullmatch(r"[0-9a-f]{64}", src[key]["sha256"]), key


def test_bundle_json_and_sizes_have_route_e(base: Built) -> None:
    assert base.tree is not None
    meta = json.loads((base.tree / "bundle.json").read_text())
    assert meta["counts"]["E"] == base.sel()["counts"]["E"]
    sizes = json.loads((base.work / "nissle_r200k" / "bundles" / "sizes.json").read_text())
    assert sizes.get("E_sixframe", 0) > 0, sizes


# ---------------------------------------------------------------- PROVENANCE, licence, README


def test_provenance_names_a4(base: Built) -> None:
    assert base.tree is not None
    text = flat((base.tree / "PROVENANCE.md").read_text())
    for needle in ("A4", "--six-frame-orfs", "diann_A4.cfg", "1918156", "DIA-NN 2.3.1"):
        assert needle in text, needle
    for key in ("a4_report", "a4_db_id_map", "a4_strain_db"):
        assert base.sel()["sources"][key]["sha256"] in text, key
    assert re.search(r"\b0 novel|called 0|0 called", text), "the genome-wide A4 run called 0"


def test_licence_covers_e_sixframe_and_the_renumbering(base: Built) -> None:
    assert base.tree is not None
    text = flat((base.tree / "LICENSE-DATA.md").read_text())
    for needle in ("`E_sixframe/`", "six-frame ORFs", "made by GIPPA, not by the depositors",
                   "renumbered to the region's coordinates", "PXD078459"):  # fmt: skip
        assert needle in text, needle


def test_readme_states_the_route_e_caveats(base: Built) -> None:
    assert base.tree is not None
    text = flat((base.tree / "README.md").read_text())
    assert "Route E: DIA-NN report, six-frame database" in text
    for needle in ("six-frame", "1,965,000", "`mayu`", "floor"):
        assert needle in text, needle


# ---------------------------------------------------------------- property: rules 3E and 9

_CLASSES = ("antisense", "cds_outframe", "intergenic", "nterm_extension")


@st.composite
def _orfs(draw: st.DrawFn) -> list[Entry]:
    out = []
    for i in range(draw(st.integers(1, 6))):
        n = 3 * draw(st.integers(20, 100))
        s0 = draw(st.integers(max(0, S2 - 400), min(E2 + 100, L2 - n)))
        out.append(Entry(f"G{i}", "orf", s0 + 1, s0 + n, draw(st.sampled_from("+-")),
                         cls=draw(st.sampled_from(_CLASSES)),
                         first_m=draw(st.booleans())))  # fmt: skip
    return out


# Each example is a whole build (about 0.5 s); no shrink phase, so a red run stays short.
@settings(max_examples=6, deadline=None, phases=(Phase.explicit, Phase.reuse, Phase.generate),
          suppress_health_check=[HealthCheck.too_slow, HealthCheck.data_too_large])  # fmt: skip
@given(orfs=_orfs(), decoy=st.lists(st.booleans(), min_size=6, max_size=6))
def test_property_rule_3e_and_9(orfs: list[Entry], decoy: list[bool]) -> None:
    """For any ORF placement (inside, crossing an end, outside; either strand; any class): a row
    is kept iff every ID is a region entry, and every kept ORF ID is the map interval shifted by
    ``-S`` on the region contig, read back by ``parse_orf_id``."""
    from gippa.novel.fdr import parse_orf_id

    ents = entries(S2, E2)
    taken = {e.db_id for e in ents.values()}
    rows = list(ROWS)
    for i, o in enumerate(orfs):
        if o.db_id in taken:
            continue
        taken.add(o.db_id)
        ents[o.key] = o
        tag = "ABCDEF"[i]
        ids = ("CDS_X", o.key) if i % 2 else (o.key,)
        rows.append(Row(f"GEN{tag}K", ids, int(decoy[i]), o.in_region(S2, E2)))
    with tempfile.TemporaryDirectory() as d:
        root = (Path(d) / "gippa_root").resolve()
        make_pilot_e(root, S2, E2, L2, ents=ents, rows=rows)
        b = build(root, Path(d), (S2, E2))
        assert b.code == 0, b.out
        got = _by_peptide(b.report())
    assert set(got) == {r.peptide for r in rows if r.kept}
    for row in rows:
        if not row.kept:
            continue
        want = ";".join(ents[k].rewritten(S2, E2) if k in ents else k for k in row.ids)
        assert {r["Protein.Ids"] for r in got[row.peptide]} == {want}
        for k in row.ids:
            e = ents.get(k)
            if e is not None and e.kind == "orf":
                assert parse_orf_id(e.rewritten(S2, E2)) == e.parsed(S2, E2)
                _, contig, s0, e0, _ = e.parsed(S2, E2)
                assert (contig, s0 + S2 + 1, e0 + S2) == (region_contig(S2, E2), e.s1, e.e1)
