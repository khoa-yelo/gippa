"""M12-01 amendment 2026-11-06 (b): route E gets its own region (work-split row 11a, red).

Spec docs/specs/M12-01-example-dataset.md, "Amendment 2026-11-06 (b)": route E reads
``NZ_CP058217.1:1400001-1600000`` (``[1400000, 1600000)``), written as the linear contig
``nissle1917_NZ_CP058217_1_1400000_1600000`` to ``genome/nissle_r200k_e.fna.gz`` (rule 1E); rules
3E, 7 (route E's part) and 9 use that region, with ``S = S_E = 1400000``; routes A-D do not change,
byte for byte (AC 25); AC 26; the ``region_e`` keys (Q-L16, accepted); route E's commands use the
``_e`` genome (AC 14 as amended); ``measure()`` gives ``n_novel_unpaired = 0`` when the FDR table is
empty or absent; AC 22 as amended; AC 24 (the tutorial compares A-D only, accepted); the README,
example page and tutorial wording. Lead decisions: Q-L14 (genome FASTA only, no annotation file),
Q-L16 (the name ``nissle_r200k_e.fna.gz`` and the key ``region_e``), Q-L17.

Seams pinned here (the spec leaves them to the implementer):

- ``build_example.START_E`` / ``END_E``: route E's region, module constants read at call time
  (``--start``/``--end`` move the shared region only); the AC 25 builds set them on the loaded
  module (``helpers.example_a4.set_route_e``);
- ``SELECTION.json > region_e``: the fields of ``region`` for route E's region;
  ``bundle.json > region_e``: as ``bundle.json > region`` (the 1-based closed string), or the
  ``SELECTION.json`` object;
- ``SELECTION.json > orf_id_rewrite``: ``contig`` = route E's contig, ``offset`` = ``S_E``.

Synthetic data only (tests/helpers/example_a4.py; R11). The build makes no statistical claim, so
there is no calibration test.
"""

from __future__ import annotations

import contextlib
import gzip
import hashlib
import importlib.util
import io
import json
import re
import tarfile
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd
import pyarrow.parquet as pq
import pytest
from helpers import example_routes as er
from helpers.example_a4 import (
    CHROM,
    END_E,
    REPORT_REL,
    REWRITTEN_KEYS,
    ROWS,
    RUNS,
    SHARED_ROWS,
    START_E,
    entries,
    make_pilot_e,
    region_contig,
    set_route_e,
)
from helpers.example_m12_01 import REPO, flat, load_routes, options
from helpers.example_route_e import route_e_misses

SCRIPT = REPO / "scripts" / "example" / "build_example.py"
START, END = 1_965_000, 2_165_000  # the shared region of routes A-D (rule 1)
NAME = "nissle_r200k"
C_E = "nissle1917_NZ_CP058217_1_1400000_1600000"  # route E's contig (rule 1E)
C_SHARED = "nissle1917_NZ_CP058217_1_1965000_2165000"
GENOME, GENOME_E = f"genome/{NAME}.fna.gz", f"genome/{NAME}_e.fna.gz"
E_REPORT = "E_sixframe/report.parquet"
A_TO_D = ("A_diann", "B_proteins", "C_cascadia", "D_generic")
# the small layout of the AC 25 builds: shared region, two places for route E, chromosome length
S2, E2 = 10_000, 40_000
E_AT_1, E_AT_2, L2 = (60_000, 90_000), (110_000, 140_000), 150_000


def _minor() -> str:
    import gippa

    m = re.match(r"(\d+)\.(\d+)", str(gippa.__version__))
    assert m is not None
    return f"{m.group(1)}.{m.group(2)}"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gippa_build_example_route_e_region", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@dataclass(frozen=True)
class Built:
    root: Path
    chrom: str
    code: int
    out: str
    work: Path
    tree: Path | None

    def path(self, rel: str) -> Path:
        assert self.tree is not None, self.out
        return self.tree / rel

    def text(self, rel: str) -> str:
        return self.path(rel).read_text(encoding="utf-8")

    def json(self, rel: str) -> dict[str, Any]:
        return json.loads(self.text(rel))  # type: ignore[no-any-return]

    def report(self) -> list[dict[str, Any]]:
        return pq.read_table(self.path(E_REPORT)).to_pylist()  # type: ignore[no-any-return]

    def by_peptide(self) -> dict[str, list[dict[str, Any]]]:
        out: dict[str, list[dict[str, Any]]] = {}
        for r in self.report():
            out.setdefault(str(r["Stripped.Sequence"]), []).append(r)
        return out


def _build(
    tmp: Path,
    shared: tuple[int, int],
    e_region: tuple[int, int],
    length: int,
    *,
    region_args: bool,
    set_e: bool,
) -> Built:
    """A pilot with routes A-D in ``shared`` and A4 in ``e_region``, built in-process."""
    root = (tmp / "gippa_root").resolve()
    chrom, _ = make_pilot_e(root, *shared, length, e_region=e_region)
    work, rel = tmp / "work", tmp / "rel"
    argv = ["--out", str(work), "--bundle-out", str(rel), "--allow-dirty"]
    if region_args:
        argv += ["--start", str(shared[0]), "--end", str(shared[1])]
    mod, buf = _load(), io.StringIO()
    if set_e:
        set_route_e(mod, *e_region)
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
    return Built(root, chrom, int(code or 0), buf.getvalue(), work, tree)


@pytest.fixture(scope="module")
def default(tmp_path_factory: pytest.TempPathFactory) -> Built:
    """Every default: routes A-D in the shared region, A4 in route E's own region; no seam set."""
    tmp = tmp_path_factory.mktemp("route_e_region")
    b = _build(tmp, (START, END), (START_E, END_E), END + 5_000, region_args=False, set_e=False)
    assert b.code == 0, b.out
    return b


def _small(tmp: Path, e_region: tuple[int, int]) -> Built:
    return _build(tmp, (S2, E2), e_region, L2, region_args=True, set_e=True)


@pytest.fixture(scope="module")
def moved(tmp_path_factory: pytest.TempPathFactory) -> tuple[Built, Built]:
    """AC 25: the same routes A-D pilot, route E's start at two places."""
    one = _small(tmp_path_factory.mktemp("e_at_1"), E_AT_1)
    two = _small(tmp_path_factory.mktemp("e_at_2"), E_AT_2)
    assert one.chrom[S2:E2] == two.chrom[S2:E2], "fixture drift: the shared region differs"
    assert one.code == 0, one.out
    assert two.code == 0, two.out
    return one, two


def _fasta(path: Path) -> list[tuple[str, str]]:
    recs: list[tuple[str, str]] = []
    for line in gzip.decompress(path.read_bytes()).decode().splitlines():
        if line.startswith(">"):
            recs.append((line[1:].split()[0], ""))
        elif recs:
            recs[-1] = (recs[-1][0], recs[-1][1] + line.strip())
    return recs


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


ENT = entries(START_E, END_E)


# ---------------------------------------------------------------- rule 1E, AC 26: the region


def test_route_e_region_constants() -> None:
    mod = _load()
    assert (mod.START_E, mod.END_E) == (1_400_000, 1_600_000)
    assert (mod.START, mod.END, mod.NAME) == (START, END, NAME), "the shared region is unchanged"


def test_fixture_contig_is_the_spec_contig() -> None:
    assert region_contig(START_E, END_E) == C_E


def test_route_e_genome_is_its_own_region(default: Built) -> None:
    """AC 26: one record, ID ``nissle1917_NZ_CP058217_1_1400000_1600000``, equal to the source
    chromosome's ``[1400000, 1600000)``."""
    path = default.path(GENOME_E)
    assert path.is_file(), f"{GENOME_E} missing from the bundle"
    recs = _fasta(path)
    assert [r[0] for r in recs] == [C_E]
    assert recs[0][1] == default.chrom[START_E:END_E].upper()


def test_shared_genome_is_unchanged(default: Built) -> None:
    recs = _fasta(default.path(GENOME))
    assert [r[0] for r in recs] == [C_SHARED]
    assert recs[0][1] == default.chrom[START:END].upper()


def test_bundle_ships_no_annotation_file(default: Built) -> None:
    """Q-L14 (lead): route E ships the genome FASTA only."""
    assert default.tree is not None
    genome = sorted(p.name for p in (default.tree / "genome").iterdir())
    assert genome == [f"{NAME}.fna.gz", f"{NAME}_e.fna.gz"]
    assert sorted(p.name for p in (default.tree / "E_sixframe").iterdir()) == ["report.parquet"]


# ---------------------------------------------------------------- the region_e keys (Q-L16)

REGION_E = {
    "contig": CHROM,
    "start_0based": 1_400_000,
    "end_exclusive": 1_600_000,
    "length_bp": 200_000,
    "one_based_closed": "NZ_CP058217.1:1400001-1600000",
    "output_contig_id": C_E,
    "topology": "linear (no tag)",
}


def test_selection_region_e(default: Built) -> None:
    sel = default.json("SELECTION.json")
    assert sel.get("region_e") == REGION_E
    assert sel["region"]["start_0based"] == START and sel["region"]["output_contig_id"] == C_SHARED
    assert set(sel["region_e"]) == set(sel["region"]), "same fields as region"


def test_bundle_json_region_e(default: Built) -> None:
    meta = default.json("bundle.json")
    assert meta["region"] == f"{CHROM}:{START + 1}-{END}", "the shared region is unchanged"
    assert "region_e" in meta, sorted(meta)
    assert meta["region_e"] in (REGION_E["one_based_closed"], REGION_E), meta["region_e"]


def test_orf_id_rewrite_names_route_e_contig_and_offset(default: Built) -> None:
    rw = default.json("SELECTION.json")["orf_id_rewrite"]
    assert (rw["contig"], rw["offset"]) == (C_E, START_E)
    assert rw["distinct_ids"] == len(REWRITTEN_KEYS)


def test_route_e_genome_is_listed_and_checksummed(default: Built) -> None:
    outputs = default.json("SELECTION.json")["outputs"]
    assert GENOME_E in outputs, sorted(outputs)
    assert outputs[GENOME_E]["sha256"] == _sha(default.path(GENOME_E))
    sums = dict(
        reversed(ln.split(maxsplit=1)) for ln in default.text("CHECKSUMS.sha256").splitlines()
    )
    assert sums.get(GENOME_E) == _sha(default.path(GENOME_E)), "CHECKSUMS.sha256 lacks it"


# ---------------------------------------------------------------- rule 9 with S = S_E

# The spec's rule 9 examples, literally (S = 1400000, C = route E's contig): (input Protein.Ids,
# the row's peptide in the fixture, expected Protein.Ids or None when dropped).
SPEC_EXAMPLES = [
    ("CDS_X", "CDSTARGETK", "CDS_X"),
    ("GIPPAORF_antisense_NZ_CP058217.1_1400101_1400400_m", "ASDECOYK",
     f"GIPPAORF_antisense_{C_E}_101_400_m"),
    ("CDS_X;GIPPAORF_cds_outframe_NZ_CP058217.1_1435002_1435301_p", "MIXTARGETK",
     f"CDS_X;GIPPAORF_cds_outframe_{C_E}_35002_35301_p"),
    ("GIPPAORF_intergenic_NZ_CP058217.1_1599900_1600120_p", "XENDK", None),
    ("CDS_X;GIPPAORF_antisense_NZ_CP058217.1_10_309_m", "FARORFK", None),
    ("CDS_Y;GIPPAORF_antisense_NZ_CP058217.1_1400101_1400400_m", "SHAREDCDSK", None),
    ("GIPPAORF_antisense_NZ_CP058217.1_1965101_1965400_m", "SHAREDORFK", None),
    ("Cont_P1;GIPPAORF_antisense_NZ_CP058217.1_1400101_1400400_m", "CONTORFK", None),
    ("GIPPAORF_intergenic_NZ_CP058217.1_1400001_1400300_p", "FIRSTBASEK",
     f"GIPPAORF_intergenic_{C_E}_1_300_p"),
    ("GIPPAORF_intergenic_NZ_CP058217.1_1599701_1600000_p", "LASTBASEK",
     f"GIPPAORF_intergenic_{C_E}_199701_200000_p"),
]  # fmt: skip


@pytest.mark.parametrize(("ids", "peptide", "want"), SPEC_EXAMPLES,
                         ids=[x[1] for x in SPEC_EXAMPLES])  # fmt: skip
def test_spec_rule_9_examples_at_route_e_region(
    default: Built, ids: str, peptide: str, want: str | None
) -> None:
    src = {r["Protein.Ids"] for r in pq.read_table(default.root / REPORT_REL).to_pylist()
           if r["Stripped.Sequence"] == peptide}  # fmt: skip
    assert src == {ids}, "fixture drift: the row does not carry the spec's input"
    got = default.by_peptide().get(peptide)
    if want is None:
        assert got is None, f"{ids} must be dropped"
    else:
        assert got is not None, f"{ids} must be kept"
        assert {r["Protein.Ids"] for r in got} == {want}


def test_spec_parse_example_at_route_e_region() -> None:
    from gippa.novel.fdr import parse_orf_id

    assert parse_orf_id(f"GIPPAORF_antisense_{C_E}_101_400_m") == ("antisense", C_E, 100, 400, "-")


def test_every_orf_id_names_route_e_contig_with_offset_s_e(default: Built) -> None:
    """AC 26: every rewritten ORF ID names route E's contig and uses ``S = S_E``; read back by
    ``parse_orf_id`` it is the map interval minus 1,400,000; nothing names the shared contig."""
    from gippa.novel.fdr import parse_orf_id

    want = {ENT[k].rewritten(START_E, END_E): ENT[k] for k in REWRITTEN_KEYS}
    seen = set()
    for r in default.report():
        for col in ("Protein.Ids", "Protein.Group", "Protein.Names"):
            for pid in str(r[col]).split(";"):
                if not pid.startswith("GIPPAORF_"):
                    continue
                assert pid in want, (col, pid)
                cls, contig, s0, e0, strand = parse_orf_id(pid)
                e = want[pid]
                assert (cls, contig, strand) == (e.cls, C_E, e.strand), pid
                assert (s0 + START_E + 1, e0 + START_E) == (e.s1, e.e1), pid
                seen.add(pid)
        for v in r.values():
            if isinstance(v, str):
                assert C_SHARED not in v and f"_{CHROM}_" not in v, v
    assert seen == set(want)


# ---------------------------------------------------------------- rule 3E on route E's region


@pytest.mark.parametrize("row", SHARED_ROWS, ids=lambda r: r.peptide)
def test_rule_3e_drops_shared_region_entries(default: Built, row: Any) -> None:
    """AC 26: no kept row names an entry outside route E's region, even one in the shared region."""
    assert row.peptide not in default.by_peptide(), row.ids


def test_route_e_counts_on_its_region(default: Built) -> None:
    e = default.json("SELECTION.json")["counts"]["E"]
    kept = [r for r in ROWS if r.kept]
    assert e["rows_total"] == (len(ROWS) + len(SHARED_ROWS)) * len(RUNS)
    assert e["rows_kept"] == len(kept) * len(RUNS)
    assert e["decoy_rows_kept"] == sum(r.decoy for r in kept) * len(RUNS)
    n_orf = sum(1 for k, v in ENT.items() if v.kind == "orf" and v.in_region(START_E, END_E))
    assert e["region_orf_entries"] == n_orf


# ---------------------------------------------------------------- AC 25: A-D byte for byte


def _a_to_d_files(b: Built) -> dict[str, str]:
    assert b.tree is not None, b.out
    files = {p.relative_to(b.tree).as_posix(): _sha(p) for d in A_TO_D
             for p in sorted((b.tree / d).iterdir())}  # fmt: skip
    files[GENOME] = _sha(b.tree / GENOME)
    return files


def test_a_to_d_byte_identical_when_route_e_moves(moved: tuple[Built, Built]) -> None:
    one, two = moved
    a, b = _a_to_d_files(one), _a_to_d_files(two)
    assert len(a) == 8 and a == b, {k: (a.get(k), b.get(k)) for k in a.keys() | b.keys()
                                    if a.get(k) != b.get(k)}  # fmt: skip


def test_a_to_d_counts_unchanged_when_route_e_moves(moved: tuple[Built, Built]) -> None:
    one, two = (m.json("SELECTION.json") for m in moved)
    for key in ("A", "B", "C", "D", "region_proteins"):
        assert one["counts"][key] == two["counts"][key], key
    assert one["region"] == two["region"]


def test_route_e_genome_follows_its_start(moved: tuple[Built, Built]) -> None:
    """AC 26 for any ``S_E``: ``genome/<name>_e.fna.gz`` is ``[S_E, END_E)`` of the chromosome,
    named after that interval, and the rewrite's contig and offset follow it."""
    for b, (s, e) in zip(moved, (E_AT_1, E_AT_2), strict=True):
        assert b.json("SELECTION.json")["region_e"]["start_0based"] == s
        assert _fasta(b.path(GENOME_E)) == [(region_contig(s, e), b.chrom[s:e].upper())]
        rw = b.json("SELECTION.json")["orf_id_rewrite"]
        assert (rw["contig"], rw["offset"]) == (region_contig(s, e), s)
        assert {r["Stripped.Sequence"] for r in b.report()} == {r.peptide for r in ROWS if r.kept}


# ---------------------------------------------------------------- README and PROVENANCE


def test_readme_route_e_region_wording(default: Built) -> None:
    text = flat(default.text("README.md"))
    for needle in (GENOME, GENOME_E, "add 1,400,000", "route E has its own"):
        assert needle in text, needle
    assert "add 1,965,000" not in text
    assert re.search(r"routes A(\u2013|-| to )D share one region", text), "A-D share one region"
    assert "The same region is given through five input routes" not in text


def test_readme_intro_names_route_e_genome(default: Built) -> None:
    text = flat(default.text("README.md"))
    assert f"`NZ_CP058217.1:{START_E + 1}-{END_E}`, in `{GENOME_E}`" in text


def test_readme_counts_orf_entries_of_route_e_region(default: Built) -> None:
    e = default.json("SELECTION.json")["counts"]["E"]
    named, in_region = e["orf_entries_named"], e["region_orf_entries"]
    want = f"{named:,} ORF entries (of {in_region:,} in route E's region)"
    assert want in flat(default.text("README.md"))


def test_readme_explains_why_route_e_has_its_own_region(default: Built) -> None:
    text = flat(default.text("README.md"))
    assert "own 200 kb region" in text and "global q-value" in text


def test_provenance_names_route_e_region(default: Built) -> None:
    text = flat(default.text("PROVENANCE.md"))
    for needle in ("NZ_CP058217.1:1400001-1600000", "[1400000, 1600000)", C_E,
                   _sha(default.path(GENOME_E))):  # fmt: skip
        assert needle in text, needle


# ---------------------------------------------------------------- AC 14: the commands


def test_route_e_commands_use_the_e_genome() -> None:
    cmds = load_routes()["E"]["commands"]
    for c in cmds[1:]:
        assert options(c)["--genome"] == ["{dir}/genome/nissle_r200k_e.fna.gz"], c


@pytest.mark.parametrize("route", "ABCD")
def test_a_to_d_keep_the_shared_genome(route: str) -> None:
    for c in load_routes()[route]["commands"]:
        if "--genome" in options(c):
            assert options(c)["--genome"] == ["{dir}/genome/nissle_r200k.fna.gz"], c


def test_build_writes_every_genome_the_routes_name(default: Built) -> None:
    named = {p for r in load_routes().values() for c in r["commands"]
             for p in re.findall(r"\{dir\}/(genome/\S+)", c)}  # fmt: skip
    assert named == {GENOME, GENOME_E}
    assert all(default.path(p).is_file() for p in named)


# ---------------------------------------------------------------- measure(): unpaired 0


def _run_dir(tmp: Path, *, fdr: str) -> Path:
    (tmp / "tables").mkdir(parents=True)
    (tmp / "novel").mkdir()
    man = {"status": "completed", "detection_tier": "mayu", "layers_missing": [], "warnings": [],
           "novel_orfs": {"classes": []}}  # fmt: skip
    (tmp / "run_manifest.json").write_text(json.dumps(man))
    (tmp / "tables" / "novel_orfs.tsv").write_text("# gippa novel_orfs\norf_id\tcalled\tclass\n")
    if fdr == "empty":
        pd.DataFrame({"orf_id": pd.Series([], dtype=str), "reason": pd.Series([], dtype=str)}
                     ).to_parquet(tmp / "novel" / "novel_class_fdr.parquet")  # fmt: skip
    return tmp


@pytest.mark.parametrize("fdr", ["empty", "absent"])
def test_measure_unpaired_is_zero_without_fdr_rows(tmp_path: Path, fdr: str) -> None:
    got = er.measure("E", _run_dir(tmp_path, fdr=fdr))
    assert got["n_novel_unpaired"] == 0 and type(got["n_novel_unpaired"]) is int
    assert got["n_novel_candidates"] == 0


# ---------------------------------------------------------------- AC 22 (amended)

_OK = {"detection_tier": "mayu", "n_novel_called": 0, "n_novel_candidates": 7,
       "n_novel_unpaired": 0, "novel_classes": [{"regime": "picked_decoy"}]}  # fmt: skip


@pytest.mark.parametrize(
    ("change", "field"),
    [({}, None),
     ({"n_novel_candidates": 0, "novel_classes": []}, "n_novel_candidates"),
     ({"n_novel_unpaired": None}, "n_novel_unpaired"),
     ({"n_novel_unpaired": 7}, "n_novel_unpaired"),
     ({"n_novel_unpaired": 0.0}, "n_novel_unpaired"),
     ({"novel_classes": []}, "novel_classes")],
)  # fmt: skip
def test_ac22_route_e_misses(change: dict[str, Any], field: str | None) -> None:
    assert [m.split(":")[0] for m in route_e_misses(_OK | change)] == ([field] if field else [])


# ---------------------------------------------------------------- docs

DOCS = REPO / "docs"
DIR, OUT = "gippa-example", "gippa-example/runs"
TITLE_E = "Route E: DIA-NN report, six-frame database"
FENCE = re.compile(r"^```[^\n]*\n(.*?)^```", re.S | re.M)


def _doc(rel: str) -> str:
    return (DOCS / rel).read_text(encoding="utf-8")


def _fenced(text: str) -> list[str]:
    return [" ".join(ln.split()) for b in FENCE.findall(text)
            for ln in b.replace("\\\n", " ").splitlines() if ln.startswith("gippa ")]  # fmt: skip


def _section(text: str, heading: str) -> str:
    m = re.search(rf"^{re.escape(heading)}[ \t]*$", text, re.M)
    assert m, heading
    nxt = re.compile(r"^#{1,3} ", re.M).search(text, m.end())
    return text[m.end() : nxt.start() if nxt else None]


def test_page_route_e_section_wording() -> None:
    body = _section(_doc("get-started/example.md"), f"### {TITLE_E}")
    cmds = [c for c in _fenced(body) if c.split()[1] in ("preflight", "run")]
    assert len(cmds) == 2
    assert all(f"--genome {DIR}/genome/nissle_r200k_e.fna.gz " in c for c in cmds), cmds
    text = flat(body)
    for needle in ("add 1,400,000", "53 ORF entries (of 5,609 in route E's region)",
                   "own 200 kb region", "global q-value"):  # fmt: skip
        assert needle.lower() in text.lower(), needle
    assert "1,965,000" not in text and "35 ORF entries" not in text


def test_page_routes_a_to_d_keep_the_shared_genome() -> None:
    cmds = [c for c in _fenced(_doc("get-started/example.md")) if "--genome" in c]
    shared = [c for c in cmds if f"{OUT}/E/" not in c]
    assert len(shared) == 8, shared
    assert all(f"--genome {DIR}/genome/nissle_r200k.fna.gz " in c for c in shared)


def test_page_links_the_tutorial_by_its_new_title() -> None:
    text = flat(_doc("get-started/example.md"))
    assert "Nissle regions, five routes](../tutorials/example-routes.md)" in text
    assert "One Nissle region, five routes" not in text


def test_tutorial_title_and_intro() -> None:
    text = _doc("tutorials/example-routes.md")
    assert text.splitlines()[0] == "# Tutorial: Nissle regions, five routes"
    intro = flat(text.split("\n## ", 1)[0])
    assert re.search(r"route E runs on its own (200 kb )?region", intro, re.I), intro


def test_tutorial_route_e_commands_use_the_e_genome() -> None:
    cmds = [c for c in _fenced(_doc("tutorials/example-routes.md"))
            if f"{OUT}/E/" in c and c.split()[1] in ("preflight", "run")]  # fmt: skip
    assert len(cmds) == 2 and all("genome/nissle_r200k_e.fna.gz" in c for c in cmds), cmds


def test_tutorial_orf_entry_sentence() -> None:
    text = flat(_doc("tutorials/example-routes.md"))
    assert "53 ORF entries (of 5,609 in route E's region)" in text
    assert "5,678" not in text


def test_tutorial_compares_a_to_d_only() -> None:
    """AC 24 (lead: accepted): exactly one ``gippa compare``, runs A-D, labels A-D, ``--outdir``,
    no ``{out}/E``; no row-level content."""
    text = _doc("tutorials/example-routes.md")
    cmds = [c for c in _fenced(text) if c.split()[1] == "compare"]
    assert len(cmds) == 1, cmds
    words = cmds[0].split()
    assert words[2:6] == [f"{OUT}/{r}" for r in "ABCD"], cmds[0]
    assert words[words.index("--labels") + 1] == "A,B,C,D"
    assert "--outdir" in words and f"{OUT}/E" not in cmds[0]
    assert "GIPPAORF_" not in text


# ---------------------------------------------------------------- sizing tars and guards (review)


def _tar_genomes(b: Built, route: str) -> set[str]:
    tar = b.work / NAME / "bundles" / f"{NAME}_{route}.tar.gz"
    with tarfile.open(tar) as t:
        return {m.split("/", 1)[-1] for m in t.getnames() if "/genome/" in f"/{m}"}


@pytest.mark.parametrize("route", ["A_diann", "B_proteins", "C_cascadia", "D_generic"])
def test_sizing_tar_a_to_d_hold_the_shared_genome_only(default: Built, route: str) -> None:
    assert _tar_genomes(default, route) == {GENOME}


def test_sizing_tar_e_holds_its_own_genome_only(default: Built) -> None:
    assert _tar_genomes(default, "E_sixframe") == {GENOME_E}


def _rerun(pilot: Built, tmp: Path, e_region: tuple[int, int]) -> tuple[int, str]:
    """``main`` on an existing small pilot, route E's seam set to ``e_region``."""
    argv = ["--out", str(tmp / "work"), "--bundle-out", str(tmp / "rel"), "--allow-dirty",
            "--start", str(S2), "--end", str(E2)]  # fmt: skip
    mod, buf = _load(), io.StringIO()
    set_route_e(mod, *e_region)
    with (
        pytest.MonkeyPatch.context() as mp,
        contextlib.redirect_stdout(buf),
        contextlib.redirect_stderr(buf),
    ):
        mp.setenv("GIPPA_ROOT", str(pilot.root))
        mp.delenv("SLURM_JOB_ID", raising=False)
        try:
            code = mod.main(argv)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
    return int(code or 0), buf.getvalue()


@pytest.mark.parametrize(
    ("e_region", "message"),
    [((90_000, 60_000), "is empty or negative"), ((120_000, L2 + 1), "is beyond")],
)
def test_route_e_region_guards(
    moved: tuple[Built, Built], tmp_path: Path, e_region: tuple[int, int], message: str
) -> None:
    code, out = _rerun(moved[0], tmp_path, e_region)
    assert code != 0, out
    assert message in out and "route E" in out, out
