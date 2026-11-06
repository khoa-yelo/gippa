"""M12-01 route E, review round 1 of #574: ``scripts/example/build_example.py`` hardening.

Spec docs/specs/M12-01-example-dataset.md (amendment 2026-11-06): rule 3E (kinds cds, sorf, orf
only), rule 8 and PROVENANCE.md (the A4 commands "with their paths shortened to base names", R12),
rule 9 (checks; "any other string column that holds one"; ``distinct_ids``), and the route E caveat
("the kept rows name 35 ORF entries (of 5,678 in the region)", lead fix e8aa8230). Per finding:

1. PROVENANCE path scrubbing: a file argument keeps its base name, a directory-only argument
   becomes ``<dir>``; no user name, ``/home``, ``/scratch`` or ``~/`` is left. The leak guard
   fails the build (exit 1) when the login name (``getpass.getuser()``) is in a command file.
2. A report column stored as ``dictionary<string>`` or ``string_view`` is rewritten or the build
   fails; it never exits 0 with chromosome ORF IDs left. ``distinct_ids`` = the distinct region
   ORF IDs in the written report.
3. An out-of-region ORF only in ``Protein.Group`` fails naming it; a region ORF only there is
   rewritten.
4. An ``nterm_extension`` map row with the wrong strand fails.
5. A record one residue longer than the translation fails.
6. A ``kind=orf`` row whose ID does not fit ``GIPPAORF_{class}_{contig}_{s1}_{e1}_{p|m}`` fails
   with a ``build_example:`` message naming it (no traceback).
7. Map kinds ``patch`` and ``contaminant`` are not region entries: rows naming them are dropped.
8. The README and the example page state the kept and region ORF-entry counts, not "a few
   dozen".

Synthetic A4 only (tests/helpers/example_a4.py); the user is monkeypatched to ``alice``. No
statistical claim, so no calibration test.
"""

from __future__ import annotations

import contextlib
import getpass
import importlib.util
import io
import json
import re
import tarfile
import traceback
from dataclasses import dataclass
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow.parquet as pq
import pytest
from helpers.example_a4 import (
    CHROM,
    REWRITTEN_KEYS,
    ROWS,
    Entry,
    Row,
    entries,
    make_pilot_e,
)
from helpers.example_m12_01 import REPO, flat

SCRIPT = REPO / "scripts" / "example" / "build_example.py"
S, E, L = 10_000, 40_000, 50_000
E_REPORT = "E_sixframe/report.parquet"
USER = "alice"
HOME = Path(f"/home/users/{USER}")
_CHROM_ORF = re.compile(rf"GIPPAORF_\w+?_{re.escape(CHROM)}_\d+_\d+_[pm]")

STEP1 = (
    f"diann --f /home/users/{USER}/x/report.parquet --temp /scratch/users/{USER}/ "
    f"--dir /home/users/{USER} --fasta=/a/b/db.faa --threads 8\n"
)
STEP2 = f"apptainer exec -B/home/{USER}/x:/data diann.sif --out ~/y --cfg ../z --qvalue 0.01\n"
EXTRAS = {
    "A4/command_step1.txt": STEP1,
    "A4/command_step2.txt": STEP2,
    "jobs/diann_A4.cfg": "--threads 8\n--qvalue 0.01\n",
    "db/makedb_manifest.json": '{"min_orf_aa": 20, "translation_table": 11}\n',
    "A4/diann_version.txt": "DIA-NN 2.3.1 Academia\n",
}


def _minor() -> str:
    import gippa

    m = re.match(r"(\d+)\.(\d+)", str(gippa.__version__))
    assert m is not None
    return f"{m.group(1)}.{m.group(2)}"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("gippa_build_example_route_e_r1", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@dataclass(frozen=True)
class Built:
    code: int
    out: str
    tree: Path | None
    archive: Path

    def text(self, name: str) -> str:
        assert self.tree is not None, self.out
        return (self.tree / name).read_text(encoding="utf-8")

    def sel(self) -> dict[str, Any]:
        return json.loads(self.text("SELECTION.json"))  # type: ignore[no-any-return]

    def report(self) -> list[dict[str, Any]]:
        assert self.tree is not None, self.out
        return pq.read_table(self.tree / E_REPORT).to_pylist()  # type: ignore[no-any-return]


def build(tmp: Path, **a4: Any) -> Built:
    """A small-region pilot with ``write_a4(**a4)``, built in-process. An uncaught exception is
    exit 99 with its traceback in ``out`` (finding 6 wants a ``build_example:`` message)."""
    root = (tmp / "gippa_root").resolve()
    make_pilot_e(root, S, E, L, **a4)
    rel = tmp / "rel"
    argv = ["--start", str(S), "--end", str(E), "--out", str(tmp / "work"), "--bundle-out",
            str(rel), "--allow-dirty"]  # fmt: skip
    buf = io.StringIO()
    with (
        pytest.MonkeyPatch.context() as mp,
        contextlib.redirect_stdout(buf),
        contextlib.redirect_stderr(buf),
    ):
        mp.setenv("GIPPA_ROOT", str(root))
        mp.delenv("SLURM_JOB_ID", raising=False)
        try:
            code = _load().main(argv)
        except SystemExit as exc:
            code = exc.code if isinstance(exc.code, int) else 1
        except Exception:  # reported; the test asserts the exit code
            buf.write(traceback.format_exc())
            code = 99
    archive = rel / f"gippa-example-{_minor()}.tar.gz"
    tree = None
    if archive.is_file():
        with tarfile.open(archive) as tar:
            tar.extractall(tmp / "unpacked", filter="data")
        tree = tmp / "unpacked" / f"gippa-example-{_minor()}"
    return Built(int(code or 0), buf.getvalue(), tree, archive)


@pytest.fixture
def as_alice(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(getpass, "getuser", lambda: USER)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: HOME))


def _fails_naming(b: Built, *needles: str) -> None:
    assert b.code == 1, b.out
    assert "build_example:" in b.out and "Traceback" not in b.out, b.out
    assert any(n in b.out for n in needles), (needles, b.out)
    assert not b.archive.exists()


def _orf_ids(rows: list[dict[str, Any]]) -> set[str]:
    return {p for r in rows for v in r.values() if isinstance(v, str) for p in v.split(";")
            if p.startswith("GIPPAORF_")}  # fmt: skip


# ---------------------------------------------------------------- 1. PROVENANCE scrubbing


@pytest.fixture
def scrubbed(tmp_path: Path, as_alice: None) -> Built:
    b = build(tmp_path, extras=EXTRAS)
    assert b.code == 0, b.out
    return b


@pytest.mark.parametrize("needle", [USER, "/home", "/scratch", "~/"])
def test_provenance_holds_no_personal_path(scrubbed: Built, needle: str) -> None:
    assert needle not in scrubbed.text("PROVENANCE.md")


def test_provenance_commands_keep_file_base_names(scrubbed: Built) -> None:
    text = flat(scrubbed.text("PROVENANCE.md"))
    assert re.search(r"(?<![\w.])report\.parquet", text), "the --f file keeps its base name"
    assert re.search(r"(?<![\w.])db\.faa", text), "the --fasta= file keeps its base name"
    assert "--threads 8" in text and "--qvalue 0.01" in text, "other arguments are kept"


def test_provenance_directory_arguments_become_a_placeholder(scrubbed: Built) -> None:
    text = flat(scrubbed.text("PROVENANCE.md"))
    assert "--temp <dir>" in text
    assert "--dir <dir>" in text


def test_leak_guard_fails_on_the_login_name(tmp_path: Path, as_alice: None) -> None:
    extras = {**EXTRAS, "A4/command_step2.txt": f"diann --label {USER} --qvalue 0.01\n"}
    b = build(tmp_path, extras=extras)
    assert b.code == 1, b.out
    assert "build_example:" in b.out and "Traceback" not in b.out, b.out
    assert not b.archive.exists()


def test_leak_guard_control_other_login(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The same command file builds when the login name is someone else's."""
    monkeypatch.setattr(getpass, "getuser", lambda: "bob")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: Path("/home/users/bob")))
    extras = {**EXTRAS, "A4/command_step2.txt": "diann --label carol --qvalue 0.01\n"}
    b = build(tmp_path, extras=extras)
    assert b.code == 0, b.out


# ---------------------------------------------------------------- 2. no-op rewrite


@pytest.mark.parametrize(
    "types",
    [{"Protein.Ids": "dictionary"}, {"Protein.Group": "dictionary"},
     {"Protein.Ids": "string_view", "Protein.Names": "string_view"}],
    ids=["ids_dictionary", "group_dictionary", "string_view"],
)  # fmt: skip
def test_non_plain_string_columns_are_rewritten_or_fail(
    tmp_path: Path, types: dict[str, str]
) -> None:
    b = build(tmp_path, column_types=types)
    assert b.code in (0, 1), b.out
    if b.code == 1:
        assert "build_example:" in b.out and not b.archive.exists(), b.out
        return
    rows = b.report()
    left = sorted({v for r in rows for v in r.values() if isinstance(v, str)
                   and _CHROM_ORF.search(v)})  # fmt: skip
    assert left == [], f"exit 0 with chromosome ORF IDs left: {left[:3]}"
    ents = entries(S, E)
    assert _orf_ids(rows) == {ents[k].rewritten(S, E) for k in REWRITTEN_KEYS}
    assert b.sel()["orf_id_rewrite"]["distinct_ids"] == len(REWRITTEN_KEYS)


def test_distinct_ids_counts_the_region_orf_ids_written(tmp_path: Path) -> None:
    b = build(tmp_path)
    assert b.code == 0, b.out
    assert b.sel()["orf_id_rewrite"]["distinct_ids"] == len(_orf_ids(b.report()))


# ---------------------------------------------------------------- 3. Protein.Group only


def _with_group_orf(ents: dict[str, Entry], key: str) -> tuple[Row, ...]:
    return (*ROWS, Row("GROUPONLYK", ("CDS_X",), 0, True, ("CDS_X", key)))


def test_out_of_region_orf_only_in_the_group_fails(tmp_path: Path) -> None:
    ents = entries(S, E)
    b = build(tmp_path, rows=_with_group_orf(ents, "O_FAR"))
    _fails_naming(b, ents["O_FAR"].db_id)


def test_region_orf_only_in_the_group_is_rewritten(tmp_path: Path) -> None:
    ents = entries(S, E)
    ents["O_GRP"] = Entry("O_GRP", "orf", S + 20_001, S + 20_300, "+", cls="intergenic")
    b = build(tmp_path, ents=ents, rows=_with_group_orf(ents, "O_GRP"))
    assert b.code == 0, b.out
    got = {r["Protein.Group"] for r in b.report() if r["Stripped.Sequence"] == "GROUPONLYK"}
    assert got == {f"CDS_X;{ents['O_GRP'].rewritten(S, E)}"}
    assert b.sel()["orf_id_rewrite"]["distinct_ids"] == len(REWRITTEN_KEYS) + 1


# ---------------------------------------------------------------- 4-6. rule 9 checks


def test_nterm_extension_with_the_wrong_map_strand_fails(tmp_path: Path) -> None:
    e = entries(S, E)["O_NT"]
    _fails_naming(build(tmp_path, map_strand="O_NT"), e.db_id, e.rewritten(S, E))


def test_record_one_residue_longer_than_the_translation_fails(tmp_path: Path) -> None:
    e = entries(S, E)["O_OUT"]
    _fails_naming(build(tmp_path, long_record="O_OUT"), e.db_id, e.rewritten(S, E))


@pytest.mark.parametrize(
    "bad_id",
    [f"GIPPAORF_intergenic_{CHROM}_x_{S + 800}_p", f"GIPPAORF_intergenic_{CHROM}_{S + 501}_p",
     f"GIPPAORF_intergenic_{CHROM}_{S + 501}_{S + 800}_q"],
    ids=["non_numeric", "missing_end", "bad_strand"],
)  # fmt: skip
def test_malformed_orf_id_fails_with_a_message(tmp_path: Path, bad_id: str) -> None:
    ents = entries(S, E)
    ents["O_BAD"] = Entry("O_BAD", "orf", S + 501, S + 800, "+", cls="intergenic", name=bad_id)
    rows = (*ROWS, Row("BADIDK", ("O_BAD",), 0, True))
    _fails_naming(build(tmp_path, ents=ents, rows=rows), bad_id)


# ---------------------------------------------------------------- 7. kinds patch, contaminant


@pytest.mark.parametrize(("kind", "name"), [("patch", "PATCH_1"), ("contaminant", "CONT_P9")])
def test_patch_and_contaminant_entries_are_not_region_entries(
    tmp_path: Path, kind: str, name: str
) -> None:
    ents = entries(S, E)
    ents["K"] = Entry("K", kind, S + 8_001, S + 8_300, "+", name=name)
    rows = (*ROWS, Row("KINDTGTK", ("K",), 0, False), Row("KINDMIXK", ("CDS_X", "K"), 0, False))
    b = build(tmp_path, ents=ents, rows=rows)
    assert b.code == 0, b.out
    peps = {r["Stripped.Sequence"] for r in b.report()}
    assert not peps & {"KINDTGTK", "KINDMIXK"}
    assert peps == {r.peptide for r in ROWS if r.kept}


# ---------------------------------------------------------------- 8. the ORF-entry caveat


def test_readme_counts_the_kept_and_region_orf_entries(tmp_path: Path) -> None:
    """Five region ORF entries are named by kept rows; three more lie in the region unnamed."""
    ents = entries(S, E)
    for i in range(3):
        a = S + 22_001 + 600 * i
        ents[f"O_U{i}"] = Entry(f"O_U{i}", "orf", a, a + 299, "-+"[i % 2], cls="antisense")
    b = build(tmp_path, ents=ents)
    assert b.code == 0, b.out
    text = flat(b.text("README.md"))
    n_region = sum(1 for e in ents.values() if e.kind == "orf" and e.in_region(S, E))
    assert n_region == len(REWRITTEN_KEYS) + 3
    assert f"{len(REWRITTEN_KEYS)} ORF entries (of {n_region:,} in the region)" in text, text
    assert "few dozen" not in text


def test_example_page_counts_the_kept_and_region_orf_entries() -> None:
    text = flat((REPO / "docs" / "get-started" / "example.md").read_text(encoding="utf-8"))
    assert "35 ORF entries (of 5,678 in the region)" in text
    assert "few dozen" not in text
