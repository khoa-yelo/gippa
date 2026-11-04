"""M12-03 (GIPPA 1.1.0): the ``function`` module runs GO and COG when KEGG is not installed.

Khoa's decision (2026-11-04): KEGG is free for academic use only, so non-academic users install
with ``gippa setup --no-kegg``. GO (CC BY 4.0) and COG (NCBI) carry no such restriction, so
``function`` must still run them and skip only the KEGG parts.

Contract pinned here (test-author reading of the lead's brief):

1. **Preflight.** ``function`` on, GO and COG installed, KEGG not: no error. Exactly one warning,
   rule id :data:`RULE` (``environment.kegg``, beside ``environment.tools`` and
   ``environment.resources``; input ``environment``, exit code 0), whose message says KEGG is not
   installed, KEGG pathway analyses are skipped, and GO and COG still run. ``environment.tools``
   no longer lists ``kegg`` as missing. GO or COG missing is still an ``environment.tools`` error
   (exit 2). ``function`` off, or KEGG installed: no :data:`RULE` check at all. This holds for the
   injected probe (:func:`environment_checks`) and for the real probe (``run_preflight`` reading
   ``db_manifest.json``) under both ``mode="run"`` and ``mode="preflight"``; ``preflight.json``
   still validates against its schema.
2. **Workflow.** Pinned by workflow/tests/main_function_no_kegg.nf.test: ``FUNCTION_STEP`` runs
   once with the GO and COG items only (``functionResources()`` stages installed items only) and
   publishes the three tables.
3. **Outputs (one rule).** The three function tables keep their full header; they carry no
   ``collection == "kegg"`` row; every other table is untouched. ``enrichment_detection`` and
   ``enrichment_gsea`` equal the KEGG run's tables with the KEGG rows removed;
   ``category_summary`` equals them on every column but the bootstrap CI (``ci_low``,
   ``ci_high``: see the design question in the report). The step's run flags carry
   ``layer_missing:kegg``; ``run_manifest.json`` records ``layers_missing: [..., "kegg"]``,
   ``warnings: [..., "preflight:environment.kegg"]`` (the why), and no KEGG licence note.
4. **Report.** The function section says "KEGG not installed: pathway enrichment skipped", with
   no academic-use sentence and no KEGG retrieval date, also when only ``enrichment_gsea`` exists.
5. **No regression.** KEGG installed: the same environment checks as before, KEGG rows present,
   the KEGG caveat with its retrieval date in the report.
6. **``gippa setup --no-kegg``** never tells the user that ``function`` is unavailable; the docs
   no longer tell ``--no-kegg`` users to drop ``function``.
"""

from __future__ import annotations

import json
import logging
import re
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.config import GippaConfig, resolve_config
from gippa.exitcodes import ExitCode
from gippa.preflight.checks import Check
from gippa.preflight.environment import GIB, EnvProbe, environment_checks, requirements

REPO = Path(__file__).resolve().parents[2]
RULE = "environment.kegg"
FUNCTION_DBS = ("go", "kegg", "cog")
KEGG_ITEMS = ("kegg_ko_pathway", "kegg_pathway_names")
SETUP_MANIFEST = REPO / "tests" / "data" / "db_manifest" / "db_manifest.json"
PREFLIGHT_SCHEMA = REPO / "src" / "gippa" / "schemas" / "preflight.schema.json"
LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"
REPORT_NOTE = "KEGG not installed: pathway enrichment skipped"


# ================================================================ 1. preflight, injected probe


def _found(modules: Sequence[str], missing: Sequence[str] = ()) -> dict[str, str | None]:
    """Every requirement of ``modules`` (and every function DB) found, then ``missing`` absent."""
    reqs = requirements(list(modules), "targeted")
    found: dict[str, str | None] = {r.name: r.min_version or "1.0" for r in reqs}
    for name in FUNCTION_DBS:
        found.setdefault(name, "1.0")
    for name in missing:
        found.pop(name, None)
    return found


def _checks(modules: Sequence[str], missing: Sequence[str] = ()) -> dict[str, Check]:
    reqs = requirements(list(modules), "targeted")
    probe = EnvProbe(found=_found(modules, missing), free_disk_bytes=500 * GIB,
                     total_ram_bytes=64 * GIB)  # fmt: skip
    out = environment_checks(reqs, probe, 0, 0)
    ids = [c.rule_id for c in out]
    assert len(ids) == len(set(ids)), f"rule reported twice: {ids}"
    for c in out:
        json.dumps(c.value)
    return {c.rule_id: c for c in out}


def _assert_kegg_warning(c: Check) -> None:
    assert c.level == "warn", c
    assert c.input == "environment", c
    assert c.exit_code == 0, c
    assert c.rule, "empty human rule"
    msg = c.message.lower()
    assert "kegg" in msg and "not installed" in msg, c.message
    assert "skip" in msg, f"the warning must say KEGG analyses are skipped: {c.message!r}"
    assert re.search(r"\bgo\b", msg) and "cog" in msg, (
        f"the warning must say GO and COG still run: {c.message!r}"
    )
    assert "kegg" in json.dumps(c.value).lower(), c.value


def test_kegg_missing_with_function_is_one_warning_not_an_error() -> None:
    checks = _checks(["core", "function"], missing=["kegg"])
    assert RULE in checks, f"no {RULE} check when KEGG is missing: {sorted(checks)}"
    _assert_kegg_warning(checks[RULE])
    tools = checks["environment.tools"]
    assert tools.level == "ok", f"KEGG missing must not fail environment.tools: {tools}"
    assert tools.exit_code == 0
    assert [c for c in checks.values() if c.level == "error"] == []
    assert [c.rule_id for c in checks.values() if c.level == "warn"] == [RULE]


@pytest.mark.parametrize("db", ["go", "cog"])
def test_go_or_cog_missing_is_still_an_environment_error(db: str) -> None:
    checks = _checks(["core", "function"], missing=[db])
    tools = checks["environment.tools"]
    assert tools.level == "error" and tools.exit_code == 2, tools
    assert db in tools.value["missing"], tools.value  # type: ignore[index,call-overload]
    assert RULE not in checks, f"KEGG is installed, yet {RULE}: {checks.get(RULE)}"


def test_go_and_kegg_missing_error_on_go_and_warn_on_kegg() -> None:
    checks = _checks(["core", "function"], missing=["go", "kegg"])
    tools = checks["environment.tools"]
    assert tools.level == "error" and tools.exit_code == 2, tools
    missing = tools.value["missing"]  # type: ignore[index,call-overload]
    assert "go" in missing and "kegg" not in missing, missing
    assert RULE in checks, sorted(checks)
    _assert_kegg_warning(checks[RULE])


def test_kegg_missing_with_function_off_says_nothing() -> None:
    checks = _checks(["core"], missing=["kegg"])
    assert RULE not in checks, checks.get(RULE)
    assert checks["environment.tools"].level == "ok"


def test_kegg_installed_keeps_todays_two_checks() -> None:
    """Item 5: with KEGG installed nothing changes (``environment.tools`` ok, no new check)."""
    checks = _checks(["core", "mobile", "bgc", "function"])
    assert set(checks) == {"environment.tools", "environment.resources"}, sorted(checks)
    assert checks["environment.tools"].level == "ok"


@settings(max_examples=40, deadline=None)
@given(
    missing=st.sets(st.sampled_from(FUNCTION_DBS)),
    function_on=st.booleans(),
    others=st.sets(st.sampled_from(["mobile", "bgc"])),
)
def test_property_only_go_and_cog_are_required_kegg_only_warns(
    missing: set[str], function_on: bool, others: set[str]
) -> None:
    """For every subset of {go, kegg, cog} missing: an error iff function is on and GO or COG is
    missing; the KEGG warning iff function is on and KEGG is missing; ``kegg`` never in the
    error's ``missing`` list."""
    modules = ["core", *sorted(others), *(["function"] if function_on else [])]
    checks = _checks(modules, missing=sorted(missing))
    tools = checks["environment.tools"]
    want_error = function_on and bool(missing & {"go", "cog"})
    assert (tools.level == "error") is want_error, (modules, missing, tools)
    if tools.level == "error":
        assert "kegg" not in tools.value.get("missing", []), tools.value  # type: ignore[union-attr]
    want_warn = function_on and "kegg" in missing
    assert (RULE in checks) is want_warn, (modules, missing, sorted(checks))
    if want_warn:
        _assert_kegg_warning(checks[RULE])


# ================================================================ 1. preflight, real probe


def _db_without_kegg(root: Path, drop: Sequence[str] = KEGG_ITEMS) -> Path:
    """The real-format setup manifest (tests/data/db_manifest) without the ``drop`` items."""
    doc = json.loads(SETUP_MANIFEST.read_text(encoding="utf-8"))
    doc["items"] = [i for i in doc["items"] if i["name"] not in set(drop)]
    db = root / "db"
    db.mkdir(parents=True, exist_ok=True)
    (db / "db_manifest.json").write_text(json.dumps(doc), encoding="utf-8")
    return db


@pytest.fixture
def isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    repo = tmp_path / "repo"
    repo.mkdir()
    monkeypatch.setenv("GIPPA_REPO", str(repo))
    monkeypatch.setenv("GIPPA_ROOT", str(repo))
    monkeypatch.delenv("GIPPA_DB_DIR", raising=False)
    return tmp_path


def _inputs(root: Path) -> Any:
    from gippa.preflight.run import PreflightInputs

    g = root / "g.fa"
    g.write_text(">chr1 circular=true\n" + ("ACGGTCAT" * 18_750) + "\n", encoding="utf-8")
    p = root / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n>S1\nMKVLAAGIVGR\n", encoding="utf-8")
    pep = root / "peptides.tsv"
    pep.write_text("peptide\tprotein_ids\nPEPTIDEK\tP1\n", encoding="utf-8")
    return PreflightInputs(genome=g, proteins=p, peptides=pep)


def _cfg(modules: str = "[core,function]") -> GippaConfig:
    return resolve_config(overrides=[*LOW, f"modules={modules}"])


def _real_preflight(root: Path, mode: str, *, drop: Sequence[str] = KEGG_ITEMS,
                    modules: str = "[core,function]") -> Any:  # fmt: skip
    from gippa.preflight.run import run_preflight

    db = _db_without_kegg(root, drop)
    report, _ = run_preflight(_inputs(root), _cfg(modules), None, mode,  # type: ignore[arg-type]
                              root / "OUT", db)  # fmt: skip
    return report


@pytest.mark.parametrize("mode", ["run", "preflight"])
def test_real_probe_without_kegg_passes_with_one_kegg_warning(isolated: Path, mode: str) -> None:
    from gippa.preflight.run import exit_code

    report = _real_preflight(isolated, mode)
    env = {c.rule_id: c for c in report.checks if c.input == "environment"}
    assert env["environment.tools"].level in ("ok", "info"), env["environment.tools"]
    errors = [c for c in report.checks if c.level == "error"]
    assert errors == [], errors
    kegg = [c for c in report.checks if c.rule_id == RULE]
    assert len(kegg) == 1, f"expected one {RULE} check, got {[c.rule_id for c in report.checks]}"
    _assert_kegg_warning(kegg[0])
    assert exit_code(report, strict=False) == ExitCode.SUCCESS


def test_real_probe_with_kegg_has_no_kegg_warning(isolated: Path) -> None:
    report = _real_preflight(isolated, "run", drop=())
    assert not [c for c in report.checks if c.rule_id == RULE]
    assert not [c for c in report.checks if c.level == "error"], report.checks


@pytest.mark.parametrize("item", ["go_basic", "cog2020_fun"])
def test_real_probe_without_go_or_cog_still_exits_2(isolated: Path, item: str) -> None:
    from gippa.preflight.run import exit_code

    report = _real_preflight(isolated, "run", drop=(*KEGG_ITEMS, item))
    tools = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert len(tools) == 1 and tools[0].level == "error", tools
    assert exit_code(report, strict=False) == ExitCode.ENVIRONMENT


def test_real_probe_without_kegg_function_off_has_no_kegg_warning(isolated: Path) -> None:
    report = _real_preflight(isolated, "run", modules="[core]")
    assert not [c for c in report.checks if c.rule_id == RULE]


def test_preflight_json_with_the_kegg_warning_validates(isolated: Path) -> None:
    import jsonschema

    from gippa.preflight.report import to_json_dict

    doc = to_json_dict(_real_preflight(isolated, "run"))
    schema = json.loads(PREFLIGHT_SCHEMA.read_text(encoding="utf-8"))
    jsonschema.Draft202012Validator(schema).validate(doc)
    rows = [c for c in doc["checks"] if c["rule_id"] == RULE]  # type: ignore[index,union-attr]
    assert len(rows) == 1 and rows[0]["level"] == "warn", rows  # type: ignore[index]


# ================================================================ 3. outputs of the step


@pytest.fixture(autouse=True, scope="module")
def _simple_gsea_engine() -> Iterator[None]:
    """fgsea_simple, never R (ADR-0038 §4): unit tests use the in-process engine."""
    import gippa.steps.function as fn
    from gippa.function.gsea import SIMPLE_ENGINE

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(fn, "_LazyRFgsea", lambda log: SIMPLE_ENGINE)
        yield


TABLES = ("category_summary", "enrichment_detection", "enrichment_gsea")
CI_COLUMNS = ("ci_low", "ci_high")


def _step_cfg() -> GippaConfig:
    return GippaConfig.model_validate({"seed": 7, "stats": {"bootstrap": 100},
                                       "function": {"gsea": {"n_perm_simple": 200}}})  # fmt: skip


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


class _StepRun:
    def __init__(self, out: Path, records: list[logging.LogRecord]) -> None:
        self.out = out
        self.records = records

    def path(self, key: str) -> Path:
        return self.out / "tables" / f"{key}.tsv"

    def frame(self, key: str) -> pd.DataFrame:
        return pd.read_csv(self.path(key), sep="\t", dtype=str, keep_default_na=False)

    def header(self, key: str) -> str:
        return self.path(key).read_text(encoding="utf-8").split("\n", 1)[0]

    def run_flags(self) -> list[str]:
        found = [r.data for r in self.records  # type: ignore[attr-defined]
                 if getattr(r, "event", None) == "function.run_flags"]  # fmt: skip
        assert len(found) == 1, found
        return [str(f) for f in found[0]["flags"]]


class _World:
    """The ~600-gene world of tests/unit/steps/test_function_step.py (fitted detectability,
    synthetic GO/KEGG/COG DB from helpers.function_db), run with and without KEGG."""

    def __init__(self, root: Path) -> None:
        from helpers.function_db import write_function_db
        from helpers.function_step_data import (
            world_pre,
            write_gene_table,
            write_replicons,
            write_tu_folds,
        )

        from gippa.genetable.covariates import COVARIATE_COLUMNS
        from gippa.genetable.detectability import fit_detectability
        from gippa.io.genes_tsv import read_genes_parquet
        from gippa.rng import for_step
        from gippa.steps import load_step_modules
        from gippa.steps import registry as reg

        load_step_modules()
        self.root = root
        _, xrefs = world_pre()
        self.db = write_function_db(root / "db", xrefs)
        pre, _ = world_pre(root / "db")
        pre_path = write_gene_table(root / "genetable" / "genes_pre.parquet", pre)
        cfg = _step_cfg()
        self.gene_table = root / "genetable" / "genes.parquet"
        reg.run_step("detectability", cfg, {"genes_pre": pre_path},
                     {"gene_table": self.gene_table,
                      "detectability_summary": root / "genetable" / "det.parquet"})  # fmt: skip
        p = read_genes_parquet(pre_path)
        idx = pd.Index(p["locus_tag"].astype(str))
        det = pd.Series(p["detected"].astype(bool).to_numpy(), index=idx, dtype="boolean")
        cov = p[["locus_tag", *COVARIATE_COLUMNS]].astype(
            dict.fromkeys(COVARIATE_COLUMNS, "float64")
        )
        fit = fit_detectability(det, cov, pd.Series(p["tu_id"].to_numpy(), index=idx),
                                cfg.detectability, for_step("detectability", cfg.seed))  # fmt: skip
        self.tu_folds = write_tu_folds(root / "genetable" / "tu_folds.parquet",
                                       dict(fit.fold_of_tu))  # fmt: skip
        self.replicons = write_replicons(root / "replicons.parquet")

    def run(self, name: str, items: Sequence[str]) -> _StepRun:
        from gippa.steps import registry as reg

        ins = {"gene_table": self.gene_table, "replicons": self.replicons,
               "tu_folds": self.tu_folds} | self.db.inputs(items)  # fmt: skip
        out = self.root / name
        outputs = {k: out / "tables" / f"{k}.tsv" for k in TABLES}
        rec = _Recorder()
        log = logging.getLogger("gippa.steps.function")
        old = log.level
        log.addHandler(rec)
        log.setLevel(logging.DEBUG)
        try:
            reg.run_step("function", _step_cfg(), dict(ins), outputs)
        finally:
            log.removeHandler(rec)
            log.setLevel(old)
        return _StepRun(out, rec.records)


ALL_ITEMS = ("go_basic", "goslim_generic", "kegg_ko_pathway", "kegg_pathway_names", "cog2020_fun")
NO_KEGG_ITEMS = tuple(i for i in ALL_ITEMS if i not in KEGG_ITEMS)


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> _World:
    return _World(tmp_path_factory.mktemp("m12_03_world"))


@pytest.fixture(scope="module")
def with_kegg(world: _World) -> _StepRun:
    return world.run("with_kegg", ALL_ITEMS)


@pytest.fixture(scope="module")
def without_kegg(world: _World) -> _StepRun:
    return world.run("without_kegg", NO_KEGG_ITEMS)


def _non_kegg(df: pd.DataFrame) -> pd.DataFrame:
    return df[df["collection"] != "kegg"].reset_index(drop=True)


def test_with_kegg_the_world_has_kegg_go_and_cog_rows(with_kegg: _StepRun) -> None:
    """Item 5 sanity: the comparison below is not vacuous."""
    for key in TABLES:
        cols = set(with_kegg.frame(key)["collection"])
        assert "kegg" in cols, (key, cols)
        assert "cog" in cols, (key, cols)
    assert {"go_bp", "go_mf", "go_cc"} & set(with_kegg.frame("enrichment_detection")["collection"])
    assert "layer_missing:kegg" not in with_kegg.run_flags()


def test_without_kegg_the_tables_keep_their_header_and_drop_kegg_rows(
    with_kegg: _StepRun, without_kegg: _StepRun
) -> None:
    for key in TABLES:
        assert without_kegg.header(key) == with_kegg.header(key), key
        df = without_kegg.frame(key)
        assert "kegg" not in set(df["collection"]), (key, set(df["collection"]))
        assert len(df) > 0, f"{key}: GO and COG rows must remain"


def test_without_kegg_go_and_cog_still_run(without_kegg: _StepRun) -> None:
    flags = without_kegg.run_flags()
    assert "layer_missing:kegg" in flags, flags
    assert "layer_missing:go" not in flags and "layer_missing:cog" not in flags, flags
    det = set(without_kegg.frame("enrichment_detection")["collection"])
    assert "cog" in det and {"go_bp", "go_mf", "go_cc"} & det, det
    assert {"cog", "go_slim"} <= set(without_kegg.frame("category_summary")["collection"])


@pytest.mark.parametrize("key", ["enrichment_detection", "enrichment_gsea"])
def test_without_kegg_enrichment_equals_the_kegg_run_minus_kegg_rows(
    with_kegg: _StepRun, without_kegg: _StepRun, key: str
) -> None:
    """One rule: skipping KEGG removes the KEGG rows and changes nothing else (q is BH within a
    collection; the ORA draws nothing; every GSEA job seeds its own generator)."""
    pd.testing.assert_frame_equal(without_kegg.frame(key), _non_kegg(with_kegg.frame(key)))


def test_without_kegg_category_summary_equals_the_kegg_run_minus_kegg_rows_but_the_ci(
    with_kegg: _StepRun, without_kegg: _StepRun
) -> None:
    keep = [c for c in without_kegg.frame("category_summary").columns if c not in CI_COLUMNS]
    pd.testing.assert_frame_equal(without_kegg.frame("category_summary")[keep],
                                  _non_kegg(with_kegg.frame("category_summary"))[keep])  # fmt: skip


# ================================================================ 3. the run manifest


def test_manifest_records_kegg_skipped_and_why(isolated: Path) -> None:
    from helpers.manifest_data import event, sources, spec, validate

    from gippa.manifest import KEGG_NOTE, build_manifest
    from gippa.preflight.report import to_json_dict

    pre = to_json_dict(_real_preflight(isolated, "run"))
    db_doc = json.loads((isolated / "db" / "db_manifest.json").read_text(encoding="utf-8"))
    flags = event("function", "function.run_flags", {"flags": ["layer_missing:kegg"]})
    warn = event("function", "function.run_flags", {}, level="WARNING")
    m = build_manifest(sources(
        cfg=_cfg(), preflight=pre, db_manifest=db_doc, events=[flags, warn],
        registered=[spec("gene_table"), spec("function", module="function")],
        local_steps=[{"name": n, "status": "completed", "wall_time_s": 1.0, "cpu_time_s": 1.0}
                     for n in ("gene_table", "function")],
    ))  # fmt: skip
    validate(m)
    assert "kegg" in m["layers_missing"], m["layers_missing"]
    assert "layer_missing:kegg" in m["flags"], m["flags"]
    assert f"preflight:{RULE}" in m["warnings"], (
        f"the manifest must say why KEGG was skipped: {m['warnings']}"
    )
    assert KEGG_NOTE not in m["licence_notes"], m["licence_notes"]
    assert "function" not in m["modules_off"], m["modules_off"]


# ================================================================ 4. the report


def _report_out(tmp_path: Path, *, summary: bool = True) -> Path:
    from helpers.report_content import write_category_summary, write_gsea
    from helpers.report_fixtures import run_out

    out = run_out(tmp_path)
    if summary:
        write_category_summary(out, {"cog": 0.62, "go_slim": 0.55})
    write_gsea(out, [None])
    return out


def _function_text(out: Path) -> str:
    from helpers.report_content import page, section_text

    return section_text(page(out), "function")


def _no_kegg_manifest(out: Path) -> None:
    from helpers.report_content import manifest, write_manifest

    doc = manifest()
    doc["databases"] = [d for d in doc["databases"] if not str(d["name"]).startswith("kegg")]
    write_manifest(out, doc)


@pytest.mark.parametrize("summary", [True, False], ids=["with_summary", "gsea_only"])
def test_report_says_kegg_not_installed_pathway_enrichment_skipped(
    tmp_path: Path, summary: bool
) -> None:
    from helpers.report_content import write_flags

    out = _report_out(tmp_path, summary=summary)
    write_flags(out, "function", ["layer_missing:kegg"])
    _no_kegg_manifest(out)
    text = _function_text(out)
    assert REPORT_NOTE.lower() in text.lower(), text
    assert text.lower().count("kegg not installed") == 1, text
    assert "academic use only" not in text, text
    assert "KEGG retrieval date" not in text, text


def test_report_with_kegg_keeps_the_kegg_caveat(tmp_path: Path) -> None:
    """Item 5: KEGG installed, the academic-use caveat and its retrieval date stay."""
    from helpers.report_content import KEGG_DATE, write_manifest

    out = _report_out(tmp_path)
    write_manifest(out)
    text = _function_text(out)
    assert "academic use only" in text, text
    assert f"KEGG retrieval date: {KEGG_DATE}" in text, text
    assert "kegg not installed" not in text.lower(), text
    assert "pathway enrichment skipped" not in text.lower(), text


# ================================================================ 6. setup --no-kegg and docs


_UNAVAILABLE = re.compile(
    r"function[^\n]*(unavailable|not available|cannot run|can't run|disabled|requires? kegg)"
    r"|(--modules|modules) core,mobile,bgc(?!,function)",
    re.IGNORECASE,
)


@pytest.fixture
def fake_setup(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[str]:
    """``gippa setup`` with every side effect faked: installs are recorded, never run."""
    import shutil
    import socket
    import subprocess
    from collections import namedtuple

    import gippa.setupdb as sdb
    from gippa.setupdb import run as srun
    from gippa.setupdb.items import DbRecord

    installed: list[str] = []

    def fake_install(item: Any, db_dir: Path, _runner: Any, _fetcher: Any) -> DbRecord:
        installed.append(item.name)
        return DbRecord(name=item.name, version=item.version, url=item.url,
                        checksum=item.checksum or "sha256:" + "0" * 64, size_bytes=1,
                        retrieved="2026-11-04", path=item.path, status="ok",
                        message="downloaded")  # fmt: skip

    def deny(*a: Any, **_k: Any) -> Any:
        raise FileNotFoundError(f"blocked in this test: {a!r}")

    def deny_net(*a: Any, **_k: Any) -> Any:
        raise AssertionError(f"gippa setup touched the network: {a!r}")

    usage = namedtuple("usage", "total used free")
    monkeypatch.setattr(srun, "install", fake_install)
    monkeypatch.setattr(sdb.HttpsFetcher, "fetch", deny_net)
    monkeypatch.setattr(subprocess, "run", deny)
    monkeypatch.setattr(subprocess, "Popen", deny)
    monkeypatch.setattr(socket, "create_connection", deny_net)
    monkeypatch.setattr(socket.socket, "connect", deny_net)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: usage(1 << 50, 0, 1 << 50))
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR", "SLURM_CLUSTER_NAME"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg_config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg_data"))
    return installed


@pytest.mark.parametrize("extra", [[], ["--minimal"]], ids=["default", "minimal"])
def test_setup_no_kegg_never_says_function_is_unavailable(
    tmp_path: Path, fake_setup: list[str], extra: list[str]
) -> None:
    from click.testing import CliRunner

    from gippa.cli import cli

    argv = ["setup", "--db-dir", str(tmp_path / "db"), "--cache-dir", str(tmp_path / "cache"),
            "--no-images", "--no-nextflow", "--accept-download-plan", "--no-kegg",
            *extra]  # fmt: skip
    r = CliRunner().invoke(cli, argv, catch_exceptions=False)
    assert not set(KEGG_ITEMS) & set(fake_setup), fake_setup
    assert {"go_basic", "cog2020_fun"} <= set(fake_setup), fake_setup
    bad = [ln for ln in r.output.splitlines() if _UNAVAILABLE.search(ln)]
    assert bad == [], f"setup --no-kegg tells users function is unavailable: {bad}"


DOCS = ("docs/get-started/setup-databases.md", "docs/get-started/quickstart.md",
        "docs/get-started/manual-databases.md", "docs/get-started/install.md",
        "docs/user-guide/run.md", "docs/user-guide/licences.md", "README.md")  # fmt: skip
_DOC_REQUIRES_KEGG = (
    "Preflight requires all of GO, KEGG and COG",
    "both KEGG\ntables",
    "both KEGG tables",
)


@pytest.mark.parametrize("doc", DOCS)
def test_docs_do_not_tell_no_kegg_users_to_drop_function(doc: str) -> None:
    text = (REPO / doc).read_text(encoding="utf-8")
    no_kegg = re.compile(r"--no-kegg|without kegg", re.IGNORECASE)
    drop_function = re.compile(r"core,mobile,bgc(?!,function)")
    bad = [ln for ln in text.splitlines() if no_kegg.search(ln) and drop_function.search(ln)]
    assert bad == [], f"{doc}: --no-kegg users are told to run without function: {bad}"
    flat = " ".join(text.split())
    for phrase in _DOC_REQUIRES_KEGG:
        assert " ".join(phrase.split()) not in flat, (
            f"{doc}: says function requires KEGG ({phrase!r})"
        )
