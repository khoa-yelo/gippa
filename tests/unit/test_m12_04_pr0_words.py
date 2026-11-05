"""M12-04 PR 0: ``gippa.viz.words``, the one plain-text dictionary for flags, skip reasons,
feature names and tiers (docs/specs/M12-04-outputs-revision.md § PR 0, AC 0.2, 0.3; § 1 rule 2).

Interface assumed beyond the spec text:
- A parameterised token (``layer_missing:<layer>``, ``missing_in:<run>``, ``module_off:<m>``,
  ``input_absent:<key>``) "is looked up by its prefix and then filled in". A dictionary therefore
  covers it when it holds the exact token, the prefix with the colon (``"layer_missing:"``) or the
  bare prefix (``"layer_missing"``). The test accepts any of the three keys.
- ``plain()`` covers both ``FLAG_TEXT`` and ``SKIP_TEXT`` (PR 2: "flags and skip reasons in prose
  go through ``plain()``"); its "logged warning" goes through ``logging`` at level WARNING.
- ``TIER_TEXT`` values must pass AC 0.3 (no underscore), so they cannot start with the raw tier
  token the PR 1 table shows ("decoy_picked: ..."); the caller adds the tier name. The tests only
  check the plain-word content of each sentence, and any ``{alpha}``-style placeholder is allowed.
"""

from __future__ import annotations

import importlib
import json
import logging
import pkgutil
import re
from collections.abc import Mapping
from pathlib import Path

import pytest
from hypothesis import given
from hypothesis import strategies as st

ROOT = Path(__file__).resolve().parents[2]
SCHEMAS = ROOT / "src" / "gippa" / "schemas" / "tables"

# GIPPA_SPEC.md §9.4, scope Run (layer_missing:<layer> is expanded below).
SPEC_RUN_FLAGS = ("anti_sd_default", "c_assumed", "bakta_db_light", "operon_source_distance_model")
# Every layer name src/gippa emits as layer_missing:<layer> on main (2026-11-04).
LAYERS = (
    "abundance", "antismash", "codon", "cog", "features", "genomad", "go", "hits_tolerant",
    "kegg", "ostir", "promoters", "rbs", "regions", "terminators", "tx_rates",
)  # fmt: skip
TIERS = ("decoy_picked", "mayu", "two_peptide", "protein_list")
FORBIDDEN = re.compile(r"ADR-\d|§|GIPPA_SPEC|_")


def _words():  # type: ignore[no-untyped-def]
    return importlib.import_module("gippa.viz.words")


def _alternatives(pattern: str) -> list[str]:
    """The alternatives of a flat ``(?:a|b)(?:;(?:a|b))*`` flags pattern."""
    m = re.match(r"\(\?:([^()]*)\)", pattern)
    assert m, f"unexpected flags pattern {pattern!r}"
    return m.group(1).split("|")


def _schema_flag_tokens() -> list[tuple[str, str]]:
    """``(table, token)`` for every token of every flags pattern / enum in the table schemas.
    ``operon_flags.flag`` (renamed ``flags`` by PR 11b, Q2) holds §9.4 Operons tokens too."""
    out: set[tuple[str, str]] = set()
    for path in sorted(SCHEMAS.glob("*.schema.json")):
        schema = json.loads(path.read_text(encoding="utf-8"))
        for field in schema.get("fields", []):
            if field.get("name") not in ("flags", "flag"):
                continue
            cons = field.get("constraints") or {}
            tokens = [str(v) for v in cons.get("enum", [])]
            if not tokens and "pattern" in cons:
                tokens = _alternatives(str(cons["pattern"]))
            out.update((str(schema["name"]), t) for t in tokens)
    return sorted(out)


SCHEMA_TOKENS = _schema_flag_tokens()


def _concrete(token: str) -> str:
    """A pattern alternative with a regex parameter (``missing_in:[^;]+``) as a real token."""
    return re.sub(r"\[\^;\]\+$", "run2", token)


def _has_entry(table: Mapping[str, str], token: str) -> bool:
    if token in table:
        return True
    if ":" in token:
        prefix = token.split(":", 1)[0]
        return f"{prefix}:" in table or prefix in table
    return False


def _warnings(caplog: pytest.LogCaptureFixture) -> list[logging.LogRecord]:
    return [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_schema_scan_found_the_vocabulary() -> None:
    """Guard: the parametrisation below is not silently empty."""
    tokens = {t for _, t in SCHEMA_TOKENS}
    assert {"partial", "mass_equivalent", "rbs_supported", "unexplained_discordance"} <= tokens
    assert len(SCHEMA_TOKENS) > 80


# --- AC 0.2: every flag token has a sentence -----------------------------------------------------


@pytest.mark.parametrize(
    ("table", "token"), SCHEMA_TOKENS, ids=[f"{a}-{b}" for a, b in SCHEMA_TOKENS]
)
def test_every_schema_flag_has_flag_text(table: str, token: str) -> None:
    words = _words()
    concrete = _concrete(token)
    assert _has_entry(words.FLAG_TEXT, concrete), f"{table}: {token!r} has no FLAG_TEXT entry"


@pytest.mark.parametrize(
    ("table", "token"), SCHEMA_TOKENS, ids=[f"{a}-{b}" for a, b in SCHEMA_TOKENS]
)
def test_plain_turns_every_schema_flag_into_words(
    table: str, token: str, caplog: pytest.LogCaptureFixture
) -> None:
    words = _words()
    concrete = _concrete(token)
    caplog.set_level(logging.WARNING)
    text = words.plain(concrete)
    assert isinstance(text, str) and text.strip()
    assert text != concrete, f"{table}: plain({concrete!r}) returned the raw token"
    assert not _warnings(caplog), f"{table}: plain({concrete!r}) warned"


RUN_TOKENS = [*SPEC_RUN_FLAGS, *(f"layer_missing:{layer}" for layer in LAYERS)]


@pytest.mark.parametrize("token", RUN_TOKENS)
def test_every_run_flag_has_flag_text(token: str, caplog: pytest.LogCaptureFixture) -> None:
    words = _words()
    assert _has_entry(words.FLAG_TEXT, token), f"run flag {token!r} has no FLAG_TEXT entry"
    caplog.set_level(logging.WARNING)
    text = words.plain(token)
    assert text != token and text.strip()
    assert not _warnings(caplog)


def test_flag_text_spec_examples() -> None:
    words = _words()
    assert words.plain("layer_missing:codon") == "codon usage not computed"
    assert words.plain("c_assumed") == "de novo precision assumed, not measured"


def test_parameterised_flag_is_filled_in(caplog: pytest.LogCaptureFixture) -> None:
    """Edge case: ``layer_missing:<layer>`` is looked up by its prefix and then filled in, so an
    unlisted layer still gets a sentence that names it."""
    words = _words()
    caplog.set_level(logging.WARNING)
    text = words.plain("layer_missing:zzlayer")
    assert "zzlayer" in text
    assert text != "layer_missing:zzlayer"
    assert "layer_missing" not in text
    text = words.plain("missing_in:runB")
    assert "runB" in text and "missing_in" not in text


# --- SKIP_TEXT -----------------------------------------------------------------------------------


def test_skip_text_spec_examples() -> None:
    words = _words()
    assert words.plain("module_off:function") == "the function module was not run"
    assert words.plain("input_absent:tables/variants.tsv") == "no variant table in this run"


@pytest.mark.parametrize("prefix", ["module_off", "layer_missing", "input_absent"])
def test_skip_text_covers_every_declared_skip_prefix(prefix: str) -> None:
    words = _words()
    assert _has_entry(words.SKIP_TEXT, f"{prefix}:x") or _has_entry(
        words.FLAG_TEXT, f"{prefix}:x"
    ), prefix


def _registered_skip_reasons() -> list[str]:
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base

    load_plot_modules()
    reasons: set[str] = set()
    for spec in base.registered().values():
        if spec.module != "core":
            reasons.add(f"module_off:{spec.module}")
        reasons.update(f"layer_missing:{layer}" for layer in spec.layers)
        reasons.update(f"input_absent:{key}" for key in spec.needs)
    return sorted(reasons)


def test_every_skip_reason_a_registered_plot_can_return_is_plain(
    caplog: pytest.LogCaptureFixture,
) -> None:
    words = _words()
    reasons = _registered_skip_reasons()
    assert any(r.startswith("input_absent:") for r in reasons)
    caplog.set_level(logging.WARNING)
    raw = [r for r in reasons if words.plain(r) == r]
    assert not raw, f"skip reasons without a sentence: {raw}"
    assert not _warnings(caplog)


# --- FEATURE_LABEL and TIER_TEXT -----------------------------------------------------------------


@pytest.mark.parametrize(
    ("column", "label"),
    [
        ("tir_ostir", "translation initiation rate (OSTIR)"),
        ("cpb", "codon pair bias"),
        ("at_2_10", "AT content, codons 2\u201310"),
        ("terminator_conf_downstream", "downstream terminator confidence"),
    ],
)
def test_feature_label_spec_examples(column: str, label: str) -> None:
    assert _words().FEATURE_LABEL[column] == label


def test_feature_labels_are_human_names() -> None:
    """A "human name" carries no raw column token (§ 1 rule 2)."""
    bad = {k: v for k, v in _words().FEATURE_LABEL.items() if FORBIDDEN.search(v)}
    assert not bad


def test_tier_text_covers_every_tier() -> None:
    assert set(TIERS) <= set(_words().TIER_TEXT)


@pytest.mark.parametrize(
    ("tier", "phrases"),
    [
        ("decoy_picked", ("gene-level q-value", "target\u2013decoy competition")),
        ("mayu", ("accepted peptide", "Mayu")),
        ("two_peptide", ("≥ 2 distinct accepted peptides", "no FDR is estimated")),
        ("protein_list", ("protein in your list", "no FDR is estimated")),
    ],
)
def test_tier_text_sentences(tier: str, phrases: tuple[str, ...]) -> None:
    text = _words().TIER_TEXT[tier]
    assert "a gene counts as detected when" in text
    for phrase in phrases:
        assert phrase in text, (tier, phrase)


# --- AC 0.3: no internal references or raw tokens in the values ---------------------------------


@pytest.mark.parametrize("name", ["FLAG_TEXT", "SKIP_TEXT", "TIER_TEXT"])
def test_values_hold_no_internal_reference_or_raw_token(name: str) -> None:
    table = getattr(_words(), name)
    assert table, f"{name} is empty"
    bad = {k: v for k, v in table.items() if FORBIDDEN.search(v)}
    assert not bad, f"{name} values match {FORBIDDEN.pattern}: {bad}"


@pytest.mark.parametrize("name", ["FLAG_TEXT", "SKIP_TEXT", "TIER_TEXT", "FEATURE_LABEL"])
def test_tables_are_str_to_nonempty_str(name: str) -> None:
    table = getattr(_words(), name)
    assert all(isinstance(k, str) and isinstance(v, str) and v.strip() for k, v in table.items())


# --- plain(): unknown tokens, never raises -------------------------------------------------------


def test_plain_unknown_token_returns_it_with_a_warning(caplog: pytest.LogCaptureFixture) -> None:
    caplog.set_level(logging.WARNING)
    assert _words().plain("no_such_flag_zz") == "no_such_flag_zz"
    assert any("no_such_flag_zz" in r.getMessage() for r in _warnings(caplog))


@given(st.text(max_size=40))
def test_plain_property_never_raises_and_returns_text(token: str) -> None:
    out = _words().plain(token)
    assert isinstance(out, str)


def test_words_module_is_pure() -> None:
    """``gippa.viz.words`` is a pure module (no file I/O)."""
    import gippa.viz

    names = {m.name for m in pkgutil.iter_modules(gippa.viz.__path__)}
    assert "words" in names
    src = (ROOT / "src" / "gippa" / "viz" / "words.py").read_text(encoding="utf-8")
    assert "open(" not in src and "read_text" not in src
