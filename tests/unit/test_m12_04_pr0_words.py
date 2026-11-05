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


# Parameters of the run-scope templates (``uncalibrated:<key>`` keys of the calibration block).
TEMPLATE_PARAMS = {
    "layer_missing:": LAYERS,
    "uncalibrated:": ("operon_distance", "uniop_pair_threshold", "terminatornet_yes"),
}


def _run_tokens() -> list[str]:
    """Review M1: every token of ``gippa.manifest.RUN_FLAGS``, templates expanded."""
    from gippa.manifest import RUN_FLAGS

    out: set[str] = set(SPEC_RUN_FLAGS)
    for token in RUN_FLAGS:
        if token.endswith(":"):
            params = TEMPLATE_PARAMS.get(token, ("some_param",))
            out.update(f"{token}{p}" for p in params)
        else:
            out.add(token)
    return sorted(out)


RUN_TOKENS = _run_tokens()


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


# --- review round 1 (PR #549) ------------------------------------------------------------------
# Lead decisions H1, M1-M3, M7, TIER_NAME, feature_label and the template mutation gaps.
# Interface assumed: ``plain(token, table=None)``; ``TIER_NAME`` and ``feature_label(col)`` live
# in gippa.viz.words.


def test_cai_reference_says_the_gene_is_in_the_reference_set() -> None:
    """H1: ``cai_reference`` marks a ribosomal protein in the codon-adaptation reference set,
    left out of the model fit; it is not "outside" the set."""
    text = _words().plain("cai_reference")
    assert "reference set" in text
    assert "outside" not in text
    assert "fit" in text


# --- M1: tokens the determinants, null and de novo precision modules emit ----------------------

_EMITTING_MODULES = (
    "gippa.determinants.models",
    "gippa.determinants.null",
    "gippa.determinants.attribution",
    "gippa.steps.determinants",
    "gippa.mapping.denovo_precision",
)
_TOKEN = re.compile(r"[a-z][a-z0-9_]*(:[a-z0-9_]*)?")


def _emitted_flags(path: Path) -> set[str]:
    """Flag tokens a module writes, found by an AST scan of its flag contexts: assignments to a
    ``*flag*`` name, ``"flags"`` dict entries, arguments of ``*flag*`` calls, arguments bound to a
    ``*flag*`` parameter or dataclass field of the same module, ``x.append(...)`` on a list the
    module ``";".join``s, and the returns of ``*flag*`` functions. Module-level string constants
    are resolved by name."""
    import ast

    tree = ast.parse(path.read_text(encoding="utf-8"))
    consts: dict[str, str] = {}
    for node in tree.body:
        if isinstance(node, (ast.Assign, ast.AnnAssign)) and isinstance(node.value, ast.Constant):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for t in targets:
                if isinstance(t, ast.Name) and isinstance(node.value.value, str):
                    consts[t.id] = node.value.value

    def ext(e: ast.AST | None) -> set[str]:
        if e is None:
            return set()
        if isinstance(e, ast.Constant) and isinstance(e.value, str):
            return {p for p in e.value.split(";") if _TOKEN.fullmatch(p)}
        if isinstance(e, ast.Name):
            v = consts.get(e.id)
            return ext(ast.Constant(v)) if v is not None else set()
        if isinstance(e, ast.IfExp):
            return ext(e.body) | ext(e.orelse)
        if isinstance(e, (ast.List, ast.Tuple, ast.Set)):
            return set().union(*(ext(x) for x in e.elts)) if e.elts else set()
        if isinstance(e, ast.BoolOp):
            return set().union(*(ext(x) for x in e.values))
        if (isinstance(e, ast.Call) and isinstance(e.func, ast.Attribute)
                and e.func.attr == "join" and e.args):  # fmt: skip
            return ext(e.args[0])
        return set()

    def name_of(f: ast.AST) -> str:
        return f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else ""

    # parameters (functions) and fields (classes) by position and name
    params: dict[str, list[str]] = {}
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            params[node.name] = [a.arg for a in [*node.args.posonlyargs, *node.args.args]]
        elif isinstance(node, ast.ClassDef):
            fields = [s.target for s in node.body if isinstance(s, ast.AnnAssign)]
            params[node.name] = [t.id for t in fields if isinstance(t, ast.Name)]
    joined = {name_of(n.args[0]) for n in ast.walk(tree)
              if isinstance(n, ast.Call) and name_of(n.func) == "join" and n.args}  # fmt: skip

    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Assign, ast.AnnAssign)):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            if any(isinstance(t, ast.Name) and "flag" in t.id.lower() for t in targets):
                out |= ext(node.value)
        elif isinstance(node, ast.Dict):
            for k, v in zip(node.keys, node.values, strict=True):
                if isinstance(k, ast.Constant) and k.value == "flags":
                    out |= ext(v)
        elif isinstance(node, ast.Call):
            fname = name_of(node.func)
            if "flag" in fname.lower() and fname != "_flag":  # ``_flag``: a bool-column helper
                for a in node.args:
                    out |= ext(a)
            receiver = name_of(node.func.value) if isinstance(node.func, ast.Attribute) else ""
            if fname == "append" and (receiver in joined or "flag" in receiver) and node.args:
                out |= ext(node.args[0])
            sig = params.get(fname, [])
            if fname == "append":
                sig = []
            for i, a in enumerate(node.args):
                if i < len(sig) and "flag" in sig[i].lower():
                    out |= ext(a)
            for kw in node.keywords:
                if kw.arg and "flag" in kw.arg.lower():
                    out |= ext(kw.value)
        elif isinstance(node, ast.FunctionDef) and "flags" in node.name.lower():
            for sub in ast.walk(node):
                if isinstance(sub, ast.Return):
                    out |= ext(sub.value)
    return out


def _module_flag_tokens() -> list[str]:
    out: set[str] = set()
    for mod in _EMITTING_MODULES:
        path = ROOT / "src" / Path(*mod.split(".")).with_suffix(".py")
        out |= _emitted_flags(path)
    # the documented lists of the two determinants tables (their flags field has no pattern)
    for name in ("determinants_model", "determinants_attribution"):
        schema = json.loads((SCHEMAS / f"{name}.schema.json").read_text(encoding="utf-8"))
        for field in schema["fields"]:
            if field["name"] == "flags":
                listed = str(field.get("description", "")).split(":", 1)[-1]
                out.update(t.strip() for t in listed.split(",") if _TOKEN.fullmatch(t.strip()))
    return sorted(out)


MODULE_TOKENS = _module_flag_tokens()


def test_module_flag_scan_found_the_tokens() -> None:
    """Guard: the scan sees the tokens known to be emitted (2026-11-04), so it is not silently
    empty; a new token joins the parametrisation below on its own."""
    known = {
        "too_few_rows", "unestimated", "null_degenerate", "null_uncalibrated", "ci_unestimated",
        "abundance_not_explained_by_sequence", "cai_reference", "c_assumed", "unreported",
        "layer_missing:abundance",
    }  # fmt: skip
    assert known <= set(MODULE_TOKENS), sorted(known - set(MODULE_TOKENS))


@pytest.mark.parametrize("token", MODULE_TOKENS)
def test_every_module_emitted_flag_has_flag_text(
    token: str, caplog: pytest.LogCaptureFixture
) -> None:
    words = _words()
    assert _has_entry(words.FLAG_TEXT, token), f"{token!r} has no FLAG_TEXT entry"
    caplog.set_level(logging.WARNING)
    text = words.plain(token)
    assert text != token and text.strip()
    assert not _warnings(caplog)


@pytest.mark.parametrize("token", ["too_few_rows", "null_uncalibrated", "unreported"])
def test_review_m1_tokens_have_exact_entries(token: str) -> None:
    text = _words().FLAG_TEXT[token]
    assert text.strip() and not FORBIDDEN.search(text)


# --- M2: corrected meanings and the table scope ------------------------------------------------


def test_tolerant_evidence_means_has_tolerant_evidence() -> None:
    text = _words().plain("tolerant_evidence")
    assert "mismatch-tolerant" in text
    assert "only" not in text


def test_ci_few_units_does_not_claim_a_wrong_coverage() -> None:
    text = _words().plain("ci_few_units")
    assert "not a 95" not in text
    assert "unreliable" in text or "not computed" in text


def test_initiator_m_mentions_a_non_atg_start() -> None:
    text = _words().plain("initiator_m")
    assert "ATG" in text and "start codon" in text
    assert re.search(r"non-ATG|not ATG|other than ATG", text), text


def test_abundance_imputed_pair_scope() -> None:
    words = _words()
    generic = words.plain("abundance_imputed")
    scoped = words.plain("abundance_imputed", table="pairs")
    assert "pair" not in generic
    assert "a gene of the pair" in scoped
    assert scoped != generic


def test_null_degenerate_is_table_dependent() -> None:
    words = _words()
    generic = words.plain("null_degenerate")
    scoped = words.plain("null_degenerate", table="determinants_model")
    assert scoped != generic
    assert "detection" not in scoped  # the determinants null permutes TUs, not detections
    assert not FORBIDDEN.search(scoped)


@pytest.mark.parametrize("token", ["partial", "c_assumed", "layer_missing:codon", "no_such_zz"])
@pytest.mark.parametrize("table", [None, "pairs", "determinants_model", "genes", "zz_table"])
def test_plain_scope_falls_back_to_the_generic_sentence(token: str, table: str | None) -> None:
    words = _words()
    assert words.plain(token, table=table) == words.plain(token)
    assert words.plain(token, None) == words.plain(token)


@given(st.text(max_size=30), st.one_of(st.none(), st.text(max_size=20)))
def test_plain_property_with_scope_never_raises(token: str, table: str | None) -> None:
    assert isinstance(_words().plain(token, table=table), str)


# --- M3: templates never leak an underscore or a placeholder -----------------------------------


@pytest.mark.parametrize(
    "token",
    ["uncalibrated:operon_distance", "uncalibrated:uniop_pair_threshold",
     "uncalibrated:terminatornet_yes"],
)  # fmt: skip
def test_uncalibrated_keys_have_exact_entries(token: str) -> None:
    words = _words()
    assert token in words.FLAG_TEXT
    assert words.plain(token) == words.FLAG_TEXT[token]
    assert "_" not in words.plain(token)


def test_uncalibrated_operon_distance_names_it_in_words() -> None:
    assert "operon distance" in _words().plain("uncalibrated:operon_distance")


@pytest.mark.parametrize(
    ("token", "phrase"),
    [
        ("layer_missing:new_layer", "new layer"),
        ("uncalibrated:some_new_key", "some new key"),
        ("module_off:new_module", "new module"),
    ],
)
def test_template_fill_puts_no_underscore_into_prose(token: str, phrase: str) -> None:
    text = _words().plain(token)
    assert phrase in text
    assert "_" not in text


@pytest.mark.parametrize("token", ["layer_missing:", "uncalibrated:", "missing_in:",
                                   "module_off:", "input_absent:"])  # fmt: skip
def test_bare_template_key_never_returns_a_placeholder(token: str) -> None:
    text = _words().plain(token)
    assert "{" not in text and "}" not in text
    assert "_" not in text


# ``missing_in:<run>`` is left out: its parameter is a user's run name, shown verbatim.
@given(st.sampled_from(["layer_missing:", "uncalibrated:", "module_off:"]),
       st.text(alphabet="abcxyz_", max_size=12))  # fmt: skip
def test_template_property_no_underscore_no_placeholder(prefix: str, param: str) -> None:
    text = _words().plain(prefix + param)
    assert "_" not in text and "{}" not in text


# --- M7: no value embeds an underscored flag token ---------------------------------------------


def _underscored_tokens() -> set[str]:
    words = _words()
    keys = {k.rstrip(":") for k in words.FLAG_TEXT} | {k.rstrip(":") for k in words.SKIP_TEXT}
    keys |= {t for _, t in SCHEMA_TOKENS} | set(MODULE_TOKENS) | set(RUN_TOKENS)
    return {k for k in keys if "_" in k}


@pytest.mark.parametrize("name", ["FLAG_TEXT", "SKIP_TEXT", "TIER_TEXT", "FEATURE_LABEL",
                                  "TIER_NAME"])  # fmt: skip
def test_no_value_contains_an_underscored_flag_token(name: str) -> None:
    table = getattr(_words(), name)
    tokens = _underscored_tokens()
    assert "abundance_imputed" in tokens
    bad = {k: v for k, v in table.items() for t in tokens if t in v}
    assert not bad, f"{name}: {bad}"


# --- TIER_NAME ---------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("tier", "name"),
    [("decoy_picked", "decoy-picked"), ("two_peptide", "two-peptide"), ("mayu", "Mayu"),
     ("protein_list", "protein list")],
)  # fmt: skip
def test_tier_name(tier: str, name: str) -> None:
    assert _words().TIER_NAME[tier] == name


def test_tier_name_covers_tier_text_without_underscores() -> None:
    words = _words()
    assert set(words.TIER_NAME) == set(words.TIER_TEXT) >= set(TIERS)
    assert not [v for v in words.TIER_NAME.values() if FORBIDDEN.search(v) or not v.strip()]


# --- feature_label -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("column", "label"),
    [
        ("start_codon_ATG", "start codon: ATG"),
        ("start_codon_other", "start codon: other"),
        ("replicon_type_plasmid", "replicon type: plasmid"),
        ("nterm_residue_M", "N-terminal residue: M"),
        ("tir_ostir", "translation initiation rate (OSTIR)"),
        ("cpb", "codon pair bias"),
    ],
)
def test_feature_label_examples(column: str, label: str) -> None:
    assert _words().feature_label(column) == label


def _all_feature_columns() -> list[str]:
    from gippa.determinants import features

    cols = set(features.FEATURE_COLUMNS) | set(features.IDENTITY_COLUMNS)
    for feat, levels in features.LEVELS.items():
        cols.update(f"{feat}_{lv}" for lv in levels)
    return sorted(cols)


@pytest.mark.parametrize("column", _all_feature_columns())
def test_feature_label_every_column_is_plain(column: str) -> None:
    label = _words().feature_label(column)
    assert label.strip() and label != column
    assert not FORBIDDEN.search(label), label


def test_feature_label_uses_feature_label_for_plain_features() -> None:
    words = _words()
    for col, label in words.FEATURE_LABEL.items():
        assert words.feature_label(col) == label


@given(st.text(max_size=30))
def test_feature_label_property_never_an_underscore(column: str) -> None:
    out = _words().feature_label(column)
    assert isinstance(out, str) and "_" not in out
