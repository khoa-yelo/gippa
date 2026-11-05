"""M12-04 PR 10 vocabulary (docs/specs/M12-04-outputs-revision.md "PR 10" Docs, § 5.4; sign-off
C1, C9): ``putative_start`` joins the §9.4 Variants vocabulary (the flag order of
``gippa.variants.flags``, the writer's ``VARIANT_FLAGS``, the ``variants`` schema ``flags`` pattern
and the generated ``docs/reference/flags.md``), has plain words in ``FLAG_TEXT``, and is **not**
exempt from ``id_mixups``. ``docs/methods/variants.md`` gets a "Putative starts" note."""

from __future__ import annotations

import json
import re
from pathlib import Path

FLAG = "putative_start"
REPO = Path(__file__).resolve().parents[3]


def test_flag_order_holds_putative_start_after_database_derived() -> None:
    """One §9.4 order: the writer's (``vcf_records``) and the flags module's are the same."""
    import gippa.variants.vcf_records as vr
    from gippa.variants.flags import FLAG_ORDER

    assert FLAG_ORDER.count(FLAG) == 1, FLAG_ORDER
    assert FLAG_ORDER.index("database_derived") < FLAG_ORDER.index(FLAG)
    assert tuple(vr._FLAG_ORDER) == tuple(FLAG_ORDER)


def test_writer_vocabulary_accepts_the_flag() -> None:
    from gippa.variants.vcf_records import VARIANT_FLAGS

    assert FLAG in VARIANT_FLAGS


def test_variants_schema_flags_pattern_accepts_the_spec_example() -> None:
    schema = json.loads(
        (REPO / "src/gippa/schemas/tables/variants.schema.json").read_text(encoding="utf-8")
    )
    (field,) = [f for f in schema["fields"] if f["name"] == "flags"]
    pattern = field["constraints"]["pattern"]
    assert re.fullmatch(pattern, "snv;strain_difference;database_derived;putative_start")
    assert re.fullmatch(pattern, FLAG)


def test_reference_flags_page_lists_it_under_variants() -> None:
    text = (REPO / "docs/reference/flags.md").read_text(encoding="utf-8")
    section = text.split("## `variants`", 1)[1].split("\n## ", 1)[0]
    assert f"- `{FLAG}`" in section, section


def test_flag_text_has_plain_words() -> None:
    from gippa.viz.words import FLAG_TEXT, plain

    assert FLAG in FLAG_TEXT
    words = FLAG_TEXT[FLAG]
    assert words and "_" not in words and "§" not in words, words
    assert plain(FLAG) == words


def test_id_mixups_counts_a_putative_start_row() -> None:
    """C1 on the pure function: a same-db peptide on a ``putative_start`` row is a mix-up;
    ``EXEMPT_FLAGS`` stays ``{start_codon, stop_codon}`` while the flag is in the vocabulary."""
    import pandas as pd

    from gippa.variants.qc import EXEMPT_FLAGS, id_mixups
    from gippa.variants.vcf_records import VARIANT_FLAGS

    assert FLAG in VARIANT_FLAGS
    assert frozenset({"start_codon", "stop_codon"}) == EXEMPT_FLAGS

    sites = pd.DataFrame({"variant_id": ["VAR_1", "VAR_2"],
                          "flags": [("snv", FLAG), ("snv", "start_codon")]})  # fmt: skip
    support = pd.DataFrame({"variant_id": ["VAR_1", "VAR_2"], "peptide": ["MAAK", "MGGK"],
                            "source": ["database"] * 2, "search_db": ["same"] * 2})  # fmt: skip
    assert id_mixups(sites, support) == (1, 2)


def test_methods_note_on_putative_starts() -> None:
    text = (REPO / "docs/methods/variants.md").read_text(encoding="utf-8")
    assert re.search(r"^#+ Putative starts", text, flags=re.MULTILINE), "no Putative starts note"
    assert FLAG in text
