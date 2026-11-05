"""M12-04 PR 12 (O-01): step intermediates under ``OUT/intermediate/``; published paths unchanged
(docs/specs/M12-04-outputs-revision.md > PR 12, AC 12.1, 12.2, 12.4; ADR-0075 decision 3).

Interface pinned here (test author, from the PR 12 text):
- the published top level of ``OUT/`` is ``README.txt``, ``run_manifest.json``,
  ``resolved_config.yaml``, ``preflight/``, ``annotation/``, ``mapping/``, ``tables/``, ``plots/``,
  ``report/``, ``schemas/``, ``logs/``, ``intermediate/`` and ``work/``;
- ``gippa.runner.local.PUBLISHED_FILES`` keeps the Nextflow ``key=path`` names (pinned against
  ``workflow/`` by tests/unit/test_local_publish_names.py); ``gippa.runner.local.publish`` copies
  every entry whose first path component is not a published name to ``intermediate/<path>``;
  ``work/…`` entries stay unpublished;
- ``gippa.io.run_view.PATHS`` gives those keys as ``intermediate/<1.0 path>``; ``_FALLBACK`` holds
  the 1.0 path, so ``RunView.from_outdir`` still reads a 1.0 tree. When both exist (a ``--resume``
  into a 1.0 directory leaves the old copies in place) the ``intermediate/`` copy wins.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

from gippa.io.run_view import _FALLBACK, PATHS, RunView
from gippa.runner.local import PUBLISHED, PUBLISHED_FILES, SUPPORT_OUTPUTS, publish

#: AC 12.1: the only top-level entries of OUT/ after a run
TOP_LEVEL = frozenset({
    "README.txt", "run_manifest.json", "resolved_config.yaml", "preflight", "annotation",
    "mapping", "tables", "plots", "report", "schemas", "logs", "intermediate", "work",
})  # fmt: skip
PUBLISHED_HEADS = frozenset({"preflight", "annotation", "mapping", "tables", "plots", "report",
                             "schemas", "logs"})  # fmt: skip

#: the 1.0 OUT-relative paths of the step intermediates PR 12 moves (main@v1.0.0)
MOVED_1_0: dict[str, str] = {
    "detectability_summary": "genetable/detection_summary_parts/detectability.parquet",
    "gene_context": "context/gene_context.parquet",
    "gene_table": "genetable/genes.parquet",
    "genes_pre": "genetable/genes_pre.parquet",
    "island_diagnostics": "context/island_diagnostics.parquet",
    "novel_class_fdr": "novel/novel_class_fdr.parquet",
    "tier_summary": "genetable/detection_summary_parts/tier.parquet",
    "tu_folds": "genetable/tu_folds.parquet",
}

#: AC 12.2: published result paths of 1.0 (main@v1.0.0), unchanged in 1.1
PUBLISHED_1_0: dict[str, str] = {
    "autocorrelation": "tables/autocorrelation.tsv",
    "category_summary": "tables/category_summary.tsv",
    "detection_summary": "tables/detection_summary.tsv",
    "determinants_attribution": "tables/determinants_attribution.tsv",
    "determinants_features": "tables/determinants_features.tsv",
    "determinants_model": "tables/determinants_model.tsv",
    "enrichment_detection": "tables/enrichment_detection.tsv",
    "enrichment_gsea": "tables/enrichment_gsea.tsv",
    "genes_tsv": "tables/genes.tsv",
    "islands": "tables/islands.tsv",
    "novel_orfs": "tables/novel_orfs.tsv",
    "novel_orfs_faa": "tables/novel_orfs.faa",
    "novel_orfs_gff3": "tables/novel_orfs.gff3",
    "annotation_patch": "tables/annotation_patch.gff3",
    "operon_flags": "tables/operon_flags.tsv",
    "operons": "tables/operons.tsv",
    "pairs": "tables/pairs.tsv",
    "regions_stats": "tables/regions_stats.tsv",
    "pileup": "tables/pileup.tsv.gz",
    "sensitivity": "tables/sensitivity.tsv",
    "variants": "tables/variants.tsv",
    "variants_vcf": "tables/variants.vcf",
    "preflight": "preflight/preflight.json",
    "peptide_hits": "mapping/peptide_hits.tsv",
    "protein_hits": "mapping/protein_hits.tsv",
    "unmapped": "mapping/unmapped.tsv",
    "mapping_fdr": "mapping/mapping_fdr.tsv",
    "tus_tsv": "annotation/tus.tsv",
    "regions_tsv": "annotation/regions.tsv",
    "promoters_bed": "annotation/promoters.bed",
    "terminators_bed": "annotation/terminators.bed",
    "events": "logs/events.jsonl",
    "resolved_config": "resolved_config.yaml",
    "run_manifest": "run_manifest.json",
}


def _frame(tag: str) -> pd.DataFrame:
    return pd.DataFrame({"locus_tag": [f"{tag}_1", f"{tag}_2"], "value": [1.0, 2.0]})


def _write_parquet(path: Path, tag: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    _frame(tag).to_parquet(path, index=False)


def _top(out: Path) -> set[str]:
    return {p.name for p in out.iterdir() if not p.name.startswith(".")}


def _work(root: Path, keys: list[str]) -> Path:
    """A local-executor work dir holding ``workdir / key`` for ``keys`` (files), every
    ``PUBLISHED`` directory and every analysis support directory."""
    work = root / "work"
    work.mkdir(parents=True)
    for key in keys:
        (work / key).write_text(f"{key}\n")
    for key in PUBLISHED:
        (work / key).mkdir()
        (work / key / "x.tsv").write_text("x\n")
    for key in SUPPORT_OUTPUTS:
        (work / key / "support").mkdir(parents=True)
        (work / key / "support" / f"{key}.tsv").write_text("s\n")
    return work


# ---------------------------------------------------------------- PATHS (reads)


@pytest.mark.parametrize("key", sorted(MOVED_1_0))
def test_paths_of_moved_keys_are_under_intermediate(key: str) -> None:
    assert PATHS[key] == f"intermediate/{MOVED_1_0[key]}", (key, PATHS[key])


@pytest.mark.parametrize("key", sorted(MOVED_1_0))
def test_fallback_holds_the_1_0_location(key: str) -> None:
    assert MOVED_1_0[key] in _FALLBACK.get(key, ()), (key, _FALLBACK.get(key))


@pytest.mark.parametrize("key", sorted(PUBLISHED_1_0))
def test_published_paths_are_unchanged_from_1_0(key: str) -> None:
    """AC 12.2 at the key level: every published result keeps its 1.0 path."""
    assert PATHS[key] == PUBLISHED_1_0[key]


def test_every_path_is_published_intermediate_or_work() -> None:
    """AC 12.1 at the key level: no key reads a top-level location outside the published list."""
    bad = {k: p for k, p in PATHS.items() if p.split("/")[0] not in TOP_LEVEL}
    assert not bad, f"keys outside the published top level: {bad}"


def test_resolved_config_keeps_its_work_fallback() -> None:
    assert "work/resolved_config.yaml" in _FALLBACK["resolved_config"]


# ---------------------------------------------------------------- RunView.from_outdir


@pytest.mark.parametrize("key", sorted(MOVED_1_0))
def test_from_outdir_reads_the_1_1_layout(tmp_path: Path, key: str) -> None:
    out = tmp_path / "OUT"
    (out / "tables").mkdir(parents=True)
    _write_parquet(out / "intermediate" / MOVED_1_0[key], key)
    view = RunView.from_outdir(out)
    assert view.path(key) == out / "intermediate" / MOVED_1_0[key]
    pd.testing.assert_frame_equal(view.table(key), _frame(key))


@pytest.mark.parametrize("key", sorted(MOVED_1_0))
def test_from_outdir_reads_a_1_0_tree_through_the_fallbacks(tmp_path: Path, key: str) -> None:
    """AC 12.4: intermediates at the top level (1.0 layout) are still found."""
    out = tmp_path / "OUT"
    (out / "tables").mkdir(parents=True)
    _write_parquet(out / MOVED_1_0[key], key)
    view = RunView.from_outdir(out)
    assert view.path(key) == out / MOVED_1_0[key]
    pd.testing.assert_frame_equal(view.table(key), _frame(key))


@pytest.mark.parametrize("key", ["gene_table", "gene_context", "island_diagnostics"])
def test_intermediate_copy_wins_over_a_legacy_copy(tmp_path: Path, key: str) -> None:
    """``--resume`` into a 1.0 directory leaves the old files in place; the new run's
    ``intermediate/`` copy must be the one read."""
    out = tmp_path / "OUT"
    (out / "tables").mkdir(parents=True)
    _write_parquet(out / MOVED_1_0[key], "legacy")
    _write_parquet(out / "intermediate" / MOVED_1_0[key], "new")
    view = RunView.from_outdir(out)
    assert view.path(key) == out / "intermediate" / MOVED_1_0[key]
    pd.testing.assert_frame_equal(view.table(key), _frame("new"))


def test_from_inputs_still_finds_a_moved_key_by_file_name(tmp_path: Path) -> None:
    """A step's input map names the staged file (``genes.parquet``) or its directory."""
    staged = tmp_path / "in1" / "genes.parquet"
    _write_parquet(staged, "gt")
    assert RunView.from_inputs({"gene_table": staged}).path("gene_table") == staged
    assert RunView.from_inputs({"gene_table": staged.parent}).path("gene_table") == staged


# ---------------------------------------------------------------- local publish (AC 12.1)


def test_local_publish_moves_intermediates_and_keeps_results(tmp_path: Path) -> None:
    keys = sorted(PUBLISHED_FILES)
    work, out = _work(tmp_path, keys), tmp_path / "OUT"
    done = publish(work, out)
    assert _top(out) <= TOP_LEVEL, sorted(_top(out) - TOP_LEVEL)
    for key, rel in MOVED_1_0.items():
        assert (out / "intermediate" / rel).read_text() == f"{key}\n", rel
        assert not (out / rel).exists(), f"{rel} still published at the top level"
        assert out / "intermediate" / rel in done
    for key, rel in PUBLISHED_FILES.items():
        if rel.split("/")[0] in PUBLISHED_HEADS:
            assert (out / rel).read_text() == f"{key}\n", rel
    for sub in PUBLISHED.values():
        assert (out / sub / "x.tsv").is_file(), sub
    for key in SUPPORT_OUTPUTS:
        assert (out / "plots" / "support" / f"{key}.tsv").is_file(), key
    assert not list(out.rglob("variants_qc*")), "work/variants_qc is never published"


def test_local_publish_writes_no_top_level_genetable_context_or_novel(tmp_path: Path) -> None:
    work, out = _work(tmp_path, sorted(MOVED_1_0)), tmp_path / "OUT"
    publish(work, out)
    for legacy in ("genetable", "context", "novel"):
        assert not (out / legacy).exists(), f"OUT/{legacy}/ published at the top level"
        assert (out / "intermediate" / legacy).is_dir(), f"no OUT/intermediate/{legacy}/"


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(st.sets(st.sampled_from(sorted(PUBLISHED_FILES))))
def test_local_publish_top_level_is_always_a_published_name(
    tmp_path_factory: pytest.TempPathFactory, keys: set[str]
) -> None:
    """Property (AC 12.1): whatever subset of step outputs exists, ``publish`` writes only the
    published top-level names, and every moved key lands under ``intermediate/``."""
    root = tmp_path_factory.mktemp("pub")
    work, out = _work(root, sorted(keys)), root / "OUT"
    publish(work, out)
    assert _top(out) <= TOP_LEVEL, sorted(_top(out) - TOP_LEVEL)
    for key in keys & set(MOVED_1_0):
        assert (out / "intermediate" / MOVED_1_0[key]).is_file(), key


def test_published_files_keep_the_nextflow_names() -> None:
    """The ``intermediate/`` prefix is a publish rule, not a rename of the task output names
    (tests/unit/test_local_publish_names.py pins those against ``workflow/``)."""
    for key, rel in MOVED_1_0.items():
        assert PUBLISHED_FILES[key] == rel, key
