"""Installed ``gippa setup`` items that pipeline steps read (GIPPA_SPEC.md §7.9.1)."""

from __future__ import annotations

from pathlib import Path

__all__ = ["FUNCTION_RESOURCES", "KEGG_RESOURCES", "cog_table", "function_resources", "resource"]

#: the setup items the ``function`` step reads (M8-06 item 7), in ``read_resources`` key form
FUNCTION_RESOURCES: tuple[str, ...] = ("go_basic", "goslim_generic", "kegg_ko_pathway",
                                       "kegg_pathway_names", "cog2020_fun")  # fmt: skip


#: KEGG counts as installed only when both items are (M12-03); a lone item is not staged
KEGG_RESOURCES: tuple[str, ...] = ("kegg_ko_pathway", "kegg_pathway_names")


def function_resources(db_dir: Path | None) -> dict[str, Path]:
    """The installed ``FUNCTION_RESOURCES`` to stage, without a lone KEGG item (M12-03)."""
    found = {n: p for n in FUNCTION_RESOURCES if (p := resource(db_dir, n)) is not None}
    if not all(n in found for n in KEGG_RESOURCES):
        found = {n: p for n, p in found.items() if n not in KEGG_RESOURCES}
    return found


def resource(db_dir: Path | None, name: str) -> Path | None:
    """``db_dir / SOURCES[name].path`` when ``db_dir``'s setup manifest lists ``name`` and the
    file exists (M4-15 decision 5, generalised in M8-06 item 7), else ``None``."""
    from gippa.preflight.environment import read_manifest
    from gippa.setupdb.items import SOURCES

    if db_dir is None:
        return None
    path = db_dir / SOURCES[name].path
    return path if name in (read_manifest(str(db_dir)) or {}) and path.is_file() else None


def cog_table(db_dir: Path | None) -> Path | None:
    """The COG 2020 ``cog-20.def.tab`` (``db_sources.yaml`` ``cog2020_def``) when installed
    (M4-15 decisions 3, 17, 22), else ``None``."""
    return resource(db_dir, "cog2020_def")
