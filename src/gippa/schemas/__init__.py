"""JSON Schemas shipped with GIPPA: run manifest and resolved config.

``OUTPUT_SCHEMA_VERSION`` (M11-14h; ADR-0060 §6) is the ``MAJOR.MINOR`` version of the run's
outputs (the tables and the manifest), written to ``run_manifest.json`` and carried by
``datapackage.json``. A new column, table or manifest key bumps the minor; a removed, renamed or
retyped one bumps the major. It is apart from ``schema_version``, the config format (Q-L13).
``1.0``: the first published output schema (spec §4, Q-L13, ADR-0060 §6). ``1.1``: the
manifest's ``annotation.translation_table`` (M11-16). ``1.2``: the
``peptide_hits`` flag ``initiator_m`` (M11-20). ``1.3``: the 1.1 outputs revision, which keeps
major 1 by a one-time exception (ADR-0075, ``scripts/release/schema_exceptions.json``): the
``genes.tsv`` column groups, the ``gippa:run_level`` field property and ``detectability_constant``
as a run flag (M12-04 PR 11a), and the column renames (M12-04 PR 11b).
"""

from __future__ import annotations

__all__ = ["OUTPUT_SCHEMA_VERSION"]

#: keep equal to ``output_schema_version`` in ``datapackage.json``
OUTPUT_SCHEMA_VERSION = "1.3"
