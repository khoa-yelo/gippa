"""M12-01 PR 4: the bundle pin (docs/specs/M12-01-example-dataset.md section 1 > Pin, AC 13).

- ``DbItem.size_bytes`` is a new optional key of the ``db_sources.yaml`` loader (default 0,
  ignored by setup); the loader rejects two ``example_*`` rows with the same ``version``.
- Every ``example_<X>_<Y>`` row of ``db_sources.yaml`` has the shape the spec shows, and exactly
  one row of ``validation/example/bundles.tsv`` with the same version, sha256, size and URL, and
  vice versa (CI-safe, no network). With no bundle rows yet the row checks pass on the empty set;
  the pin PR (PR 7) adds the "current minor has a row" assertion.
"""

from __future__ import annotations

import csv
import dataclasses
import re
from pathlib import Path

import pytest
import yaml
from helpers.example_m12_01 import ASSET_HOST, REPO, url_of

SOURCES_YAML = REPO / "src" / "gippa" / "data" / "db_sources.yaml"
BUNDLES_TSV = REPO / "validation" / "example" / "bundles.tsv"
COLUMNS = ["version", "file", "size_bytes", "sha256", "url", "release_tag", "build_commit",
           "build_job", "expected_commit", "asset_digest"]  # fmt: skip
HEX64 = re.compile(r"[0-9a-f]{64}")


def _raw_rows() -> dict[str, dict[str, object]]:
    raw = yaml.safe_load(SOURCES_YAML.read_text(encoding="utf-8"))
    return {k: v for k, v in raw.items() if k.startswith("example_")}


def _bundles() -> list[dict[str, str]]:
    with BUNDLES_TSV.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _hex(digest: str) -> str:
    return digest.removeprefix("sha256:")


# ---------------------------------------------------------------- the loader


def test_dbitem_has_optional_size_bytes() -> None:
    from gippa.setupdb.items import DbItem

    fields = {f.name: f for f in dataclasses.fields(DbItem)}
    assert "size_bytes" in fields
    assert fields["size_bytes"].default == 0


def _load(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, extra: str) -> object:
    from gippa.setupdb import items

    d = tmp_path / "data"
    d.mkdir()
    (d / "db_sources.yaml").write_text(SOURCES_YAML.read_text(encoding="utf-8") + extra,
                                       encoding="utf-8")  # fmt: skip
    monkeypatch.setattr(items, "files", lambda _pkg: d)
    return items._load_sources()


def _row(version: str, size: str = "") -> str:
    x, y = version.split(".")
    return (f"example_{x}_{y}:\n  version: \"{version}\"\n  url: {url_of(version)}\n"
            f"  checksum: sha256:{'a' * 64}\n{size}  path: example/gippa-example-{version}.tar.gz\n"
            f"  tool: null\n  check: sha256\n  size_gib: 1\n"
            f"  redirect_hosts: [{ASSET_HOST}]\n")  # fmt: skip


def test_loader_reads_size_bytes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    sources = _load(tmp_path, monkeypatch, _row("1.1", "  size_bytes: 480123\n"))
    row = sources["example_1_1"]  # type: ignore[index]
    assert row.size_bytes == 480123
    assert row.redirect_hosts == (ASSET_HOST,)
    assert all(r.size_bytes == 0 for k, r in sources.items() if k != "example_1_1")  # type: ignore[attr-defined]


def test_loader_size_bytes_defaults_to_zero(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    sources = _load(tmp_path, monkeypatch, _row("1.1"))
    assert sources["example_1_1"].size_bytes == 0  # type: ignore[index]


def test_loader_rejects_two_example_rows_with_one_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Section 6: two pins with an equal version are a ValueError at load time."""
    twin = _row("1.1", "  size_bytes: 1\n").replace("example_1_1:", "example_1_1_again:")
    with pytest.raises(ValueError, match=r"1\.1"):
        _load(tmp_path, monkeypatch, _row("1.1", "  size_bytes: 1\n") + twin)


def test_loader_still_accepts_equal_versions_of_other_rows() -> None:
    """``bakta_light`` and ``bakta_full`` share version 6.0: the rule is for example rows only."""
    from gippa.setupdb.items import SOURCES

    assert SOURCES["bakta_light"].version == SOURCES["bakta_full"].version


# ---------------------------------------------------------------- the pin rows (AC 13)


@pytest.mark.parametrize("name", sorted(_raw_rows()))
def test_pin_row_shape(name: str) -> None:
    from gippa.setupdb.items import SOURCES

    row = _raw_rows()[name]
    v = str(row["version"])
    assert re.fullmatch(r"\d+\.\d+", v) and name == "example_" + v.replace(".", "_"), (name, v)
    assert row["url"] == url_of(v)
    assert str(row["checksum"]).startswith("sha256:") and HEX64.fullmatch(
        _hex(str(row["checksum"]))
    )
    assert row["check"] == "sha256" and row["tool"] is None
    assert row["path"] == f"example/gippa-example-{v}.tar.gz"
    assert row["redirect_hosts"] == [ASSET_HOST]
    assert isinstance(row["size_bytes"], int) and row["size_bytes"] > 0
    assert SOURCES[name].size_bytes == row["size_bytes"]


def test_bundles_tsv_exists_with_its_columns() -> None:
    """``validation/example/bundles.tsv`` (the release asset manifest; PR 3 adds the header)."""
    assert BUNDLES_TSV.is_file(), "validation/example/bundles.tsv is missing"
    head = BUNDLES_TSV.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert head == COLUMNS


def test_pins_and_bundles_tsv_are_equal_sets() -> None:
    pins = {(str(r["version"]), _hex(str(r["checksum"])), int(r["size_bytes"]), str(r["url"]))  # type: ignore[call-overload]
            for r in _raw_rows().values()}  # fmt: skip
    rows = _bundles()
    tsv = [(b["version"], _hex(b["sha256"]), int(b["size_bytes"]), b["url"]) for b in rows]
    assert len(tsv) == len(set(tsv)), "duplicate bundles.tsv rows"
    assert set(tsv) == pins


def test_bundles_tsv_rows_are_consistent() -> None:
    for b in _bundles():
        v = b["version"]
        assert b["file"] == f"gippa-example-{v}.tar.gz", b
        assert b["url"] == url_of(v) and b["release_tag"] == f"v{v}.0", b
        assert HEX64.fullmatch(_hex(b["sha256"])), b
        if b["asset_digest"]:
            assert _hex(b["asset_digest"]) == _hex(b["sha256"]), b
