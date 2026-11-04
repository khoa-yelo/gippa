"""M12-01 PR 4: ``gippa example fetch`` (docs/specs/M12-01-example-dataset.md section 1, AC 1-12).

The network is mocked: a fake fetcher (the ``Fetcher`` protocol) in place of ``HttpsFetcher``
serves a tiny synthetic bundle built in ``tmp_path``, the selected pin row is monkeypatched to that
archive's url, sha256 and size, ``http.client.HTTPSConnection`` raises, and sockets and DNS are
blocked (``no_network``). The real ``HttpsFetcher`` is driven through the fake ``http.client`` of
``helpers.uniop_14u1`` only (the redirect to the release-asset host, a 404, a timeout).
Seams: see ``tests/helpers/example_m12_01.py``.
"""

from __future__ import annotations

import ast
import errno
import itertools
import os
import re
import shlex
import sys
import tarfile
import tempfile
from collections.abc import Callable, Iterator
from pathlib import Path

import pytest
from helpers.example_m12_01 import (
    ASSET_HOST,
    BOGUS_BUNDLE_COMMAND,
    EXAMPLE_PY,
    Bundle,
    FetchLog,
    Member,
    block_https,
    build_bundle,
    bundle_files,
    checksums_text,
    fill,
    flat,
    gippa,
    install_fetcher,
    install_pins,
    listing,
    load_routes,
    pin_row,
    raise_,
    serve,
    set_installed,
    sha256_hex,
    tree_members,
    url_of,
    verify_tree,
    write_tar,
)
from helpers.uniop_14u1 import Route, install_fake_net
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

pytestmark = pytest.mark.usefixtures("no_network")

SETUP_LINE = "gippa setup --no-genomad"


@pytest.fixture
def world(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> tuple[Path, Bundle]:
    """Installed GIPPA 1.1.0, one pin ``example_1_1`` = a synthetic bundle; ``work`` is the cwd."""
    work = tmp_path.resolve() / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    block_https(monkeypatch)
    bundle = build_bundle(tmp_path.resolve() / "src", "1.1")
    set_installed(monkeypatch, "1.1.0")
    install_pins(monkeypatch, [pin_row(bundle)])
    return work, bundle


def _fetch(*args: str) -> object:
    return gippa(["example", "fetch", *args])


def _assert_fetched(res: object, out_dir: Path, bundle: Bundle) -> None:
    assert res.code == 0, res.err  # type: ignore[attr-defined]
    verify_tree(out_dir, bundle.files)
    lines = [ln for ln in res.out.splitlines() if ln.strip()]  # type: ignore[attr-defined]
    assert lines[0] == f"gippa example {bundle.version}: {out_dir}", lines[0]
    text = flat(res.out)  # type: ignore[attr-defined]
    for rid, route in load_routes().items():
        assert route["title"] in res.out, rid  # type: ignore[attr-defined]
        for cmd in route["commands"]:
            assert fill(cmd, str(out_dir)) in text, (rid, cmd)
    assert SETUP_LINE in lines[-1] and "container" in lines[-1], lines[-1]
    assert BOGUS_BUNDLE_COMMAND not in res.out  # type: ignore[attr-defined]  # yaml, not bundle


def _no_leftovers(work: Path, before: list[str]) -> None:
    assert listing(work) == before, set(listing(work)) ^ set(before)


# ---------------------------------------------------------------- AC 1, 4, 6, 6b: success


def test_fetch_into_absent_outdir(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 1: download into a staging dir next to DIR, verify, unpack, rename; DIR is printed
    absolute although ``--outdir`` is relative."""
    work, bundle = world
    log = install_fetcher(monkeypatch, serve(bundle.data))
    res = _fetch("--outdir", "ex")
    _assert_fetched(res, work / "ex", bundle)
    assert [u for u, _ in log.calls] == [bundle.url]
    staged = log.calls[0][1]
    assert staged.parent.parent == work and staged.parent.name.startswith("ex.partial-"), staged
    assert listing(work) == ["ex", *[f"ex/{p}" for p in listing(work / "ex")]]  # no staging/archive
    assert not list((work / "ex").glob("*.tar.gz"))


def test_default_outdir_is_gippa_example(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    _assert_fetched(_fetch(), work / "gippa-example", bundle)


def test_parents_of_outdir_are_created(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    _assert_fetched(_fetch("--outdir", "a/b/ex"), work / "a" / "b" / "ex", bundle)


def test_from_file_never_touches_the_network(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 4: the fetcher is never constructed; FILE is unchanged (bytes and mtime) and stays."""
    work, bundle = world
    log = install_fetcher(monkeypatch, None)
    hosts = block_https(monkeypatch)
    before = (bundle.path.read_bytes(), bundle.path.stat().st_mtime_ns)
    res = _fetch("--from-file", str(bundle.path), "--outdir", "ex")
    _assert_fetched(res, work / "ex", bundle)
    assert log.constructed == 0 and hosts == []
    assert (bundle.path.read_bytes(), bundle.path.stat().st_mtime_ns) == before


def test_from_file_gives_the_same_tree_as_a_download(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    assert _fetch("--outdir", "net").code == 0
    assert _fetch("--from-file", str(bundle.path), "--outdir", "off").code == 0
    net = {p: (work / "net" / p).read_bytes() for p in listing(work / "net")
           if (work / "net" / p).is_file()}  # fmt: skip
    off = {p: (work / "off" / p).read_bytes() for p in listing(work / "off")
           if (work / "off" / p).is_file()}  # fmt: skip
    assert net == off


def test_existing_empty_outdir_is_filled(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 6."""
    work, bundle = world
    (work / "ex").mkdir()
    install_fetcher(monkeypatch, serve(bundle.data))
    _assert_fetched(_fetch("--outdir", "ex"), work / "ex", bundle)
    assert not [p for p in work.iterdir() if p.name != "ex"]


def test_refetch_is_idempotent_and_offline(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 6b: a previous fetch of the same bundle: ``already fetched``, commands, exit 0, no
    network, nothing changed."""
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    assert _fetch("--outdir", "ex").code == 0
    snap = {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}
    log = install_fetcher(monkeypatch, None)
    res = _fetch("--outdir", "ex")
    assert res.code == 0, res.err
    assert "already fetched" in res.out + res.err
    for cmd in load_routes()["A"]["commands"]:
        assert fill(cmd, str(work / "ex")) in flat(res.out)
    assert log.constructed == 0
    assert snap == {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}


def test_refetch_of_an_altered_fetch_is_refused(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    assert _fetch("--outdir", "ex").code == 0
    (work / "ex" / "D_generic" / "peptides.tsv").write_bytes(b"edited by the user\n")
    install_fetcher(monkeypatch, None)
    res = _fetch("--outdir", "ex")
    assert res.code == 1
    assert "exists and is not this example; use --outdir or remove it" in res.err


def test_shared_options_are_accepted(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Section 1: ``--config``, ``--set`` and ``--log-level`` are accepted and ignored."""
    work, bundle = world
    cfg = work.parent / "c.yaml"
    cfg.write_text("{}\n", encoding="utf-8")
    install_fetcher(monkeypatch, serve(bundle.data))
    res = _fetch("--outdir", "ex", "--log-level", "DEBUG", "--set", "threads=2",
                 "--config", str(cfg))  # fmt: skip
    _assert_fetched(res, work / "ex", bundle)


# ---------------------------------------------------------------- AC 6c, DIR problems: exit 1


@pytest.mark.parametrize("kind", ["file_inside", "is_a_file", "other_bundle"])
def test_outdir_that_is_not_this_example(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """AC 6c: exit 1, nothing written, no network."""
    work, bundle = world
    target = work / "ex"
    if kind == "file_inside":
        target.mkdir()
        (target / "notes.txt").write_text("mine\n", encoding="utf-8")
    elif kind == "is_a_file":
        target.write_text("mine\n", encoding="utf-8")
    else:  # a verified fetch of bundle 1.2 is still "not this example" for 1.1
        other = build_bundle(work.parent / "other", "1.2")
        install_pins(monkeypatch, [pin_row(bundle), pin_row(other)])
        install_fetcher(monkeypatch, None)
        assert _fetch("--from-file", str(other.path), "--version", "1.2",
                      "--outdir", "ex").code == 0  # fmt: skip
    before = {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}
    log = install_fetcher(monkeypatch, serve(bundle.data))
    res = _fetch("--outdir", "ex")
    assert res.code == 1, res.err
    assert f"{target.name} exists and is not this example; use --outdir or remove it" in res.err
    assert log.calls == []
    assert before == {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}


def test_outdir_that_cannot_be_created(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    (work / "plain").write_text("a file, not a directory\n", encoding="utf-8")
    install_fetcher(monkeypatch, serve(bundle.data))
    before = listing(work)
    res = _fetch("--outdir", "plain/ex")
    assert res.code == 1, res.err
    _no_leftovers(work, before)


# ---------------------------------------------------------------- AC 2, 3, 3b, 10: exit 2


def _fails_clean(work: Path, res: object, before: list[str], code: int = 2) -> str:
    assert res.code == code, (res.code, res.err)  # type: ignore[attr-defined]
    _no_leftovers(work, before)
    return res.err  # type: ignore[attr-defined, no-any-return]


def test_download_with_another_sha256(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 2: same size, other bytes: ``checksum mismatch … pinned sha256:<pin>``."""
    work, bundle = world
    data = bytearray(bundle.data)
    data[-9] ^= 0xFF  # inside the gzip trailer (CRC): same size, other sha256
    install_fetcher(monkeypatch, serve(bytes(data)))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "checksum mismatch" in err and f"pinned {bundle.sha256}" in err, err


def test_download_with_another_size(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data + b"\0"))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "size" in err and str(bundle.size) in err, err


@pytest.mark.parametrize("reason, extra", [
    ("timed out", ""),
    ("HTTP 404 for " + url_of("1.1"), "the bundle may not be published yet (GIPPA 1.1 release)"),
])  # fmt: skip
def test_no_network_message(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, reason: str, extra: str
) -> None:
    """AC 3, 3b: the URL, the ``--from-file`` command, exit 2, one attempt, no DIR."""
    from gippa.setupdb.run import SourceError

    work, bundle = world
    log = install_fetcher(monkeypatch, raise_(SourceError("github.com", reason)))
    out_dir = work / "ex"
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", str(out_dir)), before)
    assert f"gippa example: cannot download {bundle.url} (github.com: {reason})" in err, err
    assert "gippa example: on a machine with internet access, download that URL" in err
    assert (f"gippa example fetch --from-file gippa-example-1.1.tar.gz --outdir {out_dir}"
            in flat(err)), err  # fmt: skip
    assert extra in err
    if not extra:
        assert "may not be published yet" not in err
    assert [u for u, _ in log.calls] == [bundle.url]  # no retry, no other host


@pytest.mark.parametrize("code", [errno.ENOSPC, errno.EDQUOT])
def test_disk_full_while_downloading(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, code: int
) -> None:
    """AC 10: exit 2 with the disk-full advice; nothing left."""
    work, _ = world
    install_fetcher(monkeypatch, raise_(OSError(code, os.strerror(code))))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert re.search(r"disk|quota", err, re.I) and "--outdir" in err, err


def test_disk_full_while_unpacking(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Edge case: exit 2; staging removed; DIR not created."""
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))

    def full(*_a: object, **_k: object) -> None:
        raise OSError(errno.ENOSPC, "No space left on device")

    monkeypatch.setattr(tarfile.TarFile, "extractall", full)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert re.search(r"disk|quota", err, re.I), err


def test_failure_in_existing_empty_outdir_leaves_it_empty(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The failure happens after two entries have moved into the user's empty DIR (the third
    move into DIR hits a full disk): DIR is left empty, staging removed, exit 2."""
    work, bundle = world
    target = work / "ex"
    target.mkdir()
    install_fetcher(monkeypatch, serve(bundle.data))
    moved: list[str] = []
    real_rename, real_replace = os.rename, os.replace

    def failing(real: Callable[..., None]) -> Callable[..., None]:
        def move(src: str | os.PathLike[str], dst: str | os.PathLike[str], **k: object) -> None:
            if Path(dst).parent == target:
                if len(moved) == 2:
                    raise OSError(errno.ENOSPC, "No space left on device")
                moved.append(Path(dst).name)
            real(src, dst, **k)

        return move

    monkeypatch.setattr(os, "rename", failing(real_rename))
    monkeypatch.setattr(os, "replace", failing(real_replace))
    before = listing(work)
    assert before == ["ex"]
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert len(moved) == 2, (moved, err)  # the failure came after entries had moved into DIR
    assert list(target.iterdir()) == []
    assert re.search(r"disk|quota", err, re.I), err


def test_created_parents_are_removed_on_failure(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A failed fetch into ``new/deeper/ex`` removes ``new/`` (it created it) and nothing else."""
    work, _ = world
    (work / "keep").mkdir()
    install_fetcher(monkeypatch, serve(b"garbage"))
    before = listing(work)
    assert before == ["keep"]
    _fails_clean(work, _fetch("--outdir", "new/deeper/ex"), before)
    assert not (work / "new").exists() and (work / "keep").is_dir()


def test_unexpected_exception_is_exit_5(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Exit 5 with the standard internal-error message; staging removed."""
    work, _ = world
    install_fetcher(monkeypatch, raise_(RuntimeError("boom")))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before, code=5)
    assert "internal error (a bug)" in err, err


# ---------------------------------------------------------------- AC 8, 8b, 9: bad archives


def _pinned(monkeypatch: pytest.MonkeyPatch, bundle: Bundle) -> FetchLog:
    install_pins(monkeypatch, [pin_row(bundle)])
    return install_fetcher(monkeypatch, serve(bundle.data))


@pytest.mark.parametrize("member", [
    Member("../x", b"escaped\n"),
    Member("/tmp/gippa-m12-01-absolute-member", b"escaped\n"),
    Member("gippa-example-1.1/evil", kind="sym", link="../../x"),
])  # fmt: skip
def test_unsafe_member(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, member: Member
) -> None:
    """AC 8: exit 2 ``unsafe member …``; nothing outside staging; staging removed."""
    work, _ = world
    bad = build_bundle(work.parent / "bad", "1.1", members=[member])
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "unsafe member" in err, err
    assert not (work.parent / "x").exists() and not (work / "x").exists()
    assert not Path("/tmp/gippa-m12-01-absolute-member").exists()


@pytest.mark.parametrize("kind", ["two_dirs", "top_file", "wrong_name"])
def test_malformed_top_level(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """AC 8b: exactly one top-level directory ``gippa-example-<X.Y>/``, else exit 2."""
    work, _ = world
    extra = {"two_dirs": [Member("other", kind="dir"), Member("other/f.txt", b"x\n")],
             "top_file": [Member("stray.txt", b"x\n")], "wrong_name": []}[kind]  # fmt: skip
    top = "gippa-example-1.2" if kind == "wrong_name" else None
    bad = build_bundle(work.parent / "bad", "1.1", members=extra, top=top)
    _pinned(monkeypatch, bad)
    before = listing(work)
    _fails_clean(work, _fetch("--outdir", "ex"), before)


@pytest.mark.parametrize("kind, name", [
    ("altered", "D_generic/peptides.tsv"),
    ("extra", "A_diann/unlisted.tsv"),
    ("missing", "C_cascadia/WT_2.ssl"),
])  # fmt: skip
def test_content_checksums(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, kind: str, name: str
) -> None:
    """AC 9: a member altered after ``CHECKSUMS.sha256``, an unlisted file or a listed file that
    is absent: exit 2 naming the file."""
    work, _ = world
    kw = {"altered": {"alter": {name: b"changed\n"}}, "extra": {"extra": {name: b"x\n"}},
          "missing": {"drop": [name]}}[kind]  # fmt: skip
    bad = build_bundle(work.parent / "bad", "1.1", **kw)  # type: ignore[arg-type]
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert name in err, err


@settings(max_examples=20, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(data=st.data())
def test_property_any_tampered_file_is_named(
    world: tuple[Path, Bundle], data: st.DataObject
) -> None:
    """Property: tampering with any one listed file gives exit 2 naming it, and leaves nothing."""
    work, bundle = world
    name = data.draw(st.sampled_from(sorted(set(bundle.files) - {"CHECKSUMS.sha256"})))
    blob = data.draw(st.binary(min_size=0, max_size=40).filter(lambda b: b != bundle.files[name]))
    with pytest.MonkeyPatch.context() as mp, tempfile.TemporaryDirectory(dir=work.parent) as d:
        bad = build_bundle(Path(d), "1.1", alter={name: blob})
        _pinned(mp, bad)
        before = listing(work)
        err = _fails_clean(work, _fetch("--outdir", "ex"), before)
        assert name in err, (name, err)


# ---------------------------------------------------------------- AC 5, 5b: --from-file input


def test_from_file_of_another_pin(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 5: exit 1 ``this is example bundle 1.2; use --version 1.2``."""
    work, bundle = world
    other = build_bundle(work.parent / "other", "1.2")
    set_installed(monkeypatch, "1.2.0")
    install_pins(monkeypatch, [pin_row(bundle), pin_row(other)])
    install_fetcher(monkeypatch, None)
    before = listing(work)
    res = _fetch("--from-file", str(other.path), "--version", "1.1", "--outdir", "ex")
    err = _fails_clean(work, res, before, code=1)
    assert "this is example bundle 1.2; use --version 1.2" in err, err


@pytest.mark.parametrize("kind", ["empty", "directory", "missing", "wrong_bytes"])
def test_from_file_bad_input(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """AC 5b: exit 1 (empty or other bytes: a size or sha256 message; missing: Click's error)."""
    work, bundle = world
    f = work.parent / "in.tar.gz"
    if kind == "empty":
        f.write_bytes(b"")
    elif kind == "directory":
        f.mkdir()
    elif kind == "wrong_bytes":
        data = bytearray(bundle.data)
        data[-9] ^= 0xFF
        f.write_bytes(bytes(data))
    install_fetcher(monkeypatch, None)
    before = listing(work)
    err = _fails_clean(work, _fetch("--from-file", str(f), "--outdir", "ex"), before, code=1)
    if kind == "empty":
        assert "size" in err, err
    if kind == "wrong_bytes":
        assert "sha256" in err or "checksum" in err, err


# ---------------------------------------------------------------- AC 7-7e: version selection


def _pins(monkeypatch: pytest.MonkeyPatch, *versions: str) -> None:
    install_pins(monkeypatch, [pin_row(version=v) for v in versions])


def _asked(monkeypatch: pytest.MonkeyPatch, *args: str) -> tuple[object, list[str]]:
    """Run fetch with a fetcher that records the URL and fails (the selection is in the URL)."""
    from gippa.setupdb.run import SourceError

    log = install_fetcher(monkeypatch, raise_(SourceError("github.com", "timed out")))
    res = _fetch("--outdir", "ex", *args)
    return res, [u for u, _ in log.calls]


def test_unknown_version_lists_the_pins(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 7."""
    _pins(monkeypatch, "1.0", "1.1")
    res, urls = _asked(monkeypatch, "--version", "9.9")
    assert res.code == 1 and urls == []  # type: ignore[attr-defined]
    assert "1.0" in res.err and "1.1" in res.err  # type: ignore[attr-defined]


@pytest.mark.parametrize("bad", ["abc", "1", "1.x", ""])
def test_malformed_version_is_exit_1(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, bad: str
) -> None:
    res, urls = _asked(monkeypatch, "--version", bad)
    assert res.code == 1 and urls == []  # type: ignore[attr-defined]


def test_older_pin_is_used_with_a_warning(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 7b: installed 1.2.0, only 1.1 pinned: 1.1, the warning on stderr, exit 0."""
    work, bundle = world
    set_installed(monkeypatch, "1.2.0")
    install_fetcher(monkeypatch, serve(bundle.data))
    res = _fetch("--outdir", "ex")
    _assert_fetched(res, work / "ex", bundle)
    warn = (r"gippa example: no bundle for 1\.2(\.0)?; using 1\.1 \(expected numbers were made "
            r"with GIPPA 1\.1 and may differ, see Upgrade\)")  # fmt: skip
    assert re.search(warn, flat(res.err)), res.err
    assert "no bundle for" not in res.out


def test_explicit_older_version_warns_too(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    _pins(monkeypatch, "1.1", "1.2")
    set_installed(monkeypatch, "1.2.0")
    res, urls = _asked(monkeypatch, "--version", "1.1")
    assert urls == [url_of("1.1")]
    assert "no bundle for" in res.err and "using 1.1" in res.err  # type: ignore[attr-defined]


@pytest.mark.parametrize("installed", ["2.0.0", "2.3.1"])
def test_other_major_is_exit_2(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, installed: str
) -> None:
    """AC 7c."""
    set_installed(monkeypatch, installed)
    res, urls = _asked(monkeypatch)
    assert res.code == 2 and urls == []  # type: ignore[attr-defined]
    assert "no example bundle for GIPPA 2.x" in res.err  # type: ignore[attr-defined]


def test_no_pins_at_all_is_exit_2(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Before the pin PR: a 1.x install exits 2 ``no example bundle for GIPPA 1.x``."""
    install_pins(monkeypatch, [])
    res, urls = _asked(monkeypatch)
    assert res.code == 2 and urls == []  # type: ignore[attr-defined]
    assert "no example bundle for GIPPA 1.x" in res.err  # type: ignore[attr-defined]


@pytest.mark.parametrize("installed", ["1.1.0rc1", "1.1.0.dev3", "1.1.0", "1.1.7"])
def test_dev_rc_and_patch_count_as_their_minor(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, installed: str
) -> None:
    """AC 7d: uses 1.1, no warning."""
    work, bundle = world
    set_installed(monkeypatch, installed)
    install_fetcher(monkeypatch, serve(bundle.data))
    res = _fetch("--outdir", "ex")
    _assert_fetched(res, work / "ex", bundle)
    assert "no bundle for" not in res.err


def test_patch_version_option_is_read_as_minor(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """AC 7e: ``--version 1.1.4`` is ``1.1``."""
    _pins(monkeypatch, "1.0", "1.1")
    res, urls = _asked(monkeypatch, "--version", "1.1.4")
    assert urls == [url_of("1.1")], res.err  # type: ignore[attr-defined]


_VERSIONS = st.tuples(st.integers(1, 2), st.integers(0, 5))


@settings(max_examples=40, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(pins=st.sets(_VERSIONS, max_size=5), inst=st.tuples(st.integers(1, 3), st.integers(0, 6)),
       patch=st.integers(0, 3), suffix=st.sampled_from(["", "rc1", ".dev3"]))  # fmt: skip
def test_property_default_selection(
    world: tuple[Path, Bundle],
    pins: set[tuple[int, int]],
    inst: tuple[int, int],
    patch: int,
    suffix: str,
) -> None:
    """Property (section 1 step 1): the installed minor; else the newest older minor of the same
    major with the warning; else exit 2 ``no example bundle for GIPPA <MAJOR>.x``."""
    same = sorted(m for M, m in pins if inst[0] == M)
    want = inst[1] if inst[1] in same else max((m for m in same if m < inst[1]), default=None)
    with pytest.MonkeyPatch.context() as mp:
        set_installed(mp, f"{inst[0]}.{inst[1]}.{patch}{suffix}")
        _pins(mp, *(f"{M}.{m}" for M, m in sorted(pins)))
        res, urls = _asked(mp)
    err = res.err  # type: ignore[attr-defined]
    if want is None:
        assert res.code == 2 and urls == [], err  # type: ignore[attr-defined]
        assert f"no example bundle for GIPPA {inst[0]}.x" in err
    else:
        assert urls == [url_of(f"{inst[0]}.{want}")], err
        assert ("no bundle for" in err) == (want != inst[1]), err


@settings(max_examples=30, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(pins=st.sets(st.integers(0, 5), min_size=1, max_size=4), inst=st.integers(0, 5),
       req=st.integers(0, 6))  # fmt: skip
def test_property_explicit_version(
    world: tuple[Path, Bundle], pins: set[int], inst: int, req: int
) -> None:
    """Property: ``--version 1.<r>`` takes that pin (warning iff older than the installed minor),
    else exit 1 listing every pinned version."""
    with pytest.MonkeyPatch.context() as mp:
        set_installed(mp, f"1.{inst}.0")
        _pins(mp, *(f"1.{m}" for m in sorted(pins)))
        res, urls = _asked(mp, "--version", f"1.{req}")
    err = res.err  # type: ignore[attr-defined]
    if req in pins:
        assert urls == [url_of(f"1.{req}")], err
        assert ("no bundle for" in err) == (req < inst), err
    else:
        assert res.code == 1 and urls == [], err  # type: ignore[attr-defined]
        assert all(f"1.{m}" in err for m in pins), err


# ---------------------------------------------------------------- the real HttpsFetcher, faked net


def _real_fetcher_world(world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch) -> object:
    from gippa.setupdb.run import HttpsFetcher

    for mod in ("gippa.setupdb.run", "gippa.setupdb", "gippa.example"):
        if mod in sys.modules:
            monkeypatch.setattr(sys.modules[mod], "HttpsFetcher", HttpsFetcher, raising=False)
    return install_fake_net(monkeypatch)


def test_real_fetcher_follows_the_release_asset_redirect(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The pin row's ``redirect_hosts`` let ``HttpsFetcher`` follow GitHub's 302 to the asset
    host, unchanged (section 1 > Pin)."""
    work, bundle = world
    net = _real_fetcher_world(world, monkeypatch)
    asset = f"https://{ASSET_HOST}/github-production-release-asset/1/gippa-example-1.1"
    path = bundle.url.removeprefix("https://github.com")
    net.routes[("github.com", path)] = Route(302, {"Location": asset})  # type: ignore[attr-defined]
    net.routes[(ASSET_HOST, asset.split(ASSET_HOST)[1])] = Route(200, {}, bundle.data)  # type: ignore[attr-defined]
    _assert_fetched(_fetch("--outdir", "ex"), work / "ex", bundle)
    assert [h for h, _ in net.contacted] == ["github.com", ASSET_HOST]  # type: ignore[attr-defined]


def test_real_fetcher_404_says_not_published(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, _ = world
    _real_fetcher_world(world, monkeypatch)  # every unknown route is a 404
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "may not be published yet" in err and "--from-file" in err, err


# ---------------------------------------------------------------- AC 11: reuse, no second copy


def _imports(path: Path) -> tuple[set[str], set[tuple[str, str]]]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    mods: set[str] = set()
    names: set[tuple[str, str]] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods |= {a.name for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module:
            mods.add(node.module)
            names |= {(node.module, a.name) for a in node.names}
    return mods, names


def test_example_module_reuses_the_setup_downloader() -> None:
    """AC 11: ``src/gippa/example.py`` imports ``HttpsFetcher``, ``file_digest`` and the safe
    unpack from ``gippa.setupdb``, and no HTTPS or archive library itself."""
    assert EXAMPLE_PY.is_file(), "src/gippa/example.py does not exist"
    mods, names = _imports(EXAMPLE_PY)
    banned = {"http.client", "http", "urllib.request", "tarfile", "zipfile", "requests", "urllib3"}
    assert not mods & banned, mods & banned
    got = {n for m, n in names if m == "gippa.setupdb" or m.startswith("gippa.setupdb.")}
    assert "HttpsFetcher" in got and "file_digest" in got, got
    assert got & {"unpack_archive", "_unpack"}, got


def test_unpack_archive_is_the_public_alias() -> None:
    """Q-L3 (default yes): ``unpack_archive`` is public and ``_unpack`` stays the same function."""
    from gippa.setupdb import run

    assert run.unpack_archive is run._unpack
    assert "unpack_archive" in run.__all__


# ---------------------------------------------------------------- AC 12: setup never installs it


def test_setup_plan_never_holds_an_example_row(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """AC 12: ``items.plan()`` for every ``SetupOptions`` combination has no ``example_*`` item."""
    from gippa.setupdb import items

    install_pins(monkeypatch, [pin_row(version="1.1"), pin_row(version="1.0")])
    flags = ("genomad", "antismash", "promotech_models", "check", "go", "kegg", "cog",
             "contaminants", "uniop", "minimal", "nextflow")  # fmt: skip
    for bakta, values in itertools.product(
        ("light", "full", "skip"), itertools.product((False, True), repeat=len(flags))
    ):
        opts = items.SetupOptions(tmp_path, bakta, **dict(zip(flags, values, strict=True)))  # type: ignore[arg-type]
        assert not [i.name for i in items.plan(opts) if i.name.startswith("example_")], opts


# ---------------------------------------------------------------- review #542: links, duplicates


def _checksums_plus(version: str, *lines: str) -> bytes:
    files = bundle_files(version)
    return checksums_text(files) + "".join(f"{ln}\n" for ln in lines).encode()


def test_listed_in_tree_symlink_is_refused(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """A symlink listed in ``CHECKSUMS.sha256`` that points at a file with the listed content is
    still not a regular file: exit 2, nothing left."""
    work, _ = world
    data = bundle_files("1.1")["D_generic/peptides.tsv"]
    sums = _checksums_plus("1.1", f"{sha256_hex(data)}  D_generic/link.tsv")
    link = Member("gippa-example-1.1/D_generic/link.tsv", kind="sym", link="peptides.tsv")
    bad = build_bundle(work.parent / "bad", "1.1", alter={"CHECKSUMS.sha256": sums},
                       members=[link])  # fmt: skip
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "link.tsv" in err, err


def test_unlisted_in_tree_symlink_is_refused(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, _ = world
    link = Member("gippa-example-1.1/A_diann/alias.parquet", kind="sym", link="report.parquet")
    bad = build_bundle(work.parent / "bad", "1.1", members=[link])
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "alias.parquet" in err, err


@pytest.mark.parametrize("with_files", [False, True])
def test_top_level_symlink_is_refused(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, with_files: bool
) -> None:
    """``gippa-example-1.1 -> .`` (alone, or with the bundle files at the archive root, so that
    following the link would show a valid tree) is a malformed archive: exit 2."""
    work, _ = world
    files = bundle_files("1.1")
    listed = {**files, "CHECKSUMS.sha256": checksums_text(files)}
    members = [Member("gippa-example-1.1", kind="sym", link=".")]
    if with_files:
        members += [m for m in tree_members(".", listed) if m.name != "."]
        members = [Member(m.name.removeprefix("./"), m.data, m.kind, m.link) for m in members]
    path = write_tar(work.parent / "bad" / "gippa-example-1.1.tar.gz", members)
    data = path.read_bytes()
    bad = Bundle("1.1", path, f"sha256:{sha256_hex(data)}", len(data), listed)
    _pinned(monkeypatch, bad)
    before = listing(work)
    _fails_clean(work, _fetch("--outdir", "ex"), before)


def test_duplicate_name_in_checksums_is_refused(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, _ = world
    data = bundle_files("1.1")["README.md"]
    sums = _checksums_plus("1.1", f"{sha256_hex(data)}  README.md")
    bad = build_bundle(work.parent / "bad", "1.1", alter={"CHECKSUMS.sha256": sums})
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "CHECKSUMS.sha256" in err, err


# ---------------------------------------------------------------- review #542 M1: runs/ in DIR


def test_refetch_with_route_outputs_in_runs_is_a_noop(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """M1: the printed routes write to ``{dir}/runs``; a re-fetch after a route ran is still
    ``already fetched`` (exit 0, offline) and leaves ``runs/`` untouched."""
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    assert _fetch("--outdir", "ex").code == 0
    out = work / "ex" / "runs" / "A" / "peptides.tsv"
    out.parent.mkdir(parents=True)
    out.write_bytes(b"peptide\trun\nPEPTIDE\tWT_1\n")
    snap = {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}
    log = install_fetcher(monkeypatch, None)
    hosts = block_https(monkeypatch)
    res = _fetch("--outdir", "ex")
    assert res.code == 0, res.err
    assert "already fetched" in res.out + res.err
    assert log.constructed == 0 and hosts == []
    assert snap == {p: (work / p).read_bytes() for p in listing(work) if (work / p).is_file()}


def test_refetch_with_runs_still_refuses_an_edited_listed_file(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    assert _fetch("--outdir", "ex").code == 0
    (work / "ex" / "runs" / "A").mkdir(parents=True)
    (work / "ex" / "runs" / "A" / "peptides.tsv").write_bytes(b"peptide\trun\n")
    (work / "ex" / "D_generic" / "peptides.tsv").write_bytes(b"edited by the user\n")
    install_fetcher(monkeypatch, None)
    res = _fetch("--outdir", "ex")
    assert res.code == 1, res.err
    assert "exists and is not this example; use --outdir or remove it" in res.err


# ---------------------------------------------------------------- review #542 M2: not writable


@pytest.fixture
def read_only() -> Iterator[Callable[[Path], Path]]:
    """``chmod 555`` on a directory, restored at teardown; skipped as root."""
    if hasattr(os, "geteuid") and os.geteuid() == 0:
        pytest.skip("root ignores directory permissions")
    done: list[Path] = []

    def make(d: Path) -> Path:
        d.chmod(0o555)
        done.append(d)
        return d

    yield make
    for d in done:
        d.chmod(0o755)


def test_outdir_parent_not_writable_is_exit_1(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, read_only: Callable[[Path], Path]
) -> None:
    """M2: the parent of ``--outdir`` exists but cannot be written: exit 1, "choose another
    --outdir", nothing left (no staging, no DIR), no download."""
    work, bundle = world
    (work / "ro").mkdir()
    read_only(work / "ro")
    log = install_fetcher(monkeypatch, serve(bundle.data))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ro/ex"), before, code=1)
    assert "choose another --outdir" in err, err
    assert log.calls == []


def test_existing_empty_outdir_not_writable_is_exit_1(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, read_only: Callable[[Path], Path]
) -> None:
    work, bundle = world
    (work / "ex").mkdir()
    read_only(work / "ex")
    log = install_fetcher(monkeypatch, serve(bundle.data))
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before, code=1)
    assert "choose another --outdir" in err, err
    assert list((work / "ex").iterdir()) == []
    assert log.calls == []


# ---------------------------------------------------------------- review #542 L4, L5: quoting


_ODD = "my  ex\tdir"  # two consecutive spaces and a tab
_TOKEN = "@@GIPPA-DIR@@"


def test_printed_commands_round_trip_through_shlex(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """L4: each printed command, split by ``shlex``, is the argv of its ``example_routes.yaml``
    command with DIR (two spaces and a tab in its name) as one unchanged word part."""
    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))
    out_dir = work / _ODD
    res = _fetch("--outdir", _ODD)
    assert res.code == 0, res.err
    verify_tree(out_dir, bundle.files)
    printed = [ln for ln in res.out.splitlines()[1:] if ln.startswith("gippa ")]
    want = [[w.replace(_TOKEN, str(out_dir)) for w in shlex.split(fill(cmd, _TOKEN))]
            for route in load_routes().values() for cmd in route["commands"]]  # fmt: skip
    assert [shlex.split(ln) for ln in printed] == want


def test_no_network_hint_quotes_the_outdir(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """L5: the ``--from-file`` command of the no-network message is shell-quoted."""
    from gippa.setupdb.run import SourceError

    work, _ = world
    install_fetcher(monkeypatch, raise_(SourceError("github.com", "timed out")))
    out_dir = work / _ODD
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", str(out_dir)), before)
    start = err.index("gippa example fetch --from-file")
    hint = err[start:].splitlines()[0]
    want = ["gippa", "example", "fetch", "--from-file", "gippa-example-1.1.tar.gz"]
    assert shlex.split(hint) == [*want, "--outdir", str(out_dir)], hint


# ---------------------------------------------------------------- review #542 L1: bugs are exit 5


def test_programming_error_in_unpack_is_exit_5(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """L1: a ``TypeError`` raised inside the unpacker is a bug (exit 5, the internal-error
    message), not a "malformed archive" (exit 2)."""
    import gippa.example

    work, bundle = world
    install_fetcher(monkeypatch, serve(bundle.data))

    def bug(*_a: object, **_k: object) -> None:
        raise TypeError("unpack_archive() got an unexpected keyword argument")

    monkeypatch.setattr(gippa.example, "unpack_archive", bug)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before, code=5)
    assert "malformed archive" not in err, err
    assert "internal error (a bug)" in err, err


@pytest.mark.parametrize("kind", ["truncated", "not_a_tar"])
def test_corrupt_archive_is_still_exit_2(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch, kind: str
) -> None:
    """L1: a real corrupt ``.tar.gz`` whose bytes match the pin is a malformed archive, exit 2."""
    import gzip

    work, bundle = world
    data = (bundle.data[: len(bundle.data) // 2] if kind == "truncated"
            else gzip.compress(b"this is not a tar archive\n" * 40, mtime=0))  # fmt: skip
    path = work.parent / "corrupt" / "gippa-example-1.1.tar.gz"
    path.parent.mkdir()
    path.write_bytes(data)
    bad = Bundle("1.1", path, f"sha256:{sha256_hex(data)}", len(data), {})
    _pinned(monkeypatch, bad)
    before = listing(work)
    err = _fails_clean(work, _fetch("--outdir", "ex"), before)
    assert "internal error" not in err, err


# ---------------------------------------------------------------- review #542 L3: other major


def test_explicit_version_of_another_major_is_exit_2(
    world: tuple[Path, Bundle], monkeypatch: pytest.MonkeyPatch
) -> None:
    """L3, spec section 1 > Pin ("every older bundle row it can still run (same major)") and
    section 6 ("another major: exit 2"): ``--version 1.1`` on GIPPA 2.0 is exit 2, without the
    same-major fallback warning and without a download."""
    _pins(monkeypatch, "1.1")
    set_installed(monkeypatch, "2.0.0")
    res, urls = _asked(monkeypatch, "--version", "1.1")
    err = res.err  # type: ignore[attr-defined]
    assert res.code == 2 and urls == [], err  # type: ignore[attr-defined]
    assert "no bundle for" not in err and "using 1.1" not in err, err
    assert "GIPPA 2" in err, err


def test_select_pin_never_warns_across_majors() -> None:
    """L3 (pure): ``select_pin`` with installed 2.0 and ``--version 1.1`` raises exit 2 and never
    builds the "no bundle for 2.0; using 1.1" warning."""
    from gippa.example import ExampleError, select_pin

    row = pin_row(version="1.1")
    with pytest.raises(ExampleError) as info:
        select_pin({(1, 1): row}, "2.0.0", (1, 1))
    assert int(info.value.code) == 2
    assert not any("using 1.1" in ln for ln in info.value.lines)
