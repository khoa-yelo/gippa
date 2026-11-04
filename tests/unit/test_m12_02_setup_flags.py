"""M12-02 (GIPPA 1.1.0): ``gippa setup --accept-download-plan`` and ``--install-uniop-from-git``.

Khoa 2026-11-03: "be more specific than --yes, and have sth like --install-uniop-from-git", and
"explain the flags for user as well". Goal: ONE non-interactive full setup,
``gippa setup --db-dir D --profile apptainer --accept-download-plan --install-uniop-from-git``
installs the databases, pulls the images and builds UniOP, with no prompt of any kind.

Contract pinned here (the lead's revised brief):

1. ``--accept-download-plan`` replaces ``--yes``: it skips the ``Proceed?`` confirmation of the
   download plan and nothing else. It does not accept the UniOP terms (``setup_uniop`` still gets
   ``accepted=False``; on a TTY the UniOP prompt still comes).
2. ``--install-uniop-from-git`` replaces ``--accept-uniop-terms`` (and only that flag): it accepts
   the UniOP terms, so setup never prompts for them. It does not imply ``--uniop``: the databases
   are planned and downloaded, the images pulled, and UniOP is built when the config needs it
   (``setup_uniop(..., accepted=True, explicit=False)``). The UniOP build fetches from the authors'
   GitHub repository at the pinned commit and checks the sha256; the one-line terms notice is still
   printed once, before the fetch.
3. ``--uniop`` is not deprecated: visible in ``--help``, "install UniOP only" (Q-K27), no
   deprecation line. ``--uniop --install-uniop-from-git`` is the old ``--uniop
   --accept-uniop-terms``.
4. Usage errors (exit 2, the message names both flags, nothing is set up):
   ``--install-uniop-from-git --no-uniop``; ``--install-uniop-from-git --no-images`` without
   ``--uniop`` (with ``--uniop`` it is allowed).
5. Deprecated, hidden from ``--help``, still working, one deprecation line each (a line containing
   "deprecat" naming the old and the new flag): ``--yes`` -> ``--accept-download-plan``,
   ``--accept-uniop-terms`` -> ``--install-uniop-from-git``. ``--uniop`` without the terms: the
   prompt on a TTY; without one, exit 2 with "terms not accepted" naming
   ``--install-uniop-from-git``.
6. ``--help`` explains both new flags (key phrases, not exact text).
7. Preflight's hint: ``gippa setup --uniop --install-uniop-from-git --profile <p> --db-dir <d>``.
8. ``--check`` and ``--minimal`` behave as before (``--minimal`` keeps the ``core,function``
   images; UniOP is ``core``, so its step runs, without consent unless the flag is given).

Two kinds of tests: flag-mapping tests record what ``setup`` hands to ``install_images``,
``run_setup`` and ``setup_uniop`` (all faked, as in test_m11_15c2_setup_defaults.py); behaviour
tests run the real UniOP step against the fake network, docker and archive of the M11-14u1/u2
helpers. No test touches the network, a real tool, ``$HOME`` or ``$GIPPA_DB_DIR``.
"""

from __future__ import annotations

import dataclasses
import os
import re
import shutil
import socket
import subprocess
from collections import namedtuple
from pathlib import Path
from typing import Any

import pytest
from click.testing import CliRunner
from helpers.uniop_14u1 import (
    PKG,
    SKIP_DBS,
    FakeNet,
    docker_calls,
    fixture_archive,
    install_fake_net,
    isolate,
    manifest_items,
    pin_uniop_source,
    sha256_hex,
)
from helpers.uniop_14u2 import NOTICE, PROMPT, install_fake_runtimes
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

NEW_PLAN = "--accept-download-plan"
NEW_UNIOP = "--install-uniop-from-git"
OLD_YES = "--yes"
UNIOP = "--uniop"
OLD_TERMS = "--accept-uniop-terms"
DEPRECATED = (OLD_YES, OLD_TERMS)
PROCEED = "Proceed?"
TERMS_NOT_ACCEPTED = "terms not accepted"
MARK = "<<FIRST-FETCH>>"

_Usage = namedtuple("_Usage", "total used free")
HUGE = _Usage(1 << 50, 0, 1 << 50)

pytestmark = pytest.mark.usefixtures("no_network")


def _flag(name: str) -> re.Pattern[str]:
    """``name`` as a whole option token (``--uniop`` does not match ``--no-uniop``)."""
    return re.compile(rf"(?<![\w-]){re.escape(name)}(?![\w-])")


def _deprecations(out: str) -> list[str]:
    return [ln for ln in out.splitlines() if "deprecat" in ln.lower()]


def _tty(monkeypatch: pytest.MonkeyPatch, is_tty: bool = True) -> None:
    """A terminal on stdin for both the plan prompt and the UniOP prompt (CliRunner has none)."""
    from gippa.setupdb import space, uniop_build

    monkeypatch.setattr(uniop_build, "stdin_is_tty", lambda: is_tty)
    monkeypatch.setattr(space, "stdin_is_tty", lambda: is_tty, raising=False)


def _invoke(args: list[str], answer: str | None = None) -> tuple[int, str]:
    """``gippa setup <args>``; returns the exit code and the combined stdout + stderr."""
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["setup", *args], input=answer)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res.exit_code, res.output


# ======================================================== flag mapping (all side effects faked)


@dataclasses.dataclass
class Calls:
    installed: list[str] = dataclasses.field(default_factory=list)
    images: list[list[str]] = dataclasses.field(default_factory=list)
    uniop: list[dict[str, Any]] = dataclasses.field(default_factory=list)
    planned: list[set[str]] = dataclasses.field(default_factory=list)
    checks: list[bool] = dataclasses.field(default_factory=list)

    def state(self) -> tuple[Any, ...]:
        return (list(self.uniop), list(self.images), list(self.installed), list(self.planned),
                list(self.checks))  # fmt: skip

    def clear(self) -> None:
        for lst in (self.uniop, self.images, self.installed, self.planned, self.checks):
            lst.clear()


@pytest.fixture
def calls(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Calls:
    import gippa.setupdb as sdb
    from gippa.exitcodes import ExitCode
    from gippa.setupdb import images as img
    from gippa.setupdb import run as srun
    from gippa.setupdb import uniop_build
    from gippa.setupdb.items import DbRecord

    c = Calls()

    def fake_install(item: Any, db_dir: Path, _runner: Any, _fetcher: Any) -> DbRecord:
        c.installed.append(item.name)
        return DbRecord(name=item.name, version=item.version, url=item.url,
                        checksum=item.checksum or "sha256:" + "0" * 64, size_bytes=1,
                        retrieved="2026-11-03", path=item.path, status="ok",
                        message="downloaded")  # fmt: skip

    def fake_images(tools: list[str], *_a: Any, **_k: Any) -> ExitCode:
        c.images.append(list(tools))
        return ExitCode.SUCCESS

    def fake_uniop(profile: str, *_a: Any, accepted: bool, explicit: bool, **_k: Any) -> ExitCode:
        c.uniop.append({"profile": profile, "accepted": accepted, "explicit": explicit})
        return ExitCode.SUCCESS

    real_run_setup = sdb.run_setup

    def recording_run_setup(opts: Any, *a: Any, **k: Any) -> Any:
        c.planned.append({i.name for i in sdb.plan(opts)})
        c.checks.append(bool(opts.check))
        return real_run_setup(opts, *a, **k)

    def deny(*a: Any, **_k: Any) -> Any:
        raise FileNotFoundError(f"blocked in this test: {a!r}")

    def deny_net(*a: Any, **_k: Any) -> Any:
        raise AssertionError(f"gippa setup touched the network: {a!r}")

    monkeypatch.setattr(srun, "install", fake_install)
    monkeypatch.setattr(img, "install_images", fake_images)
    monkeypatch.setattr(uniop_build, "setup_uniop", fake_uniop)
    monkeypatch.setattr(sdb, "run_setup", recording_run_setup)
    monkeypatch.setattr(sdb.HttpsFetcher, "fetch", deny_net)
    monkeypatch.setattr(subprocess, "run", deny)
    monkeypatch.setattr(subprocess, "Popen", deny)
    monkeypatch.setattr(socket, "create_connection", deny_net)
    monkeypatch.setattr(socket.socket, "connect", deny_net)
    monkeypatch.setattr(shutil, "disk_usage", lambda _p: HUGE)
    for var in ("GIPPA_DB_DIR", "GIPPA_ROOT", "GIPPA_CONTAINER_RUNTIME", "NXF_APPTAINER_CACHEDIR",
                "APPTAINER_CACHEDIR", "SINGULARITY_CACHEDIR", "SLURM_CLUSTER_NAME"):  # fmt: skip
        monkeypatch.delenv(var, raising=False)
    monkeypatch.setenv("SLURM_JOB_ID", "123")
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg_config"))
    monkeypatch.setenv("XDG_DATA_HOME", str(tmp_path / "xdg_data"))
    _tty(monkeypatch, False)
    return c


def _run(tmp_path: Path, *args: str, answer: str | None = None, db: str = "db") -> tuple[int, str]:
    return _invoke(["--db-dir", str(tmp_path / db), "--cache-dir", str(tmp_path / "cache"),
                    "--profile", "docker", *args], answer)  # fmt: skip


def _no_unknown_option(out: str) -> None:
    assert "no such option" not in out.lower(), out


# ---------------------------------------------------------------- 1. --accept-download-plan


def test_accept_download_plan_skips_proceed_on_a_tty(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _tty(monkeypatch)
    code, out = _run(tmp_path, *SKIP_DBS, "--no-uniop", NEW_PLAN, answer="n\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert PROCEED not in out, out
    assert len(calls.images) == 1 and calls.images[0], "the planned images are pulled"


def test_without_accept_download_plan_a_tty_is_still_asked(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Control for the test above: the same command without the flag asks, and "n" stops it."""
    _tty(monkeypatch)
    code, out = _run(tmp_path, *SKIP_DBS, "--no-uniop", answer="n\n")
    assert code == 0, out
    assert PROCEED in out
    assert calls.images == []


@pytest.mark.parametrize("tty", [True, False], ids=["tty", "no-tty"])
def test_accept_download_plan_does_not_accept_the_uniop_terms(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tty: bool
) -> None:
    _tty(monkeypatch, tty)
    code, out = _run(tmp_path, *SKIP_DBS, NEW_PLAN, answer="n\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": False, "explicit": False}]


def test_accept_download_plan_prints_no_deprecation(calls: Calls, tmp_path: Path) -> None:
    code, out = _run(tmp_path, *SKIP_DBS, NEW_PLAN)
    _no_unknown_option(out)
    assert code == 0, out
    assert _deprecations(out) == [], out


# ---------------------------------------------------------------- 2. the one full setup


def test_full_setup_in_one_command_without_any_prompt(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Khoa's goal, on a TTY: databases planned and downloaded, images pulled, UniOP built with the
    terms accepted (the normal ``needed`` path, not ``--uniop``), and nothing asked."""
    _tty(monkeypatch)
    code, out = _run(tmp_path, NEW_PLAN, NEW_UNIOP, answer="n\nn\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert calls.planned and calls.planned[0] - {"uniop_source"}, "databases are planned"
    assert calls.installed, "databases are downloaded"
    assert len(calls.images) == 1 and calls.images[0], "images are pulled"
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": False}]
    assert PROCEED not in out and PROMPT not in out, out
    assert "[y/n]" not in out.lower(), "no prompt of any kind"
    assert _deprecations(out) == [], out


@pytest.mark.parametrize("tty", [True, False], ids=["tty", "no-tty"])
def test_install_uniop_from_git_does_not_imply_uniop(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tty: bool
) -> None:
    """Without ``--accept-download-plan`` the plan is still confirmed (default yes on Enter)."""
    _tty(monkeypatch, tty)
    code, out = _run(tmp_path, NEW_UNIOP, answer="\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert calls.installed, "databases are downloaded"
    assert len(calls.images) == 1, "images are pulled"
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": False}]
    assert PROMPT not in out
    assert _deprecations(out) == [], out


def test_install_uniop_from_git_is_the_old_accept_uniop_terms(calls: Calls, tmp_path: Path) -> None:
    code_new, out_new = _run(tmp_path, NEW_UNIOP, db="new")
    _no_unknown_option(out_new)
    new = calls.state()
    calls.clear()
    code_old, out_old = _run(tmp_path, OLD_TERMS, db="old")
    assert code_new == code_old == 0, (out_new, out_old)
    assert new == calls.state(), (new, calls.state())


@pytest.mark.parametrize("tty", [True, False], ids=["tty", "no-tty"])
@pytest.mark.parametrize("order", [(UNIOP, NEW_UNIOP), (NEW_UNIOP, UNIOP)], ids=["u-g", "g-u"])
def test_uniop_with_install_uniop_from_git_installs_uniop_only(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, tty: bool,
    order: tuple[str, str],
) -> None:  # fmt: skip
    """Q-K27 with the terms accepted: explicit, accepted, no DB, no image, no prompt."""
    _tty(monkeypatch, tty)
    code, out = _run(tmp_path, *order, answer="n\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": True}]
    assert calls.images == [] and calls.installed == []
    assert all(p <= {"uniop_source"} for p in calls.planned), calls.planned
    assert PROCEED not in out and PROMPT not in out, out
    assert _deprecations(out) == [], out


def test_uniop_with_install_uniop_from_git_is_the_old_pair(calls: Calls, tmp_path: Path) -> None:
    code_new, out_new = _run(tmp_path, UNIOP, NEW_UNIOP, db="new")
    _no_unknown_option(out_new)
    new = calls.state()
    calls.clear()
    code_old, out_old = _run(tmp_path, UNIOP, OLD_TERMS, db="old")
    assert code_new == code_old == 0, (out_new, out_old)
    assert new == calls.state(), (new, calls.state())


# ---------------------------------------------------------------- 3. usage errors


@pytest.mark.parametrize("other", ["--no-uniop", "--no-images"])
@pytest.mark.parametrize("first", [True, False], ids=["flag-first", "other-first"])
def test_install_uniop_from_git_conflicts(
    calls: Calls, tmp_path: Path, other: str, first: bool
) -> None:
    args = (NEW_UNIOP, other) if first else (other, NEW_UNIOP)
    code, out = _run(tmp_path, *args)
    _no_unknown_option(out)
    assert code == 2, out
    assert _flag(NEW_UNIOP).search(out) and _flag(other).search(out), out
    assert calls.uniop == [] and calls.images == [] and calls.installed == []
    assert not (tmp_path / "db" / "db_manifest.json").exists()


def test_install_uniop_from_git_no_images_is_allowed_with_uniop(
    calls: Calls, tmp_path: Path
) -> None:
    code, out = _run(tmp_path, UNIOP, NEW_UNIOP, "--no-images")
    _no_unknown_option(out)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": True}]
    assert calls.images == []


def test_install_uniop_from_git_no_uniop_is_an_error_even_with_uniop(
    calls: Calls, tmp_path: Path
) -> None:
    code, out = _run(tmp_path, UNIOP, NEW_UNIOP, "--no-uniop")
    assert code == 2, out
    assert calls.uniop == []


_HARMLESS = ["--no-go", "--no-kegg", "--no-cog", NEW_PLAN, "--bakta=skip", "--no-genomad"]


@settings(max_examples=25, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])  # fmt: skip
@given(
    st.sampled_from(["--no-uniop", "--no-images"]).flatmap(
        lambda other: st.lists(st.sampled_from(_HARMLESS), unique=True).flatmap(
            lambda extra: st.permutations([NEW_UNIOP, other, *extra])
        )
    )
)
def test_conflict_is_a_usage_error_in_any_order(
    calls: Calls, tmp_path: Path, argv: list[str]
) -> None:
    """Property: whatever else is passed and in whatever order, the pair is exit 2, the message
    names both flags and nothing is set up."""
    other = "--no-uniop" if "--no-uniop" in argv else "--no-images"
    code, out = _run(tmp_path, *argv)
    assert code == 2, (argv, out)
    assert _flag(NEW_UNIOP).search(out) and _flag(other).search(out), (argv, out)
    assert calls.uniop == [] and calls.images == [] and calls.installed == []


# ---------------------------------------------------------------- 4. deprecated flags still work


def test_yes_still_skips_proceed_and_warns_once(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _tty(monkeypatch)
    code, out = _run(tmp_path, *SKIP_DBS, "--no-uniop", OLD_YES, answer="n\n")
    assert code == 0, out
    assert PROCEED not in out
    assert len(calls.images) == 1
    dep = _deprecations(out)
    assert len(dep) == 1, out
    assert _flag(OLD_YES).search(dep[0]) and _flag(NEW_PLAN).search(dep[0]), dep


def test_yes_does_not_accept_the_uniop_terms(calls: Calls, tmp_path: Path) -> None:
    code, out = _run(tmp_path, *SKIP_DBS, OLD_YES)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": False, "explicit": False}]


def test_accept_uniop_terms_alone_under_images_still_accepts_and_warns_once(
    calls: Calls, tmp_path: Path
) -> None:
    code, out = _run(tmp_path, *SKIP_DBS, OLD_TERMS)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": False}]
    assert len(calls.images) == 1, "the other images are still pulled"
    dep = _deprecations(out)
    assert len(dep) == 1, out
    assert _flag(OLD_TERMS).search(dep[0]) and _flag(NEW_UNIOP).search(dep[0]), dep


def test_uniop_with_accept_uniop_terms_warns_once_about_the_terms_flag_only(
    calls: Calls, tmp_path: Path
) -> None:
    code, out = _run(tmp_path, UNIOP, OLD_TERMS)
    assert code == 0, out
    assert calls.uniop == [{"profile": "docker", "accepted": True, "explicit": True}]
    dep = _deprecations(out)
    assert len(dep) == 1, out
    assert _flag(OLD_TERMS).search(dep[0]) and _flag(NEW_UNIOP).search(dep[0]), dep
    assert not _flag(UNIOP).search(dep[0]), "--uniop is not deprecated"


def test_uniop_alone_is_not_deprecated_and_still_explicit_without_terms(
    calls: Calls, tmp_path: Path
) -> None:
    _code, out = _run(tmp_path, UNIOP)
    assert calls.uniop == [{"profile": "docker", "accepted": False, "explicit": True}], out
    assert calls.images == [] and calls.installed == []
    assert _deprecations(out) == [], out


def test_uniop_and_no_uniop_still_exclude_each_other(calls: Calls, tmp_path: Path) -> None:
    code, out = _run(tmp_path, UNIOP, "--no-uniop")
    assert code == 2, out
    assert calls.uniop == []


# ---------------------------------------------------------------- 5. --help


def _help() -> str:
    from gippa.cli import cli

    res = CliRunner().invoke(cli, ["setup", "--help"], terminal_width=200)
    assert res.exit_code == 0, res.output
    return res.output


@pytest.mark.parametrize("flag", [NEW_PLAN, NEW_UNIOP, UNIOP, "--no-uniop"])
def test_help_lists_the_flags(flag: str) -> None:
    text = _help()
    assert _flag(flag).search(text), f"{flag} not in --help:\n{text}"


@pytest.mark.parametrize("flag", DEPRECATED)
def test_help_hides_the_deprecated_flags(flag: str) -> None:
    text = _help()
    assert not _flag(flag).search(text), f"{flag} still in --help:\n{text}"


# ---------------------------------------------------------------- 6. --check and --minimal


@pytest.mark.parametrize("extra", [(), (UNIOP,)], ids=["plain", "uniop"])
def test_check_with_install_uniop_from_git_is_check_with_the_old_flag(
    calls: Calls, tmp_path: Path, extra: tuple[str, ...]
) -> None:
    """``--check`` only verifies: no build and no pull, exactly as with ``--accept-uniop-terms``."""
    code_new, out_new = _run(tmp_path, "--check", *extra, NEW_UNIOP, db="new")
    _no_unknown_option(out_new)
    new = calls.state()
    calls.clear()
    code_old, out_old = _run(tmp_path, "--check", *extra, OLD_TERMS, db="old")
    old = calls.state()
    assert code_new == code_old, (out_new, out_old)
    assert new == old, (new, old)
    assert new[0] == [] and new[1] == [] and new[2] == [], "--check builds and pulls nothing"


def test_check_with_accept_download_plan_downloads_nothing(calls: Calls, tmp_path: Path) -> None:
    _code, out = _run(tmp_path, "--check", NEW_PLAN)
    _no_unknown_option(out)
    assert calls.installed == [] and calls.images == [] and calls.uniop == []
    assert all(calls.checks), calls.checks


def test_minimal_with_accept_download_plan_is_minimal_with_yes(
    calls: Calls, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--minimal`` keeps its DBs and images, and (images on, UniOP is ``core``) still runs the
    UniOP step without consent; ``--accept-download-plan`` only removes the question."""
    _tty(monkeypatch)
    code_new, out_new = _run(tmp_path, "--minimal", NEW_PLAN, answer="n\n", db="new")
    _no_unknown_option(out_new)
    new = calls.state()
    calls.clear()
    code_old, out_old = _run(tmp_path, "--minimal", OLD_YES, answer="n\n", db="old")
    assert code_new == code_old == 0, (out_new, out_old)
    assert new == calls.state(), (new, calls.state())
    assert PROCEED not in out_new
    assert new[0] == [{"profile": "docker", "accepted": False, "explicit": False}]
    assert "genomad" not in new[2] and "antismash" not in new[2], new[2]


def test_minimal_with_install_uniop_from_git_is_minimal_with_accept_uniop_terms(
    calls: Calls, tmp_path: Path
) -> None:
    code_new, out_new = _run(tmp_path, "--minimal", NEW_UNIOP, db="new")
    _no_unknown_option(out_new)
    new = calls.state()
    calls.clear()
    code_old, out_old = _run(tmp_path, "--minimal", OLD_TERMS, db="old")
    assert code_new == code_old == 0, (out_new, out_old)
    assert new == calls.state(), (new, calls.state())
    assert new[0] == [{"profile": "docker", "accepted": True, "explicit": False}]
    assert "genomad" not in new[2] and "antismash" not in new[2], new[2]


# ---------------------------------------------------------------- 7. the preflight hint


def test_preflight_hint_adds_uniop_only_with_the_new_flag() -> None:
    from gippa.preflight.environment import uniop_setup_command

    want = f"{PKG} setup --uniop --install-uniop-from-git --profile apptainer --db-dir /data/db"
    assert uniop_setup_command("apptainer", "/data/db") == want


def test_preflight_hint_placeholders() -> None:
    from gippa.preflight.environment import uniop_setup_command

    want = f"{PKG} setup --uniop --install-uniop-from-git --profile <profile> --db-dir <db_dir>"
    assert uniop_setup_command() == want


# ======================================================== behaviour (real UniOP step, fakes)


@pytest.fixture
def fake_host(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Any:
    isolate(tmp_path, monkeypatch)
    install_fake_runtimes(tmp_path, monkeypatch)
    _tty(monkeypatch, False)
    old = os.umask(0o022)
    yield
    os.umask(old)


@pytest.fixture
def net(fake_host: None, monkeypatch: pytest.MonkeyPatch) -> FakeNet:
    return install_fake_net(monkeypatch)


@pytest.fixture
def archive(net: FakeNet, monkeypatch: pytest.MonkeyPatch) -> bytes:
    data = fixture_archive()
    pin_uniop_source(monkeypatch, sha256_hex(data))
    net.serve_github(data)
    return data


@pytest.fixture
def marked(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Each fetch prints ``MARK`` on stdout first (then fetches as before)."""
    import gippa.setupdb as sdb

    real = sdb.HttpsFetcher.fetch
    urls: list[str] = []

    def fetch(self: Any, url: str, dest: Path) -> None:
        print(MARK, flush=True)
        urls.append(url)
        real(self, url, dest)

    monkeypatch.setattr(sdb.HttpsFetcher, "fetch", fetch)
    return urls


def _full(db: Path, *extra: str) -> list[str]:
    """A full setup (images on, UniOP needed) with the DB downloads opted out (no fake for them)."""
    return ["--db-dir", str(db), *SKIP_DBS, "--profile", "docker", *extra]


def _uniop_only(db: Path, *extra: str) -> list[str]:
    return _full(db, UNIOP, NEW_UNIOP, *extra)


@pytest.mark.parametrize("tty", [True, False], ids=["tty", "no-tty"])
def test_full_setup_pulls_images_and_builds_uniop_without_a_prompt(
    archive: bytes, marked: list[str], tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    tty: bool,
) -> None:  # fmt: skip
    _tty(monkeypatch, tty)
    db = tmp_path / "db"
    code, out = _invoke(["--cache-dir", str(tmp_path / "cache"),
                         *_full(db, NEW_PLAN, NEW_UNIOP)], answer="n\nn\n")  # fmt: skip
    _no_unknown_option(out)
    assert code == 0, out
    assert PROCEED not in out and PROMPT not in out, out
    items = manifest_items(db)
    assert items["image:bakta"]["status"] == "ok", "the images are pulled"
    assert items["image:uniop"]["status"] == "ok", items["image:uniop"]
    assert items["image:uniop"]["local_build"]["terms"] == "flag"
    assert len(docker_calls(tmp_path, "build")) == 1
    assert out.count(NOTICE) == 1, out
    assert marked, "the archive was fetched"
    assert out.index(NOTICE) < out.index(MARK), "the notice comes before anything is fetched"


def test_uniop_install_uniop_from_git_builds_uniop_only(
    archive: bytes, net: FakeNet, marked: list[str], tmp_path: Path
) -> None:
    db = tmp_path / "db"
    code, out = _invoke(["--cache-dir", str(tmp_path / "cache"), *_uniop_only(db)])
    _no_unknown_option(out)
    assert code == 0, out
    assert len(docker_calls(tmp_path, "build")) == 1
    assert docker_calls(tmp_path, "pull") == [], "no image pulled under --uniop"
    assert {h for h, _ in net.contacted} <= {"github.com", "codeload.github.com"}, net.contacted
    items = manifest_items(db)
    assert set(items) <= {"image:uniop", "uniop_source"}, set(items)
    assert items["image:uniop"]["status"] == "ok"
    assert items["image:uniop"]["local_build"]["terms"] == "flag"
    assert out.count(NOTICE) == 1, out
    assert out.index(NOTICE) < out.index(MARK)
    assert PROMPT not in out and TERMS_NOT_ACCEPTED not in out


def test_uniop_install_uniop_from_git_never_prompts_on_a_tty(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_uniop_only(db), answer="n\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert PROMPT not in out and PROCEED not in out, out
    assert manifest_items(db)["image:uniop"]["local_build"]["terms"] == "flag"


@pytest.mark.parametrize("uniop_only", [True, False], ids=["uniop", "full"])
def test_install_uniop_from_git_checks_the_sha256(
    net: FakeNet, archive: bytes, tmp_path: Path, uniop_only: bool
) -> None:
    bad = bytearray(archive)
    bad[len(bad) // 2] ^= 0x01
    net.serve_github(bytes(bad))
    db = tmp_path / "db"
    code, out = _invoke(_uniop_only(db) if uniop_only else _full(db, NEW_UNIOP))
    _no_unknown_option(out)
    assert code == 2, out
    assert "uniop: archive sha256" in " ".join(out.split()), out
    assert docker_calls(tmp_path, "build") == [], "nothing is built from a mismatching archive"


@pytest.mark.parametrize("other", ["--no-uniop", "--no-images"])
def test_conflicts_fetch_nothing(archive: bytes, net: FakeNet, tmp_path: Path, other: str) -> None:
    code, out = _invoke(_full(tmp_path / "db", NEW_UNIOP, other))
    assert code == 2, out
    assert _flag(NEW_UNIOP).search(out) and _flag(other).search(out), out
    assert net.contacted == [] and docker_calls(tmp_path) == []
    assert NOTICE not in out


def test_uniop_without_terms_tells_to_use_the_new_flag(
    archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    """``--uniop`` alone, no TTY: exit 2 as before, the message names ``--install-uniop-from-git``
    (and no longer sends the user to the deprecated flag)."""
    db = tmp_path / "db"
    code, out = _invoke(_full(db, UNIOP))
    assert code == 2, out
    assert net.contacted == [] and docker_calls(tmp_path, "build") == []
    lines = [ln for ln in out.splitlines() if TERMS_NOT_ACCEPTED in ln]
    assert lines, out
    assert all(_flag(NEW_UNIOP).search(ln) for ln in lines), lines
    assert not any(_flag(OLD_TERMS).search(ln) for ln in lines), lines
    msg = manifest_items(db)["image:uniop"]["message"]
    assert TERMS_NOT_ACCEPTED in msg and _flag(NEW_UNIOP).search(msg), msg


def test_images_without_terms_warning_names_the_new_flag(
    archive: bytes, net: FakeNet, tmp_path: Path
) -> None:
    """Under ``--images`` the same message is a warning (exit 0), with the new flag."""
    code, out = _invoke(_full(tmp_path / "db", "--images"))
    assert code == 0, out
    lines = [ln for ln in out.splitlines() if TERMS_NOT_ACCEPTED in ln]
    assert lines and all(_flag(NEW_UNIOP).search(ln) for ln in lines), out
    assert not any(_flag(OLD_TERMS).search(ln) for ln in lines), lines


def test_uniop_without_terms_on_a_tty_still_prompts(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_full(db, UNIOP), answer="y\n")
    assert code == 0, out
    assert PROMPT in out
    assert manifest_items(db)["image:uniop"]["local_build"]["terms"] == "prompt"


def test_accept_download_plan_on_a_tty_still_asks_the_uniop_terms(
    archive: bytes, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """``--accept-download-plan`` on a TTY: no ``Proceed?``, but the UniOP prompt still comes
    (its answer is the only line on stdin) and a "n" builds nothing."""
    _tty(monkeypatch)
    db = tmp_path / "db"
    code, out = _invoke(_full(db, "--images", NEW_PLAN), answer="n\n")
    _no_unknown_option(out)
    assert code == 0, out
    assert PROCEED not in out
    assert PROMPT in out, out
    assert docker_calls(tmp_path, "build") == []
    items = manifest_items(db)
    assert items["image:uniop"]["status"] == "absent"
    assert items["image:bakta"]["status"] == "ok"


# ---------------------------------------------------------------- the help explains the flags
# Khoa 2026-11-03: "explain the flags for user as well". Key phrases only, not the exact text.


def _option_help(flag: str) -> str:
    """The ``--help`` text of ``flag`` (lower case, whitespace collapsed), as the user reads it."""
    from gippa.cli import cli

    cmd = cli.commands["setup"]  # type: ignore[attr-defined]
    found = [p for p in cmd.params if flag in getattr(p, "opts", [])]
    assert found, f"gippa setup has no option {flag}"
    opt = found[0]
    assert not getattr(opt, "hidden", False), f"{flag} is hidden from --help"
    text = " ".join((opt.help or "").split()).lower()
    assert text, f"{flag} has no help text"
    assert " ".join(text.split()[:4]) in " ".join(_help().split()).lower(), "shown in --help"
    return text


def test_help_accept_download_plan_names_the_download_plan() -> None:
    text = _option_help(NEW_PLAN)
    assert "plan" in text and "download" in text, text


def test_help_accept_download_plan_says_it_does_not_accept_licence_terms() -> None:
    text = _option_help(NEW_PLAN)
    assert re.search(r"\b(not|n't|never|nor)\b.*\b(licen[cs]es?|terms)\b", text), text


def test_help_install_uniop_from_git_names_the_authors_github_repository() -> None:
    text = _option_help(NEW_UNIOP)
    assert "github" in text, text
    assert "author" in text, text


def test_help_install_uniop_from_git_names_the_pinned_commit_or_sha256() -> None:
    text = _option_help(NEW_UNIOP)
    assert re.search(r"\bpinned\b|\bcommit\b|sha-?256|3ed14dfc", text), text


def test_help_install_uniop_from_git_says_passing_it_accepts_the_terms() -> None:
    text = _option_help(NEW_UNIOP)
    assert re.search(r"accept", text) and "terms" in text, text
    assert re.search(r"\bno licen[cs]e\b|licen[cs]e file", text), text


def test_help_install_uniop_from_git_says_it_builds_during_the_normal_setup() -> None:
    text = _option_help(NEW_UNIOP)
    assert re.search(
        r"\b(normal|full|regular|whole|rest of( the)?|with the|along with( the)?|during( the)?)"
        r"\b[^.;]*\b(setup|databases?|images)\b",
        text,
    ), text


def test_help_install_uniop_from_git_says_add_uniop_for_uniop_only() -> None:
    text = _option_help(NEW_UNIOP)
    assert _flag(UNIOP).search(text), f"--uniop not named: {text}"
    assert re.search(r"--uniop\b[^.;]*\bonly\b|\bonly\b[^.;]*--uniop\b", text), text


def test_help_uniop_still_means_uniop_only() -> None:
    text = _option_help(UNIOP)
    assert "uniop" in text and "deprecat" not in text, text
