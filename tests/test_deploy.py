"""The deploy kit, checked from the repo rather than from the droplet.

Two things worth testing without a server: that the shell scripts parse, and
that the values duplicated across Python, nginx, gunicorn and systemd actually
agree. Every one of those pairs is annotated "MUST stay in sync" in both files,
which is a convention -- this is the part that enforces it.
"""

from __future__ import annotations

import re
import subprocess

import pytest

from conftest import REPO_ROOT

DEPLOY = REPO_ROOT / "deploy"


def _read(name: str) -> str:
    return (DEPLOY / name).read_text()


# ---------------------------------------------------------------------------
#  The kit is complete and runnable
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", [
    "deploy.sh", "provision.sh", "gunicorn.conf.py",
    "flexappeal-web.service", "flexappeal-scratch-clean.service",
    "flexappeal-scratch-clean.timer", "nginx-flexappeal.conf",
    "nginx-flexappeal-limits.conf", ".env.example",
])
def test_the_kit_is_complete(name):
    assert (DEPLOY / name).is_file(), f"deploy/{name} is missing"


@pytest.mark.parametrize("script", ["deploy.sh", "provision.sh"])
def test_shell_scripts_parse(script):
    """`bash -n` catches the quoting mistakes that only bite at deploy time."""
    result = subprocess.run(["bash", "-n", str(DEPLOY / script)],
                            capture_output=True, text=True)
    assert result.returncode == 0, f"{script} has a syntax error:\n{result.stderr}"


@pytest.mark.parametrize("script", ["deploy.sh", "provision.sh"])
def test_shell_scripts_are_strict(script):
    assert "set -euo pipefail" in _read(script), \
        f"{script} would carry on after a failed step"


# ---------------------------------------------------------------------------
#  Cross-file invariants
# ---------------------------------------------------------------------------


def test_upload_cap_matches_nginx():
    """nginx rejects an oversized body before Python sees it; Flask catches the
    rest. A mismatch means one of them produces an unhelpful error."""
    from flexappeal.webapp import MAX_CONTENT_LENGTH

    nginx_mb = int(re.search(r"client_max_body_size (\d+)m",
                             _read("nginx-flexappeal.conf")).group(1))
    assert MAX_CONTENT_LENGTH // (1024 * 1024) == nginx_mb


def test_gunicorn_and_nginx_timeouts_match():
    gunicorn = int(re.search(r"^timeout = (\d+)", _read("gunicorn.conf.py"),
                             re.M).group(1))
    nginx = int(re.search(r"proxy_read_timeout (\d+)s",
                          _read("nginx-flexappeal.conf")).group(1))
    assert gunicorn == nginx


def test_the_port_is_consistent_everywhere():
    """8004 on a droplet that also runs four other apps."""
    gunicorn = re.search(r"127\.0\.0\.1:(\d+)", _read("gunicorn.conf.py")).group(1)
    env = re.search(r"BIND_ADDR=127\.0\.0\.1:(\d+)", _read(".env.example")).group(1)
    assert gunicorn == env == "8004"


def test_the_port_does_not_collide_with_a_sibling_app():
    taken = {"8000": "AlphaFraud", "8001": "chem_sage",
             "8002": "chatPDB", "8003": "BoltzMaker"}
    port = re.search(r"127\.0\.0\.1:(\d+)", _read("gunicorn.conf.py")).group(1)
    assert port not in taken, f"port {port} already belongs to {taken.get(port)}"


def test_the_analysis_budget_fits_inside_the_cgroup_limit():
    """MemoryMax is the thing that actually holds when the budget's estimate is
    optimistic; the budget has to be chosen to fit inside it."""
    from flexappeal.analysis import BUDGET_ATOM_FRAMES

    memory_mb = int(re.search(r"MemoryMax=(\d+)M",
                              _read("flexappeal-web.service")).group(1))
    # 12 bytes per atom-frame for coordinates, and MDTraj's atom_slice and
    # superpose each take a copy -- so allow for several times the raw array.
    coordinates_mb = BUDGET_ATOM_FRAMES * 12 / 1e6
    assert coordinates_mb * 4 < memory_mb, (
        f"a maximum-size job is {coordinates_mb:.0f} MB of coordinates; with "
        f"MDTraj's copies that does not fit in {memory_mb} MB"
    )


def test_the_service_user_matches_across_the_kit():
    for name in ("flexappeal-web.service", "flexappeal-scratch-clean.service"):
        assert "User=flexappeal" in _read(name)
    assert "chown flexappeal:flexappeal" in _read("deploy.sh")
    assert 'useradd --system' in _read("provision.sh")


def test_static_path_matches_where_the_files_actually_are():
    alias = re.search(r"alias ([^;]+);", _read("nginx-flexappeal.conf")).group(1)
    assert alias == "/opt/flexappeal/flexappeal/static/"
    # And that path must exist relative to the repo, or nginx serves nothing.
    assert (REPO_ROOT / "flexappeal" / "static").is_dir()
    assert (REPO_ROOT / "flexappeal" / "static" / "vendor" / "plotly.min.js").is_file()


# ---------------------------------------------------------------------------
#  Lessons the sibling projects paid for
# ---------------------------------------------------------------------------


def test_certbot_runs_unconditionally():
    """Skipping certbot when a certificate exists leaves the re-templated vhost
    with no TLS block, and nginx then answers HTTPS from another app's cert.
    boltzmaker.mdeller.com once served AlphaFraud's."""
    text = _read("provision.sh")
    assert "certbot --nginx" in text
    assert "for attempt in 1 2 3" in text, "certbot's own timer can hold the lock"


def test_http2_is_patched_in():
    """certbot does not enable HTTP/2 on nginx 1.24."""
    assert "http2" in _read("provision.sh")


def test_deploy_chowns_on_every_run_not_only_provisioning():
    """rsync from a Mac preserves 0600 dellboy:staff, and the service user then
    fails at runtime rather than at deploy time."""
    text = _read("deploy.sh")
    assert "chown flexappeal:flexappeal" in text
    assert "-prune" in text, "the venv and live scratch must be pruned from the chown"


def test_deploy_excludes_local_only_directories():
    text = _read("deploy.sh")
    for excluded in (".pixi/", ".venv/", "web_scratch/", ".env"):
        assert f"--exclude '{excluded}'" in text, f"{excluded} would be synced"


def test_deploy_does_not_exclude_the_vendored_javascript():
    """nginx serves it from disk, so it genuinely has to be on the server."""
    assert "--exclude 'flexappeal/static" not in _read("deploy.sh")


def test_provision_warns_when_dns_is_missing():
    """certbot proves domain control over HTTP; without the A record it fails."""
    assert "does not resolve" in _read("provision.sh")


def test_provision_chowns_before_building_the_venv():
    """rsync runs as root, so /opt/flexappeal arrives root-owned.

    Building the virtual environment as the service user then fails with a bare
    "Permission denied: '/opt/flexappeal/.venv'", which says nothing about
    ownership being the cause. Caught on the first real provisioning run.
    """
    text = _read("provision.sh")
    first_chown = text.index("chown -R")
    venv_build = text.index("python3 -m venv")
    assert first_chown < venv_build, \
        "the venv is built before ownership is fixed; it will fail as the service user"


def test_provision_chowns_again_at_the_end():
    """The steps between create files as root, so once is not enough."""
    assert _read("provision.sh").count("chown -R") >= 2


def test_http2_patch_matches_both_listen_directives():
    """The regex must patch the IPv4 line as well as the IPv6 one.

    In nginx 1.24 `ssl` and `http2` are protocol options on the listening
    socket, not per-server settings. One vhost declaring 0.0.0.0:443
    differently from its neighbours makes nginx warn

        [warn] protocol options redefined for 0.0.0.0:443

    and honour whichever server block was parsed first. The original pattern
    put the space inside the optional group, so it required "listen443" when
    that group was absent and silently patched only IPv6.
    """
    import re

    pattern = re.search(r're\.sub\(r"([^"]+)"', _read("provision.sh")).group(1)
    compiled = re.compile(pattern.encode().decode("unicode_escape"))

    for line in ("    listen [::]:443 ssl; # managed by Certbot",
                 "    listen 443 ssl; # managed by Certbot"):
        assert compiled.search(line), f"the http2 patch would miss: {line.strip()}"


def test_http2_patch_is_idempotent():
    """provision.sh is meant to be re-run; a second pass must not double-apply."""
    import re

    pattern = re.search(r're\.sub\(r"([^"]+)"', _read("provision.sh")).group(1)
    compiled = re.compile(pattern.encode().decode("unicode_escape"))
    assert not compiled.search("    listen 443 ssl http2; # managed by Certbot")


def test_the_droplet_can_import_the_cli_it_spawns():
    """The /reanalyse route Popen's FlexAppeal.py, so its imports are runtime deps.

    This is easy to get wrong in exactly one direction: the CLI's dependencies
    look like developer tooling, so they get left out of requirements.txt, and
    the failure surfaces only as a re-analysis job that never completes -- the
    subprocess dies on ImportError with stdout and stderr both at DEVNULL.
    """
    import ast
    import pathlib

    root = pathlib.Path(__file__).resolve().parent.parent
    requirements = (root / "requirements.txt").read_text()

    # What FlexAppeal.py pulls in transitively via flexappeal/console.py.
    console = ast.parse((root / "flexappeal" / "console.py").read_text())
    third_party = {
        node.module.split(".")[0]
        for node in ast.walk(console)
        if isinstance(node, ast.ImportFrom) and node.module and node.level == 0
    } | {
        alias.name.split(".")[0]
        for node in ast.walk(console)
        if isinstance(node, ast.Import) for alias in node.names
    }
    stdlib = {"shutil", "sys", "__future__", "os", "json", "pathlib"}

    for module in third_party - stdlib:
        assert module in requirements, (
            f"flexappeal/console.py imports {module!r}, which the CLI loads at "
            f"module scope, but requirements.txt does not install it")


def test_the_example_artefacts_are_not_gitignored():
    """The Example tab needs files that .gitignore's broad rules would exclude.

    `*.fxa` and `*.command` are ignored globally -- they are large and
    regenerable everywhere except here, where the committed copies are the
    entire point of the page. Negations bring them back, and a negation that
    does not match is invisible on a Mac: the filesystem is case-insensitive, so
    `!examples/...` appears to work locally while matching nothing on the Linux
    droplet. The result would be an Example tab that renders in development and
    shows "not built yet" in production.
    """
    import pathlib
    import subprocess

    root = pathlib.Path(__file__).resolve().parent.parent
    if not (root / ".git").exists():
        pytest.skip("not a git checkout")

    required = [
        "Examples/lysozyme_10ns/config.json",
        "Examples/lysozyme_10ns/flexappeal_lysozyme_10ns.command",
        "Examples/lysozyme_10ns/output/lysozyme_10ns.fxa",
    ]
    ignored = subprocess.run(
        ["git", "check-ignore", "--no-index", *required],
        cwd=root, capture_output=True, text=True).stdout.split()
    assert not ignored, (
        f".gitignore excludes files the Example tab requires: {ignored}. "
        f"Check the negation rules match the directory's real case.")


def test_the_example_directory_case_matches_what_the_code_expects():
    """Linux is case-sensitive; this Mac is not, so only a check catches it."""
    import pathlib

    from flexappeal.webapp import EXAMPLE_DIR

    root = pathlib.Path(__file__).resolve().parent.parent
    on_disk = [p.name for p in root.iterdir() if p.name.lower() == "examples"]
    assert on_disk == ["Examples"], f"expected Examples/, found {on_disk}"
    assert EXAMPLE_DIR.parent.name == "Examples"
