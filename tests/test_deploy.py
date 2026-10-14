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


def test_static_requests_are_still_logged():
    """`access_log off;` in the /static/ block pins this app's hit count at zero.

    mdeller.com's launcher counts a visit for each app by watching its access log
    for one request that only a rendering browser makes -- for FlexAppeal that is
    /static/app.js, because a scanner fetches the HTML and stops. Turning logging
    off for that location hides those requests, and the count then reads zero
    forever with nothing anywhere to say why. This block did have it, and a
    provision run silently undid the fix once already.
    """
    directives = [line for line in _read("nginx-flexappeal.conf").splitlines()
                  if not line.lstrip().startswith("#")]
    assert not any("access_log off" in line for line in directives), \
        "static requests must stay in the log; mdeller.com counts hits from them"


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


def test_the_droplet_installs_what_the_run_pages_measure_with():
    """A fallback that is never exercised locally is a fallback nobody sees.

    `measure.detect_equilibration` prefers pymbar and falls back to this
    repository's own implementation of the same rule, saying which answered in
    the method string. That is honest and it is not enough on its own: with
    pymbar in the development environment and not on the droplet, the site and
    a local run printed different numbers for the same file, 3.6 ns against
    3.2, and only the method string said why.
    """
    import pathlib

    root = pathlib.Path(__file__).resolve().parent.parent
    requirements = (root / "requirements.txt").read_text()
    assert "pymbar" in requirements, (
        "the run pages detect equilibration at render time, so the droplet "
        "needs what they detect it with")


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

    # Derived from the constant the app actually reads, so changing which run
    # the Example tab shows cannot leave this test guarding the old one.
    from flexappeal.webapp import EXAMPLE_NAME

    required = [
        f"examples/{EXAMPLE_NAME}/config.json",
        f"examples/{EXAMPLE_NAME}/flexappeal_{EXAMPLE_NAME}.command",
        f"examples/{EXAMPLE_NAME}/output/{EXAMPLE_NAME}.fxa",
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
    assert on_disk == ["examples"], f"expected examples/, found {on_disk}"
    assert EXAMPLE_DIR.parent.name == "examples"


def test_the_deploy_reports_its_payload_rather_than_claiming_it():
    """The comment said "under 5 MB" and it passed 5 MB during the build.

    A results file carries a trajectory and those grow with the run, so any
    figure written into a comment is a figure that will be wrong later. There
    is plenty of room, section 6.3 capping a per-system web bundle at 150 MB,
    which is exactly why nobody would notice the comment going stale.
    """
    script = _read("deploy.sh")
    assert "du -ch web_runs" in script
    assert "to send" in script
    # The comment still quotes the old wording, which is the point of it.
    # What must be gone is the claim.
    assert "are under 5 MB" not in script, "the claim should be gone"


def test_the_deploy_still_sends_what_the_pages_read():
    """What the app serves must not join the exclude list.

    `web_runs/` is gitignored, so a tidier `--filter=':- .gitignore'` would
    ship a droplet whose Runs tab holds the committed example and nothing else.
    The reference-state runs are four results files in there, and D035 made
    them the only thing the pages read: the explorer database is written and
    kept, and no longer served.
    """
    script = _read("deploy.sh")
    # Code, not comments: the script explains at length why the gitignore
    # filter would be wrong, so the string is in the file on purpose.
    code = [line for line in script.splitlines()
            if not line.lstrip().startswith("#")]
    joined = " ".join(code)
    for served in ("web_runs/", "examples/"):
        assert f"--exclude '{served}'" not in joined, (
            f"{served} must still be sent")
    assert "--filter=':- .gitignore'" not in joined, (
        "it would ship a droplet with an empty run library")
    for excluded in ("data/md/", "data/cache/", "data/refs/"):
        assert f"--exclude '{excluded}'" in joined, f"{excluded} should not be sent"


def test_the_deploy_cannot_delete_a_run_somebody_kept():
    """`web_runs/` has two writers and `--delete` only knew about one.

    The four reference-state runs are packed here and travel with the code;
    anything a visitor keeps from the Analysis tab is written there by the
    server. Without a protect filter the next deploy would remove the second
    kind, and both sides were empty until the day the four arrived, so nothing
    would have noticed.
    """
    script = _read("deploy.sh")
    code = " ".join(line for line in script.splitlines()
                    if not line.lstrip().startswith("#"))
    assert "--delete" in code, "the sync is still meant to prune the tree"
    assert "--filter='protect web_runs/**'" in code


def test_the_rate_limit_counts_writes_rather_than_paths():
    """`/analysis/` serves the viewer as well as accepting an upload.

    The zone was keyed on the address alone and applied to the location, so a
    reader opening three run pages in a minute spent the whole allowance on the
    GETs that fetch a structure and a trajectory, and the viewer came back 503.
    It happened on the live site the first time the four reference-state runs
    were deployed into the Runs library, where there are four pages to click
    between and each needs two or three of those requests.

    nginx does not account a request whose key is empty, so mapping every read
    to an empty key limits the uploads and leaves the viewer alone.
    """
    import re

    limits = _read("nginx-flexappeal-limits.conf")
    assert "map $request_method $flexappeal_write" in limits
    assert "limit_req_zone $flexappeal_write zone=flexappeal_writes" in limits
    assert "limit_req_zone $binary_remote_addr" not in limits, (
        "keyed on the address alone, every read counts against the upload budget")
    for method in ("POST", "PUT"):
        # Whitespace-insensitive: the file aligns the map's columns, and a test
        # that counts spaces breaks the next time it is tidied.
        assert re.search(rf"\n\s*{method}\s+\$binary_remote_addr;", limits), method
    assert 'default ""' in limits, "a read must map to an empty key"

    # And the locations still reference the zone, which is the other half: a
    # zone nothing uses limits nothing.
    site = _read("nginx-flexappeal.conf")
    assert site.count("limit_req zone=flexappeal_writes") == 2
    assert "flexappeal_upload" not in site, (
        "the zone was renamed because nginx refuses a key change on reload")


def test_the_workers_preload_the_application():
    """One copy of the modules shared across the workers rather than one each.

    On a box with 3.8 GB between five applications, preloading puts them in
    pages the workers share after the fork. It was added while chasing a 12.7
    second first request and did not fix it, which is 4bv: that was the
    binding-site measurement and is cached on disk now. Kept for the memory,
    which is the reason that survived the measurement.
    """
    config = _read("gunicorn.conf.py")
    assert "preload_app = True" in config
