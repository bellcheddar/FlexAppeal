"""Web application tests, driven through Flask's test client.

Covers the whole Prepare flow: acquire a structure, configure against what was
found in it, submit. Also covers the ways a request can be hostile or malformed,
because those paths are the ones nobody exercises by hand.

No network. The one test that would fetch from the RCSB is marked and skipped by
default -- everything else uses the checked-in fixtures.
"""

from __future__ import annotations

import io
import json
import pathlib

import pytest

from flexappeal import bundle, options as opts
from flexappeal.webapp import create_app

from conftest import FIXTURES


@pytest.fixture
def app(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXAPPEAL_SCRATCH_ROOT", str(tmp_path / "scratch"))
    application = create_app()
    application.config.update(TESTING=True)
    return application


@pytest.fixture
def client(app):
    return app.test_client()


def _upload(client, fixture="1aki.pdb"):
    """Run step one and return (response, session token)."""
    data = {
        "input_source": "upload",
        "input_file": (io.BytesIO((FIXTURES / fixture).read_bytes()), fixture),
    }
    response = client.post("/prepare/structure", data=data,
                           content_type="multipart/form-data")
    html = response.get_data(as_text=True)
    token = ""
    marker = 'name="_token" value="'
    if marker in html:
        token = html.split(marker, 1)[1].split('"', 1)[0]
    return response, token, html


# ---------------------------------------------------------------------------
#  Basics
# ---------------------------------------------------------------------------


def test_healthz(client):
    response = client.get("/healthz")
    assert response.status_code == 200
    assert response.get_json()["status"] == "ok"


def test_landing_page_renders(client):
    html = client.get("/").get_data(as_text=True)
    assert html.count("md-source-tab") >= 4  # one per input source
    assert "Load structure" in html


def test_unknown_page_is_a_friendly_404(client):
    response = client.get("/no-such-page")
    assert response.status_code == 404
    assert "Traceback" not in response.get_data(as_text=True)


# ---------------------------------------------------------------------------
#  Step one: acquiring a structure
# ---------------------------------------------------------------------------


def test_upload_renders_every_option(client):
    response, token, html = _upload(client)
    assert response.status_code == 200
    assert token, "the session token was not rendered into the form"
    # Every option in the registry must reach the form, or the registry has
    # stopped being the single source of truth for the UI.
    assert html.count("data-option=") == len(opts.OPTIONS)
    assert html.count('class="md-group"') == len(opts.GROUPS)


def test_upload_reports_the_structure(client):
    _, _, html = _upload(client)
    assert "129 res" in html                      # lysozyme's chain length
    assert html.count("md-chip-ss") == 4          # its four disulfides
    assert "LYSOZYME" in html.upper()


def test_upload_seeds_the_job_name_from_the_filename(client):
    _, _, html = _upload(client)
    assert 'name="job_name" value="1aki"' in html


def test_upload_preselects_every_chain(client):
    _, _, html = _upload(client)
    assert 'name="chains" value="A"\n                   checked' in html \
        or ('value="A"' in html and "checked" in html)


def test_upload_of_a_non_structure_is_a_friendly_error(client):
    data = {
        "input_source": "upload",
        "input_file": (io.BytesIO(b"just some notes\n" * 100), "notes.txt"),
    }
    response = client.post("/prepare/structure", data=data,
                           content_type="multipart/form-data")
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "md-issue-error" in html
    assert "Traceback" not in html


def test_upload_with_no_file_is_a_friendly_error(client):
    response = client.post("/prepare/structure", data={"input_source": "upload"},
                           content_type="multipart/form-data")
    assert "choose a PDB or mmCIF file" in response.get_data(as_text=True)


def test_bad_pdb_id_is_rejected_before_any_network_call(client):
    response = client.post("/prepare/structure",
                           data={"input_source": "rcsb", "pdb_id": "ZZZZ"},
                           content_type="multipart/form-data")
    html = response.get_data(as_text=True)
    assert "four characters beginning with a digit" in html


def test_haemoglobin_offers_its_haem(client):
    _, _, html = _upload(client, "4hhb.pdb")
    assert 'value="HEM"' in html
    assert "haem b" in html


# ---------------------------------------------------------------------------
#  Session tokens
# ---------------------------------------------------------------------------


def test_build_without_a_token_is_rejected(client):
    response = client.post("/prepare/build", data={"_token": ""})
    assert response.status_code == 400


def test_build_with_an_unknown_token_is_rejected(client):
    response = client.post("/prepare/build", data={"_token": "not-a-real-session"})
    assert response.status_code == 400
    assert "expired" in response.get_data(as_text=True)


@pytest.mark.parametrize("hostile", [
    "../../../etc",
    "..%2f..%2fetc",
    "/etc/passwd",
    "a" * 200,
])
def test_token_cannot_escape_the_scratch_directory(client, hostile):
    """The token names a path and comes from the browser: the classic traversal shape."""
    response = client.post("/prepare/build", data={"_token": hostile})
    assert response.status_code == 400
    assert "Traceback" not in response.get_data(as_text=True)


# ---------------------------------------------------------------------------
#  Step two: building
# ---------------------------------------------------------------------------


def _form_from_defaults(token, **overrides):
    """Build a form payload the way a browser would submit one."""
    cfg = opts.defaults() | overrides
    data = {"_token": token, "input_source": "upload"}
    for opt in opts.OPTIONS:
        value = cfg[opt.id]
        if opt.widget == "checkbox":
            if value:
                data[opt.id] = "1"          # unticked boxes are simply absent
        elif opt.widget == "multiselect":
            data[opt.id] = list(value)
        elif opt.widget == "file":
            continue
        else:
            data[opt.id] = str(value)
    return data


def test_valid_build_returns_a_runnable_bundle(client):
    _, token, _ = _upload(client)
    response = client.post("/prepare/build",
                           data=_form_from_defaults(token, job_name="lysozyme_test"))
    assert response.status_code == 200
    assert "flexappeal_lysozyme_test.command" in response.headers["Content-Disposition"]
    # Not application/x-sh: browsers and mail gateways get helpful about shell
    # scripts, and a base64 payload has to arrive byte for byte.
    assert response.mimetype == "application/octet-stream"

    content = response.get_data()
    assert content.startswith(b"#!/usr/bin/env bash")

    files = bundle.unpack(content)
    assert {"run.py", "analyse.py", "pixi.toml", "config.json", "README.md"} <= set(files)

    payload = json.loads(files["config.json"])
    assert payload["config"]["job_name"] == "lysozyme_test"
    # Nothing internal should leak into the bundle's configuration.
    assert not any(k.startswith("_") for k in payload["config"])

    # The structure that reaches the bundle must be the one that was uploaded.
    assert files["input.pdb"] == (FIXTURES / "1aki.pdb").read_bytes()


def test_the_delivered_bundle_contains_valid_python(client):
    """The end of the pipeline the user actually touches."""
    import py_compile
    import tempfile

    _, token, _ = _upload(client)
    response = client.post("/prepare/build", data=_form_from_defaults(token))
    files = bundle.unpack(response.get_data())

    with tempfile.TemporaryDirectory() as tmp:
        for name in ("run.py", "analyse.py"):
            path = pathlib.Path(tmp) / name
            path.write_bytes(files[name])
            py_compile.compile(str(path), doraise=True)


def test_invalid_build_re_renders_the_form_with_errors(client):
    _, token, _ = _upload(client)
    response = client.post(
        "/prepare/build",
        data=_form_from_defaults(token, padding=0.5, nonbonded_cutoff=1.0),
    )
    html = response.get_data(as_text=True)
    assert response.status_code == 200
    assert "md-issue-error" in html
    assert "periodic image" in html
    # The form must come back populated rather than reset.
    assert html.count("data-option=") == len(opts.OPTIONS)


def test_build_rejects_a_job_name_that_is_a_path(client):
    _, token, _ = _upload(client)
    response = client.post("/prepare/build",
                           data=_form_from_defaults(token, job_name="../../escape"))
    assert response.mimetype != "application/json"
    assert "md-issue-error" in response.get_data(as_text=True)


def test_multiselect_values_survive_the_round_trip(client):
    _, token, _ = _upload(client, "4hhb.pdb")
    data = _form_from_defaults(token)
    data["chains"] = ["A", "B"]
    response = client.post("/prepare/build", data=data)
    payload = json.loads(bundle.unpack(response.get_data())["config.json"])
    assert payload["config"]["chains"] == ["A", "B"], \
        "a multiselect posted repeatedly must not collapse to its first value"


# ---------------------------------------------------------------------------
#  The estimate API
# ---------------------------------------------------------------------------


def test_estimate_returns_derived_numbers(client):
    response = client.post("/api/estimate", json=opts.defaults() | {
        "_estimated_atoms": 25000, "_solute_atoms": 2500,
    })
    assert response.status_code == 200
    data = response.get_json()
    assert data["derived"]["production_steps"] > 0
    assert data["wall"]["basis"] == "estimated"
    assert "traj_size_human" in data["derived"]


def test_estimate_reports_active_fields(client):
    """The active list is what the browser uses to show and hide fields."""
    membrane = client.post("/api/estimate", json=opts.defaults() | {"use_membrane": True})
    soluble = client.post("/api/estimate", json=opts.defaults() | {"use_membrane": False})
    assert "lipid_type" in membrane.get_json()["active"]
    assert "lipid_type" not in soluble.get_json()["active"]


def test_estimate_reports_validation_errors(client):
    response = client.post("/api/estimate", json=opts.defaults() | {"padding": 0.5})
    errors = response.get_json()["errors"]
    assert any(e["option"] == "padding" for e in errors)


@pytest.mark.parametrize("body", ["3", "[1,2,3]", '"a string"', "null", "true"])
def test_estimate_rejects_json_that_is_not_an_object(client, body):
    """json.loads("3") is an int, and .get on it is an AttributeError.

    A bare scalar or list is valid JSON, so the endpoint has to check the type
    rather than assume a dict came back.
    """
    response = client.post("/api/estimate", data=body,
                           content_type="application/json")
    assert response.status_code == 400
    assert "Traceback" not in response.get_data(as_text=True)


def test_estimate_survives_complete_nonsense(client):
    response = client.post("/api/estimate", json={
        "timestep": "not a number", "replicates": ["a", "list"], "bogus_key": 1,
    })
    assert response.status_code == 200
    assert "Traceback" not in response.get_data(as_text=True)


# ---------------------------------------------------------------------------
#  Scratch hygiene
# ---------------------------------------------------------------------------


def test_sessions_land_in_scratch_and_sweep_away(app, client, tmp_path):
    import os
    import time

    from flexappeal.webapp import sweep_scratch

    _, token, _ = _upload(client)
    root = tmp_path / "scratch"
    assert (root / token).is_dir()
    assert (root / token / "report.json").is_file()

    # Nothing expired yet.
    assert sweep_scratch(root, ttl=3600) == 0

    old = time.time() - 10_000
    os.utime(root / token, (old, old))
    assert sweep_scratch(root, ttl=3600) == 1
    assert not (root / token).exists()


# ---------------------------------------------------------------------------
#  Analysis tab
# ---------------------------------------------------------------------------


def _analyse(client, payload=None):
    data = {"fxa_file": (io.BytesIO(payload if payload is not None
                                    else (FIXTURES / "lysozyme.fxa").read_bytes()),
                         "results.fxa")}
    return client.post("/analysis", data=data, content_type="multipart/form-data")


def test_analysis_page_offers_an_upload(client):
    html = client.get("/analysis").get_data(as_text=True)
    assert "fxa_file" in html
    assert "md-drop" in html


def test_uploading_a_real_fxa_renders_every_panel(client):
    html = _analyse(client).get_data(as_text=True)
    assert "data-figure=" in html
    # RMSD, Rg, RMSF, SASA, secondary structure, DSSP, contacts, PCA, clusters,
    # plus the convergence multiples.
    assert html.count("data-figure=") >= 10
    assert 'class="md-tile"' in html
    assert "molstar-viewer" in html


def test_analysis_serves_the_viewer_assets(client):
    html = _analyse(client).get_data(as_text=True)
    token = html.split("/analysis/", 1)[1].split("/", 1)[0]

    structure = client.get(f"/analysis/{token}/structure")
    assert structure.status_code == 200
    assert structure.get_data().startswith((b"ATOM", b"REMARK", b"CRYST", b"MODEL", b"HEADER"))

    trajectory = client.get(f"/analysis/{token}/trajectory")
    assert trajectory.status_code == 200
    assert len(trajectory.get_data()) > 1000


def test_viewer_assets_reject_a_hostile_token(client):
    for hostile in ("../../etc", "/etc/passwd", "a" * 200):
        response = client.get(f"/analysis/{hostile}/structure")
        assert response.status_code in (400, 404)
        assert "Traceback" not in response.get_data(as_text=True)


def test_uploading_rubbish_is_a_friendly_error(client):
    html = _analyse(client, b"not a results file at all" * 50).get_data(as_text=True)
    assert "md-issue-error" in html
    assert "Traceback" not in html


def test_uploading_nothing_is_a_friendly_error(client):
    response = client.post("/analysis", data={}, content_type="multipart/form-data")
    assert "Choose a .fxa" in response.get_data(as_text=True)


def test_analysis_page_includes_a_table_view(client):
    """Identity is never colour-alone: every plotted fact is also readable as text."""
    html = _analyse(client).get_data(as_text=True)
    assert '<table class="md-table">' in html
    assert "Metrics computed" in html


def test_analysis_vendors_its_javascript(client):
    """No CDN: the page must work with no third-party network access."""
    html = _analyse(client).get_data(as_text=True)
    assert "/static/vendor/plotly.min.js" in html
    assert "/static/vendor/molstar.js" in html
    assert "cdn." not in html.split("<footer")[0]


# ---------------------------------------------------------------------------
#  Re-analysis routes
# ---------------------------------------------------------------------------


def _analysis_token(client):
    html = _analyse(client).get_data(as_text=True)
    return html.split("/analysis/", 1)[1].split("/", 1)[0]


def test_the_analysis_page_offers_reanalysis(client):
    html = _analyse(client).get_data(as_text=True)
    assert "reanalyse-card" in html
    assert 're-metric' in html
    # The resid/resSeq distinction is a genuine trap: resid is a zero-based
    # index, so `resid 195` is rarely residue 195.
    assert "resSeq" in html


def test_reanalysis_rejects_a_bad_request_without_starting_work(client):
    token = _analysis_token(client)
    response = client.post(f"/analysis/{token}/reanalyse",
                           json={"metrics": ["nonsense"], "selection": "protein"})
    assert response.status_code == 400
    assert "nonsense" in response.get_json()["message"]


def test_reanalysis_rejects_a_hostile_selection(client):
    token = _analysis_token(client)
    response = client.post(
        f"/analysis/{token}/reanalyse",
        json={"metrics": ["rmsd"], "selection": "__import__('os').system('id')"})
    assert response.status_code == 400


def test_reanalysis_refuses_a_hostile_token(client):
    for hostile in ("../../etc", "/etc/passwd", "a" * 200):
        response = client.post(f"/analysis/{hostile}/reanalyse",
                               json={"metrics": ["rmsd"], "selection": "protein"})
        assert response.status_code in (400, 404)
        assert "Traceback" not in response.get_data(as_text=True)


def test_reanalysis_status_is_idle_before_anything_runs(client):
    token = _analysis_token(client)
    response = client.get(f"/analysis/{token}/reanalyse/status")
    assert response.get_json()["status"] == "idle"


def test_reanalysis_is_limited_to_one_job(client, app):
    """The second request must be told to wait rather than starting work."""
    from flexappeal import analysis as analysis_module

    token = _analysis_token(client)
    analysis_module.acquire_lock(app.config["SCRATCH_ROOT"])
    try:
        response = client.post(f"/analysis/{token}/reanalyse",
                               json={"metrics": ["rgyr"], "selection": "protein"})
        assert response.status_code == 429
        assert response.get_json()["status"] == "busy"
    finally:
        analysis_module.release_lock(app.config["SCRATCH_ROOT"])


def test_reanalysis_status_reports_a_finished_job(client, app):
    """The status route reads what the detached worker wrote."""
    import json as json_module

    from flexappeal import analysis as analysis_module

    token = _analysis_token(client)
    session = pathlib.Path(app.config["SCRATCH_ROOT"]) / token
    request_path = session / "reanalyse_request.json"
    request_path.write_text(json_module.dumps(
        {"metrics": ["rgyr"], "selection": "protein"}))

    analysis_module.acquire_lock(app.config["SCRATCH_ROOT"])
    analysis_module.run_to_file(session / "results.fxa", request_path,
                                session / "reanalyse_result.json")

    payload = client.get(f"/analysis/{token}/reanalyse/status").get_json()
    assert payload["status"] == "ready"
    assert payload["metrics"]["rgyr_nm"]
    # The status route builds the figures, so the browser only draws.
    assert any(p["id"] == "re-rgyr" for p in payload["figures"])


# ---------------------------------------------------------------------------
#  Browser-side validity of the number inputs
# ---------------------------------------------------------------------------


def _number_inputs(html):
    """Every rendered number input as a dict of its attributes."""
    import re

    found = {}
    for tag in re.findall(r'<input type="number"[^>]*>', html):
        name = re.search(r'name="([^"]+)"', tag)
        if not name:
            continue
        attrs = dict(re.findall(r'(\w[\w-]*)="([^"]*)"', tag))
        found[name.group(1)] = attrs
    return found


def test_the_browser_accepts_every_default(client):
    """HTML5 validates a number input as `value == min + n*step`, anchored at min.

    A step that does not divide (default - min) makes the browser reject the
    form's own default with "Enter a valid value". Nine options did exactly
    that: production_duration had min 0.001 and step 1, so both the 100 ns
    default and a typed 10 ns were refused, and the form could not be submitted
    at all.
    """
    _, _, html = _upload(client)
    offenders = []
    for name, attrs in _number_inputs(html).items():
        step = attrs.get("step", "any")
        if step == "any":
            continue
        low = float(attrs.get("min", 0))
        n = (float(attrs["value"]) - low) / float(step)
        if abs(n - round(n)) > 1e-9:
            offenders.append(f"{name} (min={low}, step={step}, value={attrs['value']})")
    assert not offenders, "the browser would reject these defaults: " + "; ".join(offenders)


def test_float_inputs_do_not_constrain_to_a_grid(client):
    """schema.py accepts any value between minimum and maximum.

    The browser must not be stricter than the server, or a legitimate value is
    refused before it can even be submitted.
    """
    _, _, html = _upload(client)
    inputs = _number_inputs(html)
    for opt in opts.OPTIONS:
        if opt.widget == "number" and opt.id in inputs:
            assert inputs[opt.id].get("step") == "any", (
                f"{opt.id} pins the browser to a grid, but schema.py accepts any "
                f"value in range"
            )


def test_integer_inputs_step_by_one_from_a_whole_number(client):
    """step=1 is only exact if the minimum is itself an integer."""
    _, _, html = _upload(client)
    inputs = _number_inputs(html)
    for opt in opts.OPTIONS:
        if opt.widget == "int" and opt.id in inputs:
            assert inputs[opt.id].get("step") == "1"
            low = float(inputs[opt.id].get("min", 0))
            assert low == int(low), f"{opt.id} has a fractional minimum with step=1"


@pytest.mark.parametrize("value", ["10", "0.5", "137.4", "1", "100"])
def test_a_typed_production_time_is_accepted_end_to_end(client, value):
    """The reported bug: typing 10 into Production and being told it is invalid."""
    _, token, _ = _upload(client)
    response = client.post("/prepare/build",
                           data=_form_from_defaults(token, production_duration=value))
    assert response.mimetype == "application/octet-stream", \
        f"production_duration={value} was rejected: " \
        f"{response.get_data(as_text=True)[:200]}"


def test_the_page_does_not_promise_double_click(client):
    """A file downloaded over HTTP has no executable bit, so Finder refuses it:

        The file "flexappeal_1EVS.command" could not be executed because you do
        not have appropriate access privileges.

    The page told people to double-click it anyway.
    """
    _, _, html = _upload(client)

    # Not a ban on the word: the page should still explain that double-clicking
    # fails. What it must never do is offer it as the way to run the bundle, so
    # every mention has to sit in a sentence that says it does not work.
    import re

    for sentence in re.split(r"(?<=[.!?])\s+", re.sub(r"<[^>]+>", " ", html)):
        if "ouble-click" in sentence:
            assert re.search(r"\bfail|will not work|cannot|does not work", sentence), (
                f"double-click is presented as a method, not a warning: "
                f"{' '.join(sentence.split())[:140]}"
            )


def test_the_page_gives_the_commands_that_actually_work(client):
    _, _, html = _upload(client)
    assert "chmod +x" in html
    assert "cd ~/Downloads" in html
    assert "not optional" in html, "the chmod has to read as required, not advisory"


def test_the_commands_name_the_file_that_will_be_downloaded(client):
    """The block is meant to be copied verbatim, so the filename must be real."""
    from flexappeal import bundle as bundle_mod

    _, token, html = _upload(client)
    # The page seeds the job name from the uploaded filename.
    assert 'name="job_name" value="1aki"' in html
    expected = bundle_mod.bundle_filename({"job_name": "1aki"})

    response = client.post("/prepare/build", data=_form_from_defaults(token, job_name="1aki"))
    assert expected in response.headers["Content-Disposition"]
    # And the page's JS builds the same name from the same rule.
    assert "'flexappeal_' + safe + '.command'" in (
        (pathlib.Path(__file__).parent.parent / "flexappeal" / "static" / "app.js").read_text())


# ---------------------------------------------------------------------------
#  Banner
# ---------------------------------------------------------------------------


def _banner_lines():
    """The banner art, unescaped.

    It is stored HTML-escaped because the slant font draws x as ">  <", and an
    unescaped "</" makes the browser close the <pre> and swallow the rest.
    """
    import html as html_mod

    from flexappeal.webapp import PACKAGE_ROOT

    markup = (PACKAGE_ROOT / "templates" / "_banner.html").read_text()
    art = markup.split("<pre>", 1)[1].split("</pre>", 1)[0]
    return html_mod.unescape(art).split("\n")


def test_the_banner_art_is_html_escaped():
    """slant renders x as ">  <"; unescaped, the browser eats the banner."""
    from flexappeal.webapp import PACKAGE_ROOT

    markup = (PACKAGE_ROOT / "templates" / "_banner.html").read_text()
    art = markup.split("<pre>", 1)[1].split("</pre>", 1)[0]
    assert "<" not in art and ">" not in art, \
        "raw angle brackets in the art will be parsed as markup"


def test_both_landing_pages_show_the_banner(client):
    for page in ("/", "/analysis"):
        html = client.get(page).get_data(as_text=True)
        assert "md-banner" in html, f"{page} has no banner"


def test_the_banner_matches_the_one_the_bundle_prints(client):
    """The site and the downloaded .command should look like one tool."""
    from flexappeal.webapp import PACKAGE_ROOT

    bootstrap = (PACKAGE_ROOT / "runtime" / "bootstrap.sh.j2").read_text()
    shell_art = bootstrap.split("cat <<'BANNER'\n", 1)[1].split("\nBANNER", 1)[0]
    assert shell_art.split("\n") == _banner_lines(), \
        "the page banner and the shell banner have drifted apart"


def test_the_banner_art_is_well_formed():
    """Structural sanity, independent of which figlet font is in use.

    Deliberately not pinned to one font's dimensions: the first version of this
    test hard-coded smslant's five lines and 46 columns, and broke the moment
    the font changed for a legibility reason rather than a correctness one.
    """
    lines = _banner_lines()

    assert 4 <= len(lines) <= 8, f"{len(lines)} lines is not a figlet banner"
    assert max(len(l) for l in lines) <= 60, "too wide for an 80-column terminal"
    assert all(l == l.rstrip() for l in lines), "trailing whitespace in the art"
    assert not any("\t" in l for l in lines), "a tab will not align anywhere"

    # Flush left: at least one row must reach column zero, or the whole block
    # is floating on indentation nobody intended.
    #
    # Note what this deliberately does NOT assert. The previous version required
    # *exactly* one row at column zero, which is true of an italic face like
    # slant (each row steps left by one) and false of an upright one like
    # standard (most rows start at zero). That is the second time this test
    # encoded the incumbent font's geometry as if it were a correctness
    # property and then failed on a deliberate font change. Shape is a design
    # decision; only structural sanity belongs here.
    assert any(l and not l.startswith(" ") for l in lines), \
        "no row reaches column 0 -- the art is floating on stray indentation"

    # And the glyph rows must actually carry glyphs.
    assert all(l.strip() for l in lines[:-1]), "a blank row inside the art"


def test_the_banner_is_hidden_from_screen_readers(client):
    html = client.get("/").get_data(as_text=True)
    banner = html.split('class="md-banner"', 1)[1][:400]
    assert 'aria-hidden="true"' in html.split("md-banner", 1)[0][-80:] or \
           'aria-hidden' in html.split('class="md-banner"', 1)[0][-80:] or \
           'aria-hidden="true"' in banner or 'aria-hidden' in html
    assert "md-visually-hidden" in html, "no text alternative for the art"


def test_the_tab_row_is_centred(client):
    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    assert ".md-tabs-inner { justify-content: center; }" in css \
        or "justify-content: center" in css.split(".md-tabs-inner")[1][:200]


def test_html_is_never_heuristically_cached(client):
    """The counterpart to the ?v=mtime on static assets.

    Static files are immutable for a year, so the only thing that can deliver a
    changed stylesheet is a fresh copy of the HTML that names it. With no
    Cache-Control at all the browser invents its own freshness lifetime and
    holds the page -- and with it the stale asset URL.
    """
    for path in ("/", "/analysis"):
        r = client.get(path)
        assert r.status_code == 200
        assert "no-cache" in r.headers["Cache-Control"], path


def test_static_assets_still_carry_a_cache_busting_version(client):
    """no-cache on HTML must not have been applied to static files too."""
    body = client.get("/").get_data(as_text=True)
    assert "brand.css?v=" in body


def test_the_banner_font_size_is_always_a_whole_number_of_pixels():
    """Monospace art only stays in column at an integer font-size.

    A vw-derived size is fractional at almost every window width (13.05px at
    900px, 14.5px at 1000px). The character advance is then fractional too, the
    browser rounds each glyph origin to the device pixel grid independently,
    and the error accumulates differently per row -- so the columns drift and
    the strokes stop meeting. It is invisible at 1200px, where the old clamp
    happened to land on a clean 16px, which is why it survived several rounds
    of screenshot checks.
    """
    import re

    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    block = css[css.index(".md-banner pre"):]

    for size in re.findall(r"font-size:\s*([^;]+);", block[:block.index("}")] ):
        assert re.fullmatch(r"\d+px", size.strip()), \
            f"banner font-size {size.strip()!r} is not a whole number of pixels"

    # ...including the step-down breakpoint.
    for size in re.findall(r"\.md-banner pre\s*{\s*font-size:\s*([^;]+);", css):
        assert re.fullmatch(r"\d+px", size.strip()), size


def test_the_banner_pre_overrides_the_inherited_centring():
    """.md-banner centres the block; the <pre> must not inherit that.

    text-align inherits. With display:inline-block the box is sized by its
    longest line, but each line is then centred inside it independently, so
    every row shorter than the longest slides right by half the difference.
    The banner's two descender rows are 40 columns against 55, which pushed the
    Ps' legs 7.5 columns right into the following letters -- the long-running
    "legs are shifted right" bug, which survived six font changes because it
    was never in the art.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    block = css[css.index(".md-banner pre"):]
    block = block[:block.index("}")]
    assert "text-align: left" in block, \
        ".md-banner pre must reset the centring it inherits from .md-banner"


def test_banner_rows_are_not_all_the_same_length():
    """Guards the premise of the test above.

    If every row were the longest, centring would be a no-op and the reset
    would look like dead code to a future reader. It is not: the descender rows
    are genuinely shorter.
    """
    lines = _banner_lines()
    assert len(set(len(l) for l in lines)) > 1


def test_the_banner_is_never_hidden_on_narrow_screens():
    """It used to be display:none below 560px; a phone should still see it.

    That rule was written for 54-column art that genuinely could not scale down
    far enough. smslant is 45 columns and reaches a 280px viewport at 8px, so
    hiding it is no longer the trade.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    for chunk in css.split(".md-banner")[1:]:
        head = chunk[:chunk.index("}")] if "}" in chunk else chunk
        assert "display: none" not in head, "the banner must stay visible on mobile"


def test_every_banner_font_size_step_fits_its_viewport():
    """The narrow end of each breakpoint must still fit on screen.

    Character advance is 0.6021em for the banner's monospace stack (measured in
    Chrome, not assumed). Available width is the viewport less .md-main's
    padding, which is 16px a side below 768px. Each rule is valid down to the
    width where the next one takes over, so it is that width that has to fit.
    """
    import re

    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    steps = [(int(w), int(f)) for w, f in re.findall(
        r"@media \(max-width: (\d+)px\) \{ \.md-banner pre \{ font-size:\s*(\d+)px", css)]
    assert steps, "the breakpoint ladder has gone missing"

    cols = max(len(l) for l in _banner_lines())
    # Where the next rule takes over. The last step has no successor, so it is
    # bounded at 253px -- the narrowest viewport the ladder is designed for, and
    # far below any shipping phone. Under that, .md-banner's overflow-x: auto is
    # the fallback rather than a broken layout.
    lower = [s[0] + 1 for s in steps[1:]] + [253]

    for (bp, size), floor in zip(steps, lower):
        width = cols * 0.6021 * size
        available = floor - 32
        assert width <= available, (
            f"at {floor}px the {size}px step needs {width:.0f}px "
            f"but only {available}px is available")


# ---------------------------------------------------------------------------
#  The worked example
# ---------------------------------------------------------------------------


def _example_built():
    from flexappeal.webapp import EXAMPLE_FXA

    return EXAMPLE_FXA.is_file()


def test_the_example_tab_is_in_the_navigation(client):
    """Three tabs, and the third points at the example."""
    html = client.get("/").get_data(as_text=True)
    assert html.count("md-tab-num") == 3
    assert "/example" in html


def test_the_example_page_answers_either_way(client):
    """Built or not, it must be a page rather than a 404 or a traceback.

    The results file is large enough that a shallow or partial clone can be
    missing it, and a stack trace is a poor way to find that out.
    """
    response = client.get("/example")
    assert response.status_code == 200
    assert "Traceback" not in response.get_data(as_text=True)


@pytest.mark.skipif(not _example_built(), reason="the example run has not been built")
def test_the_example_renders_the_real_analysis(client):
    """The page must be the Analysis view, not a description of it.

    The whole point of committing a real run is that this page shares a code
    path with an upload, so it cannot quietly drift from the product. If these
    stop appearing, the page has become a brochure.
    """
    html = client.get("/example").get_data(as_text=True)
    assert "md-plot" in html, "no Plotly panels"
    assert "molstar-viewer" in html, "no structure viewer"
    assert "A complete run, start to finish" in html
    assert "Every setting this run used" in html


@pytest.mark.skipif(not _example_built(), reason="the example run has not been built")
def test_the_example_lists_every_option_from_the_registry(client):
    """The tables are generated, not written out, so they cannot fall behind.

    A hand-maintained table of 115 options would be stale within a release.
    """
    from flexappeal import options as opts

    html = client.get("/example").get_data(as_text=True)
    missing = [o.label for o in opts.OPTIONS if o.label not in html]
    assert not missing, f"options absent from the example tables: {missing[:5]}"


@pytest.mark.skipif(not _example_built(), reason="the example run has not been built")
def test_the_example_run_actually_went_the_distance(client):
    """Guards against committing a truncated or aborted run as the reference."""
    from flexappeal import fxa
    from flexappeal.webapp import EXAMPLE_FXA

    results = fxa.load(EXAMPLE_FXA.read_bytes(), verify_checksums=False)
    config = results.manifest.get("config", {})
    assert config.get("production_duration") == 10.0, "the example is meant to be 10 ns"

    times = results.metrics.get("time_ns") or []
    assert times, "no time axis in the results"
    assert times[-1] >= 9.5, f"the run stopped at {times[-1]:.2f} ns, short of 10"


@pytest.mark.skipif(not _example_built(), reason="the example run has not been built")
def test_the_terminal_captures_are_real_output_not_placeholders(client):
    """Each capture must carry styling from the recording and its own classes.

    rich numbers its CSS classes from zero in every export, so several captures
    on one page collide unless namespaced -- and a class may not begin with a
    digit, which an earlier version of the generator got wrong, silently
    dropping every colour.
    """
    import re

    from flexappeal.webapp import EXAMPLE_DIR

    captures = sorted((EXAMPLE_DIR / "screenshots").glob("*.html"))
    assert captures, "no terminal captures were generated"

    # Both halves of what a user does, grouped by the phase in the filename.
    phases = {c.name.split("-", 1)[0] for c in captures}
    assert phases == {"run", "analyse"}, f"expected both phases, found {phases}"

    # These go on a public page; the run summary prints its output directory.
    for capture in captures:
        assert "/Users/" not in capture.read_text(), \
            f"{capture.name} leaks a home directory onto a public page"
    seen = set()
    for capture in captures:
        text = capture.read_text()
        classes = set(re.findall(r"\.([A-Za-z][\w-]*)-r\d+", text))
        assert classes, f"{capture.name} has no styling"
        for name in classes:
            assert name[0].isalpha(), f"{name!r} is not a valid CSS class"
        assert not (classes & seen), f"{capture.name} reuses another capture's classes"
        seen |= classes


@pytest.mark.skipif(not _example_built(), reason="the example run has not been built")
def test_the_example_survives_its_session_being_swept(client, app):
    """The scratch sweeper must not be able to break this page.

    The example caches a session token so every visitor is not given their own
    copy of a five-megabyte results file. That cache outlives the directory it
    names -- the sweeper removes anything older than the TTL -- so the second
    visit after a sweep has to recreate it silently. It used to abort with 400
    instead, which meant the page worked for four hours after a deploy and then
    stopped, on a timer, with nothing in the logs pointing at the cause.
    """
    import shutil

    from flexappeal import webapp

    assert client.get("/example").status_code == 200

    root = pathlib.Path(app.config["SCRATCH_ROOT"])
    for entry in root.iterdir():
        shutil.rmtree(entry, ignore_errors=True)
    assert webapp._EXAMPLE_TOKEN, "the token should still be cached"

    assert client.get("/example").status_code == 200, "did not recover from a sweep"


def test_the_viewer_autoplays_the_trajectory_for_fifty_seconds():
    """The structure panel starts playing on load, one pass in 50 s.

    Mol*'s parameter is durationInS and it clamps to 1-120, so this is seconds
    rather than milliseconds -- passing 50000 would silently clamp to 120 and
    the loop would take two minutes.

    The animation is looked up by name rather than imported: the viewer bundle
    does not export AnimateModelIndex on the global, and a version that renames
    it should leave a static structure with working manual controls rather than
    throwing inside the load chain and losing the viewer.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    js = (PACKAGE_ROOT / "static" / "analysis.js").read_text()
    assert "PLAYBACK_SECONDS = 50" in js
    assert "built-in.animate-model-index" in js
    assert "durationInS: PLAYBACK_SECONDS" in js, "the duration must be in seconds"
    assert "durationInMs" not in js, "durationInMs is the wrong parameter and clamps"

    body = js[js.index("function autoplay"):]
    body = body[: body.index("\n  function ")]
    assert "prefers-reduced-motion" in body, \
        "a scene that starts moving on its own must honour reduced motion"
    assert "catch" in body, "a failed lookup must not take the viewer down with it"


def test_the_vendored_molstar_still_registers_the_animation():
    """Guards the name this depends on across a Mol* upgrade.

    If a new bundle renames it the page degrades quietly to a static structure,
    which is the right behaviour but the wrong thing to discover in production.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    bundle = (PACKAGE_ROOT / "static" / "vendor" / "molstar.js").read_text(errors="replace")
    assert "built-in.animate-model-index" in bundle
    assert "durationInS" in bundle, "the duration parameter has been renamed"


def test_the_top_panel_is_full_width_on_every_tab():
    """The opening card fills the content column, as it always did on mobile.

    It used to carry max-width: 760px with auto margins, which only ever took
    effect on desktop -- below 768px the cap never binds. So the two layouts
    disagreed for no reason anyone chose. The outside padding comes from
    .md-main (24px, 16px under 768px), not from the card, so removing the cap
    cannot push it against the viewport edge.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    rules = [line for line in css.splitlines()
             if line.strip().startswith(".md-intro") and "{" in line]
    assert rules, ".md-intro has gone missing entirely"
    for rule in rules:
        assert "max-width" not in rule, f"the top panel is still capped: {rule.strip()}"

    # .md-main is what supplies the outside padding and the overall cap.
    assert ".md-main { max-width: 1400px; margin: 0 auto; padding: 24px; }" in css


def test_every_tab_uses_that_same_top_panel(client):
    """Prepare, Analysis and Example must not drift apart on this."""
    for path in ("/", "/analysis", "/example"):
        html = client.get(path).get_data(as_text=True)
        assert "md-card md-intro" in html, f"{path} has no md-intro top panel"


# ---------------------------------------------------------------------------
#  Trajectory clip
# ---------------------------------------------------------------------------


def _video_tag(html):
    return "<video" + html.split("<video", 1)[1].split(">", 1)[0]


def test_every_tab_shows_the_clip_beside_its_top_panel(client):
    """The 75/25 split is the same on all three, from one included partial."""
    for path in ("/", "/analysis", "/example"):
        html = client.get(path).get_data(as_text=True)
        assert 'class="md-split"' in html, f"{path} does not split its top panel"
        assert "md-clip-card" in html, f"{path} has no clip beside the top panel"


def test_the_split_is_three_to_one(client):
    """75/25 as fractions, so the gap comes out of the grid, not the columns."""
    from flexappeal.webapp import PACKAGE_ROOT

    css = (PACKAGE_ROOT / "static" / "brand.css").read_text()
    assert "grid-template-columns: 3fr 1fr;" in css


def test_the_clip_is_not_fetched_with_the_page(client):
    """Decorative, so it must not compete with the content for bandwidth.

    The markup carries a poster and no source at all; clip.js attaches the real
    file after window.load. A plain src= or a <source> here would put nearly a
    megabyte of video back on the critical path, which is invisible on a fast
    connection and the whole cost of the panel on a slow one.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    for path in ("/", "/analysis", "/example"):
        html = client.get(path).get_data(as_text=True)
        tag = _video_tag(html)
        assert "<source" not in html.split("</video", 1)[0].split("<video", 1)[1]
        assert " src=" not in tag, f"{path} loads the clip eagerly: {tag}"
        assert "data-src=" in tag and "poster=" in tag

    clip_js = (PACKAGE_ROOT / "static" / "clip.js").read_text()
    assert "window.addEventListener('load'" in clip_js
    assert "prefers-reduced-motion" in clip_js, "reduced motion must skip the download"


def test_the_clip_stays_small():
    """Re-encoding it is easy; re-encoding it carelessly is easier.

    811 KB at the time of writing: 480x552, 20 fps, H.264 with -tune animation,
    which beat both VP9 and AV1 at matched SSIM on this flat cel-style artwork.
    The budget is here so a future re-export cannot quietly land a 15 MB file
    on every page of the site.
    """
    from flexappeal.webapp import PACKAGE_ROOT

    clip = PACKAGE_ROOT / "static" / "lysozyme.mp4"
    poster = PACKAGE_ROOT / "static" / "lysozyme-poster.webp"
    assert clip.stat().st_size < 900_000, f"the clip has grown to {clip.stat().st_size:,} bytes"
    assert poster.stat().st_size < 40_000, "the poster is meant to paint instantly"
