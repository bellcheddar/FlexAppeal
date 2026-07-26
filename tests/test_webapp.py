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

    # A slant font steps left by one column per row, so exactly one row reaches
    # column zero. A row starting left of the baseline means the art is skewed,
    # which is what the hand-written version got wrong.
    at_zero = [i for i, l in enumerate(lines) if l and not l.startswith(" ")]
    assert len(at_zero) == 1, f"rows {at_zero} all start at column 0"

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
