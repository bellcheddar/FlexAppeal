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
