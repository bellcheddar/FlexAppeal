"""Section 10's page QC: a real browser at both viewports, and axe-core.

Separate from the other test files because these are the only tests that need a
browser, and a suite that cannot run without one is a suite people stop running.
Everything here skips cleanly when Chromium is not installed, so `pixi run test`
works on a fresh checkout and `pixi run playwright-browser` is what turns these
on.

What this catches that a server-side check cannot: a JavaScript exception, a
plot that throws on an empty series, a contrast ratio, focus order, and an ARIA
relationship that only resolves once the page is running. Those are most of the
reason section 10 asks for a browser at all.

The two viewports are the spec's: 1440x900 and 390x844.
"""

from __future__ import annotations

import json
import socket
import threading
from pathlib import Path

import pytest

VIEWPORTS = {"desktop": {"width": 1440, "height": 900},
             "mobile": {"width": 390, "height": 844}}

# Section 10 wants zero serious or critical violations. Moderate and minor are
# reported rather than failed, because a page can carry a defensible minor
# finding and a test that fails on those gets disabled.
BLOCKING_IMPACTS = {"serious", "critical"}

# Console noise that is not a page fault. Kept deliberately short: the point of
# the check is that it fails on real errors, and a long allowlist is how that
# stops being true.
IGNORED_CONSOLE = (
    # Mol* reports WebGL unavailability in headless Chromium, which is a
    # property of the test environment rather than of the page.
    "WebGL",
    "webgl",
)


def _free_port() -> int:
    with socket.socket() as sock:
        sock.bind(("127.0.0.1", 0))
        return int(sock.getsockname()[1])


@pytest.fixture(scope="module")
def chromium():
    """A browser, or a skip if one is not installed."""
    playwright = pytest.importorskip("playwright.sync_api")
    with playwright.sync_playwright() as driver:
        try:
            browser = driver.chromium.launch()
        except Exception as exc:                                # noqa: BLE001
            pytest.skip(f"Chromium is not installed: "
                        f"run `pixi run playwright-browser` ({str(exc)[:60]})")
        yield browser
        browser.close()


@pytest.fixture(scope="module")
def site(tmp_path_factory):
    """The real app on a real port, in a thread.

    A live server rather than Flask's test client, because the whole point of
    these tests is the part that runs in a browser, and the test client never
    serves a page to one.
    """
    from werkzeug.serving import make_server

    from flexappeal.refstates import store
    from flexappeal.webapp import create_app

    root = tmp_path_factory.mktemp("site")
    database = root / "refstates.db"
    payload = {
        "name": "fixture", "kind": "kinase",
        "core": {"n_residues": 43, "source": "fixture core"},
        "motif_residues": {"xDFG": [380, 381, 382, 383],
                           "hinge": [316, 317, 318]},
        "ensembles": {"active": {"codes": ["2F4J"]},
                      "inactive": {"codes": ["3PYY"]}},
        "frames": [
            {"time_ns": t / 10.0, "n_conflicts": 0, "motifs": [
                {"motif": "xDFG", "d_act_a": 1.0, "d_inact_a": 4.0,
                 "s": 0.6, "label": "active"},
                {"motif": "hinge", "d_act_a": 2.0, "d_inact_a": 2.0,
                 "s": 0.0, "label": "intermediate"}]}
            for t in range(12)],
        "notes": [],
    }
    store.ingest_analysis(database, payload, system_id="abl1", run_id="abl1_apo",
                          family="kinase", environment="water",
                          run_fields={"start_state": "DFG-in", "status": "complete",
                                      "ns_per_day": 213.0})

    # Reference geometry, so the state map has points and actually renders.
    # Without it the panel is correctly hidden, and a test that then asserts
    # the map drew is testing the fixture rather than the page.
    with store.connect(database) as connection:
        for code, state, alphac, lys in (("2F4J", "in", 6.83, 13.63),
                                         ("2GQG", "in", 6.64, 14.04),
                                         ("2HYY", "out", 13.82, 11.14),
                                         ("3PYY", "out", 13.10, 9.64)):
            store.put_reference(
                connection, f"abl1:{code}", "abl1", pdb=code, chain="A",
                source="KLIFS", source_state=state, our_state=state,
                metrics={"dfg_d_alphac": alphac, "dfg_d_lys": lys})

    # The run's analysis file, which the DSSP ribbon and the viewer read.
    analysis = root / "analysis"
    analysis.mkdir(exist_ok=True)
    (analysis / "abl1_apo.json").write_text(json.dumps({
        "dssp": {"n_frames": 12, "residues": [316, 317, 318, 380, 381],
                 "frames": ["HHHEE"] * 12, "motifs": {}}}))

    app = create_app()
    app.config["REFSTATES_DB"] = str(database)

    import flexappeal.refstates.web as web

    web.analysis_file = lambda run_id: analysis / f"{run_id}.json"
    port = _free_port()
    server = make_server("127.0.0.1", port, app, threaded=True)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{port}"
    server.shutdown()
    thread.join(timeout=5)


PAGES = [
    ("prepare", "/"),
    ("analysis", "/analysis"),
    ("example", "/example"),
    ("runs", "/runs"),
    ("refstates", "/refstates/"),
    ("refstates-about", "/refstates/about"),
    ("refstates-run", "/refstates/run/abl1_apo"),
    ("refstates-references", "/refstates/references/abl1"),
]


def _open(browser, site: str, path: str, viewport: dict):
    """One page, with every console error and page error collected."""
    context = browser.new_context(viewport=viewport)
    page = context.new_page()
    problems: list[str] = []
    page.on("console", lambda message: problems.append(
        f"console {message.type}: {message.text}")
        if message.type == "error" else None)
    page.on("pageerror", lambda error: problems.append(f"pageerror: {error}"))
    response = page.goto(f"{site}{path}", wait_until="networkidle")
    return context, page, response, problems


@pytest.mark.parametrize("name,path", PAGES)
@pytest.mark.parametrize("viewport", list(VIEWPORTS), ids=list(VIEWPORTS))
def test_every_page_loads_without_console_errors(chromium, site, name, path,
                                                 viewport, tmp_path_factory):
    context, page, response, problems = _open(
        chromium, site, path, VIEWPORTS[viewport])
    try:
        assert response is not None and response.status == 200, \
            f"{path} returned {response.status if response else 'nothing'}"
        real = [p for p in problems
                if not any(ignore in p for ignore in IGNORED_CONSOLE)]
        assert not real, f"{path} at {viewport}: {real}"

        shots = Path(tmp_path_factory.getbasetemp()) / "screenshots"
        shots.mkdir(exist_ok=True)
        page.screenshot(path=str(shots / f"{name}-{viewport}.png"),
                        full_page=True)
    finally:
        context.close()


@pytest.mark.parametrize("name,path", PAGES)
def test_no_serious_or_critical_accessibility_violations(chromium, site, name,
                                                         path):
    """axe-core at the desktop viewport. Section 10 wants zero serious or critical.

    Run at one viewport rather than both: axe's serious rules here are
    properties of the markup, which does not change between the two, and a
    second run doubles the time for findings that are the same.
    """
    from axe_playwright_python.sync_playwright import Axe

    context, page, response, _problems = _open(
        chromium, site, path, VIEWPORTS["desktop"])
    try:
        assert response is not None and response.status == 200
        results = Axe().run(page)
        blocking = [v for v in results.response.get("violations", [])
                    if v.get("impact") in BLOCKING_IMPACTS]
        if blocking:
            detail = json.dumps(
                [{"id": v["id"], "impact": v["impact"],
                  "help": v.get("help"),
                  "nodes": [n.get("target") for n in v.get("nodes", [])][:3]}
                 for v in blocking], indent=2)
            pytest.fail(f"{path}: {len(blocking)} blocking violations\n{detail}")
    finally:
        context.close()


def test_the_reference_state_plots_actually_draw(chromium, site):
    """The three Plotly panels have to produce a plot, not an empty div.

    A plot that throws leaves the div in place and the page looks fine, which is
    exactly the failure a status-code check cannot see. Asserted by counting the
    SVG Plotly injects rather than by looking for the div.
    """
    context, page, response, problems = _open(
        chromium, site, "/refstates/run/abl1_apo", VIEWPORTS["desktop"])
    try:
        assert response is not None and response.status == 200
        page.wait_for_selector("#state-timeline .main-svg", timeout=15000)
        page.wait_for_selector("#state-map .main-svg", timeout=15000)
        for panel in ("#state-timeline", "#dssp-ribbon", "#state-map"):
            count = page.locator(f"{panel} .main-svg").count()
            assert count > 0, f"{panel} drew nothing"
        real = [p for p in problems
                if not any(ignore in p for ignore in IGNORED_CONSOLE)]
        assert not real, real
    finally:
        context.close()
