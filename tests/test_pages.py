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


# Mol*'s symptom when the GPU it is sharing is busy. The projection matrix
# inverts to NaN, the NaN reaches a typed array length as -Infinity, and the
# page logs an uncaught RangeError that has nothing to do with the page.
CONTENDED_WEBGL = "Invalid typed array length: -Infinity"


def _production_holds_the_gpu() -> bool:
    """Whether an MD run is on the device right now.

    `tests/test_refstates.py` already skips its resume check on this, for the
    same machine and the same reason. Here it decides whether a single Mol*
    symptom is worth failing on, because OpenMM and headless Chromium are
    asking the one GPU for work at the same time and the browser is the one
    that loses.

    It is checked rather than assumed, and it produces a skip rather than a
    filter in `IGNORED_CONSOLE`. A filter is permanent and silent: it would
    have hidden the next real RangeError for as long as it stayed there, and
    nothing would have said these checks had stopped covering the viewer. A
    skip says so, every run, and goes away on its own when the queue drains.
    """
    import subprocess

    return subprocess.run(["pgrep", "-f", "tools/(queue|run)_md\\.py"],
                          capture_output=True).returncode == 0


def _assert_console_clean(problems, where: str) -> None:
    real = [p for p in problems
            if not any(ignore in p for ignore in IGNORED_CONSOLE)]
    contended = [p for p in real if CONTENDED_WEBGL in p]
    if contended and len(contended) == len(real) and _production_holds_the_gpu():
        pytest.skip(f"{where}: a production run holds the GPU, so Mol* could "
                    f"not get a usable WebGL context. Everything else on the "
                    f"page was clean.")
    assert not real, f"{where}: {real}"


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


def _write_bundle(directory, stem: str, n_residues: int, n_frames: int,
                  *, resname: str = "ALA", first_residue: int = 1) -> None:
    """A viewer bundle a real Mol* can actually parse.

    Generated here rather than committed, and from nothing but a residue
    count, so it carries no structure from anyone's work.

    It replaces a placeholder PDB line and 512 zero bytes called an XTC. Those
    had never been parsed: the page resolved its viewer directory itself while
    the fixture patched only `analysis_file`, so every browser test had been
    loading the repository's own production trajectories. On a clean checkout
    there would have been nothing there to load, and the console-error check
    would have been asserting against a viewer that never started.

    `first_residue` exists because the numbering has to agree with the motifs
    the same fixture declares. Numbered from 1 against motifs in the 300s, the
    cartoon's motif paint selects no atoms, and the minimum over an empty
    selection reaches Mol* as a typed array of length -Infinity. That is a
    fixture that disagrees with itself, not a page fault, and it is the shape
    of thing these tests exist to catch, so it is worth the argument.
    """
    import mdtraj as md
    import numpy as np

    topology = md.Topology()
    chain = topology.add_chain()
    for index in range(n_residues):
        residue = topology.add_residue(resname, chain,
                                       resSeq=first_residue + index)
        topology.add_atom("CA", md.element.carbon, residue)
    # A gentle helix, so the thing on screen is a shape rather than a line,
    # and a small per-frame wobble so successive frames differ.
    turn = np.arange(n_residues) * 1.745
    base = np.stack([np.cos(turn) * 0.23,
                     np.sin(turn) * 0.23,
                     np.arange(n_residues) * 0.15], axis=1)
    xyz = np.stack([base + 0.01 * frame for frame in range(n_frames)])
    traj = md.Trajectory(xyz.astype(np.float32), topology)
    traj[0].save_pdb(str(directory / f"{stem}.pdb"))
    traj.save_xtc(str(directory / f"{stem}.xtc"))


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

    # A membrane row as well as the water one. Section 10 asks for every page
    # at both viewports with no console errors and no serious axe violation,
    # and until this fixture had a receptor the only run page ever opened in a
    # browser was a kinase in water. The bilayer table, its three verdict
    # badges and the lipid toggle are markup and script that exist only on a
    # membrane row, so none of them had been through either check.
    receptor = {
        "name": "fixture-gpcr", "kind": "gpcr",
        "core": {"n_residues": 180, "source": "fixture core"},
        "motif_residues": {"NPxxY": [322, 323, 324, 325, 326],
                           "TM6": [267, 268, 269, 270]},
        "ensembles": {"active": {"codes": ["3P0G"]},
                      "inactive": {"codes": ["2RH1"]}},
        "frames": [
            {"time_ns": t / 10.0, "n_conflicts": 0, "motifs": [
                {"motif": "NPxxY", "d_act_a": 0.4, "d_inact_a": 3.1,
                 "s": 0.77, "label": "active"},
                {"motif": "TM6", "d_act_a": 1.5, "d_inact_a": 3.4,
                 "s": 0.39, "label": "active"}]}
            for t in range(12)],
        "notes": [],
    }
    store.ingest_analysis(database, receptor, system_id="b2ar",
                          run_id="b2ar_apo", family="gpcr",
                          environment="membrane",
                          run_fields={"start_state": "active",
                                      "status": "complete",
                                      "ns_per_day": 59.7})

    # One chunk per verdict, so the badge macro is exercised in all three of
    # its branches rather than only the happy one. 67.0 is inside the expected
    # window, 69.5 is outside it but still acceptable, and 80.0 raises G6.
    with store.connect(database) as connection:
        store.put_membrane(connection, "b2ar_apo", [
            {"start_ns": 0.0, "apl_a2": 67.0, "thickness_a": 38.9,
             "tilt_deg": 13.6, "core_waters": 54.4},
            {"start_ns": 10.0, "apl_a2": 69.5, "thickness_a": 39.4,
             "tilt_deg": 15.4, "core_waters": 49.4},
            {"start_ns": 20.0, "apl_a2": 80.0, "thickness_a": 42.5,
             "tilt_deg": 9.6, "core_waters": 28.1},
        ])

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
    # Frames carrying the Dunbrack distances, so the state map has a
    # trajectory path to draw and the click test is testing the page rather
    # than the fixture. A drift across the plane, not a constant, so a click in
    # the middle lands somewhere distinguishable from either end.
    (analysis / "abl1_apo.json").write_text(json.dumps({
        "dssp": {"n_frames": 12, "residues": [316, 317, 318, 380, 381],
                 "frames": ["HHHEE"] * 12, "motifs": {}},
        "viewer": {"topology": "abl1_apo.pdb", "trajectory": "abl1_apo.xtc",
                   "n_frames": 6, "stride": 10, "frame_ns": 1.0,
                   "bytes": 1024},
        "frames": [
            {"time_ns": t / 2.0, "measurements": {"dfg_plane": {
                "status": "ok", "d_alphac_a": 6.5 + t * 0.05,
                "d_lys_a": 14.0 - t * 0.04}}}
            for t in range(12)],
    }))
    # The viewer panel only renders when both files exist, and the slider only
    # when the bundle has more than one frame.
    # Spanning both of this fixture's motifs, 316 to 318 and 380 to 383.
    _write_bundle(analysis, "abl1_apo", n_residues=80, n_frames=6,
                  first_residue=310)

    # The receptor's own bundle. `viewer.lipids` is what draws the bilayer
    # toggle, and the page only offers it when both lipid files are present,
    # so both are written.
    (analysis / "b2ar_apo.json").write_text(json.dumps({
        "dssp": {"n_frames": 12, "residues": [267, 268, 322, 323, 324],
                 "frames": ["HHHHH"] * 12, "motifs": {}},
        "viewer": {"topology": "b2ar_apo.pdb", "trajectory": "b2ar_apo.xtc",
                   "n_frames": 6, "stride": 10, "frame_ns": 1.0,
                   "bytes": 1024,
                   "lipids": {"topology": "b2ar_apo.lipids.pdb",
                              "trajectory": "b2ar_apo.lipids.xtc",
                              "n_atoms": 476, "n_frames": 6}},
        "frames": [{"time_ns": t / 2.0, "measurements": {}}
                   for t in range(12)],
    }))
    # Spanning this fixture's motifs, 267 to 270 and 322 to 326.
    _write_bundle(analysis, "b2ar_apo", n_residues=75, n_frames=6,
                  first_residue=260)
    # The lipid bundle is phosphorus atoms, which is what the real one holds
    # and what the bilayer toggle draws as spheres.
    _write_bundle(analysis, "b2ar_apo.lipids", n_residues=24, n_frames=6,
                  resname="POP")

    app = create_app()
    app.config["REFSTATES_DB"] = str(database)
    # Configuration rather than a monkeypatch. Patching `web.analysis_file`
    # redirected the run's JSON and nothing else, because the two routes that
    # serve a viewer bundle resolved the directory themselves, so a fixture
    # run's viewer was always looked for in the repository's own data and a
    # patch left behind leaked into every later test in the session. One key,
    # read the same way `REFSTATES_DB` is, covers all three.
    app.config["REFSTATES_ANALYSIS_DIR"] = str(analysis)
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
    ("refstates-run-membrane", "/refstates/run/b2ar_apo"),
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
        _assert_console_clean(problems, f"{path} at {viewport}")

        shots = Path(tmp_path_factory.getbasetemp()) / "screenshots"
        shots.mkdir(exist_ok=True)
        destination = str(shots / f"{name}-{viewport}.png")
        # The run page has grown tall enough that Chromium refuses a full-page
        # capture, so a viewport shot is the fallback. Both can fail: under a
        # loaded machine the capture itself returns "Unable to capture
        # screenshot" from the browser protocol, and with a 50 ns membrane run
        # on the GPU that happened to three of these sixteen while every one of
        # them passed in isolation.
        #
        # So the capture is best effort, which the previous comment already
        # argued for and then did not do. This test's assertion is that the
        # page loads without console errors. A missing artefact is worth
        # saying out loud and is not a reason to call a clean page broken,
        # which is the tail wagging the dog one level further down.
        captured, reason = None, "capture produced an empty file"
        for full_page in (True, False):
            try:
                page.screenshot(path=destination, full_page=full_page)
                captured = Path(destination)
                break
            except Exception as exc:                            # noqa: BLE001
                reason = str(exc)[:80]
        if captured is None or captured.stat().st_size <= 1000:
            print(f"\n  no screenshot for {name}-{viewport}: {reason}")
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
        _assert_console_clean(problems, "/refstates/run/abl1_apo")
    finally:
        context.close()


# ===========================================================================
#  Section 6.3's performance targets
# ===========================================================================

# Section 6.3: first paint under 2 s, frame scrub under 100 ms, per-system web
# bundle at or under 150 MB.
FIRST_PAINT_BUDGET_MS = 2000
SCRUB_BUDGET_MS = 100
BUNDLE_BUDGET_BYTES = 150 * 1024 * 1024


def test_the_explorer_paints_within_its_budget(chromium, site):
    """First contentful paint under 2 s, from the browser's own timing.

    Measured with `paint` entries rather than by timing the navigation call,
    because what section 6.3 asks about is when a reader sees something, not
    when the last asset finished arriving. The viewer and the plots load after
    first paint by design, so a slow Mol* does not count against this.
    """
    context = chromium.new_context(viewport=VIEWPORTS["desktop"])
    page = context.new_page()
    try:
        page.goto(f"{site}/refstates/run/abl1_apo", wait_until="domcontentloaded")
        # Waited for rather than read immediately. The paint entry is not always
        # recorded by the time DOMContentLoaded fires, and reading it at that
        # instant made this test fail intermittently on a page that paints in
        # about 400 ms. The wait does not inflate the measurement: the entry
        # carries its own startTime from the navigation.
        page.wait_for_function(
            """() => performance.getEntriesByType('paint')
                        .some(e => e.name === 'first-contentful-paint')""",
            timeout=FIRST_PAINT_BUDGET_MS * 3)
        paint = page.evaluate(
            """() => {
                const entries = performance.getEntriesByType('paint');
                const fcp = entries.find(e => e.name === 'first-contentful-paint');
                return fcp ? fcp.startTime : null;
            }""")
        assert paint is not None, "the browser reported no paint timing"
        assert paint < FIRST_PAINT_BUDGET_MS, (
            f"first contentful paint {paint:.0f} ms, budget "
            f"{FIRST_PAINT_BUDGET_MS} ms")
    finally:
        context.close()


def test_scrubbing_a_frame_stays_within_its_budget(chromium, site):
    """Under 100 ms per frame, measured on the series the timeline draws.

    Section 6.3 ties the budget to precomputed per-frame JSON, so what is timed
    is a frame change against data already in the page rather than a fetch. A
    scrub that had to go back to the server could not meet 100 ms and the
    design is that it does not have to.
    """
    context = chromium.new_context(viewport=VIEWPORTS["desktop"])
    page = context.new_page()
    try:
        page.goto(f"{site}/refstates/run/abl1_apo", wait_until="networkidle")
        page.wait_for_selector("#state-timeline .main-svg", timeout=15000)
        # Relayout the timeline to a single frame's window, twenty times, and
        # take the worst. Plotly's relayout is the same work a scrub does.
        worst = page.evaluate(
            """async () => {
                const node = document.getElementById('state-timeline');
                let worst = 0;
                for (let i = 0; i < 20; i++) {
                    const start = performance.now();
                    await Plotly.relayout(node, {'xaxis.range': [i * 0.1, i * 0.1 + 0.4]});
                    worst = Math.max(worst, performance.now() - start);
                }
                return worst;
            }""")
        assert worst is not None
        assert worst < SCRUB_BUDGET_MS, (
            f"worst frame change {worst:.0f} ms, budget {SCRUB_BUDGET_MS} ms")
    finally:
        context.close()


def test_the_web_bundle_is_within_its_budget():
    """Per-system bundle at or under 150 MB, measured on what is actually served.

    Everything the explorer fetches for one run: the page, the three JSON
    endpoints, the viewer topology and coordinates. Skipped rather than failed
    when a run has not been analysed, because the size is a property of the
    data and not of the code.
    """
    analysis = Path(__file__).resolve().parents[1] / "data" / "analysis"
    topology = analysis / "abl1_apo.pdb"
    if not topology.is_file():
        pytest.skip("no analysed run on disk to measure")

    served = sorted(p for p in analysis.glob("abl1_apo.*"))
    total = sum(p.stat().st_size for p in served)
    detail = ", ".join(f"{p.name} {p.stat().st_size / 1e6:.1f} MB" for p in served)
    assert total <= BUNDLE_BUDGET_BYTES, (
        f"bundle {total / 1e6:.1f} MB over the "
        f"{BUNDLE_BUDGET_BYTES / 1e6:.0f} MB budget: {detail}")


def test_clicking_the_state_map_moves_the_viewer_to_that_frame(chromium, site):
    """Section 6.2's "clicking a point jumps the viewer to that frame".

    Only a browser can check this, and it caught two faults a markup check
    passed: `const runId` declared below the fetch that used it, which threw a
    ReferenceError the fetch's own catch swallowed into "no reference map"; and
    a duration of n rather than (n - 1) strides, which made a frame at 4.2 ns
    read as 5.0.
    """
    context = chromium.new_context(viewport=VIEWPORTS["desktop"])
    page = context.new_page()
    try:
        page.goto(f"{site}/refstates/run/abl1_apo", wait_until="networkidle")
        page.wait_for_selector("#state-map .main-svg", timeout=20000)

        names = page.evaluate(
            "()=>document.getElementById('state-map').data.map(t=>t.name)")
        assert "this run" in names, f"the trajectory path is not drawn: {names}"

        if page.locator("#scrub").count() == 0:
            pytest.skip("this run has too few viewer frames for a slider")

        # The midpoint of the path, so neither end can pass by accident.
        clicked = page.evaluate("""()=>{
            const n = document.getElementById('state-map');
            const t = n.data.findIndex(d => d.name === 'this run');
            const i = Math.round((n.data[t].x.length - 1) / 2);
            n.emit('plotly_click', {points: [{curveNumber: t, pointNumber: i,
                customdata: n.data[t].customdata[i]}]});
            return n.data[t].customdata[i];
        }""")
        assert clicked, "the path carries no frame index"

        # The slider's value, not the readout. The readout is written by the
        # handler analysis.js attaches once the trajectory has loaded, so
        # asserting on it tests Mol* rather than this click handler, and it
        # fails whenever the structure cannot load. What the handler promises is
        # to move the slider and fire its input event; that is what is checked.
        moved = int(page.input_value("#scrub"))
        total = page.evaluate("""()=>{
            const n = document.getElementById('state-map');
            const t = n.data.findIndex(d => d.name === 'this run');
            return n.data[t].x.length;
        }""")
        frames = int(page.get_attribute("#scrub", "max")) + 1
        expected = round((int(clicked[0]) / (total - 1)) * (frames - 1))
        assert moved == expected, (
            f"clicked frame {clicked[0]} of {total}; slider went to {moved} "
            f"and should have gone to {expected}")
        assert moved > 0, "a click in the middle of the path left the slider at zero"
    finally:
        context.close()
