"""Shared pytest configuration.

Puts the repo root on sys.path so `import flexappeal` works without installing
the package, matching BoltzMaker's tests/conftest.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
EXAMPLES = REPO_ROOT / "examples"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(autouse=True)
def _keep_the_locale():
    """Put the C locale back after every test.

    **OpenMM changes it.** Constructing an `app.Simulation` loads the platform
    plugins, and that leaves `locale.getpreferredencoding()` reporting
    `US-ASCII` where it reported `UTF-8` a line earlier. Measured, not assumed:
    importing OpenMM does not do it, deserialising a system does not do it, and
    `app.Simulation(topology, system, integrator)` does.

    Everything downstream that reads a file without naming an encoding then
    uses ASCII and raises on the first non-ASCII byte. This repository is full
    of them, because its own prose contains β2AR, Cα and Å, so a single OpenMM
    test poisoned every later test in the session that read a source file or a
    template: 62 failures across the suite from one test, all of them
    `UnicodeDecodeError`.

    It stayed hidden for the whole build because the test that builds a
    `Simulation` skips while a production run holds the GPU, and a production
    run held the GPU continuously from the moment that test existed until the
    queue drained. The first suite run after Tier 2 finished was the first to
    execute it.

    Restoring the locale here rather than in that one test, because the next
    test to touch OpenMM would reintroduce it and the failure appears somewhere
    else entirely.
    """
    import locale

    saved = locale.setlocale(locale.LC_ALL)
    try:
        yield
    finally:
        try:
            locale.setlocale(locale.LC_ALL, saved)
        except locale.Error:                                    # noqa: PERF203
            pass
