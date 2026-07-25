"""Shared pytest configuration.

Puts the repo root on sys.path so `import flexappeal` works without installing
the package, matching BoltzMaker's tests/conftest.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
EXAMPLES = REPO_ROOT / "examples"

if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
