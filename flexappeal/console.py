"""The one place the toolchain's terminal output is styled.

FlexAppeal prints from three places -- this CLI, the run script inside a
generated bundle, and the analysis script beside it -- and they should look like
one program. The palette is the same #1e73be brand blue used by brand.css and
the shell installers, expressed here as rich styles.

The generated scripts cannot import this module: they run inside the bundle's
own pixi environment, which has no FlexAppeal package in it. They therefore
carry their own copy of THEME (see runtime/run.py.j2). Keep the two in step --
the test suite asserts they have not drifted.
"""

from __future__ import annotations

import shutil
import sys

from rich.console import Console
from rich.theme import Theme

# Truecolor, matching :root in brand.css exactly.
BRAND = {
    "primary": "#1e73be",
    "green": "#00d084",
    "amber": "#fcb900",
    "red": "#d63638",
    "muted": "#6b7c93",
}

THEME = Theme({
    "info": BRAND["primary"],
    "ok": BRAND["green"],
    "warn": BRAND["amber"],
    "err": BRAND["red"],
    "muted": BRAND["muted"],
    "step": BRAND["primary"],
    "heading": f"bold {BRAND['primary']}",
    "value": "bold",
    "unit": BRAND["muted"],
    # Progress bars. "bar.back" is the unfilled track; leaving it default makes
    # the bar invisible against a light terminal background.
    "bar.complete": BRAND["primary"],
    "bar.finished": BRAND["green"],
    "bar.pulse": BRAND["primary"],
    "progress.percentage": BRAND["muted"],
    "progress.elapsed": BRAND["muted"],
    "progress.remaining": BRAND["muted"],
})


def make_console(stderr: bool = False) -> Console:
    """A Console that behaves in a pipe as well as at a prompt.

    soft_wrap=False lets rich wrap to the real terminal width; when there is no
    terminal it falls back to 100 columns rather than rich's default 80, because
    the option tables the CLI prints are wider than that and truncating them in
    a redirected file is worse than a long line.
    """
    return Console(
        theme=THEME,
        stderr=stderr,
        width=None if sys.stdout.isatty() else max(100, shutil.get_terminal_size((100, 24)).columns),
    )


console = make_console()
err_console = make_console(stderr=True)


def info(msg: str) -> None:
    console.print(f"[info]ℹ[/info] {msg}")


def ok(msg: str) -> None:
    console.print(f"[ok]✓[/ok] {msg}")


def step(msg: str) -> None:
    console.print(f"[step]→[/step] {msg}")


def warn(msg: str) -> None:
    console.print(f"[warn]⚠[/warn] {msg}")


def fail(msg: str) -> None:
    err_console.print(f"[err]✗[/err] {msg}")


def human_bytes(n: float) -> str:
    """1.5 GB, 900 MB, 12 kB -- three characters of number, no more.

    Decimal units, not binary: this is reported next to disk-space figures the
    operating system also quotes decimally, and the difference between GiB and
    GB in a progress readout is noise nobody acts on.
    """
    for unit in ("B", "kB", "MB", "GB", "TB"):
        if abs(n) < 1000 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1000.0
    return f"{n:.1f} TB"


def human_duration(seconds: float | None) -> str:
    """A duration a human reads at a glance: 12s, 4m30s, 2h05m."""
    if seconds is None or seconds != seconds or seconds < 0:
        return "--"
    seconds = int(seconds)
    if seconds < 60:
        return f"{seconds}s"
    if seconds < 3600:
        return f"{seconds // 60}m{seconds % 60:02d}s"
    return f"{seconds // 3600}h{(seconds % 3600) // 60:02d}m"
