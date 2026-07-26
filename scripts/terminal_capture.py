#!/usr/bin/env python3
"""Turn a captured terminal session into embeddable HTML for the Example tab.

The input is what `script(1)` recorded while a bundle ran: real bytes from a
real run, escape sequences and all. The output is HTML, produced by rich's own
exporter, so the colours and glyphs on the page are the ones that were on the
terminal rather than an artist's impression of them.

Why not a PNG screen grab. A grab needs a windowed terminal and a human to
press the button: it cannot run unattended and cannot be regenerated when the
run changes. This reads the log the run already produced, so the page can never
show a screenshot of an older version of the tool.

Why HTML and not rich's SVG export, which was the first attempt. The SVG
exporter emits a positioned element per styled run of characters, which came to
2 MB for a twenty-six line frame -- ten megabytes across the set, for something
that is fundamentally coloured text. The HTML export of the same frames is
around a fortieth of that, stays selectable and searchable in the browser, and
is read out sensibly by a screen reader. Both are equally faithful: the colours
come from the recorded bytes either way.

Class-based CSS, not inline styles, for the same reason -- a style attribute on
every character cell was itself a megabyte a frame. rich numbers its classes
.r1, .r2, ... from zero in each export, so several captures on one page would
collide; each fragment's classes are prefixed with its own name.

The interesting frames are the ones mid-progress. A progress bar redraws over
itself, so a naive capture of the final state shows every bar at 100% and none
of the live readouts -- exactly the thing worth showing. `--at` picks frames by
the fraction of the way through each stage instead.

    python scripts/terminal_capture.py Examples/lysozyme_10ns/terminal.log \\
        --analysis-log Examples/lysozyme_10ns/analyse.log \\
        --out Examples/lysozyme_10ns/screenshots

Both halves of what a user does are captured, because both are what they will
see: run.py produces the trajectory and analyse.py turns it into the results
file. Names carry the phase and then a number -- run-01-startup,
analyse-02-packing -- so sorting within a phase gives the order things
happened, and the page can group them without inferring anything from the
numbering.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from rich.console import Console
from rich.text import Text

# `script` writes the pty's raw output: carriage returns separate the redraws of
# a live region, and these two sequences are terminal control rather than
# styling, which rich's ANSI decoder does not consume.
_CURSOR = re.compile(r"\x1b\[\?25[lh]")
_ERASE = re.compile(r"\x1b\[\d*[JK]")


_SPINNER = "\u280b\u2819\u2839\u2838\u283c\u2834\u2826\u2827\u2807\u280f"


def frames(log: str) -> list[str]:
    """Every state the terminal displayed, one per entry, in order.

    Both carriage returns and newlines end a frame. Under a pty, `script`
    records a live region's redraws as separate lines anyway, so the two
    separators are equivalent here -- an earlier version of this split on them
    in sequence, treating one as "lines" and the other as "redraws within a
    line", and matched nothing because every line turned out to hold exactly
    one redraw.
    """
    normalised = log.replace("\r\n", "\n").replace("\r", "\n")
    out = []
    for raw in normalised.split("\n"):
        cleaned = _ERASE.sub("", _CURSOR.sub("", raw))
        # script(1) echoes the terminating EOF as "^D" and leaves the odd
        # backspace behind. They belong to the recording, not the run.
        cleaned = cleaned.replace("^D", "").replace("\x08", "")
        if cleaned.strip():
            out.append(cleaned)
    return out


def _plain(s: str) -> str:
    return re.sub(r"\x1b\[[0-9;]*[A-Za-z]", "", s).strip()


def shorten_home(text: str) -> str:
    """Replace an absolute home directory with ~.

    The run summary prints the output directory, which is somebody's home
    directory, and these captures go on a public page. This is the only edit
    made to the recorded bytes anywhere in this script, and the page says so --
    an undisclosed edit would make its claim that these are the real bytes
    false, which is worse than the path being there.

    Trimming instead was the first thought, but the path sits mid-table in the
    summary between the wall time and the peak memory rows, so cutting there
    would lose the memory figures that are the reason the summary is shown.
    """
    home = str(Path.home())
    return text.replace(home, "~") if home not in ("/", "") else text


def _is_live(frame: str) -> bool:
    """Is this frame one redraw of a bar or spinner, rather than settled text?"""
    text = _plain(frame)
    return bool(text) and (text[0] in _SPINNER or "\u2501" in text)


def settled(frames_: list[str]) -> list[str]:
    """The log with each live region collapsed to its final redraw.

    Reads as the user would remember it: the lines that stayed on screen, and
    one line per stage that was animating.
    """
    out = []
    for frame in frames_:
        if _is_live(frame) and out and _is_live(out[-1]):
            out[-1] = frame          # same live region, later state
        else:
            out.append(frame)
    return out


def stage_capture(frames_: list[str], pattern: str, at: float,
                  context: int = 3) -> list[str] | None:
    """Settled context, then one live line caught `at` of the way through.

    `at` indexes the matching redraws, so 0.5 is genuinely mid-stage. The
    finished bar would sit at 100% with its readouts gone, which is the one
    state not worth showing.
    """
    matching = [i for i, f in enumerate(frames_)
                if _is_live(f) and re.search(pattern, _plain(f))]
    if not matching:
        return None
    chosen = matching[min(int(len(matching) * at), len(matching) - 1)]
    before = [f for f in frames_[:matching[0]] if not _is_live(f)][-context:]
    return before + [frames_[chosen]]


_BODY = re.compile(r"<pre[^>]*>(.*)</pre>", re.S)


def render(lines: list[str], path: Path, title: str, width: int = 120) -> None:
    """Write one capture as an HTML fragment the page can drop straight in."""
    import io as _io

    console = Console(record=True, width=width, file=_io.StringIO())
    for line in lines:
        # from_ansi is what makes this a capture rather than a re-creation: the
        # styling comes out of the recorded bytes, not from restyling the text.
        console.print(Text.from_ansi(shorten_home(line)), overflow="fold")

    # Just the styled body. rich's export wraps it in a whole document with its
    # own <html> and background, which cannot be nested inside a page.
    full = console.export_html(inline_styles=False,
                               code_format="<style>{stylesheet}</style><pre>{code}</pre>")

    # Namespace this fragment's classes so several can share a page.
    # "t" first, deliberately: rich numbers its classes .r1, .r2 from zero in
    # every export, so a page with several captures needs them namespaced -- but
    # a CSS class may not begin with a digit, and these files are 01-, 02-, ...
    # ".01-startup-r1" is silently invalid and every colour is simply dropped.
    prefix = "t" + re.sub(r"[^a-z0-9]+", "-", path.stem.lower()).strip("-")
    full = re.sub(r"\.r(\d+)\b", rf".{prefix}-r\1", full)
    full = re.sub(r'class="r(\d+)"', rf'class="{prefix}-r\1"', full)

    style = re.search(r"<style>(.*?)</style>", full, re.S)
    body = _BODY.search(full)
    path.write_text(
        f"<figure class=\"md-term\">\n"
        f"  <figcaption>{title}</figcaption>\n"
        f"  <style>{style.group(1) if style else ''}</style>\n"
        f"  <pre class=\"md-term-body\">{body.group(1) if body else full}</pre>\n"
        f"</figure>\n")


def capture_analysis(log: Path, out: Path, width: int) -> None:
    """The second half: analyse.py turning the trajectory into a results file.

    Matched by position rather than by label. The run has named stages that can
    be searched for -- "minimising", "producing" -- but the analysis bar is
    relabelled with whichever metric is in flight, and those depend on what the
    user asked for. So the mid-analysis frame is simply the middle live frame,
    whatever it happens to be computing.
    """
    all_frames = frames(log.read_text(errors="replace"))
    if not all_frames:
        return
    flat = settled(all_frames)

    # r"\S" matches any non-blank frame; stage_capture has already narrowed to
    # the live ones, so this is "halfway through the bar" and nothing more.
    middle = stage_capture(all_frames, r"\S", 0.5)
    if middle:
        render(middle, out / "analyse-01-metrics.html",
               "Analysing: one bar across every metric, relabelled as it goes", width)

    # Ends at the "results file(s) ready" tick rather than running to the end of
    # the log. The two lines after it echo the results file's absolute path,
    # which is somebody's home directory and has no business on a public page.
    # Trimmed rather than rewritten: the section claims these are the bytes the
    # scripts wrote, and editing them would make that claim false, whereas
    # showing a shorter window of them does not.
    for i, line in enumerate(flat):
        if "packing" in _plain(line):
            tail = flat[max(0, i - 1):]
            for j, later in enumerate(tail):
                if "results file" in _plain(later):
                    tail = tail[: j + 1]
                    break
            render(tail, out / "analyse-02-packing.html",
                   "Packing the results file the Analysis tab reads", width)
            break


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("log", type=Path, help="the run capture (script(1) output)")
    ap.add_argument("--analysis-log", type=Path, default=None,
                    help="the analyse.py capture, if there is one")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--width", type=int, default=120)
    args = ap.parse_args()

    if not args.log.is_file():
        print(f"no such log: {args.log}", file=sys.stderr)
        return 1
    args.out.mkdir(parents=True, exist_ok=True)

    all_frames = frames(args.log.read_text(errors="replace"))
    if not all_frames:
        print("the log contained no output", file=sys.stderr)
        return 1
    flat = settled(all_frames)

    # Startup: everything up to and including the platform decision.
    head = []
    for line in flat:
        head.append(line)
        if "using" in _plain(line) and "ns/day" in _plain(line):
            break
    render(head, args.out / "run-01-startup.html",
           "Startup, and the platform benchmark that picks OpenCL", args.width)

    for name, pattern, at, title in [
        ("02-minimising", r"minimising", 0.6,
         "Energy minimisation: a spinner, not a bar -- OpenMM's iteration count is not monotonic"),
        ("03-heating", r"heating", 0.6, "Heating, 50 K to 310 K"),
        ("04-equilibrating", r"equilibrat", 0.5,
         "Equilibration, with the positional restraints released in stages"),
        ("05-producing", r"producing", 0.5,
         "Production: ns/day, resident memory, free memory and swap, all live"),
    ]:
        capture = stage_capture(all_frames, pattern, at)
        if capture is None:
            print(f"  (nothing matched {pattern!r}; skipping {name})")
            continue
        render(capture, args.out / f"run-{name}.html", title, args.width)

    # The closing summary.
    for i, line in enumerate(flat):
        if "all runs complete" in _plain(line):
            render(flat[max(0, i - 2):], args.out / "run-06-summary.html",
                   "The closing summary", args.width)
            break

    if args.analysis_log and args.analysis_log.is_file():
        capture_analysis(args.analysis_log, args.out, args.width)

    written = sorted(p.name for p in args.out.glob("*.html"))
    print(f"wrote {len(written)} screenshots to {args.out}:")
    for w in written:
        print(f"  {w}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
