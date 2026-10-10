#!/usr/bin/env python3
"""Raise every ANSI colour in the committed terminal clips that fails WCAG AA.

The clips are `rich` HTML exports of real run output, and several of the
palette's colours fall under the 4.5 contrast AA asks for against the terminal
background `#14161f`:

    #3a3a3a  grey    1.59
    #800080  purple  1.91
    #008000  green   3.51
    #008080  teal    3.78
    #0087af  cyan    4.37

The first version of this script fixed only the two axe happened to name in its
first report, which was the wrong shape of fix: axe samples nodes, so patching
what it reported left three colours that it had not got to yet. Every colour in
every clip is now checked against the threshold and the script fails if one is
missed, so a new recording cannot quietly introduce another.

Only the colour changes; the text of each recording is untouched, so what the
clip says about the run is exactly what it said before. Replacements keep the
hue and clear the threshold with margin, so a re-record at the same palette will
not drift back under it.

Re-runnable: an already-corrected file is left alone, and the script reports
what it changed rather than working silently.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CLIPS = ROOT / "examples" / "lysozyme_benzene_100ns" / "screenshots"

# The clip background, from the exported stylesheet.
BACKGROUND = "#14161f"
AA_NORMAL_TEXT = 4.5

# Each replacement holds the hue and clears the threshold with room: the
# computed ratios are 5.72, 5.77, 6.00, 5.61 and 6.57.
REPLACEMENTS = {
    "#3a3a3a": "#919191",
    "#800080": "#cc6fcc",
    "#008000": "#2fab2f",
    "#008080": "#00a0a0",
    "#0087af": "#2aa8d4",
}


def _luminance(colour: str) -> float:
    digits = colour.lstrip("#")
    channels = [int(digits[i:i + 2], 16) / 255 for i in (0, 2, 4)]
    linear = [(v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4)
              for v in channels]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def contrast(one: str, two: str) -> float:
    first, second = _luminance(one), _luminance(two)
    high, low = max(first, second), min(first, second)
    return (high + 0.05) / (low + 0.05)


def main() -> int:
    if not CLIPS.is_dir():
        print(f"no clips at {CLIPS}")
        return 1

    clips = sorted(CLIPS.glob("*.html"))
    touched = 0
    for path in clips:
        text = original = path.read_text()
        for old, new in REPLACEMENTS.items():
            text = text.replace(old, new)
        if text != original:
            counts = {old: original.count(old) for old in REPLACEMENTS
                      if old in original}
            path.write_text(text)
            touched += 1
            detail = ", ".join(f"{old} x{n}" for old, n in counts.items())
            print(f"  {path.name}: {detail}")
    print(f"{touched} of {len(clips)} clips changed")

    # Every colour that survives, not just the ones that were replaced. This is
    # what the first version of this script lacked, and it is why three failing
    # colours stayed in the clips after it ran.
    remaining: dict[str, float] = {}
    for path in clips:
        for colour in re.findall(r"color:\s*(#[0-9a-fA-F]{6})", path.read_text()):
            ratio = contrast(colour, BACKGROUND)
            if ratio < AA_NORMAL_TEXT:
                remaining[colour.lower()] = ratio
    if remaining:
        print(f"\nstill under {AA_NORMAL_TEXT} against {BACKGROUND}:")
        for colour, ratio in sorted(remaining.items(), key=lambda kv: kv[1]):
            print(f"  {colour}  {ratio:.2f}")
        return 1
    print(f"every clip colour clears {AA_NORMAL_TEXT} against {BACKGROUND}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
