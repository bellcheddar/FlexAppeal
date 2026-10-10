#!/usr/bin/env python3
"""Generate the favicon set from the FlexAppeal mark.

Section 11 lists a favicon set as a deliverable and the site had none: no
`favicon.ico`, no `apple-touch-icon`, and no icon link in `base.html` at all, so
every tab showed the browser's default document glyph.

Source is `FlexAppeal_icon_300.png`, the 300 px raster of the brand mark. A
raster rather than the SVG because the SVG carries a large C2PA manifest and
there is no SVG rasteriser in this environment, while the 300 px PNG downsamples
cleanly to every size needed.

Re-runnable, and it reports each file's size so a silently empty output is
visible.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "FlexAppeal_icon_300.png"
OUT = ROOT / "flexappeal" / "static" / "icons"

# The conventional set. `favicon.ico` carries three sizes in one file because
# that is what a browser asks for at `/favicon.ico` when nothing else is
# declared, and some still do regardless of the link tags.
ICO_SIZES = [(16, 16), (32, 32), (48, 48)]
PNG_SIZES = {
    "favicon-16x16.png": 16,
    "favicon-32x32.png": 32,
    # 180 is what iOS asks for, and it is not padded here: the mark already has
    # its own margin, and adding more makes it small on a home screen.
    "apple-touch-icon.png": 180,
    "android-chrome-192x192.png": 192,
    "android-chrome-512x512.png": 512,
}

MANIFEST = """{
  "name": "FlexAppeal",
  "short_name": "FlexAppeal",
  "icons": [
    { "src": "/static/icons/android-chrome-192x192.png", "sizes": "192x192", "type": "image/png" },
    { "src": "/static/icons/android-chrome-512x512.png", "sizes": "512x512", "type": "image/png" }
  ],
  "theme_color": "#1e73be",
  "background_color": "#ffffff",
  "display": "browser"
}
"""


def main() -> int:
    from PIL import Image

    if not SOURCE.is_file():
        print(f"no source mark at {SOURCE.name}")
        return 1

    OUT.mkdir(parents=True, exist_ok=True)
    source = Image.open(SOURCE).convert("RGBA")
    print(f"source {SOURCE.name}: {source.width}x{source.height}")

    for name, size in PNG_SIZES.items():
        # LANCZOS both ways: 512 is an upscale from a 300 px source, and a
        # nearest-neighbour upscale on a mark with curves looks broken at the
        # size an Android splash screen uses it.
        resized = source.resize((size, size), Image.LANCZOS)
        path = OUT / name
        resized.save(path, "PNG", optimize=True)
        print(f"  {name:32s} {size:>4}px  {path.stat().st_size:>7,} bytes")

    ico = OUT / "favicon.ico"
    source.save(ico, "ICO", sizes=ICO_SIZES)
    print(f"  {'favicon.ico':32s} {'multi':>4}    {ico.stat().st_size:>7,} bytes"
          f"  ({', '.join(f'{w}x{h}' for w, h in ICO_SIZES)})")

    manifest = OUT / "site.webmanifest"
    manifest.write_text(MANIFEST)
    print(f"  {'site.webmanifest':32s} {'':>4}    {manifest.stat().st_size:>7,} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
