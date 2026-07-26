#!/usr/bin/env python3
"""Key a flat background out of a video and emit it with a real alpha channel.

Written for the lysozyme animation on the Example tab, whose source is H.264
yuv420p -- a pixel format with no alpha at all -- on a flat warm off-white. The
page needs the site background to show through, so the transparency has to be
created rather than preserved.

Why not a colour key. The cartoon's coil regions are white fills inside black
outlines, and measured against this source they sit at (249.5, 248.7, 245.2)
while the background is (250, 249, 245). Nothing distinguishes them by colour;
a key on either punches holes straight through the protein. What separates them
is connectivity: the background reaches the frame border and the loop fills do
not, because the outline encloses them. So the mask is a flood fill from the
edges, and the enclosed whites stay opaque.

Edges are un-premultiplied rather than simply feathered. An antialiased pixel is
a blend of the subject and the old background, so carrying it over unchanged
paints a warm halo on a page that is not warm off-white. Given the known
background B and the computed alpha a, the true foreground is
(observed - (1-a)B) / a, which composites cleanly over anything.

Output is an animated WebP, which is the only web format this machine can
actually produce with an alpha channel. That is not the obvious choice and it
was not the first one; the alternatives were measured:

    libvpx-vp9  -pix_fmt yuva420p   alpha silently dropped
    libvpx      -pix_fmt yuva420p   alpha silently dropped
    libaom-av1  -pix_fmt yuva420p   alpha silently dropped
    hevc_videotoolbox  -pix_fmt bgra   alpha preserved, but Safari renders it
                                       and Chrome does not
    Pillow's WebP writer             this build has no animation support

ffmpeg 8.1.2 lists yuva420p among libvpx-vp9's supported formats and then
writes yuv420p regardless -- verified by round-tripping a half-transparent test
frame, not by reading the container metadata, which reports the colour stream's
format either way and looks the same whether the alpha survived or not.

Animated WebP is supported by Chrome 32+, Firefox 65+, Safari 14+ and Edge 18+,
carries a real alpha channel, and goes in an <img>, so it loops on its own with
no video element, no JavaScript and no controls to suppress.

On size: frame rate is the only lever that matters. Measured on this clip at
440px, dropping quality from 40 to 20 saved 11% while halving the frame rate
saved 47%, because WebP differences consecutive frames but has no motion
compensation, and the sharp outlines and alpha edges change everywhere as the
molecule turns. Scaling down after keying makes it *larger*, counter-intuitively
-- resampling smears the flat zero-alpha background into a field of slightly
different values, and the differencing has nothing left to exploit. Scale before
the key, which is why the fps and scale filters are on the decoder.

Usage:
    python scripts/key_video.py in.mp4 --out-dir examples/lysozyme_10ns/media \\
        --name lysozyme --width 500
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage

# How far a pixel may sit from the background colour and still be considered
# part of the background region for the flood fill. Generous: the outlines are
# near-black, hundreds of units away, so there is no risk of leaking through
# them into the enclosed fills.
CONNECT_TOLERANCE = 40

# Below this distance from the background colour a pixel is background, full
# stop. The source is lossy H.264, so its "flat" background is not flat: block
# noise moves pixels a few units either way, and without a dead zone those come
# out slightly opaque and print faint rectangular ghosts of the codec's macro
# blocks across the page. Visible immediately when composited over a saturated
# colour, and easy to miss over a near-white one.
NOISE_FLOOR = 12

# Alpha ramps from transparent to opaque between NOISE_FLOOR and here.
# Antialiased edge pixels land in this band and come out partly transparent.
EDGE_SOFTNESS = 30


def background_colour(frame: np.ndarray) -> np.ndarray:
    """The most common colour around the frame's border."""
    border = np.concatenate([
        frame[0, :], frame[-1, :], frame[:, 0], frame[:, -1],
    ])
    colours, counts = np.unique(border.reshape(-1, 3), axis=0, return_counts=True)
    return colours[counts.argmax()].astype(int)


def key_frame(frame: np.ndarray, bg: np.ndarray) -> np.ndarray:
    """One RGB frame in, one RGBA frame out."""
    distance = np.abs(frame.astype(int) - bg).max(axis=2)

    # Connectivity, not colour: only the region that reaches the border.
    labels, _ = ndimage.label(distance <= CONNECT_TOLERANCE)
    edge_labels = np.unique(np.concatenate([
        labels[0, :], labels[-1, :], labels[:, 0], labels[:, -1],
    ]))
    edge_labels = edge_labels[edge_labels != 0]
    outside = np.isin(labels, edge_labels)

    alpha = np.ones(distance.shape, dtype=np.float32)
    ramp = np.clip((distance - NOISE_FLOOR) / (EDGE_SOFTNESS - NOISE_FLOOR), 0.0, 1.0)
    alpha[outside] = ramp[outside]

    # Un-premultiply against the old background so edges carry no halo.
    a = alpha[..., None]
    rgb = np.where(
        a > 0.004,
        (frame.astype(np.float32) - (1.0 - a) * bg.astype(np.float32)) / np.maximum(a, 0.004),
        frame.astype(np.float32),
    )
    rgb = np.clip(rgb, 0, 255)

    return np.dstack([rgb, alpha * 255.0]).astype(np.uint8)


def probe(path: Path) -> tuple[int, int, str]:
    out = subprocess.run(
        ["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries",
         "stream=width,height,r_frame_rate", "-of", "csv=p=0:s=x", str(path)],
        capture_output=True, text=True, check=True).stdout.strip()
    w, h, rate = out.split("x")
    return int(w), int(h), rate


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--name", default="clip")
    ap.add_argument("--width", type=int, default=500, help="output width; height follows")
    ap.add_argument("--fps", type=float, default=8.0,
                    help="output frame rate. This is the only real size lever: WebP has "
                         "frame differencing but no motion compensation, so bytes scale "
                         "with frame count while quality barely moves them")
    ap.add_argument("--quality", type=int, default=40, help="WebP quality, 0-100")
    ap.add_argument("--poster-at", type=float, default=2.0, help="seconds for the poster frame")
    ap.add_argument("--keep-frames", type=Path, default=None,
                    help="write the keyed PNGs here and leave them, to retune the encode")
    args = ap.parse_args()

    if not args.source.is_file():
        print(f"no such file: {args.source}", file=sys.stderr)
        return 1
    args.out_dir.mkdir(parents=True, exist_ok=True)

    width, height, rate = probe(args.source)
    out_w = args.width - (args.width % 2)
    out_h = int(round(height * out_w / width))
    out_h -= out_h % 2
    print(f"  source {width}x{height} @ {rate}  ->  {out_w}x{out_h}")

    decode = subprocess.Popen(
        ["ffmpeg", "-v", "error", "-i", str(args.source),
         "-vf", f"fps={args.fps},scale={out_w}:{out_h}:flags=lanczos",
         "-f", "rawvideo", "-pix_fmt", "rgb24", "-"],
        stdout=subprocess.PIPE)

    frames_dir = Path(args.keep_frames) if args.keep_frames else Path(
        tempfile.mkdtemp(prefix="keyed-"))
    frames_dir.mkdir(parents=True, exist_ok=True)
    frame_bytes = out_w * out_h * 3
    bg = None
    poster_index = int(args.poster_at * args.fps)
    poster = None
    count = 0
    while True:
        raw = decode.stdout.read(frame_bytes)
        if len(raw) < frame_bytes:
            break
        frame = np.frombuffer(raw, np.uint8).reshape(out_h, out_w, 3)
        if bg is None:
            bg = background_colour(frame)
            print(f"  background keyed: rgb{tuple(int(c) for c in bg)}")
        rgba = key_frame(frame, bg)
        if count == poster_index:
            poster = rgba.copy()
        Image.fromarray(rgba, "RGBA").save(frames_dir / f"{count:05d}.png")
        count += 1
        if count % 200 == 0:
            print(f"  {count} frames", flush=True)

    decode.wait()
    print(f"  {count} frames keyed")

    if poster is not None:
        Image.fromarray(poster, "RGBA").save(args.out_dir / f"{args.name}-poster.png")

    webp = args.out_dir / f"{args.name}.webp"
    frame_ms = int(round(1000 / args.fps))
    subprocess.run(
        # -lossy is not the default: img2webp encodes every frame losslessly
        # unless told otherwise, which for 758 frames of flat artwork came to
        # 48 MB. -mixed lets it fall back to lossless per frame where that is
        # actually smaller, which for large flat areas it sometimes is.
        ["img2webp", "-loop", "0", "-d", str(frame_ms),
         "-lossy", "-q", str(args.quality), "-m", "6",
         *[str(p) for p in sorted(frames_dir.glob("*.png"))],
         "-o", str(webp)],
        check=True, capture_output=True)
    if not args.keep_frames:
        shutil.rmtree(frames_dir, ignore_errors=True)

    for path in (webp, args.out_dir / f"{args.name}-poster.png"):
        if path.is_file():
            print(f"  {path.name:24s} {path.stat().st_size / 1e6:6.2f} MB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
