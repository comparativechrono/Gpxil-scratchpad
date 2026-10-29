"""Albrecht-Buehler's GA-complex image (PLoS ONE 2009, Figs 1, 3, 4), redrawn on GRCh38: one row per pure GA-sequence
of a chromosome, its 400 upstream bases then its first 200, aligned at its start, in his greys (A black, C dark grey,
T light grey, G white: gpxil's gpxi palette). Rows are sorted by the run of A that starts the GA-sequence (its poly(A)
segment; a run of A just before G/A bases is part of the GA-sequence), longest at the top, as in his Fig. 4.

    python3 stack.py REF_DIR [CHROM]       -> figures/stack_CHROM.png (one pixel a base), figures/stack_CHROM_labelled.png
"""

import os
import sys

import numpy as np
from PIL import Image, ImageDraw

import ga
from gpxiltools import codes as C

UP, IN = 400, 200
GREY = {C.BASE_A: 0, C.BASE_C: 85, C.BASE_G: 255, C.BASE_T: 170}


def poly_a(seq: np.ndarray, at: int, tolerance: int = 2) -> int:
    """Length of the A-rich run from ``at``: A's, allowing up to ``tolerance`` other bases, stopping at the next;
    trailing non-A bases don't count."""
    n, other = 0, 0
    for x in seq[at:].tolist():
        if x != C.BASE_A:
            other += 1
            if other > tolerance:
                break
        n += 1
    while n and seq[at + n - 1] != C.BASE_A:
        n -= 1
    return n


def rows_of(base: np.ndarray, runs, kind: str = "GA"):
    """The stack's rows (UP bases before the GA-sequence's start, IN from it, on its purine strand) and poly(A)s."""
    rows, pa = [], []
    for s, e, k in runs:
        if k != kind:
            continue
        if kind == "GA":
            if s < UP or s + IN > base.size:
                continue
            row = base[s - UP : s + IN]
        else:
            if e + UP > base.size or e - IN < 0:
                continue
            row = ga.COMP[base[e - IN : e + UP][::-1]]
        rows.append(row)
        pa.append(poly_a(row, UP))
    return np.array(rows, dtype=np.uint8), np.array(pa)


def main(ref_dir: str, ch: str = "chr1") -> None:
    runs = ga.read_runs()[ch]
    with ga.reference(ref_dir, ch) as g:
        base = g[ch].read("ref") & 0x0F
    rows, pa = rows_of(base, runs)
    order = np.argsort(-pa, kind="stable")
    lut = np.full(16, 230, dtype=np.uint8)
    for k, v in GREY.items():
        lut[k] = v
    img = lut[rows[order]]
    os.makedirs(ga.FIGURES, exist_ok=True)
    Image.fromarray(img, mode="L").save(os.path.join(ga.FIGURES, f"stack_{ch}.png"))
    # a labelled view, 880 rows high
    h = 880
    view = Image.fromarray(img, mode="L").convert("RGB").resize((img.shape[1], h), Image.NEAREST)
    W = img.shape[1] + 230
    out = Image.new("RGB", (W, h + 70), (255, 255, 255))
    d = ImageDraw.Draw(out)
    out.paste(view, (10, 50))
    d.text((10, 6), f"GRCh38 {ch}: its {len(rows):,} pure GA-sequences, one per row: 400 bases upstream, then their first 200,", fill=(0, 0, 0))
    d.text((10, 20), "aligned at the GA-sequence's start. A black, C dark grey, T light grey, G white (gpxi)", fill=(0, 0, 0))
    d.line([(10 + UP, 44), (10 + UP, 48)], fill=(200, 0, 0), width=2)
    d.text((10 + UP - 70, 32), "GA-sequence starts", fill=(200, 0, 0))
    d.text((10, 50 + h + 4), "<- 400 bases upstream", fill=(0, 0, 0))
    k, n = 0, len(rows)
    for label, f in ((">= 15", lambda a: a >= 15), ("5-14", lambda a: (a >= 5) & (a < 15)), ("< 5", lambda a: a < 5)):
        m = int(f(pa).sum())
        y1, y2 = 50 + int(h * k / n) + 2, 50 + int(h * (k + m) / n) - 3
        x = img.shape[1] + 20
        d.line([(x, y1), (x, y2)], fill=(0, 0, 0), width=3)
        d.text((x + 8, (y1 + y2) // 2 - 12), f"poly(A) {label} bases", fill=(0, 0, 0))
        d.text((x + 8, (y1 + y2) // 2 + 2), f"{m:,} rows", fill=(80, 80, 80))
        k += m
    d.text((img.shape[1] + 28, 60), "Alu-like stripes upstream", fill=(200, 0, 0))
    out.save(os.path.join(ga.FIGURES, f"stack_{ch}_labelled.png"))
    print(f"{ch}: {len(rows):,} rows; poly(A) >= 15: {(pa >= 15).sum():,}, 5-14: {((pa >= 5) & (pa < 15)).sum():,}, < 5: {(pa < 5).sum():,}")


if __name__ == "__main__":
    if len(sys.argv) not in (2, 3):
        sys.exit(__doc__)
    main(*sys.argv[1:])
