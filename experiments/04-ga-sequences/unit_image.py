"""One pure GA-sequence as a genome pixel image wrapped at its repeat period: GRCh38, then HG002's two haplotypes, each
unit of PERIOD bases a column, read left to right; HG002's substitutions outlined in red, deleted bases pale yellow.
Greys: A black, C dark grey, G white, T light grey (gpxi).

    python3 unit_image.py REF_DIR GPXIL_FILE CHROM START END PERIOD OUT.png      (START, END: 1-based, inclusive)
"""

import sys

import numpy as np
from PIL import Image, ImageDraw

from gpxiltools import GpxilFile
from gpxiltools import codes as C

GREY = {C.BASE_A: (0, 0, 0), C.BASE_C: (85, 85, 85), C.BASE_G: (255, 255, 255), C.BASE_T: (170, 170, 170)}
DEL, UNK = (255, 240, 170), (170, 200, 235)


def main(ref_dir, path, ch, start, end, period, out) -> None:
    s, e, period = int(start) - 1, int(end), int(period)
    with GpxilFile(f"{ref_dir.rstrip('/')}/{ch}.gpxil.tif") as gr, GpxilFile(path) as g:
        ref = gr[ch].read("ref", s, e) & 0x0F
        gt = g[ch].read("gt", s, e)
    haps = [C.normalize_states(C.hap1(gt)), C.normalize_states(C.hap2(gt))]
    px, gap, top, margin = 12, 10, 16, 6
    n = e - s
    cols = (n + period - 1) // period
    panels = [("GRCh38", None), ("HG002 haplotype 1", haps[0]), ("HG002 haplotype 2", haps[1])]
    img = Image.new("RGB", (cols * px + 2 * margin, len(panels) * (period * px + gap + top) + 20), (255, 255, 255))
    d = ImageDraw.Draw(img)
    y0 = 0
    for label, h in panels:
        d.text((margin, y0 + 2), label, fill=(0, 0, 0))
        for i in range(n):
            x, y = margin + (i // period) * px, y0 + top + (i % period) * px
            st = 0 if h is None else int(h[i])
            col = GREY.get(int(ref[i]), UNK) if st == C.HAP_REF else GREY[st] if st in C.HAP_SUB_CODES else DEL if st == C.HAP_DEL else UNK
            d.rectangle([x, y, x + px - 1, y + px - 1], fill=col, outline=(200, 200, 200))
            if h is not None and st in C.HAP_SUB_CODES:
                d.rectangle([x, y, x + px - 1, y + px - 1], outline=(220, 0, 0), width=2)
        y0 += period * px + gap + top
    d.text((margin, y0), f"{ch}:{int(start):,}-{int(end):,}, {period} bases a column", fill=(80, 80, 80))
    img.save(out)
    print(out, img.size, "substitutions:", int(sum(np.isin(h, C.HAP_SUB_CODES).sum() for h in haps)))


if __name__ == "__main__":
    if len(sys.argv) != 8:
        sys.exit(__doc__)
    main(*sys.argv[1:])
