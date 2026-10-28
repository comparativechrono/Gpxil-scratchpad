"""Every record that spec 0.2's rule for overlapping calls (section 6.5) puts on the other haplotype, on every
chromosome of a VCF, and why: on the haplotype its written order gives it, the first base its allele occupies is one
that an earlier call deletes ("deleted"), or that an earlier call's allele starts at ("same position"), or neither.

    python3 -I moves.py VCF BED REFERENCE_DIR > moves.tsv

REFERENCE_DIR holds the Gpxil browser's chr1.gpxil.tif ... chr22.gpxil.tif (the reference, one file a chromosome).
Needs gpxiltools 0.6 or later. results-0.2/moves.tsv has HG002's, HG003's and HG004's v4.2.1 benchmarks.
"""

from __future__ import annotations

import os
import sys

from gpxiltools import codes as C
from gpxiltools import vcf as V
from gpxiltools.container import GpxilFile


class Watched(V.PlaneBuilder):
    """gpxiltools' builder, noting each record its rule moves."""

    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.moves: list[tuple] = []
        self.starts: list[set[int]] = [set(), set()]

    def _place_overlapping(self, alleles, call, p0, span_end):
        out = super()._place_overlapping(alleles, call, p0, span_end)
        if out != alleles:
            why = []
            for h, a in enumerate(alleles):  # the allele that didn't fit where the written order put it
                c = self._occupies(a, call, p0)
                if c is None or self._fits(h, c):
                    continue
                if (int(self.w[h][c[0]]) & 0x0F) == C.HAP_DEL:
                    why.append("deleted")
                elif c[0] in self.starts[h]:
                    why.append("same position")
                else:
                    why.append("other")
            self.moves.append((call.pos, call.ref, ",".join(call.alts), call.gt, "+".join(why)))
        return out

    def _apply_allele(self, h, a, call, p0, m, span_end):
        out = super()._apply_allele(h, a, call, p0, m, span_end)
        c = self._occupies(a, call, p0)
        if c is not None:
            self.starts[h].add(c[0])
        return out


def main(vcf: str, bed: str, refdir: str) -> None:
    sample = V.read_header(vcf).samples[0]
    with V.TempSplit(vcf, 0) as split:
        for chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
            path = os.path.join(refdir, f"{chrom}.gpxil.tif")
            if chrom not in split.index or not os.path.exists(path):
                continue
            with GpxilFile(path) as g:
                ref = g[chrom].read("ref")
            b = Watched(ref, callable_mask=V.load_bed(bed, chrom, ref.shape[0]), contig=chrom)
            for call in split.calls(chrom):
                b.add(call)
            for m in b.moves:
                print(sample, chrom, *m, sep="\t", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:])
