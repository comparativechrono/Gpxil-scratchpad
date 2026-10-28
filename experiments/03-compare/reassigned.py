"""Where a VCF encoded by gpxil 0.1's rules and by 0.2's differ: the bases, the insertions that moved to the other
haplotype, and the VCF records around each, to see that 0.2 placed them as the records mean.

    python3 -I reassigned.py OLD.gpxil.tif NEW.gpxil.tif VCF.gz [CONTIG]

(needs pysam and gpxiltools; results-0.2/reassigned.txt has HG002's, HG003's and HG004's chromosome 17)
"""

from __future__ import annotations

import sys

import numpy as np
import pysam
from gpxiltools import codes as C
from gpxiltools.container import GpxilFile


def genotypes(gt: np.ndarray, lo: int, hi: int) -> str:
    return " ".join(C.describe_gt_byte(int(v)) for v in gt[lo:hi])


def records(vcf, sample: str, contig: str, lo: int, hi: int) -> list[str]:
    out = []
    for rec in vcf.fetch(contig, max(0, lo), hi):
        s = rec.samples[sample]
        gt = ("|" if s.phased else "/").join("." if g is None else str(g) for g in s["GT"])
        out.append(f"  {rec.pos:,} {rec.ref}>{','.join(rec.alts)} {gt}")
    return out


def main(old_path: str, new_path: str, vcf_path: str, contig: str = "chr17") -> None:
    vcf = pysam.VariantFile(vcf_path)
    sample = list(vcf.header.samples)[0]
    with GpxilFile(old_path) as a, GpxilFile(new_path) as b:
        ca, cb = a[contig], b[contig]
        gta, gtb = ca.read("gt"), cb.read("gt")
        diff = np.flatnonzero(gta != gtb)
        flags = int(np.count_nonzero(ca.read("flags") != cb.read("flags")))
        sa = {(e.offset, e.hap, e.kind, e.data) for e in ca.sidecar()}
        sb = {(e.offset, e.hap, e.kind, e.data) for e in cb.sidecar()}
        moved = sorted({x[0] for x in sa ^ sb})
        enc = cb.entry.get("encoding", {})
        print(f"{contig}: {diff.size} gt values and {flags} flags differ; {len(moved)} insertions moved")
        print(f"overlapping_calls_reassigned {enc.get('overlapping_calls_reassigned')}, conflicts {enc.get('conflicts')}")
        for run in np.split(diff, np.flatnonzero(np.diff(diff) > 30) + 1) if diff.size else []:
            lo, hi = int(run[0]), int(run[-1]) + 1
            print(f"\nbases {lo + 1:,}-{hi:,}")
            print("  0.1:", genotypes(gta, lo - 1, hi + 1))
            print("  0.2:", genotypes(gtb, lo - 1, hi + 1))
            print("\n".join(records(vcf, sample, contig, lo - 30, hi + 30)))
        for o in moved:
            before = sorted(x[1:] for x in sa if x[0] == o)
            after = sorted(x[1:] for x in sb if x[0] == o)
            print(f"\ninsertion after {o + 1:,}: 0.1 {before}, 0.2 {after}")
            print("  0.1:", genotypes(gta, o - 3, o + 4))
            print("\n".join(records(vcf, sample, contig, o - 30, o + 30)))


if __name__ == "__main__":
    main(*sys.argv[1:])
