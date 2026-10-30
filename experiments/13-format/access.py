"""How fast a region comes out of a gpxil, on this computer and over the web.

Regions of 1 kb, 100 kb and 1 Mb, placed at random (seeded) on GRCh38's chromosomes 1-22 in proportion to their
lengths, are read:

- the reference: from the GRCh38 gpxil (Contig.sequence), and from a bgzipped FASTA of the same sequence with its .fai
  and .gzi indexes (pysam's faidx); the two are checked to give the same bases;
- HG002's calls: from its gpxil (the gt and flags planes and the insertions of both haplotypes), and the VCF's records
  with their genotypes (tabix, pysam);
- HG002's two haplotypes as sequences, from its gpxil beside the reference gpxil (haplotype_sequence).

The files are opened once and read region after region; each kind of read is summed up by the median and the 90th
percentile over its regions. Opening is timed apart. Then some regions are read from the Gpxil browser's per-chromosome
GRCh38 files on GitHub Pages with HTTP range requests, each from a newly opened file: how many requests, how many bytes,
and how long. The local files were in the operating system's cache.

    python3 access.py REF_GPXIL SAMPLE_GPXIL VCF WORK   -> results/access.json

REF_GPXIL: GRCh38.gpxil.tif. SAMPLE_GPXIL: roundtrip.py's HG002_v5.0q.gpxil.tif. VCF: HG002_GRCh38_v5.0q_smvar.vcf.gz
(with its .tbi). WORK: where the bgzipped FASTA is written (and left, for sizes.py).
"""

from __future__ import annotations

import json
import os
import platform
import sys
import time

import gpxiltools
import numpy as np
import pysam
from gpxiltools.container import GpxilFile
from gpxiltools.export import haplotype_sequence, iter_reference

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
CHROMS = [f"chr{i}" for i in range(1, 23)]
SIZES = (1_000, 100_000, 1_000_000)
LOCAL_REGIONS = 100  #: regions of each size read from the local files
REMOTE_REGIONS = 5  #: and from the web
WEB = "https://comparativechrono.github.io/Gpxil-browser/data/grch38/gpxil/{}.gpxil.tif"
FASTA = "GRCh38.fa.bgz"


def fasta_blocks(name: str, pieces, width: int = 60):
    """One FASTA record, ``width`` bases a line, in blocks of many lines (gpxiltools' write_fasta_stream writes a line
    at a time, which is slow through BGZF)."""
    yield f">{name}\n".encode()
    carry = b""
    for piece in pieces:
        buf = carry + piece
        cut = len(buf) - len(buf) % width
        if cut:
            a = np.frombuffer(buf, dtype=np.uint8, count=cut).reshape(-1, width)
            yield np.hstack([a, np.full((a.shape[0], 1), ord("\n"), dtype=np.uint8)]).tobytes()
        carry = buf[cut:]
    if carry:
        yield carry + b"\n"


def bgzipped_fasta(ref: str, work: str) -> str:
    """The reference gpxil's sequences as a bgzipped FASTA (60 bases a line), indexed with .fai and .gzi."""
    path = os.path.join(work, FASTA)
    if os.path.exists(path + ".gzi") and os.path.exists(path + ".fai"):
        return path
    t = time.perf_counter()
    with GpxilFile(ref) as g:
        out = pysam.BGZFile(path, "wb")
        for name in g.contigs:
            for block in fasta_blocks(name, iter_reference(g[name], chunk=60 << 17)):
                out.write(block)
        out.close()
    pysam.faidx(path)
    print(f"wrote {path} in {time.perf_counter() - t:.0f} s", flush=True)
    return path


def regions(ref: str, seed: int = 13) -> dict[int, list[tuple[str, int, int]]]:
    with GpxilFile(ref) as g:
        lengths = np.array([g[c].length for c in CHROMS], dtype=np.float64)
    rng = np.random.default_rng(seed)
    out = {}
    for size in SIZES:
        picked = rng.choice(len(CHROMS), size=LOCAL_REGIONS, p=lengths / lengths.sum())
        out[size] = [
            (CHROMS[i], int(s), int(s) + size) for i, s in ((i, rng.integers(0, lengths[i] - size)) for i in picked)
        ]
    return out


def summary(times: list[float]) -> dict:
    a = np.array(times) * 1000
    return {
        "median_ms": round(float(np.median(a)), 3),
        "p90_ms": round(float(np.percentile(a, 90)), 3),
        "n": len(times),
    }


def timed(fn, regs) -> list[float]:
    out = []
    for r in regs:
        t = time.perf_counter()
        fn(*r)
        out.append(time.perf_counter() - t)
    return out


def local(ref: str, sample: str, vcf: str, fasta: str, regs: dict) -> dict:
    res: dict = {"open_ms": {}}
    for name, opener in (
        ("reference gpxil", lambda: GpxilFile(ref)),
        ("bgzipped FASTA", lambda: pysam.FastaFile(fasta)),
        ("sample gpxil, with the reference gpxil", lambda: GpxilFile(sample, reference=ref)),
        ("VCF and its index", lambda: pysam.VariantFile(vcf)),
    ):
        ts = []
        for _ in range(5):
            t = time.perf_counter()
            f = opener()
            ts.append(time.perf_counter() - t)
            f.close()
        res["open_ms"][name] = round(float(np.median(ts)) * 1000, 2)

    with (
        GpxilFile(ref) as g,
        pysam.FastaFile(fasta) as fa,
        GpxilFile(sample, reference=ref) as s,
        pysam.VariantFile(vcf) as vf,
    ):
        # the two references give the same bases
        for c, a, b in regs[SIZES[0]][:20] + regs[SIZES[-1]][:3]:
            assert g[c].sequence(a, b) == fa.fetch(c, a, b), (c, a, b)
        # (each file read through once beforehand, so the first regions aren't the ones that find the pages)
        for c, a, b in regs[SIZES[0]]:
            g[c].sequence(a, b), fa.fetch(c, a, b), s[c].read("gt", a, b), list(vf.fetch(c, a, b))

        def calls_gpxil(c, a, b):
            x = s[c]
            return x.read("gt", a, b), x.read("flags", a, b), x.insertions(1, a, b), x.insertions(2, a, b)

        def calls_vcf(c, a, b):
            return [r.samples[0]["GT"] for r in vf.fetch(c, a, b)]

        def haplotypes(c, a, b):
            return haplotype_sequence(s[c], 1, a, b), haplotype_sequence(s[c], 2, a, b)

        kinds = {
            "reference: gpxil": lambda c, a, b: g[c].sequence(a, b),
            "reference: bgzipped FASTA (faidx)": lambda c, a, b: fa.fetch(c, a, b),
            "calls: gpxil planes and insertions": calls_gpxil,
            "calls: VCF records with genotypes (tabix)": calls_vcf,
            "haplotypes: gpxil, both as sequences": haplotypes,
        }
        res["reads"] = {}
        res["records_per_region"] = {}
        for size in SIZES:
            res["reads"][str(size)] = {k: summary(timed(fn, regs[size])) for k, fn in kinds.items()}
            res["records_per_region"][str(size)] = float(np.median([sum(1 for _ in vf.fetch(*r)) for r in regs[size]]))
    return res


def remote(regs: dict) -> dict:
    out = {}
    t = time.perf_counter()
    g = GpxilFile(WEB.format("chr1"))
    out["open"] = {"requests": g._fh.requests, "bytes": g._fh.fetched, "ms": round((time.perf_counter() - t) * 1000, 1)}
    g.close()
    for size in SIZES:
        rows = []
        for c, a, b in regs[size][:REMOTE_REGIONS]:
            t = time.perf_counter()
            with GpxilFile(WEB.format(c)) as g:
                seq = g[c].sequence(a, b)
                rows.append(
                    {
                        "region": f"{c}:{a + 1}-{b}",
                        "requests": g._fh.requests,
                        "bytes": g._fh.fetched,
                        "file_bytes": g._fh.size,
                        "ms": round((time.perf_counter() - t) * 1000, 1),
                        "bases": len(seq),
                    }
                )
        out[str(size)] = {
            "regions": rows,
            "median_bytes": float(np.median([r["bytes"] for r in rows])),
            "median_requests": float(np.median([r["requests"] for r in rows])),
            "median_ms": float(np.median([r["ms"] for r in rows])),
        }
    return out


def main() -> None:
    if len(sys.argv) != 5:
        raise SystemExit(__doc__)
    ref, sample, vcf, work = sys.argv[1:]
    fasta = bgzipped_fasta(ref, work)
    regs = regions(ref)
    res = {
        "gpxiltools": gpxiltools.__version__,
        "pysam": pysam.__version__,
        "python": platform.python_version(),
        "cpu": platform.machine(),
        "cores": os.cpu_count(),
        "local": local(ref, sample, vcf, fasta, regs),
    }
    try:
        res["web"] = remote(regs)
    except OSError as e:  # (no network: the local numbers stand alone)
        res["web"] = {"error": str(e)}
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "access.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
