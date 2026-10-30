"""What gpxil files take on disk, set against the files they stand for (the whole genome).

The reference, GRCh38 (Ensembl release 116, the soft-masked primary assembly):
- the gpxil, one file (GRCh38.gpxil.tif): its ref plane's stored bytes, and the rest (IFDs, tags, metadata);
- the Gpxil browser's per-chromosome gpxil files, together;
- FASTA, uncompressed (60 bases a line: computed); Ensembl's FASTA.gz as distributed (gpxiltools' benchmark records
  its size); the bgzipped FASTA access.py wrote, with its .fai and .gzi; and UCSC's .2bit, its exact size worked out
  from the sequence (two bits a base, and blocks for the runs of N and of soft-masked bases). A .2bit keeps no IUPAC
  code but N: the others become N, and how many there are is counted.

The sample, HG002 (GIAB's v5.0q small variants):
- the VCF.gz as distributed, with its .tbi, and the benchmark regions' BED, which says where the VCF's silence means
  the reference and where it means nothing is known;
- the same records as BCF (bcftools view -Ob) with a .csi index;
- the sample-only gpxil (roundtrip.py's), which holds both: each plane's stored bytes, the sidecar's, and the rest.

    python3 sizes.py REF_GPXIL BROWSER_DIR SAMPLE_GPXIL VCF BED WORK   -> results/sizes.json

BROWSER_DIR: the Gpxil browser's site/data/grch38/gpxil. WORK: where access.py left GRCh38.fa.bgz, and the BCF is
written (and removed).
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter

import gpxiltools
import numpy as np
from gpxiltools import codes as C
from gpxiltools.container import SIDECAR_TAG, GpxilFile
from pysam import bcftools

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
ENSEMBL_FASTA_GZ = 942_484_168  #: Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz (gpxiltools' benchmarks/results)
ACGTN = [C.BASE_A, C.BASE_C, C.BASE_G, C.BASE_T, C.BASE_N]


def stored(path: str) -> dict:
    """A gpxil file's bytes: each plane's strips, the sidecars, and the rest."""
    per = Counter()
    with GpxilFile(path) as g:
        tf = g._tf
        for e in g.metadata["contigs"]:
            for plane, idx in (e.get("pages") or {}).items():
                page = tf.pages[idx]
                per[plane] += int(sum(page.databytecounts))
                if plane == "gt":
                    tag = page.tags.get(SIDECAR_TAG)
                    if tag is not None:
                        per["sidecar"] += len(bytes(tag.value))
        n = len(g.metadata["contigs"])
    total = os.path.getsize(path)
    return {"bytes": total, "contigs": n, "planes": dict(per), "other": total - sum(per.values())}


def runs(mask: np.ndarray) -> int:
    d = np.diff(mask.view(np.int8), prepend=np.int8(0), append=np.int8(0))
    return int(np.count_nonzero(d == 1))


def reference_baselines(ref: str) -> dict:
    """The exact size of a .2bit and of an uncompressed FASTA of the reference gpxil's sequences."""
    names = []
    twobit = 16  # signature, version, sequence count, reserved
    fasta = 0
    bases = n_bases = masked = iupac = 0
    with GpxilFile(ref) as g:
        for name in g.contigs:
            codes = g[name].read("ref")
            n = codes.shape[0]
            low = np.bitwise_and(codes, 0x0F)
            not_acgt = ~np.isin(low, [C.BASE_A, C.BASE_C, C.BASE_G, C.BASE_T])
            nb = runs(not_acgt)  # (a .2bit writes every non-ACGT base as N)
            mb = runs((codes & C.SOFTMASK) != 0)
            twobit += 1 + len(name) + 4  # index entry: name length, name, offset
            twobit += 4 + 4 + 8 * nb + 4 + 8 * mb + 4 + (n + 3) // 4
            fasta += len(f">{name}\n") + n + -(-n // 60)
            bases += n
            n_bases += int(np.count_nonzero(low == C.BASE_N))
            iupac += int(np.count_nonzero(~np.isin(low, ACGTN)))
            masked += int(np.count_nonzero(codes & C.SOFTMASK))
            names.append(name)
            del codes, low, not_acgt
    return {
        "sequences": len(names),
        "bases": bases,
        "N": n_bases,
        "other_iupac": iupac,
        "soft_masked": masked,
        "fasta_uncompressed_bytes": fasta,
        "twobit_bytes": twobit,
    }


def bcf(vcf: str, work: str) -> dict:
    out = os.path.join(work, "HG002_v5.0q.bcf")
    bcftools.view("-Ob", "-o", out, vcf, catch_stdout=False)
    bcftools.index(out, catch_stdout=False)
    res = {"bytes": os.path.getsize(out), "csi_bytes": os.path.getsize(out + ".csi")}
    os.unlink(out)
    os.unlink(out + ".csi")
    return res


def main() -> None:
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    ref, browser, sample, vcf, bed, work = sys.argv[1:]
    fa = os.path.join(work, "GRCh38.fa.bgz")
    files = sorted(f for f in os.listdir(browser) if f.endswith(".gpxil.tif"))
    res = {
        "gpxiltools": gpxiltools.__version__,
        "reference": {
            "gpxil": stored(ref),
            "browser_files": {
                "files": len(files),
                "bytes": sum(os.path.getsize(os.path.join(browser, f)) for f in files),
            },
            "ensembl_fasta_gz_bytes": ENSEMBL_FASTA_GZ,
            "bgzipped_fasta": {
                "bytes": os.path.getsize(fa),
                "fai_bytes": os.path.getsize(fa + ".fai"),
                "gzi_bytes": os.path.getsize(fa + ".gzi"),
            }
            if os.path.exists(fa + ".gzi")
            else None,
            **reference_baselines(ref),
        },
        "sample": {
            "vcf_gz_bytes": os.path.getsize(vcf),
            "tbi_bytes": os.path.getsize(vcf + ".tbi"),
            "bed_bytes": os.path.getsize(bed),
            "bcf": bcf(vcf, work),
            "gpxil": stored(sample),
        },
    }
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "sizes.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
