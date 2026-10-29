"""Shared code of experiment 05: the Platinum Pedigree (CEPH 1463) and where its files are.

The family (Porubsky et al. 2025; github.com/Platinum-Pedigree-Consortium): four grandparents (G1), their children
NA12877 and NA12878 (G2), eight G3 children of those two, and the G4 children of two of them, NA12879 (with 200080)
and NA12886 (with 200100). Five members are in controlled access only (NA12883, NA12884, NA12887, 200103, 200105),
and 200081 has no calls; the 22 others are in the open joint calls. Every open child with both parents open is a
trio: 16 of them.
"""

from __future__ import annotations

import gzip
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(HERE, "figures")
EXP04 = os.path.join(os.path.dirname(HERE), "04-ga-sequences")
sys.path.insert(0, EXP04)

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
BUCKET = "https://platinum-pedigree-data.s3.amazonaws.com"
HIFI = f"{BUCKET}/variants/small_variants/GRCh38/CEPH1463.GRCh38.deepvariant.glnexus.oa.vcf.gz"
ONT = f"{BUCKET}/variants/small_variants/GRCh38/CEPH1463.GRCh38.ont-clair3.glnexus.oa.vcf.gz"
ILLUMINA = f"{BUCKET}/variants/small_variants/GRCh38/CEPH1463.GRCh38.illumina-dragen.oa.vcf.gz"
HQ = f"{BUCKET}/variants/small_variant_truthset/GRCh38/hq_regions_final.bed.gz"


def dipcall(sample: str) -> str:
    """The assembly-based calls (dipcall, GRCh38) of a sample."""
    name = sample if sample.startswith("NA") else f"K1463_{sample}"
    return f"{BUCKET}/variants/assembly-based/dipcall/GRCh38/{name}.dip.vcf.gz"


# child: (father, mother), for the open trios
TRIOS = {
    "NA12877": ("NA12889", "NA12890"),
    "NA12878": ("NA12891", "NA12892"),
    **{c: ("NA12877", "NA12878") for c in ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")},
    **{c: ("200080", "NA12879") for c in ("200082", "200084", "200085", "200086", "200087")},
    **{c: ("NA12886", "200100") for c in ("200101", "200102", "200104", "200106")},
}
GENERATION = {"NA12877": "G2", "NA12878": "G2", "NA12879": "G3", "NA12881": "G3", "NA12882": "G3", "NA12885": "G3",
              "NA12886": "G3"}
GENERATION.update({c: "G4" for c in TRIOS if c not in GENERATION})


def children(p: str) -> list[str]:
    return [c for c, fm in TRIOS.items() if p in fm]


def descendants(p: str) -> set[str]:
    out = set()
    for c in children(p):
        out |= {c} | descendants(c)
    return out


def siblings(c: str) -> list[str]:
    return [s for s, fm in TRIOS.items() if s != c and fm == TRIOS.get(c)]


def read_runs(kinds=("GA", "TC")) -> dict:
    """Experiment 04's pure GA- and TC-sequences of GRCh38: {chrom: [(start, end, kind)]}, sorted."""
    runs: dict = {}
    with gzip.open(os.path.join(EXP04, "results", "ga-sequences.GRCh38.bed.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c, s, e, k = line.split("\t")[:4]
            if k in kinds:
                runs.setdefault(c, []).append((int(s), int(e), k))
    for v in runs.values():
        v.sort()
    return runs
