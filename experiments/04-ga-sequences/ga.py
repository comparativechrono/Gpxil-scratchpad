"""Shared code of experiment 04: pure GA-sequences of GRCh38, the gpxil files they are read from, and the call sets
of HG002 that are tested in them.

A pure GA-sequence (Albrecht-Buehler 2008, 2009) is a run of more than 50 bases that are only G and A on the + strand;
a pure TC-sequence, the same of C and T, is a pure GA-sequence of the other strand. Both are read here on their purine
strand, so that a TC-sequence's bases, and a variant in it, are complemented and taken in the opposite direction.
"""

from __future__ import annotations

import gzip
import os

import numpy as np

from gpxiltools import GpxilFile
from gpxiltools import codes as C
from gpxiltools.vcf import load_bed

AUTOSOMES = [f"chr{i}" for i in range(1, 23)]
CHROMS = AUTOSOMES + ["chrX", "chrY"]
PURE = 51  #: more than 50 bases (Albrecht-Buehler's pure GA-sequences)
COMMON = 21  #: more than 20 bases (his "common" GA-sequences)
FLANK = 400  #: the flanks of a GA-complex, either side (Albrecht-Buehler 2009)
CHUNK = 1 << 24

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(HERE, "figures")
RUNS_BED = os.path.join(RESULTS, "ga-sequences.GRCh38.bed.gz")

#: base code -> 0..3 for A, C, G, T, so that 3 - i is the complement; -1 for anything else
IDX4 = np.full(256, -1, dtype=np.int64)
IDX4[C.BASE_A], IDX4[C.BASE_C], IDX4[C.BASE_G], IDX4[C.BASE_T] = 0, 1, 2, 3
#: base code -> complement code
COMP = np.zeros(16, dtype=np.uint8)
for _a, _b in ((C.BASE_A, C.BASE_T), (C.BASE_C, C.BASE_G), (C.BASE_G, C.BASE_C), (C.BASE_T, C.BASE_A), (C.BASE_N, C.BASE_N)):
    COMP[_a] = _b
LETTER = {C.BASE_A: "A", C.BASE_C: "C", C.BASE_G: "G", C.BASE_T: "T"}


def runs(mask: np.ndarray, min_len: int) -> tuple[np.ndarray, np.ndarray]:
    """Starts and ends (0-based, end excluded) of the runs of True at least ``min_len`` long."""
    d = np.diff(np.concatenate(([0], mask.view(np.int8), [0])))
    s, e = np.flatnonzero(d == 1), np.flatnonzero(d == -1)
    keep = (e - s) >= min_len
    return s[keep], e[keep]


def mask_of(n: int, starts, ends) -> np.ndarray:
    m = np.zeros(n, dtype=bool)
    for s, e in zip(np.asarray(starts).tolist(), np.asarray(ends).tolist()):
        m[max(s, 0) : min(e, n)] = True
    return m


def read_runs(path: str = RUNS_BED) -> dict[str, list[tuple[int, int, str]]]:
    """The GA- and TC-sequences, by chromosome: (start, end, kind)."""
    out: dict[str, list[tuple[int, int, str]]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            ch, s, e, kind = line.split("\t")[:4]
            out.setdefault(ch, []).append((int(s), int(e), kind))
    return out


def reference(ref_dir: str, ch: str) -> GpxilFile:
    """The Gpxil browser's gpxil file of a GRCh38 chromosome (a directory, or its URL)."""
    return GpxilFile(f"{ref_dir.rstrip('/')}/{ch}.gpxil.tif")


class CallSet:
    """HG002's calls, one gpxil file a chromosome (HG002_chrN.gpxil.tif, encoded by encode_hg002.sh), and the BED
    file of the regions its benchmark calls with confidence."""

    def __init__(self, name: str, gpxil_dir: str, bed: str):
        self.name, self.dir, self.bed = name, gpxil_dir, bed

    def path(self, ch: str) -> str:
        return os.path.join(self.dir, f"HG002_{ch}.gpxil.tif")

    def open(self, ch: str) -> GpxilFile:
        return GpxilFile(self.path(ch))

    def callable(self, ch: str, n: int) -> np.ndarray:
        m = load_bed(self.bed, ch, n)
        return np.zeros(n, dtype=bool) if m is None else m


def haplotypes(gt: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return C.normalize_states(C.hap1(gt)), C.normalize_states(C.hap2(gt))


def variant_bases(h1, h2, flags) -> np.ndarray:
    """Bases where either haplotype differs from the reference: a substitution, a deleted base, a COMPLEX allele, or
    an insertion after the base."""
    other = lambda h: (h != C.HAP_REF) & (h != C.HAP_UNKNOWN) & (h != C.HAP_ABSENT)  # noqa: E731
    return other(h1) | other(h2) | ((flags & (C.FLAG_INS1 | C.FLAG_INS2)) != 0)


def purine_strand(codes: np.ndarray, kind: str) -> np.ndarray:
    """A GA-sequence's bases as they are; a TC-sequence's, reverse-complemented."""
    return codes if kind == "GA" else COMP[codes[::-1]]


def write_json(name: str, data) -> str:
    import json

    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, name)
    with open(path, "w") as fh:
        json.dump(data, fh, indent=1)
        fh.write("\n")
    return path
