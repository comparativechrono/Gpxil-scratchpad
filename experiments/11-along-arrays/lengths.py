"""Where along a repeat array its alleles differ in length from the reference: the Platinum Pedigree parents' alleles
(TRGT) and HG002's (GIAB v5.0q, read from its gpxil files), at the TE-tail arrays of tails.py, read outward from the
element; and, for every catalogue locus the family's calls cover, whether the difference lies in the locus's longest
stretch in frame.

An allele that differs from the reference in length differs somewhere in a stretch of the reference that the two
sequences can't narrow further: the bases between their longest shared start and longest shared end. In a pure array
that is the whole array; where an impurity interrupts it, the stretch stops there. The stretch is given as a share of
the array from its proximal end (next to the element). A difference is localised when its stretch is at most half
the array, proximal when the stretch lies in the half next to the element, distal when in the far half. Placed at
random along the array, a localised stretch would be proximal as often as distal, and would touch the longest stretch
in frame as often as the expected share computed for each one says.

    python3 -I lengths.py REF_DIR HG002_DIR TRGT_VCF FEATURES     -> results/lengths.json, results/windows.tsv.gz

REF_DIR: the Gpxil browser's GRCh38 gpxil files. HG002_DIR: HG002_chrN.gpxil.tif from GIAB's v5.0q small-variant VCF,
written without the reference (experiment 04's encode_hg002.sh). TRGT_VCF: the Platinum Pedigree's pedigree-consistent
tandem repeats (TRGT; its open-access samples), ceph_1463_tandem_repeats.oa.vcf.gz. FEATURES: experiment 06's locus measures (for loci beyond the tails).
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(EXP, "04-ga-sequences"))

import ga  # noqa: E402  (experiment 04)
from gpxiltools import codes as GC  # noqa: E402
from gpxiltools.container import open_gpxil  # noqa: E402
from tails import AGES, FAMILIES, age_label, in_frame_longest  # noqa: E402

RESULTS = os.path.join(HERE, "results")
PARENTS = ("NA12877", "NA12878")  #: the G2 couple: the family's four founding haplotypes in the open-access calls
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def window(r: str, a: str) -> tuple[int, int, bool]:
    """The stretch [lo, hi) of the reference r where its difference from allele a can lie, and whether one insertion or
    deletion explains it: then the span of the places it could be put; otherwise the bases between the longest shared
    start and end (where, for a difference that also changes bases, the stretch can be narrower than where the change
    in length could lie)."""
    n, m = len(r), len(a)
    P = 0
    while P < min(n, m) and r[P] == a[P]:
        P += 1
    S = 0
    while S < min(n, m) and r[n - 1 - S] == a[m - 1 - S]:
        S += 1
    if n > m:  # D bases fewer: deleted bases r[i:i+D] for i from n - S - D to P
        D = n - m
        lo, last = max(0, n - S - D), min(P, n - D)
        if lo <= last:
            return lo, last + D, True
        return min(P, n - S), max(P, n - S), False
    if m > n:  # inserted between r[i-1] and r[i] for i from n - S to P
        lo, hi = max(0, n - S), min(P, n)
        if lo <= hi:
            return lo, hi, True
        return min(P, n - S), max(P, n - S), False
    return min(P, n - S), max(P, n - S), False


def oriented(lo: int, hi: int, n: int, strand: str) -> tuple[int, int]:
    return (lo, hi) if strand == "+" else (n - hi, n - lo)


def family(vcf: str, tails, others: dict) -> tuple[list[tuple], dict]:
    """The parents' alleles that differ in length from the reference, at the tails (oriented) and at every locus
    (unoriented, with its longest stretch in frame). others: trid -> (shortest motif, start, end), for every locus."""
    tl = {t: i for i, t in enumerate(tails.trid)}
    start = tails.start.to_numpy()
    out, allloci, skipped = [], [], Counter()
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                ip = [samples.index(p) for p in PARENTS]
                continue
            trid = f[7].split("TRID=")[1].split(";")[0]
            k = tl.get(trid)
            p = others.get(trid)
            if k is None and p is None:
                skipped["not in the catalogue's autosomal loci"] += 1
                continue
            s0, e0 = (int(tails.start.iat[k]), int(tails.end.iat[k])) if k is not None else p[1:]
            ref = f[3]
            alle = [ref] + (f[4].split(",") if f[4] != "." else [])
            pad = int(f[1]) == s0  # (POS, 1-based, is the base before the locus: a padding base leads every allele)
            if pad:
                alle = [x[1:] for x in alle]
            ref = alle[0]
            if len(ref) != e0 - s0:
                skipped["reference allele not the locus"] += 1
                continue
            # each of the parents' haplotypes counts once (two haplotypes with one allele count twice)
            haps = []
            for j in ip:
                for g in f[9 + j].split(":")[0].replace("|", "/").split("/"):
                    if g not in (".", "0") and len(alle[int(g)]) != len(ref):
                        haps.append(alle[int(g)])
            if not haps:
                continue
            n = len(ref)
            if p is not None:
                bl, bh = in_frame_longest(ref.upper(), int(p[0]))
                for a in haps:
                    lo, hi, single = window(ref, a)
                    allloci.append((n, len(a) - n, lo, hi, bl, bh, int(single)))
            if k is not None:
                st = tails.te_strand.iat[k]
                for a in haps:
                    lo, hi, single = window(ref, a)
                    lo, hi = oriented(lo, hi, n, st)
                    out.append(("family", trid, k, len(a) - n, lo, hi, int(single)))
    return out, {"skipped": dict(skipped), "all_loci": allloci}


SUB_LETTER = {GC.BASE_A: "A", GC.BASE_C: "C", GC.BASE_G: "G", GC.BASE_T: "T"}


def haplotype(ref: bytes, states: np.ndarray, ins_flag: np.ndarray, inserted: dict, hap: int, a: int, b: int) -> str:
    """Haplotype ``hap``'s sequence over the reference's [a, b), from its states (the gt plane's half), its insertion
    flags and the inserted sequences (the sidecar's, decoded once a chromosome): as gpxiltools' haplotype_sequence
    builds it (unknown and complex bases N; deleted bases dropped; an insertion after the base it follows)."""
    out = []
    for i in range(a, b):
        x = int(states[i])
        if x == GC.HAP_REF:
            out.append(chr(ref[i]))
        elif x in SUB_LETTER:
            out.append(SUB_LETTER[x])
        elif x not in (GC.HAP_DEL, GC.HAP_ABSENT):
            out.append("N")
        if ins_flag[i]:
            out.append(inserted.get((hap, i), "N").upper())
    return "".join(out)


def hg002(ref_dir: str, hg_dir: str, tails) -> tuple[list[tuple], dict]:
    """HG002's haplotypes that differ in length from the reference over a tail array, among the arrays the benchmark
    calls on both haplotypes (and the base either side)."""
    out, counts = [], Counter()
    t0 = time.time()
    for c in ga.AUTOSOMES:
        idx = np.flatnonzero(tails.chrom.to_numpy() == c)
        if not len(idx):
            continue
        with open_gpxil(os.path.join(hg_dir, f"HG002_{c}.gpxil.tif"), reference=os.path.join(ref_dir, f"{c}.gpxil.tif")) as h:
            hc = h.find(c)
            n_c = hc.length
            gt = hc.read("gt", 0, n_c)
            h1, h2 = GC.hap1(gt), GC.hap2(gt)
            del gt
            fl = hc.read("flags", 0, n_c)
            i1, i2 = (fl & GC.FLAG_INS1) != 0, (fl & GC.FLAG_INS2) != 0
            del fl
            inserted = {(e.hap, e.offset): e.data for e in hc.sidecar() if e.kind == "INS"}
            ref = bytes(GC.REF_TO_ASCII[hc.read("ref", 0, n_c) & 0x0F])
        unknown = np.flatnonzero((h1 == GC.HAP_UNKNOWN) | (h2 == GC.HAP_UNKNOWN))
        dels = np.flatnonzero((h1 == GC.HAP_DEL) | (h2 == GC.HAP_DEL))
        ins = np.flatnonzero(i1 | i2)
        st_, en_ = tails.start.to_numpy()[idx], tails.end.to_numpy()[idx]
        a_, b_ = np.maximum(st_ - 1, 0), np.minimum(en_ + 1, n_c)
        called = np.searchsorted(unknown, a_, "left") == np.searchsorted(unknown, b_, "left")
        # a deleted base inside the array, or an insertion after one of its bases or after the base before it
        changed = ((np.searchsorted(dels, st_, "left") < np.searchsorted(dels, en_, "left"))
                   | (np.searchsorted(ins, a_, "left") < np.searchsorted(ins, en_, "left")))
        counts["arrays not called on both haplotypes"] += int((~called).sum())
        counts["arrays called"] += int(called.sum())
        for k, s, e, a, b in zip(idx[called & changed], st_[called & changed], en_[called & changed], a_[called & changed], b_[called & changed]):
            r = ref[a:b].decode()
            n = int(e - s)
            st = tails.te_strand.iat[k]
            for hap, states, flag in ((1, h1, i1), (2, h2, i2)):
                x = haplotype(ref, states, flag, inserted, hap, int(a), int(b))
                if len(x) == len(r):
                    continue
                lo, hi, single = window(r, x)
                lo, hi = min(max(lo - int(s - a), 0), n), min(max(hi - int(s - a), 0), n)
                lo, hi = oriented(lo, hi, n, st)
                out.append(("HG002", tails.trid.iat[k], int(k), len(x) - len(r), lo, hi, int(single)))
                counts["haplotypes differing in length"] += 1
        del h1, h2, i1, i2, ref, inserted
        print(c, dict(counts), f"[{time.time() - t0:.0f}s]", flush=True)
    return out, dict(counts)


def expected_touch(n: np.ndarray, lo: np.ndarray, hi: np.ndarray, bl: np.ndarray, bh: np.ndarray) -> np.ndarray:
    """For each stretch [lo, hi) of an array of n bases, the chance that a stretch as wide, put at one of the n - w + 1
    places along the array at random, touches [bl, bh) (overlaps it or meets one of its ends, as touches() counts)."""
    w = hi - lo
    span = np.maximum(n - w, 0)
    a, b = np.maximum(0, bl - w), np.minimum(span, bh)
    return np.clip(b - a + 1, 0, None) / (span + 1)


def touches(lo, hi, bl, bh):
    return (lo <= bh) & (hi >= bl)


def halves(n, lo, hi):
    """Near: the stretch lies in the half next to the element; far: in the other half (a stretch of no width at the very
    middle is in neither)."""
    near = (hi <= n / 2) & (lo < n / 2)
    far = (lo >= n / 2) & (hi > n / 2)
    return near, far


def counts(loc, near, far, m) -> dict:
    x = loc & m
    d = {"localised": int(x.sum()), "near": int((x & near).sum()), "far": int((x & far).sum())}
    d["near_to_far"] = round(d["near"] / d["far"], 3) if d["far"] else None
    return d


def summarise(w) -> dict:
    n = w.n.to_numpy()
    lo, hi = w.lo.to_numpy(), w.hi.to_numpy()
    bl, bh = w.best_lo.to_numpy(), w.best_hi.to_numpy()
    width = (hi - lo) / n
    loc = width <= 0.5
    near, far = halves(n, lo, hi)
    touch = touches(lo, hi, bl, bh)
    best_near, best_far = bh <= n / 2, bl >= n / 2
    every = np.ones(len(w), dtype=bool)
    d = {"differences": int(len(w)), "arrays": int(w.trid.nunique()),
         "gains": int((w["diff"] > 0).sum()), "losses": int((w["diff"] < 0).sum()),
         "single_events": int(w.single.sum()),
         "stretch_share_of_array_quartiles": [round(float(x), 4) for x in np.quantile(width, [0.25, 0.5, 0.75])] if len(w) else None,
         **counts(loc, near, far, every)}
    d["across_the_middle"] = d["localised"] - d["near"] - d["far"]
    if loc.any():
        d["touch_longest_in_frame"] = int((loc & touch).sum())
        d["touch_longest_in_frame_expected_at_random"] = round(float(expected_touch(n[loc], lo[loc], hi[loc], bl[loc], bh[loc]).sum()), 1)
        d["longest in frame wholly in the near half"] = counts(loc, near, far, best_near)
        d["longest in frame wholly in the far half"] = counts(loc, near, far, best_far)
        d["longest in frame across the middle"] = counts(loc, near, far, ~best_near & ~best_far)
        d["not touching the longest in frame"] = counts(loc, near, far, ~touch)
        d["single events"] = counts(loc, near, far, w.single.to_numpy() == 1)
        d["gains"] = counts(loc, near, far, w["diff"].to_numpy() > 0)
        d["losses"] = counts(loc, near, far, w["diff"].to_numpy() < 0)
        d["by_te_strand"] = {s: counts(loc, near, far, w.te_strand.to_numpy() == s) for s in "+-"}
        # (by which half of the array is richer in A in the reference: the purest stretch is one measure of purity)
        an, af = w.a_near.to_numpy(), w.a_far.to_numpy()
        d["near half richer in A"] = counts(loc, near, far, an > af + 0.05)
        d["halves within 0.05 in A"] = counts(loc, near, far, np.abs(an - af) <= 0.05)
        d["far half richer in A"] = counts(loc, near, far, af > an + 0.05)
        a, b = d["longest in frame wholly in the near half"], d["longest in frame wholly in the far half"]
        if a["near"] and a["far"] and b["near"] and b["far"]:
            # (were changes to follow the purest stretch alone, the near-to-far ratio where it lies near would be the
            # far-to-near ratio where it lies far; the geometric mean of their quotient is what the near end adds)
            d["near_end_excess_beyond_the_purest_stretch"] = round(float(np.sqrt(a["near_to_far"] * b["near_to_far"])), 3)
    return d


def main(ref_dir: str, hg_dir: str, vcf: str, features: str) -> None:
    import pandas as pd
    tails = pd.read_csv(os.path.join(RESULTS, "tail_loci.tsv.gz"), sep="\t")
    print("tail arrays", len(tails), flush=True)
    with gzip.open(vcf, "rt") as fh:  # (the loci the family's calls cover: only those need their motifs)
        called = {line.split("TRID=")[1].split(";")[0] for line in fh if not line.startswith("#")}
    others = {}
    for f in pd.read_csv(features, sep="\t", usecols=["chrom", "start", "end", "trid", "motif_min"],
                         dtype={"chrom": "category", "motif_min": "int32"}, chunksize=1_000_000):
        f = f[f.chrom.isin(ga.AUTOSOMES) & f.trid.isin(called)]
        others.update(zip(f.trid, zip(f.motif_min, f.start, f.end)))
    del called
    fam_rows, fam_info = family(vcf, tails, others)
    print("family differences at tails", len(fam_rows), fam_info["skipped"], flush=True)
    hg_rows, hg_counts = hg002(ref_dir, hg_dir, tails)
    w = pd.DataFrame(fam_rows + hg_rows, columns=["source", "trid", "k", "diff", "lo", "hi", "single"])
    for col in ("fam", "te_strand", "age", "div", "best_lo", "best_hi"):
        w[col] = tails[col].to_numpy()[w.k.to_numpy()]
    # each half's share of A in the reference (the mean of its tenths, read from the element)
    w["a_near"] = tails[[f"a{i}" for i in range(5)]].to_numpy().mean(axis=1)[w.k.to_numpy()].round(4)
    w["a_far"] = tails[[f"a{i}" for i in range(5, 10)]].to_numpy().mean(axis=1)[w.k.to_numpy()].round(4)
    w["n"] = (tails.end - tails.start).to_numpy()[w.k.to_numpy()]
    w = w.drop(columns="k")
    os.makedirs(RESULTS, exist_ok=True)
    w.to_csv(os.path.join(RESULTS, "windows.tsv.gz"), sep="\t", index=False)
    out = {"parents": list(PARENTS), "family_skipped": fam_info["skipped"], "hg002": hg_counts, "tails": {}}
    for src in ("family", "HG002"):
        for fam in FAMILIES:
            x = w[(w.source == src) & (w.fam == fam)]
            out["tails"].setdefault(src, {})[fam] = {"all": summarise(x)}
            for age in [age_label(a + 0.01) for a in AGES[:-1]]:
                y = x[x.age == age]
                if len(y):
                    out["tails"][src][fam][age] = summarise(y)
    a = np.array(fam_info["all_loci"], dtype=np.int64).reshape(-1, 7)
    n, lo, hi, bl, bh = a[:, 0], a[:, 2], a[:, 3], a[:, 4], a[:, 5]
    loc = (hi - lo) / n <= 0.5
    touch = touches(lo, hi, bl, bh)
    out["all_loci_family"] = {"differences": int(len(a)), "localised": int(loc.sum()),
                              "touch_longest_in_frame": int((loc & touch).sum()),
                              "expected_at_random": round(float(expected_touch(n[loc], lo[loc], hi[loc], bl[loc], bh[loc]).sum()), 1),
                              "longest_in_frame_mean_share_of_locus": round(float(((bh - bl) / n)[loc].mean()), 4)}
    with open(os.path.join(RESULTS, "lengths.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1)[:6000])


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:5])
