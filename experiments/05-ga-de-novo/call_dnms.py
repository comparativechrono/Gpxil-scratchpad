"""Candidate de novo mutations of the 16 open trios of CEPH 1463, genome-wide on the autosomes, from the family's
HiFi joint calls (DeepVariant and GLnexus).

An allele is a candidate de novo mutation of a child when:

- the child carries it (het or hom), with GQ >= 20, DP >= 10 and at least a quarter of its reads;
- both parents are called without it, each with GQ >= 20, DP >= 10 and no read of it;
- none of the three has more than 2.5 times its median depth (on chr20), as in collapsed repeats;
- no one else in the family carries it but the child's descendants and siblings (and theirs): a sibling who does
  shares it from a parent's mosaic germline; a grandparent, a spouse or a cousin who does means it was inherited.

For each, the reads and calls of the trio, which siblings carry it, and how many of the child's children do.

    python3 -I call_dnms.py VCF [WORKERS]       -> results/candidates.tsv.gz, results/candidates.json

VCF: the joint calls, CEPH1463.GRCh38.deepvariant.glnexus.oa.vcf.gz, downloaded with its index (1.6 GB; ped.HIFI).
"""

import gzip
import json
import os
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

import ped  # noqa: E402

MIN_GQ, MIN_DP, MIN_AB, MAX_DP_FACTOR = 20, 10, 0.25, 2.5
COLS = ["chrom", "pos", "ref", "alt", "child", "generation", "father", "mother", "kind", "length_change",
        "child_gt", "child_gq", "child_dp", "child_alt_reads", "father_gq", "father_dp", "mother_gq", "mother_dp",
        "siblings_with", "children_with", "children_called"]


def kind_of(ref: str, alt: str) -> tuple[str, int]:
    if len(ref) == len(alt) == 1:
        return "SNV", 0
    if len(ref) == len(alt):
        return "MNV", 0
    d = len(alt) - len(ref)
    return ("insertion" if d > 0 else "deletion"), d


def median_depths(path: str, ch: str = "chr20", n: int = 200_000) -> dict:
    """Each sample's median depth over the first n records of a chromosome where it is called."""
    import cyvcf2

    vcf = cyvcf2.VCF(path)
    dps = []
    for i, r in enumerate(vcf(ch)):
        if i >= n:
            break
        dp = r.format("DP")
        if dp is not None:
            dps.append(np.where(r.genotype.array()[:, 0] >= 0, dp[:, 0], -1))
    a = np.array(dps)
    return {s: float(np.median(a[a[:, i] >= 0, i])) for i, s in enumerate(vcf.samples)}


def scan(args):
    """One chromosome: the candidates, and the records looked at."""
    path, ch, medians = args
    import cyvcf2

    vcf = cyvcf2.VCF(path)
    idx = {s: i for i, s in enumerate(vcf.samples)}
    trios = [(c, idx[c], idx[f], idx[m]) for c, (f, m) in ped.TRIOS.items()]
    max_dp = np.array([MAX_DP_FACTOR * medians[s] for s in vcf.samples])
    allowed = {}  # (the family members who may carry a child's de novo allele)
    for c, (f, m) in ped.TRIOS.items():
        keep = {c} | ped.descendants(c)
        for s in ped.siblings(c):
            keep |= {s} | ped.descendants(s)
        allowed[c] = np.array([s not in keep for s in vcf.samples])  # (True: must not carry it)
    out, records = [], 0
    for r in vcf(ch):
        records += 1
        g = r.genotype.array()[:, :2]
        n_alt = len(r.ALT)
        if n_alt == 0:
            continue
        ad = r.format("AD")
        gq = r.format("GQ")
        dp = r.format("DP")
        if ad is None or gq is None or dp is None:
            continue
        gq, dp = gq[:, 0], dp[:, 0]
        for k in range(1, n_alt + 1):
            has = (g[:, 0] == k) | (g[:, 1] == k)
            if not has.any():
                continue
            for c, ci, fi, mi in trios:
                if not has[ci] or has[fi] or has[mi]:
                    continue
                if (g[[fi, mi]] < 0).any() or (g[ci] < 0).any():
                    continue
                if min(gq[ci], gq[fi], gq[mi]) < MIN_GQ or min(dp[ci], dp[fi], dp[mi]) < MIN_DP:
                    continue
                if (dp[[ci, fi, mi]] > max_dp[[ci, fi, mi]]).any():
                    continue
                if k >= ad.shape[1] or ad[ci, k] < 0 or ad[fi, k] != 0 or ad[mi, k] != 0:
                    continue
                tot = int(ad[ci].clip(0).sum())
                if tot == 0 or ad[ci, k] / tot < MIN_AB:
                    continue
                if has[allowed[c]].any():
                    continue
                kids = [idx[x] for x in ped.children(c)]
                kids_called = [x for x in kids if (g[x] >= 0).all()]
                sibs = [s for s in ped.siblings(c) if has[idx[s]]]
                alt = r.ALT[k - 1]
                kind, d = kind_of(r.REF, alt)
                out.append([ch, r.POS, r.REF, alt, c, ped.GENERATION[c], vcf.samples[fi], vcf.samples[mi], kind, d,
                            "/".join(str(x) for x in g[ci]), int(gq[ci]), int(dp[ci]), int(ad[ci, k]),
                            int(gq[fi]), int(dp[fi]), int(gq[mi]), int(dp[mi]), ",".join(sibs),
                            sum(bool(has[x]) for x in kids_called), len(kids_called)])
    return ch, out, records


def main(path: str, workers: int = 2) -> None:
    t0 = time.time()
    medians = median_depths(path)
    print("median depths:", medians, flush=True)
    rows, records = [], {}
    with Pool(workers) as pool:
        for ch, out, n in pool.imap_unordered(scan, [(path, ch, medians) for ch in ped.AUTOSOMES]):
            rows += out
            records[ch] = n
            print(f"{ch}: {n:,} records, {len(out)} candidates [{time.time() - t0:.0f}s]", flush=True)
    order = {c: i for i, c in enumerate(ped.AUTOSOMES)}
    rows.sort(key=lambda x: (order[x[0]], x[1], x[4]))
    os.makedirs(ped.RESULTS, exist_ok=True)
    with gzip.open(os.path.join(ped.RESULTS, "candidates.tsv.gz"), "wt") as fh:
        fh.write("\t".join(COLS) + "\n")
        for x in rows:
            fh.write("\t".join(map(str, x)) + "\n")
    by_child = {}
    for x in rows:
        b = by_child.setdefault(x[4], {"SNV": 0, "MNV": 0, "insertion": 0, "deletion": 0})
        b[x[8]] += 1
    summary = {"vcf": os.path.basename(path), "records": sum(records.values()), "candidates": len(rows),
               "filters": {"min_gq": MIN_GQ, "min_dp": MIN_DP, "min_child_allele_fraction": MIN_AB, "parent_alt_reads": 0,
                           "max_dp_factor": MAX_DP_FACTOR},
               "median_depths": medians,
               "by_child": by_child}
    with open(os.path.join(ped.RESULTS, "candidates.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    print(json.dumps(summary["by_child"]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 2)
