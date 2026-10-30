"""Substitutions by distance from the nearest tandem repeat, observed against expected, from shadow.py's counts.

The expected count at each distance is its callable bases times the rate of their trinucleotides, the rate being HG002's
substitutions of that kind per callable base 201 bases or more from any locus (the far bases). New substitutions are
set against all bases outside low-complexity regions (not only HG002's callable ones), with their own rate 201 to 1,000
bases from any locus: farther out lie large regions where the two studies called almost none.

For the loci where HG002 varies on one haplotype only, a heterozygous substitution beside them is in cis when it is on
the haplotype that varies in the locus. Its share in cis far from the locus measures how the two haplotypes' histories
differ there; if part of the excess beside the locus arose with no regard to the haplotypes (a share in cis of 1/2),
the share beside it falls by that part.

    python3 -I distances.py          -> results/distances.json     (reads results/shadow.npz)
"""

from __future__ import annotations

import json
import os

import numpy as np
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
CLASSES = ["below 1e-6", "1e-6 to 1e-5", "1e-5 to 1e-4", "1e-4 or more"]
FAR = slice(7, 10)  #: 201 bases or more
FAR_NEW = slice(7, 9)  #: 201 to 1,000 bases: farther out lie large regions where the studies called almost no new mutations
MID = [4, 5, 6, 7]  #: 21 to 500 bases: the reference for the loci's own surroundings


def poisson_ci(k: int, level: float = 0.95) -> tuple[float, float]:
    a = 1 - level
    return (chi2.ppf(a / 2, 2 * k) / 2 if k > 0 else 0.0), chi2.ppf(1 - a / 2, 2 * k + 2) / 2


def oe(o: float, e: float) -> dict:
    lo, hi = poisson_ci(int(round(o)))
    return {"observed": int(round(o)), "expected": round(float(e), 2), "ratio": round(o / e, 4) if e else None,
            "ci": [round(lo / e, 4), round(hi / e, 4)] if e else None}


def main() -> None:
    z = np.load(os.path.join(RESULTS, "shadow.npz"))
    edges = z["edges"]
    nb = len(edges) - 1
    bins = [f"{edges[i] + 1}-{edges[i + 1]}" if i < nb - 1 else f"{edges[i] + 1}+" for i in range(nb)]
    a_call = z["A_call"]
    out = {"bins": bins, "classes": CLASSES, "callable_bases": int(a_call.sum()), "bases": int(z["A_all"].sum())}
    base_far = a_call[:, :, FAR, :].sum(axis=(0, 1, 2))
    for name, x in (("heterozygous", z["A_het"]), ("heterozygous, isolated", z["A_iso"]), ("on either haplotype", z["A_snv"])):
        r_ctx = (x[:, :, FAR].sum(axis=(0, 1, 2)) / np.maximum(base_far, 1)[:, None]).sum(axis=1)
        res = {}

        def add(key, xs, bs):
            res[key] = {bins[b]: oe(xs[..., b, :, :].sum(), (bs[..., b, :].reshape(-1, 32).sum(axis=0) * r_ctx).sum()) for b in range(nb)}

        add("all loci", x, a_call)
        for k, cl in enumerate(CLASSES):
            add(f"nearest locus {cl}", x[k], a_call[k])
        for st, sel in (("HG002 the same in the locus", [0]), ("HG002 varies in the locus", [1, 2, 3]),
                        ("HG002 varies on one haplotype", [1, 2]), ("HG002 varies on both haplotypes", [3])):
            add(st, x[:, sel], a_call[:, sel])
            for k, cl in enumerate(CLASSES):
                add(f"{st}, {cl}", x[k, sel], a_call[k, sel])
        out[name] = res
    # B: the distance to the nearest locus of each class, whatever lies between
    for name, x in (("heterozygous, by class alone", z["B_het"]), ("heterozygous, isolated, by class alone", z["B_iso"])):
        res = {}
        for k, cl in enumerate(CLASSES):
            r_ctx = x[k, FAR].sum(axis=(0, 2)) / np.maximum(z["B_call"][k, FAR, :].sum(axis=0), 1)
            res[cl] = {bins[b]: oe(x[k, b].sum(), (z["B_call"][k, b] * r_ctx).sum()) for b in range(nb)}
        out[name] = res
    # cis or trans, beside loci where HG002 varies on one haplotype
    ct = z["CT"]
    cis = {}
    for iso, lab in ((slice(0, 2), "all"), (slice(1, 2), "isolated")):
        for k in [None] + list(range(len(CLASSES))):
            kk = slice(None) if k is None else k
            row = {}
            for b in range(nb):
                c, t = int(ct[kk, b, iso, 1].sum()), int(ct[kk, b, iso, 0].sum())
                if c + t:
                    p = c / (c + t)
                    se = np.sqrt(p * (1 - p) / (c + t))
                    row[bins[b]] = {"cis": c, "trans": t, "share_in_cis": round(p, 3), "ci": [round(p - 1.96 * se, 3), round(p + 1.96 * se, 3)]}
            cis[f"{lab}, {'all classes' if k is None else CLASSES[k]}"] = row
    out["cis_or_trans"] = cis
    # the excess beside one-sided loci, and the share in cis it would leave if it took no account of the haplotypes
    blind = {}
    for iso, lab, x in ((slice(0, 2), "all", z["A_het"]), (slice(1, 2), "isolated", z["A_iso"])):
        r_ctx = (x[:, :, FAR].sum(axis=(0, 1, 2)) / np.maximum(base_far, 1)[:, None]).sum(axis=1)
        for k in [None] + list(range(len(CLASSES))):
            kk = slice(None) if k is None else k

            def o_e(bb):
                xs, bs = x[kk][..., [1, 2], :, :, :], a_call[kk][..., [1, 2], :, :]
                o = xs[..., bb, :, :].sum()
                e = (bs[..., bb, :].reshape(-1, 32).sum(axis=0) * r_ctx).sum()
                return o, e

            o0, e0 = o_e(0)
            om = sum(o_e(b)[0] for b in MID)
            em = sum(o_e(b)[1] for b in MID)
            rel = (o0 / e0) / (om / em)
            c0, t0 = ct[kk, 0, iso, 1].sum(), ct[kk, 0, iso, 0].sum()
            cm, tm = ct[kk, MID, iso, 1].sum(), ct[kk, MID, iso, 0].sum()
            f = 1 - 1 / rel
            far_share = cm / (cm + tm)
            blind[f"{lab}, {'all classes' if k is None else CLASSES[k]}"] = {
                "excess at 1-2 against 21-500": round(float(rel), 3), "share of 1-2 that is excess": round(float(f), 3),
                "share in cis at 21-500": round(float(far_share), 3), "share in cis at 1-2": round(float(c0 / (c0 + t0)), 3),
                "substitutions at 1-2": int(c0 + t0),
                "share in cis at 1-2 if the excess took no account of haplotypes": round(float((1 - f) * far_share + f / 2), 3) if f > 0 else None,
                "the excess's own share in cis": round(float((c0 - (1 - f) * far_share * (c0 + t0)) / (f * (c0 + t0))), 3) if f > 0 else None,
                "at most this much of the excess in cis by the same event, were the rest blind":
                    round(float(((c0 - (1 - f) * far_share * (c0 + t0)) / (f * (c0 + t0)) - 0.5) / 0.5), 3) if f > 0 else None}
    out["cis_if_the_excess_were_blind"] = blind
    # new substitutions, against all bases outside low-complexity regions
    a_all = z["A_new"]
    sets = [str(s) for s in z["sets"]]
    dn = {}
    for lab, sel in (("all", list(range(len(sets)))), ("Sasani", [i for i, s in enumerate(sets) if s.startswith("Sasani")]),
                     ("Platinum", [i for i, s in enumerate(sets) if s.startswith("Platinum")])):
        d = z["D_A"][sel].sum(axis=0)
        mu = d[:, FAR_NEW].sum(axis=(0, 1)).sum(axis=1) / np.maximum(a_all[:, FAR_NEW].sum(axis=(0, 1)), 1)
        res = {"read": int(d.sum()), "all loci": {bins[b]: oe(d[:, b].sum(), (a_all[:, b].sum(axis=0) * mu).sum()) for b in range(nb)}}
        for k, cl in enumerate(CLASSES):
            res[f"nearest locus {cl}"] = {bins[b]: oe(d[k, b].sum(), (a_all[k, b] * mu).sum()) for b in range(nb)}
            res[f"nearest locus {cl}, within 20 bases"] = oe(d[k, :4].sum(), (a_all[k, :4].sum(axis=0) * mu).sum())
            res[f"nearest locus {cl}, 3 to 20 bases"] = oe(d[k, 1:4].sum(), (a_all[k, 1:4].sum(axis=0) * mu).sum())
        dn[lab] = res
    out["new substitutions"] = dn
    with open(os.path.join(RESULTS, "distances.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for name in ("heterozygous", "heterozygous, isolated"):
        print(name)
        for k, v in out[name].items():
            if k.startswith("HG002") and "," in k:
                continue
            print(f"  {k:34s}", "  ".join(f"{b}:{x['ratio']}" for b, x in v.items()))
    print(json.dumps(blind, indent=0)[:2000])


if __name__ == "__main__":
    main()
