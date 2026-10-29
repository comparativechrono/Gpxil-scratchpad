"""Do the substitutions in pure GA-sequences copy a neighbouring repeat unit?

Inside a pure GA-sequence, read on its purine strand, every base has one possible transition (A<->G) and two
transversions, which would put a pyrimidine into it. For each base in the regions a benchmark calls with confidence,
the k bases around it after its transition either occur elsewhere in the same GA-sequence (the transition "copies a
unit") or don't. If HG002's transitions came about by copying (slippage, unequal exchange or gene conversion between
the copies of a motif), they would sit where a transition copies a unit more often than the share of such bases.
Reported for k = 4, 5 and 6 (the window centred on the base, left of centre when k is even), for all of HG002's
transitions there and separately for isolated ones (no other variant within 20 bases) and clustered ones.

    python3 unit_copy.py REF_DIR NAME GPXIL_DIR BED          -> results/unit_copy_NAME.json
"""

import sys
import time

import numpy as np

import ga
from gpxiltools import codes as C

KS = (4, 5, 6)
NEAR = 20
TI = {"A": "G", "G": "A"}


def main(ref_dir: str, name: str, gpxil_dir: str, bed: str) -> None:
    from scipy.stats import binomtest

    cs = ga.CallSet(name, gpxil_dir, bed)
    runs = ga.read_runs()
    res = {k: {"sites": 0, "sites_copy": 0, **{g: {"ti": 0, "ti_copy": 0, "tv": 0} for g in ("all", "isolated", "clustered")}} for k in KS}
    t0 = time.time()
    for ch in ga.AUTOSOMES:
        with ga.reference(ref_dir, ch) as gr, cs.open(ch) as gh:
            rc, hc = gr[ch], gh[ch]
            n = rc.length
            cal = cs.callable(ch, n)
            for s, e, kind in runs.get(ch, []):
                a, b = max(0, s - NEAR), min(n, e + NEAR)
                gt, fl = hc.read("gt", a, b), hc.read("flags", a, b)
                h1, h2 = ga.haplotypes(gt)
                var = ga.variant_bases(h1, h2, fl)
                csum = np.concatenate(([0], np.cumsum(var, dtype=np.int64)))
                idx = np.arange(s, e) - a
                near = csum[np.minimum(idx + NEAR + 1, b - a)] - csum[np.maximum(idx - NEAR, 0)]
                s1, s2 = np.isin(h1[idx], C.HAP_SUB_CODES), np.isin(h2[idx], C.HAP_SUB_CODES)
                snv, alt = s1 | s2, np.where(s1, h1[idx], h2[idx])
                c = cal[s:e]
                base = rc.read("ref", s, e) & 0x0F
                if kind == "TC":  # (on the purine strand)
                    base = ga.COMP[base[::-1]]
                    snv, alt, near, c = snv[::-1], ga.COMP[alt[::-1]], near[::-1], c[::-1]
                seq = "".join(ga.LETTER[int(x)] for x in base)
                for k in KS:
                    left = (k - 1) // 2
                    words = {seq[j : j + k] for j in range(len(seq) - k + 1)}
                    r = res[k]
                    for j in range(left, len(seq) - (k - 1 - left)):
                        if not c[j]:
                            continue
                        w = seq[j - left : j - left + k]
                        copy = (w[:left] + TI[seq[j]] + w[left + 1 :]) in words
                        r["sites"] += 1
                        r["sites_copy"] += copy
                        if not snv[j]:
                            continue
                        is_ti = ga.LETTER.get(int(alt[j])) == TI[seq[j]]
                        for g in ("all", "isolated" if near[j] == 1 else "clustered"):
                            if is_ti:
                                r[g]["ti"] += 1
                                r[g]["ti_copy"] += copy
                            else:
                                r[g]["tv"] += 1
        print(f"{ch} [{time.time() - t0:.0f}s]", flush=True)
    print(f"\n{name}")
    for k, r in res.items():
        p = r["sites_copy"] / max(r["sites"], 1)
        r["expected_share"] = round(p, 4)
        line = [f"k={k}: a transition copies a unit at {100 * p:.1f} % of {r['sites']:,} bases"]
        for g in ("all", "isolated", "clustered"):
            x = r[g]
            if x["ti"]:
                bt = binomtest(x["ti_copy"], x["ti"], p)
                x["share"], x["p"] = round(x["ti_copy"] / x["ti"], 4), float(bt.pvalue)
                x["transitions"] = round(x["ti"] / (x["ti"] + x["tv"]), 4)
                line.append(f"{g} {x['ti_copy']}/{x['ti']} = {100 * x['share']:.1f} % (p = {bt.pvalue:.2g}; transitions {100 * x['transitions']:.0f} % of {x['ti'] + x['tv']})")
        print("; ".join(line))
    ga.write_json(f"unit_copy_{name}.json", {"call_set": name, "k": {str(k): v for k, v in res.items()}})


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
