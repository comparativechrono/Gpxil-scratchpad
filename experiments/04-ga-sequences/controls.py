"""Two controls for unit copying (unit_copy.py).

1. **Context.** Some 7-base contexts mutate faster than others, and the bases where a transition would copy a unit
   may be such contexts. Here each base of a GA-sequence is weighted by HG002's rate of transitions at the middle
   base of its 7-base context (read on the purine strand), measured over the callable bases outside GA-sequences and
   their flanks. The expected share of copying transitions is then the rate-weighted share of copying bases.
2. **Neighbours.** If transitions copy a neighbouring unit of their own GA-sequence, the copied 6 bases must be in
   that GA-sequence, not just anywhere in GA-sequences. Here the same test is made with the 6-base words of another
   GA-sequence, drawn at random among those of nearly the same length (within 20 %): the excess should vanish.

    python3 controls.py REF_DIR NAME GPXIL_DIR BED          -> results/controls_NAME.json
"""

import os
import random
import sys
import time

import numpy as np

import ga
from gpxiltools import codes as C

K = 6
TI = {"A": "G", "G": "A"}
W7 = 7  # (context length)
NEAR = 20  # (an isolated substitution has no other variant within this many bases)
CHUNK = 1 << 22


def context_ids(i4: np.ndarray) -> np.ndarray:
    """7-base context ids of the middle bases of a + strand array of 0-3 codes (-1 where any base isn't A/C/G/T),
    on the strand where the middle base is a purine: len(i4) - 6 ids."""
    n = len(i4) - (W7 - 1)
    fwd = np.zeros(n, dtype=np.int64)
    rev = np.zeros(n, dtype=np.int64)
    bad = np.zeros(n, dtype=bool)
    for k in range(W7):
        x = i4[k : k + n]
        bad |= x < 0
        fwd = fwd * 4 + x
        rev = rev * 4 + (3 - i4[W7 - 1 - k : W7 - 1 - k + n])  # (the other strand, read 5' to 3')
    mid = i4[W7 // 2 : W7 // 2 + n]
    ids = np.where((mid == 1) | (mid == 3), rev, fwd)
    ids[bad] = -1
    return ids


def rates(ref_dir: str, cs: ga.CallSet, runs) -> tuple[np.ndarray, np.ndarray]:
    """Bases and transitions (at the middle base) of each 7-base context, over callable bases outside GA-sequences
    and their flanks."""
    bases = np.zeros(4**W7, dtype=np.int64)
    trans = np.zeros(4**W7, dtype=np.int64)
    t0 = time.time()
    for ch in ga.AUTOSOMES:
        with ga.reference(ref_dir, ch) as gr, cs.open(ch) as gh:
            rc, hc = gr[ch], gh[ch]
            n = rc.length
            cal = cs.callable(ch, n)
            rr = runs.get(ch, [])
            away = ~ga.mask_of(n, [s - ga.FLANK for s, _, _ in rr], [e + ga.FLANK for _, e, _ in rr])
            for lo in range(3, n - 3, CHUNK):
                hi = min(n - 3, lo + CHUNK)
                base = rc.read("ref", lo - 3, hi + 3) & 0x0F
                ids = context_ids(ga.IDX4[base])
                gt = hc.read("gt", lo, hi)
                h1, h2 = ga.haplotypes(gt)
                b = base[3:-3]
                s1, s2 = np.isin(h1, C.HAP_SUB_CODES), np.isin(h2, C.HAP_SUB_CODES)
                alt = np.where(s1, h1, h2)
                ti = (s1 | s2) & (((b | alt) == (C.BASE_A | C.BASE_G)) | ((b | alt) == (C.BASE_C | C.BASE_T)))
                m = cal[lo:hi] & away[lo:hi] & (ids >= 0)
                bases += np.bincount(ids[m], minlength=4**W7)
                trans += np.bincount(ids[m & ti], minlength=4**W7)
        print(f"{ch} [{time.time() - t0:.0f}s]", flush=True)
    return bases, trans


def main(ref_dir: str, name: str, gpxil_dir: str, bed: str) -> None:
    from scipy.stats import binomtest

    cs = ga.CallSet(name, gpxil_dir, bed)
    runs = ga.read_runs()
    cache = os.path.join(ga.RESULTS, f"context7_{name}.npz")  # (the genome-wide pass, kept: it takes minutes)
    if os.path.exists(cache):
        z = np.load(cache)
        nb, nt = z["bases"], z["transitions"]
    else:
        nb, nt = rates(ref_dir, cs, runs)
        np.savez_compressed(cache, bases=nb, transitions=nt)
    rate = nt / np.maximum(nb, 1)
    # GA-sequences, on their purine strand: words, each base's context id, callable, and HG002's transitions there
    seqs = []
    for ch in ga.AUTOSOMES:
        with ga.reference(ref_dir, ch) as gr, cs.open(ch) as gh:
            rc, hc = gr[ch], gh[ch]
            cal = cs.callable(ch, rc.length)
            for s, e, kind in runs.get(ch, []):
                if s < 3 or e + 3 > rc.length:
                    continue
                base = rc.read("ref", s - 3, e + 3) & 0x0F
                a, b = max(0, s - NEAR), min(rc.length, e + NEAR)
                gt_w, fl_w = hc.read("gt", a, b), hc.read("flags", a, b)
                g1, g2 = ga.haplotypes(gt_w)
                var = ga.variant_bases(g1, g2, fl_w)
                csum = np.concatenate(([0], np.cumsum(var, dtype=np.int64)))
                idx = np.arange(s, e) - a
                iso = (csum[np.minimum(idx + NEAR + 1, b - a)] - csum[np.maximum(idx - NEAR, 0)]) == 1
                h1, h2 = g1[idx], g2[idx]
                s1, s2 = np.isin(h1, C.HAP_SUB_CODES), np.isin(h2, C.HAP_SUB_CODES)
                snv, alt = s1 | s2, np.where(s1, h1, h2)
                c = cal[s:e]
                if kind == "TC":
                    base = ga.COMP[base[::-1]]
                    snv, alt, c, iso = snv[::-1], ga.COMP[alt[::-1]], c[::-1], iso[::-1]
                ids = context_ids(ga.IDX4[base])
                word = "".join(ga.LETTER[int(x)] for x in base[3:-3])
                ti = np.array([snv[j] and ga.LETTER.get(int(alt[j])) == TI[word[j]] for j in range(len(word))], dtype=bool)
                seqs.append((word, ids, c, ti, iso))
    left = (K - 1) // 2

    def tally(own: bool, rng=None, by_len=None):
        out = {"sites": 0, "copy": 0, "w_sites": 0.0, "w_copy": 0.0, "ti": 0, "ti_copy": 0, "ti_iso": 0, "ti_iso_copy": 0}
        for i, (word, ids, c, ti, iso) in enumerate(seqs):
            src = word
            if not own:
                pool = by_len[len(word)]
                j = rng.choice(pool) if len(pool) > 1 else i
                while j == i and len(pool) > 1:
                    j = rng.choice(pool)
                src = seqs[j][0]
            words = {src[k : k + K] for k in range(len(src) - K + 1)}
            for j in range(left, len(word) - (K - 1 - left)):
                if not c[j] or ids[j] < 0:
                    continue
                w = word[j - left : j - left + K]
                cp = (w[:left] + TI[word[j]] + w[left + 1 :]) in words
                r = rate[ids[j]]
                out["sites"] += 1
                out["copy"] += cp
                out["w_sites"] += r
                out["w_copy"] += r * cp
                if ti[j]:
                    out["ti"] += 1
                    out["ti_copy"] += cp
                    if iso[j]:
                        out["ti_iso"] += 1
                        out["ti_iso_copy"] += cp
        return out

    lengths = [len(w) for w, *_ in seqs]
    by_len = {}
    order = sorted(range(len(seqs)), key=lambda i: lengths[i])
    for L in set(lengths):
        lo_, hi_ = np.searchsorted([lengths[i] for i in order], [0.8 * L, 1.2 * L + 1])
        by_len[L] = [order[k] for k in range(lo_, hi_)]
    res = {}
    for label, own in (("own", True), ("other", False)):
        x = tally(own, random.Random(11), by_len)
        p_sites = x["copy"] / x["sites"]
        p_ctx = x["w_copy"] / x["w_sites"]
        x["share_sites"], x["share_context"] = round(p_sites, 4), round(p_ctx, 4)
        x["observed"] = round(x["ti_copy"] / max(x["ti"], 1), 4)
        x["p_vs_sites"] = float(binomtest(x["ti_copy"], x["ti"], p_sites).pvalue) if x["ti"] else None
        x["p_vs_context"] = float(binomtest(x["ti_copy"], x["ti"], p_ctx).pvalue) if x["ti"] else None
        x["observed_isolated"] = round(x["ti_iso_copy"] / max(x["ti_iso"], 1), 4)
        x["p_isolated_vs_context"] = float(binomtest(x["ti_iso_copy"], x["ti_iso"], p_ctx).pvalue) if x["ti_iso"] else None
        x["w_sites"], x["w_copy"] = round(x["w_sites"], 3), round(x["w_copy"], 3)
        res[label] = x
        print(f"{name}, words of the {'same' if own else 'another'} GA-sequence: transitions copying {100 * x['observed']:.1f} % "
              f"of {x['ti']:,}; expected {100 * p_sites:.1f} % from the bases (p = {x['p_vs_sites']:.2g}), "
              f"{100 * p_ctx:.1f} % weighted by 7-base context rates (p = {x['p_vs_context']:.2g}); isolated transitions "
              f"{100 * x['observed_isolated']:.1f} % of {x['ti_iso']:,} (p = {x['p_isolated_vs_context']:.2g} against the context share)")
    ga.write_json(f"controls_{name}.json", {"call_set": name, "k": K, "context": W7, **res})


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
