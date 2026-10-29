"""HG002's variants inside the pure GA/TC-sequences of GRCh38, in their flanks (400 bases either side, outside any of
them) and in all bases, within the regions a benchmark calls with confidence.

Albrecht-Buehler (2008) found that pure GA-sequences "contained many times fewer point mutations than one should
expect from the genome at large", counting from the sequences themselves. Here the point mutations are a person's own.
A GA-sequence holds only A:T pairs and G:C pairs outside CpG (it has no C on its strand), and CpG pairs change ten
times faster, so the counts are kept by base pair class (A:T, G:C outside CpG, CpG), and each region's substitutions
are compared with what its trinucleotide contexts lead one to expect, from the rate of each substitution in each
context across all callable bases. A substitution is "isolated" when no other variant (a substitution, a deleted
base, an insertion or a COMPLEX allele, on either haplotype) lies within 20 bases; clustered ones may be a change of a
repeat written as several substitutions.

    python3 variants.py REF_DIR NAME GPXIL_DIR BED

NAME names the call set (v4.2.1, v5.0q); GPXIL_DIR holds its HG002_chrN.gpxil.tif (encode_hg002.sh), BED its regions.
Writes results/variants_NAME.json, and prints the comparison.
"""

import sys
import time

import numpy as np

import ga
from gpxiltools import codes as C

REGIONS = ("inside", "flank", "all")
CLASSES = ("AT", "GC", "CpG")
NEAR = 20
NAMES = "ACGT"


def chunk_counts(base, h1, h2, flags, cal, inside, flank, lo, hi, tot) -> None:
    """Counts over [lo, hi). The arrays hold NEAR bases more on each side (index NEAR is lo); bases beyond the
    chromosome are 0."""
    p = NEAR
    sl = slice(p, p + hi - lo)
    b = base[sl]
    at = (b == C.BASE_A) | (b == C.BASE_T)
    gc = (b == C.BASE_C) | (b == C.BASE_G)
    cg = (base[:-1] == C.BASE_C) & (base[1:] == C.BASE_G)
    cpg = cg[p : p + hi - lo] | cg[p - 1 : p - 1 + hi - lo]
    cls = {"AT": at, "GC": gc & ~cpg, "CpG": cpg}
    s1, s2 = np.isin(h1[sl], C.HAP_SUB_CODES), np.isin(h2[sl], C.HAP_SUB_CODES)
    snv = s1 | s2
    hom = s1 & s2 & (h1[sl] == h2[sl])
    alt = np.where(s1, h1[sl], h2[sl])
    ti = snv & (((b | alt) == (C.BASE_A | C.BASE_G)) | ((b | alt) == (C.BASE_C | C.BASE_T)))
    dele = (h1 == C.HAP_DEL) | (h2 == C.HAP_DEL)
    indel = (dele[sl] & ~dele[p - 1 : p - 1 + hi - lo]) | ((flags[sl] & (C.FLAG_INS1 | C.FLAG_INS2)) != 0)
    var = ga.variant_bases(h1, h2, flags)
    cs = np.concatenate(([0], np.cumsum(var, dtype=np.int64)))
    idx = np.arange(hi - lo) + p
    iso = snv & ((cs[idx + NEAR + 1] - cs[idx - NEAR]) == 1)  # (the substitution is the only variant base near it)
    real = (b != C.BASE_N) & (b != 0)
    # trinucleotide contexts, read on the strand where the middle base is a purine; the substitution likewise
    i4 = ga.IDX4[base[p - 1 : p + hi - lo + 1]]
    prv, mid, nxt = i4[:-2], i4[1:-1], i4[2:]
    good = (prv >= 0) & (mid >= 0) & (nxt >= 0)
    flip = (mid == 1) | (mid == 3)
    ctx = np.where(good, np.where(flip, 3 - nxt, prv) * 16 + np.where(flip, 3 - mid, mid) * 4 + np.where(flip, 3 - prv, nxt), -1)
    a4 = ga.IDX4[alt]
    alt_i = np.where(snv & (a4 >= 0), np.where(flip, 3 - a4, a4), 0)
    for rg, rmask in (("inside", inside[sl]), ("flank", flank[sl]), ("all", None)):
        m = cal[sl] if rmask is None else (rmask & cal[sl])
        d = tot.setdefault(rg, {"bp": 0, "bp_all": 0, "indels": 0, "classes": {}})
        d["bp"] += int(m.sum())
        d["bp_all"] += int(real.sum() if rmask is None else (rmask & real).sum())
        d["indels"] += int((indel & m).sum())
        for k, cm in cls.items():
            mm = m & cm
            x = d["classes"].setdefault(k, {"bp": 0, "snv": 0, "het": 0, "hom": 0, "ti": 0, "isolated": 0})
            x["bp"] += int(mm.sum())
            x["snv"] += int((snv & mm).sum())
            x["hom"] += int((hom & mm).sum())
            x["het"] += int((snv & ~hom & mm).sum())
            x["ti"] += int((ti & mm).sum())
            x["isolated"] += int((iso & mm).sum())
        cx = d.setdefault("ctx", {"bp": np.zeros(64, np.int64), "snv": np.zeros((64, 4), np.int64), "iso": np.zeros((64, 4), np.int64)})
        ok = m & (ctx >= 0)
        cx["bp"] += np.bincount(ctx[ok], minlength=64)
        for key, sel in (("snv", snv & ok), ("iso", iso & ok)):
            cx[key] += np.bincount(ctx[sel] * 4 + alt_i[sel], minlength=256).reshape(64, 4)


def count(ref_dir: str, cs: ga.CallSet) -> dict:
    runs = ga.read_runs()
    tot: dict = {}
    t0 = time.time()
    for ch in ga.AUTOSOMES:
        with ga.reference(ref_dir, ch) as gr, cs.open(ch) as gh:
            rc, hc = gr[ch], gh[ch]
            n = rc.length
            cal = cs.callable(ch, n)
            rr = runs.get(ch, [])
            s = np.array([x[0] for x in rr], dtype=np.int64)
            e = np.array([x[1] for x in rr], dtype=np.int64)
            inside = ga.mask_of(n, s, e)
            flank = ga.mask_of(n, s - ga.FLANK, e + ga.FLANK) & ~inside
            for lo in range(0, n, ga.CHUNK):
                hi = min(n, lo + ga.CHUNK)
                a, b = max(0, lo - NEAR), min(n, hi + NEAR)

                def padded(x, fill=0):
                    left = np.full(NEAR - (lo - a), fill, dtype=x.dtype)
                    right = np.full(NEAR - (b - hi), fill, dtype=x.dtype)
                    return np.concatenate((left, x, right))

                base = padded(rc.read("ref", a, b) & 0x0F)
                gt = padded(hc.read("gt", a, b), C.GT_UNKNOWN)
                h1, h2 = ga.haplotypes(gt)
                chunk_counts(base, h1, h2, padded(hc.read("flags", a, b)), padded(cal[a:b]), padded(inside[a:b]),
                             padded(flank[a:b]), lo, hi, tot)
        print(f"{ch}: GA/TC bases called {int((inside & cal).sum()):,} of {int(inside.sum()):,} [{time.time() - t0:.0f}s]", flush=True)
    for d in tot.values():
        d["ctx"] = {k: v.tolist() for k, v in d["ctx"].items()}
    return tot


def compare(tot: dict) -> dict:
    """Observed against expected substitutions, inside and in the flanks, from the rate of each substitution in each
    context over all callable bases: all substitutions, and isolated ones alone."""
    from scipy.stats import binomtest, poisson

    out = {}
    for key in ("snv", "iso"):
        bp_all, snv_all = np.array(tot["all"]["ctx"]["bp"]), np.array(tot["all"]["ctx"][key])
        rate = snv_all / np.maximum(bp_all, 1)[:, None]
        for rg in ("inside", "flank"):
            bp, obs = np.array(tot[rg]["ctx"]["bp"]), np.array(tot[rg]["ctx"][key])
            exp = bp[:, None] * rate
            o, x = int(obs.sum()), float(exp.sum())
            row = {"observed": o, "expected": round(x, 1), "ratio": round(o / x, 3) if x else None,
                   "p": float(min(1.0, 2 * min(poisson.cdf(o, x), poisson.sf(o - 1, x)))) if x else None, "types": {}}
            for mb in (0, 2):  # (A and G in the middle)
                for ab in range(4):
                    if ab == mb:
                        continue
                    sel = [c for c in range(64) if (c // 4) % 4 == mb]
                    oo, ee = int(obs[sel, ab].sum()), float(exp[sel, ab].sum())
                    share = binomtest(oo, o, ee / x).pvalue if o and x else None
                    row["types"][f"{NAMES[mb]}>{NAMES[ab]}"] = {"observed": oo, "expected": round(ee, 1),
                                                                 "ratio": round(oo / ee, 3) if ee else None, "share_p": share}
            out[f"{key} {rg}"] = row
    return out


def report(name: str, tot: dict, cmp: dict) -> None:
    from scipy.stats import binomtest

    print(f"\n{name}")
    for rg in REGIONS:
        d = tot[rg]
        print(f"{rg:7s} called {d['bp']:>13,} of {d['bp_all']:>13,} bp ({100 * d['bp'] / max(d['bp_all'], 1):.1f} %); "
              f"indels {d['indels']:,} = {1000 * d['indels'] / max(d['bp'], 1):.3f}/kb")
        for k in CLASSES:
            x = d["classes"][k]
            print(f"   {k:4s} {x['bp']:>13,} bp  SNV {x['snv']:>9,} = {1000 * x['snv'] / max(x['bp'], 1):.3f}/kb  "
                  f"(isolated {x['isolated']:,}; het {x['het']:,}, hom {x['hom']:,}; transitions {100 * x['ti'] / max(x['snv'], 1):.0f} %)")
    a, b = tot["inside"], tot["all"]
    pr = b["indels"] / b["bp"]
    res = binomtest(a["indels"], a["bp"], pr)
    ci = res.proportion_ci()
    print(f"indels inside / all: {a['indels'] / a['bp'] / pr:.2f} (95 % CI {ci.low / pr:.2f}-{ci.high / pr:.2f})")
    for key, row in cmp.items():
        print(f"{key}: observed {row['observed']:,}, expected from contexts {row['expected']:,} (ratio {row['ratio']}, p = {row['p']:.2g})")
        print("   " + "; ".join(f"{t} {v['observed']}/{v['expected']} ({v['ratio']})" for t, v in row["types"].items()))


def main(ref_dir: str, name: str, gpxil_dir: str, bed: str) -> None:
    tot = count(ref_dir, ga.CallSet(name, gpxil_dir, bed))
    cmp = compare(tot)
    ga.write_json(f"variants_{name}.json", {"call_set": name, "counts": tot, "expected": cmp})
    report(name, tot, cmp)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
