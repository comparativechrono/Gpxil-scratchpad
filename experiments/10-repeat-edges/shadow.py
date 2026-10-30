"""Do tandem repeats cast a mutational shadow on the bases beside them? The counts.

HG002's gpxil planes (GIAB's v5.0q small-variant benchmark, both haplotypes), read beside GRCh38's reference plane:
every autosomal base outside the Platinum TR catalogue's 7.3 million loci, by its distance to the nearest locus and
that locus's instability, and the substitutions HG002 carries there. Instability is the rate of new TR mutations that
experiment 06's Poisson model (length, motif size, purity and class, fitted to the 440 new TR mutations of eight
children) predicts for the locus, per locus and transmission, in four classes: below 1e-6, 1e-6 to 1e-5, 1e-5 to
1e-4, and 1e-4 or more.

* A: each base is given its nearest locus (any class); whether HG002 varies inside that locus, and on which haplotype;
* B: for each class, each base's distance to the nearest locus of that class.

A base counts when it is not N and has a known trinucleotide; it is callable when it is also in the benchmark's
regions, known on both haplotypes and in no deletion or complex allele. A substitution is a haplotype state that is a
base other than the reference's; heterozygous when the other haplotype is the reference's; isolated when no other
variant (any haplotype state, or an insertion) lies within 10 bases. Counts are kept by distance, class and
trinucleotide (read with a pyrimidine in the middle) and new base. The new substitutions of the CEPH families (Sasani
et al. 2019) and the Platinum Pedigree (Porubsky et al. 2025), read as experiment 07 reads them, are placed the same
way. Sasani et al. called no new mutations in low-complexity regions, so for new substitutions both the bases and the
substitutions are counted outside LCR-hs38 (Li 2014), which stands in for their GRCh37 mask.

    python3 -I shadow.py REF_DIR HG002_DIR BED FEATURES DNM_DIR CHAIN37 TABLE.xlsx CHAIN_T2T LCR
        -> results/shadow.npz, results/shadow_counts.json

HG002_DIR: HG002_chrN.gpxil.tif from v5.0q, encoded as experiment 04's encode_hg002.sh does. BED: v5.0q's benchmark
regions. FEATURES: the locus measures experiment 06 writes with --features. LCR: LCR-hs38.bed.gz (lh3/varcmp).
"""

from __future__ import annotations

import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(EXP, "04-ga-sequences"))
sys.path.insert(0, os.path.join(EXP, "05-ga-de-novo"))
sys.path.insert(0, os.path.join(EXP, "07-copied-units"))

import ga  # noqa: E402  (experiment 04)
from gpxiltools.container import open_gpxil  # noqa: E402

RESULTS = os.path.join(HERE, "results")
STRATA = os.path.join(EXP, "06-repeat-mutability", "results", "strata.tsv.gz")
EDGES = np.array([0, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 10**12])  #: distance bins: 1-2, 3-5, ..., 501-1000, 1001+
NB = len(EDGES) - 1
CLS = np.array([0, 1e-6, 1e-5, 1e-4, np.inf])  #: instability classes, new TR mutations per locus and transmission
NC = len(CLS) - 1
NS = 4  #: HG002 inside the nearest locus: the same as the reference on both haplotypes, varies on hap1, on hap2, on both
ISO = 10
CHUNK = 8_000_000
L06 = [0, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 200, 250, 350, 500, 1000, 10**7]  #: experiment 06's strata
M06 = [0, 1, 2, 3, 4, 5, 6, 12, 10**5]
P06 = [-1.0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.999, 1.0]

IDX = np.full(256, 4, dtype=np.int8)  #: reference code -> A 0, C 1, G 2, T 3, else 4
for _i, _b in enumerate((1, 2, 4, 8)):
    IDX[_b] = IDX[_b | 0x10] = _i
HAP = np.full(16, 4, dtype=np.int8)  #: haplotype state -> the base it holds, else 4
HAP[[1, 2, 4, 8]] = [0, 1, 2, 3]


def label(edges: list, i: int, digits: int = 0) -> str:
    a, b = edges[i], edges[i + 1]
    if digits:
        return f"{max(a, 0):.{digits}f}-{b:.{digits}f}"
    return f"{a + 1}-{b}" if b < edges[-1] else f"{a + 1}+"


def strata_rates() -> dict:
    """(length, motif, purity, class) labels of experiment 06's strata -> the rate its class model predicts."""
    import pandas as pd
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    st = pd.read_csv(STRATA, sep="\t", dtype=str)
    for k in ("loci", "bases", "mutations"):
        st[k] = st[k].astype(int)
    st = st.rename(columns={"class": "cls"})
    fit = smf.glm("mutations ~ C(length) + C(motif) + C(purity) + C(cls)", data=st, family=sm.families.Poisson(),
                  offset=np.log(st.loci * 8)).fit()
    mu = fit.predict(st, offset=np.zeros(len(st)))
    print(f"experiment 06's class model: deviance {fit.deviance:.1f} on {int(fit.df_resid)}", flush=True)
    return {(a, b, c, d): float(x) for a, b, c, d, x in zip(st.length, st.motif, st.purity, st.cls, mu)}


def read_loci(features: str, rates: dict) -> dict:
    """chrom -> (starts, ends, instability class, predicted rate, shortest motif), in start order."""
    import pandas as pd

    f = pd.read_csv(features, sep="\t", usecols=["chrom", "start", "end", "length", "motif_min", "purity", "class"])
    li = np.searchsorted(L06, f.length.values, "left") - 1
    mi = np.searchsorted(M06, f.motif_min.values, "left") - 1
    pi = np.searchsorted(P06, f.purity.values, "left") - 1
    lab_l = [label(L06, i) for i in range(len(L06) - 1)]
    lab_m = [label(M06, i) for i in range(len(M06) - 1)]
    lab_p = [label(P06, i, 3) for i in range(len(P06) - 1)]
    mu = np.array([rates.get((lab_l[a], lab_m[b], lab_p[c], d), np.nan) for a, b, c, d in zip(li, mi, pi, f["class"].values)])
    print("loci without a rate:", int(np.isnan(mu).sum()), flush=True)
    mu = np.nan_to_num(mu, nan=np.nanmedian(mu))
    cls = np.searchsorted(CLS, mu, "right") - 1
    out = {}
    for c, d in f.assign(mu=mu, cls=cls).groupby("chrom"):
        d = d.sort_values("start")
        out[c] = (d.start.values.astype(np.int64), d.end.values.astype(np.int64), d.cls.values.astype(np.int8), d.mu.values,
                  d.motif_min.values.astype(np.int64))
    return out


def read_bed(path: str) -> dict:
    import gzip

    iv = defaultdict(list)
    with (gzip.open(path, "rt") if path.endswith(".gz") else open(path)) as fh:
        for line in fh:
            c, s, e = line.split("\t")[:3]
            iv[c].append((int(s), int(e)))
    return {c: (np.array([x[0] for x in v]), np.array([x[1] for x in v])) for c, v in iv.items()}


def interval_mask(a: int, b: int, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
    m = np.zeros(b - a, dtype=bool)
    i0 = np.searchsorted(ends, a, "right")
    i1 = np.searchsorted(starts, b, "left")
    for s, e in zip(starts[i0:i1], ends[i0:i1]):
        m[max(s, a) - a:min(e, b) - a] = True
    return m


def nearest(q: np.ndarray, s: np.ndarray, e: np.ndarray):
    """For positions outside the loci: the distance to the nearest locus (1 = right beside it) and its index."""
    if not len(s):
        return np.full(len(q), 10**12), np.zeros(len(q), dtype=np.int64)
    i = np.searchsorted(s, q, "right")
    left = i - 1
    dl = np.where(left >= 0, q - e[np.maximum(left, 0)] + 1, 10**12)
    dr = np.where(i < len(s), s[np.minimum(i, len(s) - 1)] - q, 10**12)
    use_left = dl <= dr
    return np.where(use_left, dl, dr), np.where(use_left, left, np.minimum(i, len(s) - 1))


def contexts(x: np.ndarray):
    """x: base indexes of [a-1, b+1). For [a, b): the trinucleotide read with a pyrimidine in the middle (0-31; -1 when
    any base is unknown), and whether it was read on the other strand."""
    l, c, r = x[:-2].astype(np.int16), x[1:-1].astype(np.int16), x[2:].astype(np.int16)
    ok = (l < 4) & (c < 4) & (r < 4)
    flip = (c == 0) | (c == 2)
    l2, c2, r2 = np.where(flip, 3 - r, l), np.where(flip, 3 - c, c), np.where(flip, 3 - l, r)
    return np.where(ok, (l2 * 2 + (c2 == 3)) * 4 + r2, -1), flip


def dnms(ref_dir: str, dnm_dir: str, chain37: str, table: str, chain_t2t: str) -> list[dict]:
    """The new substitutions of both studies; the Platinum Pedigree's that Sasani et al. also called (the CEPH 1463
    family is in both) are counted once, as Sasani et al.'s, as in experiment 07."""
    import new_mutations as nm  # experiment 07

    snvs, _ = nm.sasani(dnm_dir, chain37)
    seen = {(d["chrom"], d["pos"], d["alt"]) for d in snvs}
    p, _ = nm.platinum(table, chain_t2t, ref_dir)
    return snvs + [d for d in p if (d["chrom"], d["pos"], d["alt"]) not in seen]


def main(ref_dir: str, hg_dir: str, bed_path: str, features: str, dnm_dir: str, chain37: str, table: str, chain_t2t: str,
         lcr_path: str) -> None:
    t0 = time.time()
    loci = read_loci(features, strata_rates())
    bed = read_bed(bed_path)
    lcr = read_bed(lcr_path)
    dn = dnms(ref_dir, dnm_dir, chain37, table, chain_t2t)
    dn_by = defaultdict(list)
    for d in dn:
        dn_by[d["chrom"]].append(d)
    sets = sorted({d["set"] for d in dn})
    print(f"{len(dn)} new substitutions in {sets} [{time.time() - t0:.0f}s]", flush=True)

    shp_a, shp_b = (NC, NS, NB, 32), (NC, NB, 32)
    a_all, a_call = np.zeros(shp_a, np.int64), np.zeros(shp_a, np.int64)
    a_snv, a_het, a_iso = (np.zeros(shp_a + (4,), np.int64) for _ in range(3))
    b_all, b_call = np.zeros(shp_b, np.int64), np.zeros(shp_b, np.int64)
    b_het, b_iso = np.zeros(shp_b + (4,), np.int64), np.zeros(shp_b + (4,), np.int64)
    cis_trans = np.zeros((NC, NB, 2, 2), np.int64)  # loci varying on one haplotype: [class, bin, isolated, cis]
    d_a = np.zeros((len(sets), NC, NB, 32, 4), np.int64)
    d_b = np.zeros((len(sets), NC, NB, 32, 4), np.int64)
    a_new, b_new = np.zeros(shp_b, np.int64), np.zeros(shp_b, np.int64)  # bases outside low-complexity regions
    seen = Counter()

    for c in ga.AUTOSOMES:
        s_, e_, cls_, _, _ = loci[c]
        per = [(s_[cls_ == k], e_[cls_ == k]) for k in range(NC)]
        bs, be = bed.get(c, (np.array([], int), np.array([], int)))
        ls, le = lcr.get(c, (np.array([], int), np.array([], int)))
        with ga.reference(ref_dir, c) as g, open_gpxil(os.path.join(hg_dir, f"HG002_{c}.gpxil.tif")) as h:
            ref, hc = g[c], h[c]
            n = hc.length
            # pass 1: where each haplotype varies, for the loci's states and for isolation
            v1, v2, unphased = [], [], 0
            for a in range(0, n, CHUNK):
                b = min(n, a + CHUNK)
                gt, fl = hc.read("gt", a, b), hc.read("flags", a, b)
                h1, h2 = gt >> 4, gt & 15
                v1.append(np.flatnonzero(((h1 != 0) & (h1 != 15)) | ((fl & 2) != 0)) + a)
                v2.append(np.flatnonzero(((h2 != 0) & (h2 != 15)) | ((fl & 4) != 0)) + a)
                unphased += int(((fl & 16) != 0).sum()) - int(((fl & 17) == 17).sum())
            v1, v2 = np.concatenate(v1), np.concatenate(v2)
            seen["variant bases not phased"] += unphased
            allv = np.union1d(v1, v2)
            state = np.zeros(len(s_), dtype=np.int8)
            for bit, v in ((1, v1), (2, v2)):
                li = np.searchsorted(s_, v, "right") - 1
                ok = (li >= 0) & (v < e_[np.maximum(li, 0)])
                state[np.unique(li[ok])] |= bit
            seen["loci where HG002 varies"] += int((state > 0).sum())
            seen["loci where HG002 varies on one haplotype"] += int(((state == 1) | (state == 2)).sum())
            # pass 2: the bases
            for a in range(0, n, CHUNK):
                b = min(n, a + CHUNK)
                xr = IDX[ref.read("ref", max(a - 1, 0), min(b + 1, n))]
                if a == 0:
                    xr = np.concatenate([[4], xr])
                if b == n:
                    xr = np.concatenate([xr, [4]])
                ctx, flip = contexts(xr)
                xc = xr[1:-1]
                base = (ctx >= 0) & ~interval_mask(a, b, s_, e_)
                q = np.flatnonzero(base) + a
                d, li = nearest(q, s_, e_)
                db = np.searchsorted(EDGES, d, "left") - 1
                k_a, s_a, cq = cls_[li], state[li], ctx[q - a]
                idx_a = ((k_a.astype(np.int64) * NS + s_a) * NB + db) * 32 + cq
                a_all += np.bincount(idx_a, minlength=a_all.size).reshape(shp_a)
                gt = hc.read("gt", a, b)
                h1, h2 = gt >> 4, gt & 15
                known = (h1 != 15) & (h2 != 15) & (h1 != 12) & (h2 != 12) & (h1 != 13) & (h2 != 13)
                cb = (interval_mask(a, b, bs, be) & known)[q - a]
                nl = ~interval_mask(a, b, ls, le)[q - a]
                a_new += np.bincount(((k_a.astype(np.int64) * NB + db) * 32 + cq)[nl], minlength=a_new.size).reshape(shp_b)
                a_call += np.bincount(idx_a[cb], minlength=a_all.size).reshape(shp_a)
                db_k = []
                for k in range(NC):
                    dk, _ = nearest(q, *per[k])
                    dbk = np.searchsorted(EDGES, dk, "left") - 1
                    db_k.append(dbk)
                    idx_b = (k * NB + dbk.astype(np.int64)) * 32 + cq
                    b_all += np.bincount(idx_b, minlength=b_all.size).reshape(shp_b)
                    b_call += np.bincount(idx_b[cb], minlength=b_all.size).reshape(shp_b)
                    b_new += np.bincount(idx_b[nl], minlength=b_all.size).reshape(shp_b)
                a1, a2 = HAP[h1[q - a]], HAP[h2[q - a]]
                fq, rq = flip[q - a], xc[q - a]
                for al in (a1, a2):
                    sub = cb & (al < 4) & (al != rq)
                    np.add.at(a_snv, (k_a[sub], s_a[sub], db[sub], cq[sub], np.where(fq, 3 - al, al)[sub]), 1)
                het1 = cb & (a1 < 4) & (a1 != rq) & (h2[q - a] == 0)
                het2 = cb & (a2 < 4) & (a2 != rq) & (h1[q - a] == 0)
                het = het1 | het2
                alh = np.where(het1, a1, a2)
                alt = np.where(fq, 3 - alh, alh)
                pos = q[het]
                kk = np.searchsorted(allv, pos)
                prev = np.where(kk > 0, allv[np.maximum(kk - 1, 0)], -10**9)
                nxt = np.where(kk + 1 < len(allv), allv[np.minimum(kk + 1, len(allv) - 1)], 10**12)
                iso = (pos - prev > ISO) & (nxt - pos > ISO)
                hh = (k_a[het], s_a[het], db[het], cq[het], alt[het])
                np.add.at(a_het, hh, 1)
                np.add.at(a_iso, tuple(v[iso] for v in hh), 1)
                for k in range(NC):
                    hb = (np.full(int(het.sum()), k), db_k[k][het], cq[het], alt[het])
                    np.add.at(b_het, hb, 1)
                    np.add.at(b_iso, tuple(v[iso] for v in hb), 1)
                one = (s_a[het] == 1) | (s_a[het] == 2)
                cis = np.where(het1[het], 1, 2) == s_a[het]
                np.add.at(cis_trans, (k_a[het][one], db[het][one], iso[one].astype(int), cis[one].astype(int)), 1)
                seen["bases"] += int(base.sum())
                seen["callable bases"] += int(cb.sum())
                seen["heterozygous substitutions"] += int(het.sum())
            for dd in dn_by.get(c, []):
                p = dd["pos"]
                if not 0 < p < n - 1:
                    continue
                xr = IDX[ref.read("ref", p - 1, p + 2)]
                ctx, flip = contexts(xr)
                if ctx[0] < 0 or "ACGT"[xr[1]] != dd["ref"]:
                    seen["new substitutions: reference differs, or N"] += 1
                    continue
                li0 = np.searchsorted(s_, p, "right") - 1
                if li0 >= 0 and p < e_[li0]:
                    seen["new substitutions: in a locus"] += 1
                    continue
                lj = np.searchsorted(ls, p, "right") - 1
                if lj >= 0 and p < le[lj]:
                    seen["new substitutions: in a low-complexity region"] += 1
                    continue
                d, li = nearest(np.array([p]), s_, e_)
                al = "ACGT".index(dd["alt"])
                alt = 3 - al if flip[0] else al
                si = sets.index(dd["set"])
                d_a[si, cls_[li[0]], int(np.searchsorted(EDGES, d[0], "left") - 1), ctx[0], alt] += 1
                for k in range(NC):
                    dk, _ = nearest(np.array([p]), *per[k])
                    d_b[si, k, int(np.searchsorted(EDGES, dk[0], "left") - 1), ctx[0], alt] += 1
                seen["new substitutions: placed"] += 1
        print(c, dict(seen), f"[{time.time() - t0:.0f}s]", flush=True)
    os.makedirs(RESULTS, exist_ok=True)
    np.savez_compressed(os.path.join(RESULTS, "shadow.npz"), A_all=a_all, A_call=a_call, A_snv=a_snv, A_het=a_het, A_iso=a_iso, A_new=a_new, B_new=b_new,
                        B_all=b_all, B_call=b_call, B_het=b_het, B_iso=b_iso, CT=cis_trans, D_A=d_a, D_B=d_b, sets=np.array(sets),
                        edges=EDGES, cls=CLS)
    with open(os.path.join(RESULTS, "shadow_counts.json"), "w") as fh:
        json.dump(dict(seen) | {"new substitutions read": len(dn), "sets": sets}, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    main(*sys.argv[1:10])
