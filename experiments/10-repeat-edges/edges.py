"""The bases right beside a repeat: do their substitutions extend it?

For each catalogue locus, the first bases beside it on either side, up to its shortest motif's length p and at most 4
(base 1 is the one touching it): the base that would continue the repeat there is the one p bases back inside it. A
substitution extends the repeat when its new base is that one; a base whose reference already continues the repeat is
counted apart. Bases inside another locus, N, or (for HG002) outside the benchmark's callable bases are left out.

* HG002's heterozygous substitutions there (v5.0q), and those isolated (no other variant within 10 bases), against the
  spectrum of its heterozygous substitutions at the far bases (shadow.npz): all substitutions, those that extend and
  those that don't, by the locus's instability class;
* the new substitutions of the CEPH families and the Platinum Pedigree, the same way, against their own spectrum 201 to
  1,000 bases from any locus, bases and substitutions both outside low-complexity regions (LCR-hs38), where Sasani et
  al. called none;
* the bases beside pure two-letter runs (experiment 09's census): a substitution extends the run when its new base is
  in the run's alphabet.

    python3 -I edges.py REF_DIR HG002_DIR BED FEATURES DNM_DIR CHAIN37 TABLE.xlsx CHAIN_T2T LCR     -> results/edges.json
"""

from __future__ import annotations

import gzip
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import shadow as SH  # noqa: E402  (this experiment's counts and helpers)
from distances import poisson_ci  # noqa: E402
import ga  # noqa: E402  (experiment 04, on the path through shadow)
from gpxiltools.container import open_gpxil  # noqa: E402

RUNS = os.path.join(SH.EXP, "09-alphabets", "results", "runs.tsv.gz")
DMAX = 4
FAR = slice(7, 10)
LETTERS = {"A/G": "AG", "C/T": "CT", "A/T": "AT", "G/C": "GC", "A/C": "AC", "G/T": "GT"}


def ratio(o: float, e: float) -> dict:
    lo, hi = poisson_ci(int(round(o)))
    return {"observed": int(round(o)), "expected": round(float(e), 2), "ratio": round(o / e, 4) if e else None,
            "ci": [round(lo / e, 4), round(hi / e, 4)] if e else None}


def ctx_of(x: np.ndarray, pos: np.ndarray):
    l, m, r = x[pos - 1].astype(np.int16), x[pos].astype(np.int16), x[pos + 1].astype(np.int16)
    good = (l < 4) & (m < 4) & (r < 4)
    flip = (m == 0) | (m == 2)
    ctx = (np.where(flip, 3 - r, l) * 2 + (np.where(flip, 3 - m, m) == 3)) * 4 + np.where(flip, 3 - l, r)
    return np.where(good, ctx, 0), flip, good, m


def summarise(t: np.ndarray) -> dict:
    """t: bases, expected substitutions, expected extending ones, observed, observed extending, the observed ones'
    summed chance of extending (given a substitution there)."""
    n, e, e_ext, o, o_ext, p_ext = t
    return {"bases": int(n), "all": ratio(o, e), "extending": ratio(o_ext, e_ext), "others": ratio(o - o_ext, e - e_ext),
            "share extending": round(o_ext / o, 4) if o else None, "share extending by the spectrum, where they are": round(p_ext / o, 4) if o else None,
            "q, the spectrum's chance of extending at these bases": round(e_ext / e, 4) if e else None}


def main(ref_dir, hg_dir, bed_path, features, dnm_dir, chain37, table, chain_t2t, lcr_path) -> None:
    z = np.load(os.path.join(HERE, "results", "shadow.npz"))
    base_far = z["A_call"][:, :, FAR, :].sum(axis=(0, 1, 2))
    rate = z["A_het"][:, :, FAR].sum(axis=(0, 1, 2)) / np.maximum(base_far, 1)[:, None]  # HG002: [ctx, new base]
    rate_iso = z["A_iso"][:, :, FAR].sum(axis=(0, 1, 2)) / np.maximum(base_far, 1)[:, None]
    a_all = z["A_new"]  # (bases outside low-complexity regions)
    far_new = slice(7, 9)  # (201 to 1,000 bases: see distances.py)
    mu = z["D_A"].sum(axis=0)[:, far_new].sum(axis=(0, 1)) / np.maximum(a_all[:, far_new].sum(axis=(0, 1)), 1)[:, None]  # new substitutions
    sets = [str(x) for x in z["sets"]]
    plat = [i for i, x in enumerate(sets) if x.startswith("Platinum")]
    share = {"Platinum": float(z["D_A"][plat][:, :, far_new].sum() / z["D_A"][:, :, far_new].sum())}
    share["Sasani"] = 1 - share["Platinum"]
    loci = SH.read_loci(features, SH.strata_rates())
    bed = SH.read_bed(bed_path)
    lcr = SH.read_bed(lcr_path)
    dn_by = defaultdict(list)
    for d in SH.dnms(ref_dir, dnm_dir, chain37, table, chain_t2t):
        dn_by[d["chrom"]].append(d)
    runs = defaultdict(list)
    with gzip.open(RUNS, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["chrom"] in ga.AUTOSOMES:
                runs[r["chrom"]].append((int(r["start"]), int(r["end"]), r["alphabet"], LETTERS[r["kind"]], r["origin"].endswith(" tail")))
    hg = defaultdict(lambda: np.zeros(6))  # (class, base, reference continues, isolated)
    new = defaultdict(lambda: np.zeros(6))  # (class, base)
    run_hg = defaultdict(lambda: np.zeros(6))  # (alphabet, tail, isolated)
    run_new = defaultdict(lambda: np.zeros(6))  # alphabet
    by_study = defaultdict(lambda: np.zeros(4))  # (study, class) at base 1: expected, expected extending, observed, extending
    for c in ga.AUTOSOMES:
        s_, e_, cls_, _, p_ = loci[c]
        bs, be = bed[c]
        ls, le = lcr.get(c, (np.array([], int), np.array([], int)))

        def in_lcr(pos: np.ndarray) -> np.ndarray:
            j = np.searchsorted(ls, pos, "right") - 1
            return (j >= 0) & (pos < le[np.maximum(j, 0)]) if len(ls) else np.zeros(len(pos), dtype=bool)

        with ga.reference(ref_dir, c) as g:
            x = SH.IDX[g[c].read("ref")]
        with open_gpxil(os.path.join(hg_dir, f"HG002_{c}.gpxil.tif")) as h:
            gt = h[c].read("gt")
            fl = h[c].read("flags")
        n = len(x)
        h1, h2 = gt >> 4, gt & 15
        del gt
        anyv = np.flatnonzero(((h1 != 0) & (h1 != 15)) | ((h2 != 0) & (h2 != 15)) | ((fl & 6) != 0))
        del fl
        dn_pos = defaultdict(list)
        for d in dn_by.get(c, []):
            dn_pos[d["pos"]].append(d)
        dn_arr = np.array(sorted(dn_pos), dtype=np.int64)

        def hg002(pos: np.ndarray):
            """callable, and the heterozygous substitution there (new base, or 4), and whether isolated"""
            bi = np.searchsorted(bs, pos, "right") - 1
            call = (bi >= 0) & (pos < be[np.maximum(bi, 0)])
            g1, g2 = h1[pos], h2[pos]
            call &= (g1 != 15) & (g2 != 15) & (g1 != 12) & (g2 != 12) & (g1 != 13) & (g2 != 13)
            a1, a2, m = SH.HAP[g1], SH.HAP[g2], x[pos]
            het1 = (a1 < 4) & (a1 != m) & (g2 == 0)
            het2 = (a2 < 4) & (a2 != m) & (g1 == 0)
            alt = np.where(het1, a1, np.where(het2, a2, 4))
            k = np.searchsorted(anyv, pos)
            prev = np.where(k > 0, anyv[np.maximum(k - 1, 0)], -10**9)
            nxt = np.where(k + 1 < len(anyv), anyv[np.minimum(k + 1, len(anyv) - 1)], 10**12)
            return call, alt, (pos - prev > SH.ISO) & (nxt - pos > SH.ISO)

        for side in ("left", "right"):
            for d in range(1, DMAX + 1):
                ok = (p_ >= d) & ((e_ - s_) >= p_)
                pos = (s_ - d) if side == "left" else (e_ - 1 + d)
                cont = (pos + p_) if side == "left" else (pos - p_)
                ok &= (pos >= 1) & (pos < n - 1)
                pos, cont, cls = pos[ok], cont[ok], cls_[ok]
                li = np.searchsorted(s_, pos, "right") - 1
                inside = (li >= 0) & (pos < e_[np.maximum(li, 0)])
                ctx, flip, good, m = ctx_of(x, pos)
                cb = x[cont].astype(np.int16)
                keep = ~inside & good & (cb < 4)
                refc = m == cb
                cbc = np.where(flip, 3 - np.minimum(cb, 3), np.minimum(cb, 3))
                call, alt, iso = hg002(pos)
                for iso_only, rt in ((0, rate), (1, rate_iso)):
                    r_all, r_ext = rt.sum(axis=1)[ctx], rt[ctx, cbc] * ~refc
                    p_ext = r_ext / np.maximum(r_all, 1e-12)
                    sub = call & (alt < 4) & (iso if iso_only else True)
                    ext = sub & (alt == cb)
                    for k in range(SH.NC):
                        for rc in (0, 1):
                            sel = keep & call & (cls == k) & (refc == bool(rc))
                            t = hg[(k, d, rc, iso_only)]
                            t += [sel.sum(), r_all[sel].sum(), r_ext[sel].sum(), (sub & sel).sum(), (ext & sel).sum(), p_ext[sub & sel].sum()]
                # new substitutions, where the reference breaks the repeat, outside low-complexity regions
                r_all, r_ext = mu.sum(axis=1)[ctx], mu[ctx, cbc]
                open_ = keep & ~refc & ~in_lcr(pos)
                for k in range(SH.NC):
                    sel = open_ & (cls == k)
                    t = new[(k, d)]
                    t[0] += sel.sum()
                    t[1] += r_all[sel].sum()
                    t[2] += r_ext[sel].sum()
                    if d == 1:
                        for study in ("Sasani", "Platinum"):
                            by_study[(study, k)][:2] += [r_all[sel].sum() * share[study], r_ext[sel].sum() * share[study]]
                for i in np.flatnonzero(open_ & np.isin(pos, dn_arr)):
                    for dd in dn_pos[int(pos[i])]:
                        if "ACGT"[m[i]] != dd["ref"]:
                            continue
                        t = new[(int(cls[i]), d)]
                        t[3] += 1
                        t[4] += "ACGT".index(dd["alt"]) == cb[i]
                        t[5] += r_ext[i] / r_all[i]
                        if d == 1:
                            by_study[(dd["set"].split()[0], int(cls[i]))][2:] += [1, "ACGT".index(dd["alt"]) == cb[i]]
        # the pure runs' edges
        rr = runs.get(c, [])
        if rr:
            pos = np.array([p for s, e, *_ in rr for p in (s - 1, e)])
            info = [(al, letters, tail) for s, e, al, letters, tail in rr for _ in (0, 1)]
            ok = (pos >= 1) & (pos < n - 1)
            lcr_run = in_lcr(pos)
            ctx, flip, good, m = ctx_of(x, np.clip(pos, 1, n - 2))
            call, alt, iso = hg002(np.clip(pos, 1, n - 2))
            for i in np.flatnonzero(ok & good):
                al, letters, tail = info[i]
                ext_b = [b for b in ("ACGT".index(ch) for ch in letters) if b != m[i]]
                ext_c = [3 - b if flip[i] else b for b in ext_b]
                for iso_only, rt in ((0, rate), (1, rate_iso)):
                    if not call[i]:
                        continue
                    t = run_hg[(al, tail, iso_only)]
                    r_all, r_ext = rt[ctx[i]].sum(), rt[ctx[i], ext_c].sum()
                    t[:3] += [1, r_all, r_ext]
                    if alt[i] < 4 and (iso[i] or not iso_only):
                        t[3] += 1
                        t[4] += alt[i] in ext_b
                        t[5] += r_ext / r_all
                if lcr_run[i]:
                    continue
                t = run_new[al]
                r_all, r_ext = mu[ctx[i]].sum(), mu[ctx[i], ext_c].sum()
                t[:3] += [1, r_all, r_ext]
                for dd in dn_pos.get(int(pos[i]), []):
                    if "ACGT"[m[i]] == dd["ref"]:
                        t[3] += 1
                        t[4] += "ACGT".index(dd["alt"]) in ext_b
                        t[5] += r_ext / r_all
        print(c, flush=True)
        del x, h1, h2
    classes = ["below 1e-6", "1e-6 to 1e-5", "1e-5 to 1e-4", "1e-4 or more"]
    out = {"loci, HG002": {}, "loci, new substitutions": {}, "runs, HG002": {}, "runs, new substitutions": {}}
    for (k, d, rc, iso_only), t in sorted(hg.items()):
        key = f"{classes[k]}, base {d}, reference {'continues' if rc else 'breaks'} the repeat{', isolated' if iso_only else ''}"
        out["loci, HG002"][key] = summarise(t)
    for (k, d), t in sorted(new.items()):
        out["loci, new substitutions"][f"{classes[k]}, base {d}"] = summarise(t)
    tot = sum(new[(k, 1)] for k in (1, 2, 3))
    out["loci, new substitutions"]["1e-6 or more, base 1"] = summarise(tot)
    for (al, tail, iso_only), t in sorted(run_hg.items()):
        out["runs, HG002"][f"{al}, {'tails' if tail else 'others'}{', isolated' if iso_only else ''}"] = summarise(t)
    for al in sorted({k[0] for k in run_hg}):
        for iso_only in (0, 1):
            t = run_hg[(al, True, iso_only)] + run_hg[(al, False, iso_only)]
            out["runs, HG002"][f"{al}{', isolated' if iso_only else ''}"] = summarise(t)
    for al, t in sorted(run_new.items()):
        out["runs, new substitutions"][al] = summarise(t)
    out["loci, new substitutions, by study, base 1"] = {
        f"{study}, {classes[k] if k is not None else 'all classes'}": {
            "expected": round(float(v[0]), 2), "observed": int(v[2]), "expected extending": round(float(v[1]), 2), "extending": int(v[3])}
        for study in ("Sasani", "Platinum") for k in list(range(SH.NC)) + [None]
        for v in [by_study[(study, k)] if k is not None else sum(by_study[(study, kk)] for kk in range(SH.NC))]}
    out["new substitutions' share by study, 201 to 1,000 bases from any locus"] = {k: round(v, 4) for k, v in share.items()}
    with open(os.path.join(HERE, "results", "edges.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for sec, v in out.items():
        print(sec)
        for k, r in v.items():
            if "base 1" in k or sec.startswith("runs"):
                print(f"  {k}: all {r['all']['ratio']} ({r['all']['observed']}), extending {r['extending']['ratio']} ({r['extending']['observed']} vs {r['extending']['expected']}), others {r['others']['ratio']}")


if __name__ == "__main__":
    if len(sys.argv) != 10:
        sys.exit(__doc__)
    main(*sys.argv[1:10])
