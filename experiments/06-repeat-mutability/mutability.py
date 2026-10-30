"""Are pure GA-sequences more mutable than other tandem repeats as long, with units as short, and as pure?

Every autosomal locus of the tandem-repeat catalogue that TRGT genotyped in the Platinum Pedigree (Platinum TRs v1.0,
Dashnow 2024), measured on GRCh38's gpxil files: its length, its motifs' sizes, how periodic it is, what bases it is
made of, and how much of it lies in pure GA-sequences (experiment 04); and the new TR mutations that Porubsky et al.
2025 (Supplementary Table 10) called in eight G3 children and that pass their filters, by locus.

- **purity**: the share of the locus's bases that equal the base one unit on, for the best of its motifs' sizes (1 for
  a perfect repeat, about 0.25 for random sequence);
- **class**: a GA-sequence when at least half of the locus lies in GA-sequences; otherwise purine-rich (90 % or more of
  its bases A and G, or C and T), AT- or CG-rich, AC- or GT-rich (the same, for those pairs), or other;
- the loci are counted by length, the size of their shortest motif, purity and class, and those counts are all that the
  statistics use: rates per 1,000 loci a transmission; GA-sequences' mutations against those of other loci in the same
  strata (indirect standardisation), and a Poisson model of length, motif size, purity and class.

    python3 -I mutability.py CATALOG TABLE.xlsx REF_DIR [--features FILE]    -> results/strata.tsv.gz, results/mutability.json

CATALOG: human_GRCh38_no_alt_analysis_set.platinumTRs-v1.0.trgt.bed.gz (doi:10.5281/zenodo.13178746, CC BY 4.0).
TABLE.xlsx: Supplementary Table 10 of Porubsky et al. 2025. REF_DIR: the Gpxil browser's GRCh38 gpxil files.
--features also writes every locus's measures (7.3 million lines, 170 MB).
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "05-ga-de-novo"))

import ped  # noqa: E402  (puts experiment 04's ga module on the path)
import ga  # noqa: E402
from analyse import poisson_ci  # noqa: E402
from gpxiltools import codes as GC  # noqa: E402
from paper_calls import Runs, sheet  # noqa: E402

RESULTS = os.path.join(HERE, "results")
T = 8  #: transmissions: the eight G3 children
LENGTHS = [0, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 200, 250, 350, 500, 1000, 10**7]
MOTIFS = [0, 1, 2, 3, 4, 5, 6, 12, 10**5]
PURITY = [-1.0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.999, 1.0]
CLASSES = ("GA-sequence", "purine-rich", "AT- or CG-rich", "AC- or GT-rich", "other")
IDX = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate((GC.BASE_A, GC.BASE_C, GC.BASE_G, GC.BASE_T)):
    IDX[_b] = IDX[_b | GC.SOFTMASK] = _i


def bin_of(edges: list, x: float) -> int:
    """i where edges[i] < x <= edges[i + 1]"""
    lo, hi = 0, len(edges) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x > edges[mid]:
            lo = mid
        else:
            hi = mid
    return lo


def label(edges: list, i: int, digits: int = 0) -> str:
    a, b = edges[i], edges[i + 1]
    if digits:
        return f"{max(a, 0):.{digits}f}-{b:.{digits}f}"
    return f"{a + 1}-{b}" if b < edges[-1] else f"{a + 1}+"


def purity(seq: np.ndarray, sizes: list[int]) -> float:
    """The share of bases equal to the base one unit on, for the best of the motif sizes."""
    best = 0.0
    for k in sorted(set(sizes)):
        if k < len(seq):
            a, b = seq[:-k], seq[k:]
            good = (a < 4) & (b < 4)
            if good.any():
                best = max(best, float((a[good] == b[good]).mean()))
    return best


def class_of(n: dict, acgt: int, in_ga: bool) -> str:
    if in_ga:
        return CLASSES[0]
    for name, a, b in (("purine-rich", "R", "Y"), ("AT- or CG-rich", "W", "S"), ("AC- or GT-rich", "M", "K")):
        if max(n[a], n[b]) >= 0.9 * acgt:
            return name
    return CLASSES[-1]


def measure(catalog: str, table: str, ref_dir: str, features: str | None) -> tuple[dict, dict]:
    import openpyxl

    wb = openpyxl.load_workbook(table, read_only=True, data_only=True)
    muts = Counter(r["trid"] for r in sheet(wb, "TRs (GRCh38)", 1) if r["pass_all_filters"] == "YES" and r["chrom"] in ped.AUTOSOMES)
    runs = Runs()
    loci: dict = {c: [] for c in ped.AUTOSOMES}
    with gzip.open(catalog, "rt") as fh:
        for line in fh:
            c, s, e, info = line.rstrip("\n").split("\t")[:4]
            if c in loci:
                f = dict(x.split("=", 1) for x in info.split(";"))
                loci[c].append((int(s), int(e), f["ID"], [len(m) for m in f["MOTIFS"].split(",") if m]))
    strata: dict = defaultdict(lambda: [0, 0, 0])  # (length, motif, purity, class) -> loci, bases, mutations
    seen = Counter()
    fo = gzip.open(features, "wt") if features else None
    if fo:
        fo.write("chrom\tstart\tend\ttrid\tlength\tmotif_min\tpurity\tclass\tga_bases\tmutations\n")
    t0 = time.time()
    for c in ped.AUTOSOMES:
        with ga.reference(ref_dir, c) as g:
            x = IDX[g[c].read("ref")]
        rr = runs.runs.get(c, [])
        in_ga = ga.mask_of(len(x), [s for s, _, _ in rr], [e for _, e, _ in rr])
        lines = []
        for s, e, tid, sizes in loci[c]:
            seq = x[s:e]
            a, cc, gg, t, _ = np.bincount(seq, minlength=5).tolist()
            acgt = a + cc + gg + t
            if not acgt:
                seen["no bases"] += 1
                continue
            n = {"R": a + gg, "Y": cc + t, "W": a + t, "S": cc + gg, "M": a + cc, "K": gg + t}
            ga_bases = int(in_ga[s:e].sum())
            p = purity(seq, sizes)
            cls = class_of(n, acgt, ga_bases >= (e - s) / 2)
            m = muts.get(tid, 0)
            key = (bin_of(LENGTHS, e - s), bin_of(MOTIFS, min(sizes)), bin_of(PURITY, p), CLASSES.index(cls))
            v = strata[key]
            v[0] += 1
            v[1] += e - s
            v[2] += m
            seen["loci"] += 1
            seen["mutations"] += m
            if fo:
                lines.append(f"{c}\t{s}\t{e}\t{tid}\t{e - s}\t{min(sizes)}\t{p:.4f}\t{cls}\t{ga_bases}\t{m}")
        if fo:
            fo.write("\n".join(lines) + "\n")
        del x, in_ga
        print(f"{c}: {len(loci[c]):,} loci [{time.time() - t0:.0f}s]", flush=True)
    if fo:
        fo.close()
    seen["mutations_in_table"] = sum(muts.values())
    seen["loci_with_a_mutation"] = len(muts)
    return dict(strata), dict(seen)


def write_strata(strata: dict, path: str) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("length\tmotif\tpurity\tclass\tloci\tbases\tmutations\n")
        for (lb, kb, pb, cl), (n, b, m) in sorted(strata.items()):
            fh.write(f"{label(LENGTHS, lb)}\t{label(MOTIFS, kb)}\t{label(PURITY, pb, 3)}\t{CLASSES[cl]}\t{n}\t{b}\t{m}\n")


def read_strata(path: str):
    import pandas as pd

    df = pd.read_csv(path, sep="\t")
    lab = {name: [label(edges, i, d) for i in range(len(edges) - 1)] for name, edges, d in
           (("length", LENGTHS, 0), ("motif", MOTIFS, 0), ("purity", PURITY, 3))}
    for name, labels in lab.items():
        df[name] = pd.Categorical(df[name], categories=labels, ordered=True)
    df["ga"] = df["class"] == CLASSES[0]
    return df


def per1000(m: float, n: float) -> float | None:
    return round(1000 * m / n / T, 4) if n else None


def stats(df) -> dict:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf

    out: dict = {"loci": int(df.loci.sum()), "mutations": int(df.mutations.sum()), "transmissions": T}
    by = df.groupby("class").agg(loci=("loci", "sum"), bases=("bases", "sum"), mutations=("mutations", "sum"))
    out["by_class"] = {c: {"loci": int(r.loci), "bases": int(r.bases), "mutations": int(r.mutations),
                           "per_1000_loci_per_transmission": per1000(r.mutations, r.loci)} for c, r in by.iterrows()}
    # indirect standardisation: GA-sequences against the other loci of their strata
    std = {}
    for keys in (["length"], ["length", "motif"], ["length", "motif", "purity"]):
        ref = df[~df.ga].groupby(keys, observed=True).agg(n=("loci", "sum"), m=("mutations", "sum"))
        ref["rate"] = ref.m / ref.n
        g = df[df.ga].groupby(keys, observed=True).agg(n=("loci", "sum"), m=("mutations", "sum")).join(ref["rate"])
        ok = g.rate.notna() & (g.n > 0)
        obs, exp = int(g[ok].m.sum()), float((g[ok].n * g[ok].rate).sum())
        lo, hi = poisson_ci(obs)
        std[" x ".join(keys)] = {"ga_loci": int(g[ok].n.sum()), "of": int(g.n.sum()), "observed": obs, "expected": round(exp, 1),
                                 "ratio": round(obs / exp, 2), "ci": [round(lo / exp, 2), round(hi / exp, 2)]}
    out["ga_against_other_loci_of_their_strata"] = std
    # Poisson models on the strata counts (the same likelihood as per locus, for these categorical terms)
    d = df[df.loci > 0].copy()
    d["cls"] = d["class"]
    for col in ("length", "motif", "purity"):
        d[col] = d[col].cat.remove_unused_categories()
    models = {}
    for name, formula in (("ga", "mutations ~ C(length) + C(motif) + C(purity) + ga"),
                          ("class", "mutations ~ C(length) + C(motif) + C(purity) + C(cls, Treatment('other'))")):
        d["ga"] = d["class"].eq(CLASSES[0]).astype(int)
        fit = smf.glm(formula, data=d, family=sm.families.Poisson(), offset=np.log(d.loci * T)).fit()
        terms = {}
        for p in fit.params.index:
            if p == "ga" or p.startswith("C(cls"):
                b, se = fit.params[p], fit.bse[p]
                terms[p.replace("C(cls, Treatment('other'))[T.", "").rstrip("]")] = {
                    "rate_ratio": round(float(np.exp(b)), 2), "ci": [round(float(np.exp(b - 1.96 * se)), 2), round(float(np.exp(b + 1.96 * se)), 2)],
                    "p": float(f"{fit.pvalues[p]:.2g}")}
        models[name] = {"formula": formula, "terms": terms, "deviance": round(float(fit.deviance), 1), "df_resid": int(fit.df_resid)}
    out["poisson_models"] = models
    # rates by purity within length bands, GA-sequences against the rest (motifs of 6 bases or fewer)
    bands = {"51-100": ["51-60", "61-75", "76-100"], "101-250": ["101-125", "126-150", "151-200", "201-250"],
             "251+": ["251-350", "351-500", "501-1000", "1001+"]}
    small = df[df.motif.isin(["1-1", "2-2", "3-3", "4-4", "5-5", "6-6"])]
    by_purity = {}
    for band, labs in bands.items():
        sub = small[small.length.isin(labs)]
        rows = []
        for p, grp in sub.groupby("purity", observed=True):
            ga_, ot = grp[grp.ga], grp[~grp.ga]
            row = {"purity": str(p)}
            for tag, x in (("ga", ga_), ("other", ot)):
                n, m = int(x.loci.sum()), int(x.mutations.sum())
                lo, hi = poisson_ci(m)
                row[tag] = {"loci": n, "mutations": m, "per_1000": per1000(m, n),
                            "ci": [per1000(lo, n), per1000(hi, n)] if n else None}
            rows.append(row)
        by_purity[band] = rows
    out["by_purity_in_length_bands"] = by_purity
    # the long, pure repeats: what they are made of
    long_pure = {}
    for min_len, min_pur in ((101, 0.85), (101, 0.8), (101, 0.9), (76, 0.85), (151, 0.85)):
        lab_len = [label(LENGTHS, i) for i in range(len(LENGTHS) - 1) if LENGTHS[i] + 1 >= min_len]
        lab_pur = [label(PURITY, i, 3) for i in range(len(PURITY) - 1) if PURITY[i] >= min_pur - 1e-9]
        sub = small[small.length.isin(lab_len) & small.purity.isin(lab_pur)]
        tot_n, tot_m = int(sub.loci.sum()), int(sub.mutations.sum())
        rows = {}
        for c, grp in sub.groupby("class"):
            n, m = int(grp.loci.sum()), int(grp.mutations.sum())
            rows[c] = {"loci": n, "share_of_loci": round(n / tot_n, 3), "mutations": m,
                       "share_of_mutations": round(m / tot_m, 3) if tot_m else None, "per_1000": per1000(m, n)}
        long_pure[f"length >= {min_len}, motif <= 6, purity >= {min_pur}"] = {"loci": tot_n, "mutations": tot_m, "by_class": rows}
    out["long_pure_repeats"] = long_pure
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("catalog")
    ap.add_argument("table")
    ap.add_argument("ref_dir")
    ap.add_argument("--features", help="also write every locus's measures here (.tsv.gz)")
    ap.add_argument("--strata", help="skip the measuring: read these strata counts instead")
    a = ap.parse_args()
    os.makedirs(RESULTS, exist_ok=True)
    path = os.path.join(RESULTS, "strata.tsv.gz")
    seen = {}
    if not a.strata:
        strata, seen = measure(a.catalog, a.table, a.ref_dir, a.features)
        write_strata(strata, path)
    out = {"source": {"catalog": os.path.basename(a.catalog), "table": "Porubsky et al. 2025, Supplementary Table 10"}, "counted": seen}
    out.update(stats(read_strata(a.strata or path)))
    with open(os.path.join(RESULTS, "mutability.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in ("counted", "by_class", "ga_against_other_loci_of_their_strata", "poisson_models")}, indent=1))


if __name__ == "__main__":
    main()
