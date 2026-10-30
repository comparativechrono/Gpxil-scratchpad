"""Do repeats whose two strands differ mutate differently depending on which strand a gene's transcript carries?

A repeat whose strands differ, a poly(A) or a GA-sequence, sits in a gene one way or the other: its transcript carries
the A-rich (or purine-rich, or G-rich) strand, and the polymerase reads the other, or the reverse. Every autosomal
locus of the catalogue (Platinum TRs v1.0; experiment 06's base counts on the + strand) is put in a group by how its
strands differ:

- poly(A): one strand 90 % A or more (the other T);
- purine: one strand 90 % purines or more, not a poly(A) (GA-sequences, (AAAG)n, (AG)n ...);
- A-rich A/T: 90 % A or T, one strand at least three-quarters A, not a poly(A) ((AAAT)n ...);
- G-rich: one strand at least 60 % G and at most 10 % C, not in the groups above;
- and, as a control without a strand, A/T alternating: 90 % A or T, neither strand more than 60 % A ((AT)n ...);

and by where it lies: in the introns of genes on one strand only (Ensembl 116's protein-coding and lncRNA genes; loci
within 100 bases of an exon of them, where splice sites' signals are, or in genes on both strands, are left out),
sense (the transcript carries the A-, purine- or G-rich strand) or antisense; or between genes. For the control,
sense means a gene on the + strand.

Three measures of how much a locus mutates:

- the family: whether either parent of the Platinum Pedigree's G3 children (NA12877, NA12878) has a haplotype whose
  allele differs in length from the reference (TRGT's open-access calls; loci they don't list: none does);
- HG002: whether either haplotype has an insertion or deletion in the locus, among loci GIAB's v5.0q calls on both
  haplotypes (with the base either side), read from its gpxil files;
- new TR mutations: the Platinum Pedigree's, in eight G3 children (experiment 06).

Sense is compared with antisense within experiment 06's strata of length, motif size and purity: a Mantel-Haenszel odds
ratio for the family and HG002, a Mantel-Haenszel rate ratio for new mutations; and each with the group's loci between
genes. Alus and L1s in introns lean to one orientation, so the comparison is made again within the loci at Alu and
L1 tails (experiment 11's) and the others.

    python3 -I strands.py GTF CATALOG_FEATURES TRGT_VCF REF_DIR HG002_DIR     -> results/strands.json, results/loci.tsv.gz

GTF: Ensembl's Homo_sapiens.GRCh38.116.gtf.gz. CATALOG_FEATURES: experiment 06's measures with base counts (R, Y, W,
S, M, K, acgt). TRGT_VCF: ceph_1463_tandem_repeats.oa.vcf.gz. REF_DIR, HG002_DIR: as in experiment 10.
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
sys.path.insert(0, os.path.join(EXP, "04-ga-sequences"))
sys.path.insert(0, os.path.join(EXP, "06-repeat-mutability"))

import ga  # noqa: E402  (experiment 04)
import mutability as M06  # noqa: E402  (experiment 06: the strata)
from gpxiltools import codes as GC  # noqa: E402
from gpxiltools.container import open_gpxil  # noqa: E402

RESULTS = os.path.join(HERE, "results")
TAILS = os.path.join(EXP, "11-along-arrays", "results", "tail_loci.tsv.gz")  #: the loci at Alu and L1 tails (experiment 11)
BIOTYPES = ("protein_coding", "lncRNA")
PARENTS = ("NA12877", "NA12878")
EXON_MARGIN = 100  #: loci this close to an exon are left out (splice sites' polypyrimidine tracts are under selection)
GROUPS = ("poly(A)", "purine", "A-rich A/T", "G-rich", "A/T alternating")
WHERE = ("sense", "antisense", "between genes")


def genes(gtf: str) -> tuple[dict, dict]:
    """chrom -> strand -> merged gene bodies (starts, ends), 0-based; chrom -> merged exons of those genes."""
    body = defaultdict(lambda: defaultdict(list))
    exon = defaultdict(list)
    keep = set()
    with gzip.open(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t", 9)
            if f[2] not in ("gene", "exon"):
                continue
            c = "chr" + f[0]
            if c not in ga.AUTOSOMES:
                continue
            attrs = f[8]
            gid = attrs.split('gene_id "', 1)[1].split('"', 1)[0]
            if f[2] == "gene":
                bt = attrs.split('gene_biotype "', 1)[1].split('"', 1)[0] if 'gene_biotype "' in attrs else ""
                if bt in BIOTYPES:
                    keep.add(gid)
                    body[c][f[6]].append((int(f[3]) - 1, int(f[4])))
            elif gid in keep:  # (a gene line comes before its exons in Ensembl's GTF)
                exon[c].append((int(f[3]) - 1, int(f[4])))

    def merged(v):
        v.sort()
        s, e = [], []
        for a, b in v:
            if s and a <= e[-1]:
                e[-1] = max(e[-1], b)
            else:
                s.append(a)
                e.append(b)
        return np.array(s, dtype=np.int64), np.array(e, dtype=np.int64)

    return ({c: {s: merged(v) for s, v in d.items()} for c, d in body.items()}, {c: merged(v) for c, v in exon.items()})


def overlaps(s: np.ndarray, e: np.ndarray, iv: tuple) -> np.ndarray:
    """Whether each [s, e) overlaps the merged intervals iv."""
    if iv is None or not len(iv[0]):
        return np.zeros(len(s), dtype=bool)
    st, en = iv
    i = np.searchsorted(en, s, "right")  # the first interval ending after s
    return (i < len(st)) & (st[np.minimum(i, len(st) - 1)] < e)


def family_polymorphic(vcf: str) -> set:
    """The loci where a parent's haplotype differs in length from the reference."""
    out = set()
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = f[9:]
                ip = [samples.index(p) for p in PARENTS]
                continue
            alle = [f[3]] + (f[4].split(",") if f[4] != "." else [])
            n = len(f[3])
            for j in ip:
                if any(g not in (".", "0") and len(alle[int(g)]) != n for g in f[9 + j].split(":")[0].replace("|", "/").split("/")):
                    out.add(f[7].split("TRID=")[1].split(";")[0])
                    break
    return out


def hg002_indels(ref_dir: str, hg_dir: str, loci) -> tuple[np.ndarray, np.ndarray]:
    """Per locus: called on both haplotypes over [start - 1, end + 1); an insertion or deletion inside."""
    called = np.zeros(len(loci), dtype=bool)
    indel = np.zeros(len(loci), dtype=bool)
    chrom = loci.chrom.to_numpy()
    t0 = time.time()
    for c in ga.AUTOSOMES:
        idx = np.flatnonzero(chrom == c)
        if not len(idx):
            continue
        with open_gpxil(os.path.join(hg_dir, f"HG002_{c}.gpxil.tif"), reference=os.path.join(ref_dir, f"{c}.gpxil.tif")) as h:
            hc = h.find(c)
            n = hc.length
            gt = hc.read("gt", 0, n)
            h1, h2 = GC.hap1(gt), GC.hap2(gt)
            del gt
            unk = np.flatnonzero((h1 == GC.HAP_UNKNOWN) | (h2 == GC.HAP_UNKNOWN))
            dels = (h1 == GC.HAP_DEL) | (h2 == GC.HAP_DEL)
            del h1, h2
            fl = hc.read("flags", 0, n)
            ins = np.flatnonzero((fl & (GC.FLAG_INS1 | GC.FLAG_INS2)) != 0)
            dl = np.flatnonzero(dels)
            del fl, dels
        s = loci.start.to_numpy()[idx]
        e = loci.end.to_numpy()[idx]
        a, b = np.maximum(s - 1, 0), np.minimum(e + 1, n)
        called[idx] = np.searchsorted(unk, a, "left") == np.searchsorted(unk, b, "left")
        # a deleted base inside the locus; an insertion after one of its bases or after the base before it
        indel[idx] = ((np.searchsorted(dl, s, "left") < np.searchsorted(dl, e, "left"))
                      | (np.searchsorted(ins, a, "left") < np.searchsorted(ins, e, "left")))
        print(c, int(called[idx].sum()), int((indel[idx] & called[idx]).sum()), f"[{time.time() - t0:.0f}s]", flush=True)
    return called, indel


def mh_odds(y: np.ndarray, x: np.ndarray, strata: np.ndarray) -> dict:
    """Mantel-Haenszel odds ratio of y (bool) for x (bool: True the exposed group) over strata, with its 95 % interval
    (Robins, Breslow and Greenland)."""
    k = strata.max() + 1
    a = np.bincount(strata[x & y], minlength=k).astype(float)
    b = np.bincount(strata[x & ~y], minlength=k).astype(float)
    c = np.bincount(strata[~x & y], minlength=k).astype(float)
    d = np.bincount(strata[~x & ~y], minlength=k).astype(float)
    n = a + b + c + d
    ok = n > 1
    a, b, c, d, n = a[ok], b[ok], c[ok], d[ok], n[ok]
    R, S = a * d / n, b * c / n
    if R.sum() == 0 or S.sum() == 0:
        return {"odds_ratio": None}
    P, Q = (a + d) / n, (b + c) / n
    v = ((P * R).sum() / (2 * R.sum() ** 2) + ((P * S + Q * R).sum()) / (2 * R.sum() * S.sum()) + (Q * S).sum() / (2 * S.sum() ** 2))
    orr = R.sum() / S.sum()
    se = np.sqrt(v)
    return {"odds_ratio": round(float(orr), 3), "ci": [round(float(orr * np.exp(-1.96 * se)), 3), round(float(orr * np.exp(1.96 * se)), 3)],
            "z": round(float(np.log(orr) / se), 2)}


def mh_rate(m: np.ndarray, x: np.ndarray, strata: np.ndarray) -> dict:
    """Mantel-Haenszel rate ratio of counts m for x over strata (each locus the same exposure), with a 95 % interval
    from the variance of the exposed group's count given each stratum's total (binomial, as Greenland and Robins)."""
    k = strata.max() + 1
    m1 = np.bincount(strata[x], weights=m[x], minlength=k)
    m0 = np.bincount(strata[~x], weights=m[~x], minlength=k)
    n1 = np.bincount(strata[x], minlength=k).astype(float)
    n0 = np.bincount(strata[~x], minlength=k).astype(float)
    n = n1 + n0
    ok = (n1 > 0) & (n0 > 0)
    m1, m0, n1, n0, n = m1[ok], m0[ok], n1[ok], n0[ok], n[ok]
    num, den = (m1 * n0 / n).sum(), (m0 * n1 / n).sum()
    out = {"exposed": int(m1.sum()), "other": int(m0.sum())}
    if num == 0 or den == 0:
        return {**out, "rate_ratio": None}
    rr = num / den
    v = ((n1 * n0 * (m1 + m0)) / n ** 2).sum() / (num * den)
    se = np.sqrt(v)
    return {**out, "rate_ratio": round(float(rr), 3),
            "ci": [round(float(rr * np.exp(-1.96 * se)), 3), round(float(rr * np.exp(1.96 * se)), 3)], "z": round(float(np.log(rr) / se), 2)}


def main(gtf: str, features: str, vcf: str, ref_dir: str, hg_dir: str) -> None:
    import pandas as pd
    cols = ["chrom", "start", "end", "trid", "motif_min", "share", "R", "Y", "W", "S", "M", "K", "acgt", "mutations"]
    dt = {k: "int32" for k in cols if k not in ("chrom", "trid", "share", "start", "end")}
    dt.update({"chrom": "category", "share": "float32", "start": "int64", "end": "int64"})
    f = pd.read_csv(features, sep="\t", usecols=cols, dtype=dt)
    f = f[f.chrom.isin(ga.AUTOSOMES)].reset_index(drop=True)
    acgt = f.acgt.to_numpy().astype(float)
    A = (f.R + f.W + f.M - f.acgt).to_numpy() / 2
    C, G, T = f.M.to_numpy() - A, f.R.to_numpy() - A, f.W.to_numpy() - A
    R, Y, W = f.R.to_numpy(), f.Y.to_numpy(), f.W.to_numpy()
    with np.errstate(divide="ignore", invalid="ignore"):
        poly_a = (A >= 0.9 * acgt) | (T >= 0.9 * acgt)
        purine = ~poly_a & ((R >= 0.9 * acgt) | (Y >= 0.9 * acgt))
        a_rich = ~poly_a & ~purine & (W >= 0.9 * acgt) & ((A >= 0.75 * W) | (T >= 0.75 * W))
        g_rich = ~poly_a & ~purine & ~a_rich & (((G >= 0.6 * acgt) & (C <= 0.1 * acgt)) | ((C >= 0.6 * acgt) & (G <= 0.1 * acgt)))
        at_alt = (W >= 0.9 * acgt) & (A <= 0.6 * W) & (T <= 0.6 * W)
    group = np.full(len(f), "", dtype=object)
    rich_plus = np.zeros(len(f), dtype=bool)  # the rich strand is the + strand
    for name, m, plus in (("poly(A)", poly_a, A >= T), ("purine", purine, R >= Y), ("A-rich A/T", a_rich, A >= T),
                          ("G-rich", g_rich, G >= C), ("A/T alternating", at_alt, np.ones(len(f), dtype=bool))):
        group[m] = name
        rich_plus[m] = plus[m]
    print({g: int((group == g).sum()) for g in GROUPS}, flush=True)
    body, exons = genes(gtf)
    in_p = np.zeros(len(f), dtype=bool)
    in_m = np.zeros(len(f), dtype=bool)
    exonic = np.zeros(len(f), dtype=bool)
    chrom = f.chrom.to_numpy()
    for c in ga.AUTOSOMES:
        idx = np.flatnonzero(chrom == c)
        s, e = f.start.to_numpy()[idx], f.end.to_numpy()[idx]
        in_p[idx] = overlaps(s, e, body.get(c, {}).get("+"))
        in_m[idx] = overlaps(s, e, body.get(c, {}).get("-"))
        exonic[idx] = overlaps(s - EXON_MARGIN, e + EXON_MARGIN, exons.get(c))
    where = np.full(len(f), "", dtype=object)
    one = in_p ^ in_m
    gene_plus = in_p & ~in_m
    intronic = one & ~exonic
    where[intronic & (gene_plus == rich_plus)] = "sense"
    where[intronic & (gene_plus != rich_plus)] = "antisense"
    where[~in_p & ~in_m] = "between genes"
    print("where", Counter(where.tolist()), flush=True)
    poly = family_polymorphic(vcf)
    fam = f.trid.isin(poly).to_numpy()
    print("family: loci where a parent's allele differs in length", int(fam.sum()), flush=True)
    called, indel = hg002_indels(ref_dir, hg_dir, f)
    length = (f.end - f.start).to_numpy()
    strata = ((np.searchsorted(M06.LENGTHS, length, "left") - 1) * 1000
              + (np.searchsorted(M06.MOTIFS, f.motif_min.to_numpy(), "left") - 1) * 10
              + (np.searchsorted(M06.PURITY, f.share.to_numpy(), "left") - 1))
    _, sid = np.unique(strata, return_inverse=True)
    muts = f.mutations.to_numpy().astype(float)
    tails = pd.read_csv(TAILS, sep="\t", usecols=["trid", "fam"])
    tail_of = np.full(len(f), "elsewhere", dtype=object)
    for fam_name, label in (("Alu", "Alu tail"), ("L1", "L1 tail")):
        tail_of[f.trid.isin(set(tails.trid[tails.fam == fam_name])).to_numpy()] = label
    del tails
    out = {"biotypes": list(BIOTYPES), "parents": list(PARENTS), "loci": int(len(f)), "groups": {}}
    for g in GROUPS:
        gm = group == g
        d = {"loci": int(gm.sum())}
        for w in WHERE:
            m = gm & (where == w)
            d[w] = {"loci": int(m.sum()), "family_polymorphic": int(fam[m].sum()), "family_share": round(float(fam[m].mean()), 4) if m.any() else None,
                    "hg002_called": int(called[m].sum()), "hg002_indel": int((indel & called)[m].sum()),
                    "hg002_share": round(float(indel[m & called].mean()), 4) if (m & called).any() else None,
                    "new_tr_mutations": int(muts[m].sum()), "mean_length": round(float(length[m].mean()), 1) if m.any() else None}
        for label, x_name, y_name in (("sense vs antisense", "sense", "antisense"), ("sense vs between genes", "sense", "between genes"),
                                      ("antisense vs between genes", "antisense", "between genes")):
            m = gm & ((where == x_name) | (where == y_name))
            x = where[m] == x_name
            mc = m & called
            d[label] = {"family": mh_odds(fam[m], x, sid[m]),
                        "hg002": mh_odds(indel[mc], where[mc] == x_name, sid[mc]),
                        "new_tr_mutations": mh_rate(muts[m], x, sid[m])}
        # (the same comparison within the loci at Alu tails, at L1 tails and elsewhere: Alus and L1s in introns lean to one
        # orientation, so their tails could carry the difference between sense and antisense)
        d["sense vs antisense, by origin"] = {}
        for origin in ("Alu tail", "L1 tail", "elsewhere"):
            m = gm & ((where == "sense") | (where == "antisense")) & (tail_of == origin)
            mc = m & called
            d["sense vs antisense, by origin"][origin] = {
                "sense": int((m & (where == "sense")).sum()), "antisense": int((m & (where == "antisense")).sum()),
                "family": mh_odds(fam[m], where[m] == "sense", sid[m]), "hg002": mh_odds(indel[mc], where[mc] == "sense", sid[mc])}
        out["groups"][g] = d
        print(g, json.dumps({k: d[k] for k in ("sense", "antisense", "sense vs antisense")}), flush=True)
    os.makedirs(RESULTS, exist_ok=True)
    keep = group != ""
    pd.DataFrame({"trid": f.trid[keep], "group": group[keep], "rich_strand": np.where(rich_plus[keep], "+", "-"), "where": where[keep],
                  "family_polymorphic": fam[keep].astype(int), "hg002_called": called[keep].astype(int),
                  "hg002_indel": (indel & called)[keep].astype(int), "mutations": f.mutations[keep]}).to_csv(
        os.path.join(RESULTS, "loci.tsv.gz"), sep="\t", index=False)
    with open(os.path.join(RESULTS, "strands.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:6])
