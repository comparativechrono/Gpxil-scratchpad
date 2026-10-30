"""Repeat arrays that begin at a transposable element's 3' end, read outward from the element: how pure each part of
them is in the reference, by the element's age.

Arrays: the catalogue's autosomal loci (Platinum TRs v1.0; experiment 06's measures of them) that begin within 30 bases
of an Alu's or an L1's 3' end, on its strand, the rule by which experiments 04 and 09 give a run its origin. Each is
read on the element's strand, from the end next to the element (proximal) to the far end (distal). Its age is the
element's divergence from its consensus in RepeatMasker's annotation, as in experiment 09.

For each array, on GRCh38's gpxil files: its length, the share of A in each tenth of it, and its longest stretch in
frame (each base the same as the base one shortest motif on; for a poly(A), its longest run of A): where it lies and
how long it is. For the pure two-letter runs of experiment 09 that begin at an Alu's 3' end, the share of the run's
second letter (G in a purine run, T in an A/T run, C in an A/C run) in each tenth. And the new TR mutations of the
Platinum Pedigree at TE tails, against the catalogue's other loci of their strata (experiment 06), as in experiment 09.

    python3 -I tails.py REF_DIR RMSK FEATURES      -> results/tail_loci.tsv.gz, results/tails.json

REF_DIR: the Gpxil browser's GRCh38 gpxil files. RMSK: UCSC's rmsk.txt.gz for hg38. FEATURES: the locus measures
experiment 06 writes with --features.
"""

from __future__ import annotations

import gzip
import json
import os
import sys
import time
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(EXP, "04-ga-sequences"))
sys.path.insert(0, os.path.join(EXP, "06-repeat-mutability"))
sys.path.insert(0, os.path.join(EXP, "07-copied-units"))

import ga  # noqa: E402  (experiment 04)
import mutability as M06  # noqa: E402  (experiment 06: the strata)
import new_mutations as nm  # noqa: E402  (experiment 07: base letters)

RESULTS = os.path.join(HERE, "results")
RUNS09 = os.path.join(EXP, "09-alphabets", "results", "runs.tsv.gz")
NEAR = 30  #: a TE's 3' end this close to an array's start makes the array its tail (experiments 04 and 09)
FAMILIES = ("Alu", "L1")
AGES = [0, 4, 8, 12, 16, 20, 100]  #: divergence bins (%), as in experiment 09
DECILES = 10
COMP = str.maketrans("ACGTN", "TGCAN")
SECOND = {"purine/pyrimidine": "G", "A/T": "T", "A/C": "C"}  #: a two-letter tail run's letter other than A


def age_label(d: float) -> str:
    i = int(np.searchsorted(AGES, d, "right") - 1)
    i = min(max(i, 0), len(AGES) - 2)
    return f"{AGES[i]}-{AGES[i + 1]} %" if AGES[i + 1] < 100 else f"{AGES[i]} %+"


def te_ends(rmsk: str) -> dict:
    """chrom -> strand -> (3' ends sorted, family, name, divergence %, consensus bases missing at the 3' end)."""
    acc = defaultdict(lambda: defaultdict(list))
    with gzip.open(rmsk, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if f[12] not in FAMILIES or f[5] not in ga.AUTOSOMES:
                continue
            plus = f[9] == "+"
            end3 = int(f[7]) if plus else int(f[6])
            missing = -int(f[15]) if plus else -int(f[13])
            acc[f[5]][f[9]].append((end3, f[12], f[10], int(f[2]) / 10, missing))
    out = {}
    for c, d in acc.items():
        out[c] = {}
        for s, v in d.items():
            v.sort()
            out[c][s] = (np.array([x[0] for x in v]), np.array([x[1] for x in v]), np.array([x[2] for x in v]),
                         np.array([x[3] for x in v]), np.array([x[4] for x in v]))
    return out


def match_tails(loci, ends: dict):
    """For each locus, the TE whose 3' end lies within NEAR of the locus's start on the TE's strand (for a TE on '-',
    of the locus's end): the nearest, else none. Returns per-locus arrays (index into loci order)."""
    n = len(loci)
    fam = np.full(n, "", dtype=object)
    strand = np.full(n, "", dtype=object)
    name = np.full(n, "", dtype=object)
    div = np.full(n, np.nan)
    missing = np.full(n, -1)
    dist = np.full(n, NEAR + 1)
    chrom = loci.chrom.to_numpy()
    for c in ga.AUTOSOMES:
        sel = np.flatnonzero(chrom == c)
        if not len(sel) or c not in ends:
            continue
        for s, key in (("+", "start"), ("-", "end")):
            if s not in ends[c]:
                continue
            pos, f_, nm_, dv, ms = ends[c][s]
            q = loci[key].to_numpy()[sel]
            i = np.searchsorted(pos, q)
            for j in (i - 1, i):
                ok = (j >= 0) & (j < len(pos))
                jj = np.clip(j, 0, len(pos) - 1)
                d = np.where(ok, np.abs(q - pos[jj]), NEAR + 1)
                better = d < dist[sel]
                k = sel[better]
                dist[k] = d[better]
                fam[k], strand[k], name[k], div[k], missing[k] = f_[jj[better]], s, nm_[jj[better]], dv[jj[better]], ms[jj[better]]
    tail = dist <= NEAR
    return tail, fam, strand, name, div, missing, dist


def in_frame_longest(seq: str, p: int) -> tuple[int, int]:
    """The longest stretch [a, b) in which every base equals the base p on (the stretch's last p bases end it): for a
    poly(A), its longest run of one base. (0, 0) for an empty sequence."""
    n = len(seq)
    if n == 0:
        return 0, 0
    p = max(1, min(p, n))
    best = (0, min(p, n))
    a = 0
    for i in range(n - p):
        if seq[i] != seq[i + p]:
            if i + p - a > best[1] - best[0]:
                best = (a, i + p)
            a = i + 1
    if n - a > best[1] - best[0]:
        best = (a, n)
    return best


def measure(ref_dir: str, loci) -> dict:
    """Per array (on the TE's strand, proximal first): A share by tenth, the longest stretch in frame."""
    out = {k: np.zeros((len(loci), DECILES)) if k == "a" else np.zeros(len(loci), dtype=np.int64)
           for k in ("a", "best_lo", "best_hi")}
    t0 = time.time()
    for c in ga.AUTOSOMES:
        idx = np.flatnonzero(loci.chrom.to_numpy() == c)
        if not len(idx):
            continue
        with ga.reference(ref_dir, c) as g:
            codes = g[c].read("ref")
        let = nm.LET[codes].tobytes().decode().upper()
        del codes
        st, en = loci.start.to_numpy()[idx], loci.end.to_numpy()[idx]
        sd, mm = loci.te_strand.to_numpy()[idx], loci.motif_min.to_numpy()[idx]
        for k, s, e, d, m in zip(idx, st, en, sd, mm):
            seq = let[s:e]
            if d == "-":
                seq = seq.translate(COMP)[::-1]
            n = len(seq)
            a = np.frombuffer(seq.encode(), dtype=np.uint8) == ord("A")
            dec = np.minimum((np.arange(n) * DECILES) // n, DECILES - 1)
            out["a"][k] = np.bincount(dec, weights=a, minlength=DECILES) / np.bincount(dec, minlength=DECILES)
            out["best_lo"][k], out["best_hi"][k] = in_frame_longest(seq, int(m))
        del let
        print(c, len(idx), f"[{time.time() - t0:.0f}s]", flush=True)
    return out


def second_letters(ref_dir: str) -> dict:
    """Experiment 09's pure runs that begin at an Alu's 3' end: the share of the second letter by tenth, by age."""
    import pandas as pd
    runs = pd.read_csv(RUNS09, sep="\t")
    runs = runs[(runs.origin == "Alu tail") & runs.chrom.isin(ga.AUTOSOMES) & runs.alphabet.isin(list(SECOND))].copy()
    runs["div"] = runs.te.str.split(":").str[1].astype(float) / 10
    acc = defaultdict(lambda: {"runs": 0, "bases": 0, "second": 0, "by_tenth": np.zeros(DECILES), "tenth_bases": np.zeros(DECILES)})
    for c, d in runs.groupby("chrom"):
        with ga.reference(ref_dir, c) as g:
            codes = g[c].read("ref")
        let = nm.LET[codes].tobytes().decode().upper()
        del codes
        for s, e, st, al, dv in zip(d.start, d.end, d.a_strand, d.alphabet, d["div"]):
            seq = let[s:e]
            if st == "-":
                seq = seq.translate(COMP)[::-1]
            n = len(seq)
            m = np.frombuffer(seq.encode(), dtype=np.uint8) == ord(SECOND[al])
            dec = np.minimum((np.arange(n) * DECILES) // n, DECILES - 1)
            for key in (al, f"{al}|{age_label(dv)}"):
                x = acc[key]
                x["runs"] += 1
                x["bases"] += n
                x["second"] += int(m.sum())
                x["by_tenth"] += np.bincount(dec, weights=m, minlength=DECILES)
                x["tenth_bases"] += np.bincount(dec, minlength=DECILES)
    out = {}
    for key, x in sorted(acc.items()):
        al, _, age = key.partition("|")
        out.setdefault(al, {"second_letter": SECOND[al], "by_age": {}})
        d = {"runs": x["runs"], "mean_length": round(x["bases"] / x["runs"], 1), "second_share": round(x["second"] / x["bases"], 6),
             "by_tenth": [round(v, 6) for v in (x["by_tenth"] / x["tenth_bases"]).tolist()]}
        if age:
            out[al]["by_age"][age] = d
        else:
            out[al]["all"] = d
    return out


def mutations_vs_strata(loci_all, tail: np.ndarray, fam: np.ndarray) -> dict:
    """New TR mutations at TE-tail loci against the catalogue's other loci of the same length, motif size and purity
    (experiment 06's strata; indirect standardisation), as experiment 09 does for runs."""
    length = (loci_all.end - loci_all.start).to_numpy()
    strata = ((np.searchsorted(M06.LENGTHS, length, "left") - 1) * 1000
              + (np.searchsorted(M06.MOTIFS, loci_all.motif_min.to_numpy(), "left") - 1) * 10
              + (np.searchsorted(M06.PURITY, loci_all.purity.to_numpy(), "left") - 1))
    _, sid = np.unique(strata, return_inverse=True)
    muts = loci_all.mutations.to_numpy()
    k = sid.max() + 1
    out = {"loci": int(len(muts)), "mutations": int(muts.sum())}
    for name in FAMILIES + ("either",):
        inside = tail & ((fam == name) if name != "either" else True)
        n_i = np.bincount(sid[inside], minlength=k)
        n_o = np.bincount(sid[~inside], minlength=k)
        m_o = np.bincount(sid[~inside], weights=muts[~inside], minlength=k)
        ok = n_o > 0
        exp = float((n_i[ok] * m_o[ok] / n_o[ok]).sum())
        m_in = int(muts[inside].sum())
        lo, hi = M06.poisson_ci(m_in)
        out[f"{name} tails"] = {"loci": int(inside.sum()), "share_of_loci": round(float(inside.mean()), 4), "mutations": m_in,
                                "share_of_mutations": round(m_in / max(int(muts.sum()), 1), 4),
                                "expected_from_strata": round(exp, 1), "ratio": round(m_in / exp, 2) if exp else None,
                                "ci": [round(lo / exp, 2), round(hi / exp, 2)] if exp else None}
        print("mutations", name, out[f"{name} tails"], flush=True)
    return out


def summarise(t, meas: dict) -> dict:
    n = (t.end - t.start).to_numpy()
    blen = meas["best_hi"] - meas["best_lo"]
    centre = (meas["best_lo"] + meas["best_hi"]) / 2
    out = {}
    for f in FAMILIES:
        for age in ["all"] + [age_label(a + 0.01) for a in AGES[:-1]]:
            m = (t.fam.to_numpy() == f) & ((t.age.to_numpy() == age) if age != "all" else True)
            if not m.any():
                continue
            a = meas["a"][m]
            d = {"arrays": int(m.sum()), "mean_length": round(float(n[m].mean()), 2),
                 "a_share": round(float((a.mean(1)).mean()), 6),
                 "a_by_tenth": [round(float(v), 6) for v in a.mean(0)],
                 "longest_in_frame_share_of_array": round(float((blen[m] / n[m]).mean()), 6),
                 "longest_in_frame_mean_bases": round(float(blen[m].mean()), 2),
                 "longest_in_frame_centre_relative": round(float((centre[m] / n[m]).mean()), 6),
                 "longest_in_frame_centre_bases_from_element": round(float(centre[m].mean()), 2),
                 "longest_in_frame_in_proximal_half": round(float((centre[m] < n[m] / 2).mean()), 6),
                 "poly_a_loci": int((t.motif_min.to_numpy()[m] == 1).sum())}
            out.setdefault(f, {})[age] = d
    return out


def main(ref_dir: str, rmsk: str, features: str) -> None:
    import pandas as pd
    cols = ["chrom", "start", "end", "trid", "motif_min", "purity", "class", "mutations"]
    loci = pd.read_csv(features, sep="\t", usecols=cols, dtype={"chrom": "category", "start": "int64", "end": "int64",
                       "motif_min": "int32", "purity": "float32", "class": "category", "mutations": "int16"})
    loci = loci[loci.chrom.isin(ga.AUTOSOMES)].reset_index(drop=True)
    print("loci", len(loci), flush=True)
    tail, fam, strand, name, div, missing, dist = match_tails(loci, te_ends(rmsk))
    out = {"near": NEAR, "loci": int(len(loci)), "tails": {f: int((tail & (fam == f)).sum()) for f in FAMILIES}}
    out["new_tr_mutations"] = mutations_vs_strata(loci, tail, fam)
    t = loci[tail].copy()
    t["fam"], t["te_strand"], t["te_name"], t["div"] = fam[tail], strand[tail], name[tail], div[tail]
    t["te_end_missing"], t["te_distance"] = missing[tail], dist[tail]
    t["age"] = [age_label(x) for x in t["div"]]
    t = t.reset_index(drop=True)
    meas = measure(ref_dir, t)
    for i in range(DECILES):
        t[f"a{i}"] = meas["a"][:, i].round(3)
    t["best_lo"], t["best_hi"] = meas["best_lo"], meas["best_hi"]
    os.makedirs(RESULTS, exist_ok=True)
    t.drop(columns=["class"]).to_csv(os.path.join(RESULTS, "tail_loci.tsv.gz"), sep="\t", index=False)
    out["arrays"] = summarise(t, meas)
    out["two_letter_tail_runs"] = second_letters(ref_dir)
    with open(os.path.join(RESULTS, "tails.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out["arrays"])[:2500])


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
