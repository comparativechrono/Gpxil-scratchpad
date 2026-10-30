"""Pure runs in each of the genome's two-letter alphabets: how many, made of what, where from, how mutable.

A pure GA-sequence (experiment 04) is a run of more than 50 bases drawn from a two-letter alphabet, purines, on one
strand. The four bases make three such pairs of complementary alphabets, and each is a bit test on gpxil's 4-bit
reference codes, which gpxiltools' run finder takes as a set of letters:

* purine/pyrimidine: A/G, and C/T, its other strand (Albrecht-Buehler's GA- and TC-sequences);
* A/T and G/C, each its own other strand;
* A/C, and G/T, its other strand.

For every run of more than 50 bases on GRCh38 (all chromosomes; the statistics, the autosomes): its length, period and
purity (gpxiltools' period finder), the commonest base's share, its unit; its origin in UCSC's RepeatMasker annotation
of hg38, read on its A-rich strand as experiment 04 did (a transposable element's 3' end within 30 bases of the run's
start: a tail; else the class covering half of it or more); the share of its bases inside HG002's benchmarks (v4.2.1
and v5.0q); and the new tandem-repeat mutations of the Platinum Pedigree in catalogue loci at least half inside such
runs, against loci of the same length, motif size and purity (experiment 06's strata). The autosomes' dinucleotide
counts give each alphabet's expected number of runs under a first-order Markov model of the sequence.

    python3 -I alphabets.py REF_DIR RMSK FEATURES BED_V4.2.1 BED_V5.0Q      -> results/alphabets.json, results/runs.tsv.gz

REF_DIR: the Gpxil browser's GRCh38 gpxil files. RMSK: UCSC's rmsk.txt.gz for hg38. FEATURES: the locus measures
experiment 06 writes with --features. The BEDs: GIAB's HG002 benchmark regions on GRCh38.
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
sys.path.insert(0, os.path.join(EXP, "07-copied-units"))

import ga  # noqa: E402  (experiment 04)
import mutability as M06  # noqa: E402  (experiment 06: the strata)
import new_mutations as nm  # noqa: E402  (experiment 07: canonical motifs)
from gpxiltools import codes as GC  # noqa: E402
from gpxiltools.period import best_period  # noqa: E402
from gpxiltools.runs import iter_runs  # noqa: E402
from gpxiltools.vcf import load_bed  # noqa: E402

RESULTS = os.path.join(HERE, "results")
#: alphabet -> (the letters given to the run finder, both strands); a run on '-' is the complementary set's on '+'
ALPHABETS = {"purine/pyrimidine": ("AG", True), "A/T": ("AT", False), "G/C": ("GC", False), "A/C": ("AC", True)}
KINDS = {"purine/pyrimidine": ("A/G", "C/T"), "A/C": ("A/C", "G/T")}
MIN = ga.PURE  #: more than 50 bases
NEAR = 30  #: a TE's 3' end this close to a run's start makes the run its tail (experiment 04)
TE_CLASSES = {"SINE", "LINE", "LTR", "DNA", "Retroposon"}
IDX = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate((GC.BASE_A, GC.BASE_C, GC.BASE_G, GC.BASE_T)):
    IDX[_b] = IDX[_b | GC.SOFTMASK] = _i


def te_family(cls: str, fam: str) -> str:
    if fam in ("Alu", "L1"):
        return fam
    if fam == "SVA" or cls == "Retroposon":
        return "SVA"
    return {"SINE": "other SINE", "LINE": "other LINE", "LTR": "LTR"}.get(cls, "DNA")


def read_rmsk(path: str) -> dict:
    """chrom -> arrays of RepeatMasker's annotation, by start: start, end, strand, name, class, family, milliDiv."""
    out = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[5]].append((int(f[6]), int(f[7]), f[9], f[10], f[11], f[12], int(f[2])))
    res = {}
    for c, v in out.items():
        v.sort()
        res[c] = {"start": np.array([x[0] for x in v]), "end": np.array([x[1] for x in v]), "strand": [x[2] for x in v],
                  "name": [x[3] for x in v], "cls": [x[4] for x in v], "fam": [x[5] for x in v], "div": [x[6] for x in v]}
    return res


def origin(rm: dict, c: str, s: int, e: int, a_plus: bool) -> tuple[str, str | None]:
    """A TE whose 3' end lies within NEAR of the run's start, on the run's A-rich strand: '<family> tail' and its
    'name:milliDiv'; else the class of the annotation covering half of the run or more, or 'none'."""
    r = rm.get(c)
    if r is None:
        return "none", None
    st, en = r["start"], r["end"]
    best = None
    i0 = np.searchsorted(st, s - 20000, "left")
    i1 = np.searchsorted(st, e + NEAR, "right")
    for i in range(i0, i1):
        if r["cls"][i] not in TE_CLASSES:
            continue
        if a_plus and r["strand"][i] == "+" and abs(int(en[i]) - s) <= NEAR:
            d = abs(int(en[i]) - s)
        elif not a_plus and r["strand"][i] == "-" and abs(int(st[i]) - e) <= NEAR:
            d = abs(int(st[i]) - e)
        else:
            continue
        if best is None or d < best[0]:
            best = (d, te_family(r["cls"][i], r["fam"][i]), f'{r["name"][i]}:{r["div"][i]}')
    if best:
        return best[1] + " tail", best[2]
    cover = Counter()
    for i in range(i0, i1):
        a, b = max(s, int(st[i])), min(e, int(en[i]))
        if b > a:
            k = r["cls"][i]
            cover[te_family(k, r["fam"][i]) + " body" if k in TE_CLASSES else k] += b - a
    if cover:
        k, v = cover.most_common(1)[0]
        if v >= (e - s) / 2:
            return k, None
    return "none", None


def unit_of(seq: str, p: int) -> str:
    """The run's commonest p-base word, in its canonical form (the least rotation of it or its reverse complement)."""
    words = Counter(seq[i:i + p] for i in range(0, len(seq) - p + 1))
    w = words.most_common(1)[0][0] if words else seq[:p]
    return nm.canon(w)[0] if set(w) <= set("ACGT") else w


def markov_expectation(di: np.ndarray, length: int = MIN) -> dict:
    """Expected maximal runs of `length` bases or more in each alphabet, in sequence of the autosomes' first-order
    Markov chain (from their dinucleotide counts), and with independent bases of their composition."""
    D = di[:4, :4].astype(float)
    pi = D.sum(1) / D.sum()
    P = D / D.sum(1, keepdims=True)
    n = D.sum()
    out = {}
    for name, letters in (("A/G", "AG"), ("C/T", "CT"), ("A/T", "AT"), ("G/C", "GC"), ("A/C", "AC"), ("G/T", "GT")):
        X = ["ACGT".index(b) for b in letters]
        notx = [i for i in range(4) if i not in X]
        start = np.array([sum(pi[j] * P[j, i] for j in notx) for i in X])
        v = start.copy()
        for _ in range(length - 1):
            v = v @ P[np.ix_(X, X)]
        px = pi[X].sum()
        out[name] = {"markov": float(f"{n * v.sum():.3g}"), "independent": float(f"{n * (1 - px) * px ** length:.3g}")}
    return out


def census(ref_dir: str, rm: dict) -> tuple[list[dict], dict, np.ndarray]:
    runs, masks = [], {}
    di = np.zeros(25, dtype=np.int64)
    t0 = time.time()
    for c in ga.CHROMS:
        with ga.reference(ref_dir, c) as g:
            codes = g[c].read("ref")
        if c in ga.AUTOSOMES:
            for a in range(0, len(codes) - 1, 1 << 24):
                x = IDX[codes[a:a + (1 << 24) + 1]].astype(np.int16)
                di += np.bincount(x[:-1] * 5 + x[1:], minlength=25)
        letters = nm.LET[codes].tobytes().decode()
        del codes
        with ga.reference(ref_dir, c) as g:
            for name, (of, both) in ALPHABETS.items():
                iv = []
                for r in iter_runs(g, of=of, min_length=MIN, both_strands=both):
                    s, e = r.start, r.end
                    seq = letters[s:e]
                    cnt = Counter(seq)
                    if name in KINDS:
                        kind = KINDS[name][0] if r.strand == "+" else KINDS[name][1]
                        a_plus = r.strand == "+"
                    else:
                        kind = name
                        a_plus = cnt["A"] >= cnt["T"] if of == "AT" else cnt["G"] >= cnt["C"]
                    per, share, _ = best_period(np.array([nm.CODE.get(x, 15) for x in seq], dtype=np.uint8), 12)
                    o, te = origin(rm, c, s, e, a_plus)
                    runs.append({"alphabet": name, "kind": kind, "chrom": c, "start": s, "end": e, "length": e - s,
                                 "period": per, "purity": round(share, 3), "major": round(max(cnt.values()) / (e - s), 3),
                                 "unit": unit_of(seq, per), "a_strand": "+" if a_plus else "-", "origin": o, "te": te or ""})
                    iv.append((s, e))
                masks[(name, c)] = sorted(iv)
        del letters
        print(c, len(runs), f"[{time.time() - t0:.0f}s]", flush=True)
    return runs, masks, di.reshape(5, 5)


def coverage(ref_dir: str, beds: dict, masks: dict) -> dict:
    """The share of each alphabet's run bases inside each benchmark's regions; of the genome's, over its bases that
    are not N."""
    cov = {}
    for b, bed in beds.items():
        tot = Counter()
        for c in ga.AUTOSOMES:
            with ga.reference(ref_dir, c) as g:
                known = IDX[g[c].read("ref")] < 4
            n = len(known)
            m = load_bed(bed, c, n)
            for name in ALPHABETS:
                for s, e in masks[(name, c)]:
                    tot[(name, "bases")] += e - s
                    tot[(name, "called")] += int(m[s:e].sum()) if m is not None else 0
            tot[("genome", "bases")] += int(known.sum())
            tot[("genome", "called")] += int((m & known).sum()) if m is not None else 0
        cov[b] = {k: round(tot[(k, "called")] / max(tot[(k, "bases")], 1), 4) for k in list(ALPHABETS) + ["genome"]}
        print("coverage", b, cov[b], flush=True)
    return cov


def tr_mutations(features: str, masks: dict) -> dict:
    """New TR mutations in loci at least half inside an alphabet's runs, against the other loci of their strata."""
    import pandas as pd

    f = pd.read_csv(features, sep="\t", usecols=["chrom", "start", "end", "motif_min", "purity", "mutations"])
    f = f[f.chrom.isin(ga.AUTOSOMES)].sort_values(["chrom", "start"])
    length = (f.end - f.start).to_numpy()
    strata = ((np.searchsorted(M06.LENGTHS, length, "left") - 1) * 1000 + (np.searchsorted(M06.MOTIFS, f.motif_min.to_numpy(), "left") - 1) * 10
              + (np.searchsorted(M06.PURITY, f.purity.to_numpy(), "left") - 1))
    _, sid = np.unique(strata, return_inverse=True)
    muts = f.mutations.to_numpy()
    chrom = f.chrom.to_numpy()
    starts, ends = f.start.to_numpy(), f.end.to_numpy()
    out = {"loci": int(len(f)), "mutations": int(muts.sum())}
    for name in ALPHABETS:
        inside = np.zeros(len(f), dtype=bool)
        for c in ga.AUTOSOMES:
            sel = np.flatnonzero(chrom == c)
            iv = masks[(name, c)]
            if not iv or not len(sel):
                continue
            rs, re_ = np.array([x[0] for x in iv]), np.array([x[1] for x in iv])  # (an alphabet's runs don't overlap)
            st, en = starts[sel], ends[sel]
            cum = np.concatenate([[0], np.cumsum(re_ - rs)])
            i0 = np.searchsorted(re_, st, "right")  # the first run ending after the locus starts
            i1 = np.searchsorted(rs, en, "left")  # one past the last run starting before it ends
            some = i1 > i0
            a0, a1 = np.minimum(i0, len(rs) - 1), np.maximum(i1 - 1, 0)
            ov = np.where(some, cum[i1] - cum[i0] - np.clip(st - rs[a0], 0, None) - np.clip(re_[a1] - en, 0, None), 0)
            inside[sel] = ov >= (en - st) / 2
        k = sid.max() + 1
        n_i = np.bincount(sid[inside], minlength=k)
        n_o = np.bincount(sid[~inside], minlength=k)
        m_o = np.bincount(sid[~inside], weights=muts[~inside], minlength=k)
        ok = n_o > 0
        exp = float((n_i[ok] * m_o[ok] / n_o[ok]).sum())
        m_in = int(muts[inside].sum())
        lo, hi = M06.poisson_ci(m_in)
        out[name] = {"loci": int(inside.sum()), "mutations": m_in, "expected_from_strata": round(exp, 1),
                     "ratio": round(m_in / exp, 2) if exp else None, "ci": [round(lo / exp, 2), round(hi / exp, 2)] if exp else None,
                     "share_of_all_mutations": round(m_in / int(muts.sum()), 4)}
        print("mutations", name, out[name], flush=True)
    return out


def summary(runs: list[dict]) -> dict:
    res = {}
    for name in ALPHABETS:
        for where, chroms in (("all chromosomes", set(ga.CHROMS)), ("autosomes", set(ga.AUTOSOMES))):
            sel = [r for r in runs if r["alphabet"] == name and r["chrom"] in chroms]
            L = np.array([r["length"] for r in sel])
            d = {"runs": len(sel), "bases": int(L.sum()) if len(L) else 0, "by_kind": dict(Counter(r["kind"] for r in sel))}
            if where == "autosomes":
                d.update({"length_quantiles_50_90_99": [int(np.percentile(L, q)) for q in (50, 90, 99)] if len(L) else None,
                          "longest": int(L.max()) if len(L) else 0, "one_base_98": sum(r["major"] >= 0.98 for r in sel),
                          "top_units": Counter(r["unit"] for r in sel).most_common(12),
                          "periods": dict(sorted(Counter(min(r["period"], 13) for r in sel).items())),
                          "origin": dict(Counter(r["origin"] for r in sel).most_common())})
            res.setdefault(name, {})[where] = d
    return res


def main(ref_dir: str, rmsk: str, features: str, bed4: str, bed5: str) -> None:
    rm = read_rmsk(rmsk)
    runs, masks, di = census(ref_dir, rm)
    os.makedirs(RESULTS, exist_ok=True)
    cols = ["alphabet", "kind", "chrom", "start", "end", "length", "period", "purity", "major", "unit", "a_strand", "origin", "te"]
    with gzip.open(os.path.join(RESULTS, "runs.tsv.gz"), "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in runs:
            fh.write("\t".join(str(r[k]) for k in cols) + "\n")
    out = {"min_length": MIN, "census": summary(runs),
           "composition_autosomes": dict(zip("ACGT", (di[:4, :4].sum(1) / di[:4, :4].sum()).round(4).tolist())),
           "dinucleotides_autosomes": di.tolist(), "expected_runs_autosomes": markov_expectation(di)}
    print(json.dumps(out["census"])[:3000], flush=True)
    out["coverage"] = coverage(ref_dir, {"v4.2.1": bed4, "v5.0q": bed5}, masks)
    out["tr_mutations"] = tr_mutations(features, masks)
    with open(os.path.join(RESULTS, "alphabets.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:6])
