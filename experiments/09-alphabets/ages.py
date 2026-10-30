"""Which alphabet does a retrotransposon's tail stay pure in, as the element ages?

Every autosomal Alu and L1 of UCSC's RepeatMasker annotation of hg38, and the runs of alphabets.py's census whose
origin is such a tail: the element's 3' end lies within 30 bases of where the run starts, on the run's A-rich strand
(the closest element, when several do). For each alphabet: how often an element's tail is a pure run in it, by the
element's divergence from its subfamily's consensus (RepeatMasker's milliDiv, a measure of age) and by subfamily; within
the 16 largest Alu subfamilies, by thirds of their divergence; and how long the runs are, and how much of them is their
commonest base, by age. The G/C alphabet has no tails and is left out.

    python3 -I ages.py RMSK           -> results/ages.json      (reads results/runs.tsv.gz)
"""

from __future__ import annotations

import gzip
import json
import os
import sys
from collections import Counter, defaultdict

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
AUTOSOMES = {f"chr{i}" for i in range(1, 23)}
NEAR = 30
ALPHA = ("purine/pyrimidine", "A/T", "A/C")
BINS = {"Alu": [0, 40, 60, 80, 100, 120, 140, 160, 180, 200, 10**6], "L1": [0, 50, 100, 150, 200, 250, 300, 10**6]}
LARGEST = 16  #: the largest Alu subfamilies, looked at within themselves
MIN_SUBFAMILY = 5000  #: subfamilies this large go into the correlation with age


def group(family: str, name: str) -> str:
    if family == "Alu":
        for p in ("AluY", "AluS", "AluJ"):
            if name.startswith(p):
                return p
        return "FLAM/FRAM" if name.startswith(("FLAM", "FRAM")) else "other"
    for p in ("L1HS", "L1PA", "L1PB", "L1MA", "L1MB", "L1MC", "L1MD", "L1ME"):
        if name.startswith(p):
            return p
    return "other L1M" if name.startswith("L1M") else "other"


def label(bins: list, i: int) -> str:
    a, b = bins[i] / 10, bins[i + 1] / 10
    return f"{a:g}-{b:g} %" if bins[i + 1] < 10**6 else f"{a:g} %+"


def read_elements(rmsk: str) -> dict:
    """(family, chrom) -> elements in start order: start, end, strand, name, milliDiv."""
    out = defaultdict(list)
    with gzip.open(rmsk, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[5] in AUTOSOMES and f[12] in ("Alu", "L1"):
                out[(f[12], f[5])].append((int(f[6]), int(f[7]), f[9], f[10], int(f[2])))
    for v in out.values():
        v.sort()
    return out


def match(elements: dict, runs: list[dict]) -> dict:
    """(family, chrom, element index) -> the tail runs: (alphabet, length, the commonest base's share)."""
    hits = defaultdict(list)
    for fam in ("Alu", "L1"):
        for c in AUTOSOMES:
            v = elements.get((fam, c), [])
            plus = [(e, i) for i, (s, e, st, _, _) in enumerate(v) if st == "+"]
            minus = [(s, i) for i, (s, e, st, _, _) in enumerate(v) if st == "-"]
            plus.sort()
            minus.sort()
            p_end, p_idx = np.array([x[0] for x in plus]), [x[1] for x in plus]
            m_start, m_idx = np.array([x[0] for x in minus]), [x[1] for x in minus]
            for r in runs:
                if r["chrom"] != c or r["origin"] != f"{fam} tail":
                    continue
                if r["a_strand"] == "+":
                    arr, idx, x = p_end, p_idx, r["start"]
                else:
                    arr, idx, x = m_start, m_idx, r["end"]
                j0, j1 = np.searchsorted(arr, x - NEAR, "left"), np.searchsorted(arr, x + NEAR, "right")
                best = min(((abs(int(arr[j]) - x), idx[j]) for j in range(j0, j1)), default=None)
                if best:
                    hits[(fam, c, best[1])].append((r["alphabet"], r["length"], r["major"]))
    return hits


def incidence(rows: list) -> dict:
    """rows: (alphabets with a tail run) per element."""
    n = len(rows)
    d = {"elements": n}
    for a in ALPHA:
        k = sum(a in x for x in rows)
        d[a] = {"tails": k, "per_1000": round(1000 * k / n, 2) if n else None}
    d["in several alphabets"] = sum(len(x) > 1 for x in rows)
    return d


def main(rmsk: str) -> None:
    runs = []
    with gzip.open(os.path.join(RESULTS, "runs.tsv.gz"), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["chrom"] in AUTOSOMES and r["alphabet"] in ALPHA and r["origin"] in ("Alu tail", "L1 tail"):
                r["start"], r["end"], r["length"], r["major"] = int(r["start"]), int(r["end"]), int(r["length"]), float(r["major"])
                runs.append(r)
    elements = read_elements(rmsk)
    hits = match(elements, runs)
    out = {"tail_runs_matched": sum(len(v) for v in hits.values()), "tail_runs": len(runs)}
    for fam in ("Alu", "L1"):
        table = []  # (group, name, milliDiv, alphabets, runs)
        for (f, c), v in elements.items():
            if f != fam:
                continue
            for i, (s, e, st, name, div) in enumerate(v):
                h = hits.get((fam, c, i), [])
                table.append((group(fam, name), name, div, {x[0] for x in h}, h))
        res: dict = {"elements": len(table)}
        bins = BINS[fam]
        bi = [next(k for k in range(len(bins) - 1) if bins[k] <= t[2] < bins[k + 1]) for t in table]
        res["by_divergence"] = {label(bins, k): incidence([t[3] for t, b in zip(table, bi) if b == k]) for k in range(len(bins) - 1)}
        res["by_group"] = {g: incidence([t[3] for t in table if t[0] == g]) for g in sorted({t[0] for t in table})}
        # the runs themselves, by age: length and the commonest base's share
        runs_by = {}
        for k in range(len(bins) - 1):
            for a in ALPHA:
                L = [x[1] for t, b in zip(table, bi) if b == k for x in t[4] if x[0] == a]
                M = [x[2] for t, b in zip(table, bi) if b == k for x in t[4] if x[0] == a]
                if L:
                    runs_by.setdefault(a, {})[label(bins, k)] = {"runs": len(L), "mean_length": round(float(np.mean(L)), 2),
                                                                  "median_length": float(np.median(L)), "major_share": round(float(np.mean(M)), 4)}
        res["runs_by_divergence"] = runs_by
        # subfamilies
        by_name = defaultdict(list)
        for t in table:
            by_name[t[1]].append(t)
        subs = []
        for name, v in by_name.items():
            if len(v) >= MIN_SUBFAMILY:
                d = incidence([t[3] for t in v])
                d.update({"subfamily": name, "median_milliDiv": float(np.median([t[2] for t in v]))})
                subs.append(d)
        subs.sort(key=lambda d: d["median_milliDiv"])
        res["subfamilies"] = subs
        if fam == "Alu":
            x = [d["median_milliDiv"] for d in subs]
            res["subfamily_correlation_with_age"] = {
                a: {"spearman_rho": round(float(r.statistic), 3), "p": float(f"{r.pvalue:.2g}"), "subfamilies": len(subs)}
                for a in ALPHA for r in [spearmanr(x, [d[a]["per_1000"] for d in subs])]}
            within = {}
            for name in sorted(by_name, key=lambda k: -len(by_name[k]))[:LARGEST]:
                v = by_name[name]
                q = np.quantile([t[2] for t in v], [1 / 3, 2 / 3])
                edges = [-1, q[0], q[1], 10**6]
                thirds = []
                for k in range(3):
                    sel = [t for t in v if edges[k] < t[2] <= edges[k + 1]]
                    d = incidence([t[3] for t in sel])
                    d["milliDiv"] = [int(min(t[2] for t in sel)), int(max(t[2] for t in sel))]
                    for a in ALPHA:
                        rr = [x for t in sel for x in t[4] if x[0] == a]
                        d[a]["mean_length"] = round(float(np.mean([x[1] for x in rr])), 2) if rr else None
                        d[a]["major_share"] = round(float(np.mean([x[2] for x in rr])), 4) if rr else None
                    thirds.append(d)
                within[name] = {"elements": len(v), "thirds": thirds}
            res["within_subfamilies"] = within
            res["within_subfamilies_summary"] = {
                a: {"rises_first_to_last_third": sum(w["thirds"][2][a]["per_1000"] > w["thirds"][0][a]["per_1000"] for w in within.values()),
                    "longer_in_last_third": sum((w["thirds"][2][a]["mean_length"] or 0) > (w["thirds"][0][a]["mean_length"] or 0) for w in within.values()),
                    "of": len(within)} for a in ALPHA}
        out[fam] = res
        print(fam, json.dumps(res["by_divergence"])[:1500], flush=True)
    print(json.dumps(out["Alu"].get("subfamily_correlation_with_age")), json.dumps(out["Alu"].get("within_subfamilies_summary")))
    with open(os.path.join(RESULTS, "ages.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
