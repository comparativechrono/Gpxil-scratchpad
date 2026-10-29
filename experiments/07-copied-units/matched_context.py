"""Does a repeat's own content change its new substitutions, at a fixed local context?

new_mutations.py sets each new substitution against the rate of its trinucleotide context, within its own locus. A bias
of wider contexts could still make a substitution that recreates a word of the locus look like a copy. Here the context
is held fixed: every new substitution's 6-base window (the base, two before it and three after) and its new base make a
stratum, with its locus's purity, length and unit size. In every catalogue locus, each place with that window is an
opportunity for the same substitution, of one of two sorts: the change would make a 6-base word that the locus already
holds elsewhere (a copy), or it would not. If only the window mattered, a stratum's mutations would fall among its
opportunities at random; conditioned on their number, the copies among them are binomial. The Mantel-Haenszel rate
ratio says how much more often the same substitution, in the same window, happens where it copies the locus.

Strata of loci, from coarse to fine: purity (bins at 0.6, 0.8, 0.9) and length (50, 100, 250); the same and the unit
size (1 to 6, 7 and more); experiment 06's finer bins of purity (9) and length (17), the unit size and the base-
composition class. Each is also run with the sites in low-complexity regions (Li 2014, LCR-hs38) left out, mutations
and opportunities alike, where short-read calls are least reliable and fewest de novo calls were made.

    python3 -I matched_context.py REF_DIR FEATURES LCR_BED
        (reads results/new_mutations.tsv.gz, and results/in_gnomad.tsv.gz when there is one)
        -> results/matched_context.json
"""

from __future__ import annotations

import gzip
import json
import math
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import new_mutations as nm  # noqa: E402

K, LEFT = nm.K, nm.LEFT
LENGTHS = (0, 50, 100, 250, 10 ** 9)  #: coarse length bins
FINE_LENGTHS = [0, 15, 20, 25, 30, 40, 50, 60, 75, 100, 125, 150, 200, 250, 350, 500, 1000, 10 ** 7]  #: experiment 06's
FINE_PURITY = [-1.0, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.98, 0.999, 1.0]  #: experiment 06's (right-closed bins)
LAYERS = {"purity, length": (0, 1), "purity, length, unit": (0, 1, 2), "fine purity, fine length, unit, class": (3, 4, 2, 5)}


def lbin(n: int) -> int:
    return next(i for i in range(len(LENGTHS) - 1) if LENGTHS[i] <= n < LENGTHS[i + 1])


def fine(edges: list, x: float) -> int:
    """i where edges[i] < x <= edges[i + 1] (experiment 06's bins)"""
    for i in range(len(edges) - 1):
        if x <= edges[i + 1]:
            return i
    return len(edges) - 2


def stratum(purity: float, length: int, unit: int, cls: int) -> tuple:
    """A locus's attributes: coarse purity bin, coarse length bin, unit size, fine purity bin, fine length bin, class."""
    return (int(np.digitize(purity, nm.PURITY[1:-1])), lbin(length), min(int(unit), 7), fine(FINE_PURITY, purity),
            fine(FINE_LENGTHS, length), int(cls))


def read_lcr(path: str) -> dict:
    out = defaultdict(list)
    with gzip.open(path, "rt") as fh:
        for line in fh:
            c, s, e = line.split()[:3]
            out[c].append((int(s), int(e)))
    return out


def read_mutations(path: str) -> list[dict]:
    out = []
    with gzip.open(path, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            if r["also_in_sasani"] == "True":
                continue
            out.append({"set": r["set"], "chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"],
                        "purity": float(r["purity"]), "length": int(r["length"]), "unit": int(r["unit"]), "class": r["class"],
                        "copies": r["copies"] == "True", "allele_balance": float(r["allele_balance"]) if r["allele_balance"] else None,
                        "cluster": r["cluster"] == "True", "transition": r["transition"] == "True", "kind": r["kind"]})
    return out


def main(ref_dir: str, features: str, lcr_bed: str) -> None:
    muts = read_mutations(os.path.join(nm.RESULTS, "new_mutations.tsv.gz"))
    gnomad = {}
    path = os.path.join(nm.RESULTS, "in_gnomad.tsv.gz")
    if os.path.exists(path):  # (in_gnomad.py's: whether each new allele is already in gnomAD)
        with gzip.open(path, "rt") as fh:
            fh.readline()
            for line in fh:
                c, p, a, st = line.rstrip("\n").split("\t")
                gnomad[(c, int(p), a)] = st
    for m in muts:
        m["gnomad"] = gnomad.get((m["chrom"], m["pos"], m["alt"]))
    loci = nm.Loci(features)
    lcr = read_lcr(lcr_bed)
    by_chrom = defaultdict(list)
    for m in muts:
        by_chrom[m["chrom"]].append(m)
        m["sid"] = stratum(m["purity"], m["length"], m["unit"], loci.classes.index(m["class"]))
    want: dict = defaultdict(set)
    opp: dict = defaultdict(lambda: [0, 0])  # (window, new base, locus attributes, in an LCR) -> [not a copy, a copy]
    for c in nm.ped.AUTOSOMES:
        with nm.ga.reference(ref_dir, c) as g:
            s = nm.LET[g[c].read("ref")].tobytes().decode()
        mask = np.zeros(len(s), dtype=bool)
        for a, b in lcr.get(c, []):
            mask[a:b] = True
        for m in by_chrom.get(c, []):
            m["window"] = s[m["pos"] - LEFT:m["pos"] - LEFT + K]
            assert m["window"][LEFT] == m["ref"], m
            m["lcr"] = bool(mask[m["pos"]])
            if "N" not in m["window"]:
                want[m["window"]].add(m["alt"])
        del s, mask
    words_wanted = set(want)
    sids: dict = {}
    print(f"{len(muts):,} mutations, {len(words_wanted):,} windows", flush=True)
    for c in nm.ped.AUTOSOMES:
        with nm.ga.reference(ref_dir, c) as g:
            s = nm.LET[g[c].read("ref")].tobytes().decode()
        mask = np.zeros(len(s), dtype=bool)
        for a, b in lcr.get(c, []):
            mask[a:b] = True
        a, b = loci.bounds[c]
        for i in range(a, b):
            s0, e0 = int(loci.start[i]), int(loci.end[i])
            seq = s[s0:e0]
            if len(seq) < K:
                continue
            words = {seq[x:x + K] for x in range(len(seq) - K + 1)}
            hit = words & words_wanted
            if not hit:
                continue
            sid = stratum(float(loci.purity[i]), e0 - s0, int(loci.unit[i]), int(loci.cls[i]))
            sid = sids.setdefault(sid, sid)
            for x in range(len(seq) - K + 1):
                w = seq[x:x + K]
                if w in hit:
                    inl = bool(mask[s0 + x + LEFT])
                    for alt in want[w]:
                        opp[(w, alt, sid, inl)][(w[:LEFT] + alt + w[LEFT + 1:]) in words] += 1
        del s, mask
        print(f"{c}: {len(opp):,} strata", flush=True)
    # the opportunities of each layer, with and without the LCR sites
    agg = {}
    for layer, cols in LAYERS.items():
        for outside in (False, True):
            d: dict = defaultdict(lambda: [0, 0])
            for (w, alt, sid, inl), v in opp.items():
                if outside and inl:
                    continue
                k = (w, alt) + tuple(sid[i] for i in cols)
                d[k][0] += v[0]
                d[k][1] += v[1]
            agg[(layer, outside)] = d
    pure = [m for m in muts if m["purity"] >= 0.6]
    groups = {"purity 0.6+": pure, "purity 0.6+, unit 1-6": [m for m in pure if m["unit"] <= 6],
              "purity 0.6+, unit 7+": [m for m in pure if m["unit"] > 6],
              "purity 0.9+": [m for m in pure if m["purity"] >= 0.9], "purity 0.6-0.9": [m for m in pure if m["purity"] < 0.9],
              "purity under 0.6": [m for m in muts if m["purity"] < 0.6],
              "purity 0.6+, outside clusters": [m for m in pure if not m["cluster"]],
              "purity 0.6+, allele balance 0.45+": [m for m in pure if m["allele_balance"] is not None and m["allele_balance"] >= 0.45],
              "purity 0.6+, transitions": [m for m in pure if m["transition"]],
              "purity 0.6+, transversions": [m for m in pure if not m["transition"]],
              "purity 0.6+, toward the consensus": [m for m in pure if m["kind"] == "toward"],
              "purity 0.6+, away from the consensus": [m for m in pure if m["kind"] == "away"],
              "purity 0.6+, neither": [m for m in pure if m["kind"] == "neither"],
              "purity 0.6+, under three units": [m for m in pure if not m["kind"]]}
    for s in sorted({m["set"] for m in muts}):
        groups[f"purity 0.6+, {s}"] = [m for m in pure if m["set"] == s]
    if gnomad:
        groups["purity 0.6+, allele absent from gnomAD"] = [m for m in pure if m["gnomad"] == "absent"]
        groups["purity 0.6+, allele in gnomAD"] = [m for m in pure if m["gnomad"] not in (None, "absent")]
        groups["purity 0.6+, allele in gnomAD, PASS"] = [m for m in pure if (m["gnomad"] or "").startswith("PASS")]
        groups["purity 0.6+, allele in gnomAD, not PASS"] = [m for m in pure if (m["gnomad"] or "").startswith("not PASS")]
        groups["purity under 0.6, allele absent from gnomAD"] = [m for m in muts if m["purity"] < 0.6 and m["gnomad"] == "absent"]
    res = {}
    for name, sel in groups.items():
        res[name] = {}
        for layer, cols in LAYERS.items():
            res[name][layer] = mantel_haenszel(sel, agg[(layer, False)], cols, False)
            res[name][layer + ", outside LCRs"] = mantel_haenszel(sel, agg[(layer, True)], cols, True)
    summary = {"mutations": len(muts), "windows": len(words_wanted), "strata": len(opp),
               "mutations_in_lcrs": sum(m["lcr"] for m in muts), "tests": res}
    with open(os.path.join(nm.RESULTS, "matched_context.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for name, t in res.items():
        print(name, json.dumps(t))


def mantel_haenszel(sel: list[dict], opp: dict, cols: tuple, outside: bool) -> dict:
    """Copies among a group's mutations against their strata's shares of copy opportunities (conditional binomial),
    and the Mantel-Haenszel rate ratio of copying to not copying (Greenland and Robins' variance)."""
    obs: dict = defaultdict(lambda: [0, 0])
    for m in sel:
        if "N" in m["window"] or (outside and m["lcr"]):
            continue
        obs[(m["window"], m["alt"]) + tuple(m["sid"][i] for i in cols)][int(m["copies"])] += 1
    o1 = e1 = var = r = s = gr = 0.0
    informative = 0
    for sk, (oo0, oo1) in obs.items():
        n0, n1 = opp[sk]
        n = n0 + n1
        assert n0 >= oo0 and n1 >= oo1 and n, (sk, n0, n1, oo0, oo1)
        o = oo0 + oo1
        p = n1 / n
        o1 += oo1
        e1 += o * p
        var += o * p * (1 - p)
        informative += 0 < p < 1
        r += oo1 * n0 / n
        s += oo0 * n1 / n
        gr += n1 * n0 * o / n ** 2
    out = {"mutations": sum(sum(v) for v in obs.values()), "strata": len(obs), "informative": informative,
           "copies": int(o1), "expected": round(e1, 1), "ratio": round(o1 / e1, 2) if e1 else None,
           "z": round((o1 - e1) / math.sqrt(var), 1) if var else None}
    if r and s:
        rr = r / s
        se = math.sqrt(gr / (r * s))
        out.update({"rate_ratio": round(rr, 2), "ci": [round(rr * math.exp(-1.96 * se), 2), round(rr * math.exp(1.96 * se), 2)]})
    return out


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
