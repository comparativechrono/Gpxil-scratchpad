"""The copying tests on the variants of 76,215 people: gnomAD v4.1 genomes' singletons in repeats.

A singleton (allele count 1) is the nearest a population sample comes to a new mutation: its allele is almost surely
the derived one, and too recent for selection or gene conversion to have shaped it. Two sets of loci, each read on the
+ strand:

* the pure GA/TC-sequences of the autosomes, whose gnomAD PASS records experiment 04 read (its local copy,
  GNOMAD_RECORDS, also holds 10,000 random 100-base windows elsewhere: the background); period and purity from
  gpxiltools' period finder;
* 20,000 of the catalogue's other tandem repeats (purity 0.6 or more, 20 bases or more), drawn at random by
  gnomad_tr.py, which read their records the same way (TR_RECORDS); purity and unit from experiment 06's measures.

The SNVs kept are biallelic (no other record at the position) and their 6-base window lies in the locus; the
well-called ones were called in 90 % of the genomes or more.

* **Within the locus**, as new_mutations.py: each singleton's chance of copying (its new 6-base word is elsewhere in the
  locus) under the rates of the trinucleotide contexts, here the background's singletons'; and of moving toward the
  consensus unit, away from it, or neither.
* **At a fixed context**, as matched_context.py: the same 6-base window and new base, in loci of the same purity,
  length and unit, where the change would copy and where it would not; Mantel-Haenszel rate ratio.

Allele-count classes beyond singletons are given too, for the minor alleles only (AF < 0.5).

    python3 -I population.py REF_DIR GNOMAD_RECORDS [TR_RECORDS FEATURES]       -> results/population.json
"""

from __future__ import annotations

import gzip
import json
import math
import os
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "04-ga-sequences"))

import new_mutations as nm  # noqa: E402
import gnomad as G  # noqa: E402  (experiment 04)
from gpxiltools.period import best_period  # noqa: E402

K, LEFT = nm.K, nm.LEFT
GA_LENGTHS = (51, 76, 101, 151, 251, 10 ** 9)
TR_LENGTHS = (0, 50, 100, 250, 10 ** 9)
PURITY = nm.PURITY
CLASSES = ("singleton", "rare", "uncommon", "common")


def letters(codes: np.ndarray) -> str:
    return nm.LET[codes].tobytes().decode()


def bin_of(edges: tuple, x: float) -> int:
    return next(i for i in range(len(edges) - 1) if edges[i] <= x < edges[i + 1])


def biallelic_snvs(recs: list) -> list:
    """The SNV records at positions that hold no other record."""
    n = Counter(r[0] for r in recs)
    return [r for r in recs if n[r[0]] == 1 and len(r[1]) == 1 and len(r[2]) == 1 and r[1] in "ACGT" and r[2] in "ACGT"]


def background_rates(records: str):
    """The background windows' singleton rate of each (pyrimidine trinucleotide, new base)."""
    tri_n, tri_mut = Counter(), Counter()
    for (ch, s, e, kind), codes, a, recs in G.cached(records):
        if kind != "background":
            continue
        seq = letters(codes)
        for x in range(s + 1, e - 1):
            tri = seq[x - a - 1:x - a + 2]
            if "N" not in tri:
                tri_n[nm.pyr(tri, "A")[0]] += 1
        for pos, ref, alt, ac, af, an in biallelic_snvs(recs):
            if ac == 1 and s + 1 <= pos < e - 1:
                tri = seq[pos - a - 1:pos - a + 2]
                if "N" not in tri and tri[1] == ref:
                    tri_mut[nm.pyr(tri, alt)] += 1

    def rate(seq_: str, x: int, alt: str) -> float:
        tri = seq_[x - 1:x + 2]
        if "N" in tri:
            return 0.0
        t, b = nm.pyr(tri, alt)
        return tri_mut[(t, b)] / tri_n[t] if tri_n[t] else 0.0

    return rate, sum(tri_mut.values())


def analyse(loci, rate) -> dict:
    """loci: (chrom, start, end, sequence, SNV records, purity, unit, length bins, maximum period)."""
    rows = []
    opp: dict = defaultdict(lambda: [0, 0])
    for ch, s, e, seq, snvs, purity, unit, lengths, maxp in loci:
        codes = np.array([nm.CODE.get(c, 15) for c in seq], dtype=np.uint8)
        period, share, _ = best_period(codes, maxp)
        if purity is None:
            purity, unit = share, period
        cons = nm.consensus(seq, period)
        own = {seq[x:x + K] for x in range(len(seq) - K + 1)}
        sid = (bin_of(PURITY, purity), bin_of(lengths, len(seq)), min(unit, 7))
        for x in range(LEFT, len(seq) - (K - 1 - LEFT)):
            w = seq[x - LEFT:x - LEFT + K]
            if "N" in w:
                continue
            for alt in "ACGT":
                if alt != w[LEFT]:
                    opp[(w, alt) + sid][(w[:LEFT] + alt + w[LEFT + 1:]) in own] += 1
        if not snvs:
            continue
        ex = nm.expectations(seq, period, cons, {"tri": rate}, own, {})
        if not ex:
            continue
        for pos, ref, alt, ac, af, an in snvs:
            j = pos - s
            if seq[j] != ref:
                continue
            w = seq[j - LEFT:j - LEFT + K]
            if "N" in w:
                continue
            new = w[:LEFT] + alt + w[LEFT + 1:]
            c = cons[j % period] if cons else None
            kd = None if c is None else "away" if ref == c else "toward" if alt == c else "neither"
            ik = ex.get("in_kind", {}).get(kd, {}) if kd else {}
            rows.append({"class": G.af_class(ac, af), "af": af, "well_called": an >= G.WELL_CALLED, "purity": purity, "sid": sid,
                         "window": w, "alt": alt, "transition": {ref, alt} in ({"A", "G"}, {"C", "T"}), "kind": kd,
                         "copies": new in own, "expected_copies": ex["expected_copies"],
                         "near": alt in {seq[y] for y in (j - period, j + period) if 0 <= y < len(seq)},
                         "expected_near": ex["expected_near"], "expected_toward": ex.get("expected_toward"),
                         "expected_copies_in_kind": ik.get("copies")})
    single = [r for r in rows if r["class"] == "singleton"]
    groups = {"singletons": single, "singletons, well called": [r for r in single if r["well_called"]],
              "singletons, purity under 0.6": [r for r in single if r["purity"] < 0.6],
              "singletons, purity 0.6-0.9": [r for r in single if 0.6 <= r["purity"] < 0.9],
              "singletons, purity 0.9+": [r for r in single if r["purity"] >= 0.9],
              "singletons, unit 1-6": [r for r in single if r["sid"][2] <= 6],
              "singletons, unit 7+": [r for r in single if r["sid"][2] > 6],
              "singletons, transitions": [r for r in single if r["transition"]],
              "singletons, transversions": [r for r in single if not r["transition"]],
              "singletons, toward the consensus": [r for r in single if r["kind"] == "toward"],
              "singletons, away from the consensus": [r for r in single if r["kind"] == "away"]}
    res = {}
    for name, sel in groups.items():
        phased = [r for r in sel if r["kind"]]
        away = [r for r in sel if r["kind"] == "away"]
        res[name] = {"copies": nm.test(sel, "copies", "expected_copies"), "near": nm.test(sel, "near", "expected_near"),
                     "toward": nm.test([dict(r, is_toward=r["kind"] == "toward") for r in phased], "is_toward", "expected_toward"),
                     "away_copies": nm.test(away, "copies", "expected_copies_in_kind"),
                     "matched_context": mantel_haenszel(sel, opp)}
    by_class = {}
    for cl in CLASSES:
        sel = [r for r in rows if r["class"] == cl and r["af"] < 0.5]
        by_class[cl] = {"copies": nm.test(sel, "copies", "expected_copies"), "matched_context": mantel_haenszel(sel, opp)}
    return {"loci": len(loci), "snvs": len(rows), "singletons": len(single), "strata": len(opp), "tests": res,
            "minor_alleles_by_count": by_class}


def main(ref_dir: str, records: str, tr_records: str | None = None, features: str | None = None) -> None:
    G.init(ref_dir)
    rate, n_bg = background_rates(records)
    out = {"background_singletons": n_bg}
    ga_loci = []
    for (ch, s, e, kind), codes, a, recs in G.cached(records):
        if kind in ("GA", "TC"):
            seq = letters(codes)[s - a:e - a]
            ga_loci.append((ch, s, e, seq, [r for r in biallelic_snvs(recs) if s + LEFT <= r[0] < e - (K - 1 - LEFT)], None, None, GA_LENGTHS, 12))
    out["ga_sequences"] = analyse(ga_loci, rate)
    print("GA-sequences", json.dumps(out["ga_sequences"]["tests"]["singletons"]), flush=True)
    del ga_loci
    if tr_records and features:
        regions = []
        for (ch, s, e, kind), codes, a, recs in G.cached(tr_records):
            regions.append((ch, s, e, letters(codes)[s - a:e - a], [r for r in biallelic_snvs(recs) if s + LEFT <= r[0] < e - (K - 1 - LEFT)]))
        keys = {(ch, s, e) for ch, s, e, _, _ in regions}
        meas = {}
        with gzip.open(features, "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                k = (f[0], int(f[1]), int(f[2]))
                if k in keys:
                    d = dict(zip(head, f))
                    meas[k] = (float(d["purity"]), int(d["motif_min"]))
        tr_loci = [(ch, s, e, seq, snvs, meas[(ch, s, e)][0], meas[(ch, s, e)][1], TR_LENGTHS, min(max(12, meas[(ch, s, e)][1]), 100))
                   for ch, s, e, seq, snvs in regions]
        out["tandem_repeats"] = analyse(tr_loci, rate)
        print("tandem repeats", json.dumps(out["tandem_repeats"]["tests"]["singletons"]), flush=True)
    with open(os.path.join(nm.RESULTS, "population.json"), "w") as fh:
        json.dump(out, fh, indent=1)


def mantel_haenszel(sel: list[dict], opp: dict) -> dict:
    """As matched_context.py: strata of window, new base, purity, length and unit."""
    obs: dict = defaultdict(lambda: [0, 0])
    for r in sel:
        obs[(r["window"], r["alt"]) + r["sid"]][int(r["copies"])] += 1
    o1 = e1 = var = rr_n = rr_d = gr = 0.0
    for k, (oo0, oo1) in obs.items():
        n0, n1 = opp[k]
        n = n0 + n1
        o = oo0 + oo1
        p = n1 / n
        o1 += oo1
        e1 += o * p
        var += o * p * (1 - p)
        rr_n += oo1 * n0 / n
        rr_d += oo0 * n1 / n
        gr += n1 * n0 * o / n ** 2
    out = {"mutations": sum(sum(v) for v in obs.values()), "strata": len(obs), "copies": int(o1), "expected": round(e1, 1),
           "ratio": round(o1 / e1, 4) if e1 else None, "z": round((o1 - e1) / math.sqrt(var), 1) if var else None}
    if rr_n and rr_d:
        rr = rr_n / rr_d
        se = math.sqrt(gr / (rr_n * rr_d))
        out.update({"rate_ratio": round(rr, 4), "ci": [round(rr * math.exp(-1.96 * se), 4), round(rr * math.exp(1.96 * se), 4)]})
    return out


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(*sys.argv[1:5])
