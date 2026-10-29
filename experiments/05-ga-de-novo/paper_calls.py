"""Experiment 05 against the family's own calls: the de novo mutations of Porubsky et al. 2025 (Nature 643: 427-436),
Supplementary Table 10, in pure GA-sequences.

The paper called de novo tandem-repeat (TR) mutations with TRGT in eight G3 children, on GRCh38 and on T2T-CHM13,
and de novo substitutions (SNVs) and indels outside TRs in G2 and G3, on T2T-CHM13. Here:

- the TR mutations (GRCh38) at a GA-sequence: "in" one when at least half of the TR locus lies in pure GA- or
  TC-sequences (experiment 04's census), "touching" one when any of it does; each with the change's size in units
  of the locus's motifs, its parent of origin, and whether it passed the paper's filters (its rates use those that
  did) or is flagged recurrent;
- their rate a base against the paper's TR mutations over all of the autosomes (the same calls, so the same
  sensitivity), and against the paper's genome-wide rates per transmission as experiment 05 set its own calls
  against them;
- the five G3 children both this experiment and the paper call: which of the paper's TR mutations in GA-sequences
  were among experiment 05's candidates, and with what verdict, and which of experiment 05's were the paper's; and,
  from the family's joint calls, which of call_dnms.py's rules kept the paper's from its candidates;
- the TR mutations per base, and per repeat by its length, in GA-sequences against pure runs of the other pairs of
  bases and against RepeatMasker's other simple repeats: is a GA-sequence more mutable than another repeat as long?
- the SNVs and non-TR indels, lifted from T2T-CHM13 to GRCh38 with the T2T consortium's chain (each SNV's reference
  base checked against GRCh38's), those in GA-sequences against what the autosomes' rate gives them, a base and by
  trinucleotide context; and the de novo SVs.

    python3 -I paper_calls.py TABLE.xlsx CHAIN REF_DIR [--vcf VCF] [--rmsk RMSK]
                                                      -> results/paper_calls.json, results/paper_ga_events.tsv

TABLE.xlsx: Supplementary Table 10 of the paper (41586_2025_8922_MOESM12_ESM.xlsx). CHAIN: chm13v2-grch38.chain
(s3://human-pangenomics/T2T/CHM13/assemblies/chain/v1_nflo/). REF_DIR: the Gpxil browser's GRCh38 gpxil files.
VCF: the family's HiFi joint calls (ped.HIFI), to say why experiment 05's candidates missed the paper's. RMSK: UCSC's
RepeatMasker table of GRCh38 (hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz), to set GA-sequences
beside the other simple repeats.
"""

from __future__ import annotations

import bisect
import gzip
import json
import math
import os
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ped  # noqa: E402  (experiment 04's ga module is on the path after this)
import ga  # noqa: E402
from analyse import GA_BASES, PAPER, poisson_ci  # noqa: E402
from gpxiltools import codes as C  # noqa: E402


def _genome() -> int:
    """the autosomes' A/C/G/T bases (GRCh38; experiment 04's census, as analyse.py counts them)"""
    with open(os.path.join(ped.EXP04, "results", "census.json")) as fh:
        return sum(c["acgt"] for c in json.load(fh)["chromosomes"] if c["chrom"] in ped.AUTOSOMES)


GENOME = _genome()
FLANK_BASES = 0  #: set by main: the 400-base flanks' bases (flank_bases)
SHARED = ("NA12879", "NA12881", "NA12882", "NA12885", "NA12886")  #: G3 children both call
G3 = ("NA12879", "NA12881", "NA12882", "NA12883", "NA12884", "NA12885", "NA12886", "NA12887")
LETTER = {C.BASE_A: "A", C.BASE_C: "C", C.BASE_G: "G", C.BASE_T: "T"}
IDX = np.full(256, -1, dtype=np.int64)
for _i, _b in enumerate((C.BASE_A, C.BASE_C, C.BASE_G, C.BASE_T)):
    IDX[_b] = IDX[_b | C.SOFTMASK] = _i
COMP = str.maketrans("ACGT", "TGCA")


# ------------------------------------------------------------------ the table


def sheet(wb, name: str, header_row: int) -> list[dict]:
    rows = list(wb[name].iter_rows(values_only=True))
    head = [str(x) if x is not None else f"_{i}" for i, x in enumerate(rows[header_row])]
    return [dict(zip(head, r)) for r in rows[header_row + 1:] if r and r[0] is not None]


def as_int(x) -> int | None:
    try:
        return int(float(x))
    except (TypeError, ValueError):
        return None


# ------------------------------------------------------------------ GA-sequences


class Runs:
    """Experiment 04's pure GA- and TC-sequences of GRCh38, for overlaps."""

    def __init__(self):
        self.runs = ped.read_runs()
        self.starts = {c: [s for s, _, _ in v] for c, v in self.runs.items()}
        self.longest = max(e - s for v in self.runs.values() for s, e, _ in v)
        origin = {}
        with gzip.open(os.path.join(ped.EXP04, "results", "ga-sequences.origin.tsv.gz"), "rt") as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                origin[(r["chrom"], int(r["start"]))] = r["origin"]
        self.origin = origin

    def over(self, c: str, s: int, e: int) -> list[tuple[int, int, str]]:
        """The runs that overlap [s, e)."""
        v, st = self.runs.get(c, []), self.starts.get(c, [])
        i = bisect.bisect_left(st, s - self.longest)
        out = []
        while i < len(v) and v[i][0] < e:
            if v[i][1] > s:
                out.append(v[i])
            i += 1
        return out

    def bases(self, chroms) -> int:
        return sum(e - s for c in chroms for s, e, _ in self.runs.get(c, []))


def flank_bases(runs: Runs, flank: int = 400) -> int:
    """The bases within ``flank`` of a GA-sequence but in none, on the autosomes (the union, each base once)."""
    total = 0
    for c in ped.AUTOSOMES:
        v = runs.runs.get(c, [])
        merged: list[list[int]] = []
        for s, e, _ in v:
            a, b = max(0, s - flank), e + flank
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        total += sum(b - a for a, b in merged)
    return total - GA_BASES


def unit_change(size: int | None, motifs: str) -> str:
    """A TR mutation's size in units of the locus's motifs, the longest that divides it: '-1 x GAAA', '+2 x AG',
    '-13 x T' (only a one-base motif divides it), or '' when none does. Of motifs as long, the first listed."""
    if not size:
        return ""
    for m in sorted(dict.fromkeys(x.strip() for x in str(motifs).split(",") if x.strip()), key=lambda m: -len(m)):
        if abs(size) % len(m) == 0:
            return f"{size // len(m):+d} x {m}"
    return ""


def unit_of(x: dict) -> str:
    return x["units"].split(" x ")[1] if x["units"] else ""


LENGTHS = ((51, 75), (76, 100), (101, 150), (151, 250), (251, 10**7))


def place(events: list[dict], intervals: dict) -> list[tuple[dict, tuple[int, int]]]:
    """The events whose locus lies at least half in the intervals ({chrom: sorted [(start, end)]}), each with the
    interval it overlaps most."""
    longest = max((e - s for v in intervals.values() for s, e in v), default=0)
    starts = {c: [s for s, _ in v] for c, v in intervals.items()}
    out = []
    for x in events:
        v, st = intervals.get(x["chrom"], []), starts.get(x["chrom"], [])
        s, e = x["start"], x["end"]
        i = bisect.bisect_left(st, s - longest)
        cover, best = 0, None
        while i < len(v) and v[i][0] < e:
            ov = min(v[i][1], e) - max(v[i][0], s)
            if ov > 0:
                cover += ov
                if best is None or ov > best[0]:
                    best = (ov, v[i])
            i += 1
        if best and cover >= (e - s) / 2:
            out.append((x, best[1]))
    return out


def per_length(intervals: dict, placed: list, t: int = len(G3)) -> list[dict]:
    """Mutations per thousand intervals a transmission, by the interval's length (95 % interval)."""
    rows = []
    for lo, hi in LENGTHS:
        n = sum(lo <= e - s <= hi for v in intervals.values() for s, e in v)
        k = sum(lo <= e - s <= hi for _, (s, e) in placed)
        cl, ch = poisson_ci(k)
        rows.append({"length": f"{lo}-{hi}" if hi < LENGTHS[-1][1] else f"{lo}+", "intervals": n, "mutations": k,
                     "per_1000_per_transmission": round(1000 * k / n / t, 3) if n else None,
                     "ci": [round(1000 * cl / n / t, 3), round(1000 * ch / n / t, 3)] if n else None})
    return rows


def pure_runs(ref_dir: str, runs: Runs) -> dict:
    """Pure runs of more than 50 bases of the other pairs of bases on the autosomes, as GA-sequences are of A and G:
    of A and C (with those of G and T, the same on the other strand), of A and T, and of C and G, leaving out those
    that overlap a GA-sequence (a run of A's is in both). {pair: {chrom: sorted [(start, end)]}}"""
    pairs = {"AC/GT": ((C.BASE_A, C.BASE_C), (C.BASE_G, C.BASE_T)), "AT": ((C.BASE_A, C.BASE_T),), "CG": ((C.BASE_C, C.BASE_G),)}
    out: dict = {k: {} for k in pairs}
    for ch in ped.AUTOSOMES:
        with ga.reference(ref_dir, ch) as g:
            base = g[ch].read("ref") & 0x0F
        for k, sets in pairs.items():
            v = []
            for a, b in sets:
                s, e = ga.runs((base == a) | (base == b), ga.PURE)
                v += [(x, y) for x, y in zip(s.tolist(), e.tolist()) if not runs.over(ch, x, y)]
            out[k][ch] = sorted(v)
    return out


def simple_repeats(rmsk: str, runs: Runs) -> dict:
    """RepeatMasker's simple repeats and low-complexity sequence on the autosomes (UCSC's rmsk.txt.gz), merged where
    they overlap or touch, leaving out those that overlap a GA-sequence. {chrom: sorted [(start, end)]}"""
    iv: dict = {}
    with gzip.open(rmsk, "rt") as fh:
        for line in fh:
            f = line.split("\t", 13)
            if f[11] in ("Simple_repeat", "Low_complexity") and f[5] in ped.AUTOSOMES:
                iv.setdefault(f[5], []).append((int(f[6]), int(f[7])))
    out = {}
    for c, v in iv.items():
        merged: list[list[int]] = []
        for s, e in sorted(v):
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        out[c] = [(s, e) for s, e in merged if not runs.over(c, s, e)]
    return out


def repeats_compared(events: list[dict], runs: Runs, ref_dir: str, rmsk: str | None) -> dict:
    """The paper's passing TR mutations per base and per thousand a transmission, by length, in GA-sequences against
    pure runs of the other pairs of bases, and against RepeatMasker's other simple repeats: is a GA-sequence more
    mutable than another tandem repeat as long? A mutation counts in a set when at least half its locus is in it."""
    sel = [x for x in events if x["pass"] and x["chrom"] in ped.AUTOSOMES]
    sets = {"GA-sequences": {c: [(s, e) for s, e, _ in v] for c, v in runs.runs.items() if c in ped.AUTOSOMES}}
    for k, v in pure_runs(ref_dir, runs).items():
        sets[f"pure {k} runs"] = v
    if rmsk:
        sets["other simple repeats (RepeatMasker)"] = simple_repeats(rmsk, runs)
    out = {}
    for name, iv in sets.items():
        pool = sel if name == "GA-sequences" else [x for x in sel if x["where"] != "in"]
        placed = place(pool, iv)
        bases = sum(e - s for v in iv.values() for s, e in v)
        lo, hi = poisson_ci(len(placed))
        t = len(G3)
        out[name] = {"intervals": sum(len(v) for v in iv.values()), "bases": bases, "mutations": len(placed),
                     "per_mb_per_transmission": round(len(placed) / bases * 1e6 / t, 2),
                     "ci": [round(lo / bases * 1e6 / t, 2), round(hi / bases * 1e6 / t, 2)],
                     "by_length": per_length(iv, placed)}
    return out


# ------------------------------------------------------------------ the analyses


def tr_mutations(wb, runs: Runs, ours: list[dict]) -> tuple[dict, list[dict]]:
    rows = sheet(wb, "TRs (GRCh38)", 1)
    events = []
    for r in rows:
        c, s, e = r["chrom"], as_int(r["start"]), as_int(r["end"])
        hits = runs.over(c, s, e)
        cover = sum(min(b, e) - max(a, s) for a, b, _ in hits)
        main = max(hits, key=lambda h: min(h[1], e) - max(h[0], s)) if hits else None  # (the one it overlaps most)
        size = as_int(r["likely_denovo_size"])
        ev = {
            "sample": r["sample_id"], "chrom": c, "start": s, "end": e, "locus_length": e - s, "trid": r["trid"],
            "motifs": r["motifs"], "structure": r["simplified_motif_structure"], "size": size,
            "child_to_parent": as_int(r["child_to_parent_diff"]), "units": unit_change(size, r["motifs"]),
            "parent": r["inheritance"], "pass": r["pass_all_filters"] == "YES", "recurrent": r["recurrent_status"] == "recurrent",
            "transmitted": r["is_transmitted"], "manual": r["manual_validation_status"], "element": r["element_validation_status"],
            "ga_bases": cover, "where": "in" if hits and cover >= (e - s) / 2 else "touching" if hits else "",
            "ga": ";".join(f"{a}-{b}:{k}:{runs.origin.get((c, a), '')}" for a, b, k in hits),
            "origin": runs.origin.get((c, main[0]), "") if main else "",  # (experiment 04 classed the autosomes' only)
        }
        events.append(ev)
    auto = [x for x in events if x["chrom"] in ped.AUTOSOMES]
    out: dict = {"rows": len(events), "autosomal": len(auto), "children": dict(Counter(x["sample"] for x in events))}
    for label, keep in (("pass", lambda x: x["pass"]), ("all", lambda x: True)):
        base = [x for x in auto if keep(x)]
        out[label] = {"autosomal_tr_mutations": len(base)}
        for where in ("in", "touching"):
            sel = [x for x in base if x["where"] == "in" or (where == "touching" and x["where"])]
            n, N = len(sel), len(base)
            lo, hi = poisson_ci(n)
            per_base_ratio = (n / GA_BASES) / (N / GENOME)
            t = len(G3)
            expected_paper = PAPER["indel"] * t * GA_BASES / GENOME  # (experiment 05's yardstick)
            out[label][where] = {
                "mutations": n, "share_of_tr_mutations": round(n / N, 4), "per_transmission": round(n / t, 2),
                "ratio_to_autosomal_tr_mutations_per_base": round(per_base_ratio, 1),
                "ratio_ci": [round(lo / GA_BASES / (N / GENOME), 1), round(hi / GA_BASES / (N / GENOME), 1)],
                "expected_at_papers_genome_rate": round(expected_paper, 2),
                "ratio_to_papers_genome_rate": round(n / expected_paper, 1),
                "ratio_to_papers_genome_rate_ci": [round(lo / expected_paper, 1), round(hi / expected_paper, 1)],
                "by_child": dict(sorted(Counter(x["sample"] for x in sel).items())),
                "parent": dict(Counter(x["parent"] for x in sel)),
                "size": dict(sorted(Counter(x["size"] for x in sel).items(), key=lambda kv: (kv[0] is None, kv[0] or 0))),
                "four_bases": sum(abs(x["size"] or 0) == 4 for x in sel),
                "one_unit": sum(x["units"].startswith(("-1 x ", "+1 x ")) and len(unit_of(x)) > 1 for x in sel),
                "whole_units": sum(len(unit_of(x)) > 1 for x in sel),
                "one_base_units": sum(len(unit_of(x)) == 1 for x in sel),
                "losses": sum((x["size"] or 0) < 0 for x in sel), "gains": sum((x["size"] or 0) > 0 for x in sel),
                "transmitted": dict(Counter(x["transmitted"] for x in sel)),
                "structure": dict(Counter(x["structure"] for x in sel)),
                "origin": dict(Counter(x["origin"] for x in sel)),
            }
    # by the GA-sequence's length (the one it overlaps most), per thousand GA-sequences a transmission, as analyse.py
    # bins experiment 05's own; and how many are Alu tails, against the share of GA-sequences that are
    sel = [x for x in auto if x["pass"] and x["where"] == "in"]
    ga_runs = {c: [(s, e) for s, e, _ in v] for c, v in runs.runs.items() if c in ped.AUTOSOMES}
    by_length = per_length(ga_runs, place([x for x in auto if x["pass"]], ga_runs))
    lengths = [e - s for v in ga_runs.values() for s, e in v]
    out["pass"]["by_length"] = by_length
    alu_share = sum(1 for (c, _), o in runs.origin.items() if c in ped.AUTOSOMES and o == "Alu tail") / len(lengths)
    k = sum(x["origin"] == "Alu tail" for x in sel)
    from scipy.stats import binomtest

    # (Alu tails are a little longer than the other GA-sequences, and longer ones mutate more: their share, each
    # GA-sequence weighted by its length bin's rate, "weighted by length")
    def bin_rate(n):
        return next(b["per_1000_per_transmission"] for b, (lo, hi) in zip(by_length, LENGTHS) if lo <= n <= hi)
    weight = Counter()
    for (c, a), o in runs.origin.items():
        if c in ped.AUTOSOMES:
            e = next(e for s0, e, _ in runs.over(c, a, a + 1) if s0 == a)
            weight[o == "Alu tail"] += bin_rate(e - a)
    adjusted = weight[True] / (weight[True] + weight[False])
    out["pass"]["alu_tails"] = {"mutations": k, "of": len(sel), "share": round(k / len(sel), 3), "ga_sequences_share": round(alu_share, 3),
                                "p": float(f"{binomtest(k, len(sel), alu_share).pvalue:.3g}"),
                                "share_weighted_by_length": round(adjusted, 3),
                                "p_weighted_by_length": float(f"{binomtest(k, len(sel), adjusted).pvalue:.3g}")}
    # (the per-transmission rates of the paper's own TR calls, all and in GA-sequences, for the 160-times question's
    # counterpart: how much of a child's TR mutations GA-sequences take)
    pas = [x for x in auto if x["pass"]]
    out["pass"]["per_child_all"] = dict(sorted(Counter(x["sample"] for x in pas).items()))
    out["recurrent_loci_in_ga"] = sorted({f"{x['chrom']}:{x['start']}-{x['end']} {x['motifs']}" for x in auto if x["recurrent"] and x["where"] == "in"})
    return out, events


def concordance(events: list[dict], ours: list[dict]) -> dict:
    """The five G3 children both call: the paper's TR mutations at GA-sequences against experiment 05's candidates."""
    mine = [r for r in ours if r["carrier"] in SHARED]
    def ours_at(x):
        return [r for r in mine if r["carrier"] == x["sample"] and r["chrom"] == x["chrom"]
                and int(r["ga_start"]) < x["end"] and int(r["ga_end"]) > x["start"]]
    paper = [x for x in events if x["sample"] in SHARED and x["where"] == "in" and x["chrom"] in ped.AUTOSOMES]
    rows = []
    for x in paper:
        m = ours_at(x)
        x["experiment05"] = m[0]["class"] if m else "not a candidate"
        rows.append(x)
    found = Counter((x["pass"], x["experiment05"]) for x in rows)
    ours_events = [r for r in mine if r["class"] in ("de novo", "de novo, other allele unclear", "mosaic in a parent")]
    in_paper, lacking = Counter(), []
    for r in ours_events:
        hit = [x for x in events if x["sample"] == r["carrier"] and x["chrom"] == r["chrom"]
               and x["start"] < int(r["ga_end"]) and x["end"] > int(r["ga_start"])]
        key = "paper, passing" if any(x["pass"] for x in hit) else "paper, failing its filters" if hit else "not in the paper"
        in_paper[(r["class"], key)] += 1
        if key != "paper, passing" and r["class"] != "mosaic in a parent":
            lacking.append({"carrier": r["carrier"], "ga": f"{r['chrom']}:{int(r['ga_start']) + 1}-{r['ga_end']}", "class": r["class"],
                            "paper": key, "kind": r["kind"], "length": f"{r['old_length']} -> {r['new_length']}",
                            "unit": r["unit"], "children_with": f"{r['children_with']} of {r['children_read']}"})
    return {
        "children": list(SHARED),
        "paper_in_ga": len(rows), "paper_in_ga_pass": sum(x["pass"] for x in rows),
        "paper_by_experiment05_verdict": {f"{'pass' if p else 'fail'}: {v}": n for (p, v), n in sorted(found.items())},
        "experiment05_events": len(ours_events),
        "experiment05_in_paper": {f"{c}: {k}": n for (c, k), n in sorted(in_paper.items())},
        "experiment05_de_novo_not_passing_in_paper": lacking,
    }


RULES = ("the child or a parent isn't called", "a GQ under 20 in the trio", "a depth under 10 in the trio",
         "a depth over 2.5 times the median in the trio", "a parent has a read of it",
         "under a quarter of the child's reads", "another relative is called with it")


def broken_rules(g, gq, dp, ad, k: int, ci: int, fi: int, mi: int, max_dp, outside) -> list[str]:
    """Which of call_dnms.py's rules an allele k, new in a child (ci) to its parents (fi, mi), breaks, in its order."""
    from call_dnms import MIN_AB, MIN_DP, MIN_GQ

    trio = [ci, fi, mi]
    out = []
    if (g[[fi, mi]] < 0).any() or (g[ci] < 0).any():
        out.append(RULES[0])
    if min(gq[trio]) < MIN_GQ:
        out.append(RULES[1])
    if min(dp[trio]) < MIN_DP:
        out.append(RULES[2])
    if (dp[trio] > max_dp[trio]).any():
        out.append(RULES[3])
    if k >= ad.shape[1] or ad[ci, k] < 0 or ad[fi, k] != 0 or ad[mi, k] != 0:
        out.append(RULES[4])
    tot = int(ad[ci].clip(0).sum())
    if k < ad.shape[1] and (tot == 0 or ad[ci, k] / tot < MIN_AB):
        out.append(RULES[5])
    if ((g[:, 0] == k) | (g[:, 1] == k))[outside].any():
        out.append(RULES[6])
    return out


def why_missed(events: list[dict], vcf_path: str) -> dict:
    """For the paper's passing TR mutations in GA-sequences in the five shared children: what the HiFi joint calls
    hold over the locus, and which of call_dnms.py's rules kept each from experiment 05's candidates. A record holds
    a new allele when the child is called with an allele neither parent is; of the locus's new alleles, the one that
    breaks the fewest rules speaks for it, and the first rule it breaks, in call_dnms.py's order, is its reason. The
    calls' alleles are written against the reference, so a new allele is of the paper's size when it differs by that
    from an allele a parent is called with in the same record (or from the reference, when neither parent is called)."""
    from cyvcf2 import VCF

    from call_dnms import MAX_DP_FACTOR

    with open(os.path.join(ped.RESULTS, "candidates.json")) as fh:
        medians = json.load(fh)["median_depths"]
    vcf = VCF(vcf_path)
    ix = {x: i for i, x in enumerate(vcf.samples)}
    max_dp = np.array([MAX_DP_FACTOR * medians[x] for x in vcf.samples])
    first, every, rows = Counter(), Counter(), []
    for x in events:
        if not (x["sample"] in SHARED and x["where"] == "in" and x["pass"] and x["chrom"] in ped.AUTOSOMES):
            continue
        child = x["sample"]
        ci, fi, mi = (ix[p] for p in (child, *ped.TRIOS[child]))
        keep = {child} | ped.descendants(child)
        for sib in ped.siblings(child):
            keep |= {sib} | ped.descendants(sib)
        outside = np.array([p not in keep for p in vcf.samples])
        best, of_size = None, False
        recs = list(vcf(f"{x['chrom']}:{x['start'] - 5}-{x['end'] + 5}"))
        for v in recs:
            g = v.genotype.array()[:, :2]
            ad, gq, dp = v.format("AD"), v.format("GQ"), v.format("DP")
            if ad is None or gq is None or dp is None:
                continue
            lens = [0] + [len(a) - len(v.REF) for a in v.ALT]
            parents = {a for a in (*g[fi], *g[mi]) if a >= 0} or {0}
            for k in range(1, len(v.ALT) + 1):
                has = (g[:, 0] == k) | (g[:, 1] == k)
                if not has[ci] or has[fi] or has[mi]:
                    continue
                of_size |= any(lens[k] - lens[a] == x["size"] for a in parents)
                broken = broken_rules(g, gq[:, 0], dp[:, 0], ad, k, ci, fi, mi, max_dp, outside)
                if best is None or len(broken) < len(best[1]):
                    best = ({"at": f"{v.CHROM}:{v.POS}", "change": lens[k]}, broken)
        reason = ("no record over the locus" if not recs else "no new allele in the calls" if best is None
                  else best[1][0] if best[1] else "none: a candidate")
        first[reason] += 1
        for r in (best[1] if best else []):
            every[r] += 1
        rows.append({"sample": child, "locus": f"{x['chrom']}:{x['start']}-{x['end']}", "size": x["size"],
                     "experiment05": x["experiment05"], "record": best[0] if best else None, "reason": reason,
                     "broken": best[1] if best else [], "a_new_allele_of_the_papers_size": of_size})
    return {"mutations": len(rows), "first_rule_broken": dict(first.most_common()),
            "any_rule_broken": dict(every.most_common()),
            "with_a_new_allele_of_the_papers_size": sum(r["a_new_allele_of_the_papers_size"] for r in rows if r["reason"] != "none: a candidate"),
            "loci": rows}


def lift_and_read(rows: list[dict], pos_key: str, chain: str, ref_dir: str, kind: str) -> list[dict]:
    """Each row's T2T-CHM13 position on GRCh38 (0-based), its GRCh38 base and trinucleotide, and the check of its
    reference allele (an SNV's base, or an indel's first base)."""
    from pyliftover import LiftOver

    lo = LiftOver(chain)
    out = []
    by_chrom: dict = {}
    for r in rows:
        pos = as_int(r[pos_key])
        hit = lo.convert_coordinate(r["chromosome"], pos - 1) if pos else None
        x = {"sample": r["sample"], "chm13": f"{r['chromosome']}:{pos}", "kind": kind,
             "ref": str(r.get("reference allele") or r.get("ref") or ""), "alt": str(r.get("alternate allele") or r.get("alt") or ""),
             "origin": r.get("variant origin"), "parent": r.get("inheritance"), "chrom": None, "pos": None, "strand": None}
        if hit and len(hit) == 1:
            x["chrom"], x["pos"], x["strand"] = hit[0][0], int(hit[0][1]), hit[0][2]
            by_chrom.setdefault(x["chrom"], []).append(x)
        out.append(x)
    for ch, xs in by_chrom.items():
        if ch not in ga.CHROMS:
            continue
        with ga.reference(ref_dir, ch) as g:
            for x in xs:
                three = g[ch].read("ref", max(0, x["pos"] - 1), x["pos"] + 2) & 0x0F
                x["grch38"] = "".join(LETTER.get(int(b), "N") for b in three)
                base = x["grch38"][1] if len(x["grch38"]) == 3 else ""
                want = x["ref"][:1].upper()
                x["ref_matches"] = base == (want if x["strand"] == "+" else want.translate(COMP))
    return out


def contexts_of_genome(ref_dir: str, runs: Runs, chunk: int = 1 << 24) -> tuple[np.ndarray, np.ndarray]:
    """The autosomes' trinucleotides (64), all and with the middle base in a GA-sequence, a chunk at a time."""
    idx = IDX.astype(np.int8)
    allc = np.zeros(64, dtype=np.int64)
    gac = np.zeros(64, dtype=np.int64)
    for ch in ped.AUTOSOMES:
        with ga.reference(ref_dir, ch) as g:
            c = g[ch]
            for a in range(0, c.length - 2, chunk):
                b = min(c.length, a + chunk + 2)
                x = idx[c.read("ref", a, b)].astype(np.int16)
                ok = (x[:-2] >= 0) & (x[1:-1] >= 0) & (x[2:] >= 0)
                code = x[:-2] * 16 + x[1:-1] * 4 + x[2:]
                allc += np.bincount(code[ok], minlength=64)
                hits = runs.over(ch, a + 1, b - 1)  # (the middle bases: a + 1 to b - 2)
                inside = ga.mask_of(b - a - 2, [max(0, s - a - 1) for s, _, _ in hits], [min(b - a - 2, e - a - 1) for _, e, _ in hits])
                gac += np.bincount(code[ok & inside], minlength=64)
    return allc, gac


def pyrimidine(tri: str) -> str:
    return tri if tri[1] in "CT" else tri[::-1].translate(COMP)


def point_mutations(wb, runs: Runs, chain: str, ref_dir: str) -> dict:
    snv = lift_and_read(sheet(wb, "SNVs", 1), "T2T-CHM13v2.0 position", chain, ref_dir, "SNV")
    ind = lift_and_read(sheet(wb, "indels (non-TR)", 1), "T2T-CHM13v2.0 position", chain, ref_dir, "indel")
    allc, gac = contexts_of_genome(ref_dir, runs)
    names = ["ACGT"[i // 16] + "ACGT"[(i // 4) % 4] + "ACGT"[i % 4] for i in range(64)]
    genome_ctx, ga_ctx = Counter(), Counter()
    for i, n in enumerate(names):
        genome_ctx[pyrimidine(n)] += int(allc[i])
        ga_ctx[pyrimidine(n)] += int(gac[i])
    out = {"genome_trinucleotides": int(allc.sum()), "ga_trinucleotides": int(gac.sum())}
    for name, xs in (("SNV", snv), ("indel", ind)):
        lifted = [x for x in xs if x["chrom"]]
        auto = [x for x in lifted if x["chrom"] in ped.AUTOSOMES]
        good = [x for x in auto if x.get("ref_matches")] if name == "SNV" else auto
        for x in good:
            x["in_ga"] = bool(runs.over(x["chrom"], x["pos"], x["pos"] + 1))
            x["near_ga"] = bool(runs.over(x["chrom"], x["pos"] - 400, x["pos"] + 401))
            x["in_flank"] = x["near_ga"] and not x["in_ga"]
        res = {"rows": len(xs), "lifted": len(lifted), "autosomal": len(auto), "reference_base_matches": len(good) if name == "SNV" else None,
               "children": sorted({x["sample"] for x in xs})}
        for origin in ("germline", "postzygotic", "all"):
            sel = [x for x in good if origin == "all" or x["origin"] == origin]
            n_in = sum(x["in_ga"] for x in sel)
            expected = len(sel) * GA_BASES / GENOME
            flank_expected = len(sel) * FLANK_BASES / GENOME
            n_flank = sum(x["in_flank"] for x in sel)
            d = {"mutations": len(sel), "in_ga": n_in, "expected_per_base": round(expected, 2),
                 "within_400_bases": sum(x["near_ga"] for x in sel),
                 "in_flanks": n_flank, "expected_in_flanks_per_base": round(flank_expected, 2),
                 "p_flanks_at_most_observed": round(sum(math.exp(-flank_expected) * flank_expected ** k / math.factorial(k) for k in range(n_flank + 1)), 4),
                 "in_ga_list": [f"{x['sample']} {x['chrom']}:{x['pos'] + 1} {x['ref']}>{x['alt']} ({x['grch38']})" for x in sel if x["in_ga"]]}
            if name == "SNV":
                rate = Counter()
                for x in sel:
                    rate[pyrimidine(x["grch38"])] += 1
                d["expected_by_context"] = round(sum(rate[k] / genome_ctx[k] * ga_ctx[k] for k in rate if genome_ctx[k]), 2)
                lam = d["expected_by_context"]
                d["p_none_if_160_times_fewer"] = round(math.exp(-lam / 160), 3)
                d["p_none_at_genome_rate"] = round(math.exp(-lam), 3)
            res[origin] = d
        out[name] = res
    return out


def structural(wb, runs: Runs) -> dict:
    rows = sheet(wb, "SVs (GRCh38)", 0 + 1)
    hits = []
    for r in rows:
        c, s, e = r["chromosome"], as_int(r["start"]), as_int(r["end"])
        if s is None or e is None:
            continue
        if runs.over(c, s, e):
            hits.append(f"{r['sample']} {c}:{s}-{e} {r['SV type']} {r['SV size']} ({r['Source caller']}, {r['Mechanism']})")
    return {"rows": len(rows), "at_ga": hits}


def main(table: str, chain: str, ref_dir: str, vcf_path: str | None = None, rmsk: str | None = None) -> None:
    import openpyxl

    global FLANK_BASES
    wb = openpyxl.load_workbook(table, read_only=True, data_only=True)
    runs = Runs()
    assert runs.bases(ped.AUTOSOMES) == GA_BASES, runs.bases(ped.AUTOSOMES)
    FLANK_BASES = flank_bases(runs)
    with open(os.path.join(ped.RESULTS, "dnms_ga.tsv")) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ours = [dict(zip(head, line.rstrip("\n").split("\t"))) for line in fh]
    tr, events = tr_mutations(wb, runs, ours)
    conc = concordance(events, ours)
    if vcf_path:
        conc["why_not_candidates"] = why_missed(events, vcf_path)
    tr["compared"] = repeats_compared(events, runs, ref_dir, rmsk)
    points = point_mutations(wb, runs, chain, ref_dir)
    svs = structural(wb, runs)
    result = {"source": "Porubsky et al. 2025, Nature 643: 427-436, Supplementary Table 10", "ga_bases": GA_BASES, "genome_bases": GENOME,
              "tr": tr, "concordance": conc, "points": points, "svs": svs}
    with open(os.path.join(ped.RESULTS, "paper_calls.json"), "w") as fh:
        json.dump(result, fh, indent=1)
    cols = ["sample", "chrom", "start", "end", "locus_length", "motifs", "structure", "size", "units", "parent", "pass", "recurrent",
            "transmitted", "where", "ga_bases", "ga", "origin", "experiment05"]
    with open(os.path.join(ped.RESULTS, "paper_ga_events.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in sorted((x for x in events if x["where"]), key=lambda x: (x["sample"], ga.CHROMS.index(x["chrom"]), x["start"])):
            fh.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "points"}, indent=1)[:6000])
    print(json.dumps(points, indent=1)[:4000])


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("table")
    ap.add_argument("chain")
    ap.add_argument("ref_dir")
    ap.add_argument("--vcf", help="the family's HiFi joint calls")
    ap.add_argument("--rmsk", help="UCSC's RepeatMasker table of GRCh38, rmsk.txt.gz")
    a = ap.parse_args()
    main(a.table, a.chain, a.ref_dir, a.vcf, a.rmsk)
