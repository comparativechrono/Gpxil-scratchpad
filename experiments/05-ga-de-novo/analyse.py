"""What the family's de novo mutations in pure GA-sequences are, how many, and against the rest of the genome.

From the read-level verdicts (alleles.py), each candidate in a GA-sequence is put in one class:

- de novo: confirmed (the carrier's new allele in neither parent, its other allele in one), with 5 reads or more,
  in no sibling, and, when the carrier has 3 or more children with reads, passed on to one of them;
- de novo, other allele unclear: the same, but the carrier's other allele isn't seen (and it isn't passed on);
- mosaic in a parent: the new allele has a few reads in a parent (under a fifth), and is passed on, or the carrier
  has no children to tell;
- shared by siblings, not passed on (of 3 or more children), weak (under 5 reads), or not new.

A de novo allele is compared with the parental allele closest to it (edit distance) and read as one change, the
common prefix and suffix trimmed: a substitution, or bases inserted or deleted (complex when both). The parent it
most likely arose in is that allele's (when the closest alleles of the two parents are equally close, either); for
a mosaic, the parent with the reads of it. An inserted or deleted run is a tandem copy when it repeats the bases beside it; its
unit is its shortest period. A substitution is read on the GA-sequence's purine strand and tested like experiment
04's: does it make its 6 bases a 6-mer found elsewhere in the allele?

The G1 grandparents' DNA comes from cell lines only, where repeats can drift in culture, so the G2 children's
events are counted apart; the rates use the 14 trios of G3 and G4, whose parents' DNA is from blood. They are set
against the genome-wide rates of the family's paper (Porubsky et al. 2025: 74.5 de novo substitutions, and 7.4 + 65.3
indels outside and from tandem repeats, per transmission), over the autosomes' 2,756 million A/C/G/T bases (GRCh38,
experiment 04's census); GA-sequences hold 3,869,600 of them.

    python3 -I analyse.py          -> results/dnms_ga.tsv, results/analysis.json
"""

import gzip
import json
import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import alleles  # noqa: E402
import ped  # noqa: E402
from validate import period  # noqa: E402

GA_BASES = 3_869_600
PAPER = {"SNV": 74.5, "indel": 7.4 + 65.3}  # (per transmission, genome-wide; Porubsky et al. 2025)
CONFIRMED = ("confirmed", "confirmed, other allele unclear")
MIN_NEW_READS = 5
COMP = str.maketrans("ACGT", "TGCA")


def edit_distance(a: str, b: str) -> int:
    if len(a) < len(b):
        a, b = b, a
    prev = list(range(len(b) + 1))
    for i, x in enumerate(a, 1):
        cur = [i]
        for j, y in enumerate(b, 1):
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (x != y)))
        prev = cur
    return prev[-1]


def change(old: str, new: str) -> dict:
    """The one change that turns old into new, after the common prefix and suffix."""
    i = 0
    while i < min(len(old), len(new)) and old[i] == new[i]:
        i += 1
    j = 0
    while j < min(len(old), len(new)) - i and old[len(old) - 1 - j] == new[len(new) - 1 - j]:
        j += 1
    gone, came = old[i : len(old) - j], new[i : len(new) - j]
    out = {"at": i, "gone": gone, "came": came}
    if len(gone) == len(came) == 1:
        out["kind"] = "substitution"
    elif not gone or not came:
        x = came or gone
        unit = period(x)
        k, u = len(x), len(unit)
        after = i if came else i + k
        out.update({"kind": "insertion" if came else "deletion", "length": k, "unit": unit, "units": k // u,
                    "tandem_copy": (i >= k and old[i - k : i] == x) or old[after : after + k] == x,
                    "unit_beside": (i >= u and old[i - u : i] == unit) or old[after : after + u] == unit})
    else:
        out["kind"] = "complex"
        out["length"] = len(came) - len(gone)
    return out


def substitution_copy(old: str, new: str, at: int, k: int = 6) -> bool | None:
    """Whether a substitution at `at` makes the k bases around it (left of centre for even k) a k-mer of the old allele."""
    left = (k - 1) // 2
    lo = at - left
    if lo < 0 or lo + k > len(new):
        return None
    return new[lo : lo + k] in {old[x : x + k] for x in range(len(old) - k + 1)}


def poisson_ci(n: int, alpha: float = 0.05) -> tuple[float, float]:
    """The exact (Garwood) interval of a Poisson count."""
    from scipy.stats import chi2

    lo = chi2.ppf(alpha / 2, 2 * n) / 2 if n else 0.0
    return float(lo), float(chi2.ppf(1 - alpha / 2, 2 * (n + 1)) / 2)


def classify(e: dict, v: dict) -> str:
    kids_read = v["children_read"]
    passed = bool(v["children_with"])
    if v["status"] in CONFIRMED:
        if v["siblings_with"]:
            return "shared by siblings"
        if v["new_reads"] < MIN_NEW_READS:
            return "weak"
        if len(kids_read) >= 3 and not passed:
            return "not passed on"
        if v["status"] != "confirmed" and not passed:
            return "de novo, other allele unclear"
        return "de novo"
    if v["status"] == "low in a parent" and not v["siblings_with"] and v["new_reads"] >= MIN_NEW_READS:
        if len(kids_read) >= 3 and not passed:
            return "not passed on"
        return "mosaic in a parent"
    return "not new"


def describe(e: dict, v: dict, mosaic: bool = False) -> dict:
    new = v["new_allele"]
    f, m = ped.TRIOS[v["carrier"]]
    parental = []
    for p in (f, m):
        parental += [(a["seq"], p, edit_distance(a["seq"], new)) for a in alleles.call(e["reads"].get(p, {}))]
    if not parental:
        return {"kind": "unknown"}
    old, origin, dist = min(parental, key=lambda x: (x[2], x[1] != f))
    ch = change(old, new)
    tie = any(x[1] != origin and x[2] == dist for x in parental)
    ch["from"] = "either" if tie else ("father" if origin == f else "mother")
    if mosaic:
        fn, mn = v["father_new"], v["mother_new"]
        ch["from"] = "father" if fn >= 2 and mn <= 1 else ("mother" if mn >= 2 and fn <= 1 else "either")
    ch["old_length"], ch["new_length"] = len(old), len(new)
    if ch["kind"] == "substitution":
        purine = e["ga"][2] == "GA"
        o, n = (old, new) if purine else (old.translate(COMP)[::-1], new.translate(COMP)[::-1])
        at = ch["at"] if purine else len(new) - 1 - ch["at"]
        ch["purine_strand"] = f"{o[at]}>{n[at]}"
        ch["transition"] = {o[at], n[at]} in ({"A", "G"}, {"C", "T"})
        ch["copies_6mer"] = substitution_copy(o, n, at)
        ch["context"] = f"{o[max(0, at - 8) : at]}[{o[at]}>{n[at]}]{o[at + 1 : at + 9]}"
    return ch


def main() -> None:
    with gzip.open(os.path.join(ped.RESULTS, "read_alleles.json.gz"), "rt") as fh:
        events = json.load(fh)
    with open(os.path.join(ped.EXP04, "results", "census.json")) as fh:
        census = json.load(fh)
    genome = sum(c["acgt"] for c in census["chromosomes"] if c["chrom"] in ped.AUTOSOMES)
    origin = {}
    with gzip.open(os.path.join(ped.EXP04, "results", "ga-sequences.origin.tsv.gz"), "rt") as fh:
        cols = next(fh).rstrip("\n").split("\t")
        for line in fh:
            x = dict(zip(cols, line.rstrip("\n").split("\t")))
            origin[(x["chrom"], int(x["start"]))] = x
    rows = []
    for e in events:
        if e["set"] != "GA":
            continue
        for v in e["verdicts"]:
            cls = classify(e, v)
            row = {"class": cls, "chrom": e["chrom"], "ga_start": e["ga"][0], "ga_end": e["ga"][1], "ga_kind": e["ga"][2],
                   "ga_origin": origin.get((e["chrom"], e["ga"][0]), {}).get("origin", ""), "carrier": v["carrier"],
                   "generation": ped.GENERATION[v["carrier"]], "status": v["status"], "new_reads": v["new_reads"],
                   "carrier_reads": v["carrier_reads"], "father_new": v["father_new"], "father_reads": v["father_reads"],
                   "mother_new": v["mother_new"], "mother_reads": v["mother_reads"],
                   "children_with": len(v["children_with"]), "children_read": len(v["children_read"]),
                   "siblings_with": ",".join(v["siblings_with"]), "vcf": f"{e['pos']} {e['ref']}>{e['alt']}"}
            if cls in ("de novo", "de novo, other allele unclear", "mosaic in a parent"):
                ch = describe(e, v, mosaic=cls == "mosaic in a parent")
                row.update({k: ch.get(k) for k in ("kind", "from", "old_length", "new_length", "length", "gone", "came", "unit",
                                                   "units", "tandem_copy", "unit_beside", "purine_strand", "transition",
                                                   "copies_6mer", "context")})
            rows.append(row)
    cols = ["class", "chrom", "ga_start", "ga_end", "ga_kind", "ga_origin", "carrier", "generation", "status", "kind", "from",
            "old_length", "new_length", "length", "unit", "units", "tandem_copy", "unit_beside", "gone", "came", "purine_strand",
            "transition", "copies_6mer", "context", "new_reads", "carrier_reads", "father_new", "father_reads", "mother_new",
            "mother_reads", "children_with", "children_read", "siblings_with", "vcf"]
    order = {"de novo": 0, "de novo, other allele unclear": 1, "mosaic in a parent": 2, "not passed on": 3, "shared by siblings": 4,
             "weak": 5, "not new": 6}
    rows.sort(key=lambda r: (order[r["class"]], r["generation"], ped.AUTOSOMES.index(r["chrom"]), r["ga_start"]))
    with open(os.path.join(ped.RESULTS, "dnms_ga.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    blood = [r for r in rows if r["generation"] in ("G3", "G4")]
    n_trios = {"G2": 2, "G3+G4": sum(ped.GENERATION[c] in ("G3", "G4") for c in ped.TRIOS)}
    dn = [r for r in blood if r["class"] == "de novo"]
    mo = [r for r in blood if r["class"] == "mosaic in a parent"]
    kinds = Counter(("SNV" if r["kind"] == "substitution" else "indel") for r in dn)
    rates = {}
    for kind in ("SNV", "indel"):
        n = kinds.get(kind, 0)
        lo, hi = poisson_ci(n)
        t = n_trios["G3+G4"]
        per_base = n / t / GA_BASES
        genome_per_base = PAPER[kind] / genome
        expected = PAPER[kind] * t * GA_BASES / genome
        rates[kind] = {"de_novo": n, "ci": [round(lo, 2), round(hi, 2)], "per_transmission": round(n / t, 3),
                       "expected_at_genome_rate": round(expected, 2), "ratio_to_genome": round(per_base / genome_per_base, 1),
                       "ratio_ci": [round(lo / t / GA_BASES / genome_per_base, 3), round(hi / t / GA_BASES / genome_per_base, 2)],
                       "poisson_p_at_least_observed_if_160_times_fewer": (
                           round(1 - sum(math.exp(-expected / 160) * (expected / 160) ** k / math.factorial(k) for k in range(n)), 4)
                           if kind == "SNV" else None)}
    indels = [r for r in dn + mo if r["kind"] in ("insertion", "deletion", "complex")]
    # mutability by the GA-sequence's length: events (de novo and mosaic, G3 and G4) per thousand GA-sequences a
    # transmission, in bins of length
    lengths = [e - s for ch, v in ped.read_runs().items() if ch in ped.AUTOSOMES for s, e, _ in v]
    hit = Counter()
    for r in dn + mo:
        hit[r["ga_end"] - r["ga_start"]] += 1
    bins = [(51, 75), (76, 100), (101, 150), (151, 250), (251, 2000)]
    by_length = []
    for lo, hi in bins:
        loci = sum(lo <= x <= hi for x in lengths)
        ev = sum(k for x, k in hit.items() if lo <= x <= hi)
        cl, ch = poisson_ci(ev)
        t = n_trios["G3+G4"]
        by_length.append({"length": f"{lo}-{hi}" if hi < 2000 else f"{lo}+", "ga_sequences": loci, "events": ev,
                          "per_1000_loci_per_transmission": round(1000 * ev / loci / t, 3),
                          "ci": [round(1000 * cl / loci / t, 3), round(1000 * ch / loci / t, 3)]})
    # the controls: candidates elsewhere, checked the same way
    controls = Counter()
    for e in events:
        if e["set"] != "control":
            continue
        kind = "SNV" if e["kind"] in ("SNV", "MNV") else "indel"
        for v in e["verdicts"]:
            controls[f"{kind} {ped.GENERATION[v['carrier']]} | {v['status']}"] += 1
    with open(os.path.join(ped.RESULTS, "candidates.json")) as fh:
        cand = json.load(fh)
    per_trio = {k: round(sum(c[k] for c in cand["by_child"].values()) / len(cand["by_child"]), 1)
                for k in ("SNV", "MNV", "insertion", "deletion")}
    # Alu tails among the events, against their share of GA-sequences (experiment 04), one-sided
    from scipy.stats import binomtest

    alu_share = sum(x["origin"] == "Alu tail" for x in origin.values()) / len(origin)
    alu_p = round(float(binomtest(sum(r["ga_origin"] == "Alu tail" for r in dn + mo), len(dn + mo), alu_share,
                                  alternative="greater").pvalue), 3)
    out = {
        "trios": n_trios,
        "candidates": {"total": cand["candidates"], "per_trio": per_trio,
                       "in_ga_sequences": sum(1 for e in events if e["set"] == "GA" for _ in e["carriers"]),
                       "ga_events": sum(e["set"] == "GA" for e in events)},
        "controls": dict(sorted(controls.items())),
        "by_length": by_length, "ga_bases": GA_BASES, "genome_bases": genome, "paper_per_transmission": PAPER,
        "classes": {g: dict(Counter(r["class"] for r in rows if (r["generation"] == "G2") == (g == "G2"))) for g in ("G2", "G3+G4")},
        "de_novo_by_generation": dict(Counter(r["generation"] for r in rows if r["class"] == "de novo")),
        "rates_g3_g4": rates,
        "de_novo_and_mosaic_g3_g4": {
            "events": len(dn) + len(mo), "de_novo": len(dn), "mosaic_in_a_parent": len(mo),
            "kinds": dict(Counter(r["kind"] for r in dn + mo)),
            "insertions_deletions": [sum(r["kind"] == "insertion" for r in indels), sum(r["kind"] == "deletion" for r in indels)],
            "length_changes": dict(Counter((r["new_length"] - r["old_length"]) for r in indels)),
            "whole_units": sum(1 for r in indels if r["kind"] in ("insertion", "deletion") and r["unit_beside"]),
            "tandem_copies": sum(1 for r in indels if r["tandem_copy"]),
            "units": dict(Counter(r["unit"] for r in indels if r.get("unit")).most_common()),
            "one_unit_of_4": {"deletions": sum(r["kind"] == "deletion" and r["length"] == 4 and r["units"] == 1 for r in indels),
                              "insertions": sum(r["kind"] == "insertion" and r["length"] == 4 and r["units"] == 1 for r in indels)},
            "simple_by_class": {c: [sum(r["tandem_copy"] for r in indels if r["class"] == c and r["kind"] in ("insertion", "deletion")),
                                    sum(r["kind"] in ("insertion", "deletion") for r in indels if r["class"] == c)]
                                for c in ("de novo", "mosaic in a parent")},
            "from_father": sum(r["from"] == "father" for r in dn + mo), "from_mother": sum(r["from"] == "mother" for r in dn + mo),
            "from_either": sum(r["from"] == "either" for r in dn + mo),
            "ga_origins": dict(Counter(r["ga_origin"] for r in dn + mo)),
            "alu_tail_p": alu_p,
            "passed_on": [sum(r["children_with"] > 0 for r in dn + mo if r["children_read"] >= 3), sum(r["children_read"] >= 3 for r in dn + mo)]},
        "g2": {"de_novo": sum(r["class"] == "de novo" for r in rows if r["generation"] == "G2"),
               "mosaic_in_a_parent": sum(r["class"] == "mosaic in a parent" for r in rows if r["generation"] == "G2")},
        "substitutions": [r for r in rows if r.get("kind") == "substitution"],
    }
    with open(os.path.join(ped.RESULTS, "analysis.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "substitutions"}, indent=1))
    for r in out["substitutions"]:
        print({k: r[k] for k in ("class", "chrom", "ga_start", "carrier", "from", "purine_strand", "transition", "copies_6mer", "context")})


if __name__ == "__main__":
    main()
