"""Each person's alleles in a window, from the sequences of the reads that span it (validate.py), and the verdict on a
candidate de novo mutation: is the carrier's new allele in neither parent?

Alleles are called by length first, as repeats change by whole units and reads mostly err by a base: the lengths
seen are taken in order of how many reads have them, each becoming an allele unless one already taken is a base
away, in which case its reads join that one. An allele needs 3 reads and 15 % of the person's reads; its sequence is
the commonest read of its exact length, and when two sequences of that length both have 3 reads or more and a third
of its reads, they are two alleles (a substitution apart).

A carrier's allele is inherited when a parent has an allele of the same sequence, or one a base longer or shorter in
a run of one base (HiFi's commonest error: so a new allele that is a homopolymer a base longer or shorter is missed);
it is new when neither does, with at most one read of its sequence in each parent (or, for a new length, at most one
read of that length). With 2 reads or more of it in a parent but under a fifth of the parent's reads, it is "low in a
parent": a parent's mosaic, or noise, and not counted as new.

    python3 -I alleles.py           -> results/validated.tsv (again, from results/read_alleles.json.gz)
"""

import gzip
import json
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ped  # noqa: E402

MIN_READS, MIN_FRACTION, SPLIT, MAX_PARENT, MIN_PARENT_READS, LOW = 3, 0.15, 1 / 3, 1, 8, 0.2


def call(counts: dict) -> list[dict]:
    """A person's alleles from read sequences counted: [{length, seq, reads, exact}], commonest first."""
    n = sum(counts.values())
    by_len = Counter()
    for s, k in counts.items():
        by_len[len(s)] += k
    peaks: dict = {}
    for length, k in sorted(by_len.items(), key=lambda lk: (-lk[1], lk[0])):
        near = [p for p in peaks if abs(p - length) <= 1]
        if near:
            peaks[near[0]] += k
        else:
            peaks[length] = k
    out = []
    for length, k in sorted(peaks.items(), key=lambda lk: -lk[1]):
        if k < MIN_READS or k < MIN_FRACTION * n:
            continue
        exact = Counter({s: c for s, c in counts.items() if len(s) == length})
        seqs = exact.most_common()
        if not seqs:
            continue
        top = [(s, c) for s, c in seqs if c >= MIN_READS and c >= SPLIT * sum(exact.values())]
        if len(top) >= 2:
            for s, c in top[:2]:
                out.append({"length": length, "seq": s, "reads": c, "exact": c})
        else:
            out.append({"length": length, "seq": seqs[0][0], "reads": k, "exact": seqs[0][1]})
    return out


def homopolymer_step(a: str, b: str) -> bool:
    """Whether b is a with one base more or less in a run of that base (and nothing else different)."""
    if abs(len(a) - len(b)) != 1:
        return False
    long, short = (a, b) if len(a) > len(b) else (b, a)
    i = 0
    while i < len(short) and long[i] == short[i]:
        i += 1
    if long[:i] + long[i + 1 :] != short:
        return False
    x = long[i]
    return (i > 0 and long[i - 1] == x) or (i + 1 < len(long) and long[i + 1] == x)


def same(a: str, b: str) -> bool:
    return a == b or homopolymer_step(a, b)


def verdict(reads: dict, carrier: str) -> dict:
    """The carrier's alleles against its parents', and whether its children and siblings have the new one."""
    f, m = ped.TRIOS[carrier]
    kid, pf, pm = reads.get(carrier, {}), reads.get(f, {}), reads.get(m, {})
    n, nf, nm = sum(kid.values()), sum(pf.values()), sum(pm.values())
    a_kid, a_f, a_m = call(kid), call(pf), call(pm)
    parent_seqs = {a["seq"] for a in a_f + a_m}
    parent_lens = {a["length"] for a in a_f + a_m}
    sure_lens = {a["length"] for a in a_f + a_m if a["exact"] >= MIN_READS}

    def support(counts: dict, a: dict) -> int:
        """Reads of a parent for a carrier's allele: of its sequence, or, for a length no parent allele has, of its
        length."""
        if a["length"] in parent_lens:
            return counts.get(a["seq"], 0)
        return sum(c for s, c in counts.items() if len(s) == a["length"])

    new, low, old = [], [], []
    for a in a_kid:
        if any(same(a["seq"], s) for s in parent_seqs):
            old.append(a)
            continue
        if a["length"] in parent_lens and not (a["exact"] >= MIN_READS and a["length"] in sure_lens):
            old.append(a)  # (a parent has this length, and the reads are too noisy to tell sequences apart)
            continue
        sf, sm = support(pf, a), support(pm, a)
        if sf <= MAX_PARENT and sm <= MAX_PARENT:
            new.append(a)
        elif (sf <= MAX_PARENT or sf < LOW * nf) and (sm <= MAX_PARENT or sm < LOW * nm):
            low.append(a)
        else:
            old.append(a)
    if n < 2 * MIN_READS or nf < MIN_PARENT_READS or nm < MIN_PARENT_READS:
        status = "too few reads"
    elif not a_kid:
        status = "unclear"
    elif len(new) == 1 and old:
        status = "confirmed"
    elif len(new) == 1:
        status = "confirmed, other allele unclear"
    elif len(new) > 1:
        status = "two new alleles"
    elif low:
        status = "low in a parent"
    else:
        status = "in a parent"
    nu = new[0] if new else (low[0] if low else None)

    def has(sample: str) -> bool:
        if not nu:
            return False
        return any(same(a["seq"], nu["seq"]) or (a["length"] == nu["length"] and nu["length"] not in parent_lens)
                   for a in call(reads.get(sample, {})))

    kids = ped.children(carrier)
    return {"carrier": carrier, "status": status, "new_allele": nu["seq"] if nu else "", "new_reads": nu["reads"] if nu else 0,
            "carrier_reads": n, "carrier_alleles": [(a["length"], a["reads"]) for a in a_kid],
            "other_allele": old[0]["seq"] if old else "", "father_reads": nf, "father_alleles": [(a["length"], a["reads"]) for a in a_f],
            "father_new": support(pf, nu) if nu else 0, "mother_reads": nm, "mother_alleles": [(a["length"], a["reads"]) for a in a_m],
            "mother_new": support(pm, nu) if nu else 0,
            "children_with": [k for k in kids if has(k)], "children_read": [k for k in kids if sum(reads.get(k, {}).values()) >= 2 * MIN_READS],
            "siblings_with": [s for s in ped.siblings(carrier) if has(s)]}


COLS = ["set", "chrom", "pos", "ref", "alt", "kind", "length_change", "carrier", "status", "carrier_reads", "carrier_alleles",
        "new_reads", "father_reads", "father_alleles", "father_new", "mother_reads", "mother_alleles", "mother_new",
        "children_with", "children_read", "siblings_with", "window", "ga", "new_allele", "other_allele"]


def write(events: list) -> Counter:
    rows = []
    for e in events:
        e["verdicts"] = [verdict(e["reads"], c) for c in e["carriers"]]
        rows += [e | v for v in e["verdicts"]]
    with open(os.path.join(ped.RESULTS, "validated.tsv"), "w") as fh:
        fh.write("\t".join(COLS) + "\n")
        for x in rows:
            vals = []
            for c in COLS:
                v = x.get(c)
                vals.append(";".join(map(str, v)) if isinstance(v, (list, tuple)) else ("" if v is None else str(v)))
            fh.write("\t".join(vals) + "\n")
    return Counter((x["set"], x["status"]) for x in rows)


def main() -> None:
    path = os.path.join(ped.RESULTS, "read_alleles.json.gz")
    with gzip.open(path, "rt") as fh:
        events = json.load(fh)
    tab = write(events)
    with gzip.open(path, "wt") as fh:
        json.dump(events, fh)
    for k in sorted(tab):
        print(k, tab[k])


if __name__ == "__main__":
    main()
