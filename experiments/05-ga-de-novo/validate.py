"""Each candidate de novo mutation in a pure GA-sequence, checked in the family's HiFi reads, and a random sample of the
candidates elsewhere, checked the same way, to see how often a candidate is real.

For each candidate, a window over the whole repeat it sits in (the GA-sequence, or the run of the indel's unit around
it), widened until each end has 15 bases found only once nearby (anchors); then, for the carrier, its parents,
siblings, children and grandparents, the sequence between the anchors of every read that spans the window (MAPQ >= 1,
primary), counted by sequence: each person's alleles, read whole, without trusting the alignment inside the repeat.
alleles.py then calls each person's alleles and the verdict on each candidate.

    python3 -I validate.py REF_DIR [WORKERS [CONTROLS]]   -> results/read_alleles.json.gz, results/validated.tsv

REF_DIR: the Gpxil browser's GRCh38 gpxil files (for the windows). The reads are read from the family's public BAM
files a window at a time; run it in a directory of its own, where htslib keeps the indexes it downloads (about
500 MB). CONTROLS: how many candidates outside GA-sequences to check (default 150).
"""

import bisect
import gzip
import json
import os
import random
import sys
import time
from collections import Counter
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import alleles  # noqa: E402
import ped  # noqa: E402

FLANK, ANCHOR, REACH = 10, 15, 150


def bam_url(sample: str) -> str:
    return f"{ped.BUCKET}/data/hifi/mapped/GRCh38/{sample}.GRCh38.haplotagged.bam"


def load_candidates() -> list[dict]:
    rows = []
    with gzip.open(os.path.join(ped.RESULTS, "candidates.tsv.gz"), "rt") as fh:
        cols = next(fh).rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(cols, line.rstrip("\n").split("\t"))))
    return rows


def span(x: dict) -> tuple[int, int]:
    """0-based [a, b) of the reference bases a candidate changes, or between which it inserts."""
    p0 = int(x["pos"]) - 1
    if x["kind"] in ("SNV", "MNV"):
        return p0, p0 + len(x["ref"])
    if x["kind"] == "deletion":
        return p0 + 1, p0 + len(x["ref"])
    return p0, p0 + 2


class GA:
    def __init__(self):
        self.runs = ped.read_runs()
        self.starts = {c: [r[0] for r in v] for c, v in self.runs.items()}

    def hit(self, c: str, a: int, b: int):
        v = self.runs.get(c, [])
        i = bisect.bisect_right(self.starts.get(c, []), b) - 1
        for j in (i, i - 1):
            if 0 <= j < len(v) and v[j][0] < b and a < v[j][1]:
                return v[j]
        return None


def period(unit: str) -> str:
    """The shortest unit that a sequence is a whole number of copies of (AAAGAAAG: AAAG)."""
    for k in range(1, len(unit) + 1):
        if len(unit) % k == 0 and unit[:k] * (len(unit) // k) == unit:
            return unit[:k]
    return unit


def repeat_window(seq: str, offset: int, a: int, b: int, unit: str) -> tuple[int, int]:
    """Widen [a, b) over the run of a unit (the indel's inserted or deleted bases, at their shortest period) around
    it, in a reference string that starts at offset."""
    if not unit:
        return a, b
    k = len(period(unit))
    lo, hi = a - offset, b - offset
    # extend a base at a time while the sequence keeps its period k
    while lo - 1 >= 0 and lo - 1 + k < len(seq) and seq[lo - 1] == seq[lo - 1 + k]:
        lo -= 1
    while hi < len(seq) and hi - k >= 0 and seq[hi] == seq[hi - k]:
        hi += 1
    return lo + offset, hi + offset


def windows(ref_dir: str, n_controls: int) -> list[dict]:
    """The events to check: candidates in GA-sequences, and a random sample of the others; one per site and allele."""
    import ga  # (experiment 04's)

    rows = load_candidates()
    gaseq = GA()
    events: dict = {}
    for x in rows:
        a, b = span(x)
        hit = gaseq.hit(x["chrom"], a, b)
        key = (x["chrom"], x["pos"], x["ref"], x["alt"])
        e = events.setdefault(key, {"chrom": x["chrom"], "pos": int(x["pos"]), "ref": x["ref"], "alt": x["alt"], "kind": x["kind"],
                                    "length_change": int(x["length_change"]), "carriers": [], "ga": list(hit) if hit else None})
        e["carriers"].append(x["child"])
    ga_events = [e for e in events.values() if e["ga"]]
    rng = random.Random(5)
    others = [e for e in events.values() if not e["ga"]]
    controls = rng.sample(others, min(n_controls, len(others)))
    import numpy as np

    lut = np.full(256, ord("N"), dtype=np.uint8)
    for code, letter in ga.LETTER.items():
        lut[code] = ord(letter)
    order = {c: i for i, c in enumerate(ped.AUTOSOMES)}
    out = []
    refs: dict = {}
    for e in sorted(ga_events + controls, key=lambda e: (order[e["chrom"]], e["pos"])):
        c = e["chrom"]
        if c not in refs:
            with ga.reference(ref_dir, c) as g:
                refs.clear()
                refs[c] = lut[g[c].read("ref") & 0x0F].tobytes().decode()
        seq = refs[c]
        a, b = span(e)
        if e["kind"] in ("insertion", "deletion"):
            unit = e["alt"][1:] if e["kind"] == "insertion" else e["ref"][1:]
            a, b = repeat_window(seq, 0, max(a - 1, 0), b, unit)
        if e["ga"]:
            a, b = min(a, e["ga"][0]), max(b, e["ga"][1])
        a, b = anchored(seq, a, b)
        e["window"] = [a, b]
        e["window_ref"] = seq[a:b]
        e["set"] = "GA" if e["ga"] else "control"
        out.append(e)
    return out


def anchored(seq: str, a: int, b: int) -> tuple[int, int]:
    """Widen [a, b) by FLANK bases and an anchor's length, then base by base (up to REACH) until the ANCHOR bases at
    each end, outside [a - FLANK, b + FLANK), occur only once within 200 bases."""
    lo, hi = max(0, a - FLANK - ANCHOR), min(len(seq), b + FLANK + ANCHOR)
    while lo > max(0, a - REACH) and seq.count(seq[lo : lo + ANCHOR], max(0, lo - 200), min(len(seq), hi + 200)) != 1:
        lo -= 1
    while hi < min(len(seq), b + REACH) and seq.count(seq[hi - ANCHOR : hi], max(0, lo - 200), min(len(seq), hi + 200)) != 1:
        hi += 1
    return lo, hi


def people(carrier: str) -> list[str]:
    f, m = ped.TRIOS[carrier]
    out = [carrier, f, m] + ped.siblings(carrier) + ped.children(carrier)
    for p in (f, m):
        if p in ped.TRIOS:
            out += list(ped.TRIOS[p])
    return list(dict.fromkeys(out))


def read_alleles(bam, ch: str, a: int, b: int, left: str, right: str) -> dict:
    """The sequences of the reads that span [a, b), from the left anchor to the right one, counted."""
    counts: Counter = Counter()
    for r in bam.fetch(ch, a, b):
        if r.is_secondary or r.is_supplementary or r.is_unmapped or r.mapping_quality < 1:
            continue
        if r.reference_start > a or r.reference_end is None or r.reference_end < b or r.query_sequence is None:
            continue
        qa = qb = None
        for q, rp in r.get_aligned_pairs(matches_only=True):
            if qa is None and rp >= a:
                qa = q
            if rp >= b - 1:
                qb = q
                break
        if qa is None or qb is None:
            continue
        s = r.query_sequence
        i = s.find(left, max(0, qa - 60), qa + 60 + len(left))
        j = s.rfind(right, max(0, qb - 60 - len(right)), qb + 60)
        if i < 0 or j < 0 or j + len(right) <= i:
            continue
        counts[s[i : j + len(right)]] += 1
    return dict(counts)


def job(args):
    """One person's reads in every window that needs them."""
    sample, todo = args
    import pysam

    url = bam_url(sample)
    out = {}
    for attempt in range(6):
        try:
            bam = pysam.AlignmentFile(url, index_filename=url + ".bai")
            break
        except (OSError, ValueError):
            if attempt == 5:
                raise
            time.sleep(3 + 5 * attempt)
    for key, ch, a, b, left, right in sorted(todo, key=lambda t: (t[1], t[2])):
        for attempt in range(6):
            try:
                out[key] = read_alleles(bam, ch, a, b, left, right)
                break
            except (OSError, ValueError):
                if attempt == 5:
                    raise
                time.sleep(3 + 5 * attempt)
                bam = pysam.AlignmentFile(url, index_filename=url + ".bai")
    return sample, out


def main(ref_dir: str, workers: int = 4, n_controls: int = 150) -> None:
    t0 = time.time()
    events = windows(ref_dir, n_controls)
    print(f"{sum(e['set'] == 'GA' for e in events)} events in GA-sequences, {sum(e['set'] == 'control' for e in events)} controls "
          f"[{time.time() - t0:.0f}s]", flush=True)
    todo: dict = {}
    for i, e in enumerate(events):
        e["id"] = i
        need = set()
        for c in e["carriers"]:
            need |= set(people(c))
        w = e["window_ref"]
        for s in need:
            todo.setdefault(s, []).append((i, e["chrom"], e["window"][0], e["window"][1], w[:ANCHOR], w[-ANCHOR:]))
    reads: dict = {}
    with Pool(workers) as pool:
        for sample, out in pool.imap_unordered(job, list(todo.items())):
            for i, counts in out.items():
                reads.setdefault(i, {})[sample] = counts
            print(f"{sample}: {len(out)} windows [{time.time() - t0:.0f}s]", flush=True)
    for e in events:
        e["reads"] = reads.get(e["id"], {})
    tab = alleles.write(events)
    with gzip.open(os.path.join(ped.RESULTS, "read_alleles.json.gz"), "wt") as fh:
        json.dump(events, fh)
    for k in sorted(tab):
        print(k, tab[k])


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 4, int(sys.argv[3]) if len(sys.argv) > 3 else 150)
