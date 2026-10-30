"""gnomAD v4.1 genomes' PASS records for a random sample of the catalogue's tandem repeats, outside the GA-sequences.

The loci: those of experiment 06's per-locus measures on the autosomes with purity 0.6 or more and 20 bases or more,
not GA-sequences (experiment 04 read gnomAD for those), drawn at random (seed 7). The records are read from the web a
locus at a time, as experiment 04 read them, and written in experiment 04's local-copy format (a `#region` line a locus,
kind `TR`, then its records), for population.py.

    CURL_CA_BUNDLE=... python3 -I gnomad_tr.py FEATURES OUT.tsv.gz [N [WORKERS]]
"""

from __future__ import annotations

import gzip
import os
import random
import sys
import time
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "04-ga-sequences"))

import gnomad as G  # noqa: E402  (experiment 04)

AUTOSOMES = {f"chr{i}" for i in range(1, 23)}


def job(tasks):
    out = []
    for ch, s, e in sorted(tasks, key=lambda t: t[1]):
        out.append((ch, s, e, G.fetch(ch, s, e)))
    G._tabix.clear()
    return out


def main(features: str, out: str, n: int = 20000, workers: int = 8) -> None:
    loci = []
    with gzip.open(features, "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(head, line.rstrip("\n").split("\t")))
            if f["chrom"] in AUTOSOMES and float(f["purity"]) >= 0.6 and int(f["length"]) >= 20 and f["class"] != "GA-sequence":
                loci.append((f["chrom"], int(f["start"]), int(f["end"])))
    rng = random.Random(7)
    pick = sorted(rng.sample(loci, n))
    print(f"{len(loci):,} loci of the kind; {len(pick):,} drawn", flush=True)
    by = {}
    for x in pick:
        by.setdefault(x[0], []).append(x)
    pieces = [v[i:i + 50] for v in by.values() for i in range(0, len(v), 50)]
    t0 = time.time()
    done = 0
    with gzip.open(out + ".part", "wt") as w, Pool(workers) as pool:
        for batch in pool.imap_unordered(job, pieces):
            for ch, s, e, recs in batch:
                w.write(f"#region\t{ch}\t{s}\t{e}\tTR\n")
                w.writelines(f"{pos0}\t{ref}\t{alt}\t{ac}\t{af!r}\t{an}\n" for pos0, ref, alt, ac, af, an in recs)
            done += len(batch)
            if done % 2000 < len(batch):
                print(f"{done:,} of {len(pick):,} [{time.time() - t0:.0f}s]", flush=True)
    os.replace(out + ".part", out)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], *(int(x) for x in sys.argv[3:5]))
