"""A check of gnomad.py's reads from the web: a random sample of its regions read again with another reader (pysam's
VariantFile, which parses each record, where gnomad.py split the lines itself), and their PASS records counted.

    CURL_CA_BUNDLE=... python3 -I check_gnomad.py [N]          -> results/gnomad_check.json

Every count should agree with results/gnomad_counts.tsv.gz; a read cut short by the network without an error would
show here as fewer records.
"""

import gzip
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ga  # noqa: E402
from gnomad import URL  # noqa: E402

SEED = 9


def count(vf, ch: str, s: int, e: int) -> int:
    return sum(1 for r in vf.fetch(ch, s, e) if s <= r.start < e and list(r.filter.keys()) == ["PASS"])


def main(n: int = 300) -> None:
    import pysam

    rows = []
    with gzip.open(os.path.join(ga.RESULTS, "gnomad_counts.tsv.gz"), "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                ch, s, e, kind, k = line.rstrip("\n").split("\t")
                rows.append((ch, int(s), int(e), kind, int(k)))
    sample = sorted(random.Random(SEED).sample(rows, n), key=lambda r: (ga.CHROMS.index(r[0]), r[1]))
    t0 = time.time()
    agree, differ, records, vf, open_ch = 0, [], 0, None, None
    for ch, s, e, kind, k in sample:
        for attempt in range(6):
            try:
                if open_ch != ch:
                    vf, open_ch = pysam.VariantFile(URL.format(ch), index_filename=URL.format(ch) + ".tbi"), ch
                got = count(vf, ch, s, e)
                break
            except (OSError, ValueError):
                open_ch = None
                if attempt == 5:
                    raise
                time.sleep(2 + 3 * attempt)
        records += got
        if got == k:
            agree += 1
        else:
            differ.append({"chrom": ch, "start": s, "end": e, "kind": kind, "gnomad_py": k, "again": got})
    out = {"regions": len(rows), "sample": n, "seed": SEED, "agree": agree, "differ": differ, "pass_records": records,
           "background_windows": sum(r[3] == "background" for r in sample)}
    ga.write_json("gnomad_check.json", out)
    print(f"{agree} of {n} regions agree ({records:,} PASS records read again) [{time.time() - t0:.0f}s]")
    for d in differ:
        print("differs:", d)


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 300)
