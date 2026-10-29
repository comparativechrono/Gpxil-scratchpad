"""The pure GA- and TC-sequences of GRCh38 (more than 50 bases of only G/A, or only C/T, on the + strand), and the
"common" ones (more than 20 bases), per chromosome, from the Gpxil browser's gpxil files.

    python3 census.py REF_DIR

REF_DIR holds chr1.gpxil.tif ... chrY.gpxil.tif (the Gpxil browser's data/grch38/gpxil). Writes
results/ga-sequences.GRCh38.bed.gz (every pure GA- and TC-sequence: chrom, start, end, kind, length) and
results/census.json.
"""

import gzip
import sys
import time

import numpy as np

import ga
from gpxiltools import codes as C


def main(ref_dir: str) -> None:
    t0 = time.time()
    rows, lines = [], []
    for ch in ga.CHROMS:
        with ga.reference(ref_dir, ch) as g:
            base = g[ch].read("ref") & 0x0F
        acgt = np.isin(base, (C.BASE_A, C.BASE_C, C.BASE_G, C.BASE_T))
        pur = (base == C.BASE_A) | (base == C.BASE_G)
        pyr = (base == C.BASE_C) | (base == C.BASE_T)
        row = {"chrom": ch, "length": int(base.size), "acgt": int(acgt.sum())}
        for kind, m in (("GA", pur), ("TC", pyr)):
            s, e = ga.runs(m, ga.PURE)
            L = e - s
            homo = sum(np.unique(base[a:b]).size == 1 for a, b in zip(s.tolist(), e.tolist()))
            row[kind] = {"n": int(s.size), "bases": int(L.sum()), "max": int(L.max()) if L.size else 0, "homopolymers": int(homo)}
            sc, ec = ga.runs(m, ga.COMMON)
            row[kind + "_common"] = {"n": int(sc.size), "bases": int((ec - sc).sum())}
            lines += [(ch, a, b, kind) for a, b in zip(s.tolist(), e.tolist())]
        rows.append(row)
        print(f"{ch}: GA {row['GA']['n']:,} TC {row['TC']['n']:,} [{time.time() - t0:.0f}s]", flush=True)
    order = {c: i for i, c in enumerate(ga.CHROMS)}
    lines.sort(key=lambda x: (order[x[0]], x[1]))
    with gzip.open(ga.RUNS_BED, "wt") as fh:
        fh.write("#chrom\tstart\tend\tkind\tlength\n")
        for ch, a, b, kind in lines:
            fh.write(f"{ch}\t{a}\t{b}\t{kind}\t{b - a}\n")
    lengths = np.array([b - a for _, a, b, _ in lines])
    total = {k: sum(r[k]["n"] for r in rows) for k in ("GA", "TC", "GA_common", "TC_common")}
    summary = {
        "total": total,
        "bases": sum(r["GA"]["bases"] + r["TC"]["bases"] for r in rows),
        "length_quantiles": {q: float(np.percentile(lengths, q)) for q in (50, 90, 99)},
        "longest": int(lengths.max()),
        "homopolymers": sum(r["GA"]["homopolymers"] + r["TC"]["homopolymers"] for r in rows),
        "acgt_bases": sum(r["acgt"] for r in rows),
        "per_mb": {r["chrom"]: round((r["GA"]["n"] + r["TC"]["n"]) / (r["acgt"] / 1e6), 2) for r in rows},
    }
    ga.write_json("census.json", {"summary": summary, "chromosomes": rows})
    print("total:", total, "bases:", summary["bases"], "longest:", summary["longest"])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
