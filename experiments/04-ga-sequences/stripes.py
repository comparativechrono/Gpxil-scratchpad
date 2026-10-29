"""How many GA-sequences have the Alu pattern upstream? An image measure, without alignment: in the stack of
stack.py, the most frequent base in each column, among the rows whose GA-sequence starts with a poly(A) of 15 or more,
makes a "consensus row"; each row scores the fraction of its 300 upstream bases equal to it. A row passes when it
scores above the 99.9th percentile of rows shuffled within themselves. TC-sequences, the GA-sequences of the other
strand, are scored the same way on their own strand.

    python3 stripes.py REF_DIR [CHROM ...]       -> results/stripes.json (default: chr1, chr7, chrX)
"""

import sys

import numpy as np

import ga
from gpxiltools import codes as C
from stack import UP, rows_of

ACGT = np.array([C.BASE_A, C.BASE_C, C.BASE_G, C.BASE_T])


def main(ref_dir: str, chroms) -> None:
    runs = ga.read_runs()
    report = {}
    for ch in chroms:
        with ga.reference(ref_dir, ch) as g:
            base = g[ch].read("ref") & 0x0F
        for kind in ("GA", "TC"):
            rows, pa = rows_of(base, runs[ch], kind)
            up = rows[:, UP - 300 : UP]
            top = up[pa >= 15]
            cons = ACGT[[np.bincount(col, minlength=16)[ACGT].argmax() for col in top.T]]
            score = (up == cons).mean(axis=1)
            rng = np.random.default_rng(1)
            null = np.array([(rng.permutation(x) == cons).mean() for x in up])
            thr = float(np.percentile(null, 99.9))
            out = {"rows": int(len(rows)), "threshold": round(thr, 3)}
            for label, m in (("polyA>=15", pa >= 15), ("polyA 5-14", (pa >= 5) & (pa < 15)), ("polyA<5", pa < 5)):
                out[label] = {"n": int(m.sum()), "alu_like": int((score[m] > thr).sum()), "mean_score": round(float(score[m].mean()), 3)}
            report[f"{ch} {kind}"] = out
            print(ch, kind, out)
    ga.write_json("stripes.json", report)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2:] or ["chr1", "chr7", "chrX"])
