"""The Alu tails stacked as Albrecht-Buehler drew GA-sequences: one a row, lined up at the end next to the Alu, read on
the Alu's strand, youngest Alu at the top (gpxiltools' stack). And each column's bases, by the Alu's age.

    python3 -I stacks.py GENOME      -> figures/stack_alu_tails.png, results/stack_rows.npy, results/stack_profile.json

GENOME: GRCh38 as one gpxil file (the Gpxil browser's). Reads results/tail_loci.tsv.gz (tails.py).
"""

from __future__ import annotations

import dataclasses
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from gpxiltools import codes as GC  # noqa: E402
from gpxiltools.container import open_gpxil  # noqa: E402
from gpxiltools.features import Feature  # noqa: E402
from gpxiltools.stack import binned, labelled, legend_of, rendering, stack  # noqa: E402
from tails import AGES, age_label  # noqa: E402

RESULTS = os.path.join(HERE, "results")
FIGURES = os.path.join(HERE, "figures")
UP, DOWN = 40, 100  #: bases of the Alu before the array's start, and of the array and beyond from it
BINS = 600  #: rows of the binned picture


def main(genome: str) -> None:
    import pandas as pd
    t = pd.read_csv(os.path.join(RESULTS, "tail_loci.tsv.gz"), sep="\t", usecols=["chrom", "start", "end", "trid", "fam", "te_strand", "div"])
    t = t[t.fam == "Alu"]
    feats = [Feature(c, int(s), int(e), n, float(d), st) for c, s, e, n, d, st in zip(t.chrom, t.start, t.end, t.trid, t["div"], t.te_strand)]
    del t
    with open_gpxil(genome) as g:
        st = stack(g, feats, anchor="start", up=UP, down=DOWN, sort="score", reverse=True)
    print(f"{len(st.rows):,} rows; unplaced {st.unplaced}", flush=True)
    div = np.array([r.key for r in st.rows], dtype=float)
    length = np.array([r.feature.length for r in st.rows])
    nib = st.ref & 0x0F
    prof = {}
    for a in [age_label(x + 0.01) for x in AGES[:-1]]:
        m = np.array([age_label(x) == a for x in div])
        if not m.any():
            continue
        x = nib[m]
        have = (x != 0).sum(0)
        within = (np.arange(st.width)[None, :] - UP < length[m][:, None]) & (np.arange(st.width)[None, :] >= UP)
        prof[a] = {"rows": int(m.sum()),
                   **{b: [round(float(v), 4) for v in ((x == code).sum(0) / np.maximum(have, 1))]
                      for b, code in zip("ACGT", (GC.BASE_A, GC.BASE_C, GC.BASE_G, GC.BASE_T))},
                   "in_array": [round(float(v), 4) for v in within.mean(0)],
                   "a_in_array": [round(float(v), 4) for v in (((x == GC.BASE_A) & within).sum(0) / np.maximum(within.sum(0), 1))]}
    total = len(st.rows)
    every = 3  # (the picture draws every third tail, so that drawing a million rows stays small)
    st = dataclasses.replace(st, rows=st.rows[::every], ref=st.ref[::every])
    div = div[::every]
    r = rendering(st, mode="reference")
    k = max(1, len(st.rows) // BINS)
    px = binned(r, k)
    os.makedirs(FIGURES, exist_ok=True)
    np.save(os.path.join(RESULTS, "stack_rows.npy"), px)
    edges = [float(div[min(i, len(div) - 1)]) for i in range(0, len(div), k)]
    img = labelled(px, st, scale=2, title=f"{len(st.rows):,} Alu tails (every third), lined up at the end next to the Alu, youngest Alu first",
                   legend=legend_of(r), rows_note=f"each pixel row the mean of {k:,} tails")
    img.save(os.path.join(FIGURES, "stack_alu_tails.png"))
    with open(os.path.join(RESULTS, "stack_profile.json"), "w") as fh:
        json.dump({"up": UP, "down": DOWN, "tails": total, "rows_drawn": len(st.rows), "every": every, "rows_a_bin": k, "bin_divergence_start": edges,
                   "by_age": prof}, fh)
    print({a: (v["rows"], v["a_in_array"][UP:UP + 5]) for a, v in prof.items()})


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
