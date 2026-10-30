"""Which forensic STR markers are pure GA-sequences, and where did they come from?

The autosomal STRs of forensic kits are placed on GRCh38 from their sequences: lusSTR's test alleles (each marker's
longest, with its flanks) and, for SE33, which lusSTR lacks, GenBank V00481 (bases 14 to 356). A 24-base window of a
flank, at 0 to 40 bases from either end of the sequence, found exactly once in the autosomes (either strand) places a
marker; a marker whose place disagrees with its name's chromosome (D9S1122's test allele, for one) is left unplaced.
In the placed span (20 bases more each side), the marker's repeat is the longest array of whole units of its motif
(lusSTR's LUS motif; AAAG for SE33), any rotation, on one strand, gaps of up to a unit allowed; the LUS is its longest
run of whole units.

Each array is then set against experiment 04's census (the pure GA/TC-sequences of more than 50 bases) and their
origins, and against GRCh38's RepeatMasker annotation: an Alu on the array's purine-rich strand whose 3' end lies within
30 bases of the array's start (experiment 04's rule for an Alu tail, applied to the array). Liu et al. 2021's mutation
rates (19 loci; 28,313 paternity cases, Chinese Han) are set beside them.

    python3 -I markers.py REF_DIR LUSSTR_DIR V00481.fasta RMSK       -> results/markers.json, results/markers.tsv
"""

from __future__ import annotations

import csv
import glob
import gzip
import json
import os
import sys
from collections import defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
E04 = os.path.join(os.path.dirname(HERE), "04-ga-sequences")
sys.path.insert(0, E04)

import ga  # noqa: E402  (experiment 04)
from gpxiltools import codes as GC  # noqa: E402

RESULTS = os.path.join(HERE, "results")
K = 24  #: a placing window
OFFSETS = (0, 8, 16, 24, 32, 40)
NEAR = 30  #: an Alu's 3' end this close to the array's start makes it an Alu tail
FAR = 100  #: the nearest Alu 3' end on the purine strand is reported this far
CHUNK = 1 << 24
COMP = str.maketrans("ACGT", "TGCA")
IDX = np.full(256, -1, dtype=np.int8)
for _i, _b in enumerate((GC.BASE_A, GC.BASE_C, GC.BASE_G, GC.BASE_T)):
    IDX[_b] = IDX[_b | GC.SOFTMASK] = _i
LET = np.full(256, ord("N"), dtype=np.uint8)
for _b, _ch in ((GC.BASE_A, "A"), (GC.BASE_C, "C"), (GC.BASE_G, "G"), (GC.BASE_T, "T")):
    LET[_b] = LET[_b | GC.SOFTMASK] = ord(_ch)
#: the chromosome of each marker not named for it
CHROM_OF = {"CSF1PO": "chr5", "FGA": "chr4", "TH01": "chr11", "TPOX": "chr2", "VWA": "chr12", "PENTAD": "chr21", "PENTAE": "chr15",
            "SE33": "chr6"}
#: Liu et al. 2021 (Hereditas (Beijing) 43:949-961, doi:10.16288/j.yczz.21-197), Table 1 and Table 1 Supplementary: the
#: trio rate (mean of paternal and maternal), and trio mutations and trios
LIU2021 = {"CSF1PO": (2.6945e-3, 87, 22113), "D12S391": (4.315e-3, 158, 20562), "D13S317": (1.4355e-3, 51, 22113),
           "D16S539": (1.6e-3, 56, 22113), "D18S51": (3.3265e-3, 128, 22113), "D19S433": (1.728e-3, 60, 22113),
           "D1S1656": (1.9225e-3, 48, 13712), "D21S11": (2.2847e-3, 92, 22113), "D2S1338": (1.6075e-3, 64, 22113),
           "D3S1358": (1.1697e-3, 48, 22113), "D5S818": (1.8365e-3, 59, 22113), "D7S820": (1.8795e-3, 67, 22113),
           "D8S1179": (2.624e-3, 94, 22113), "FGA": (4.5905e-3, 172, 22113), "PENTAD": (1.1595e-3, 40, 20562),
           "PENTAE": (3.755e-3, 142, 20562), "TH01": (1.6306e-4, 7, 22113), "TPOX": (7.7e-5, 3, 22113), "VWA": (3.547e-3, 124, 22113)}


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def name_of(locus: str) -> str:
    return locus.strip().upper().replace(" ", "")


def sequences(lusstr: str, se33: str) -> tuple[dict, dict]:
    """The longest test allele of each autosomal marker (lusSTR), and SE33's from V00481; the motifs."""
    best: dict = {}
    for path in sorted(glob.glob(os.path.join(lusstr, "lusSTR", "tests", "data", "*.csv")) +
                       glob.glob(os.path.join(lusstr, "lusSTR", "tests", "data", "*.txt"))):
        try:
            with open(path, newline="") as fh:
                head = fh.read(4096)
                fh.seek(0)
                for r in csv.DictReader(fh, delimiter="\t" if head.count("\t") > head.count(",") else ","):
                    loc = name_of(r.get("Locus") or "")
                    seq = (r.get("Full_Sequence") or r.get("Sequence") or r.get("ForenSeq") or "").strip().upper()
                    if loc and seq and set(seq) <= set("ACGT") and len(seq) > len(best.get(loc, "")):
                        best[loc] = seq
        except (UnicodeDecodeError, csv.Error):
            pass
    with open(os.path.join(lusstr, "lusSTR", "data", "str_markers.json")) as fh:
        defs = json.load(fh)
    motifs = {name_of(m): v["LUS"] for m, v in defs.items() if v.get("LUS")}
    keep = {m: s for m, s in best.items() if m in motifs and not m.startswith(("DX", "DY", "Y-", "AMEL", "HPRT"))}
    with open(se33) as fh:
        v = "".join(x.strip() for x in fh if not x.startswith(">")).upper()
    keep["SE33"] = v[13:356]
    motifs["SE33"] = "AAAG"
    return keep, motifs


def code_of(s: str) -> int:
    v = 0
    for ch in s:
        v = v * 4 + "ACGT".index(ch)
    return v


def place(ref_dir: str, seqs: dict) -> dict:
    """Each marker's place: the first window (left side, then right; nearest the end first) found once."""
    wins = {}
    for m, s in seqs.items():
        for side in ("L", "R"):
            for o in OFFSETS:
                if o + K <= len(s) // 2:
                    w = s[o:o + K] if side == "L" else s[len(s) - o - K:len(s) - o]
                    wins[(m, side, o)] = w
    targets = defaultdict(list)
    for key, w in wins.items():
        targets[code_of(w)].append((key, "+"))
        targets[code_of(rc(w))].append((key, "-"))
    tarr = np.array(sorted(targets), dtype=np.int64)
    hits = defaultdict(list)
    for ch in ga.AUTOSOMES:
        with ga.reference(ref_dir, ch) as g:
            c = g[ch]
            for a in range(0, c.length - K + 1, CHUNK):
                b = min(c.length, a + CHUNK + K - 1)
                x = IDX[c.read("ref", a, b)]
                n = len(x) - K + 1
                h = np.zeros(n, dtype=np.int64)
                bad = np.zeros(n, dtype=bool)
                for k in range(K):
                    seg = x[k:k + n]
                    h = (h << 2) | np.maximum(seg, 0).astype(np.int64)
                    bad |= seg < 0
                for i in np.nonzero(np.isin(h, tarr) & ~bad)[0]:
                    for key, strand in targets[int(h[i])]:
                        hits[key].append((ch, a + int(i), strand))
        print(ch, flush=True)
    out = {}
    for m, s in seqs.items():
        out[m] = None
        for side in ("L", "R"):
            for o in OFFSETS:
                h = hits.get((m, side, o))
                if h is not None and len(h) == 1:
                    ch, i, strand = h[0]
                    L = len(s)
                    if side == "L":
                        s0 = i - o if strand == "+" else i + K + o - L
                    else:
                        s0 = i + K + o - L if strand == "+" else i - o
                    out[m] = {"chrom": ch, "start": s0, "end": s0 + L, "strand": strand, "anchor": f"{side}{o}"}
                    break
            if out[m]:
                break
        want = CHROM_OF.get(m) or ("chr" + m[1:m.index("S")] if m.startswith("D") and "S" in m else None)
        if out[m] and want and out[m]["chrom"] != want:
            out[m] = {"unplaced": f"a window was found once on {out[m]['chrom']}, not {want}"}
    return out


def array_of(seq: str, motif: str) -> tuple[int, int, int, str] | None:
    """The longest array of whole units of the motif's rotations on one strand, gaps of up to a unit allowed:
    (start, end, its longest run of whole units, the strand's motif)."""
    p = len(motif)
    best = None
    for m in (motif, rc(motif)):
        rot = {m[i:] + m[:i] for i in range(p)}
        hits = [i for i in range(len(seq) - p + 1) if seq[i:i + p] in rot]
        if not hits:
            continue
        groups, s, e = [], hits[0], hits[0]
        for h in hits[1:]:
            if h <= e + 2 * p:
                e = h
            else:
                groups.append((s, e + p))
                s = e = h
        groups.append((s, e + p))
        for s, e in groups:
            lus, i = 0, s
            while i <= e - p:
                if seq[i:i + p] in rot:
                    j = i
                    r = seq[i:i + p]
                    while seq[j:j + p] == r:
                        j += p
                    lus = max(lus, (j - i) // p)
                    i = j if j > i else i + 1
                else:
                    i += 1
            if best is None or e - s > best[1] - best[0]:
                best = (s, e, lus, m)
    return best


def read_alus(rmsk: str, spans: dict) -> dict:
    """RepeatMasker's Alu records near each array: by marker, (start, end, strand, name)."""
    out = defaultdict(list)
    with gzip.open(rmsk, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[12] != "Alu":
                continue
            c, s, e = f[5], int(f[6]), int(f[7])
            for m, (ch, a, b) in spans.items():
                if c == ch and s < b + 400 and e > a - 400:
                    out[m].append((s, e, f[9], f[10]))
    return out


def main(ref_dir: str, lusstr: str, se33: str, rmsk: str) -> None:
    from scipy.stats import mannwhitneyu, spearmanr

    seqs, motifs = sequences(lusstr, se33)
    placed = place(ref_dir, seqs)
    census = []
    with gzip.open(os.path.join(E04, "results", "ga-sequences.origin.tsv.gz"), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            census.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    rows = {}
    for m in sorted(seqs):
        p = placed[m]
        if not p or "unplaced" in p:
            rows[m] = {"marker": m, "placed": False, "note": (p or {}).get("unplaced", "no window found once")}
            continue
        ch, a, b = p["chrom"], p["start"] - 20, p["end"] + 20
        with ga.reference(ref_dir, ch) as g:
            span = LET[g[ch].read("ref", a, b)].tobytes().decode()
        arr = array_of(span, motifs[m])
        s0, e0 = a + arr[0], a + arr[1]
        purine = sum(ch_ in "AG" for ch_ in span[arr[0]:arr[1]]) >= (e0 - s0) / 2
        cen = [r for r in census if r["chrom"] == ch and int(r["start"]) < e0 and int(r["end"]) > s0]
        rows[m] = {"marker": m, "placed": True, "chrom": ch, "start": s0, "end": e0, "anchor": p["anchor"], "motif": motifs[m],
                   "array_length": e0 - s0, "lus_units": arr[2], "purine_strand": "+" if purine else "-",
                   "census": [{k: r[k] for k in ("start", "end", "kind", "length", "origin", "alu", "alu_class", "distance")} for r in cen]}
    alus = read_alus(rmsk, {m: (r["chrom"], r["start"], r["end"]) for m, r in rows.items() if r["placed"]})
    for m, r in rows.items():
        if not r["placed"]:
            continue
        # the array's start on its purine strand, and an Alu on that strand ending (3') near it
        if r["purine_strand"] == "+":
            near = [(abs(r["start"] - e), n) for s, e, st, n in alus[m] if st == "+" and abs(r["start"] - e) <= FAR]
        else:
            near = [(abs(s - r["end"]), n) for s, e, st, n in alus[m] if st == "-" and abs(s - r["end"]) <= FAR]
        r["nearest_alu"] = min(near)[1] if near else None
        r["alu_distance"] = min(near)[0] if near else None
        r["alu_at_start"] = r["nearest_alu"] if near and min(near)[0] <= NEAR else None
        r["ga_sequence"] = bool(r["census"])
        r["liu2021"] = LIU2021.get(m)
    os.makedirs(RESULTS, exist_ok=True)
    liu = [r for r in rows.values() if r.get("liu2021")]
    ga_rates = [r["liu2021"][0] for r in liu if r["ga_sequence"]]
    other_rates = [r["liu2021"][0] for r in liu if not r["ga_sequence"]]
    summary = {
        "markers": len(rows), "placed": sum(r["placed"] for r in rows.values()),
        "ga_sequences": sorted(m for m, r in rows.items() if r.get("ga_sequence")),
        "alu_tails": sorted(m for m, r in rows.items() if r.get("alu_at_start")),
        "liu2021": {"loci": len(liu), "ga_sequence_markers": sorted(r["marker"] for r in liu if r["ga_sequence"]),
                    "mean_rate_ga": float(np.mean(ga_rates)) if ga_rates else None, "mean_rate_other": float(np.mean(other_rates)) if other_rates else None,
                    "mann_whitney_p": float(mannwhitneyu(ga_rates, other_rates).pvalue) if ga_rates and other_rates else None,
                    "spearman_lus": [round(float(x), 3) for x in spearmanr([r["liu2021"][0] for r in liu], [r["lus_units"] for r in liu])],
                    "spearman_array_length": [round(float(x), 3) for x in spearmanr([r["liu2021"][0] for r in liu], [r["array_length"] for r in liu])]},
        "markers_table": rows}
    with open(os.path.join(RESULTS, "markers.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    cols = ["marker", "placed", "chrom", "start", "end", "motif", "array_length", "lus_units", "purine_strand", "ga_sequence",
            "census_length", "census_origin", "nearest_alu", "alu_distance", "alu_at_start", "liu2021_rate"]
    with open(os.path.join(RESULTS, "markers.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for m, r in sorted(rows.items()):
            c = r.get("census") or []
            vals = {**r, "census_length": c[0]["length"] if c else "", "census_origin": c[0]["origin"] if c else "",
                    "liu2021_rate": r["liu2021"][0] if r.get("liu2021") else ""}
            fh.write("\t".join("" if vals.get(k) is None else str(vals.get(k, "")) for k in cols) + "\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "markers_table"}, indent=1))


if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.exit(__doc__)
    main(*sys.argv[1:5])
