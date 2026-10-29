"""Where do pure GA-sequences come from? Each is put beside GRCh38's RepeatMasker annotation (UCSC's rmsk table):

- an **Alu tail**: an Alu element on the GA-sequence's purine strand whose 3' end lies within 30 bases of the
  GA-sequence's start (Alu elements end in a poly(A) tail, and a GA-sequence that starts with poly(A) begins there);
- else the repeat class it overlaps most (a simple repeat, low complexity, an L1, an SVA ...), or none.

The share of Alu tails among GA-sequences that start with a poly(A) of each length (stack.py's poly_a, as in his
image), and Alu elements a megabase of each autosome beside its GA-sequences a megabase (census.py's). Then, for Alu elements of each subfamily and of each divergence from their consensus (milliDiv, per thousand:
an age), the share that end in a pure GA-sequence, and what such tails are made of (A, and G, and their longest run of
A), and images: the tails sorted by the Alu's divergence, youngest at the top, and 150 of each of three ages.

    python3 -I alu.py REF_DIR RMSK          -> results/alu.json, results/ga-sequences.origin.tsv.gz,
                                               figures/alu_tails_by_age.png, figures/alu_tails_three_ages.png

RMSK is UCSC's hg38 rmsk.txt.gz (http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz). Run it with
python -I from a directory other than the download's.
"""

import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402
from PIL import Image, ImageDraw  # noqa: E402

import ga  # noqa: E402
from gpxiltools import codes as C  # noqa: E402
from stack import poly_a  # noqa: E402

NEAR = 30
SHOW = 120  # (bases of each tail shown in the image)


def read_rmsk(path: str):
    """RepeatMasker records of the 22 autosomes and X, Y: by chromosome, arrays of start, end, strand (+1/-1),
    milliDiv, and the indices of repName, repClass and repFamily in the lists returned."""
    from array import array

    names, classes, fams = {}, {}, {}
    keep = set(ga.CHROMS)
    cols = ("start", "end", "strand", "div", "name", "cls", "fam")
    tmp: dict[str, dict[str, array]] = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 17 or f[5] not in keep:
                continue
            t = tmp.setdefault(f[5], {c: array("q") for c in cols})
            t["start"].append(int(f[6]))
            t["end"].append(int(f[7]))
            t["strand"].append(1 if f[9] == "+" else -1)
            t["div"].append(int(f[2]))
            t["name"].append(names.setdefault(f[10], len(names)))
            t["cls"].append(classes.setdefault(f[11], len(classes)))
            t["fam"].append(fams.setdefault(f[12], len(fams)))
    out = {}
    for ch, t in tmp.items():
        a = {c: np.frombuffer(t[c], dtype=np.int64) for c in cols}
        order = np.argsort(a["start"], kind="stable")
        out[ch] = {c: a[c][order] for c in cols}
    inv = lambda d: [k for k, _ in sorted(d.items(), key=lambda kv: kv[1])]  # noqa: E731
    return out, inv(names), inv(classes), inv(fams)


def alu_class(name: str) -> str:
    for p in ("AluY", "AluS", "AluJ"):
        if name.startswith(p):
            return p
    return "other Alu"


def main(ref_dir: str, rmsk_path: str) -> None:
    rm, names, classes, fams = read_rmsk(rmsk_path)
    alu_fam = fams.index("Alu")
    runs = ga.read_runs()
    rows = []
    tails = []  # (milliDiv, subfamily, tail bases on the purine strand)
    for ch in ga.AUTOSOMES:
        r = rm[ch]
        is_alu = r["fam"] == alu_fam
        # Alu 3' ends on each strand: + strand at end, - strand at start
        ends_p = np.sort(r["end"][is_alu & (r["strand"] == 1)])
        idx_p = np.flatnonzero(is_alu & (r["strand"] == 1))[np.argsort(r["end"][is_alu & (r["strand"] == 1)], kind="stable")]
        starts_m = np.sort(r["start"][is_alu & (r["strand"] == -1)])
        idx_m = np.flatnonzero(is_alu & (r["strand"] == -1))[np.argsort(r["start"][is_alu & (r["strand"] == -1)], kind="stable")]
        maxlen = int((r["end"] - r["start"]).max())
        with ga.reference(ref_dir, ch) as g:
            base = g[ch].read("ref") & 0x0F
        for s, e, kind in runs.get(ch, []):
            five = s if kind == "GA" else e  # (the GA-sequence's 5' end on its purine strand)
            if kind == "GA":
                lo, hi = np.searchsorted(ends_p, [five - NEAR, five + NEAR + 1])
                cand = idx_p[lo:hi]
                dist = r["end"][cand] - five if cand.size else np.array([])
            else:
                lo, hi = np.searchsorted(starts_m, [five - NEAR, five + NEAR + 1])
                cand = idx_m[lo:hi]
                dist = five - r["start"][cand] if cand.size else np.array([])
            seq = ga.purine_strand(base[s:e], kind)
            n_a = int((seq == C.BASE_A).sum())
            # the longest run of A in the tail
            d = np.diff(np.concatenate(([0], (seq == C.BASE_A).view(np.int8), [0])))
            longest_a = int((np.flatnonzero(d == -1) - np.flatnonzero(d == 1)).max()) if n_a else 0
            row = {"chrom": ch, "start": s, "end": e, "kind": kind, "length": e - s, "a_fraction": round(n_a / (e - s), 3),
                   "longest_a": longest_a, "poly_a": poly_a(seq, 0)}
            if cand.size:
                k = cand[np.argmin(np.abs(dist))]
                row.update({"origin": "Alu tail", "alu": names[r["name"][k]], "alu_class": alu_class(names[r["name"][k]]),
                            "alu_div": int(r["div"][k]), "distance": int(dist[np.argmin(np.abs(dist))])})
                tails.append((int(r["div"][k]), row["alu_class"], seq[:SHOW]))
            else:
                # the repeat class overlapping the GA-sequence most
                i0 = np.searchsorted(r["start"], s - maxlen)
                i1 = np.searchsorted(r["start"], e)
                best, cover = None, 0
                for j in range(i0, i1):
                    ov = min(e, int(r["end"][j])) - max(s, int(r["start"][j]))
                    if ov > cover:
                        best, cover = j, ov
                if best is not None and cover >= 0.5 * (e - s):
                    row["origin"] = classes[r["cls"][best]] + (f" ({fams[r['fam'][best]]})" if classes[r["cls"][best]] in ("LINE", "SINE", "Retroposon") else "")
                    row["repeat"] = names[r["name"][best]]
                else:
                    row["origin"] = "none"
            rows.append(row)
        print(ch, len(rows), flush=True)
    # Alu elements of the autosomes, by subfamily class and divergence: how many end in a pure GA-sequence
    alus_by = {}
    for ch in ga.AUTOSOMES:
        r = rm[ch]
        m = r["fam"] == alu_fam
        for nm, dv in zip(r["name"][m].tolist(), r["div"][m].tolist()):
            key = (alu_class(names[nm]), min(dv // 20, 15))
            alus_by[key] = alus_by.get(key, 0) + 1
    tails_by = {}
    for dv, cl, _ in tails:
        key = (cl, min(dv // 20, 15))
        tails_by[key] = tails_by.get(key, 0) + 1
    by_class = {}
    for cl in ("AluY", "AluS", "AluJ", "other Alu"):
        n = sum(v for (c, _), v in alus_by.items() if c == cl)
        t = sum(v for (c, _), v in tails_by.items() if c == cl)
        sel = [x for x in rows if x.get("alu_class") == cl]
        by_class[cl] = {"alus": n, "with_ga_tail": t, "per_1000": round(1000 * t / max(n, 1), 3),
                        "tail_a_fraction": round(float(np.mean([x["a_fraction"] for x in sel])), 3) if sel else None,
                        "tail_length": round(float(np.median([x["length"] for x in sel])), 1) if sel else None,
                        "tail_longest_a": round(float(np.median([x["longest_a"] for x in sel])), 1) if sel else None}
    by_div = []
    for b in range(16):
        n = sum(v for (_, bb), v in alus_by.items() if bb == b)
        t = sum(v for (_, bb), v in tails_by.items() if bb == b)
        sel = [x for x in rows if "alu_div" in x and min(x["alu_div"] // 20, 15) == b]
        by_div.append({"milliDiv": f"{20 * b}-{20 * b + 19}" if b < 15 else ">=300", "alus": n, "with_ga_tail": t,
                       "per_1000": round(1000 * t / max(n, 1), 3),
                       "tail_a_fraction": round(float(np.mean([x["a_fraction"] for x in sel])), 3) if sel else None,
                       "tail_longest_a": round(float(np.median([x["longest_a"] for x in sel])), 1) if sel else None})
    origins = {}
    for x in rows:
        origins[x["origin"]] = origins.get(x["origin"], 0) + 1
    by_poly_a = {}
    for label, lo, hi in (("<5", 0, 5), ("5-14", 5, 15), (">=15", 15, 1 << 30)):
        sel = [x for x in rows if lo <= x["poly_a"] < hi]
        t = sum(x["origin"] == "Alu tail" for x in sel)
        by_poly_a[label] = {"ga_sequences": len(sel), "alu_tails": t, "share": round(t / max(len(sel), 1), 4)}
    # Alu elements a megabase (of A/C/G/T bases) of each autosome, beside its pure GA/TC-sequences a megabase
    from scipy.stats import spearmanr

    with open(os.path.join(ga.RESULTS, "census.json")) as fh:
        census = json.load(fh)
    acgt = {c["chrom"]: c["acgt"] for c in census["chromosomes"]}
    ga_per_mb = census["summary"]["per_mb"]
    alu_per_mb = {ch: round(1e6 * int((rm[ch]["fam"] == alu_fam).sum()) / acgt[ch], 1) for ch in ga.AUTOSOMES}
    rho = spearmanr([alu_per_mb[ch] for ch in ga.AUTOSOMES], [ga_per_mb[ch] for ch in ga.AUTOSOMES])
    density = {"alu_per_mb": alu_per_mb, "ga_per_mb": {ch: ga_per_mb[ch] for ch in ga.AUTOSOMES},
               "spearman_rho": round(float(rho.statistic), 3), "p": float(rho.pvalue)}
    out = {"ga_sequences": len(rows), "origins": dict(sorted(origins.items(), key=lambda kv: -kv[1])),
           "alu_tails_by_poly_a": by_poly_a, "alu_tails_by_class": by_class, "alu_tails_by_divergence": by_div,
           "density_by_autosome": density}
    ga.write_json("alu.json", out)
    with gzip.open(os.path.join(ga.RESULTS, "ga-sequences.origin.tsv.gz"), "wt") as fh:
        cols = ["chrom", "start", "end", "kind", "length", "a_fraction", "longest_a", "poly_a", "origin", "alu", "alu_class", "alu_div",
                "distance", "repeat"]
        fh.write("\t".join(cols) + "\n")
        for x in rows:
            fh.write("\t".join(str(x.get(c, "")) for c in cols) + "\n")
    print(out["origins"])
    print("Alu tails by poly(A):", by_poly_a)
    print("Alu elements and GA-sequences a megabase, by autosome:", density)
    for cl, v in by_class.items():
        print(cl, v)
    for v in by_div:
        print(v)
    image(tails)


def image(tails) -> None:
    """The Alu tails, one per row, sorted by the Alu's divergence (youngest at the top), first SHOW bases, gpxi greys."""
    tails = sorted(tails, key=lambda t: t[0])
    lut = np.full(16, 230, dtype=np.uint8)
    for k, v in {C.BASE_A: 0, C.BASE_C: 85, C.BASE_G: 255, C.BASE_T: 170}.items():
        lut[k] = v
    img = np.full((len(tails), SHOW), 230, dtype=np.uint8)
    for y, (_, _, seq) in enumerate(tails):
        img[y, : seq.size] = lut[seq]
    h = 900
    view = Image.fromarray(img, mode="L").convert("RGB").resize((SHOW * 4, h), Image.NEAREST)
    out = Image.new("RGB", (SHOW * 4 + 200, h + 60), (255, 255, 255))
    d = ImageDraw.Draw(out)
    out.paste(view, (10, 40))
    d.text((10, 6), f"The {len(tails):,} pure GA-sequences that are Alu tails, first {SHOW} bases (A black, G white),", fill=(0, 0, 0))
    d.text((10, 20), "one per row, sorted by the Alu's divergence from its consensus", fill=(0, 0, 0))
    divs = [t[0] for t in tails]
    for q in (0.0, 0.25, 0.5, 0.75, 0.999):
        y = 40 + int(h * q)
        d.text((SHOW * 4 + 20, y), f"{divs[min(int(q * len(divs)), len(divs) - 1)] / 10:.1f} % divergent", fill=(0, 0, 0))
    out.save(os.path.join(ga.FIGURES, "alu_tails_by_age.png"))
    # 150 tails drawn at random from each of three ages, side by side, three pixels a base
    import random

    rng = random.Random(7)
    groups = [("young: < 6 % divergent", lambda dv: dv < 60), ("middle: 10-14 %", lambda dv: 100 <= dv < 140),
              ("old: >= 18 %", lambda dv: dv >= 180)]
    n, sc, gap = 150, 3, 24
    panel_w = SHOW * sc
    fig = Image.new("RGB", (3 * panel_w + 2 * gap + 20, n * sc + 70), (255, 255, 255))
    d = ImageDraw.Draw(fig)
    d.text((10, 6), "Alu tails (pure GA-sequences at an Alu's 3' end), 150 at random from each age; first 120 bases, A black, G white", fill=(0, 0, 0))
    for k, (label, f) in enumerate(groups):
        pool = [t for t in tails if f(t[0])]
        pick = rng.sample(pool, min(n, len(pool)))
        arr = np.full((len(pick), SHOW), 230, dtype=np.uint8)
        for y, (_, _, seq) in enumerate(pick):
            arr[y, : seq.size] = lut[seq]
        im = Image.fromarray(arr, mode="L").convert("RGB").resize((panel_w, len(pick) * sc), Image.NEAREST)
        x0 = 10 + k * (panel_w + gap)
        fig.paste(im, (x0, 44))
        a_frac = np.mean([float((s == C.BASE_A).mean()) for _, _, s in pool])
        d.text((x0, 26), f"{label} ({len(pool):,}; {100 * a_frac:.0f} % A)", fill=(0, 0, 0))
    fig.save(os.path.join(ga.FIGURES, "alu_tails_three_ages.png"))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
