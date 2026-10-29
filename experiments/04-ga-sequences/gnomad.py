"""The GA-sequence tests on the variants of 76,215 people: gnomAD v4.1 genomes, read from the web a region at a time.

For every pure GA/TC-sequence of the 22 autosomes, and for a background of random 100-base windows elsewhere, the
PASS variants that gnomAD's genomes hold there (the sites VCF; one alternative allele a record). Then, on each
GA-sequence's purine strand:

- polymorphic substitutions per base, against what the background's rate of each substitution in each trinucleotide
  context leads one to expect; transitions against transversions;
- whether the transitions copy a neighbouring 6-base unit of the GA-sequence (unit_copy.py), by allele count:
  singletons (AC 1), rare (AC 2 to AF 0.1 %), uncommon (to 1 %) and common (1 % and more); and the same at the
  sites called in 90 % of the genomes or more (AN), where calls are most reliable;
- the same counts in each of those classes, inside and in the background, the expected ones from the background's
  rates in the same class;
- the lengths of insertions and deletions.

    CURL_CA_BUNDLE=... python3 -I gnomad.py REF_DIR [WORKERS [CHROM ...]]      -> results/gnomad.json,
                                                                                 results/gnomad_counts.tsv.gz

(CHROMs: only these, with a background in proportion, written to results/gnomad_test.json; for a trial run.)

(CURL_CA_BUNDLE, when the web is reached through a proxy with its own certificate. GNOMAD_RECORDS=FILE.tsv.gz keeps
a copy of the records read, about 40 MB, and reads them from there when it exists, instead of from the web.)
"""

import gzip
import os
import random
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np  # noqa: E402

import ga  # noqa: E402
from gpxiltools import codes as C  # noqa: E402

URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{}.vcf.bgz"
K = 6
TI = {"A": "G", "G": "A"}
BACKGROUND = 10000  # (random 100-base windows)
SEED = 4
CODE = {"A": C.BASE_A, "C": C.BASE_C, "G": C.BASE_G, "T": C.BASE_T}
CLASSES = ("singleton", "rare", "uncommon", "common")
GENOMES = 76215
WELL_CALLED = int(0.9 * 2 * GENOMES)  # (a site called in 90 % of the genomes or more: AN of 137,187 or more)
RECORDS = os.environ.get("GNOMAD_RECORDS")
_tabix: dict = {}
_ref: dict = {}
REF_DIR = None


def af_class(ac: int, af: float) -> str:
    if ac == 1:
        return "singleton"
    if af < 0.001:
        return "rare"
    if af < 0.01:
        return "uncommon"
    return "common"


def fetch(ch: str, s: int, e: int):
    """gnomAD's PASS records in [s, e): (pos0, ref, alt, ac, af, an)."""
    import pysam

    url = URL.format(ch)
    for attempt in range(6):
        try:
            if ch not in _tabix:
                _tabix[ch] = pysam.TabixFile(url, index=url + ".tbi")
            rows = list(_tabix[ch].fetch(ch, s, e))
            break
        except (OSError, ValueError):
            _tabix.pop(ch, None)
            if attempt == 5:
                raise
            time.sleep(2 + 3 * attempt)
    out = []
    for line in rows:
        f = line.split("\t", 8)
        if f[6] != "PASS":
            continue
        pos0 = int(f[1]) - 1
        if not (s <= pos0 < e):
            continue
        info = f[7]
        ac = int(info[3 : info.index(";")]) if info.startswith("AC=") else 0
        i = info.find(";AF=")
        af = float(info[i + 4 : info.index(";", i + 4)]) if i >= 0 else 0.0
        j = info.find(";AN=")
        an = int(info[j + 4 : info.index(";", j + 4)]) if j >= 0 else 0
        out.append((pos0, f[3], f[4], ac, af, an))
    return out


def reference(ch: str, s: int, e: int) -> tuple[np.ndarray, int]:
    """The + strand's bases over [s, e), cut at the chromosome's end, and its length."""
    if ch not in _ref:
        for g in _ref.values():
            g.close()
        _ref.clear()
        _ref[ch] = ga.reference(REF_DIR, ch)
    c = _ref[ch][ch]
    return c.read("ref", s, min(e, c.length)) & 0x0F, c.length


def init(ref_dir: str) -> None:
    global REF_DIR
    REF_DIR = ref_dir


def job(tasks):
    """The regions of one chromosome, in order (GA-sequences and background windows): each one's bases (+ strand, one
    more each side) and gnomAD's records. The chromosome's file is opened once, and again after a read error."""
    out = []
    for ch, s, e, kind in sorted(tasks, key=lambda t: t[1]):
        a = max(0, s - 1)
        seq, n = reference(ch, a, e + 1)
        e = min(e, n)
        out.append(((ch, s, e, kind), seq, a, fetch(ch, s, e)))
    _tabix.clear()
    return out


def with_bases(region, recs):
    """A region read from a local copy of the records, with its bases, as job() gives it."""
    ch, s, e, kind = region
    a = max(0, s - 1)
    seq, n = reference(ch, a, e + 1)
    return (ch, s, min(e, n), kind), seq, a, recs


def cached(path: str):
    """The regions and records of a local copy written by an earlier run, in its order."""
    with gzip.open(path, "rt") as fh:
        region, recs = None, []
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[0] == "#region":
                if region:
                    yield with_bases(region, recs)
                region, recs = (f[1], int(f[2]), int(f[3]), f[4]), []
            else:
                recs.append((int(f[0]), f[1], f[2], int(f[3]), float(f[4]), int(f[5])))
        if region:
            yield with_bases(region, recs)


def contexts(seq: np.ndarray, a: int, s: int, e: int):
    """Context ids (32: purine middle base) of bases s..e-1 of a + strand array that starts at a, with flip flags."""
    i4 = ga.IDX4[seq]
    out = []
    for p in range(s, e):
        j = p - a
        if j - 1 < 0 or j + 1 >= len(seq):
            out.append((-1, False))
            continue
        prv, mid, nxt = int(i4[j - 1]), int(i4[j]), int(i4[j + 1])
        if min(prv, mid, nxt) < 0:
            out.append((-1, False))
            continue
        flip = mid in (1, 3)
        out.append(((3 - nxt) * 16 + (3 - mid) * 4 + (3 - prv), True) if flip else (prv * 16 + mid * 4 + nxt, False))
    return out


def main(ref_dir: str, workers: int = 8, chroms=None) -> None:
    from scipy.stats import binomtest

    runs = ga.read_runs()
    test = bool(chroms)
    chroms = chroms or ga.AUTOSOMES
    tasks = [(ch, s, e, kind) for ch in chroms for s, e, kind in runs.get(ch, [])]
    # a background of random windows away from GA-sequences (those with an N are left out below)
    rng = random.Random(SEED)
    sizes = {}
    for ch in chroms:
        with ga.reference(ref_dir, ch) as g:
            sizes[ch] = g[ch].length
    total = sum(sizes.values())
    near = {ch: ga.mask_of(sizes[ch], [s - 500 for s, _, _ in runs.get(ch, [])], [e + 500 for _, e, _ in runs.get(ch, [])])
            for ch in chroms}
    bg = []
    n_bg = BACKGROUND if not test else max(200, int(BACKGROUND * total / 2.88e9))
    while len(bg) < n_bg:
        x = rng.randrange(total)
        for ch in chroms:
            if x < sizes[ch]:
                break
            x -= sizes[ch]
        if x + 100 > sizes[ch] or near[ch][x : x + 100].any():
            continue
        bg.append((ch, x, x + 100, "background"))
    del near
    t0 = time.time()
    def tally():
        return {"bp": 0, "snv": 0, "ti": 0, "indel": 0, "indel_len": {}, "ctx_bp": [0] * 64, "ctx_snv": [[0] * 4 for _ in range(64)],
                "by_af": {c: {"snv": 0, "ti": 0, "indel": 0, "an": 0, "ctx_snv": [[0] * 4 for _ in range(64)]} for c in CLASSES}}

    res = {"inside": tally(), "background": tally()}
    copy = {"sites": 0, "sites_copy": 0, **{c: {"ti": 0, "ti_copy": 0, "tv": 0} for c in ("singleton", "rare", "uncommon", "common", "all")},
            "well_called": {c: {"ti": 0, "ti_copy": 0, "tv": 0} for c in ("singleton", "rare", "uncommon", "common", "all")}}
    by_ch: dict = {}
    for t in bg + tasks:
        by_ch.setdefault(t[0], []).append(t)
    # (the biggest chromosomes first, a piece of about 500 regions a task, so that the workers stay busy)
    pieces = []
    for ch, ts in by_ch.items():
        ts.sort(key=lambda t: t[1])
        pieces += [ts[i : i + 500] for i in range(0, len(ts), 500)]
    pieces.sort(key=len, reverse=True)
    done = 0

    def results():
        if RECORDS and os.path.exists(RECORDS):
            init(ref_dir)
            yield from cached(RECORDS)
            return
        with Pool(workers, initializer=init, initargs=(ref_dir,)) as pool:
            for batch in pool.imap_unordered(job, pieces):
                yield from batch

    write = gzip.open(RECORDS + ".part", "wt") if RECORDS and not os.path.exists(RECORDS) else None
    counts = []
    for (ch, s, e, kind), seq, a, recs in results():
        counts.append((ch, s, e, kind, len(recs)))
        if write:
            write.write(f"#region\t{ch}\t{s}\t{e}\t{kind}\n")
            write.writelines(f"{pos0}\t{ref}\t{alt}\t{ac}\t{af!r}\t{an}\n" for pos0, ref, alt, ac, af, an in recs)
        done += 1
        if done % 2000 == 0:
            print(f"{done:,} of {len(bg) + len(tasks):,} [{time.time() - t0:.0f}s]", flush=True)
        if (seq[s - a : e - a] == C.BASE_N).any() or (seq[s - a : e - a] == 0).any():
            continue
        d = res["background" if kind == "background" else "inside"]
        d["bp"] += e - s
        ctx = contexts(seq, a, s, e)
        for c, _ in ctx:
            if c >= 0:
                d["ctx_bp"][c] += 1
        snvs = {}
        for pos0, ref, alt, ac, af, an in recs:
            if alt == "*":
                continue
            if len(ref) == 1 and len(alt) == 1 and alt in CODE:
                snvs.setdefault(pos0, []).append((alt, ac, af, an))
            else:
                d["indel"] += 1
                d["by_af"][af_class(ac, af)]["indel"] += 1
                L = len(alt) - len(ref)
                d["indel_len"][str(L)] = d["indel_len"].get(str(L), 0) + 1
        for pos0, alts in snvs.items():
            c, flip = ctx[pos0 - s]
            for alt, ac, af, an in alts:
                cl = d["by_af"][af_class(ac, af)]
                d["snv"] += 1
                cl["snv"] += 1
                cl["an"] += an
                r = ga.LETTER.get(int(seq[pos0 - a]))
                if r is None:
                    continue
                if {r, alt} in ({"A", "G"}, {"C", "T"}):
                    d["ti"] += 1
                    cl["ti"] += 1
                if c >= 0:
                    ai = int(ga.IDX4[CODE[alt]])
                    d["ctx_snv"][c][3 - ai if flip else ai] += 1
                    cl["ctx_snv"][c][3 - ai if flip else ai] += 1
        if kind == "background":
            continue
        # the unit-copy test, on the purine strand
        bases = ga.purine_strand(seq[s - a : e - a], kind)
        word = "".join(ga.LETTER[int(x)] for x in bases)
        words = {word[j : j + K] for j in range(len(word) - K + 1)}
        left = (K - 1) // 2
        site_copy = {}
        for j in range(left, len(word) - (K - 1 - left)):
            w = word[j - left : j - left + K]
            cp = (w[:left] + TI[word[j]] + w[left + 1 :]) in words
            site_copy[j] = cp
            copy["sites"] += 1
            copy["sites_copy"] += cp
        for pos0, alts in snvs.items():
            j = pos0 - s if kind == "GA" else (e - 1 - pos0)
            if j not in site_copy:
                continue
            for alt, ac, af, an in alts:
                alt_p = alt if kind == "GA" else {"A": "T", "C": "G", "G": "C", "T": "A"}[alt]
                is_ti = alt_p == TI[word[j]]
                for g in ("all", af_class(ac, af)):
                    for x in (copy[g], copy["well_called"][g]) if an >= WELL_CALLED else (copy[g],):
                        if is_ti:
                            x["ti"] += 1
                            x["ti_copy"] += site_copy[j]
                        else:
                            x["tv"] += 1
    if write:
        write.close()
        os.replace(RECORDS + ".part", RECORDS)
    if len(counts) != len(bg) + len(tasks) or {r[:2] for r in counts} != {t[:2] for t in bg + tasks}:
        sys.exit(f"{len(counts):,} regions read, of {len(bg) + len(tasks):,}: the records are from another run")
    # expected substitutions inside, from the background's rates by context
    bp_bg, snv_bg = np.array(res["background"]["ctx_bp"]), np.array(res["background"]["ctx_snv"])
    rate = snv_bg / np.maximum(bp_bg, 1)[:, None]
    exp = np.array(res["inside"]["ctx_bp"])[:, None] * rate
    obs = np.array(res["inside"]["ctx_snv"])
    names = "ACGT"

    def types_of(obs, exp):
        types = {}
        for mb in (0, 2):
            for ab in range(4):
                if ab != mb:
                    sel = [c for c in range(64) if (c // 4) % 4 == mb]
                    types[f"{names[mb]}>{names[ab]}"] = {"observed": int(obs[sel, ab].sum()), "expected": round(float(exp[sel, ab].sum()), 1),
                                                         "ratio": round(float(obs[sel, ab].sum() / max(exp[sel, ab].sum(), 1e-9)), 3)}
        return types

    types = types_of(obs, exp)
    # the same in each allele-frequency class, from the background's rates in that class; transitions (A>G, G>A at a
    # purine middle base: alternative base 2 for an A, 0 for a G) and transversions apart
    by_af = {}
    ti_mask = np.zeros((64, 4), dtype=bool)
    for c in range(64):
        mb = (c // 4) % 4
        if mb in (0, 2):
            ti_mask[c, 2 - mb] = True
    for cname in CLASSES:
        ins, bg_ = res["inside"]["by_af"][cname], res["background"]["by_af"][cname]
        rate_c = np.array(bg_["ctx_snv"]) / np.maximum(bp_bg, 1)[:, None]
        exp_c = np.array(res["inside"]["ctx_bp"])[:, None] * rate_c
        obs_c = np.array(ins["ctx_snv"])
        row = {}
        for k, d in (("inside", ins), ("background", bg_)):
            bp = res[k]["bp"]
            row[k] = {"snv": d["snv"], "ti": d["ti"], "indel": d["indel"], "snv_per_kb": round(1000 * d["snv"] / max(bp, 1), 2),
                      "transitions": round(d["ti"] / max(d["snv"], 1), 4), "indel_per_kb": round(1000 * d["indel"] / max(bp, 1), 2),
                      "mean_an_snv": round(d["an"] / max(d["snv"], 1))}
        row["expected_inside"] = {
            "all": {"observed": int(obs_c.sum()), "expected": round(float(exp_c.sum()), 1), "ratio": round(float(obs_c.sum() / max(exp_c.sum(), 1e-9)), 3)},
            "transitions": {"observed": int(obs_c[ti_mask].sum()), "expected": round(float(exp_c[ti_mask].sum()), 1),
                            "ratio": round(float(obs_c[ti_mask].sum() / max(exp_c[ti_mask].sum(), 1e-9)), 3)},
            "transversions": {"observed": int(obs_c[~ti_mask].sum()), "expected": round(float(exp_c[~ti_mask].sum()), 1),
                              "ratio": round(float(obs_c[~ti_mask].sum() / max(exp_c[~ti_mask].sum(), 1e-9)), 3)},
            "types": types_of(obs_c, exp_c)}
        row["indels_inside_vs_background"] = round(row["inside"]["indel_per_kb"] / max(row["background"]["indel_per_kb"], 1e-9), 2)
        by_af[cname] = row
    p = copy["sites_copy"] / max(copy["sites"], 1)
    copy["expected_share"] = round(p, 4)
    for g in ("all", "singleton", "rare", "uncommon", "common"):
        for x in (copy[g], copy["well_called"][g]):
            if x["ti"]:
                x["share"] = round(x["ti_copy"] / x["ti"], 4)
                x["p"] = float(binomtest(x["ti_copy"], x["ti"], p).pvalue)
                x["transitions"] = round(x["ti"] / (x["ti"] + x["tv"]), 4)
    out = {"gnomad": "v4.1 genomes, PASS", "ga_sequences": len(tasks), "background_windows": len(bg),
           "background_windows_used": res["background"]["bp"] // 100,
           "inside": {k: v for k, v in res["inside"].items() if not k.startswith("ctx") and k != "by_af"},
           "background": {k: v for k, v in res["background"].items() if not k.startswith("ctx") and k != "by_af"},
           "expected_inside": {"observed": int(obs.sum()), "expected": round(float(exp.sum()), 1), "ratio": round(float(obs.sum() / exp.sum()), 3), "types": types},
           "by_allele_frequency": by_af,
           "unit_copy": copy}

    with gzip.open(os.path.join(ga.RESULTS, "gnomad_counts.tsv.gz" if not test else "gnomad_counts_test.tsv.gz"), "wt") as fh:
        fh.write("#chrom\tstart\tend\tkind\tpass_records\n")
        for row in sorted(counts):
            fh.write("\t".join(map(str, row)) + "\n")
    ga.write_json("gnomad_test.json" if test else "gnomad.json", out)
    for k in ("inside", "background"):
        d = res[k]
        print(f"{k}: {d['bp']:,} bp; SNVs {d['snv']:,} ({1000 * d['snv'] / max(d['bp'], 1):.0f}/kb), transitions {100 * d['ti'] / max(d['snv'], 1):.0f} %; indels {d['indel']:,} ({1000 * d['indel'] / max(d['bp'], 1):.1f}/kb)")
    print("expected inside:", out["expected_inside"])
    for cname, row in by_af.items():
        print(cname, row)
    print("unit copy:", copy)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) >= 3 else 8, sys.argv[3:] or None)
