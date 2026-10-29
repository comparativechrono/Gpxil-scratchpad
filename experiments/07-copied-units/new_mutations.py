"""Do new substitutions in tandem repeats copy a neighbouring unit?

Experiment 04 found that HG002's substitutions in pure GA-sequences, and gnomAD's, mostly turn a stretch of the repeat
into one it already has, the more so the commoner the variant. Here the same question is put to new mutations, which
no selection has touched: the de novo substitutions of the CEPH three-generation families (Sasani et al. 2019, short
reads; GRCh37, lifted to GRCh38) and of the Platinum Pedigree (Porubsky et al. 2025, long reads; T2T-CHM13, lifted),
in the loci of the Platinum TRs catalogue.

A substitution at base j of a locus **copies** when the 6 bases around it (j - 2 to j + 3), after the change, occur
elsewhere in the locus. It takes the **neighbouring unit's base** when its new base is the base one period before or
after it (the period that makes the locus most regular, gpxiltools' period finder), and a **further unit's** when it is
the base two periods away. Against the locus's consensus unit (the commonest base at each place in the period) it moves
**toward** the consensus (an interruption put right), **away** from it (a new interruption), or neither.

Under the null, a new substitution in the locus falls at each base, and becomes each other base, at the rate of its
trinucleotide context: the de novo spectrum (Sasani et al.'s germline substitutions, against the autosomes'
trinucleotides). Each mutation's expectation is that rate-weighted share over its locus (or over the substitutions of
its own kind); observed counts are set against the sum of the expectations, and a normal approximation to the
Poisson-binomial gives z. The same is done with 5-base contexts: each trinucleotide's rate scaled by its 5-mer's new
mutations over those its trinucleotide predicts, shrunk toward 1 by 10 pseudo-mutations (keys ending in _5mer).

Platinum substitutions that Sasani et al. also called (the CEPH 1463 family is in both) are counted once, as Sasani's,
in the pooled tests. A substitution is in a **cluster** when another new mutation of the same child (a substitution,
an indel, or a Platinum tandem-repeat mutation) lies in its locus or within 10 bases of it.

Two controls take the 6-base words from other loci, five drawn at random for each mutation's locus, with the same
purity (bins 0.6, 0.8, 0.9) and length within 20 % (50 % when none is that close):

* **another repeat**: the same unit size (1 to 6, or 7 and more), any motif;
* **the same motif**: the same first motif in the catalogue, up to rotation and strand, read on the same strand. A
  mutational bias that came with the motif, rather than copying, would make the same words in loci of the same motif.

    python3 -I new_mutations.py REF_DIR FEATURES CATALOG DNM_DIR CHAIN37 [TABLE.xlsx CHAIN_T2T]
        -> results/new_mutations.json, results/new_mutations.tsv.gz

REF_DIR: the Gpxil browser's GRCh38 gpxil files. FEATURES: experiment 06's per-locus measures
(`mutability.py ... --features FILE`). CATALOG: the Platinum TRs catalogue (TRGT BED). DNM_DIR: ceph-dnm-manuscript's
data/ (second_gen, third_gen and post-pgcs .dnms.txt). CHAIN37: Ensembl's GRCh37_to_GRCh38.chain.gz. TABLE.xlsx and
CHAIN_T2T, optional: the Platinum Pedigree's Supplementary Table 10 and the T2T consortium's chm13v2-grch38.chain, as
experiment 05 reads them.
"""

from __future__ import annotations

import gzip
import json
import os
import random
import sys
from collections import Counter, defaultdict

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "05-ga-de-novo"))

import ped  # noqa: E402
import ga  # noqa: E402
from gpxiltools import codes as GC  # noqa: E402
from gpxiltools.period import best_period  # noqa: E402

RESULTS = os.path.join(HERE, "results")
K, LEFT = 6, 2  #: the window: 6 bases, the base third from the left
COMP = str.maketrans("ACGT", "TGCA")
PURITY = (0.0, 0.6, 0.8, 0.9, 1.01)
PICKS = 5  #: control loci for each mutation's locus
MAX_MOTIF = 50  #: longer first motifs get no same-motif control
NEAR = 10  #: another new mutation of the child's this close to a locus makes the substitution part of a cluster
SEED = 7
TAU = 10  #: pseudo-mutations that shrink each 5-mer's rate toward its trinucleotide's
LET = np.full(256, ord("N"), dtype=np.uint8)
for _b, _ch in ((GC.BASE_A, "A"), (GC.BASE_C, "C"), (GC.BASE_G, "G"), (GC.BASE_T, "T")):
    LET[_b] = LET[_b | GC.SOFTMASK] = ord(_ch)
CODE = {"A": GC.BASE_A, "C": GC.BASE_C, "G": GC.BASE_G, "T": GC.BASE_T}
KINDS = ("toward", "away", "neither")


def pyr(tri: str, alt: str) -> tuple[str, str]:
    """(trinucleotide, new base) read with a pyrimidine in the middle"""
    if tri[1] in "CT":
        return tri, alt
    return tri.translate(COMP)[::-1], alt.translate(COMP)


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def canon(m: str) -> tuple[str, int]:
    """A motif's canonical form (the least of its rotations and its reverse complement's) and the strand that gives it
    (1: the motif's, -1: the reverse complement's, 0: both)."""
    r = rc(m)
    a = min(m[i:] + m[:i] for i in range(len(m)))
    b = min(r[i:] + r[:i] for i in range(len(r)))
    return (a, 1) if a < b else (b, -1) if b < a else (a, 0)


def sasani(dnm_dir: str, chain: str) -> tuple[list[dict], list[tuple]]:
    """The new substitutions of the CEPH families on GRCh38, and every new mutation's place, (sample, chrom, start, end)."""
    from pyliftover import LiftOver

    lo = LiftOver(chain)
    snvs, places = [], []
    for name in ("second_gen", "third_gen", "post-pgcs"):
        with open(os.path.join(dnm_dir, f"{name}.dnms.txt")) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            for line in fh:
                r = dict(zip(head, line.rstrip("\n").split("\t")))
                if r["chrom"] in ("X", "Y"):
                    continue
                h = lo.convert_coordinate(r["chrom"], int(r["start"]))
                if not h or len(h) != 1:
                    continue
                chrom, pos, sample = "chr" + str(h[0][0]).removeprefix("chr"), int(h[0][1]), f"CEPH {r['new_sample_id']}"
                places.append((sample, chrom, pos, pos + max(len(r["ref"]), 1)))
                if len(r["ref"]) != 1 or len(r["alt"]) != 1:
                    continue
                ref, alt = r["ref"], r["alt"]
                if h[0][2] == "-":
                    ref, alt = ref.translate(COMP), alt.translate(COMP)
                snvs.append({"set": f"Sasani {name.replace('_', ' ').replace('-', ' ')}", "sample": sample, "chrom": chrom, "pos": pos,
                             "ref": ref, "alt": alt, "allele_balance": float(r["kid_allele_balance"]), "transmitted": None,
                             "family": r["new_family_id"]})
    return snvs, places


def platinum(table: str, chain: str, ref_dir: str) -> tuple[list[dict], list[tuple]]:
    """The Platinum Pedigree's new substitutions on GRCh38, and every new mutation's place (its indels lifted, its
    tandem-repeat mutations at their GRCh38 loci, whether they passed the paper's filters or not)."""
    import openpyxl

    from paper_calls import lift_and_read, sheet

    wb = openpyxl.load_workbook(table, read_only=True, data_only=True)
    rows = list(sheet(wb, "SNVs", 1))
    snvs, places = [], []
    for r, x in zip(rows, lift_and_read(rows, "T2T-CHM13v2.0 position", chain, ref_dir, "SNV")):
        if x["chrom"] in ped.AUTOSOMES:
            places.append((x["sample"], x["chrom"], x["pos"], x["pos"] + 1))
        if x["chrom"] in ped.AUTOSOMES and x.get("ref_matches"):
            ref, alt = x["ref"].upper(), x["alt"].upper()
            if x["strand"] == "-":
                ref, alt = ref.translate(COMP), alt.translate(COMP)
            kids = str(r.get("children with allele") or "").strip()
            ab = r.get("HiFi AB")
            snvs.append({"set": f"Platinum {x['origin']}", "sample": x["sample"], "chrom": x["chrom"], "pos": x["pos"], "ref": ref, "alt": alt,
                         "allele_balance": float(ab) if isinstance(ab, (int, float)) else None,
                         "transmitted": None if kids in ("", "NA", "None") else kids != "none"})
    ind = list(sheet(wb, "indels (non-TR)", 1))
    for x in lift_and_read(ind, "T2T-CHM13v2.0 position", chain, ref_dir, "indel"):
        if x["chrom"] in ped.AUTOSOMES:
            places.append((x["sample"], x["chrom"], x["pos"], x["pos"] + max(len(x["ref"]), 1)))
    for r in sheet(wb, "TRs (GRCh38)", 1):
        try:
            if r["chrom"] in ped.AUTOSOMES:
                places.append((r["sample_id"], r["chrom"], int(r["start"]), int(r["end"])))
        except (TypeError, ValueError):
            pass
    return snvs, places


class Loci:
    """The catalogue's loci with experiment 06's measures, as arrays in chromosome and start order."""

    def __init__(self, features: str):
        import pandas as pd

        df = pd.read_csv(features, sep="\t", usecols=["chrom", "start", "end", "motif_min", "purity", "class"],
                         dtype={"chrom": "category", "start": np.int64, "end": np.int64, "motif_min": np.int32, "purity": np.float64,
                                "class": "category"})
        self.chroms = list(ped.AUTOSOMES)
        ci = {c: i for i, c in enumerate(self.chroms)}
        lut = np.array([ci.get(c, -1) for c in df["chrom"].cat.categories], dtype=np.int16)
        chrom = lut[df["chrom"].cat.codes.to_numpy()]
        keep = chrom >= 0
        self.classes = list(df["class"].cat.categories)
        cols = {"chrom": chrom, "start": df["start"].to_numpy(), "end": df["end"].to_numpy(), "unit": df["motif_min"].to_numpy(),
                "purity": df["purity"].to_numpy(), "cls": df["class"].cat.codes.to_numpy()}
        del df
        order = np.lexsort((cols["start"][keep], cols["chrom"][keep]))
        for k, v in cols.items():
            setattr(self, k, v[keep][order])
        self.pbin = np.digitize(self.purity, PURITY[1:-1])
        self.length = self.end - self.start
        self.bounds = {c: (int(np.searchsorted(self.chrom, i, "left")), int(np.searchsorted(self.chrom, i, "right"))) for c, i in ci.items()}

    def at(self, chrom: str, pos: int) -> int:
        """The locus whose window positions hold pos (the window must fit in it), or -1."""
        if chrom not in self.bounds:
            return -1
        a, b = self.bounds[chrom]
        i = a + int(np.searchsorted(self.start[a:b], pos, "right")) - 1
        for j in range(max(a, i - 3), min(b, i + 2)):
            if self.start[j] + LEFT <= pos < self.end[j] - (K - 1 - LEFT):
                return j
        return -1

    def index(self, chrom: str, starts: np.ndarray, ends: np.ndarray) -> np.ndarray:
        """The loci with these starts and ends on chrom (-1 where none)."""
        a, b = self.bounds[chrom]
        i = np.clip(a + np.searchsorted(self.start[a:b], starts, "left"), a, max(a, b - 1))
        ok = (b > a) & (self.start[i] == starts) & (self.end[i] == ends)
        return np.where(ok, i, -1)

    def key(self, i: int) -> tuple[str, int, int]:
        return self.chroms[self.chrom[i]], int(self.start[i]), int(self.end[i])


def first_motifs(catalog: str, keys: set) -> dict:
    """The catalogue's id, first motif and number of motifs of each locus in keys ((chrom, start, end))."""
    chroms = {k[0] for k in keys}
    out = {}
    with gzip.open(catalog, "rt") as fh:
        for line in fh:
            c, s, e, info = line.rstrip("\n").split("\t")[:4]
            if c not in chroms:
                continue
            k = (c, int(s), int(e))
            if k in keys:
                f = dict(x.split("=", 1) for x in info.split(";"))
                motifs = [m.upper() for m in f["MOTIFS"].split(",") if m]
                out[k] = (f["ID"], motifs[0], len(motifs))
    return out


def same_motif_loci(catalog: str, wanted: dict) -> dict:
    """The loci whose first motif is a rotation of a wanted canonical motif or of its reverse complement, with a length
    in the wanted range: canonical motif -> [(chrom, start, end, strand)]."""
    rot = {}
    for m in wanted:
        for x in (m, rc(m)):
            for i in range(len(x)):
                r = x[i:] + x[:i]
                rot[r] = canon(r)
    out: dict = defaultdict(list)
    with gzip.open(catalog, "rt") as fh:
        for line in fh:
            c, s, e, info = line.rstrip("\n").split("\t")[:4]
            i = info.find("MOTIFS=")
            j = info.find(";", i)
            m = info[i + 7:j if j >= 0 else len(info)].split(",")[0].upper()
            hit = rot.get(m)
            if hit:
                s, e = int(s), int(e)
                lo, hi = wanted[hit[0]]
                if lo <= e - s <= hi:
                    out[hit[0]].append((c, s, e, hit[1]))
    return out


def pick(rng: random.Random, idx: np.ndarray, n: int) -> list[int]:
    if len(idx) <= n:
        return [int(x) for x in idx]
    return [int(x) for x in rng.sample(list(idx), n)]


def expectations(seq: str, period: int, cons: list | None, rates: dict, own: set, others: dict) -> dict:
    """Rate-weighted shares, over the locus's possible substitutions, of copying (into the locus's own words and each
    control's), of taking the base one and two periods away, and of each kind; and the copying shares within each kind.
    ``rates``: name -> rate(seq, x, new base); the first is the main one, and only it is applied to the controls; the
    others' keys end in _name."""
    first = next(iter(rates))
    acc = {name: [0.0, Counter(), Counter(), defaultdict(Counter)] for name in rates}  # den, num, by kind, copies by kind
    n = len(seq)
    for x in range(max(LEFT, 1), n - (K - 1 - LEFT)):
        win = seq[x - LEFT:x - LEFT + K]
        near = {seq[y] for y in (x - period, x + period) if 0 <= y < n}
        far = {seq[y] for y in (x - 2 * period, x + 2 * period) if 0 <= y < n}
        c = cons[x % period] if cons else None
        for a in "ACGT":
            if a == seq[x]:
                continue
            new = win[:LEFT] + a + win[LEFT + 1:]
            kind = None if c is None else "away" if seq[x] == c else "toward" if a == c else "neither"
            hits = {"copies": new in own, "near": a in near, "far": a in far}
            for name, words in others.items():
                hits[name] = new in words
            for name, fn in rates.items():
                w = fn(seq, x, a)
                if not w:
                    continue
                d = acc[name]
                d[0] += w
                for h, v in hits.items():
                    if v and (name == first or h in ("copies", "near", "far")):
                        d[1][h] += w
                if kind:
                    d[2][kind] += w
                    for h, v in hits.items():
                        if v and (name == first or h == "copies"):
                            d[3][kind][h] += w
    out = {}
    for name, (den, num, dk, nk) in acc.items():
        if not den:
            if name == first:
                return {}
            continue
        sfx = "" if name == first else f"_{name}"
        hs = ["copies", "near", "far"] + (list(others) if name == first else [])
        out.update({f"expected_{h}{sfx}": num[h] / den for h in hs})
        if dk:
            tot = sum(dk.values())
            for kind in KINDS:
                out[f"expected_{kind}{sfx}"] = dk[kind] / tot
            out[f"in_kind{sfx}"] = {kind: {h: nk[kind][h] / dk[kind] for h in ["copies"] + (list(others) if name == first else [])}
                                     for kind in KINDS if dk[kind]}
    return out


def pentamers(ref_dir: str, chunk: int = 1 << 24) -> np.ndarray:
    """The autosomes' 5-base words (1,024, the first base most significant)."""
    from paper_calls import IDX

    idx = IDX.astype(np.int8)
    out = np.zeros(1024, dtype=np.int64)
    for ch in ped.AUTOSOMES:
        with ga.reference(ref_dir, ch) as g:
            c = g[ch]
            for a in range(0, c.length - 4, chunk):
                b = min(c.length, a + chunk + 4)
                x = idx[c.read("ref", a, b)].astype(np.int16)
                ok = (x[:-4] >= 0) & (x[1:-3] >= 0) & (x[2:-2] >= 0) & (x[3:-1] >= 0) & (x[4:] >= 0)
                code = x[:-4] * 256 + x[1:-3] * 64 + x[2:-2] * 16 + x[3:-1] * 4 + x[4:]
                out += np.bincount(code[ok], minlength=1024)
    return out


def consensus(seq: str, period: int) -> list | None:
    """The commonest base at each place in the period, if the locus has three units or more."""
    if period < 1 or len(seq) < 3 * period:
        return None
    out = []
    for f in range(period):
        c = Counter(ch for ch in seq[f::period] if ch != "N")
        out.append(c.most_common(1)[0][0] if c else "N")
    return out


def main(ref_dir: str, features: str, catalog: str, dnm_dir: str, chain37: str, table: str | None = None, chain_t2t: str | None = None) -> None:
    rng = random.Random(SEED)
    snvs, places = sasani(dnm_dir, chain37)
    platinum_keys: set = set()
    if table and chain_t2t:
        p_snvs, p_places = platinum(table, chain_t2t, ref_dir)
        seen_sasani = {(d["chrom"], d["pos"], d["alt"]) for d in snvs}
        for d in p_snvs:
            d["also_in_sasani"] = (d["chrom"], d["pos"], d["alt"]) in seen_sasani
        platinum_keys = {(d["chrom"], d["pos"], d["alt"]) for d in p_snvs}
        snvs += p_snvs
        places += p_places
    print(f"{len(snvs):,} new substitutions read", flush=True)
    loci = Loci(features)
    print(f"{len(loci.start):,} loci", flush=True)
    rows = []
    for d in snvs:
        d["locus"] = loci.at(d["chrom"], d["pos"])
        if d["locus"] >= 0:
            rows.append(d)
    # clusters: another new mutation of the same child in the locus or within NEAR bases of it
    by_child: dict = defaultdict(list)
    for sample, c, a, b in places:
        by_child[(sample, c)].append((a, b))
    for r in rows:
        _, s0, e0 = loci.key(r["locus"])
        r["cluster"] = any(a < e0 + NEAR and b > s0 - NEAR and not (a == r["pos"] and b == r["pos"] + 1)
                           for a, b in by_child[(r["sample"], r["chrom"])])
    print(f"{len(rows):,} in loci", flush=True)
    # the mutated loci's first motifs, and loci of the same motif
    keys = {loci.key(r["locus"]) for r in rows}
    motif = first_motifs(catalog, keys)
    wanted: dict = {}
    for k in keys:
        if k in motif and len(motif[k][1]) <= MAX_MOTIF and set(motif[k][1]) <= set("ACGT"):
            cm, _ = canon(motif[k][1])
            lo, hi = wanted.get(cm, (10 ** 9, 0))
            wanted[cm] = (min(lo, int((k[2] - k[1]) * 0.5)), max(hi, int((k[2] - k[1]) * 1.5) + 1))
    same = same_motif_loci(catalog, wanted)
    groups = {}
    for cm, v in same.items():
        cs = np.array([x[0] for x in v])
        st = np.array([x[1] for x in v], dtype=np.int64)
        en = np.array([x[2] for x in v], dtype=np.int64)
        strand = np.array([x[3] for x in v], dtype=np.int8)
        idx = np.full(len(v), -1, dtype=np.int64)
        for c in set(cs.tolist()):
            if c in loci.bounds:
                m = cs == c
                idx[m] = loci.index(c, st[m], en[m])
        ok = idx >= 0
        groups[cm] = (idx[ok], strand[ok])
    print(f"same-motif candidates for {len(groups):,} motifs", flush=True)
    # the controls of each mutated locus
    pool: dict = defaultdict(lambda: np.array([], dtype=np.int64))
    unit_class = np.minimum(loci.unit, 7)
    for u in range(1, 8):
        for b in range(len(PURITY) - 1):
            idx = np.nonzero((unit_class == u) & (loci.pbin == b))[0]
            pool[(u, b)] = idx[np.argsort(loci.length[idx], kind="stable")]
    controls: dict = {}
    for k in sorted(keys):
        i = loci.index(k[0], np.array([k[1]]), np.array([k[2]]))[0]
        L, b = int(loci.length[i]), int(loci.pbin[i])
        cand = pool[(min(int(loci.unit[i]), 7), b)]
        lens = loci.length[cand]
        chosen = []
        for f in (0.2, 0.5):
            sel = cand[np.searchsorted(lens, int(L * (1 - f)), "left"):np.searchsorted(lens, int(L * (1 + f)), "right")]
            sel = sel[sel != i]
            if len(sel):
                chosen = pick(rng, sel, PICKS)
                break
        other = [(j, 1) for j in chosen]
        sm = []
        if k in motif and len(motif[k][1]) <= MAX_MOTIF and set(motif[k][1]) <= set("ACGT"):
            cm, strand = canon(motif[k][1])
            if cm in groups:
                gi, gs = groups[cm]
                for f in (0.2, 0.5):
                    m = (gi != i) & (loci.pbin[gi] == b) & (loci.length[gi] >= int(L * (1 - f))) & (loci.length[gi] <= int(L * (1 + f)))
                    if m.any():
                        sel = np.nonzero(m)[0]
                        for t in pick(rng, sel, PICKS):
                            # read the control on the strand on which its motif runs as the mutated locus's does
                            flip = strand != 0 and gs[t] != 0 and gs[t] != strand
                            sm.append((int(gi[t]), -1 if flip else 1))
                        break
        controls[k] = {"other_repeat": other, "same_motif": sm, "motif": motif.get(k, (None, None, 0))}
    # one pass over the reference: the mutations' bases, the spectrum, and the loci's sequences
    need: dict = defaultdict(set)
    for k, v in controls.items():
        need[k[0]].add((k[1], k[2]))
        for name in ("other_repeat", "same_motif"):
            for j, _ in v[name]:
                c, s, e = loci.key(j)
                need[c].add((s, e))
    by_chrom: dict = defaultdict(list)
    for d in snvs:
        by_chrom[d["chrom"]].append(d)
    spectrum, spectrum5, seen, seqs = Counter(), Counter(), Counter(), {}
    for c in ped.AUTOSOMES:
        with ga.reference(ref_dir, c) as g:
            s = LET[g[c].read("ref")].tobytes().decode()
        for d in by_chrom.get(c, []):
            p = d["pos"]
            d["ref_ok"] = s[p] == d["ref"]
            seen["checked" if d["ref_ok"] else "reference base differs"] += 1
            tri = s[p - 1:p + 2]
            if d["ref_ok"] and d["set"] in ("Sasani second gen", "Sasani third gen") and "N" not in tri:
                spectrum[pyr(tri, d["alt"])] += 1
                five = s[p - 2:p + 3]
                if "N" not in five:
                    spectrum5[pyr(five, d["alt"])] += 1
        for s0, e0 in need.get(c, ()):
            seqs[(c, s0, e0)] = s[s0:e0]
        del s
        print(f"{c} read", flush=True)
    from paper_calls import Runs, contexts_of_genome

    allc, _ = contexts_of_genome(ref_dir, Runs())
    genome = Counter()
    for i in range(64):
        genome[pyr("ACGT"[i // 16] + "ACGT"[(i // 4) % 4] + "ACGT"[i % 4], "A")[0]] += int(allc[i])

    def rate(tri: str, a: str) -> float:
        t, x = pyr(tri, a)
        return spectrum[(t, x)] / genome[t] if genome[t] else 0.0

    # 5-base contexts: the trinucleotide rate, scaled by the 5-mer's observed over expected new mutations, shrunk
    # toward 1 by TAU pseudo-mutations
    p5 = pentamers(ref_dir)
    genome5 = Counter()
    for i in range(1024):
        w5 = "".join("ACGT"[(i >> (2 * (4 - k))) & 3] for k in range(5))
        genome5[pyr(w5, "A")[0]] += int(p5[i])
    scale5 = {}
    for w5, g5 in genome5.items():
        for a in "ACGT":
            if a == w5[2]:
                continue
            e5 = g5 * rate(w5[1:4], a)
            scale5[(w5, a)] = (spectrum5[(w5, a)] + TAU) / (e5 + TAU)

    def rate3(seq: str, x: int, a: str) -> float:
        tri = seq[x - 1:x + 2]
        return 0.0 if "N" in tri else rate(tri, a)

    def rate5(seq: str, x: int, a: str) -> float:
        five = seq[x - 2:x + 3]
        if len(five) < 5 or "N" in five:
            return 0.0
        w5, b = pyr(five, a)
        return rate(five[1:4], a) * scale5[(w5, b)]

    rates = {"tri": rate3, "5mer": rate5}

    def words(seq: str) -> set:
        return {seq[x:x + K] for x in range(len(seq) - K + 1)}

    out = []
    for r in rows:
        if not r.get("ref_ok"):
            continue
        k = loci.key(r["locus"])
        i = r["locus"]
        seq = seqs[k]
        j = r["pos"] - k[1]
        ctl = controls[k]
        _, first, n_motifs = ctl["motif"]
        maxp = min(max(12, len(first or "")), 100)
        period = best_period(np.array([CODE.get(ch, 15) for ch in seq], dtype=np.uint8), maxp)[0]
        cons = consensus(seq, period)
        others = {}
        for name in ("other_repeat", "same_motif"):
            for t, (jj, strand) in enumerate(ctl[name]):
                cseq = seqs[loci.key(jj)]
                others[f"{name}_{t}"] = words(cseq if strand == 1 else rc(cseq))
        own = words(seq)
        e = expectations(seq, period, cons, rates, own, others)
        if not e:
            continue
        win = seq[j - LEFT:j - LEFT + K]
        new = win[:LEFT] + r["alt"] + win[LEFT + 1:]
        c = cons[j % period] if cons else None
        kind = None if c is None else "away" if r["ref"] == c else "toward" if r["alt"] == c else "neither"
        row = {key: r.get(key) for key in ("set", "sample", "chrom", "pos", "ref", "alt", "allele_balance", "transmitted", "cluster", "also_in_sasani")}
        row.update({"trid": ctl["motif"][0], "motif": canon(first)[0] if first and len(first) <= MAX_MOTIF and set(first) <= set("ACGT") else first,
                    "motifs": n_motifs, "unit": int(loci.unit[i]), "purity": float(loci.purity[i]), "class": loci.classes[loci.cls[i]],
                    "length": len(seq), "period": period, "transition": {r["ref"], r["alt"]} in ({"A", "G"}, {"C", "T"}), "kind": kind,
                    "copies": new in own, "expected_copies": e["expected_copies"],
                    "near": r["alt"] in {seq[y] for y in (j - period, j + period) if 0 <= y < len(seq)}, "expected_near": e["expected_near"],
                    "far": r["alt"] in {seq[y] for y in (j - 2 * period, j + 2 * period) if 0 <= y < len(seq)}, "expected_far": e["expected_far"]})
        for kk in KINDS:
            if f"expected_{kk}" in e:
                row[f"expected_{kk}"] = e[f"expected_{kk}"]
        for key in ("copies", "near", "far") + KINDS:
            if f"expected_{key}_5mer" in e:
                row[f"expected_{key}_5mer"] = e[f"expected_{key}_5mer"]
        ik5 = e.get("in_kind_5mer", {}).get(kind, {}) if kind else {}
        if kind and ik5:
            row["expected_copies_in_kind_5mer"] = ik5["copies"]
        ik = e.get("in_kind", {}).get(kind, {}) if kind else {}
        if kind and ik:
            row["expected_copies_in_kind"] = ik["copies"]
        for name in ("other_repeat", "same_motif"):
            names = [f"{name}_{t}" for t in range(len(ctl[name]))]
            if names:
                row[f"{name}_copies"] = sum(new in others[x] for x in names) / len(names)
                row[f"expected_{name}_copies"] = sum(e[f"expected_{x}"] for x in names) / len(names)
                if kind and ik:
                    row[f"expected_{name}_copies_in_kind"] = sum(ik[x] for x in names) / len(names)
        out.append(row)
    os.makedirs(RESULTS, exist_ok=True)
    cols = ["set", "sample", "chrom", "pos", "ref", "alt", "allele_balance", "transmitted", "cluster", "also_in_sasani", "trid", "motif", "motifs",
            "unit", "purity", "class", "length", "period", "transition", "kind", "expected_toward", "expected_away", "expected_neither",
            "copies", "expected_copies", "expected_copies_in_kind", "near", "expected_near", "far", "expected_far",
            "other_repeat_copies", "expected_other_repeat_copies", "expected_other_repeat_copies_in_kind",
            "same_motif_copies", "expected_same_motif_copies", "expected_same_motif_copies_in_kind",
            "expected_copies_5mer", "expected_copies_in_kind_5mer", "expected_near_5mer", "expected_far_5mer", "expected_toward_5mer",
            "expected_away_5mer", "expected_neither_5mer"]
    with gzip.open(os.path.join(RESULTS, "new_mutations.tsv.gz"), "wt") as fh:
        fh.write("\t".join(cols) + "\n")
        for x in out:
            fh.write("\t".join(fmt(x.get(c)) for c in cols) + "\n")
    summary = {"checked": dict(seen), "spectrum_mutations": sum(spectrum.values()), "spectrum_mutations_5mer": sum(spectrum5.values()),
               "pseudo_mutations": TAU, "in_loci": len(out), "platinum_family": platinum_family(snvs, out, platinum_keys),
               "with_controls": {"other_repeat": sum("other_repeat_copies" in x for x in out), "same_motif": sum("same_motif_copies" in x for x in out)},
               "tests": tests(out)}
    with open(os.path.join(RESULTS, "new_mutations.json"), "w") as fh:
        json.dump(summary, fh, indent=1)
    for name, t in summary["tests"].items():
        print(name, json.dumps(t))


def platinum_family(snvs: list[dict], out: list[dict], keys: set) -> dict:
    """The CEPH 1463 family is in both studies: Sasani et al.'s family whose new substitutions the Platinum Pedigree's
    calls hold most often, and how many of its substitutions they hold, by generation, and in the loci by copying."""
    if not keys:
        return {}
    sas = [d for d in snvs if d["set"].startswith("Sasani")]
    fam = Counter(d["family"] for d in sas if (d["chrom"], d["pos"], d["alt"]) in keys).most_common(1)[0][0]
    mine = [d for d in sas if d["family"] == fam]
    by_set: dict = {}
    for d in mine:
        v = by_set.setdefault(d["set"], [0, 0])
        v[0] += (d["chrom"], d["pos"], d["alt"]) in keys
        v[1] += 1
    samples = {d["sample"] for d in mine}
    loci: dict = {}
    for x in out:
        if x["set"].startswith("Sasani") and x["sample"] in samples:
            v = loci.setdefault("copies" if x["copies"] else "others", [0, 0])
            v[0] += (x["chrom"], x["pos"], x["alt"]) in keys
            v[1] += 1
    return {"sasani_family": fam, "held_of_all": by_set, "held_in_loci": loci}


def fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def test(sel: list[dict], obs_key: str, exp_key: str) -> dict:
    """Observed against expected: the observations are 0/1 (or a control's mean over its loci), the expectations shares."""
    sel = [x for x in sel if x.get(exp_key) is not None and x.get(obs_key) is not None]
    n = len(sel)
    if not n:
        return {"mutations": 0}
    obs = sum(float(x[obs_key]) for x in sel)
    exp = sum(x[exp_key] for x in sel)
    var = sum(x[exp_key] * (1 - x[exp_key]) for x in sel)
    z = (obs - exp) / var ** 0.5 if var else None
    return {"mutations": n, "observed": round(obs, 1), "expected": round(exp, 1), "ratio": round(obs / exp, 2) if exp else None,
            "observed_share": round(obs / n, 3), "expected_share": round(exp / n, 3), "z": round(z, 1) if z is not None else None}


def against_control(sel: list[dict], control: str, suffix: str = "") -> dict:
    """The excess of copying the locus's own words over the excess of copying a control's, mutation by mutation (the
    mutations with that control): its z ignores the two's covariance, which is positive, and so is conservative."""
    own_e = "expected_copies" + suffix
    c_o, c_e = f"{control}_copies", f"expected_{control}_copies" + suffix
    sel = [x for x in sel if x.get(own_e) is not None and x.get(c_e) is not None]
    if not sel:
        return {"mutations": 0}
    d = sum(float(x["copies"]) - x[own_e] for x in sel) - sum(x[c_o] - x[c_e] for x in sel)
    var = sum(x[own_e] * (1 - x[own_e]) + x[c_e] * (1 - x[c_e]) for x in sel)
    own = sum(float(x["copies"]) for x in sel) / sum(x[own_e] for x in sel)
    ctl = sum(x[c_o] for x in sel) / sum(x[c_e] for x in sel)
    return {"mutations": len(sel), "own_ratio": round(own, 2), "control_ratio": round(ctl, 2), "excess_difference": round(d, 1),
            "z": round(d / var ** 0.5, 1) if var else None}


def battery(sel: list[dict]) -> dict:
    """The tests of one group of mutations."""
    away = [x for x in sel if x["kind"] == "away"]
    phased = [x for x in sel if x["kind"]]
    return {
        "copies": test(sel, "copies", "expected_copies"),
        "other_repeat": test(sel, "other_repeat_copies", "expected_other_repeat_copies"),
        "same_motif": test(sel, "same_motif_copies", "expected_same_motif_copies"),
        "own_vs_same_motif": against_control(sel, "same_motif"),
        "near": test(sel, "near", "expected_near"),
        "far": test(sel, "far", "expected_far"),
        "toward": test([dict(x, is_toward=x["kind"] == "toward") for x in phased], "is_toward", "expected_toward"),
        "away": test([dict(x, is_away=x["kind"] == "away") for x in phased], "is_away", "expected_away"),
        "away_copies": test(away, "copies", "expected_copies_in_kind"),
        "away_same_motif": test(away, "same_motif_copies", "expected_same_motif_copies_in_kind"),
        "away_other_repeat": test(away, "other_repeat_copies", "expected_other_repeat_copies_in_kind"),
        "away_own_vs_same_motif": against_control(away, "same_motif", "_in_kind"),
        "5mer_copies": test(sel, "copies", "expected_copies_5mer"),
        "5mer_near": test(sel, "near", "expected_near_5mer"),
        "5mer_far": test(sel, "far", "expected_far_5mer"),
        "5mer_toward": test([dict(x, is_toward=x["kind"] == "toward") for x in phased], "is_toward", "expected_toward_5mer"),
        "5mer_away_copies": test(away, "copies", "expected_copies_in_kind_5mer"),
    }


def tests(out: list[dict]) -> dict:
    res = {}
    once = [x for x in out if not x.get("also_in_sasani")]
    for lo, hi in zip(PURITY, PURITY[1:]):
        for unit, keep in (("1-6", lambda u: u <= 6), ("7+", lambda u: u > 6)):
            res[f"purity {lo:.1f}-{min(hi, 1):.1f}, unit {unit}"] = battery([x for x in once if lo <= x["purity"] < hi and keep(x["unit"])])
    pure = [x for x in once if x["purity"] >= 0.6]
    groups = {"purity 0.6+": pure, "purity 0.6+, unit 1-6": [x for x in pure if x["unit"] <= 6],
              "purity 0.6+, unit 7+": [x for x in pure if x["unit"] > 6],
              "purity 0.6+, outside clusters": [x for x in pure if not x["cluster"]],
              "purity 0.6+, in clusters": [x for x in pure if x["cluster"]],
              "purity 0.6+, transitions": [x for x in pure if x["transition"]],
              "purity 0.6+, transversions": [x for x in pure if not x["transition"]],
              "purity 0.6+, one motif": [x for x in pure if x["motifs"] == 1],
              "GA-sequences": [x for x in once if x["class"] == "GA-sequence"]}
    for lo, hi in ((0, 0.35), (0.35, 0.45), (0.45, 1.01)):
        groups[f"purity 0.6+, allele balance {lo}-{min(hi, 1)}"] = [x for x in pure if x["allele_balance"] is not None and lo <= x["allele_balance"] < hi]
    for s in sorted({x["set"] for x in out}):
        groups[f"purity 0.6+, {s}"] = [x for x in out if x["purity"] >= 0.6 and x["set"] == s]
    groups["purity 0.6+, Platinum transmitted"] = [x for x in out if x["purity"] >= 0.6 and x["transmitted"]]
    for name, sel in groups.items():
        res[name] = battery(sel)
    return res


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.exit(__doc__)
    main(*sys.argv[1:8])
