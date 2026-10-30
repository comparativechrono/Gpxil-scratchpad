"""How often does a new substitution keep a pure run pure, in each alphabet?

A substitution inside a pure run keeps it pure when the new base is in the run's alphabet: in a purine run, a
transition (A<->G, which is C<->T on the other strand); in an A/T run, A<->T; in an A/C run, A<->C; in a G/C run, G<->C.
The new substitutions of the CEPH families (Sasani et al. 2019, lifted to GRCh38) and of the Platinum Pedigree
(Porubsky et al. 2025), as experiment 07 reads them (those both call, once), give the spectrum: at A:T pairs, the shares of T>C, T>A and T>G
(read with a pyrimidine in the middle); at C:G pairs outside CpG and in CpG, of C>T, C>A and C>G; and the rate at C:G
outside CpG against A:T, per base, from the autosomes' base counts (results/alphabets.json). Purine and A/C runs hold
no CpG (they lack C or G on one strand); G/C runs can, and are given both ways.

    python3 -I spectrum.py REF_DIR DNM_DIR CHAIN37 TABLE.xlsx CHAIN_T2T     -> results/spectrum.json
"""

from __future__ import annotations

import json
import os
import sys
from collections import Counter

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(EXP, "04-ga-sequences"))
sys.path.insert(0, os.path.join(EXP, "07-copied-units"))

import ga  # noqa: E402
import new_mutations as nm  # noqa: E402

RESULTS = os.path.join(HERE, "results")
LETTERS = "-ACMGRSVTWYHKDBN"
COMP = str.maketrans("ACGT", "TGCA")


def spectrum(snvs: list[dict], ref_dir: str) -> tuple[Counter, int]:
    """(source, pyrimidine reference base, new base, in CpG) -> new substitutions."""
    cnt, bad = Counter(), 0
    by = {}
    for d in snvs:
        by.setdefault(d["chrom"], []).append(d)
    for c, v in by.items():
        with ga.reference(ref_dir, c) as g:
            ref = g[c]
            for d in v:
                tri = "".join(LETTERS[x & 15] for x in ref.read("ref", d["pos"] - 1, d["pos"] + 2))
                if tri[1] != d["ref"] or any(b not in "ACGT" for b in tri):
                    bad += 1
                    continue
                r, a = d["ref"], d["alt"]
                if r in "AG":
                    tri, r, a = tri.translate(COMP)[::-1], r.translate(COMP), a.translate(COMP)
                cnt[(d["set"].split()[0], r, a, r == "C" and tri[2] == "G")] += 1
    return cnt, bad


def kept(cnt: Counter, source: str, bases: dict) -> dict:
    def k(r, a, cpg):
        return sum(v for (s, rr, aa, cc), v in cnt.items() if source in ("all", s) and rr == r and aa == a and cc == cpg)

    at = {a: k("T", a, False) for a in "CAG"}
    cg = {a: k("C", a, False) for a in "TAG"}
    n_at, n_cg = sum(at.values()), sum(cg.values())
    s_at = {a: at[a] / n_at for a in at}
    s_cg = {a: cg[a] / n_cg for a in cg}
    ratio = (n_cg / bases["C:G outside CpG"]) / (n_at / bases["A:T"])  # rate per base, C:G against A:T

    def two(f_at: float, keep_at: float, keep_cg: float) -> dict:
        """A run with a share f_at of its bases at A:T pairs: the share of its substitutions that keep it pure, and
        the substitutions that break it, per base, against an A:T pair's rate of all substitutions."""
        w_at, w_cg = f_at, (1 - f_at) * ratio
        keep = (w_at * keep_at + w_cg * keep_cg) / (w_at + w_cg)
        return {"kept": round(keep, 6), "broken_per_base": round((w_at + w_cg) * (1 - keep), 6)}

    cpg = {a: k("C", a, True) for a in "TAG"}
    n_cpg = sum(cpg.values())
    transitions = at["C"] + cg["T"] + cpg["T"]
    return {
        "new substitutions": {"at A:T": n_at, "at C:G outside CpG": n_cg, "at CpG": n_cpg},
        "transitions, share of all": round(transitions / (n_at + n_cg + n_cpg), 6),
        "shares at A:T (T>C, T>A, T>G)": [round(s_at[a], 6) for a in "CAG"],
        "shares at C:G outside CpG (C>T, C>A, C>G)": [round(s_cg[a], 6) for a in "TAG"],
        "shares at CpG (C>T, C>A, C>G)": [round(cpg[a] / n_cpg, 6) for a in "TAG"] if n_cpg else None,
        "rate per base, C:G outside CpG against A:T": round(ratio, 6),
        "counts": {"at A:T (T>C, T>A, T>G)": [at[a] for a in "CAG"], "at C:G outside CpG (C>T, C>A, C>G)": [cg[a] for a in "TAG"],
                   "at CpG (C>T, C>A, C>G)": [cpg[a] for a in "TAG"]},
        "runs": {"purine, 3 A to 1 G": two(0.75, s_at["C"], s_cg["T"]), "purine, 1 A to 1 G": two(0.5, s_at["C"], s_cg["T"]),
                 "A/T": two(1.0, s_at["A"], 0.0), "A/C, 3 A to 1 C": two(0.75, s_at["G"], s_cg["A"]),
                 "A/C, 1 A to 1 C": two(0.5, s_at["G"], s_cg["A"]), "G/C, outside CpG": two(0.0, 0.0, s_cg["G"]),
                 "G/C, at CpG (kept share only)": {"kept": round(cpg["G"] / n_cpg, 6)} if n_cpg else None}}


def main(ref_dir: str, dnm_dir: str, chain37: str, table: str, chain_t2t: str) -> None:
    snvs, _ = nm.sasani(dnm_dir, chain37)
    seen = {(d["chrom"], d["pos"], d["alt"]) for d in snvs}
    p, _ = nm.platinum(table, chain_t2t, ref_dir)
    snvs += [d for d in p if (d["chrom"], d["pos"], d["alt"]) not in seen]  # (the CEPH 1463 family is in both: counted once)
    cnt, bad = spectrum(snvs, ref_dir)
    with open(os.path.join(RESULTS, "alphabets.json")) as fh:
        di = np.array(json.load(fh)["dinucleotides_autosomes"])[:4, :4]
    a, c, g, t = di.sum(1)
    bases = {"A:T": int(a + t), "C:G outside CpG": int(c + g - 2 * di[1][2])}
    out = {"read": len(snvs), "left out (reference differs, or N)": bad, "autosome bases": bases}
    for source, key in (("all", "all"), ("Sasani", "Sasani"), ("Platinum", "Platinum, not also Sasani's")):
        out[key] = kept(cnt, source, bases)
        print(key, json.dumps(out[key]), flush=True)
    with open(os.path.join(RESULTS, "spectrum.json"), "w") as fh:
        json.dump(out, fh, indent=1)


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:6])
