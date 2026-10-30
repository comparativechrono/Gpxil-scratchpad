"""How much of HG002's excess of repeat-extending substitutions does the reference make by itself?

The reference decides where a repeat ends. Take a repeat whose ancestral sequence ends at base k, so that base k+1 is
the first that breaks it. A polymorphism shows as a heterozygous substitution at the reference's boundary base in three
ways, with no mutational bias at all:

  (i)   at k+1, the reference carrying the ancestral base: HG002's other allele is the new one, which continues the
        repeat with the chance q that the spectrum gives it;
  (ii)  at k+1, the reference carrying a new base that does not continue the repeat: HG002's other allele is the
        ancestral one, which doesn't either;
  (iii) at k, the reference carrying a new base: any new base breaks the repeat there, so the reference's repeat ends
        at k-1 and k is its boundary base; HG002's other allele is the ancestral one, which continues the repeat.

With derived-allele frequencies f spread as the neutral spectrum has them (density proportional to 1/f), HG002 is
heterozygous with chance 2f(1-f) and the reference carries the derived allele with chance f. Writing F for the mean of f
over heterozygous sites (1/3 under that spectrum), extending substitutions come out at (1 - F) q + F where the spectrum
expects q, a ratio of 1 - F + F/q (5/3 when q = 1/3), and the others at exactly their expected rate. A population
growing or shrinking moves F; a spectrum of density proportional to f^-a gives F = (2-a)/(4-a). A Monte Carlo of the
three cases checks the algebra. New mutations have none of this: a parent who carries the reference's base passes on a
new base that is a new mutation, wherever the repeat ends in other genomes.

    python3 -I ascertainment.py        -> results/ascertainment.json        (reads results/edges.json)
"""

from __future__ import annotations

import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
CLASSES = ["below 1e-6", "1e-6 to 1e-5", "1e-5 to 1e-4", "1e-4 or more"]


def ratio(q: float, a: float = 1.0) -> float:
    f = (2 - a) / (4 - a)
    return 1 - f + f / q


def monte_carlo(q: float, n: int = 2_000_000, two_n: int = 20_000, seed: int = 1) -> dict:
    """Polymorphisms at k and k+1 with derived-allele counts from the neutral spectrum in a sample of two_n
    haplotypes; the reference and HG002's two haplotypes drawn from them."""
    rng = np.random.default_rng(seed)
    i = np.arange(1, two_n)
    p = (1 / i) / (1 / i).sum()
    f = rng.choice(i, size=n, p=p) / two_n
    at_k = rng.random(n) < 0.5  # the polymorphism is at k (else at k+1)
    ref_derived = rng.random(n) < f
    h1, h2 = rng.random(n) < f, rng.random(n) < f
    het = h1 != h2
    continues = rng.random(n) < q  # at k+1, the new base continues the repeat
    # at the reference's boundary base: (i) k+1, reference ancestral; (ii) k+1, reference new and breaking; (iii) k, reference new
    i_ = ~at_k & ~ref_derived & het
    ii = ~at_k & ref_derived & ~continues & het
    iii = at_k & ref_derived & het
    ext = (i_ & continues).sum() + iii.sum()
    others = (i_ & ~continues).sum() + ii.sum()
    # what the spectrum expects at a site: its polymorphisms' heterozygous ones (the k+1 ones stand for any site's)
    expected = het[~at_k].sum()
    return {"extending": round(ext / (q * expected), 3), "others": round(others / ((1 - q) * expected), 3), "polymorphisms": n}


def main() -> None:
    out = {"q 1/3": {"neutral (a = 1)": round(ratio(1 / 3), 3), "a = 0.8": round(ratio(1 / 3, 0.8), 3), "a = 1.2": round(ratio(1 / 3, 1.2), 3),
                     "monte carlo, neutral": monte_carlo(1 / 3)}}
    with open(os.path.join(RESULTS, "edges.json")) as fh:
        edges = json.load(fh)
    rows = {}
    for cl in CLASSES:
        h = edges["loci, HG002"][f"{cl}, base 1, reference breaks the repeat"]
        q = h["q, the spectrum's chance of extending at these bases"]
        made = ratio(q)
        n = edges["loci, new substitutions"][f"{cl}, base 1"]
        rows[cl] = {"q, the spectrum's chance of extending": round(q, 3), "made by the reference": round(made, 4),
                    "range for a = 0.8 to 1.2": [round(ratio(q, 1.2), 3), round(ratio(q, 0.8), 3)],
                    "HG002, extending": h["extending"]["ratio"], "HG002, extending, beyond the reference's share": round(h["extending"]["ratio"] / made, 3),
                    "HG002, others": h["others"]["ratio"],
                    "new substitutions, extending": n["extending"]["ratio"], "new substitutions, extending, observed": n["extending"]["observed"],
                    "new substitutions, extending, expected": n["extending"]["expected"], "new substitutions, others": n["others"]["ratio"]}
    out["base 1, by class"] = rows
    out["monte carlo at each class's q"] = {cl: monte_carlo(r["q, the spectrum's chance of extending"]) for cl, r in rows.items()}
    with open(os.path.join(RESULTS, "ascertainment.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
