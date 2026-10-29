"""How many published short-read de novo mutations fall in pure GA-sequences: the CEPH three-generation families of
Sasani et al. (2019, eLife; github.com/quinlan-lab/ceph-dnm-manuscript, data/second_gen.dnms.txt and
third_gen.dnms.txt, GRCh37), lifted to GRCh38 (Ensembl's GRCh37_to_GRCh38.chain.gz, with pyliftover), against
experiment 04's GA-sequences, and against their 400-base flanks.

    python3 -I short_reads.py DNM_DIR CHAIN         -> results/short_reads.json
"""

import bisect
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import ped  # noqa: E402

FLANK = 400


def main(dnm_dir: str, chain: str) -> None:
    from pyliftover import LiftOver

    lo = LiftOver(chain)
    runs = ped.read_runs()
    starts = {c: [r[0] for r in v] for c, v in runs.items()}

    def near(c: str, p: int, d: int) -> bool:
        v = runs.get(c, [])
        i = bisect.bisect_right(starts.get(c, []), p + d) - 1
        return any(0 <= j < len(v) and v[j][0] - d <= p < v[j][1] + d for j in (i, i - 1))

    with open(os.path.join(ped.EXP04, "results", "census.json")) as fh:
        census = json.load(fh)
    genome = sum(c["acgt"] for c in census["chromosomes"] if c["chrom"] in ped.AUTOSOMES)
    inside_bases = sum(e - s for c, v in runs.items() if c in ped.AUTOSOMES for s, e, _ in v)
    # the flanks' bases: each GA-sequence widened by FLANK either side, merged, less the GA-sequences
    flank_bases = 0
    for c, v in runs.items():
        if c not in ped.AUTOSOMES:
            continue
        cur_s = cur_e = None
        for s, e, _ in v:
            s2, e2 = s - FLANK, e + FLANK
            if cur_e is None or s2 > cur_e:
                if cur_e is not None:
                    flank_bases += cur_e - cur_s
                cur_s, cur_e = s2, e2
            else:
                cur_e = max(cur_e, e2)
        if cur_e is not None:
            flank_bases += cur_e - cur_s
    flank_bases -= inside_bases
    out = {"genome_bases": genome, "inside_bases": inside_bases, "flank_bases": flank_bases, "sets": {}}
    for name in ("second_gen", "third_gen"):
        n = lifted = inside = flank = 0
        kinds: dict = {}
        with open(os.path.join(dnm_dir, f"{name}.dnms.txt")) as fh:
            next(fh)
            for line in fh:
                x = line.rstrip("\n").split("\t")
                if x[0] in ("X", "Y"):
                    continue
                n += 1
                r = lo.convert_coordinate(x[0], int(x[1]))
                if not r:
                    continue
                lifted += 1
                c, p = "chr" + str(r[0][0]).removeprefix("chr"), int(r[0][1])
                if near(c, p, 0):
                    inside += 1
                    kinds[x[7]] = kinds.get(x[7], 0) + 1
                elif near(c, p, FLANK):
                    flank += 1
        out["sets"][name] = {"autosomal": n, "lifted": lifted, "inside": inside, "inside_kinds": kinds,
                             "expected_inside": round(lifted * inside_bases / genome, 1), "flanks": flank,
                             "expected_flanks": round(lifted * flank_bases / genome, 1)}
    tot = {k: sum(s[k] for s in out["sets"].values()) for k in ("autosomal", "lifted", "inside", "flanks")}
    tot["expected_inside"] = round(tot["lifted"] * inside_bases / genome, 1)
    tot["expected_flanks"] = round(tot["lifted"] * flank_bases / genome, 1)
    out["total"] = tot
    with open(os.path.join(ped.RESULTS, "short_reads.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
