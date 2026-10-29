"""Are the new substitutions that copy new alleles? gnomAD v4.1 genomes at each one's place.

A true new mutation is mostly absent from a population sample: of a random 3,000 of the CEPH families' new
substitutions outside the catalogue's loci, most are not in gnomAD's 76,215 genomes, and those that are, are mostly
rare. Here every new substitution of new_mutations.py's table (the ones in the catalogue's loci) and that random 3,000
are looked up in gnomAD's sites VCFs (every record at the place, whatever its FILTER), and the within-locus tests of
new_mutations.py are split by whether the new allele is there: absent; present, and then by FILTER and frequency.

    CURL_CA_BUNDLE=... python3 -I in_gnomad.py DNM_DIR CHAIN37 RECORDS.tsv.gz      -> results/in_gnomad.json,
                                                                                    results/in_gnomad.tsv.gz

RECORDS.tsv.gz: the records read, kept (and read from there when it exists): chrom, position (0-based), then each
record's REF, ALT, FILTER, AC, AF, AN, or a line with only the place when gnomAD has none.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import random
import sys
import time
from collections import Counter
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import new_mutations as nm  # noqa: E402

URL = "https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.{}.vcf.bgz"
OTHERS = 3000
SEED = 11


def info_num(info: str, key: str, cast):
    for f in info.split(";"):
        if f.startswith(key + "="):
            try:
                return cast(f[len(key) + 1:])
            except ValueError:
                return None
    return None


def job(task):
    import pysam

    ch, places = task
    out, tb = [], None
    for pos in places:
        rows = None
        for attempt in range(6):
            try:
                if tb is None:
                    tb = pysam.TabixFile(URL.format(ch), index=URL.format(ch) + ".tbi")
                rows = list(tb.fetch(ch, pos, pos + 1))
                break
            except (OSError, ValueError):
                tb = None
                time.sleep(2 + 3 * attempt)
        if rows is None:
            raise RuntimeError(f"{ch}:{pos} could not be read")
        recs = []
        for line in rows:
            f = line.split("\t", 8)
            if int(f[1]) - 1 == pos:
                info = f[7].split("\t")[0]
                recs.append((f[3], f[4], f[6], info_num(info, "AC", int), info_num(info, "AF", float), info_num(info, "AN", int)))
        out.append((ch, pos, recs))
    return out


def read_records(path: str, places: list) -> dict:
    if not os.path.exists(path):
        by: dict = {}
        for ch, p in sorted(set(places)):
            by.setdefault(ch, []).append(p)
        tasks = [(ch, ps[i:i + 60]) for ch, ps in by.items() for i in range(0, len(ps), 60)]
        with gzip.open(path + ".part", "wt") as w, Pool(8) as pool:
            for batch in pool.imap_unordered(job, tasks):
                for ch, pos, recs in batch:
                    w.write("\t".join([ch, str(pos)] + ["|".join(str(v) for v in r) for r in recs]) + "\n")
        os.replace(path + ".part", path)
    out = {}
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            recs = []
            for x in f[2:]:
                ref, alt, filt, ac, af, an = x.split("|")
                recs.append((ref, alt, filt, None if ac == "None" else int(ac), None if af == "None" else float(af), None if an == "None" else int(an)))
            out[(f[0], int(f[1]))] = recs
    return out


def status(recs: list, alt: str) -> str:
    """The new allele in gnomAD: absent, or FILTER (PASS or not) and frequency (one genome; under 0.1 %; more)."""
    same = [r for r in recs if r[1] == alt and len(r[0]) == 1]
    if not same:
        return "absent"
    _, _, filt, ac, af, _ = same[0]
    freq = "one genome" if ac == 1 else "under 0.1 %" if (af or 0) < 0.001 else "0.1 % or more"
    return f"{'PASS' if filt == 'PASS' else 'not PASS'}, {freq}"


def test(sel: list[dict], obs: str, exp: str) -> dict:
    sel = [x for x in sel if x.get(exp) not in (None, "")]
    if not sel:
        return {"mutations": 0}
    o = sum(x[obs] == "True" for x in sel)
    e = sum(float(x[exp]) for x in sel)
    v = sum(float(x[exp]) * (1 - float(x[exp])) for x in sel)
    return {"mutations": len(sel), "observed": o, "expected": round(e, 1), "ratio": round(o / e, 2) if e else None,
            "z": round((o - e) / math.sqrt(v), 1) if v else None}


def main(dnm_dir: str, chain37: str, records: str) -> None:
    rows = []
    with gzip.open(os.path.join(nm.RESULTS, "new_mutations.tsv.gz"), "rt") as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            rows.append(dict(zip(head, line.rstrip("\n").split("\t"))))
    snvs, _ = nm.sasani(dnm_dir, chain37)
    inside = {(r["chrom"], int(r["pos"])) for r in rows}
    others = [d for d in snvs if (d["chrom"], d["pos"]) not in inside]
    random.Random(SEED).shuffle(others)
    others = others[:OTHERS]
    recs = read_records(records, sorted(inside | {(d["chrom"], d["pos"]) for d in others}))
    base = Counter(status(recs[(d["chrom"], d["pos"])], d["alt"]) for d in others)
    for r in rows:
        r["gnomad"] = status(recs[(r["chrom"], int(r["pos"]))], r["alt"])
    once = [r for r in rows if r["also_in_sasani"] != "True"]
    out = {"outside_loci_sample": {"mutations": len(others), "gnomad": dict(base)}, "in_loci": {}}
    for pname, keep in (("purity 0.6+", lambda p: p >= 0.6), ("purity under 0.6", lambda p: p < 0.6)):
        sel = [r for r in once if keep(float(r["purity"]))]
        res = {"gnomad_of_copies": dict(Counter(r["gnomad"] for r in sel if r["copies"] == "True")),
               "gnomad_of_others": dict(Counter(r["gnomad"] for r in sel if r["copies"] != "True"))}
        for sname, sk in (("all", lambda g: True), ("absent", lambda g: g == "absent"), ("present", lambda g: g != "absent"),
                          ("present, PASS", lambda g: g.startswith("PASS")), ("present, not PASS", lambda g: g.startswith("not PASS")),
                          ("absent or PASS in one genome", lambda g: g in ("absent", "PASS, one genome"))):
            s2 = [r for r in sel if sk(r["gnomad"])]
            phased = [dict(r, is_toward=str(r["kind"] == "toward")) for r in s2 if r["kind"]]
            res[sname] = {"copies": test(s2, "copies", "expected_copies"), "copies_5mer": test(s2, "copies", "expected_copies_5mer"),
                          "near": test(s2, "near", "expected_near"), "toward": test(phased, "is_toward", "expected_toward")}
        out["in_loci"][pname] = res
    with open(os.path.join(nm.RESULTS, "in_gnomad.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    with gzip.open(os.path.join(nm.RESULTS, "in_gnomad.tsv.gz"), "wt") as fh:
        fh.write("chrom\tpos\talt\tgnomad\n")
        for r in rows:
            fh.write(f"{r['chrom']}\t{r['pos']}\t{r['alt']}\t{r['gnomad']}\n")
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit(__doc__)
    main(*sys.argv[1:4])
