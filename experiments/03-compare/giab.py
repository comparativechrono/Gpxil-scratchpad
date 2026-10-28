"""03-compare: gpxil compare on GIAB's HG002, and its parents HG003 and HG004, checked against the VCF files.

Two checks, each made a second way, from the VCF records with pysam, without
gpxiltools:

1. **Two versions of HG002** (v4.2.1 and v5.0q, chromosome 17). Every stretch
   where `gpxil compare` finds the two differ base by base is spelled out
   again from each VCF's records: each file's two haplotype sequences over
   the stretch, its heterozygous records placed either way round. The stretch
   is the same, written differently, if the two files can spell the same
   pair of sequences.
2. **The trio** (v4.2.1, chromosome 17). Every Mendelian error `gpxil compare`
   finds is looked up in the three VCF and BED files; and the records are
   checked for Mendelian errors one by one, where all three files have their
   regions called with confidence.

    python3 giab.py DATA GPXIL_DIR [OUT_DIR]

DATA holds the GIAB files (giab/, giab5/, giab-trio/; see the README) and
GPXIL_DIR the gpxil files of chromosome 17: HG002_chr17.gpxil.tif (v4.2.1),
HG002_chr17_v5.0q.noref.gpxil.tif, HG003_chr17.gpxil.tif and HG004_chr17.gpxil.tif.
"""

from __future__ import annotations

import bisect
import itertools
import json
import sys
import time
from collections import Counter
from pathlib import Path

import pysam
from gpxiltools.compare import compare_pair, compare_trio, format_compare
from gpxiltools.container import GpxilFile

CONTIG = "chr17"


def bed(path: Path) -> list[tuple[int, int]]:
    out = []
    for line in open(path):
        f = line.split()
        if f and f[0] == CONTIG:
            out.append((int(f[1]), int(f[2])))
    return sorted(out)


def inside(iv, a: int, b: int | None = None) -> bool:
    """Whether [a, b) (or the base a) lies in one interval."""
    b = a + 1 if b is None else b
    i = bisect.bisect_right(iv, (a, float("inf"))) - 1
    return i >= 0 and iv[i][0] <= a and b <= iv[i][1]


def overlaps(iv, a: int, b: int) -> bool:
    i = bisect.bisect_right(iv, (b, float("inf")))
    return any(x < b and a < y for x, y in iv[max(0, i - 3) : i])


class Calls:
    """One VCF's records on chromosome 17, and its regions called with confidence."""

    def __init__(self, vcf: Path, bed_path: Path):
        self.vcf = pysam.VariantFile(str(vcf))
        self.sample = list(self.vcf.header.samples)[0]
        self.bed = bed(bed_path)

    def records(self, a: int, b: int):
        """Records overlapping [a, b), as (start, end, alleles, genotype, phased); none outside the confident regions."""
        out = []
        for rec in self.vcf.fetch(CONTIG, max(0, a - 1), b):
            start, end = rec.pos - 1, rec.pos - 1 + len(rec.ref)
            if not (start < b and a < max(end, start + 1)) or not overlaps(self.bed, start, max(end, start + 1)):
                continue
            s = rec.samples[self.sample]
            out.append((start, end, (rec.ref, *rec.alts), tuple(s["GT"]), s.phased))
        return out


def trimmed(start: int, ref: str, alt: str) -> tuple[int, int, str]:
    """An allele with the bases it shares with REF at its ends taken off: (start, end, the bases put there)."""
    while ref and alt and ref[-1] == alt[-1]:
        ref, alt = ref[:-1], alt[:-1]
    k = 0
    while k < min(len(ref), len(alt)) and ref[k] == alt[k]:
        k += 1
    return start + k, start + len(ref), alt[k:]


def spell(ref: str, r0: int, recs, hap: list[int]) -> str | None:
    """A haplotype's sequence over [r0, r0 + len(ref)) with its alleles of the records (None: two of them overlap)."""
    out, at = [], r0
    changes = []
    for (start, _, alleles, _, _), allele in zip(recs, hap):
        if allele and alleles[allele] != "*":  # (*: deleted by another record)
            changes.append(trimmed(start, alleles[0], alleles[allele]))
    for start, end, alt in sorted(changes):
        if start < at:
            return None
        out += [ref[at - r0 : start - r0], alt]
        at = end
    out.append(ref[at - r0 :])
    return "".join(out)


def pairs(ref: str, r0: int, recs) -> set:
    """Every pair of sequences the records allow, heterozygous ones placed either way round."""
    if any(None in gt for *_, gt, _ in recs):
        return {("?",)}
    choices = [[gt] if gt[0] == gt[-1] else [gt, gt[::-1]] for *_, gt, _ in recs]
    out = set()
    for pick in itertools.islice(itertools.product(*choices), 1 << 14):
        h1 = spell(ref, r0, recs, [g[0] for g in pick])
        h2 = spell(ref, r0, recs, [g[-1] for g in pick])
        if h1 is not None and h2 is not None:
            out.add(tuple(sorted((h1, h2))))
    return out


def check_pair(gfiles, calls, outd: Path, refc) -> dict:
    stretches: list = []
    t0 = time.time()
    with GpxilFile(gfiles["v4.2.1"]) as a, GpxilFile(gfiles["v5.0q"]) as b:
        res = compare_pair(a, b, stretches=stretches)
    seconds = time.time() - t0
    (outd / "pair.txt").write_text(format_compare(res) + "\n")
    (outd / "pair.json").write_text(json.dumps(res, indent=1) + "\n")
    rows, agree = [], Counter()
    for _, s1, e1, alike in stretches:
        a, b = s1 - 1, e1
        while True:  # (the records of both files over the stretch, and the stretch widened to take them in whole)
            recs = [c.records(a, b) for c in calls]
            lo = min([a] + [r[0] for rs in recs for r in rs])
            hi = max([b] + [max(r[1], r[0] + 1) for rs in recs for r in rs])
            if (lo, hi) == (a, b):
                break
            a, b = lo, hi
        ref = refc.sequence(a, b).upper()
        same = bool(pairs(ref, a, recs[0]) & pairs(ref, a, recs[1]))
        agree[(alike, same)] += 1
        rows.append(
            {
                "start": s1,
                "end": e1,
                "gpxil": "the same" if alike else "different",
                "vcf": "the same" if same else "different",
                "records": [[describe(r) for r in rs] for rs in recs],
            }
        )
    (outd / "pair-stretches.json").write_text(json.dumps(rows, indent=1) + "\n")
    return {
        "seconds": round(seconds, 1),
        "stretches": len(rows),
        "agreement": {f"gpxil {'same' if x else 'different'}, vcf {'same' if y else 'different'}": n for (x, y), n in agree.items()},
    }


def describe(r) -> str:
    start, _, alleles, gt, phased = r
    return f"{start + 1} {alleles[0]}>{','.join(alleles[1:])} {('|' if phased else '/').join(map(str, gt))}"


def genotype_at(calls: Calls, a: int, b: int):
    """A file's alleles over [a, b): from the one record that covers it exactly, REF/REF without one; None otherwise."""
    recs = calls.records(a, b)
    if not recs:
        return ("REF", "REF")
    if len(recs) == 1 and recs[0][0] == a and recs[0][1] == b:
        _, _, alleles, gt, _ = recs[0]
        return tuple(alleles[g] if g else "REF" for g in gt)
    return None


def check_trio(gfiles, trio: dict, outd: Path) -> dict:
    errors: list = []
    t0 = time.time()
    with (
        GpxilFile(gfiles["HG002"]) as c,
        GpxilFile(gfiles["HG003"]) as f,
        GpxilFile(gfiles["HG004"]) as m,
    ):
        res = compare_trio(c, f, m, errors=errors)
    seconds = time.time() - t0
    (outd / "trio.txt").write_text(format_compare(res) + "\n")
    (outd / "trio.json").write_text(json.dumps(res, indent=1) + "\n")
    # each error, looked up in the VCF and BED files
    looked = []
    for e in errors:
        p0 = e["position"] - 1
        row = {**e}
        for who, calls in trio.items():
            row[f"{who}_bed"] = inside(calls.bed, p0)
            row[f"{who}_vcf"] = [describe(r) for r in calls.records(p0, p0 + 1)]
        looked.append(row)
    (outd / "trio-errors.json").write_text(json.dumps(looked, indent=1) + "\n")
    # Mendelian errors record by record, where all three are called with confidence
    beds = [c.bed for c in trio.values()]
    sites = set()
    for calls in trio.values():
        for rec in calls.vcf.fetch(CONTIG):
            a, b = rec.pos - 1, rec.pos - 1 + len(rec.ref)
            if all(inside(iv, a, b) for iv in beds):
                sites.add((a, b))
    count = Counter()
    inconsistent = []
    for a, b in sorted(sites):
        child, father, mother = (genotype_at(c, a, b) for c in trio.values())
        if None in (child, father, mother) or len(child) != 2:
            count["not one record each (overlapping or written differently)"] += 1
            continue
        ok = (child[0] in father and child[1] in mother) or (child[1] in father and child[0] in mother)
        count["consistent" if ok else "inconsistent"] += 1
        if not ok:
            inconsistent.append({"position": a + 1, "child": child, "father": father, "mother": mother})
    # (a record's error is at its POS, gpxil's at the first base that differs, and a run of bases is one error)
    found = {e["position"] for e in errors}
    near = {p for e in errors for p in range(e["position"] - 1, e["position"] + 2)}
    by_records = {x["position"] for x in inconsistent}
    matched = {x for x in by_records if x in near}
    return {
        "seconds": round(seconds, 1),
        "gpxil_errors": len(errors),
        "gpxil_new": sum(e["new"] for e in errors),
        "errors_inside_all_beds": sum(all(r[f"{w}_bed"] for w in trio) for r in looked),
        "records": dict(count),
        "record_errors": len(inconsistent),
        "both_ways": len({e for e in found if any(abs(e - x) <= 1 for x in by_records)}),
        "only_gpxil": [e for e in looked if not any(abs(e["position"] - x) <= 1 for x in by_records)],
        "only_records": [x for x in inconsistent if x["position"] not in matched],
    }


def main(data: Path, gdir: Path, outd: Path) -> None:
    outd.mkdir(parents=True, exist_ok=True)
    gfiles = {
        "v4.2.1": gdir / "HG002_chr17.gpxil.tif",
        "v5.0q": gdir / "HG002_chr17_v5.0q.noref.gpxil.tif",
        "HG002": gdir / "HG002_chr17.gpxil.tif",
        "HG003": gdir / "HG003_chr17.gpxil.tif",
        "HG004": gdir / "HG004_chr17.gpxil.tif",
    }
    v421 = Calls(
        data / "giab/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        data / "giab/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
    )
    v5 = Calls(
        data / "giab5/HG002_GRCh38_v5.0q_smvar.vcf.gz",
        data / "giab5/HG002_GRCh38_v5.0q_smvar.benchmark.bed",
    )
    trio = {
        "child": v421,
        "father": Calls(
            data / "giab-trio/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
            data / "giab-trio/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        ),
        "mother": Calls(
            data / "giab-trio/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
            data / "giab-trio/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed",
        ),
    }
    with GpxilFile(gfiles["v4.2.1"]) as g:
        summary = {"pair": check_pair(gfiles, [v421, v5], outd, g[CONTIG])}
    summary["trio"] = check_trio(gfiles, trio, outd)
    (outd / "giab.json").write_text(json.dumps(summary, indent=1) + "\n")
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    main(
        Path(sys.argv[1]),
        Path(sys.argv[2]),
        Path(sys.argv[3]) if len(sys.argv) > 3 else Path(__file__).parent / "results",
    )
