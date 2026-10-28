"""03-compare: the same indel written at different places in a repeat. Does gpxil compare see it as the same?

Two VCF files are made on a reference of repeats: homopolymers and repeats of
2, 3 and 4 bases, from 8 to 160 bases long. In each repeat file A writes an
insertion or deletion of one or two units at one place (the left end, as
`bcftools norm` would, or the middle) and file B the same change at another
(the middle, or the right end); either heterozygous or homozygous. Some
repeats get a heterozygous substitution nearby too, phased in A on the other
haplotype and unphased in B. Controls write a different change in B (two
units for one; a deletion for an insertion; an inserted unit with a base
changed), which must not count as the same.

Both are encoded with gpxiltools and compared with `gpxil compare` at two gaps
(the default 50, and 200); every stretch where they differ base by base is
collected and matched to its repeat.

    python3 repeats.py [OUT_DIR]
"""

from __future__ import annotations

import json
import random
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

from gpxiltools.compare import compare_pair
from gpxiltools.container import GpxilFile
from gpxiltools.encode import encode

UNITS = ["A", "CA", "TTG", "AGAT"]
LENGTHS = [
    8,
    16,
    30,
    48,
    60,
    100,
    160,
]  #: bases of repeat (rounded down to whole units)
SPACER = 700  #: random bases between repeats: more than twice the larger gap tried
GAPS = (50, 200)


def flank(unit: str, rng: random.Random) -> str:
    return rng.choice([b for b in "ACGT" if b not in unit] or ["C"])


def build(rng: random.Random):
    """The reference, and the cases: one repeat each, with A's and B's records."""
    seq: list[str] = []
    cases = []

    def spacer():
        # random bases, each unlike the one before: no runs to be read as part of a repeat
        out = [rng.choice("ACGT")]
        for _ in range(SPACER - 1):
            out.append(rng.choice([b for b in "ACGT" if b != out[-1]]))
        return "".join(out)

    for unit in UNITS:
        p = len(unit)
        for length in LENGTHS:
            n = length // p
            for kind in ("ins", "del"):
                for k in (1, 2):
                    if kind == "del" and n < k + 2:
                        continue
                    for a_at, b_at in (
                        ("left", "right"),
                        ("left", "middle"),
                        ("middle", "right"),
                    ):
                        for gt in ("0/1", "1/1"):
                            for extra in ("", "snv-phased") if gt == "0/1" and k == 1 else ("",):
                                cases.append(
                                    dict(
                                        unit=unit,
                                        n=n,
                                        kind=kind,
                                        k=k,
                                        a=a_at,
                                        b=b_at,
                                        gt=gt,
                                        extra=extra,
                                        control=None,
                                    )
                                )
                    if k == 1:
                        for control in (
                            "two-units",
                            "del-for-ins" if kind == "ins" else None,
                            "changed-unit" if kind == "ins" else None,
                        ):
                            if control:
                                cases.append(
                                    dict(
                                        unit=unit,
                                        n=n,
                                        kind=kind,
                                        k=k,
                                        a="left",
                                        b="right",
                                        gt="0/1",
                                        extra="",
                                        control=control,
                                    )
                                )
    for c in cases:
        unit, n, p = c["unit"], c["n"], len(c["unit"])
        seq.append(spacer())
        left = flank(unit, rng)
        start = sum(map(len, seq))  # (the left flank's position, 0-based)
        seq.append(left + unit * n + flank(unit, rng))
        rpt = start + 1  # first base of the repeat
        c["start"], c["end"] = start - 20, rpt + n * p + 1
        c["ref_at"] = (start, left + unit * n)
    seq.append(spacer())
    ref = "".join(seq)
    for c in cases:
        c["records_a"], c["records_b"] = records(c, ref)
    return ref, cases


def placement(c, where: str, ref: str, kind=None, k=None, inserted=None):
    """One record: an insertion or deletion of k units at the left end, in the middle or at the right end."""
    unit, n, p = c["unit"], c["n"], len(c["unit"])
    kind, k = kind or c["kind"], k or c["k"]
    rpt = c["ref_at"][0] + 1
    j = {"left": 0, "middle": n // 2, "right": n if kind == "ins" else n - k}[where]
    anchor = rpt + j * p - 1  # (the base before the change)
    if kind == "ins":
        return (anchor, ref[anchor], ref[anchor] + (inserted or unit * k))
    return (anchor, ref[anchor : anchor + 1 + k * p], ref[anchor])


def records(c, ref: str):
    a = [placement(c, c["a"], ref)]
    if c["control"] == "two-units":
        b = [placement(c, c["b"], ref, k=2)]
    elif c["control"] == "del-for-ins":
        b = [placement(c, c["b"], ref, kind="del", k=1)]
    elif c["control"] == "changed-unit":
        u = c["unit"]
        changed = ("C" if u[0] != "C" else "G") + u[1:]
        b = [placement(c, c["b"], ref, inserted=changed)]
    else:
        b = [placement(c, c["b"], ref)]
    ga = gb = c["gt"]
    out_a = [(*r, ga) for r in a]
    out_b = [(*r, gb) for r in b]
    if c["extra"] == "snv-phased":  # a substitution 15 bases before the repeat: in A on the other haplotype, phased
        pos = c["ref_at"][0] - 15
        alt = next(x for x in "ACGT" if x != ref[pos])
        out_a = [
            (pos, ref[pos], alt, "0|1", pos + 1),
            *[(*r[:3], "1|0", pos + 1) for r in a],
        ]
        out_b = [(pos, ref[pos], alt, "0/1"), *out_b]
    return out_a, out_b


def write(tmp: Path, ref: str, cases, side: str) -> Path:
    vcf = tmp / f"{side}.vcf"
    head = [
        "##fileformat=VCFv4.3",
        f"##contig=<ID=chr1,length={len(ref)}>",
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">',
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{side}",
    ]
    body = []
    for c in cases:
        for pos0, r, a, gt, *ps in c[f"records_{side}"]:
            fmt, val = ("GT:PS", f"{gt}:{ps[0]}") if ps else ("GT", gt)
            body.append((pos0, f"chr1\t{pos0 + 1}\t.\t{r}\t{a}\t50\tPASS\t.\t{fmt}\t{val}"))
    vcf.write_text("\n".join(head + [x for _, x in sorted(body)]) + "\n")
    out = tmp / f"{side}.gpxil.tif"
    encode(tmp / "ref.fa", out, vcf=vcf, log=lambda m: None)
    return out


def jump(c) -> int:
    """How far apart the bases are where A and B differ, at most, from one to the next (0: written alike)."""

    def bases(rec):
        pos, ref, alt, *_ = rec
        return {pos} if len(alt) > len(ref) else set(range(pos + 1, pos + len(ref)))  # (an insertion's anchor; the deleted bases)

    diff = sorted(bases(c["records_a"][-1]) ^ bases(c["records_b"][-1]))
    return max((y - x for x, y in zip(diff, diff[1:])), default=0)


def main(out_dir: Path) -> None:
    rng = random.Random(3)
    ref, cases = build(rng)
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp / "ref.fa").write_text(">chr1\n" + "\n".join(ref[i : i + 60] for i in range(0, len(ref), 60)) + "\n")
        fa, fb = write(tmp, ref, cases, "a"), write(tmp, ref, cases, "b")
        found = {}
        with GpxilFile(fa) as a, GpxilFile(fb) as b:
            for gap in GAPS:
                st: list = []
                compare_pair(a, b, gap=gap, stretches=st)
                found[gap] = st
    by_case = {gap: defaultdict(list) for gap in GAPS}
    for gap, st in found.items():
        for _, s1, _e1, alike in st:
            hit = [i for i, c in enumerate(cases) if c["start"] <= s1 - 1 < c["end"]]
            assert len(hit) == 1, (s1, hit)
            by_case[gap][hit[0]].append(alike)
    rows = []
    for i, c in enumerate(cases):
        r = {k: c[k] for k in ("unit", "n", "kind", "k", "a", "b", "gt", "extra", "control")}
        r["repeat_bases"] = c["n"] * len(c["unit"])
        r["jump"] = jump(c)
        for gap in GAPS:
            got = by_case[gap].get(i, [])
            r[f"gap{gap}"] = "same base by base" if not got else "the same" if all(got) else "different"
        rows.append(r)
    (out_dir / "repeats.json").write_text(json.dumps(rows, indent=1) + "\n")
    report(rows, out_dir)


BINS = ((1, 10), (11, 50), (51, 100), (101, 200))


def report(rows, out_dir: Path) -> None:
    lines = ["# Indels written at different places in a repeat", ""]
    alike = [r for r in rows if not r["control"] and r["jump"] == 0]  # (the two places are one: a short repeat)
    real = [r for r in rows if not r["control"] and r["jump"] > 0]
    ctrl = [r for r in rows if r["control"]]
    assert all(r["gap50"] == r["gap200"] == "same base by base" for r in alike)
    lines += [
        f"{len(real)} repeats with the same change written at two places, and {len(ctrl)} controls with a different change "
        f"({len(alike)} more, in repeats too short for two places, are written alike).",
        "",
        "| bases between the two places | repeats | " + " | ".join(f"the same, gap {g}" for g in GAPS) + " |",
        "|---|---:|" + "---:|" * len(GAPS),
    ]
    for lo, hi in BINS:
        got = [r for r in real if lo <= r["jump"] <= hi]
        cells = [str(sum(r[f"gap{g}"] == "the same" for r in got)) for g in GAPS]
        lines.append(f"| {lo}-{hi} | {len(got)} | " + " | ".join(cells) + " |")
    assert sum(lo <= r["jump"] <= hi for r in real for lo, hi in BINS) == len(real)
    lines += ["| controls (a different change) | " + str(len(ctrl)) + " | "]
    lines[-1] += " | ".join(str(sum(r[f"gap{g}"] == "the same" for r in ctrl)) for g in GAPS) + " |"
    missed = [r for r in real if any(r[f"gap{g}"] != "the same" for g in GAPS)]
    if missed:
        lines += [
            "",
            "Not seen as the same:",
            "",
            "```",
            *[json.dumps(r) for r in missed[:20]],
            "```",
        ]
    by = Counter((r["unit"], r["extra"], r["gap50"]) for r in real)
    lines += [
        "",
        "By repeat unit, at gap 50:",
        "",
        "| unit | a phased substitution nearby | the same | different |",
        "|---|---|---:|---:|",
    ]
    for unit in UNITS:
        for extra in ("", "snv-phased"):
            lines.append(f"| {unit} | {'yes' if extra else 'no'} | {by[(unit, extra, 'the same')]} | {by[(unit, extra, 'different')]} |")
    (out_dir / "repeats.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent / "results"
    out.mkdir(parents=True, exist_ok=True)
    main(out)
