"""What a gpxil keeps of a VCF, and what comes back out of it.

GIAB's HG002 v5.0q small-variant benchmark is encoded on GRCh38 (the whole genome, sample only: the reference is the
GRCh38 gpxil's) and checked, chromosome by chromosome, three ways:

- haplotypes: each haplotype as gpxiltools reads it from the gpxil against the one bcftools consensus builds from the
  VCF, a separate implementation, compared base for base wherever the gpxil calls the base (UNKNOWN, read as N, is
  passed over). Where the two put an insertion or deletion in different places of a repeat, the stretch around it is
  compared as a whole. Where they differ, the difference is set against the records bcftools consensus left out
  because they overlap an allele it had already put on that haplotype;
- records: the VCF the gpxil exports against the VCF it came from, both normalised by bcftools norm (left-aligned,
  multiallelic records split, complex ones atomised), as (position, REF, ALT, phased genotype). '*' alleles, which
  only point at a deletion written elsewhere, and genotypes without an ALT allele are left out. A record of one that
  isn't in the other is set against the haplotype comparison: where the haplotypes agree, the two files say the same
  thing in different records;
- planes: gpxil -> VCF (with its no-call blocks) -> gpxil, and the planes compared.

    python3 roundtrip.py encode REF_GPXIL VCF BED WORK     -> results/encode.json
    python3 roundtrip.py check  REF_GPXIL VCF WORK [CHROMS] -> results/check.json
    python3 roundtrip.py summary                            -> results/check_summary.json

REF_GPXIL: GRCh38.gpxil.tif (DATA.md). VCF, BED: HG002_GRCh38_v5.0q_smvar.vcf.gz and its benchmark regions. WORK: a
folder for the big files: the sample's gpxil, and each chromosome's FASTA, haplotypes and VCFs while it is checked.
CHROMS: comma-separated (default chr1..chr22, chrX, chrY). check keeps what it has done in results/check.json, and
goes on from there when run again.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
from collections import Counter

import gpxiltools
import numpy as np
import pysam
from gpxiltools import codes as C
from gpxiltools.container import GpxilFile
from gpxiltools.export import iter_haplotype
from pysam import bcftools

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
SAMPLE = "HG002_v5.0q.gpxil.tif"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
N, STAR = ord("N"), ord("*")
REJOIN = (10_000, 50)  #: how far (bases, events) a stretch out of step is followed for the two to come back into step

# A command run in a process of its own, so that its peak memory is its own (RUSAGE_CHILDREN of the wrapper).
WRAP = (
    "import json, resource, subprocess, sys, time\n"
    "t = time.perf_counter()\n"
    "with open(sys.argv[2], 'w') as log:\n"
    "    r = subprocess.run(sys.argv[3:], stdout=log, stderr=subprocess.STDOUT)\n"
    "w = time.perf_counter() - t\n"
    "u = resource.getrusage(resource.RUSAGE_CHILDREN)\n"
    "with open(sys.argv[1], 'w') as out:\n"
    "    json.dump({'returncode': r.returncode, 'wall_s': round(w, 1), 'cpu_s': round(u.ru_utime + u.ru_stime, 1),\n"
    "               'peak_rss_mb': round(u.ru_maxrss / 1024, 1)}, out)\n"
)


def measure(cmd: list[str], log: str) -> dict:
    """Run ``cmd``, its output to ``log``; its wall and CPU time and peak memory."""
    fd, tmp = tempfile.mkstemp(suffix=".json")
    os.close(fd)
    subprocess.run([sys.executable, "-I", "-c", WRAP, tmp, log, *cmd], check=True)
    with open(tmp) as fh:
        out = json.load(fh)
    os.unlink(tmp)
    if out["returncode"] != 0:
        raise SystemExit(f"{' '.join(cmd[:2])} failed (see {log})")
    return out


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 24), b""):
            h.update(block)
    return h.hexdigest()


def save(name: str, obj) -> None:
    os.makedirs(RESULTS, exist_ok=True)
    tmp = os.path.join(RESULTS, name + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1)
    os.replace(tmp, os.path.join(RESULTS, name))


# ----------------------------------------------------------------------------------------------------------- encode


def step_encode(ref: str, vcf: str, bed: str, work: str) -> None:
    os.makedirs(work, exist_ok=True)
    out = os.path.join(work, SAMPLE)
    cmd = ["gpxil", "encode", ref, "-v", vcf, "--callable", bed, "--sex", "XY", "--no-ref", "-o", out]
    enc = measure(cmd, os.path.join(RESULTS, "encode.log"))
    val = measure(["gpxil", "validate", out], os.path.join(RESULTS, "validate.log"))
    with open(os.path.join(RESULTS, "validate.log")) as fh:
        val["says"] = fh.read().strip().splitlines()[-3:]
    totals: Counter = Counter()
    contigs = []
    with GpxilFile(out) as g:
        for e in g.metadata["contigs"]:
            contigs.append(e["name"])
            for k, v in (e.get("encoding") or {}).items():
                if isinstance(v, int) and not isinstance(v, bool):
                    totals[k] += v
            totals["callable_bases"] += e.get("callable_bases") or 0
    res = {
        "gpxiltools": gpxiltools.__version__,
        "command": " ".join(
            [
                "gpxil",
                "encode",
                "GRCh38.gpxil.tif",
                "-v",
                os.path.basename(vcf),
                "--callable",
                os.path.basename(bed),
                "--sex",
                "XY",
                "--no-ref",
                "-o",
                SAMPLE,
            ]
        ),
        "encode": enc,
        "validate": val,
        "file_bytes": os.path.getsize(out),
        "sha256": sha256(out),
        "contigs": contigs,
        "stats": dict(totals),
        "input_vcf_bytes": os.path.getsize(vcf),
    }
    save("encode.json", res)
    print(json.dumps(res, indent=1))


# ----------------------------------------------------------------------------------------------------- haplotypes


def runs(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    d = np.diff(mask.view(np.int8), prepend=np.int8(0), append=np.int8(0))
    return np.flatnonzero(d == 1), np.flatnonzero(d == -1)


def gpxil_events(states: np.ndarray, ins: dict[int, str]) -> list[tuple[int, int, int, str]]:
    """One haplotype's departures in length from the reference: runs of DEL, runs of ABSENT (a haplotype that isn't
    there, which isn't compared) and insertions, as (start, end, change in length, kind) on the reference."""
    ev = []
    for kind, code in (("del", C.HAP_DEL), ("absent", C.HAP_ABSENT)):
        a, b = runs(states == code)
        ev += [(int(s), int(e), -int(e - s), kind) for s, e in zip(a, b)]
    ev += [(a + 1, a + 1, len(s), "ins") for a, s in ins.items()]
    return ev


def chain_events(path: str) -> list[tuple[int, int, int, str]]:
    """bcftools consensus's chain (the reference is the target): its gaps, as (start, end, change in length, 'gap')."""
    ev = []
    t = None
    with open(path) as fh:
        for line in fh:
            f = line.split()
            if not f:
                continue
            if f[0] == "chain":
                t = int(f[5])
                continue
            t += int(f[0])
            if len(f) == 3:
                dt, dq = int(f[1]), int(f[2])
                ev.append((t, t + dt, dq - dt, "gap"))
                t += dt
    return ev


def clusters_of(gev, cev, n: int) -> list[tuple[int, int, int, int, bool]]:
    """Both sides' events, a base added each way (the anchor and the base after) and merged where they meet or
    overlap: (start, end, gpxil's change in length, consensus's, whether a haplotype isn't there)."""
    ev = [(max(0, s - 1), min(n, e + 1), d, 0, k == "absent") for s, e, d, k in gev]
    ev += [(max(0, s - 1), min(n, e + 1), 0, d, False) for s, e, d, k in cev]
    ev.sort()
    out: list[list] = []
    for a, b, dg, dc, ab in ev:
        if out and a <= out[-1][1]:
            o = out[-1]
            o[1], o[2], o[3], o[4] = max(o[1], b), o[2] + dg, o[3] + dc, o[4] or ab
        else:
            out.append([a, b, dg, dc, ab])
    return [tuple(o) for o in out]


def compare_haplotype(gx: bytes, cs: bytes, gev, cev, n: int) -> dict:
    """gpxil's haplotype ``gx`` against consensus's ``cs`` (both over a reference of ``n`` bases).

    Between events the two are read in step with the reference and compared base for base. An event's stretch (with
    the base either side) is compared as a whole; where the two sides change the length by different amounts there
    (one put a deletion at one end of a repeat, the other at the other end), the stretch takes in the next events until
    they agree again, within ``REJOIN``. N in either (nothing known there: UNKNOWN in the gpxil's, a base the VCF's
    reference masks in consensus's) matches anything."""
    G = np.frombuffer(gx, dtype=np.uint8)
    S = np.frombuffer(cs, dtype=np.uint8)
    cl = clusters_of(gev, cev, n)
    st = Counter()
    sites: list[tuple[int, int, bool]] = []  # [a, b) where the two differ, and whether consensus wrote a '*' there
    offg = offc = 0
    pos = 0

    def bulk(a: int, b: int) -> None:
        if b <= a:
            return
        g, c = G[a + offg : b + offg], S[a + offc : b + offc]
        assert g.shape[0] == c.shape[0] == b - a, (a, b, offg, offc)
        unk = g == N
        cun = (c == N) & ~unk
        nu, ncu = int(unk.sum()), int(cun.sum())
        st["unknown"] += nu
        st["consensus_unknown"] += ncu
        st["called"] += (b - a) - nu - ncu
        bad = np.flatnonzero((g != c) & ~unk & ~cun)
        if bad.size:
            st["differing_bases"] += int(bad.size)
            # neighbouring bases that differ are one site
            cut = np.flatnonzero(np.diff(bad) > 1) + 1
            for part in np.split(bad, cut):
                lo, hi = int(part[0]), int(part[-1]) + 1
                sites.append((a + lo, a + hi, bool(np.any(c[lo:hi] == STAR))))

    lim_b, lim_n = REJOIN
    i = 0
    while i < len(cl):
        A, B, dg, dc, ab = cl[i]
        bulk(pos, A)
        j = i
        if not ab and dg != dc:
            J, B2, dg2, dc2 = i, B, dg, dc
            while J + 1 < len(cl) and J + 1 - i <= lim_n and cl[J + 1][1] - A <= lim_b and not cl[J + 1][4]:
                J += 1
                B2, dg2, dc2 = cl[J][1], dg2 + cl[J][2], dc2 + cl[J][3]
                if dg2 == dc2:
                    break
            if dg2 == dc2:
                j, B, dg, dc = J, B2, dg2, dc2
                st["stretches_taken_together"] += 1
        if ab:
            st["absent_reference_bases"] += B - A
        else:
            g, c = G[A + offg : B + offg + dg], S[A + offc : B + offc + dc]
            unk = g == N
            st["unknown"] += int(unk.sum())
            st["called"] += int(g.shape[0] - unk.sum())
            st["event_stretches"] += 1
            if g.shape[0] != c.shape[0] or bool(np.any((g != c) & ~unk & (c != N))):
                st["differing_stretches"] += 1
                sites.append((A, B, bool(np.any(c == STAR))))
        offg += dg
        offc += dc
        pos = B
        i = j + 1
    bulk(pos, n)
    assert n + offg == G.shape[0] and n + offc == S.shape[0], "the events don't account for the lengths"
    st["gpxil_length"], st["consensus_length"] = int(G.shape[0]), int(S.shape[0])
    return {"stats": dict(st), "sites": sites}


def vcf_reference(vcf: str, chrom: str, fa: str, d: str) -> tuple[str, int]:
    """The VCF's own reference for this chromosome, for bcftools: ``fa`` with N wherever a record's REF allele has N
    and ``fa`` a base, and how many bases that masks. (GIAB's v5.0q has one such record, on chromosome 9: its
    reference masks 17,993 bases that Ensembl's has. gpxil reads N in a REF allele as any base, spec 6.7; bcftools
    consensus and norm refuse the record against the unmasked bases.)"""
    f = pysam.FastaFile(fa)
    masked = set()
    with pysam.VariantFile(vcf) as v:
        for r in v.fetch(chrom):
            ref = r.ref.upper()
            if "N" in ref:
                seq = f.fetch(chrom, r.pos - 1, r.pos - 1 + len(ref)).upper()
                masked.update(r.pos - 1 + i for i, (a, b) in enumerate(zip(ref, seq)) if a == "N" and b != "N")
    if not masked:
        f.close()
        return fa, 0
    seq = bytearray(f.fetch(chrom).encode())
    f.close()
    for p in masked:
        seq[p] = N
    out = os.path.join(d, "vcfref.fa")
    with open(out, "wb") as fh:
        fh.write(f">{chrom}\n".encode())
        fh.writelines(bytes(seq[i : i + 60]) + b"\n" for i in range(0, len(seq), 60))
    pysam.faidx(out)
    return out, len(masked)


def consensus(fa: str, vcf: str, h: int, out_fa: str, chain: str) -> dict:
    bcftools.consensus("-f", fa, "-H", str(h), "-c", chain, "-o", out_fa, vcf, catch_stdout=False)
    msg = bcftools.consensus.get_messages() or ""
    lines = [x for x in msg.splitlines() if x and not x.startswith("[W::hts_idx_load3]")]
    applied = [int(m.group(1)) for x in lines for m in [re.match(r"Applied (\d+) variants", x)] if m]
    overlaps = [
        int(m.group(1)) for x in lines for m in [re.match(r"The site \S+:(\d+) overlaps with another variant", x)] if m
    ]
    other = Counter(re.sub(r"\S+:\d+", "SITE", x)[:100] for x in lines if not x.startswith(("Applied", "The site")))
    return {"applied": applied[0] if applied else None, "overlaps": overlaps, "other_messages": dict(other)}


def read_fasta(path: str) -> bytes:
    with open(path, "rb") as fh:
        data = fh.read()
    return data.split(b"\n", 1)[1].replace(b"\n", b"").upper()


def skipped_alleles(vf: pysam.VariantFile, chrom: str, positions, h: int) -> list[tuple[int, int, str]]:
    """The records bcftools consensus left out on haplotype ``h`` at these positions (1-based), but for '*' alleles,
    which only point at a deletion already put there: [(start, end, ALT)] (0-based, end exclusive)."""
    out = []
    for pos in sorted(set(positions)):
        for r in vf.fetch(chrom, pos - 1, pos):
            if r.pos != pos:
                continue
            gt = r.samples[0].get("GT") or ()
            if len(gt) < h or gt[h - 1] in (None, 0):
                continue
            alt = r.alts[gt[h - 1] - 1]
            if alt != "*":
                out.append((r.pos - 1, r.pos - 1 + len(r.ref), alt))
    return out


class Spans:
    """Intervals [s, e) on the reference, asked whether any meets a stretch."""

    def __init__(self, spans):
        spans = sorted(spans)
        self.starts = [s for s, _ in spans]
        self.max_end = list(np.maximum.accumulate([e for _, e in spans])) if spans else []

    def meets(self, a: int, b: int) -> bool:
        """Whether any interval meets [a, b), a base either side counting: one that starts by b and ends after a - 1."""
        k = bisect.bisect_right(self.starts, b)
        return k > 0 and self.max_end[k - 1] >= a


# --------------------------------------------------------------------------------------------------------- records


def gt_text(s) -> str | None:
    gt = s.get("GT")
    if gt is None or all(a in (None, 0) for a in gt):
        return None  # (no ALT allele on either haplotype)
    return ("|" if s.phased else "/").join("." if a is None else str(a) for a in gt)


def normalised(src: str, fa: str, dst: str) -> None:
    bcftools.norm("-a", "--atom-overlaps", ".", "-m", "-any", "-f", fa, "-o", dst, src, catch_stdout=False)


def record_keys(path: str) -> tuple[Counter, dict]:
    keys: Counter = Counter()
    filt = {}
    for r in pysam.VariantFile(path):
        alt = r.alts[0] if r.alts else "."
        if alt == "*" or alt.startswith("<"):
            continue
        g = gt_text(r.samples[0])
        if g is None:
            continue
        k = (r.pos, r.ref, alt, g)
        keys[k] += 1
        filt[k] = list(r.filter) in ([], ["PASS"])
    return keys, filt


def compare_records(vcf: str, export: str, chrom: str, fa: str, vfa: str, d: str, differing: Spans) -> dict:
    src = os.path.join(d, "in.vcf")
    with pysam.VariantFile(vcf) as vin, pysam.VariantFile(src, "w", header=vin.header) as vout:
        n_in = 0
        for r in vin.fetch(chrom):
            vout.write(r)
            n_in += 1
    n_out = nocall = 0
    with pysam.VariantFile(export) as v:
        for r in v:
            n_out += 1
            nocall += bool(r.alts) and r.alts[0] == "<*>"
    a_n, b_n = os.path.join(d, "in.norm.vcf"), os.path.join(d, "out.norm.vcf")
    normalised(src, vfa, a_n)  # (the VCF against its own reference, the export against the gpxil's)
    normalised(export, fa, b_n)
    A, fa_ = record_keys(a_n)
    B, fb_ = record_keys(b_n)
    both = A & B
    kinds = Counter()
    shown = {}
    # A record of one file that the other hasn't is (a) where the haplotypes differ (above), or else the two say the
    # same haplotypes in other records: (b) the export writes an UNKNOWN base inside an allele as N (a genotype's '.'
    # is for a whole record), or (c) they are written differently (most often: the export writes the two haplotypes'
    # overlapping events as one record, which normalising doesn't split as the VCF had them).
    for side, only in (("only_in_vcf", A - B), ("only_in_export", B - A)):
        for k, m in only.items():
            pos, ref, alt, g = k
            at_difference = differing.meets(pos - 1, pos - 1 + len(ref))
            kind = (
                "at_a_haplotype_difference"
                if at_difference
                else ("unknown_bases_written_as_N" if "N" in alt else "written_differently")
            )
            kinds[(side, kind)] += m
            shown.setdefault(f"{side}/{kind}", [])
            if len(shown[f"{side}/{kind}"]) < 6:
                shown[f"{side}/{kind}"].append(f"{pos} {ref[:30]} {alt[:30]} {g}")
    return {
        "vcf_records": n_in,
        "export_records": n_out,
        "export_nocall_blocks": nocall,
        "vcf_alleles": sum(A.values()),
        "export_alleles": sum(B.values()),
        "same": sum(both.values()),
        "same_filter": sum(m for k, m in both.items() if fa_[k] == fb_[k]),
        "differing": {f"{s}/{k}": m for (s, k), m in sorted(kinds.items())},
        "examples": shown,
    }


# ----------------------------------------------------------------------------------------------------------- planes


def compare_planes(a_path: str, b_path: str, chrom: str, ref: str) -> dict:
    """The sample's gpxil against the one its export was encoded into: each haplotype's sequence (N for UNKNOWN), and
    the planes base for base. Where the export wrote the two haplotypes' overlapping events as one record, the encoder
    aligns each allele to the reference afresh (spec 6.2), which can say the same haplotype with other states."""
    out = {}
    with GpxilFile(a_path, reference=ref) as a, GpxilFile(b_path, reference=ref) as b:
        ca, cb = a[chrom], b[chrom]
        out["planes"] = [list(ca.planes), list(cb.planes)]
        out["haplotypes_same"] = {}
        for h in (1, 2):
            ha = b"".join(iter_haplotype(ca, h, keep_case=False))
            hb = b"".join(iter_haplotype(cb, h, keep_case=False))
            out["haplotypes_same"][str(h)] = ha == hb
            if ha != hb:
                k = next(i for i in range(min(len(ha), len(hb)) + 1) if i == min(len(ha), len(hb)) or ha[i] != hb[i])
                out.setdefault("haplotype_first_difference", {})[str(h)] = [k, len(ha), len(hb)]
            del ha, hb
        ga, gb = ca.read("gt"), cb.read("gt")
        bad = np.flatnonzero(ga != gb)
        out["gt_bases_differing"] = int(bad.size)
        out["gt_stretches_differing"] = int(1 + np.count_nonzero(np.diff(bad) > 10)) if bad.size else 0
        out["gt_examples"] = [[int(i) + 1, int(ga[i]), int(gb[i])] for i in bad[:10]]
        del ga, gb
        fa, fb = ca.read("flags"), cb.read("flags")
        names = {
            C.FLAG_PHASED: "PHASED",
            C.FLAG_INS1: "INS1",
            C.FLAG_INS2: "INS2",
            C.FLAG_FILTERED: "FILTERED",
            C.FLAG_VARIANT_SITE: "VARIANT_SITE",
        }
        out["flag_bases_differing"] = {nm: int(((fa ^ fb) & bit != 0).sum()) for bit, nm in names.items()}
        filt = np.flatnonzero((fa ^ fb) & C.FLAG_FILTERED)
        out["filtered_examples"] = [[int(i) + 1, int(fa[i]), int(fb[i])] for i in filt[:5]]
        del fa, fb
        sa = sorted((e.offset, e.hap, e.kind, e.data) for e in ca.sidecar())
        sb = sorted((e.offset, e.hap, e.kind, e.data) for e in cb.sidecar())
        out["sidecar_entries"] = [len(sa), len(sb)]
        out["sidecar_same"] = sa == sb
        enc = next(e for e in b.metadata["contigs"] if e["name"] == chrom).get("encoding") or {}
        out["reencode_conflicts"] = enc.get("conflicts")
    return out


# ------------------------------------------------------------------------------------------------------------ check


def check_chromosome(ref: str, vcf: str, work: str, chrom: str) -> dict:
    sample = os.path.join(work, SAMPLE)
    d = os.path.join(work, chrom)
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    res: dict = {}
    fa = os.path.join(d, "ref.fa")
    res["fasta_export"] = measure(
        ["gpxil", "fasta", ref, "--which", "ref", "--region", chrom, "-o", fa], os.path.join(d, "fasta.log")
    )
    pysam.faidx(fa)
    vfa, res["vcf_reference_masked_bases"] = vcf_reference(vcf, chrom, fa, d)
    export = os.path.join(d, "export.vcf")
    res["vcf_export"] = measure(
        ["gpxil", "vcf", sample, "--reference", ref, "--contigs", chrom, "-o", export], os.path.join(d, "vcf.log")
    )
    res["vcf_export"]["bytes"] = os.path.getsize(export)
    # haplotypes
    differing: list[tuple[int, int]] = []
    res["haplotypes"] = {}
    with GpxilFile(sample, reference=ref) as g, pysam.VariantFile(vcf) as vf:
        c = g[chrom]
        n = c.length
        gt = c.read("gt")
        for h in (1, 2):
            t = time.perf_counter()
            states = C.normalize_states(C.hap1(gt) if h == 1 else C.hap2(gt))
            gev = gpxil_events(states, c.insertions(h))
            del states
            gx = b"".join(iter_haplotype(c, h, keep_case=False))
            cons_fa, chain = os.path.join(d, f"cons{h}.fa"), os.path.join(d, f"chain{h}.txt")
            con = consensus(vfa, vcf, h, cons_fa, chain)
            cs = read_fasta(cons_fa)
            r = compare_haplotype(gx, cs, gev, chain_events(chain), n)
            del gx, cs
            left_out = skipped_alleles(vf, chrom, con["overlaps"], h)
            skipped = Spans((s, e) for s, e, _ in left_out)
            # Where the two differ, bcftools consensus has (a) written a '*' allele that no deletion of its covers into
            # the sequence as a character, or (b) left out an allele that overlaps one it had put on the haplotype,
            # where gpxil's rules (spec 6.4) keep the variant bases of both; or else (c) something else.
            kinds: Counter = Counter()
            others = []
            for a, b, star in r["sites"]:
                kind = (
                    "consensus_wrote_star"
                    if star
                    else "consensus_left_out_overlapping"
                    if skipped.meets(a, b)
                    else "other"
                )
                kinds[kind] += 1
                if kind == "other":
                    others.append([a + 1, b])
            r["differing_sites"] = dict(kinds)
            r["other_sites"] = others[:30]
            r["consensus"] = {
                "applied": con["applied"],
                "left_out_overlapping": len(con["overlaps"]),
                "left_out_not_star": len(left_out),
                "other_messages": con["other_messages"],
            }
            differing += [(a, b) for a, b, _ in r.pop("sites")]
            r["seconds"] = round(time.perf_counter() - t, 1)
            res["haplotypes"][str(h)] = r
        del gt
    # records
    t = time.perf_counter()
    res["records"] = compare_records(vcf, export, chrom, fa, vfa, d, Spans(differing))
    res["records"]["seconds"] = round(time.perf_counter() - t, 1)
    # planes: the export encoded again
    again = os.path.join(d, "again.gpxil.tif")
    res["reencode"] = measure(
        ["gpxil", "encode", ref, "-v", export, "--sex", "XY", "--no-ref", "--contigs", chrom, "-o", again],
        os.path.join(d, "reencode.log"),
    )
    res["planes"] = compare_planes(sample, again, chrom, ref)
    shutil.rmtree(d)
    return res


def step_check(ref: str, vcf: str, work: str, chroms: list[str]) -> None:
    path = os.path.join(RESULTS, "check.json")
    done = {}
    if os.path.exists(path):
        with open(path) as fh:
            done = json.load(fh).get("chromosomes", {})
    for chrom in chroms:
        if chrom in done:
            continue
        t = time.perf_counter()
        done[chrom] = check_chromosome(ref, vcf, work, chrom)
        done[chrom]["seconds"] = round(time.perf_counter() - t, 1)
        save("check.json", {"gpxiltools": gpxiltools.__version__, "pysam": pysam.__version__, "chromosomes": done})
        x = done[chrom]
        print(
            chrom,
            x["seconds"],
            "s; haplotypes differ at",
            x["haplotypes"]["1"]["differing_sites"],
            x["haplotypes"]["2"]["differing_sites"],
            "; records the same",
            x["records"]["same"],
            "of",
            x["records"]["vcf_alleles"],
            "; again: haplotypes the same",
            x["planes"]["haplotypes_same"],
            "gt bases differing",
            x["planes"]["gt_bases_differing"],
            flush=True,
        )


def step_summary() -> None:
    """check.json's chromosomes added up."""
    with open(os.path.join(RESULTS, "check.json")) as fh:
        chroms = json.load(fh)["chromosomes"]
    hap: Counter = Counter()
    kinds: Counter = Counter()
    rec: Counter = Counter()
    rec_kinds: Counter = Counter()
    planes: Counter = Counter()
    flags: Counter = Counter()
    secs: Counter = Counter()
    not_same = []
    other_sites = []
    for name, x in chroms.items():
        for h, r in x["haplotypes"].items():
            for k in (
                "called",
                "unknown",
                "consensus_unknown",
                "differing_bases",
                "differing_stretches",
                "event_stretches",
            ):
                hap[k] += r["stats"].get(k, 0)
            hap["consensus_applied"] += r["consensus"]["applied"] or 0
            hap["consensus_left_out_overlapping"] += r["consensus"]["left_out_overlapping"]
            hap["consensus_left_out_not_star"] += r["consensus"]["left_out_not_star"]
            kinds.update(r["differing_sites"])
            other_sites += [[name, h, *s] for s in r["other_sites"]]
        r = x["records"]
        for k in (
            "vcf_records",
            "export_records",
            "export_nocall_blocks",
            "vcf_alleles",
            "export_alleles",
            "same",
            "same_filter",
        ):
            rec[k] += r[k]
        rec_kinds.update(r["differing"])
        p = x["planes"]
        if not all(p["haplotypes_same"].values()):
            not_same.append(name)
        planes["gt_bases_differing"] += p["gt_bases_differing"]
        planes["gt_stretches_differing"] += p["gt_stretches_differing"]
        planes["chromosomes_sidecar_same"] += bool(p["sidecar_same"])
        planes["sidecar_entries"] += p["sidecar_entries"][0]
        planes["sidecar_entries_again"] += p["sidecar_entries"][1]
        planes["reencode_conflicts"] += p["reencode_conflicts"] or 0
        flags.update(p["flag_bases_differing"])
        for k in ("fasta_export", "vcf_export", "reencode"):
            secs[k] += x[k]["wall_s"]
        secs["vcf_export_bytes"] += x["vcf_export"]["bytes"]
        hap["vcf_reference_masked_bases"] += x.get("vcf_reference_masked_bases", 0)
    out = {
        "chromosomes": list(chroms),
        "haplotypes": dict(hap),
        "haplotype_differing_sites": dict(kinds),
        "haplotype_other_sites": other_sites,
        "records": dict(rec),
        "records_differing": dict(rec_kinds),
        "again": {**planes, "flag_bases_differing": dict(flags), "chromosomes_with_other_haplotypes": not_same},
        "wall_s": {k: round(v, 1) for k, v in secs.items() if k != "vcf_export_bytes"},
        "export_bytes": secs["vcf_export_bytes"],
        "export_peak_rss_mb": max(x["vcf_export"]["peak_rss_mb"] for x in chroms.values()),
        "reencode_peak_rss_mb": max(x["reencode"]["peak_rss_mb"] for x in chroms.values()),
    }
    save("check_summary.json", out)
    print(json.dumps(out, indent=1))


def main() -> None:
    a = sys.argv[1:]
    if a and a[0] == "encode" and len(a) == 5:
        step_encode(*a[1:])
    elif a and a[0] == "check" and len(a) in (4, 5):
        step_check(a[1], a[2], a[3], a[4].split(",") if len(a) == 5 else CHROMS)
    elif a == ["summary"]:
        step_summary()
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
