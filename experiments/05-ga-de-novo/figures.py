"""The experiment's pictures.

- figures/family_*.png: a de novo mutation in a pure GA-sequence, as pixels (gpxi greys: A black, C dark grey,
  T light grey, G white), one row per allele read whole from the HiFi reads: the father's two alleles, then the allele
  each child had from him, on the purine strand; the change outlined in orange.
- figures/rates.png: de novo mutations per base in GA-sequences, against the genome's rate, for indels and
  substitutions, in experiment 05's calls and in the family's paper's (paper_calls.py), with Albrecht-Buehler's
  "160 times fewer" for scale.
- figures/by_length.png: mutations per GA-sequence a transmission, by its length, in both.
- figures/repeats.png: the paper's TR mutations per tandem repeat a transmission, by its length, in GA-sequences
  and in other repeats as long.

    python3 figures.py
"""

import gzip
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402
from PIL import Image, ImageDraw, ImageFont  # noqa: E402

import alleles  # noqa: E402
import ped  # noqa: E402
from analyse import change, edit_distance  # noqa: E402

GREY = {"A": 0, "C": 85, "T": 170, "G": 255}
ORANGE, BLUE = (235, 104, 52), (42, 120, 214)
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
COMP = str.maketrans("ACGT", "TGCA")
EVENTS = [("chr1", 114677757, "200086", "family_substitution.png",
           "A new base: 200086's allele from the father differs from his 158-base allele by one base (G>A on the purine strand)"),
          ("chr12", 28785750, "200082", "family_unit_deletion.png",
           "A unit lost: 200082's allele from the father is his 181-base allele less one 4-base unit (GAAA)")]


def font(size: int):
    for path in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", "/usr/share/fonts/dejavu/DejaVuSans.ttf"):
        if os.path.exists(path):
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def purine(seq: str, kind: str) -> str:
    return seq if kind == "GA" else seq.translate(COMP)[::-1]


def family_image(e: dict, carrier: str, out: str, title: str) -> None:
    father, mother = ped.TRIOS[carrier]
    reads = e["reads"]
    kind = e["ga"][2]
    fa = sorted(alleles.call(reads.get(father, {})), key=lambda a: -a["length"])
    mo = {a["seq"] for a in alleles.call(reads.get(mother, {}))}
    mo_len = {len(s) for s in mo}
    v = [x for x in e["verdicts"] if x["carrier"] == carrier][0]
    new = v["new_allele"]
    rows = [(f"father {father}: {len(a['seq'])}", a["seq"], False) for a in fa]
    kids = [c for c, fm in ped.TRIOS.items() if fm == (father, mother)]
    for c in sorted(kids):
        ca = alleles.call(reads.get(c, {}))
        if c == carrier:
            rows.append((f"{c}: {len(new)}, de novo", new, True))
            continue
        pat = [a for a in ca if a["seq"] not in mo and len(a["seq"]) not in mo_len] or \
              sorted(ca, key=lambda a: min(edit_distance(a["seq"], f["seq"]) for f in fa))[:1]
        if pat:
            rows.append((f"{c}: {len(pat[0]['seq'])}", pat[0]["seq"], False))
    old = min((a["seq"] for a in fa), key=lambda s: edit_distance(s, new))
    ch = change(old, new)
    # purine strand
    rows = [(label, purine(s, kind), is_new) for label, s, is_new in rows]
    old_p, new_p = purine(old, kind), purine(new, kind)
    chp = change(old_p, new_p)
    px, rh, gap, left = 5, 12, 6, 210
    width = left + px * max(len(s) for _, s, _ in rows) + 20
    top = 58
    height = top + len(rows) * (rh + gap) + 40
    img = Image.new("RGB", (width, height), (255, 255, 255))
    d = ImageDraw.Draw(img)
    f12, f11 = font(13), font(11)
    d.text((10, 8), title, fill=(11, 11, 11), font=f12)
    d.text((10, 28), f"{e['chrom']}:{e['ga'][0] + 1:,}-{e['ga'][1]:,} ({'GA' if kind == 'GA' else 'TC'}-sequence, shown on its purine "
                     f"strand with its flanks); A black, C dark grey, T light grey, G white", fill=(82, 81, 78), font=f11)
    for r, (label, s, is_new) in enumerate(rows):
        y = top + r * (rh + gap)
        d.text((10, y), label, fill=(235, 104, 52) if is_new else (11, 11, 11), font=f11)
        for i, b in enumerate(s):
            g = GREY.get(b, 128)
            d.rectangle([left + i * px, y, left + (i + 1) * px - 1, y + rh - 1], fill=(g, g, g))
        d.rectangle([left - 1, y - 1, left + len(s) * px, y + rh], outline=(200, 199, 195))
        if is_new:
            a, b = chp["at"], chp["at"] + max(len(chp["came"]), 1)
            if chp["came"]:
                d.rectangle([left + a * px - 2, y - 3, left + b * px + 1, y + rh + 2], outline=ORANGE, width=2)
            else:  # (a deletion: mark where the bases were lost)
                d.line([left + a * px, y - 4, left + a * px, y + rh + 3], fill=ORANGE, width=3)
        elif s == old_p and chp["gone"]:
            a, b = chp["at"], chp["at"] + len(chp["gone"])
            d.rectangle([left + a * px - 2, y - 3, left + b * px + 1, y + rh + 2], outline=ORANGE, width=2)
    note = (f"the change: {chp['gone'] or '-'} > {chp['came'] or '-'}" + (f" ({ch['kind']})" if ch.get("kind") else ""))
    d.text((left, height - 26), note, fill=(82, 81, 78), font=f11)
    img.save(os.path.join(ped.FIGURES, out))
    print(out, len(rows), "rows", note)


def style(ax, grid: str = "x") -> None:
    ax.set_facecolor(SURFACE)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis=grid, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)


def rates_chart(a: dict, p: dict | None) -> None:
    r = a["rates_g3_g4"]
    t = a["trios"]["G3+G4"]
    rows = [("new indels", "experiment 05", BLUE_HEX, r["indel"]["ratio_to_genome"], r["indel"]["ratio_ci"],
             f"{r['indel']['ratio_to_genome']:g} ×   {r['indel']['de_novo']} in {t} transmissions")]
    if p:
        x = p["tr"]["pass"]["in"]
        rows.append(("new TR mutations", "the paper", ORANGE_HEX, x["ratio_to_papers_genome_rate"], x["ratio_to_papers_genome_rate_ci"],
                     f"{x['ratio_to_papers_genome_rate']:g} ×   {x['mutations']} in {len(p['tr']['children'])} transmissions"))
    rows.append(("new substitutions", "experiment 05", BLUE_HEX, r["SNV"]["ratio_to_genome"], r["SNV"]["ratio_ci"],
                 f"{r['SNV']['ratio_to_genome']:g} ×   {r['SNV']['de_novo']} in {t} transmissions"))
    if p:
        g = p["points"]["SNV"]["germline"]
        hi = 3.689 / g["expected_per_base"]  # (the exact interval of none: 0 to 3.689)
        rows.append(("new substitutions", "the paper", ORANGE_HEX, 0.0, [0.0, round(hi, 1)],
                     f"none, where {g['expected_per_base']:g} would be   ({g['mutations']} in {len(p['points']['SNV']['children'])} transmissions)"))
    left = 0.003
    fig, ax = plt.subplots(figsize=(8.4, 0.62 * len(rows) + 1.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ys = list(range(len(rows) - 1, -1, -1))
    for y, (kind, who, colour, ratio, (lo, hi), label) in zip(ys, rows):
        ax.plot([max(lo, left), hi], [y, y], color=colour, linewidth=2, solid_capstyle="round", zorder=2)
        if ratio > 0:
            ax.scatter([ratio], [y], s=64, color=colour, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        else:
            ax.scatter([left], [y], s=46, marker="<", color=colour, zorder=3, clip_on=False)
        ax.text(1.02, y, label, transform=ax.get_yaxis_transform(), va="center", fontsize=8.2, color=INK)
    ax.axvline(1, color=INK2, linewidth=1, zorder=1)
    ax.text(1.08, len(rows) - 0.45, "the genome's rate", fontsize=7.8, color=INK2, va="center")
    ax.axvline(1 / 160, color=INK2, linewidth=1, linestyle=(0, (3, 2)), zorder=1)
    ax.text(1 / 160 * 1.12, len(rows) - 0.45, "160 times fewer", fontsize=7.8, color=INK2, va="center")
    ax.set_xscale("log")
    ax.set_xlim(left, 400)
    ax.xaxis.set_major_locator(FixedLocator([0.01, 0.1, 1, 10, 100]))
    ax.xaxis.set_major_formatter(FixedFormatter(["0.01", "0.1", "1", "10", "100"]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_ylim(-0.6, len(rows) - 0.2)
    ax.set_yticks(ys)
    ax.set_yticklabels([f"{k} · {w}" for k, w, *_ in rows], fontsize=8.6, color=INK)
    ax.set_xlabel("per base, in pure GA-sequences, against the genome's rate (95 % interval)", fontsize=8.5, color=INK2)
    style(ax)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    if p:
        ax.legend(handles=[Line2D([], [], color=BLUE_HEX, marker="o", linewidth=2, label="experiment 05's calls (G3 and G4)"),
                           Line2D([], [], color=ORANGE_HEX, marker="o", linewidth=2, label="the paper's calls (G3; Porubsky et al. 2025)")],
                  loc="lower left", bbox_to_anchor=(0, 1.0), ncol=2, frameon=False, fontsize=8, labelcolor=INK2,
                  handlelength=1.6, borderaxespad=0.2)
    ax.set_title("New mutations in one family's GA-sequences", fontsize=10.5, color=INK, loc="left", pad=22 if p else 6)
    fig.tight_layout()
    fig.savefig(os.path.join(ped.FIGURES, "rates.png"), facecolor=SURFACE, bbox_inches="tight")
    print("rates.png")


def interval_series(ax, xs, rows, key, colour, offset, label, end: str | None = None) -> None:
    """A dot and its 95 % interval a length bin; ``end`` labels the last dot, in ink, beside it."""
    for x, row in zip(xs, rows):
        if row[key] is None:
            continue
        lo, hi = row["ci"]
        ax.plot([x + offset, x + offset], [lo, hi], color=colour, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([x + offset], [row[key]], s=52, color=colour, edgecolor=SURFACE, linewidth=1.5, zorder=3,
                   label=label if x == xs[0] else None)
    if end:
        ax.annotate(end, (xs[-1] + offset, rows[-1][key]), xytext=(9, 0), textcoords="offset points", va="center",
                    fontsize=7.6, color=INK2)


def length_chart(a: dict, p: dict | None) -> None:
    b = a["by_length"]
    fig, ax = plt.subplots(figsize=(6.8, 3.3), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    xs = list(range(len(b)))
    d = 0.13 if p else 0.0
    interval_series(ax, xs, b, "per_1000_loci_per_transmission", BLUE_HEX, -d,
                    f"experiment 05's calls ({a['trios']['G3+G4']} transmissions; de novo and mosaic)", "experiment 05" if p else None)
    if p:
        q = p["tr"]["pass"]["by_length"]
        interval_series(ax, xs, q, "per_1000_per_transmission", ORANGE_HEX, d,
                        f"the paper's calls ({len(p['tr']['children'])} transmissions; de novo)", "the paper")
        ax.legend(loc="upper left", frameon=False, fontsize=7.8, labelcolor=INK2, handletextpad=0.3)
        ax.set_xlim(-0.5, len(xs) - 0.5 + 0.7)
    ax.set_xticks(xs)
    ax.set_xticklabels([row["length"] for row in b], fontsize=8.5, color=INK)
    ax.set_xlabel("the GA-sequence's length (bases)", fontsize=8.5, color=INK2)
    ax.set_ylabel("mutations a transmission,\nper 1,000 GA-sequences", fontsize=8.5, color=INK2)
    ax.set_ylim(bottom=0)
    style(ax, "y")
    ax.set_title("Longer GA-sequences mutate more", fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(ped.FIGURES, "by_length.png"), facecolor=SURFACE, bbox_inches="tight")
    print("by_length.png")


def repeats_chart(p: dict) -> None:
    c = p["tr"]["compared"]
    series = [("GA-sequences", ORANGE_HEX, "GA-sequences", "GA-sequences"),
              ("other simple repeats (RepeatMasker)", AQUA_HEX, "RepeatMasker's other simple repeats", "other simple\nrepeats")]
    series = [x for x in series if x[0] in c]
    fig, ax = plt.subplots(figsize=(6.8, 3.3), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    xs = list(range(len(c["GA-sequences"]["by_length"])))
    step = 0.26
    for i, (name, colour, label, end) in enumerate(series):
        rows = c[name]["by_length"]
        interval_series(ax, xs, rows, "per_1000_per_transmission", colour, (i - (len(series) - 1) / 2) * step,
                        f"{label} ({c[name]['intervals']:,})", end)
    ax.set_xlim(-0.5, len(xs) - 0.5 + 0.8)
    ax.set_xticks(xs)
    ax.set_xticklabels([row["length"] for row in c["GA-sequences"]["by_length"]], fontsize=8.5, color=INK)
    ax.set_xlabel("the repeat's length (bases)", fontsize=8.5, color=INK2)
    ax.set_ylabel("TR mutations a transmission,\nper 1,000 repeats", fontsize=8.5, color=INK2)
    ax.set_ylim(bottom=0)
    ax.legend(loc="upper left", frameon=False, fontsize=7.8, labelcolor=INK2, handletextpad=0.3)
    style(ax, "y")
    ax.set_title("GA-sequences against other repeats as long (the paper's calls)", fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(ped.FIGURES, "repeats.png"), facecolor=SURFACE, bbox_inches="tight")
    print("repeats.png")


BLUE_HEX, ORANGE_HEX, AQUA_HEX = "#2a78d6", "#eb6834", "#1baf7a"  # (the first three slots of a palette checked for
                                                                   #  colour blindness; aqua is light, so labelled)


def main() -> None:
    os.makedirs(ped.FIGURES, exist_ok=True)
    with gzip.open(os.path.join(ped.RESULTS, "read_alleles.json.gz"), "rt") as fh:
        events = json.load(fh)
    for ch, start, carrier, out, title in EVENTS:
        e = [x for x in events if x["chrom"] == ch and x.get("ga") and x["ga"][0] == start][0]
        family_image(e, carrier, out, title)
    with open(os.path.join(ped.RESULTS, "analysis.json")) as fh:
        a = json.load(fh)
    paper = os.path.join(ped.RESULTS, "paper_calls.json")
    p = None
    if os.path.exists(paper):
        with open(paper) as fh:
            p = json.load(fh)
    rates_chart(a, p)
    length_chart(a, p)
    if p:
        repeats_chart(p)


if __name__ == "__main__":
    main()
