"""The experiment's pictures.

- figures/family_*.png: a de novo mutation in a pure GA-sequence, as pixels (gpxi greys: A black, C dark grey,
  T light grey, G white), one row per allele read whole from the HiFi reads: the father's two alleles, then the allele
  each child had from him, on the purine strand; the change outlined in orange.
- figures/rates.png: de novo mutations per base in GA-sequences, against the genome's rate, for indels and
  substitutions, with Albrecht-Buehler's "160 times fewer" for scale.
- figures/by_length.png: mutations per GA-sequence a transmission, by its length.

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


def rates_chart(a: dict) -> None:
    r = a["rates_g3_g4"]
    rows = [("de novo indels", r["indel"]), ("de novo substitutions", r["SNV"])]
    fig, ax = plt.subplots(figsize=(7.6, 2.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for y, (label, x) in zip((1, 0), rows):
        lo, hi = x["ratio_ci"]
        ax.plot([max(lo, 0.003), hi], [y, y], color=BLUE_HEX, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([max(x["ratio_to_genome"], 0.003)], [y], s=70, color=BLUE_HEX, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text(hi * 1.25, y, f"{x['ratio_to_genome']:g} × ({x['de_novo']} in {a['trios']['G3+G4']} transmissions)", va="center",
                fontsize=8.5, color=INK)
    ax.axvline(1, color=INK2, linewidth=1, zorder=1)
    ax.text(1.06, 1.45, "the genome's rate", fontsize=8, color=INK2, va="center")
    ax.axvline(1 / 160, color=ORANGE_HEX, linewidth=1, linestyle=(0, (3, 2)), zorder=1)
    ax.text(1 / 160 * 1.1, 1.45, "160 times fewer", fontsize=8, color=ORANGE_HEX, va="center")
    ax.set_xscale("log")
    ax.set_xlim(0.003, 400)
    ax.set_ylim(-0.6, 1.7)
    ax.set_yticks([1, 0])
    ax.set_yticklabels([x[0] for x in rows], fontsize=9, color=INK)
    ax.set_xlabel("per base, in pure GA-sequences, against the genome (95 % interval)", fontsize=8.5, color=INK2)
    ax.tick_params(axis="x", colors=INK2, labelsize=8)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_title("New mutations in one family's GA-sequences", fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(ped.FIGURES, "rates.png"), facecolor=SURFACE, bbox_inches="tight")
    print("rates.png")


def length_chart(a: dict) -> None:
    b = a["by_length"]
    fig, ax = plt.subplots(figsize=(6.4, 3.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    xs = list(range(len(b)))
    for x, row in zip(xs, b):
        lo, hi = row["ci"]
        ax.plot([x, x], [lo, hi], color=BLUE_HEX, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([x], [row["per_1000_loci_per_transmission"]], s=60, color=BLUE_HEX, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.text(x + 0.08, row["per_1000_loci_per_transmission"], f"{row['events']} of {row['ga_sequences']:,}", fontsize=7.5, color=INK2,
                va="center")
    ax.set_xticks(xs)
    ax.set_xticklabels([row["length"] for row in b], fontsize=8.5, color=INK)
    ax.set_xlabel("the GA-sequence's length (bases)", fontsize=8.5, color=INK2)
    ax.set_ylabel("mutations a transmission,\nper 1,000 GA-sequences", fontsize=8.5, color=INK2)
    ax.tick_params(colors=INK2, labelsize=8)
    ax.grid(axis="y", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(GRID)
    ax.set_title("Longer GA-sequences mutate more", fontsize=10.5, color=INK, loc="left")
    fig.tight_layout()
    fig.savefig(os.path.join(ped.FIGURES, "by_length.png"), facecolor=SURFACE, bbox_inches="tight")
    print("by_length.png")


BLUE_HEX, ORANGE_HEX = "#2a78d6", "#eb6834"


def main() -> None:
    os.makedirs(ped.FIGURES, exist_ok=True)
    with gzip.open(os.path.join(ped.RESULTS, "read_alleles.json.gz"), "rt") as fh:
        events = json.load(fh)
    for ch, start, carrier, out, title in EVENTS:
        e = [x for x in events if x["chrom"] == ch and x.get("ga") and x["ga"][0] == start][0]
        family_image(e, carrier, out, title)
    with open(os.path.join(ped.RESULTS, "analysis.json")) as fh:
        a = json.load(fh)
    rates_chart(a)
    length_chart(a)


if __name__ == "__main__":
    main()
