"""The experiment's picture: the Alu tails stacked at the Alu, by age; how pure each tenth of them is; and at which end
their alleles differ in length.

    python3 figures.py      -> figures/along_arrays.png
"""

import json
import os

import matplotlib
import matplotlib.patches

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
AGES = ["0-4 %", "4-8 %", "8-12 %", "12-16 %", "16-20 %", "20 %+"]
RAMP = ["#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c"]  #: young (light) to old (dark)
SOURCE = {"family": ("#2a78d6", "o", "the family's parents (TRGT)"), "HG002": ("#eb6834", "D", "HG002 (GIAB v5.0q)")}


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.8)
    ax.set_axisbelow(True)


def main() -> None:
    res = os.path.join(HERE, "results")
    px = np.load(os.path.join(res, "stack_rows.npy"))
    with open(os.path.join(res, "stack_profile.json")) as fh:
        sp = json.load(fh)
    with open(os.path.join(res, "tails.json")) as fh:
        tails = json.load(fh)
    with open(os.path.join(res, "lengths.json")) as fh:
        lengths = json.load(fh)["tails"]
    fig = plt.figure(figsize=(12, 5.0), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.1, 1.05], wspace=0.62)
    # A: the stack
    a = fig.add_subplot(gs[0])
    a.imshow(px, aspect="auto", interpolation="nearest", extent=(-sp["up"] - 0.5, sp["down"] - 0.5, len(px), 0))
    starts = np.array(sp["bin_divergence_start"])
    ticks = [int(np.searchsorted(starts, d)) for d in (4, 8, 12, 16, 20)]
    a.set_yticks(ticks)
    a.set_yticklabels([f"{d} %" for d in (4, 8, 12, 16, 20)])
    a.axvline(-0.5, color=INK, linewidth=0.8)
    a.set_xlabel("bases from the tail's start (Alu to the left)", fontsize=8, color=INK2)
    a.set_ylabel("the Alu's divergence: its age", fontsize=8, color=INK2)
    a.tick_params(colors=INK2, labelsize=7.4)
    a.set_title("A. Alu tails, lined up at the Alu", fontsize=9, color=INK, loc="left")
    from gpxiltools import codes as GC
    from gpxiltools.palettes import ref_colours
    rgb, _ = ref_colours("classic")
    keys = [matplotlib.patches.Patch(facecolor=rgb[GC.BASE_A | GC.SOFTMASK] / 255, edgecolor=GRID, label="A (in a repeat)")]
    keys += [matplotlib.patches.Patch(facecolor=rgb[code] / 255, edgecolor=GRID, label=b)
             for b, code in (("C", GC.BASE_C), ("G", GC.BASE_G), ("T", GC.BASE_T))]
    a.legend(handles=keys, loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=4, fontsize=6.8, frameon=False,
             handlelength=1, columnspacing=0.8, handletextpad=0.4)
    # B: A by tenth
    b = fig.add_subplot(gs[1])
    style(b)
    x = np.arange(10) + 0.5
    alu = tails["arrays"]["Alu"]["all arrays"]
    ends = []
    for age, colour in zip(AGES, RAMP):
        y = alu[age]["a_by_tenth"]
        b.plot(x, y, color=colour, linewidth=2, marker="o", markersize=3.5, markeredgecolor=SURFACE, markeredgewidth=0.8, zorder=3)
        ends.append([y[-1], age])
    for i in range(1, len(ends)):  # (labels kept apart, top down)
        ends[i][0] = min(ends[i][0], ends[i - 1][0] - 0.014)
    for y_, age in ends:
        b.text(10.25, y_, age, color=INK2, fontsize=7.2, va="center")
    b.set_xlim(0, 11.6)
    b.set_xticks([0, 5, 10])
    b.set_xticklabels(["next to\nthe Alu", "middle", "far\nend"])
    b.set_ylim(0.68, 1.0)
    b.grid(axis="y", color=GRID, linewidth=0.8)
    b.set_ylabel("share of A, each tenth of the array", fontsize=8, color=INK2)
    b.set_title("B. Richest in A next to the Alu, at every age", fontsize=9, color=INK, loc="left")
    # C: which end differs
    c = fig.add_subplot(gs[2])
    style(c)
    rows = [("Alu tails", "Alu", "all", None),
            ("purest stretch\nwholly in near half", "Alu", "all", "longest in frame wholly in the near half"),
            ("purest stretch\nwholly in far half", "Alu", "all", "longest in frame wholly in the far half"),
            ("L1 tails", "L1", "all", None),
            ("purest stretch\nwholly in near half", "L1", "all", "longest in frame wholly in the near half"),
            ("purest stretch\nwholly in far half", "L1", "all", "longest in frame wholly in the far half")]
    for i, (label, fam, age, part) in enumerate(rows):
        for dy, (src, (colour, marker, _)) in zip((-0.13, 0.13), SOURCE.items()):
            d = lengths[src][fam][age]
            if part:
                d = d[part]
            share = d["near"] / (d["near"] + d["far"])
            c.plot(share, i + dy, marker=marker, color=colour, markersize=6.5, markeredgecolor=SURFACE, markeredgewidth=1, zorder=3)
    c.axvline(0.5, color=INK2, linewidth=0.9, linestyle=(0, (3, 3)))
    c.text(0.51, len(rows) - 0.42, "each end alike", color=INK2, fontsize=7.2)
    c.axhline(2.5, color=GRID, linewidth=0.8)
    c.set_yticks(range(len(rows)))
    c.set_yticklabels([r[0] for r in rows], fontsize=7.8, color=INK)
    c.invert_yaxis()
    c.set_ylim(len(rows) - 0.4, -0.6)
    c.tick_params(axis="y", length=0)
    c.set_xlim(0.0, 1.0)
    c.grid(axis="x", color=GRID, linewidth=0.8)
    c.set_xlabel("of the length differences placed in one half,\nthe share in the half next to the element", fontsize=8, color=INK2)
    for src, (colour, marker, name) in SOURCE.items():
        c.plot([], [], marker=marker, color=colour, linestyle="none", markersize=6, label=name)
    c.legend(loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=2, fontsize=7.2, frameon=False)
    c.set_title("C. Length changes follow the purest stretch", fontsize=9, color=INK, loc="right")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.9, bottom=0.25)
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "along_arrays.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
