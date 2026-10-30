"""The research article's Fig 1: new TR mutations per 1,000 loci a transmission by purity, in three length bands, for
each composition class (units of 6 bases or fewer); and each class's rate ratio against mixed loci in the Poisson
model of length, unit, purity and class.

    python3 figures_by_class.py      -> figures/purity_by_class.png and .tif

Sized to PLOS's rules for figures: 7.5 inches wide, 300 dpi, an RGB TIFF with LZW compression, and text of 8 to 12
points in Arial, for which Liberation Sans, metrically the same, stands in.

(figures.py draws GA-sequences against the other repeats, as the Software article shows them.)
"""

import json
import os

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image
from scipy.stats import chi2

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
# The classes, in the model's order, with the colours of the alphabets in experiment 09's figure (purines orange, A/T
# blue, A/C green); purine-rich loci other than GA-sequences hollow, mixed loci, the model's reference, dark grey.
CLASSES = [
    ("other", "mixed", "#3d3c39", True),
    ("GA-sequence", "GA-sequences", "#eb6834", True),
    ("purine-rich", "other purine-rich, poly(A) among them", "#eb6834", False),
    ("AT- or CG-rich", "A/T- or C/G-rich, poly(A) aside", "#2a78d6", True),
    ("AC- or GT-rich", "A/C- or G/T-rich", "#1baf7a", True),
]
BANDS = {"51 to 100 bases": (51, 100), "101 to 250 bases": (101, 250), "251 bases or more": (251, 10**9)}
GROUPS = [(0.6, 0.7), (0.7, 0.8), (0.8, 0.85), (0.85, 0.9), (0.9, 1.01)]  #: purity, merged (below 0.6, hardly repeats)
MIN_LOCI = 100  #: a dot only where a group has this many loci
T = 8  #: transmissions
YMAX = 2.0
FIG_W, FIG_H, DPI = 7.5, 5.5, 300
BODY = 8
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": BODY, "axes.linewidth": 0.6})


def poisson_ci(n: int) -> tuple[float, float]:
    return (chi2.ppf(0.025, 2 * n) / 2 if n else 0.0), float(chi2.ppf(0.975, 2 * (n + 1)) / 2)


def axes_at(fig, x: float, y: float, w: float, h: float):
    """Axes placed in inches from the figure's lower left corner."""
    return fig.add_axes((x / FIG_W, y / FIG_H, w / FIG_W, h / FIG_H))


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=BODY, width=0.6, length=2.5)
    ax.set_axisbelow(True)


def marker(ax, x, y, colour, filled, s) -> None:
    ax.scatter([x], [y], s=s, color=colour if filled else SURFACE, edgecolor=colour if not filled else SURFACE,
               linewidth=1.0, zorder=3)


def main() -> None:
    res = os.path.join(HERE, "results")
    s = pd.read_csv(os.path.join(res, "strata.tsv.gz"), sep="\t")
    s["len_lo"] = s["length"].str.split("-").str[0].str.rstrip("+").astype(int)
    s["unit_lo"] = s["motif"].str.split("-").str[0].str.rstrip("+").astype(int)
    s["pur_lo"] = s["purity"].str.split("-").str[0].astype(float)
    s = s[s["unit_lo"] <= 6]  # (units of 6 bases or fewer: the bins 1-1 to 6-6)
    with open(os.path.join(res, "mutability.json")) as fh:
        terms = json.load(fh)["poisson_models"]["class"]["terms"]

    fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=DPI)
    fig.patch.set_facecolor(SURFACE)
    fig.text(0.05 / FIG_W, 5.3 / FIG_H, "A", fontsize=12, fontweight="bold", color=INK, va="baseline")
    fig.text(0.27 / FIG_W, 5.3 / FIG_H, "New TR mutations by purity, in three bands of length", fontsize=9,
             color=INK, va="baseline")
    w, gap, x0 = 2.17, 0.12, 0.72
    axes = [axes_at(fig, x0 + k * (w + gap), 3.3, w, 1.7) for k in range(3)]
    clipped = 0
    for ax, (title, (lo_len, hi_len)) in zip(axes, BANDS.items()):
        style(ax)
        band = s[(s["len_lo"] >= lo_len) & (s["len_lo"] <= hi_len)]
        xs = range(len(GROUPS))
        for k, (cls, name, colour, filled) in enumerate(CLASSES):
            off = (k - (len(CLASSES) - 1) / 2) * 0.15
            sub = band[band["class"] == cls]
            for x, (lo_p, hi_p) in zip(xs, GROUPS):
                g = sub[(sub["pur_lo"] >= lo_p) & (sub["pur_lo"] < hi_p)]
                n, m = int(g["loci"].sum()), int(g["mutations"].sum())
                if n < MIN_LOCI:
                    continue
                lo, hi = poisson_ci(m)
                rate, lo, hi = (1000 * v / n / T for v in (m, lo, hi))
                ax.plot([x + off, x + off], [lo, min(hi, YMAX)], color=colour, linewidth=1.1,
                        solid_capstyle="butt", zorder=2)
                if hi > YMAX:  # an arrowhead where the interval goes on past the top
                    clipped += 1
                    ax.plot([x + off], [YMAX], marker="^", markersize=3.2, color=colour, zorder=2, clip_on=False)
                marker(ax, x + off, rate, colour, filled, 13)
        ax.set_xticks(list(xs))
        ax.set_xticklabels([f"{a:.2f}-{b:.2f}" if b < 1.01 else f"{a:.2f}+" for a, b in GROUPS], fontsize=BODY,
                           color=INK, rotation=35, ha="right", rotation_mode="anchor")
        ax.set_xlim(-0.55, len(GROUPS) - 0.45)
        ax.set_title(title, fontsize=BODY, color=INK, loc="left", pad=3)
        ax.set_ylim(0, YMAX)
        ax.set_yticks([0, 0.5, 1, 1.5, 2])
        ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
        ax.grid(axis="y", color=GRID, linewidth=0.6)
    for ax in axes[1:]:
        ax.tick_params(labelleft=False)
    axes[0].set_ylabel("per 1,000 loci a transmission", fontsize=BODY, color=INK2)
    axes[1].set_xlabel("purity: the share of bases that equal the base one unit on", fontsize=BODY, color=INK2,
                       labelpad=4)

    # the classes, once for both panels
    handles = [
        plt.Line2D([], [], color=colour, marker="o", markersize=4.5, linewidth=1.1,
                   markerfacecolor=colour if filled else SURFACE, markeredgecolor=colour)
        for _, _, colour, filled in CLASSES
    ]
    fig.legend(handles, [c[1] for c in CLASSES], loc="lower left", bbox_to_anchor=(0.66 / FIG_W, 2.25 / FIG_H),
               ncol=3, frameon=False, fontsize=BODY, labelcolor=INK, handlelength=1.6, handletextpad=0.5,
               columnspacing=1.6)

    # B: the class terms
    fig.text(0.05 / FIG_W, 1.98 / FIG_H, "B", fontsize=12, fontweight="bold", color=INK, va="baseline")
    fig.text(0.27 / FIG_W, 1.98 / FIG_H, "Composition in the Poisson model of length, unit, purity and class",
             fontsize=9, color=INK, va="baseline")
    b = axes_at(fig, 2.6, 0.5, 2.7, 1.25)
    style(b)
    rows = [c for c in CLASSES if c[0] != "other"]
    for i, (cls, name, colour, filled) in enumerate(rows):
        t = terms[cls]
        b.plot(t["ci"], [i, i], color=colour, linewidth=1.6, solid_capstyle="butt")
        marker(b, t["rate_ratio"], i, colour, filled, 30)
        b.text(3.35, i, f"{t['rate_ratio']:.2f} ({t['ci'][0]:.2f} to {t['ci'][1]:.2f})", va="center",
               fontsize=BODY, color=INK)
    b.axvline(1, color=INK2, linewidth=0.7, linestyle=(0, (3, 3)), zorder=1)
    b.set_xscale("log")
    b.set_xticks([0.5, 1, 2])
    b.set_xticklabels(["0.5", "1", "2"])
    b.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    b.set_xlim(0.36, 3.0)
    b.set_ylim(len(rows) - 0.5, -0.5)
    b.set_yticks(range(len(rows)))
    b.set_yticklabels([r[1] for r in rows], fontsize=BODY, color=INK)
    b.tick_params(axis="y", length=0, pad=6)
    b.spines["left"].set_visible(False)
    b.grid(axis="x", color=GRID, linewidth=0.6)
    b.set_xlabel("rate ratio against mixed loci (95 %)", fontsize=BODY, color=INK2)

    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    out = os.path.join(HERE, "figures", "purity_by_class")
    fig.savefig(out + ".png", facecolor=SURFACE, dpi=DPI)
    plt.close(fig)
    with Image.open(out + ".png") as im:  # PLOS: RGB, no alpha, LZW
        im.convert("RGB").save(out + ".tif", compression="tiff_lzw", dpi=(DPI, DPI))
    print("purity_by_class.png and .tif; intervals cut at", YMAX, ":", clipped)


if __name__ == "__main__":
    main()
