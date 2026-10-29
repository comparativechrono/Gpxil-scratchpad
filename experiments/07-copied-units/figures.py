"""The experiment's picture: how much more often a substitution happens where it would copy its repeat than where the
same substitution, in the same 6-base window and the same kind of locus, would not (Mantel-Haenszel rate ratio, 95 %
interval), for the new substitutions of the families, split by whether gnomAD already holds the new allele, and for
gnomAD's own minor alleles by allele count, in the GA-sequences and in other tandem repeats.

    python3 figures.py      -> figures/copying.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator, NullLocator  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"  # (as experiments 05 and 06: GA-sequences orange, other repeats aqua)
LAYER = "purity, length, unit"
CLASSES = (("singleton", "one genome"), ("rare", "under 0.1 %"), ("uncommon", "0.1 to 1 %"), ("common", "1 % or more"))


def rows() -> list:
    with open(os.path.join(HERE, "results", "matched_context.json")) as fh:
        mc = json.load(fh)["tests"]
    with open(os.path.join(HERE, "results", "population.json")) as fh:
        pop = json.load(fh)
    out = [("New substitutions in the families (repeats of purity 0.6 or more)", None, None)]
    for key, label in (("purity 0.6+", "all"), ("purity 0.6+, allele absent from gnomAD", "new allele not in gnomAD"),
                       ("purity 0.6+, allele in gnomAD", "new allele already in gnomAD")):
        t = mc[key][LAYER]
        out.append((f"{label} ({t['mutations']:,})", t, BLUE))
    for part, title, colour in (("ga_sequences", "gnomAD's minor alleles in the GA-sequences, by frequency", ORANGE),
                                ("tandem_repeats", "gnomAD's minor alleles in 20,000 other tandem repeats", AQUA)):
        out.append((title, None, None))
        for cl, label in CLASSES:
            t = pop[part]["minor_alleles_by_count"][cl]["matched_context"]
            out.append((f"{label} ({t['mutations']:,})", t, colour))
    return out


def main() -> None:
    data = rows()
    fig, ax = plt.subplots(figsize=(8.2, 0.34 * len(data) + 0.9), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    y = len(data)
    ticks, labels = [], []
    for label, t, colour in data:
        y -= 1
        if t is None:
            ax.text(0.31, y, label, fontsize=8.6, color=INK, fontweight="bold", va="center", ha="left", zorder=4,
                    bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 1.5})
            continue
        lo, hi = t["ci"]
        ax.plot([lo, hi], [y, y], color=colour, linewidth=2, solid_capstyle="round", zorder=2)
        ax.scatter([t["rate_ratio"]], [y], s=34, color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        ax.text(hi * 1.04, y, f"{t['rate_ratio']:.2f}", fontsize=7.6, color=INK2, va="center")
        ticks.append(y)
        labels.append(label)
    ax.axvline(1, color=INK2, linewidth=0.9, zorder=1)
    ax.set_xscale("log")
    xt = [0.5, 0.75, 1, 1.5, 2, 3, 4, 6]
    ax.xaxis.set_major_locator(FixedLocator(xt))
    ax.xaxis.set_major_formatter(FixedFormatter([f"{v:g}" for v in xt]))
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlim(0.3, 9)
    ax.set_ylim(-0.7, len(data) - 0.3)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=7.8, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=INK2, labelsize=7.8)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("rate ratio: where the substitution would copy the repeat, against where it would not", fontsize=8, color=INK2)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "copying.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
