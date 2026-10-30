"""The experiment's picture: new TR mutations per 1,000 loci a transmission, by purity, in three length bands, for
GA-sequences and for the catalogue's other tandem repeats with units of 6 bases or fewer.

    python3 figures.py      -> figures/purity.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from scipy.stats import chi2  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
ORANGE, AQUA = "#eb6834", "#1baf7a"  # (as experiment 05's repeats.png: GA-sequences orange, other repeats aqua)
TITLES = {"51-100": "51 to 100 bases", "101-250": "101 to 250 bases", "251+": "251 bases or more"}
GROUPS = [(0.6, 0.7), (0.7, 0.8), (0.8, 0.85), (0.85, 0.9), (0.9, 1.01)]  #: purity, merged (below 0.6, hardly repeats)
MIN_LOCI = 100  #: a dot only where a group has this many loci
T = 8  #: transmissions
YMAX = 2.0


def poisson_ci(n: int) -> tuple[float, float]:
    return (chi2.ppf(0.025, 2 * n) / 2 if n else 0.0), float(chi2.ppf(0.975, 2 * (n + 1)) / 2)


def main() -> None:
    with open(os.path.join(HERE, "results", "mutability.json")) as fh:
        d = json.load(fh)
    bands = d["by_purity_in_length_bands"]
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, (band, rows) in zip(axes, bands.items()):
        groups = []
        for lo_edge, hi_edge in GROUPS:  # (the strata's purity bins, merged into broader ones)
            got = [r for r in rows if lo_edge <= float(r["purity"].split("-")[0]) < hi_edge]
            groups.append((f"{lo_edge:.2f}-{hi_edge:.2f}" if hi_edge < 1.01 else f"{lo_edge:.2f}+",
                           {tag: (sum(r[tag]["loci"] for r in got), sum(r[tag]["mutations"] for r in got)) for tag in ("ga", "other")}))
        xs = list(range(len(groups)))
        for tag, colour, off, name in (("ga", ORANGE, 0.14, "GA-sequences"), ("other", AQUA, -0.14, "other repeats")):
            for x, (_, g) in zip(xs, groups):
                n, m = g[tag]
                if n < MIN_LOCI:
                    continue
                lo, hi = poisson_ci(m)
                rate, lo, hi = (1000 * v / n / T for v in (m, lo, hi))
                ax.plot([x + off, x + off], [lo, min(hi, YMAX)], color=colour, linewidth=1.8, solid_capstyle="round", zorder=2)
                ax.scatter([x + off], [rate], s=30, color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=3,
                           label=name if x == xs[0] else None)
        ax.set_xticks(xs)
        ax.set_xticklabels([g[0] for g in groups], fontsize=7.2, color=INK, rotation=30)
        ax.set_title(TITLES[band], fontsize=9.5, color=INK, loc="left")
        ax.set_ylim(0, YMAX)
        ax.set_facecolor(SURFACE)
        ax.tick_params(colors=INK2, labelsize=7.5)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color(GRID)
    axes[0].set_ylabel("new TR mutations a transmission,\nper 1,000 loci (95 %)", fontsize=8.5, color=INK2)
    axes[1].set_xlabel("purity: the share of bases that equal the base one unit on", fontsize=8.5, color=INK2)
    axes[0].legend(loc="upper left", frameon=False, fontsize=8, labelcolor=INK2, handletextpad=0.3)
    fig.suptitle("GA-sequences mutate as other repeats as long and as pure do (units of 6 bases or fewer)", fontsize=10.5, color=INK, x=0.01, ha="left")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "purity.png"), facecolor=SURFACE, bbox_inches="tight")
    print("purity.png")


if __name__ == "__main__":
    main()
