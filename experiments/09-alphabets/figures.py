"""The experiment's picture: Alu tails that are pure runs in each alphabet, by the Alu's age; and the share of new
substitutions that keep a run of each alphabet pure.

    python3 figures.py      -> figures/alphabets.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
COLOUR = {"purine/pyrimidine": "#eb6834", "A/T": "#2a78d6", "A/C": "#1baf7a", "G/C": "#8d8c88"}  # (purines orange, as in 05 to 08)
NAMES = {"purine/pyrimidine": "purines (GA)", "A/T": "A/T", "A/C": "A/C", "G/C": "G/C"}


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.8)
    ax.set_axisbelow(True)


def main() -> None:
    with open(os.path.join(HERE, "results", "ages.json")) as fh:
        ages = json.load(fh)["Alu"]["by_divergence"]
    with open(os.path.join(HERE, "results", "spectrum.json")) as fh:
        spec = json.load(fh)["all"]["runs"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10, 4.2), dpi=150, gridspec_kw={"width_ratios": [1.45, 1]})
    fig.patch.set_facecolor(SURFACE)
    bins = list(ages)
    x = range(len(bins))
    style(a)
    for name in ("purine/pyrimidine", "A/T", "A/C"):
        y = [ages[k][name]["per_1000"] for k in bins]
        a.plot(x, y, color=COLOUR[name], linewidth=2, marker="o", markersize=4.5, markeredgecolor=SURFACE, markeredgewidth=1,
               label=NAMES[name], zorder=3)
        a.text(len(bins) - 0.75, y[-1], NAMES[name], color=INK2, fontsize=7.8, va="center")
    a.set_xticks(list(x))
    a.set_xticklabels([k.replace(" %", "") for k in bins], fontsize=7.4)
    a.set_xlim(-0.4, len(bins) + 0.9)
    a.set_ylim(0, None)
    a.grid(axis="y", color=GRID, linewidth=0.8)
    a.set_xlabel("Alu's divergence from its consensus (%): its age", fontsize=8, color=INK2)
    a.set_ylabel("tails that are a pure run, per 1,000 Alus", fontsize=8, color=INK2)
    a.legend(loc="upper left", fontsize=7.6, frameon=False)
    a.set_title("A. The same tails, opposite fates", fontsize=9, color=INK, loc="left")
    style(b)
    rows = [("purine/pyrimidine", "purines, 3 A to 1 G", spec["purine, 3 A to 1 G"]["kept"]),
            ("A/C", "A/C, 3 A to 1 C", spec["A/C, 3 A to 1 C"]["kept"]),
            ("G/C", "G/C", spec["G/C, outside CpG"]["kept"]),
            ("A/T", "A/T", spec["A/T"]["kept"])]
    for i, (name, lab, v) in enumerate(rows):
        b.barh(i, v, height=0.56, color=COLOUR[name], zorder=3)
        b.text(v + 0.012, i, f"{100 * v:.0f} %", va="center", fontsize=7.8, color=INK2)
    b.set_yticks(range(len(rows)))
    b.set_yticklabels([r[1] for r in rows], fontsize=8, color=INK)
    b.invert_yaxis()
    b.tick_params(axis="y", length=0)
    b.set_xlim(0, 0.75)
    b.grid(axis="x", color=GRID, linewidth=0.8)
    b.set_xlabel("share of new substitutions in it that keep it pure", fontsize=8, color=INK2)
    b.set_title("B. Transitions keep a purine run pure", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "alphabets.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
