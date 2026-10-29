"""The experiment's picture: the 19 markers of Liu et al. 2021, by mutation rate, the GA-sequences among them marked,
and those that begin at an Alu's 3' end.

    python3 figures.py      -> figures/markers.png
"""

import json
import os
from decimal import ROUND_HALF_UP, Decimal

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
ORANGE, GREY = "#eb6834", "#8d8c88"  # (GA-sequences orange, as in experiments 05 to 07)
NAMES = {"PENTAD": "Penta D", "PENTAE": "Penta E", "VWA": "vWA"}


def main() -> None:
    with open(os.path.join(HERE, "results", "markers.json")) as fh:
        rows = [r for r in json.load(fh)["markers_table"].values() if r.get("liu2021")]
    rows.sort(key=lambda r: r["liu2021"][0])
    fig, ax = plt.subplots(figsize=(6.4, 5.2), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    for y, r in enumerate(rows):
        rate = 1000 * r["liu2021"][0]
        colour = ORANGE if r["ga_sequence"] else GREY
        ax.plot([0, rate], [y, y], color=GRID, linewidth=1.2, zorder=1)
        ax.scatter([rate], [y], s=40, color=colour, edgecolor=SURFACE, linewidth=1.2, zorder=3)
        shown = f"{float(Decimal(repr(rate)).quantize(Decimal('0.01'), rounding=ROUND_HALF_UP)):.2f}"  # (3.755 shows as 3.76)
        note = shown + (f"   {r['alu_at_start']} tail" if r.get("alu_at_start") else "")
        ax.text(rate + 0.08, y, note, fontsize=7.4, color=INK2, va="center")
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([NAMES.get(r["marker"], r["marker"]) for r in rows], fontsize=8, color=INK)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=INK2, labelsize=7.8)
    ax.set_xlim(0, 6.2)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("mutations per 1,000 transmissions (Liu et al. 2021)", fontsize=8, color=INK2)
    ax.scatter([], [], s=40, color=ORANGE, label="in a pure GA-sequence")
    ax.scatter([], [], s=40, color=GREY, label="not")
    ax.legend(loc="lower right", fontsize=7.6, frameon=False)
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "markers.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
