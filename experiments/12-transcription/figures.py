"""The experiment's picture: how much more often a repeat varies in the family's parents and in HG002 when a gene's
transcript carries its A-, purine- or G-rich strand than when the polymerase reads it; and how much less genic repeats
vary than repeats between genes.

    python3 figures.py      -> figures/transcription.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
SOURCE = {"family": ("#2a78d6", "o", "the family's parents (TRGT)"), "hg002": ("#eb6834", "D", "HG002 (GIAB v5.0q)")}
GROUPS = [("poly(A)", "poly(A)"), ("purine", "purine"), ("A-rich A/T", "A-rich A/T"), ("G-rich", "G-rich"),
          ("A/T alternating", "A/T alternating\n(no strand)")]


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.8)
    ax.set_axisbelow(True)


def point(ax, x: dict, y: float, colour: str, marker: str, filled: bool = True) -> None:
    if not x.get("odds_ratio"):
        return
    lo, hi = x["ci"]
    ax.plot([lo, hi], [y, y], color=colour, linewidth=1.4, zorder=2)
    ax.plot(x["odds_ratio"], y, marker=marker, color=colour, markersize=6.5, zorder=3,
            markerfacecolor=colour if filled else SURFACE, markeredgecolor=colour if not filled else SURFACE, markeredgewidth=1.3)


def main() -> None:
    with open(os.path.join(HERE, "results", "strands.json")) as fh:
        g = json.load(fh)["groups"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.5, 4.2), dpi=150, sharey=True)
    fig.patch.set_facecolor(SURFACE)
    for ax in (a, b):
        style(ax)
        ax.set_xscale("log")
        ax.axvline(1, color=INK2, linewidth=0.9, linestyle=(0, (3, 3)))
        ax.grid(axis="x", color=GRID, linewidth=0.8, which="both")
    for i, (key, _) in enumerate(GROUPS):
        for dy, src in zip((-0.14, 0.14), SOURCE):
            colour, marker, _ = SOURCE[src]
            point(a, g[key]["sense vs antisense"][src], i + dy, colour, marker)
        for dy, src in zip((-0.22, 0.22), SOURCE):
            colour, marker, _ = SOURCE[src]
            point(b, g[key]["sense vs between genes"][src], i + dy - 0.06, colour, marker, filled=True)
            point(b, g[key]["antisense vs between genes"][src], i + dy + 0.06, colour, marker, filled=False)
    a.set_yticks(range(len(GROUPS)))
    a.set_yticklabels([lab for _, lab in GROUPS], fontsize=8, color=INK)
    a.invert_yaxis()
    a.tick_params(axis="y", length=0)
    for ax in (a, b):
        ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%.2f"))
        ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a.set_xticks([0.9, 1.0, 1.1, 1.2])
    a.set_xlim(0.84, 1.24)
    a.set_xlabel("odds of a length variant, sense over antisense\n(matched on length, motif size and purity; 95 % CI)", fontsize=8, color=INK2)
    a.set_title("A. Repeats vary more when the transcript\ncarries their A- or purine-rich strand", fontsize=9, color=INK, loc="left")
    for src, (colour, marker, name) in SOURCE.items():
        a.plot([], [], marker=marker, color=colour, linestyle="none", markersize=6, label=name)
    a.legend(loc="lower right", fontsize=7.2, frameon=False)
    b.set_xticks([0.8, 0.9, 1.0, 1.1])
    b.set_xlim(0.76, 1.12)
    b.set_xlabel("odds of a length variant, in genes over between genes\n(filled: sense; open: antisense)", fontsize=8, color=INK2)
    b.set_title("B. Repeats in genes vary less, most when the\npolymerase reads their A- or purine-rich strand", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "transcription.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
