"""The experiment's picture: HG002's heterozygous substitutions by distance from the nearest tandem repeat, for its
four classes of instability; and, at the base right beside a repeat, those that extend it and the others, in HG002
and in new mutations, with the share the reference makes by itself.

    python3 figures.py      -> figures/edges.png
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
RAMP = ["#9ec5f4", "#5598e7", "#256abf", "#0d366b"]  # instability, an ordered ramp of one hue
GREY = "#8d8c88"
CLASSES = ["below 1e-6", "1e-6 to 1e-5", "1e-5 to 1e-4", "1e-4 or more"]
SHORT = ["< 10⁻⁶", "10⁻⁶–10⁻⁵", "10⁻⁵–10⁻⁴", "≥ 10⁻⁴"]


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.8)
    ax.set_axisbelow(True)


def main() -> None:
    with open(os.path.join(HERE, "results", "distances.json")) as fh:
        dist = json.load(fh)
    with open(os.path.join(HERE, "results", "edges.json")) as fh:
        edges = json.load(fh)
    with open(os.path.join(HERE, "results", "ascertainment.json")) as fh:
        asc = json.load(fh)["base 1, by class"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(10.4, 4.3), dpi=150, gridspec_kw={"width_ratios": [1.25, 1]})
    fig.patch.set_facecolor(SURFACE)
    style(a)
    bins = dist["bins"][:-1]  # (the last, 1001+, holds few bases)
    for k, cl in enumerate(CLASSES):
        row = dist["heterozygous"][f"nearest locus {cl}"]
        y = [row[x]["ratio"] for x in bins]
        a.plot(range(len(bins)), y, color=RAMP[k], linewidth=2, marker="o", markersize=4.5, markeredgecolor=SURFACE,
               markeredgewidth=1, label=SHORT[k], zorder=3)
    a.axhline(1, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
    a.set_xticks(range(len(bins)))
    a.set_xticklabels(bins, fontsize=7.2, rotation=0)
    a.grid(axis="y", color=GRID, linewidth=0.8)
    a.set_xlabel("distance from the nearest tandem repeat (bases)", fontsize=8, color=INK2)
    a.set_ylabel("HG002's heterozygous substitutions, observed / expected", fontsize=8, color=INK2)
    a.legend(title="new TR mutations per locus\nand transmission", title_fontsize=7.4, fontsize=7.4, frameon=False, loc="upper right")
    a.set_title("A. An excess only right beside the repeat", fontsize=9, color=INK, loc="left")
    style(b)
    off = {"HG002, extending": -0.24, "HG002, others": 0.0, "new mutations, extending": 0.24}
    for k, cl in enumerate(CLASSES):
        h = edges["loci, HG002"][f"{cl}, base 1, reference breaks the repeat"]
        b.plot([k - 0.36, k - 0.12], [asc[cl]["made by the reference"]] * 2, color=INK, linewidth=1.6, zorder=2)
        b.scatter([k + off["HG002, extending"]], [h["extending"]["ratio"]], s=46, color=RAMP[k], edgecolor=SURFACE, linewidth=1, zorder=4)
        b.scatter([k + off["HG002, others"]], [h["others"]["ratio"]], s=46, color=GREY, edgecolor=SURFACE, linewidth=1, zorder=4)
        b.text(k + off["HG002, extending"] - 0.06, h["extending"]["ratio"], f"{h['extending']['ratio']:.1f}", ha="right", va="center",
               fontsize=7, color=INK2)
        n = edges["loci, new substitutions"][f"{cl}, base 1"]["extending"]
        if n["observed"]:
            x = k + off["new mutations, extending"]
            b.errorbar([x], [n["ratio"]], yerr=[[n["ratio"] - max(n["ci"][0], 0.3)], [n["ci"][1] - n["ratio"]]], fmt="D", markersize=5,
                       color=INK, markerfacecolor=SURFACE, elinewidth=1, capsize=2, zorder=5)
    b.axhline(1, color=INK2, linewidth=0.8, linestyle=(0, (3, 3)), zorder=1)
    b.set_yscale("log")
    b.set_ylim(0.5, 20)
    b.set_yticks([0.5, 1, 2, 4, 8, 16])
    b.set_yticklabels(["0.5", "1", "2", "4", "8", "16"])
    b.set_xticks(range(len(CLASSES)))
    b.set_xticklabels(SHORT, fontsize=7.6)
    b.set_xlim(-0.6, len(CLASSES) - 0.4)
    b.set_xlabel("the repeat's instability (new TR mutations per locus and transmission)", fontsize=8, color=INK2)
    b.set_ylabel("at the base beside it, observed / expected", fontsize=8, color=INK2)
    b.grid(axis="y", color=GRID, linewidth=0.8)
    from matplotlib.lines import Line2D
    handles = [Line2D([], [], marker="o", linestyle="none", markersize=7, markerfacecolor=RAMP[2], markeredgecolor=SURFACE,
                      label="HG002: substitutions extending the repeat"),
               Line2D([], [], marker="o", linestyle="none", markersize=7, markerfacecolor=GREY, markeredgecolor=SURFACE,
                      label="HG002: other substitutions"),
               Line2D([], [], color=INK, linewidth=1.6, label="extending, made by the reference alone"),
               Line2D([], [], marker="D", linestyle="none", markersize=5.5, color=INK, markerfacecolor=SURFACE,
                      label="new mutations: extending (95 % CI)")]
    b.legend(handles=handles, fontsize=7.2, frameon=False, loc="upper left")
    b.set_title("B. Substitutions beside a repeat extend it", fontsize=9, color=INK, loc="left")
    fig.tight_layout()
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "edges.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
