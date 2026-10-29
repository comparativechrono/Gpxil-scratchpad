"""The experiment's chart: how often transitions in pure GA-sequences copy a neighbouring 6-base unit, observed against
expected, for HG002 (two benchmarks; all transitions and isolated ones) and for gnomAD by allele frequency.

    python3 figures.py          -> figures/unit_copy.png (from results/unit_copy_*.json, controls_*.json, gnomad.json)
"""

import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import ga  # noqa: E402

INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#ffffff"
OBSERVED, CONTEXT = "#2a78d6", "#eb6834"  # (categorical slots 1 and 2, validated)


def load(name):
    path = os.path.join(ga.RESULTS, name)
    return json.load(open(path)) if os.path.exists(path) else None


def main() -> None:
    rows = []  # (label, expected from the bases, expected from contexts or None, observed, n)
    for cs in ("v4.2.1", "v5.0q"):
        u, c = load(f"unit_copy_{cs}.json"), load(f"controls_{cs}.json")
        if not u:
            continue
        k6 = u["k"]["6"]
        ctx = c["own"]["share_context"] if c else None
        for g, label in (("all", "all"), ("isolated", "isolated")):
            x = k6[g]
            rows.append((f"HG002 {cs}, {label} ({x['ti']:,})", k6["expected_share"], ctx, x["share"], x["ti"]))
    gn = load("gnomad.json")
    if gn:
        u = gn["unit_copy"]
        for g, label in (("singleton", "singletons"), ("rare", "AF < 0.1 %"), ("uncommon", "AF 0.1-1 %"), ("common", "AF ≥ 1 %")):
            x = u[g]
            rows.append((f"gnomAD, {label} ({x['ti']:,})", u["expected_share"], None, x["share"], x["ti"]))
    fig, ax = plt.subplots(figsize=(8.2, 0.42 * len(rows) + 1.5), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)
    # (top to bottom, with a gap between HG002's rows and gnomAD's)
    ys, y = [], 0.0
    for i, r in enumerate(rows):
        if i and r[0].split()[0] != rows[i - 1][0].split()[0]:
            y -= 0.6
        ys.append(y)
        y -= 1
    for y, (label, exp, ctx, obs, n) in zip(ys, rows):
        ax.plot([100 * exp, 100 * obs], [y, y], color=GRID, linewidth=2, zorder=1, solid_capstyle="round")
        ax.scatter([100 * exp], [y], s=64, facecolor=SURFACE, edgecolor=INK2, linewidth=1.5, zorder=3)
        if ctx is not None:
            ax.scatter([100 * ctx], [y], s=64, marker="D", facecolor=CONTEXT, edgecolor=SURFACE, linewidth=1.5, zorder=3)
        ax.scatter([100 * obs], [y], s=72, facecolor=OBSERVED, edgecolor=SURFACE, linewidth=1.5, zorder=4)
        ax.text(100 * obs + 0.8, y, f"{100 * obs:.0f} %", va="center", ha="left", fontsize=8.5, color=INK)
    ax.set_yticks(ys)
    ax.set_yticklabels([r[0] for r in rows], fontsize=9, color=INK)
    low = min(min(r[1], r[3], r[2] if r[2] is not None else 1) for r in rows)
    ax.set_xlim(min(50, 10 * int(10 * low - 0.3)), 100)
    ax.set_xlabel("transitions that copy a neighbouring 6-base unit of their GA-sequence (%)", fontsize=9, color=INK2)
    ax.tick_params(axis="x", colors=INK2, labelsize=8.5)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    handles = [
        plt.Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=SURFACE, markeredgecolor=INK2, label="expected: share of bases where a transition copies"),
        plt.Line2D([], [], marker="D", linestyle="", markersize=7, markerfacecolor=CONTEXT, markeredgecolor=SURFACE, label="expected: bases weighted by 7-base context rate"),
        plt.Line2D([], [], marker="o", linestyle="", markersize=8, markerfacecolor=OBSERVED, markeredgecolor=SURFACE, label="observed"),
    ]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.42, -0.14 - 0.25 / len(rows)), ncol=1, frameon=False, fontsize=8.5,
              labelcolor=INK)
    ax.set_title("Transitions in pure GA-sequences copy their neighbours", fontsize=11, color=INK, loc="left", pad=10)
    fig.tight_layout()
    os.makedirs(ga.FIGURES, exist_ok=True)
    out = os.path.join(ga.FIGURES, "unit_copy.png")
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    print(out, len(rows), "rows")


def alu_clock() -> None:
    """Two panels on one age axis: how many Alu elements end in a pure GA-sequence, and how much of such a tail is A."""
    a = load("alu.json")
    if not a:
        return
    rows = [r for r in a["alu_tails_by_divergence"] if r["alus"] >= 1000]
    x = [2 * i + 1 for i in range(len(rows))]  # (the middle of each 2 % step)
    per = [r["per_1000"] for r in rows]
    afr = [100 * r["tail_a_fraction"] for r in rows]
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.2), dpi=150, sharex=True)
    fig.patch.set_facecolor(SURFACE)
    for ax, ys, title, fmt in ((axes[0], per, "Alu elements that end in a pure GA-sequence, per 1,000", "{:.0f}"),
                               (axes[1], afr, "Of such a tail, A (%)", "{:.0f} %")):
        ax.set_facecolor(SURFACE)
        ax.plot(x, ys, color=OBSERVED, linewidth=2, marker="o", markersize=5, markerfacecolor=OBSERVED, markeredgecolor=SURFACE)
        for xi, yi in ((x[0], ys[0]), (x[-1], ys[-1])):
            ax.annotate(fmt.format(yi), (xi, yi), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8.5, color=INK)
        ax.set_title(title, fontsize=9.5, color=INK, loc="left")
        ax.set_xlabel("the Alu's divergence from its consensus (%): its age", fontsize=8.5, color=INK2)
        ax.tick_params(colors=INK2, labelsize=8)
        ax.grid(axis="y", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(GRID)
    axes[0].set_ylim(0, max(per) * 1.2)
    axes[1].set_ylim(50, 100)
    fig.tight_layout()
    out = os.path.join(ga.FIGURES, "alu_clock.png")
    fig.savefig(out, facecolor=SURFACE, bbox_inches="tight")
    print(out)


if __name__ == "__main__":
    main()
    alu_clock()
