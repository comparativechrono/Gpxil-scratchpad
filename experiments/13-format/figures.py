"""The experiment's picture: what the files take on disk, and how fast a region comes out of them.

python3 figures.py      -> figures/format.png
"""

import json
import os

import matplotlib
import matplotlib.ticker

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e6e5e1", "#ffffff"
GPXIL, OTHER = "#2a78d6", "#a9a8a3"
MB = 1e6


def style(ax) -> None:
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID)
    ax.tick_params(colors=INK2, labelsize=7.8)
    ax.set_axisbelow(True)


def bars(ax, rows, title, unit_label) -> None:
    """Horizontal bars, gpxil's in colour, the others in grey, each labelled with its size."""
    style(ax)
    names = [r[0] for r in rows]
    vals = [r[1] / MB for r in rows]
    colours = [GPXIL if r[2] else OTHER for r in rows]
    y = range(len(rows))
    ax.barh(y, vals, color=colours, height=0.62)
    top = max(vals)
    for i, v in zip(y, vals):
        ax.text(v + top * 0.015, i, f"{v:,.1f}", va="center", fontsize=7.4, color=INK2)
    ax.set_yticks(list(y))
    ax.set_yticklabels(names, fontsize=7.8, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, top * 1.2)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:,.0f}"))
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel(unit_label, fontsize=8, color=INK2)
    ax.set_title(title, fontsize=9, color=INK, loc="left")


def reads(ax, local, pairs, title) -> None:
    """Median time to read a region, by its size: each kind of read a line, gpxil's in colour."""
    style(ax)
    sizes = sorted(local["reads"], key=int)
    x = list(range(len(sizes)))
    for key, label, colour, dash in pairs:
        y = [local["reads"][s][key]["median_ms"] for s in sizes]
        ax.plot(
            x,
            y,
            color=colour,
            linewidth=2,
            linestyle=dash,
            marker="o",
            markersize=4,
            markeredgecolor=SURFACE,
            markeredgewidth=0.8,
        )
        ax.text(x[-1] + 0.12, y[-1], label, fontsize=7.4, color=INK2, va="center")
    ax.set_yscale("log")
    ticks = [t for t in (0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50) if ax.get_ylim()[0] <= t <= ax.get_ylim()[1]]
    ax.set_yticks(ticks)
    ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks(x)
    ax.set_xticklabels(["1 kb", "100 kb", "1 Mb"])
    ax.set_xlim(-0.3, len(x) - 1 + 0.9)
    ax.grid(axis="y", color=GRID, linewidth=0.8, which="major")
    ax.set_xlabel("size of the region read (each placed at random, 100 of each)", fontsize=8, color=INK2)
    ax.set_ylabel("median time per region (ms)", fontsize=8, color=INK2)
    ax.set_title(title, fontsize=9, color=INK, loc="left")


def main() -> None:
    res = os.path.join(HERE, "results")
    with open(os.path.join(res, "sizes.json")) as fh:
        sz = json.load(fh)
    with open(os.path.join(res, "access.json")) as fh:
        acc = json.load(fh)
    r, s = sz["reference"], sz["sample"]
    bg = r["bgzipped_fasta"]
    ref_rows = [
        ("FASTA", r["fasta_uncompressed_bytes"], False),
        ("bgzipped FASTA\n+ .fai + .gzi", bg["bytes"] + bg["fai_bytes"] + bg["gzi_bytes"], False),
        ("Ensembl's FASTA.gz", r["ensembl_fasta_gz_bytes"], False),
        ("gpxil", r["gpxil"]["bytes"], True),
        (".2bit", r["twobit_bytes"], False),
    ]
    sample_rows = [
        ("BCF + .csi", s["bcf"]["bytes"] + s["bcf"]["csi_bytes"], False),
        ("VCF.gz + .tbi\n+ regions' BED", s["vcf_gz_bytes"] + s["tbi_bytes"] + s["bed_bytes"], False),
        ("VCF.gz", s["vcf_gz_bytes"], False),
        ("gpxil, without\nthe reference", s["gpxil"]["bytes"], True),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 6.6), dpi=150)
    fig.patch.set_facecolor(SURFACE)
    bars(axes[0, 0], ref_rows, "A. GRCh38, 3.1 billion bases", "MB")
    bars(axes[0, 1], sample_rows, "B. HG002's 5.9 million calls (GIAB v5.0q)", "MB")
    loc = acc["local"]
    reads(
        axes[1, 0],
        loc,
        [("reference: gpxil", "gpxil", GPXIL, "-"), ("reference: bgzipped FASTA (faidx)", "FASTA (faidx)", OTHER, "-")],
        "C. Reading the reference",
    )
    reads(
        axes[1, 1],
        loc,
        [
            ("calls: gpxil planes and insertions", "gpxil planes", GPXIL, "-"),
            ("calls: VCF records with genotypes (tabix)", "VCF records (tabix)", OTHER, "-"),
            ("haplotypes: gpxil, both as sequences", "gpxil, both haplotypes\nas sequences", GPXIL, (0, (3, 2))),
        ],
        "D. Reading HG002's calls",
    )
    fig.subplots_adjust(left=0.13, right=0.97, top=0.94, bottom=0.08, hspace=0.42, wspace=0.42)
    os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
    fig.savefig(os.path.join(HERE, "figures", "format.png"), facecolor=SURFACE)


if __name__ == "__main__":
    main()
