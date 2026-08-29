#!/usr/bin/env python3
"""SUPERSEDED earlier-build recall/precision figure, retained for audit only.

Two panels, region-level and base-pair-level, both at the matched period range
(<=100 bp) so every point sits on one axis. BWTandem is a connected curve through
its four operating points; the competitors are single points, because each offers
one.

Numbers are transcribed from logs/score_wp1_5982266.out and re-emitted as CSV
alongside the figure so the figure and Table 1d cannot drift apart.

This script is not a manuscript input and deliberately writes only filenames
containing ``superseded``.  Figure 1 remains pending job 6129408.

Usage: make_curve_figure_superseded.py [--outdir DIR]
"""
import argparse
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# (label, region recall, region precision, bp recall, bp precision)
CURVE = [
    ("P", 54.69, 62.17, 26.81, 62.56),
    ("B", 71.90, 53.27, 31.14, 34.77),
    ("F", 79.88, 50.62, 32.72, 34.50),
    ("H", 81.60, 48.43, 33.57, 34.17),
]
POINTS = [
    ("ULTRA",  81.62, 53.66, 38.14, 61.33),
    ("tantan", 78.00, 57.66, 23.54, 70.22),
    ("TRF",    30.80, 94.73, 20.81, 93.16),
    ("TRASH",   0.91, 96.98,  5.01, 54.49),
]

BW = "#1f4e79"     # BWTandem curve
OTH = "#7f7f7f"    # competitors

# Label offsets in points, per panel. Several points sit almost on top of one
# another — B/F/H differ by 0.6 pp of bp precision, and ULTRA is 0.02 pp from H in
# region recall — so the offsets are hand-placed rather than uniform.
OFFSETS = {
    "region": {"P": (8, 6), "B": (8, 6), "F": (-6, 10), "H": (8, -4),
               "ULTRA": (10, 2), "tantan": (9, 4), "TRF": (10, -2),
               "TRASH": (10, -2)},
    "bp": {"P": (8, 6), "B": (-4, -14), "F": (6, -15), "H": (11, -6),
           "ULTRA": (-14, 10), "tantan": (9, 3), "TRF": (10, -2),
           "TRASH": (10, -2)},
}


def panel(ax, xi, yi, xlabel, ylabel, title, key):
    off = OFFSETS[key]
    xs = [r[xi] for r in CURVE]
    ys = [r[yi] for r in CURVE]
    ax.plot(xs, ys, "-o", color=BW, lw=2, ms=7, zorder=3,
            label="BWTandem (operating-point curve)")
    for lab, *v in CURVE:
        ax.annotate(lab, (v[xi - 1], v[yi - 1]), textcoords="offset points",
                    xytext=off[lab], fontsize=9, color=BW, fontweight="bold",
                    zorder=4)

    for lab, *v in POINTS:
        ax.plot(v[xi - 1], v[yi - 1], "s", color=OTH, ms=7, zorder=2)
        ax.annotate(lab, (v[xi - 1], v[yi - 1]), textcoords="offset points",
                    xytext=off[lab], fontsize=9, color="#333333", zorder=4)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, loc="left")
    ax.grid(alpha=0.25, ls=":")
    ax.set_axisbelow(True)
    # room for the labels, which otherwise run off the right edge (ULTRA) and the
    # top (TRASH)
    ax.margins(x=0.20, y=0.14)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=os.path.dirname(os.path.abspath(__file__)))
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    panel(axes[0], 1, 2, "Region recall (%)", "Region precision (%)",
          "a  Region level", "region")
    panel(axes[1], 3, 4, "Base pair recall (%)", "Base pair precision (%)",
          "b  Base pair level", "bp")
    axes[0].legend(loc="lower left", fontsize=9, frameon=False)
    fig.suptitle("Detection performance on human GRCh38 (adotto), periods ≤100 bp",
                 fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))

    for ext in ("pdf", "png"):
        p = os.path.join(a.outdir, f"figure_curve_superseded.{ext}")
        fig.savefig(p, dpi=300)
        print(f"wrote {p}")

    csv_path = os.path.join(a.outdir, "figure_curve_data_superseded.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["series", "label", "region_recall", "region_precision",
                    "bp_recall", "bp_precision"])
        for lab, rr, rp, br, bp in CURVE:
            w.writerow(["BWTandem", lab, rr, rp, br, bp])
        for lab, rr, rp, br, bp in POINTS:
            w.writerow(["competitor", lab, rr, rp, br, bp])
    print(f"wrote {csv_path}")


if __name__ == "__main__":
    main()
