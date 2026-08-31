#!/usr/bin/env python3
"""Render active Figure 1, failing safely while its CSV is pending."""
import argparse
import csv
import os

HEADER = ["series", "label", "region_recall", "region_precision",
          "bp_recall", "bp_precision"]
CURVE_LABELS = ["P", "B", "F", "H"]
POINT_LABELS = ["ULTRA", "tantan", "TRF", "TRASH"]
BW, OTH = "#1f4e79", "#7f7f7f"
OFFSETS = {
    "region": {"P": (8, 6), "B": (8, 6), "F": (-6, 10), "H": (8, -4),
               "ULTRA": (10, 2), "tantan": (9, 4), "TRF": (10, -2),
               "TRASH": (10, -2)},
    "bp": {"P": (8, 6), "B": (-4, -14), "F": (6, -15), "H": (11, -6),
           "ULTRA": (-14, 10), "tantan": (9, 3), "TRF": (10, -2),
           "TRASH": (10, -2)},
}


def fail(message):
    raise SystemExit(f"Figure 1 not rendered: {message}")


def load_rows(path):
    if not os.path.isfile(path):
        fail(f"input CSV does not exist: {path}")
    with open(path, newline="") as handle:
        rows = list(csv.reader(handle))
    if not rows:
        fail(f"input CSV is empty: {path}")
    if rows[0] != HEADER:
        if rows[0] and rows[0][0].strip().lower() == "status":
            detail = ", ".join(rows[1]) if len(rows) > 1 else "pending sentinel"
            fail(f"input is pending ({detail})")
        fail(f"unexpected CSV header: {rows[0]!r}")
    parsed = []
    for line_no, row in enumerate(rows[1:], 2):
        if len(row) != 6:
            fail(f"line {line_no} has {len(row)} fields; expected 6")
        try:
            values = [float(value) for value in row[2:]]
        except ValueError:
            fail(f"line {line_no} contains a non-numeric metric")
        if any(value < 0 or value > 100 for value in values):
            fail(f"line {line_no} contains a metric outside 0--100")
        parsed.append((row[0], row[1], *values))
    curve = {row[1]: row for row in parsed if row[0] == "BWTandem"}
    points = {row[1]: row for row in parsed if row[0] == "competitor"}
    if any(row[0] not in {"BWTandem", "competitor"} for row in parsed):
        fail("unknown series name")
    if set(curve) != set(CURVE_LABELS) or set(points) != set(POINT_LABELS):
        fail(f"expected BWTandem {CURVE_LABELS} and competitors {POINT_LABELS}")
    if len(parsed) != 8:
        fail("duplicate or extra rows are not allowed")
    return [curve[x] for x in CURVE_LABELS], [points[x] for x in POINT_LABELS]


def panel(ax, curve, points, xi, yi, xlabel, ylabel, title, key):
    ax.plot([r[xi] for r in curve], [r[yi] for r in curve], "-o", color=BW,
            lw=2, ms=7, zorder=3, label="BWTandem (operating-point curve)")
    for _, label, *values in curve:
        ax.annotate(label, (values[xi - 2], values[yi - 2]),
                    textcoords="offset points", xytext=OFFSETS[key][label],
                    fontsize=9, color=BW, fontweight="bold", zorder=4)
    for _, label, *values in points:
        ax.plot(values[xi - 2], values[yi - 2], "s", color=OTH, ms=7, zorder=2)
        ax.annotate(label, (values[xi - 2], values[yi - 2]),
                    textcoords="offset points", xytext=OFFSETS[key][label],
                    fontsize=9, color="#333333", zorder=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=11, loc="left")
    ax.grid(alpha=0.25, ls=":")
    ax.set_axisbelow(True)
    ax.margins(x=0.20, y=0.14)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(here, "figure_curve_data.csv"))
    parser.add_argument("--outdir", default=here)
    args = parser.parse_args()
    curve, points = load_rows(os.path.abspath(args.input))
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:
        fail(f"matplotlib is required after the CSV is complete ({exc})")
    os.makedirs(args.outdir, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8))
    panel(axes[0], curve, points, 2, 3, "Region recall (%)",
          "Region precision (%)", "a  Region level", "region")
    panel(axes[1], curve, points, 4, 5, "Base pair recall (%)",
          "Base pair precision (%)", "b  Base pair level", "bp")
    axes[0].legend(loc="lower left", fontsize=9, frameon=False)
    fig.suptitle("Detection performance on human GRCh38 (adotto), periods <=100 bp",
                 fontsize=12, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    for extension in ("pdf", "png"):
        path = os.path.join(args.outdir, f"figure_curve.{extension}")
        fig.savefig(path, dpi=300)
        print(f"wrote {path}")
    plt.close(fig)


if __name__ == "__main__":
    main()
