"""Shared style for the paper figures (seaborn).

Interpreter with seaborn 0.13:
    /data/gpfs/assoc/pgl/bin/conda/conda_envs/bch709_vibe_coding/bin/python
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# One colour per tool, stable across every figure.
#
# BWTandem carries the only saturated colour; every competitor is a muted slate
# so the eye lands on our tool first and the comparison stays readable in grey
# scale. The hues still differ between competitors so they remain separable, but
# their chroma is low enough that none of them competes with BWTandem for
# attention. Do not raise a competitor's saturation to "make it visible" — if a
# competitor needs emphasis in one panel, say so in that panel's caption.
BWT = "#0B5FA5"          # the accent: only BWTandem gets it
TOOL_COLORS = {
    "BWTandem": BWT,
    "ULTRA": "#B5794A",
    "tantan": "#6E8B7A",
    "TRF": "#8C7BA8",
    "TRASH": "#A9AFB8",
    "longdust": "#A79274",
    "AniAnn's": "#B08699",
}
BWT_POINT_ORDER = ["P", "B", "F", "H"]


def setup():
    sns.set_theme(style="ticks", context="paper", font_scale=1.05)
    plt.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,  # editable text in Illustrator
        "mathtext.fontset": "custom",
        "mathtext.bf": "sans:bold",
    })


def save(fig, stem):
    for ext in ("png", "pdf"):
        fig.savefig(f"{stem}.{ext}", bbox_inches="tight")
    print(f"wrote {stem}.png/.pdf")
