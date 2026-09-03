"""Shared style for the paper figures (seaborn).

Interpreter with seaborn 0.13:
    /data/gpfs/assoc/pgl/bin/conda/conda_envs/bch709_vibe_coding/bin/python
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# One colour per tool, stable across every figure.
TOOL_COLORS = {
    "BWTandem": "#1f4e79",
    "ULTRA": "#c05621",
    "tantan": "#2f855a",
    "TRF": "#6b46c1",
    "TRASH": "#718096",
    "longdust": "#b7791f",
    "AniAnn's": "#d53f8c",
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
