#!/usr/bin/env python3
"""figS1_subset — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import figstyle

figstyle.setup()

df_s1 = pd.read_csv("data/figS1_subset_sensitivity.csv")

fig = plt.figure(figsize=(12, 5))
gs = fig.add_gridspec(1, 2)

x_order = ["22 non-selection chr", "chr21+22", "all 24 chr"]
df_s1["panel"] = pd.Categorical(df_s1["panel"], categories=x_order, ordered=True)
df_s1 = df_s1.sort_values("panel")

# Panel A: Recall
ax_a = fig.add_subplot(gs[0])
ax_a.set_title("A. Matched-range region recall", fontweight="bold")
for tool in df_s1["tool"].unique():
    subset = df_s1[df_s1["tool"] == tool]
    ax_a.plot(subset["panel"], subset["regRecall"], marker="o", 
              color=figstyle.TOOL_COLORS.get(tool, "black"), label=tool, rasterized=True)
ax_a.set_ylabel("Region recall (%)")

# Panel B: Precision
ax_b = fig.add_subplot(gs[1])
ax_b.set_title("B. Matched-range region precision", fontweight="bold")
for tool in df_s1["tool"].unique():
    subset = df_s1[df_s1["tool"] == tool]
    ax_b.plot(subset["panel"], subset["regPrec"], marker="o", 
              color=figstyle.TOOL_COLORS.get(tool, "black"), label=tool, rasterized=True)
ax_b.set_ylabel("Region precision (%)")
ax_b.legend()

caption = (
    r"$\mathbf{Figure\ S1.}$ "
    "Matched-range region recall (A) and precision (B) on the 22 non-selection chromosomes, the 2 selection chromosomes, and all 24."
)
plt.figtext(0.02, -0.13, caption, ha="left", fontsize=14, wrap=True)
plt.suptitle("The human ranking is stable across chromosome subsets.", y=1.05, fontsize=14)

sns.despine()
plt.tight_layout()
figstyle.save(fig, "figS1_subset")
