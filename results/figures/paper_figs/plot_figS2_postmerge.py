#!/usr/bin/env python3
"""figS2_postmerge — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import figstyle

figstyle.setup()

df_s2 = pd.read_csv("data/figS2_postmerge.csv")

fig = plt.figure(figsize=(15, 8))
# 2 rows, 3 columns
gs = fig.add_gridspec(2, 3)

families = ["knob180", "TR-1", "CentC"]

for i, family in enumerate(families):
    subset = df_s2[df_s2["family"] == family].sort_values("merge_gap_bp")
    
    # Coverage (Row 0)
    ax_cov = fig.add_subplot(gs[0, i])
    if i == 0:
        ax_cov.set_ylabel("Coverage (%)")
    ax_cov.set_title(f"{family}", fontweight="bold")
    
    # Boundary Offset (Row 1)
    ax_off = fig.add_subplot(gs[1, i])
    if i == 0:
        ax_off.set_ylabel("Mean boundary offset (bp)")
    ax_off.set_xlabel("Merge gap (bp)")
    
    for tool in subset["tool"].unique():
        tool_data = subset[subset["tool"] == tool]
        color = figstyle.TOOL_COLORS.get(tool, "black")
        
        ax_cov.plot(tool_data["merge_gap_bp"], tool_data["coverage_pct"], marker="o", 
                    color=color, label=tool, rasterized=True)
        ax_off.plot(tool_data["merge_gap_bp"], tool_data["mean_offset_bp"], marker="s", 
                    color=color, label=tool, rasterized=True)
        
        if tool == "BWTandem":
            # Annotate calls-per-array at a few gaps (0, 1000, 10000)
            for idx, row in tool_data.iterrows():
                if row["merge_gap_bp"] in [0, 1000, 10000]:
                    ax_cov.annotate(f"{row['calls_per_array']:.1f} calls", 
                                    (row["merge_gap_bp"], row["coverage_pct"]), 
                                    xytext=(5, 5), textcoords="offset points", fontsize=8)

    ax_cov.set_xscale("symlog", linthresh=100)
    ax_off.set_xscale("symlog", linthresh=100)
    ax_off.set_yscale("log")

caption = (
    r"$\mathbf{Figure\ S2.}$ "
    "Effects of coordinate-based merging on BWTandem's fragmented satellite calls."
)
plt.figtext(0.02, -0.06, caption, ha="left", fontsize=14, wrap=True)
plt.suptitle("Merging fragments buys coverage at the price of boundaries.", y=1.006, fontsize=14)

sns.despine()
plt.tight_layout()
figstyle.save(fig, "figS2_postmerge")
