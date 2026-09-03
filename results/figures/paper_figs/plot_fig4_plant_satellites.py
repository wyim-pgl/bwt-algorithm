#!/usr/bin/env python3
"""fig4_plant_satellites — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import figstyle
import numpy as np

figstyle.setup()

df_a = pd.read_csv("data/fig4a_colcen.csv")
df_b = pd.read_csv("data/fig4b_maize_coverage.csv")
df_c = pd.read_csv("data/fig4c_band_filter_delta.csv")

fig = plt.figure(figsize=(18, 5))
gs = fig.add_gridspec(1, 3)

# Panel A: Col-CEN (two aligned dot plots, one for coverage, one for recall)
# We can use a single panel with twin y-axes or two subplots. Let's use two subplots side-by-side or a single with different x.
# Actually, "two aligned dot plots" means tool on y-axis, and two x-axes (or two panels side by side).
ax_a = fig.add_subplot(gs[0])
ax_a.set_title("A. Arabidopsis Col-CEN", fontweight="bold")

# Assign markers based on mode
markers = {"de novo": "o", "array-level (2026)": "^", "interval-only (2026)": "s"}
df_a_sorted = df_a.sort_values("cen_coverage_pct", ascending=False)

for mode, m in markers.items():
    subset = df_a_sorted[df_a_sorted["mode"] == mode]
    if not subset.empty:
        ax_a.scatter(subset["cen_coverage_pct"], subset["tool"], marker=m, 
                     color=[figstyle.TOOL_COLORS.get(t, "black") for t in subset["tool"]],
                     label=mode, s=100, rasterized=True)
ax_a.set_xlabel("Centromere coverage (%)")
leg_a = ax_a.legend(title="Mode")
for handle in getattr(leg_a, 'legend_handles', getattr(leg_a, 'legendHandles', [])):
    handle.set_color('black')

# Add subtle separator lines between tools
num_tools_a = len(df_a_sorted["tool"].unique())
for y_pos in range(num_tools_a - 1):
    ax_a.axhline(y_pos + 0.5, color='gray', linestyle='--', alpha=0.3)

# Panel B: Maize unfiltered coverage
ax_b = fig.add_subplot(gs[1])
ax_b.set_title("B. Maize (unfiltered satellite coverage)", fontweight="bold")
sns.stripplot(data=df_b, x="coverage_pct", y="family", hue="tool", ax=ax_b,
              palette=figstyle.TOOL_COLORS, dodge=True, size=8, jitter=False, rasterized=True)
ax_b.set_xlabel("Unfiltered coverage (%)")
ax_b.set_ylabel("")
ax_b.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# Add subtle separator lines between families
num_families = len(df_b["family"].unique())
for y_pos in range(num_families - 1):
    ax_b.axhline(y_pos + 0.5, color='gray', linestyle='--', alpha=0.3)

# Panel C: Band filter delta (slope plot)
ax_c = fig.add_subplot(gs[2])
ax_c.set_title("C. Maize (coverage lost to period band)", fontweight="bold")

families = df_c["family"].unique()
tools = df_c["tool"].unique()

# x = 0 is unfiltered (0 delta), x = 1 is banded (delta)
ax_c.set_xlim(-0.2, 1.2)
ax_c.set_xticks([0, 1])
ax_c.set_xticklabels(["Unfiltered", "Banded"])
ax_c.set_ylabel("Coverage change (pp)")

for i, family in enumerate(families):
    # offset each family slightly so they don't overlap completely
    offset = (i - 1) * 0.1
    for tool in tools:
        subset = df_c[(df_c["family"] == family) & (df_c["tool"] == tool)]
        if not subset.empty:
            delta = subset.iloc[0]["coverage_delta_pp"]
            x_vals = [0 + offset, 1 + offset]
            y_vals = [0, delta]
            lw = 1.5
            zo = 2
            if tool == "tantan":
                lw = 3
                zo = 2
            elif tool == "TRF":
                lw = 1.5
                zo = 3
                
            ax_c.plot(x_vals, y_vals, marker="o", color=figstyle.TOOL_COLORS.get(tool, "black"), linewidth=lw, zorder=zo, rasterized=True)
            if tool == "BWTandem":
                ax_c.annotate(f"{delta:.1f}", (1 + offset, delta), xytext=(5, 0), textcoords="offset points", va="center", fontsize=8, bbox=dict(boxstyle="square,pad=0.1", facecolor="white", edgecolor="none", alpha=0.7))

ax_c.annotate("* 2026 tools (longdust, AniAnn's) lack a banded arm and are omitted here.",
              xy=(0, -0.2), xycoords='axes fraction', ha='left', fontsize=8)

caption = (
    r"$\mathbf{Figure\ 4.}$ "
    "(A) Centromere coverage and CEN180 monomer recall in Arabidopsis Col-CEN, for the six de novo and 2026 tools; "
    "Table 2 also lists mreps, NCRF and the two TRASH rows, which are omitted here. "
    "(B) Unfiltered coverage of curated maize satellite arrays. "
    "(C) Coverage lost when applying a period band."
)
plt.figtext(0.01, -0.1, caption, ha="left", fontsize=14, wrap=True)
plt.suptitle("Leading-group coverage on plant satellites, unstable period assignment under banding.", y=1.03, fontsize=14)

sns.despine()
plt.tight_layout(w_pad=0.2)

# Manually move Panel B to the left without affecting A and C
pos_b = ax_b.get_position()
ax_b.set_position([pos_b.x0 - 0.03, pos_b.y0, pos_b.width, pos_b.height])

figstyle.save(fig, "fig4_plant_satellites")
