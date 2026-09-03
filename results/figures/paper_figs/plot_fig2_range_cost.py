#!/usr/bin/env python3
"""fig2_range_cost — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import figstyle

figstyle.setup()

df_a = pd.read_csv("data/fig2a_paired_runs.csv")
df_b = pd.read_csv("data/fig2b_maize_scaling.csv")
df_c = pd.read_csv("data/fig2c_human_cost.csv")

fig = plt.figure(figsize=(18, 5))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 1])

# Panel A: Paired slope plot
ax_a = fig.add_subplot(gs[0])
ax_a.set_title("A. Human (paired runs)", fontweight="bold")
ax_a.set_ylabel("Runtime (h)")
ax_a.set_xlabel("Requested/reportable\nmaximum period (bp)")
ax_a.set_xticks([100, 2000])

for rep in df_a["replicate"].unique():
    subset = df_a[df_a["replicate"] == rep].sort_values("max_period_bp")
    ax_a.plot(subset["max_period_bp"], subset["runtime_h"], marker="o", 
              color=figstyle.TOOL_COLORS["BWTandem"], rasterized=True)
    
ax_a.annotate("1.30×", (2000, 8.61), textcoords="offset points", xytext=(5, 5), va="bottom")
ax_a.annotate("1.30×", (2000, 8.576), textcoords="offset points", xytext=(5, -5), va="top")
ax_a.annotate("1.41×", (2000, 7.746), textcoords="offset points", xytext=(5, -2), va="center")
ax_a.set_xlim(50, 2500) # give room for annotations

# Panel B: Maize scaling
ax_b = fig.add_subplot(gs[1])
ax_b.set_title("B. Maize (scaling)", fontweight="bold")
ax_b.set_xlabel("Requested/reportable\nmaximum period (bp)")
ax_b.set_ylabel("Runtime (h)")
ax_b.set_xscale("log")
ax_b.set_yscale("log")

for tool in df_b["tool"].unique():
    subset = df_b[df_b["tool"] == tool].sort_values("max_period_bp")
    ax_b.plot(subset["max_period_bp"], subset["runtime_h"], marker="o", 
              color=figstyle.TOOL_COLORS.get(tool, "black"), label=tool, rasterized=True)
ax_b.legend()
ax_b.annotate("observed scaling on maize runs,\nnot a cross-tool mechanism comparison",
              xy=(0.05, 0.95), xycoords='axes fraction', va='top', fontsize=8)

# Panel C1: Core-hours
ax_c1 = fig.add_subplot(gs[2])
ax_c1.set_title("C. Human (cost: core-hours)", fontweight="bold")
ax_c1.set_xlabel("Core-hours")
sns.stripplot(data=df_c, x="core_hours", y="tool", hue="tool", legend=False, ax=ax_c1, 
              palette=figstyle.TOOL_COLORS, size=8, jitter=False, rasterized=True)
for idx, row in df_c.reset_index().iterrows():
    ha_val = "left" if row["tool"] in ["longdust", "AniAnn's"] else "center"
    ax_c1.annotate(row["max_period_bp"], (row["core_hours"], idx), 
                   xytext=(0, -10), textcoords="offset points", ha=ha_val, fontsize=8)
    
# Panel C2: Peak memory
ax_c2 = fig.add_subplot(gs[3], sharey=ax_c1)
ax_c2.set_title("Human (cost: peak memory)", fontweight="bold")
ax_c2.set_xlabel("Peak memory (GB)")
sns.stripplot(data=df_c, x="memory_gb", y="tool", hue="tool", legend=False, ax=ax_c2, 
              palette=figstyle.TOOL_COLORS, size=8, jitter=False, rasterized=True)
for idx, row in df_c.reset_index().iterrows():
    ha_val = "left" if row["tool"] in ["longdust", "AniAnn's"] else "center"
    ax_c2.annotate(row["max_period_bp"], (row["memory_gb"], idx), 
                   xytext=(0, -10), textcoords="offset points", ha=ha_val, fontsize=8)
ax_c2.set_ylabel("")
ax_c2.tick_params(labelleft=False)

# Footnote for C
ax_c1.annotate("*Runs are not range-matched. Competitor FASTA scope ~5% larger.",
               xy=(0, -0.3), xycoords='axes fraction', ha='left', fontsize=8)

caption = (
    r"$\mathbf{Figure\ 2.}$ "
    "(A) Computational cost of extending BWTandem's maximum period from 100 bp to 2,000 bp on GRCh38 across 3 replicates. "
    "(B) Runtime scaling of ULTRA and TRF on maize Mo17 as the maximum period increases (log–log scale). "
    "(C) Total computational resources (core-hours and peak memory) for human whole-genome processing. Points are labeled by their period cap."
)
plt.figtext(0.009, -0.03, caption, ha="left", fontsize=14, wrap=True)
plt.suptitle("Widening the requested period range 20× costs 1.3–1.4×.", y=1.01, fontsize=14)

sns.despine()
plt.tight_layout()
figstyle.save(fig, "fig2_range_cost")
