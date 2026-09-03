#!/usr/bin/env python3
"""fig1_accuracy_tradeoff — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import figstyle

figstyle.setup()

df_ab = pd.read_csv("data/fig1ab_pr_points.csv")
df_c = pd.read_csv("data/fig1c_overlap_rules.csv")

fig, axes = plt.subplots(1, 3, figsize=(16, 5))

def get_color(label):
    base_tool = label.split()[0]
    return figstyle.TOOL_COLORS.get(base_tool, figstyle.TOOL_COLORS.get(label, "#000000"))

# Panel A: Region-level precision vs recall
ax = axes[0]
ax.set_title("A. Region-level", fontweight="bold")
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_xlabel("Region recall (%)")
ax.set_ylabel("Region precision (%)")

bwt = df_ab[df_ab["series"] == "BWTandem"].set_index("label")
bwt = bwt.loc[figstyle.BWT_POINT_ORDER]
ax.plot(bwt["region_recall"], bwt["region_precision"], marker="o", 
        color=figstyle.TOOL_COLORS["BWTandem"], label="BWTandem", rasterized=True)
for idx, row in bwt.iterrows():
    ax.annotate(idx, (row["region_recall"], row["region_precision"]), 
                xytext=(5, 5), textcoords="offset points")

comps = df_ab[df_ab["series"].str.startswith("competitor")]
for idx, row in comps.iterrows():
    color = get_color(row["label"])
    ax.scatter(row["region_recall"], row["region_precision"], color=color, 
               label=row["label"], rasterized=True)
    
    anno_text = row["label"]
    if row["label"] == "AniAnn's":
        anno_text += " (window ≥1 kb)"
        
    xytext = (5, 5)
    ha = "left"
    if "longdust" in row["label"]:
        if "-k8" in row["label"]:
            # Panel A: longdust -k8w20000
            # Go up and right so it avoids the axis and goes above P
            xytext = (5, 25)
            ha = "left"
        else:
            # longdust is lower and to the right
            xytext = (5, -15)
            ha = "left"
            
    ax.annotate(anno_text, (row["region_recall"], row["region_precision"]),
                xytext=xytext, textcoords="offset points", ha=ha)

# Panel B: Base-pair precision vs recall
ax = axes[1]
ax.set_title("B. Base-pair", fontweight="bold")
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.set_xlabel("Base-pair recall (%)")
ax.set_ylabel("Base-pair precision (%)")

ax.plot(bwt["bp_recall"], bwt["bp_precision"], marker="o", 
        color=figstyle.TOOL_COLORS["BWTandem"], label="BWTandem", rasterized=True)
for idx, row in bwt.iterrows():
    ax.annotate(idx, (row["bp_recall"], row["bp_precision"]), 
                xytext=(5, 5), textcoords="offset points")

for idx, row in comps.iterrows():
    color = get_color(row["label"])
    ax.scatter(row["bp_recall"], row["bp_precision"], color=color, 
               label=row["label"], rasterized=True)
    anno_text = row["label"]
    if row["label"] == "AniAnn's":
        anno_text += " (window ≥1 kb)"
        
    xytext = (5, 5)
    ha = "left"
    if "longdust" in row["label"]:
        if "-k8" in row["label"]:
            # Panel B: longdust -k8w20000 is lower and to the left of longdust
            # Go right and down to empty space
            xytext = (5, -15)
            ha = "left"
        else:
            # longdust is higher and to the right. Go left and up.
            xytext = (-5, 5)
            ha = "right"
            
    ax.annotate(anno_text, (row["bp_recall"], row["bp_precision"]),
                xytext=xytext, textcoords="offset points", ha=ha)

# Panel C: Matched-range region recall vs overlap rule
ax = axes[2]
ax.set_title("C. Matched-range region recall", fontweight="bold")
ax.set_xlabel("Overlap rule")
ax.set_ylabel("Region recall (%)")
ax.set_ylim(0, 100)

rules_order = ["one-base", "reciprocal 0.25", "reciprocal 0.50"]
df_c["rule"] = pd.Categorical(df_c["rule"], categories=rules_order, ordered=True)
df_c = df_c.sort_values("rule")

for tool in df_c["tool"].unique():
    tool_data = df_c[df_c["tool"] == tool]
    color = get_color(tool)
    
    alpha = 1.0
    lw = 2.0
    ls = "-"
    if tool in ["TRF", "TRASH"]:
        alpha = 0.3
        lw = 1.5
    
    if tool == "BWTandem":
        ls = "--" 
        
    ax.plot(tool_data["rule"], tool_data["regRecall"], marker="o",
            color=color, alpha=alpha, linewidth=lw, linestyle=ls, rasterized=True)

# Custom legend for the whole figure
legend_elements = [
    mlines.Line2D([0], [0], color=figstyle.TOOL_COLORS["BWTandem"], marker='o', label='BWTandem (Panels A/B)'),
    mlines.Line2D([0], [0], color=figstyle.TOOL_COLORS["BWTandem"], marker='o', linestyle='--', label='BWTandem (post-hoc ≤100 bp)'),
]
for tool in df_c["tool"].unique():
    if tool != "BWTandem":
        color = get_color(tool)
        alpha = 0.3 if tool in ["TRF", "TRASH"] else 1.0
        legend_elements.append(mlines.Line2D([0], [0], color=color, marker='o', linestyle='', alpha=alpha, label=tool))

ax.legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left")

caption = (
    r"$\mathbf{Figure\ 1.}$ "
    "(A) Region-level and (B) base-pair precision–recall on GRCh38 against the adotto catalog (one-base overlap). BWTandem's four configurations "
    "(P, B, F, H) are connected lines; competing tools are single points. (C) Region recall evaluated only for periods ≤100 bp across three overlap rules. "
    "ULTRA performs best, with BWTandem second. BWTandem data in (C) uses full-range output filtered post-hoc to ≤100 bp to ensure accurate boundary "
    "assessment. Reciprocal rules here stress-test boundary precision against a merged catalog rather than absolute ranking. The 2026 tools are excluded "
    "from (C) as they lack period bounds. Panels A and B plot the four native --max-period 100 operating points (Table 1d), where F reads "
    "79.88% region recall; panel C plots the post-hoc filtered arm, which is what Table 1b reports, at 78.87%."
)
plt.figtext(0.06, -0.22, caption, ha="left", fontsize=14, wrap=True)

plt.figtext(0.55, 1.05, "Shared-range accuracy: BWTandem is second to ULTRA under every overlap rule.", ha="center", fontsize=14)
sns.despine()
figstyle.save(fig, "fig1_accuracy_tradeoff")
