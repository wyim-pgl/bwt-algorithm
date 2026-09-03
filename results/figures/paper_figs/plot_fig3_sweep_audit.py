#!/usr/bin/env python3
"""fig3_sweep_audit — see FIGURE_PLAN.md for the panel spec, data provenance and the
caption draft. Data: data/*.csv (regenerate with prep_data.py)."""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import figstyle

figstyle.setup()

df_a = pd.read_csv("data/fig3a_idsweep.csv")
df_b = pd.read_csv("data/fig3b_audit.csv")

fig = plt.figure(figsize=(12, 5))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2])

# Panel A: Identity sweep
ax_a = fig.add_subplot(gs[0])
ax_a.set_title("A. Identity sweep", fontweight="bold")

x_order = ["off", "0.80", "0.76", "0.72", "0.68"]
df_a["catchall_identity"] = pd.Categorical(df_a["catchall_identity"], categories=x_order, ordered=True)
df_a = df_a.sort_values("catchall_identity")

ax_a.plot(df_a["catchall_identity"], df_a["regRecall"], marker="o", 
          label="Region recall", color=figstyle.TOOL_COLORS["BWTandem"], rasterized=True)
ax_a.plot(df_a["catchall_identity"], df_a["regPrec"], marker="s", linestyle="--",
          label="Region precision", color="dimgrey", rasterized=True)

ax_a.set_xlabel("Catch-all identity")
ax_a.set_ylabel("Metric (%)")
ax_a.set_ylim(0, 100)
ax_a.legend()

# Panel B: Audit verdicts
ax_b = fig.add_subplot(gs[1])
ax_b.set_title("B. Audit of 400 BWTandem-only calls", fontweight="bold")

# Prepare data for stacked bar
periods = df_b["period_stratum_bp"].tolist()
supported = df_b["SUPPORTED"].values
unsupported = df_b["UNSUPPORTED"].values
unsure = df_b["UNSURE"].values

# Calculate overall
tot_sup = supported.sum()
tot_unsup = unsupported.sum()
tot_unsure = unsure.sum()

categories = periods + ["Overall"]
val_sup = list(supported) + [tot_sup]
val_unsup = list(unsupported) + [tot_unsup]
val_unsure = list(unsure) + [tot_unsure]

colors = {"SUPPORTED": "#2ca02c", "UNSUPPORTED": "#d62728", "UNSURE": "#7f7f7f"}

x = range(len(categories))
# Make overall bar stand out a bit
bar_width = [0.6]*len(periods) + [0.8]
x_pos = list(range(len(periods))) + [len(periods) + 0.5]

p1 = ax_b.bar(x_pos, val_sup, width=bar_width, color=colors["SUPPORTED"], label="SUPPORTED", rasterized=True)
p2 = ax_b.bar(x_pos, val_unsup, bottom=val_sup, width=bar_width, color=colors["UNSUPPORTED"], label="UNSUPPORTED", rasterized=True)
bottom2 = [s + u for s, u in zip(val_sup, val_unsup)]
p3 = ax_b.bar(x_pos, val_unsure, bottom=bottom2, width=bar_width, color=colors["UNSURE"], label="UNSURE", rasterized=True)

ax_b.set_xticks(x_pos)
ax_b.set_xticklabels(categories, rotation=45, ha="right")
ax_b.set_ylabel("Count")

# Add text annotations on overall bar
for i, b in enumerate([val_sup[-1], val_unsup[-1], val_unsure[-1]]):
    if b > 0:
        y_val = sum([val_sup[-1], val_unsup[-1], val_unsure[-1]][:i]) + b/2
        ax_b.text(x_pos[-1], y_val, str(b), ha="center", va="center", color="white", fontweight="bold")

ax_b.legend(bbox_to_anchor=(1.05, 1), loc="upper left")

caption = (
    r"$\mathbf{Figure\ 3.}$ "
    "(A) Precision–recall trade-off for the catch-all identity sweep (Supplementary Table S3). "
    "(B) Blinded manual audit of 400 BWTandem-exclusive calls (100 per period stratum). "
    "Results show 4 supported, 346 unsupported, and 50 unsure calls (Wilson 95% CI: 0.4–2.9% for definitive verdicts)."
)
plt.figtext(0.02, -0.15, caption, ha="left", fontsize=14, wrap=True)
plt.suptitle("Recall-favouring settings buy recall with precision, and the unmatched calls are predominantly unsupported.", y=1.05, fontsize=12)

sns.despine()
plt.tight_layout()
figstyle.save(fig, "fig3_sweep_audit")
