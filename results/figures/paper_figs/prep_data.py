#!/usr/bin/env python3
"""Extract every CSV the paper figures need from the deposited evidence.

Deterministic: re-running overwrites data/*.csv with identical content as long
as the deposited inputs are unchanged. Sources are the repository's own
results/ files wherever one exists; the few values that live only in
manuscript tables (maize competitor cost points, Col-CEN Table 2) are
hard-coded here with their manuscript line provenance in comments, per
FIGURE_PLAN.md.
"""
import csv
import json
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
REGEN = os.path.join(REPO, "results", "regen")
DATA = os.path.join(HERE, "data")
os.makedirs(DATA, exist_ok=True)


def write(name, header, rows):
    path = os.path.join(DATA, name)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"wrote {name}: {len(rows)} rows")


def matched_range(path):
    """Parse the MATCHED RANGE block of a score_table1 report."""
    out = {}
    with open(path) as fh:
        block = False
        for line in fh:
            if line.startswith("########## MATCHED RANGE"):
                block = True
                continue
            if block:
                if line.startswith("##########") or not line.strip():
                    if out:
                        break
                    continue
                parts = line.split()
                if parts[0] == "tool":
                    continue
                out[parts[0]] = dict(
                    calls=int(parts[1]), regRecall=float(parts[2]),
                    regPrec=float(parts[3]), bpRecall=float(parts[4]),
                    bpPrec=float(parts[5]))
    return out


def baseline(path):
    """Parse the BASELINE block of a score_table1 report."""
    out = {}
    with open(path) as fh:
        block = False
        for line in fh:
            if line.startswith("########## BASELINE"):
                block = True
                continue
            if block:
                if line.startswith("##########"):
                    break
                parts = line.split()
                if not parts or parts[0] == "tool":
                    continue
                out[parts[0]] = dict(
                    calls=int(parts[1]), regRecall=float(parts[2]),
                    regPrec=float(parts[3]), bpRecall=float(parts[4]),
                    bpPrec=float(parts[5]))
    return out


# ---- Fig 1C: matched-range region recall under three overlap rules --------
# NOTE (FIGURE_PLAN.md): the BWTandem series here is the regenerated
# full-range output post-hoc filtered to <=100 bp, NOT the native H run.
rules = [("one-base", "recip_none.txt"),
         ("reciprocal 0.25", "recip_0.25.txt"),
         ("reciprocal 0.50", "recip_0.50.txt")]
rows = []
for rule, fname in rules:
    mr = matched_range(os.path.join(REGEN, fname))
    for tool, v in mr.items():
        rows.append([rule, tool, v["regRecall"], v["regPrec"]])
write("fig1c_overlap_rules.csv", ["rule", "tool", "regRecall", "regPrec"], rows)

# ---- Fig 2A: BWTandem paired range-cost replicates ------------------------
rows = []
with open(os.path.join(REPO, "results", "manifest.tsv")) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        if r["table"] != "range-rep":
            continue
        m = re.match(r"p(\d+)-r(\d)", r["row"])
        h, mn, s = (int(x) for x in r["elapsed"].split(":"))
        rows.append([f"r{m.group(2)}", int(m.group(1)),
                     round(h + mn / 60 + s / 3600, 3)])
write("fig2a_paired_runs.csv", ["replicate", "max_period_bp", "runtime_h"], rows)

# ---- Fig 2B: maize competitor period-vs-runtime points --------------------
# Manuscript Tables 3A/3C/3B rows (lines ~202-203, 232-233, 277-278):
# observed points, one genome (maize), not a scaling law.
rows = [
    ["ULTRA", 6, 0.23], ["ULTRA", 200, 12.9], ["ULTRA", 500, 34.7],
    ["TRF", 6, 5.2], ["TRF", 200, 5.5], ["TRF", 500, 5.5],
]
write("fig2b_maize_scaling.csv", ["tool", "max_period_bp", "runtime_h"], rows)

# ---- Fig 2C: human whole-genome cost (not range-matched) ------------------
# Table 1a / abstract: BWTandem 12.6 h x 2 thr = 25.3 core-h, 28.08 GB sacct;
# ULTRA 29.8 h x 2 = 59.6 core-h, 1.68 GB; TRF 33.7 h x 1, 1.45 GB (GNU time).
rows = [
    ["BWTandem", 2000, 12.6, 25.3, 28.08],
    ["ULTRA", 100, 29.8, 59.6, 1.68],
    ["TRF", 500, 33.7, 33.7, 1.45],
]
write("fig2c_human_cost.csv",
      ["tool", "max_period_bp", "wall_h", "core_hours", "memory_gb"], rows)

# ---- Fig 3A: identity sweep (full range) ----------------------------------
bl = baseline(os.path.join(REGEN, "score_table1_idsweep.txt"))
arm_order = ["BWT-id-off", "BWT-id-0.80", "BWT-id-0.76",
             "BWT-id-0.72", "BWT-id-0.68"]
labels = ["off", "0.80", "0.76", "0.72", "0.68"]
rows = [[lab, bl[a]["calls"], bl[a]["regRecall"], bl[a]["regPrec"]]
        for lab, a in zip(labels, arm_order)]
write("fig3a_idsweep.csv",
      ["catchall_identity", "calls", "regRecall", "regPrec"], rows)

# ---- Fig 3B: audit verdicts by stratum ------------------------------------
rows = [
    ["1-6", 3, 59, 38], ["7-20", 1, 89, 10],
    ["21-100", 0, 99, 1], ["101-2000", 0, 99, 1],
]  # results/audit11/aggregate_reviewer2_20260831.txt
write("fig3b_audit.csv",
      ["period_stratum_bp", "SUPPORTED", "UNSUPPORTED", "UNSURE"], rows)

# ---- Fig 4A: Col-CEN coverage and monomer recall (manuscript Table 2) -----
rows = [
    ["BWTandem", 84.54, 99.72, "de novo"],
    ["ULTRA", 84.44, 99.80, "de novo"],
    ["TRF", 84.39, 97.55, "de novo"],
    ["tantan", 81.50, 99.24, "de novo"],
]
write("fig4a_colcen.csv",
      ["tool", "cen_coverage_pct", "cen180_monomer_recall_pct", "mode"], rows)

# ---- Fig 4B/4C: maize unfiltered coverage + band-filter deltas ------------
# Values from results/regen/table3bc_replacement.md (regenerated Table 3B/3C
# rescoring, hard-validated against the old-BED manuscript cells first).
fallback_cov = [
    ["BWTandem", "knob180", 79.79], ["BWTandem", "TR-1", 50.34], ["BWTandem", "CentC", 58.55],
    ["ULTRA", "knob180", 78.71], ["ULTRA", "TR-1", 36.36], ["ULTRA", "CentC", 57.73],
    ["TRF", "knob180", 80.01], ["TRF", "TR-1", 45.77], ["TRF", "CentC", 58.50],
    ["tantan", "knob180", 71.92], ["tantan", "TR-1", 22.74], ["tantan", "CentC", 56.51],
]
write("fig4b_maize_coverage.csv", ["tool", "family", "coverage_pct"], fallback_cov)
fallback_delta = [
    ["BWTandem", "knob180", -15.59], ["BWTandem", "TR-1", -33.01], ["BWTandem", "CentC", -17.84],
    ["ULTRA", "knob180", -1.20], ["ULTRA", "TR-1", -2.82], ["ULTRA", "CentC", -1.89],
    ["TRF", "knob180", -0.63], ["TRF", "TR-1", -1.54], ["TRF", "CentC", -0.87],
    ["tantan", "knob180", -0.51], ["tantan", "TR-1", -1.38], ["tantan", "CentC", -0.91],
]
write("fig4c_band_filter_delta.csv",
      ["tool", "family", "coverage_delta_pp"], fallback_delta)

# ---- Fig S1: chromosome-subset sensitivity --------------------------------
panels = [("heldout_22chrom", "22 non-selection chr"),
          ("heldout_select_chr21_22", "chr21+22"),
          ("heldout_all24", "all 24 chr")]
rows = []
for fname, label in panels:
    mr = matched_range(os.path.join(REGEN, f"{fname}.txt"))
    for tool in ("BWTandem", "ULTRA", "tantan"):
        rows.append([label, tool, mr[tool]["regRecall"], mr[tool]["regPrec"]])
write("figS1_subset_sensitivity.csv",
      ["panel", "tool", "regRecall", "regPrec"], rows)

# ---- Fig S2: maize coordinate post-merge trade-off ------------------------
extra = json.load(open(os.path.join(REGEN, "maize_extra_evidence.json")))
pm = extra.get("coordinate_postmerge") or extra.get("postmerge") or {}
rows = []


def walk(node, ctx):
    if isinstance(node, dict):
        if {"coverage_percent", "mean_offset_bp"} <= set(node):
            rows.append([ctx.get("family"), ctx.get("tool"), ctx.get("gap"),
                         node["coverage_percent"], node["mean_offset_bp"],
                         node.get("calls_per_array")])
        else:
            for k, v in node.items():
                nctx = dict(ctx)
                ks = str(k)
                if ks.lstrip("-").isdigit():
                    nctx["gap"] = int(ks)
                elif ks in ("knob180", "TR-1", "TR1", "CentC"):
                    nctx["family"] = ks
                else:
                    nctx.setdefault("tool", ks)
                walk(v, nctx)


walk(pm, {})
if rows:
    write("figS2_postmerge.csv",
          ["family", "tool", "merge_gap_bp", "coverage_pct",
           "mean_offset_bp", "calls_per_array"], rows)
else:
    print("WARNING: coordinate_postmerge not found in maize_extra_evidence.json"
          " — inspect its schema and adapt walk() (keys: %s)" % sorted(extra)[:10])
