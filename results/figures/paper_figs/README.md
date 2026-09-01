# Paper figures workspace (seaborn)

Working directory for the manuscript's display items. **Implementing the
plots is handed off — start from `HANDOFF_FILIP.md`** (English: per-figure
titles, panel specs, data-file mapping, caption drafts, ground rules). The
full adversarial design document behind it is `FIGURE_PLAN.md` (Korean,
2026-08-31), including the list of figures that must NOT be made.

## Layout

- `prep_data.py` — regenerates every `data/*.csv` from the deposited evidence
  (`results/regen/`, `results/manifest.tsv`, `results/audit11/`). The few
  values that exist only as manuscript-table cells (maize competitor cost
  points, Col-CEN Table 2) are hard-coded there with provenance comments.
- `data/` — one CSV per panel input:
  | file | feeds | source |
  |---|---|---|
  | fig1ab_pr_points.csv | Fig 1A/1B | results/figures/figure_curve_data.csv |
  | fig1c_overlap_rules.csv | Fig 1C | recip_{none,0.25,0.50}.txt MATCHED RANGE |
  | fig2a_paired_runs.csv | Fig 2A | manifest.tsv range-rep rows |
  | fig2b_maize_scaling.csv | Fig 2B | manuscript Tables 3A/3B/3C cost cells |
  | fig2c_human_cost.csv | Fig 2C | Table 1a cost cells |
  | fig3a_idsweep.csv | Fig 3A | score_table1_idsweep.txt BASELINE |
  | fig3b_audit.csv | Fig 3B | audit11/aggregate_reviewer2_20260831.txt |
  | fig4a_colcen.csv | Fig 4A | manuscript Table 2 |
  | fig4b_maize_coverage.csv | Fig 4B | table3bc_replacement.md |
  | fig4c_band_filter_delta.csv | Fig 4C | table3bc_replacement.md |
  | figS1_subset_sensitivity.csv | Fig S1 | heldout_*.txt MATCHED RANGE |
  | figS2_postmerge.csv | Fig S2 | maize_extra_evidence.json coordinate_postmerge |
- `figstyle.py` — shared seaborn theme; one stable colour per tool.
- `plot_fig*.py` — one script per figure, implemented in the order
  Fig 2 → 1 → 3 → 4 → S1 → S2 (App Note survivors first).

## Interpreter

seaborn 0.13 lives in:

```
/data/gpfs/assoc/pgl/bin/conda/conda_envs/bch709_vibe_coding/bin/python
```

## Non-negotiable honesty constraints (from FIGURE_PLAN.md)

- Fig 1C's BWTandem series is the full-range output post-hoc filtered to
  ≤100 bp, not the native H run — plot it as its own labelled series.
- Fig 2 axis label is "requested/reportable maximum period", never
  "searched range"; panels A/C are human, B is maize — no shared fit line.
- Never truncate a recall axis to manufacture the 81.60-vs-81.62 "tie".
- The audit panel shows raw stratified verdict counts; no population
  extrapolation, no pie charts.
