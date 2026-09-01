# Figure handoff brief (for Filip)

Everything needed to draw the manuscript's display items is in this folder:
one CSV per panel under `data/` (regenerate any time with `prep_data.py` —
they are derived from the deposited evidence, never hand-typed), a shared
seaborn style in `figstyle.py` (one stable colour per tool — please keep
it), and stub scripts `plot_fig*.py` to fill in. The full adversarial
design document is `FIGURE_PLAN.md` (Korean; this brief is the English
operational summary and is self-sufficient). Interpreter with seaborn 0.13:
`/data/gpfs/assoc/pgl/bin/conda/conda_envs/bch709_vibe_coding/bin/python`.

**Output**: PNG + PDF, 300 dpi, editable text (figstyle sets
`pdf.fonttype=42`). Implementation order: **Fig 2 → 1 → 3 → 4 → S1 → S2**
(the two Application-Note survivors first).

**Ground rules (non-negotiable, from the adversarial review):**
- Never truncate a recall axis to manufacture the 81.60-vs-81.62 "tie".
- Fig 2's x-axis label is **"requested/reportable maximum period (bp)"**,
  never "searched range"; panels A/C are human, B is maize — no shared fit
  line or scaling exponent across panels.
- Fig 1C's BWTandem series is the full-range output post-hoc filtered to
  ≤100 bp, not the native H run — label it as its own series; do not reuse
  the P/B/F/H connecting line there.
- The audit panel (3B) shows raw stratified verdict counts only — no
  population extrapolation, no pie charts.
- See `FIGURE_PLAN.md` §3 for seven figure types we must NOT make.

---

## Figure 1 — Accuracy–precision trade-off on the shared search range
**Title:** *Shared-range accuracy: BWTandem is second to ULTRA under every
overlap rule.*
3 panels.

- **1A** Region-level precision (y) vs recall (x) scatter, human GRCh38 vs
  adotto. Data `data/fig1ab_pr_points.csv` (`region_recall`,
  `region_precision`). Connect BWTandem P→B→F→H with a line
  (`figstyle.BWT_POINT_ORDER`), competitors as labelled points. Axes 0–100,
  do not zoom.
- **1B** Same, base-pair metrics (`bp_recall`, `bp_precision`).
- **1C** Matched-range region recall vs overlap rule. Data
  `data/fig1c_overlap_rules.csv`; x categorical: one-base → reciprocal
  0.25 → reciprocal 0.50. Emphasise ULTRA/BWTandem/tantan; draw TRF/TRASH
  faded. Annotate that the BWTandem series is "regenerated full-range
  output, post-hoc ≤100 bp".
- **Caption draft:** "Region-level (A) and base-pair (B) precision–recall
  on GRCh38 against the adotto v1.2.1 catalog (one-base overlap rule,
  Table 1a pipeline). BWTandem's four operating points (P, B, F, H) are
  connected; competitors are single points. (C) Matched-range (periods
  ≤100 bp) region recall under three overlap rules. ULTRA ranks first and
  BWTandem second under every rule; tantan's near-tie with BWTandem at the
  one-base rule does not survive reciprocal overlap. The BWTandem series
  in (C) is the regenerated full-range output filtered post hoc to
  ≤100 bp, not the native period-100 runs of (A–B). Absolute values under
  the reciprocal rules are a stress test of boundary conventions against a
  merged, slop-padded catalog, not a ranking metric."

## Figure 2 — Range–cost behaviour and whole-genome compute cost
**Title:** *Widening the requested period range 20× costs 1.3–1.4×.*
3 panels (A/C human, B maize — say so in the panel titles).

- **2A** Paired slope plot, BWTandem human replicates. Data
  `data/fig2a_paired_runs.csv`: three replicate pairs 100→2,000 bp;
  annotate each pair's ratio (1.30×, 1.30×, 1.41×). y runtime (h),
  x requested/reportable maximum period.
- **2B** ULTRA and TRF runtime vs maximum period on maize, log–log. Data
  `data/fig2b_maize_scaling.csv`. Two connected series; annotate "observed
  scaling on maize runs, not a cross-tool mechanism comparison".
- **2C** Two small horizontal dot plots: core-hours and peak memory (GB)
  for the human runs, BWTandem/ULTRA/TRF only. Data
  `data/fig2c_human_cost.csv`. Label each point with its period cap
  (2,000 / 100 / 500 bp). Footnote: runs are not range-matched and the
  competitor FASTA scope was ~5% larger.
- **Caption draft:** "(A) Paired BWTandem runs on GRCh38: raising the
  requested/reportable maximum period from 100 to 2,000 bp costs
  1.30–1.41× in wall clock across three replicate pairs. (B) Observed
  runtime versus maximum period for ULTRA and TRF on maize Mo17 (log–log);
  the two tools scale differently, and no cross-tool mechanism comparison
  is implied. (C) Measured cost of the human whole-genome runs
  (core-hours; sacct peak memory for BWTandem, GNU-time for competitors),
  each point labelled with its period cap — the runs are not
  range-matched. BWTandem's 28.08 GB peak is the cost of the design."

## Figure 3 — Sensitivity levers and the audit of unsupported unique calls
**Title:** *Recall-favouring settings buy recall with precision, and the
unmatched calls are predominantly unsupported.*
2 panels.

- **3A** Identity sweep: region recall and region precision (two clearly
  distinct series) vs categorical x = off, 0.80, 0.76, 0.72, 0.68. Data
  `data/fig3a_idsweep.csv`.
- **3B** 100% stacked bars of audit verdicts (SUPPORTED / UNSUPPORTED /
  UNSURE) per period stratum, plus the overall 4/346/50 as a separate
  right-hand summary. Data `data/fig3b_audit.csv`. The Wilson 95% CI
  (0.4–2.9%) applies to the overall definitive-verdict rate only — do not
  attach it to strata.
- **Caption draft:** "(A) The catch-all identity sweep trades region
  precision for recall (Supplementary Table S3). (B) Blinded single-reader
  audit of 400 BWTandem-only calls (100 per reported-period stratum):
  4 supported / 346 unsupported / 50 unsure overall; the two populations
  in (A) and (B) are related limitations, not a causal chain."

## Figure 4 — Plant satellites: coverage, boundary and period assignment
**Title:** *Leading-group coverage on plant satellites, unstable period
assignment under banding.*
3 panels.

- **4A** Col-CEN: centromere coverage and CEN180 monomer recall as two
  aligned dot plots. Data `data/fig4a_colcen.csv`.
- **4B** Maize unfiltered coverage per family (knob180, TR-1, CentC),
  grouped dot plot. Data `data/fig4b_maize_coverage.csv`.
- **4C** Coverage change when the period band is applied, slope plot
  unfiltered→banded per family and tool. Data
  `data/fig4c_band_filter_delta.csv` (deltas in pp; before-values are 4B).
  Use the file's values, not manuscript roundings, throughout.
- **Caption draft:** "(A) Arabidopsis Col-CEN centromere coverage and
  CEN180 monomer recall. (B) Unfiltered coverage of the curated maize
  satellite arrays. (C) Coverage lost when a period band is imposed:
  BWTandem loses 15.6–33.0 pp where competitors lose ≤2.8 pp, because its
  satellite period assignments are less stable — the trade-off named in
  the limitations."

## Figure S1 — Chromosome-subset sensitivity analysis
**Title:** *The human ranking is stable across chromosome subsets.*
2 panels (A recall, B precision), x = 22 non-selection chr / chr21+22 /
all 24; lines BWTandem, ULTRA, tantan. Data
`data/figS1_subset_sensitivity.csv`. Use the term "chromosome-subset
sensitivity analysis", NOT "held-out validation".
**Caption draft:** "Matched-range region recall (A) and precision (B) on
the 22 chromosomes not used for configuration selection, the two selection
chromosomes, and all 24. Every tool moves the same way; the selection
chromosomes are slightly easier for all tools."

## Figure S2 — Coordinate post-merge trade-off on maize satellites
**Title:** *Merging fragments buys coverage at the price of boundaries.*
3 columns (knob180, TR-1, CentC) × 2 rows (coverage %; mean boundary
offset, log y). x = merge gap 0–10,000 bp (symlog). Data
`data/figS2_postmerge.csv`; annotate calls-per-array at a few gaps. No
dual y-axes.
**Caption draft:** "Coordinate-only merging of BWTandem's fragmented
satellite calls: coverage rises toward 100% with the merge gap while the
mean boundary offset degrades by orders of magnitude, so 'the whole array
was detected' would misstate what a merged interval represents."
