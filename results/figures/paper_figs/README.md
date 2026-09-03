# Paper figures — delivered 2026-09-03

All six figures are implemented and rendered. `rendered/` holds the PNG and PDF
of each; the `plot_*.py` scripts that produced them are here, no longer stubs,
and read their inputs from `data/*.csv`.

| Figure | Renders | Status against the manuscript at `e232bcf` |
|---|---|---|
| Fig 1 accuracy trade-off | ✅ | current |
| Fig 2 range cost | ✅ | **three defects, see below** |
| Fig 3 sweep + audit | ✅ | current |
| Fig 4 plant satellites | ✅ | current |
| Fig S1 chromosome subset | ✅ | current |
| Fig S2 post-merge | ✅ | current |

## Fig 1 anticipated a decision we made the same day

Panel C plots BWTandem's **post-hoc ≤100 bp** arm, labelled as such in the legend
and stated in the caption, rather than the native `--max-period 100` rerun. That
is the arm Table 1b now uses, for the reason recorded as `quarantine.md` §6.18:
restricting BWTandem by re-running while restricting TRF by filtering is not a
shared range. The figure and the table agree, and they were decided
independently.

Panels A and B plot the P/B/F/H operating-point curve, which is Table 1d and
stays native — so F reads 79.88% there and 78.87% in Table 1b. That is not a
contradiction, but it is the kind of thing a referee asks about, and the caption
should say which panel is which arm before submission.

## Fig 2 — do not ship as rendered

1. **Panel A is superseded data.** The 1.30 / 1.30 / 1.41 ratios come from
   `data/fig2a_paired_runs.csv`, whose replicates all ran at commit `07ad6fa` —
   the only commit where Tier 3's search window was fixed at 100 bp–100 kb
   regardless of `--max-period` (`quarantine.md` §6.17). SLURM jobs
   6147698/99/700 are re-measuring the same pairs at `0363d8b`, where the
   ceiling bounds Tier 3 as well. **The panel, its title claim
   ("20× costs 1.3–1.4×") and the CSV must be regenerated from those results.**
   The x-axis phrase "requested/reportable maximum period" is the `07ad6fa`
   framing and needs the same revision.
2. **Axis label says "Peak memory (GB)".** Every memory figure in this paper is
   a kibibyte count divided by 1024², relabelled GiB throughout the manuscript in
   `63ccd07`. The axis and the `prep_data.py` comments still say GB.
3. **Footnote says "Competitor FASTA scope ~5% larger".** Measured 2026-09-03:
   3,209,286,105 bases against 3,088,269,832, i.e. **3.92%**, or 3.80% counting
   unambiguous bases only.

## Fig 4 panel A omits three tools without saying so

Table 2 has eleven rows; panel A plots six. TRASH, NCRF and mreps are absent and
the caption does not say they were dropped or why. Today's Table 2 correction —
the template row rebuilt as the true 397-region template-only union — therefore
does not touch this figure, but the omission should be stated.

## Regeneration

`prep_data.py` builds `data/*.csv` from the deposited evidence; each `plot_*.py`
reads one or more of those and writes to `rendered/`. Regenerate Fig 2 after the
C-1 jobs land, then re-check the three items above.
