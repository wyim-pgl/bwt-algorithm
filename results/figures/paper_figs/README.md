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
| Fig 5 array structure | ✅ | new 2026-09-03, not yet placed in the manuscript |

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

## Fig 2 — one defect left, and it is the data

Two of the three defects are fixed as of 2026-09-03: the memory axis reads GiB,
and the range-matching caveat moved into the caption with the measured 3.92%
(a separate footnote collided with the caption once the caption grew).

1. **Panel A is superseded data — still open.** The 1.30 / 1.30 / 1.41 ratios come from
   `data/fig2a_paired_runs.csv`, whose replicates all ran at commit `07ad6fa` —
   the only commit where Tier 3's search window was fixed at 100 bp–100 kb
   regardless of `--max-period` (`quarantine.md` §6.17). SLURM jobs
   6147698/99/700 are re-measuring the same pairs at `0363d8b`, where the
   ceiling bounds Tier 3 as well. **The panel, its title claim
   ("20× costs 1.3–1.4×") and the CSV must be regenerated from those results.**
   The x-axis phrase "requested/reportable maximum period" is the `07ad6fa`
   framing and needs the same revision.
2. ~~Axis said "Peak memory (GB)"~~ — now GiB.
3. ~~Footnote said "~5% larger"~~ — now 3.92%, in the caption.

The caption states which build panel A came from, so the figure is honest as it
stands; it is the claim in the title that the re-measurement may move.

## Fig 4 panel A omits three tools — now said in the caption

Table 2 has eleven rows; panel A plots six. TRASH, NCRF and mreps are absent, and
the caption now names them. Today's Table 2 correction —
the template row rebuilt as the true 397-region template-only union — therefore
does not touch this figure, but the omission should be stated.

## Regeneration

`prep_data.py` builds `data/*.csv` from the deposited evidence; each `plot_*.py`
reads one or more of those and writes to `rendered/`. Regenerate Fig 2 after the
C-1 jobs land, then re-check the three items above.


## Colour (2026-09-03)

`figstyle.py` gives BWTandem the only saturated colour and mutes every competitor,
so the eye lands on our tool first and the figures still read in grey scale. The
competitor hues stay distinct from each other — ULTRA and TRF had to be pulled
apart again after the first pass made them nearly identical in Fig 2B, where they
are the only two series. Do not raise a competitor's saturation to make it
visible; if one needs emphasis in a panel, say so in that panel's caption.

Rendered with `/data/gpfs/assoc/pgl/bin/conda/conda_envs/anianns/bin/python`
(matplotlib 3.10.9, seaborn 0.13.2). The `bwtandem` env has no matplotlib.


## Fig 5 — array self-similarity (new, 2026-09-03)

One satellite array per species, drawn two ways: the pairwise window-identity
matrix as a 45-degree-rotated triangle, and a self dot plot over the same
interval. Regions are BWTandem's own calls:

| Species | Region | Reported motif |
|---|---|--:|
| Human GRCh38 | chr1:122,257,803–122,317,803 (60 kb) | 1,866 bp |
| Arabidopsis Col-CEN | Chr4:4,985,644–5,045,644 (60 kb) | 178 bp |
| Maize Mo17 | CM039157.1:53,652,014–53,678,133 (26 kb) | 155 bp |

Windows are compared by shared canonical 13-mers (containment against the smaller
window) rather than by alignment, so one indel costs a few k-mers instead of
shifting everything downstream. Containment, not Jaccard, so a window straddling
an array edge is not penalised for being half empty.

What it shows: the human array carries regular bands at a spacing consistent with
its 1,866 bp reported motif — the multi-monomer higher-order repeat — while
CEN180 and CentC show finer, less regular structure. The dot plots carry the same
information as line ladders and keep the monomer period legible.

Two things to know before reusing it. The dot-plot sampling stride is calibrated
per panel to a fixed point budget, so the strides differ between rows and are
printed in each title; do not compare point density across panels. And the
identity colour scale is fixed at 0–1 across all three, which is why the human
panel reads darker — its 13-mer containment between distant windows really is
lower, not a scaling artefact.

It is not referenced by the manuscript yet.
