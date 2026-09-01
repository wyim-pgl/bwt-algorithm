# Application Note conversion plan (2026-09-01)

> ⚠️ **Superseded as the default direction (same day)**: the venue analysis
> (docs/2026-09-01-venue-and-positioning-analysis.md) recommends a FULL
> paper at GigaScience or Bioinformatics — the contribution's center of
> mass is the provable benchmark, which ~2,600 App Note words cannot hold
> once longdust/AniAnn's enter the comparison. This plan is retained as a
> fallback design only.

Decision context: shared-range accuracy is honestly non-leading (ULTRA
first by ~1pp), so the publishable claim is **capability** — one pass over
periods 1 bp–2,000 bp at a 1.30–1.41× cost for a 20× wider range, with the
only tractable >500 bp coverage — not superiority. The Application Note
format sells exactly that; the full-paper track stays viable at venues that
value the benchmark methodology itself (venue analysis pending; this plan
is venue-agnostic and parameterized where formats differ).

## What survives in the main text (~1,100–1,300 words + 2 display items)

1. **One paragraph of problem + claim.** Existing period-capped detectors
   force per-range re-runs; the wide-range single pass is the contribution.
   State ULTRA-first plainly in the same breath as the 1pp gap and the 20×
   range — the abstract already does this; compress it.
2. **Figure 1 = current Fig 2 (range–cost).** Panels: paired 100→2,000 bp
   runs (1.30–1.41×); maize ULTRA/TRF scaling; human core-hours/memory
   dot plots. This is the value proposition. If the running ULTRA `-p 2000`
   measurement (job 6145581) lands before submission, its wall clock (or
   documented non-completion) becomes the killer annotation on panel C.
3. **Figure 2 = current Fig 1 (shared-range accuracy).** P/B/F/H curve +
   competitors + overlap-rule robustness. Keeps us honest and preempts the
   "you hid the ranking" review.
4. **One results paragraph** citing: full-range region recall 80.53% /
   held-out 80.48%; matched-range BWT 78.87 vs ULTRA 81.62; the >100 bp
   stratum (3.43% of catalog, structurally absent from range-capped
   tools); Col-CEN 99.72% monomer recall in 40 min; the specificity audit
   in one sentence (4/400 supported — unmatched calls are predominantly
   over-calls).
5. **Availability/repro paragraph**: pip/conda install, deposited BEDs,
   hashed manifests, per-cell provenance.

## What moves to the supplement (all of it already exists)

- Tables 1a–1e, 2, 3A–3C and their long caveat prose (currently the bulk
  of the manuscript body).
- S1 (tool versions), S2 (unique-call characterisation), S3 (identity
  sweep), S4 (strict one-to-one panels).
- Figures 3 (sweep + audit), 4 (plant satellites), S1 (chromosome-subset),
  S2 (post-merge trade-off).
- The full Methods 2.1–2.2 text; the main text keeps two-sentence
  summaries with pointers.

## What gets cut nowhere (non-negotiables)

The retractions and limitations stay wherever their claims appear: the
audit result, the 28.08 GB memory, fragmentation, ULTRA-first under every
overlap rule, the S4 corrections trail. An Application Note with a public
evidence repo cannot afford a single suppressed cell.

## Mechanical work list (when the author picks the venue)

1. Freeze the full manuscript as `manuscript_full.md` (supplement source).
2. Write the ~1,300-word main text per the outline above; port the
   two figures once Filip's implementations land (HANDOFF_FILIP.md).
3. Re-cut the abstract to the venue limit; keep "non-leading shared-range
   accuracy" verbatim.
4. Update `results/manifest.tsv` table ids if table numbering changes
   (S-numbers stay stable if the supplement keeps current numbering — do
   that).
5. Citation/format conversion via the venue's template.

## Open inputs

- Venue analysis (running): format limits per venue, precedent papers,
  reviewer-pool risk.
- ULTRA `-p 2000` (job 6145581): the last cell of the range-cost story.
- Title decision (5 candidates pending) — should name the capability, not
  a superiority claim.
