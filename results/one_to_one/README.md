# Strict one-to-one, boundary-aware scoring (issue #15, Supplementary Table S4)

Every accuracy figure in Tables 1a–1e uses the one-base-pair region rule, and
the permissive matcher in `tests/test_ground_truth.py` (many-to-one,
integer-multiple/±20% period pooling) remains the named historical/regression
metric. The JSONs here are the strict counterpart, produced by
`scripts/scoring/score_one_to_one.py`:

- **One-to-one assignment** — maximum-cardinality bipartite matching
  (Hopcroft–Karp) between truth regions and predictions; each record
  participates in at most one match. Among equal-cardinality matchings the
  pairing is arbitrary (overlap is not a secondary objective), so the
  boundary/period/copy statistics carry that disclosed indeterminacy.
- **Eligibility** — reciprocal overlap: the intersection must cover ≥50% of
  the truth region *and* ≥50% of the prediction.
- **Boundary error** — |Δstart| and |Δend| per matched pair (median, p90,
  exact-boundary fraction).
- **Period agreement** — exact, integer-multiple, ±20%, and outside, reported
  separately instead of pooled.
- **Copy-number error** — median relative error over pairs where both sides
  carry a copy count (> 0 on the truth side).
- **Stratification** — by truth period band (1–6, 7–20, 21–100, 101–2000 bp).

## Provenance

SLURM job 6145846 (cpu-51, 2026-09-01), repo commit `4a1d368`, 0 dirty
entries; the job log records the full SHA-256 of the scorer and of every
input BED. Inputs are the Table 1a human baselines exactly as listed in
`results/manifest.tsv` (BWTandem deposited BED, ULTRA, tantan, TRF), scored
**at each tool's full published range** — unlike the reciprocal-overlap
sensitivity analysis under `results/regen/`, which restricted BWTandem and
TRF to the matched range (periods ≤100 bp). ULTRA and tantan call nothing
above period 100, so their figures are identical under either restriction;
BWTandem's 0.50-reciprocal region figure is 19.52% here versus 17.61%
matched-range, and TRF's 11.04% versus 9.80%.

- `one_to_one_<tool>_r50.json` — truth = adotto v1.2.1 region GT
  (`adotto_primary.bed`, 1,784,804 regions, coordinates only, chr1–22XY).
  Region/boundary metrics; period and copy metrics are structurally absent
  (the region GT carries no motif).
- `one_to_one_<tool>_annot_r50.json` — truth = the same regions annotated
  with the motif/copies of each region's top-score catalog annotation
  (`scripts/scoring/derive_adotto_annotated_truth.py`, a disclosed
  simplification for compound regions; all 1,784,804 rows carry a motif).
  Adds period agreement, copy-number error, and period-band strata. TRF has
  no row in this arm: its BED's column 4 is the full array sequence, so a
  length-derived period would be meaningless.

## Headline numbers (reciprocal 50%, one-to-one)

| Tool | Sens. (%) | Prec. (%) | |Δstart| med (bp) | Period exact (%) | Copy med rel err (%) | 1–6 (%) | 7–20 (%) | 21–100 (%) | 101–2000 (%) |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| BWTandem | 19.52 | 8.68 | 25 | 54.40 | 45.0 | 17.29 | 15.25 | 42.50 | **41.76** |
| ULTRA | 31.05 | 16.56 | 20 | 49.85 | 78.8 | 25.38 | 34.02 | 64.84 | 7.91 |
| tantan | 6.37 | 3.28 | 25 | 63.56 | 94.3 | 4.48 | 4.77 | 25.76 | 1.28 |
| TRF | 11.04 | 19.49 | 25 | — | — | — | — | — | — |

The absolute values are low for every tool because the adotto regions are
merged and slop-padded, so 50% reciprocal overlap against a single call is
unattainable by construction for much of the catalog (the same caveat as the
reciprocal-overlap sensitivity analysis in Methods). What the strict metric
adds over that analysis: the ULTRA-first ordering persists under one-to-one
assignment; BWTandem leads every tool in the 101–2000 bp stratum (41.76%
versus ULTRA's 7.91%, whose default range cannot reach it); and among the
matched pairs BWTandem has the lowest copy-number error (45% median relative
error versus ULTRA's 78.8% and tantan's 94.3%) and a higher period-exact
rate than ULTRA (54.4% versus 49.9%; tantan's 63.6% is computed over a
matched set one fifth the size).
