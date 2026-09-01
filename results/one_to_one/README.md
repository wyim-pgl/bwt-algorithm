# Strict one-to-one, boundary-aware scoring (issue #15, Supplementary Table S4)

Every accuracy figure in Tables 1a–1e uses the one-base-pair region rule, and
the permissive matcher in `tests/test_ground_truth.py` (many-to-one,
integer-multiple/±20% period pooling) remains the named historical/regression
metric. The JSONs here are the strict counterpart, produced by
`scripts/scoring/score_one_to_one.py`.

## Metric definitions

- **One-to-one assignment** — maximum-cardinality bipartite matching
  (Hopcroft–Karp) between truth records and predictions; each record
  participates in at most one match. Among equal-cardinality matchings the
  pairing is arbitrary (overlap is not a secondary objective), so the
  boundary/period/copy statistics carry that disclosed indeterminacy —
  the sensitivity/precision values themselves do not (maximum cardinality
  is unique).
- **Eligibility** — reciprocal overlap: the intersection must cover ≥50% of
  the truth record *and* ≥50% of the prediction.
- **Sensitivity** = matched / truth records; **precision** = matched /
  predictions (after the chromosome scope filter below). These are
  record-count fractions, distinct from Table 1a's ≥1 bp region-overlap
  fractions. Because the truth records do not overlap one another, the 50%
  reciprocal requirement makes per-truth matching essentially injective
  already: **one-to-one assignment leaves sensitivity identical to the
  many-to-one `bedtools -f 0.50 -r` baseline** (`results/regen/
  recip_0.50.txt`, full-range block) **and changes only precision**, by
  charging every additional call a tool makes on the same truth record.
- **Boundary error** — |Δstart| and |Δend| per matched pair (median, p90,
  exact-boundary fraction).
- **Period-length agreement** — exact, integer-multiple, and within-20%
  (smaller-period-relative, matching the permissive matcher's rule),
  reported separately. This compares period *lengths*, not motif sequences.
- **Copy-number error** — median relative error over pairs where both sides
  carry a copy count. **Only BWTandem's BED carries a copy count in column
  5**; the converted ULTRA/tantan/TRF/TRASH baselines carry the *period*
  there (`convert_to_bed.py`), so the metric is reported for BWTandem alone
  (`--pred-col5 period` disables it rather than mis-scoring it).
- **Strata** — sensitivity by **truth-side period band** (1–6, 7–20,
  21–100, 101–2,000 bp of the truth annotation's primitive period). This is
  a different quantity from Table 1c's prediction-side period strata.

## The two truth sets

- `*_r50.json` — **region truth**: `adotto_primary.bed` (1,784,804 regions,
  coordinates only, chr1–22XY). The catalog's regions are merged and
  slop-padded (~25 bp beyond the underlying annotation at both ends), which
  both makes 50% reciprocal overlap against a single call unattainable for
  much of the catalog and imposes a ~25 bp floor on boundary error — the
  region-arm boundary medians measure that padding, not tool accuracy.
- `*_annot_r50.json` — **annotation truth**: each region's top-score
  catalog annotation (deterministic tie-break: score, span, leftmost) at
  the *annotation's own coordinates*, with the motif primitive-reduced and
  the copy count rescaled by the reduction factor
  (`derive_adotto_annotated_truth.py`; a disclosed simplification for
  compound regions — a tool correctly reporting one of a region's *other*
  motifs scores as a period disagreement). Tighter coordinates make
  matching easier and boundary error meaningful: the best tools sit at
  0–2 bp median |Δstart| here versus the padded 20–25 bp of the region arm.

## Provenance

SLURM job **6146229** (cpu-51, 2026-09-01), repo commit `43543da`, 0 dirty
entries; the job log records the full SHA-256 of the scorer and of every
input BED. Inputs are the five Table 1a human baselines exactly as listed in
`results/manifest.tsv` (mreps is not scored on human; see Methods). Every
arm ran with `--truth-chroms-only`: predictions on sequences absent from the
chr1–22XY truth (ULTRA 130,506; tantan 149,706; TRF 48,568; TRASH 198;
BWTandem 0 — it ran on the primary FASTA) are dropped before precision is
computed, so denominators are comparable across tools run on different
FASTA scopes. ULTRA, tantan and TRASH used `--pred-col5 period`; TRF used
`--pred-col5 period --pred-motif-is-sequence` (its column 4 is the full
array sequence, so prediction-period metrics are disabled while truth-based
strata still fill). Ranges are those of each tool's Table 1a run (ULTRA
default ≤100 bp, tantan default ~100 bp window, TRF ≤500 bp, BWTandem
1–2,000 bp) — *not* each tool's maximal supported range; the tantan
2,000 bp re-run of Table 1c is not included here.

### Supersession note

`superseded-v1/` holds the first deposited measurement (job 6145846, commit
`4a1d368`). It is superseded because (a) it read the ULTRA/tantan column 5
as a copy count when it is the period, so its "copy-number error" for those
tools compared a period against a copy count; (b) it counted scaffold calls
in competitor precision denominators; (c) its annotated arm used the padded
region coordinates and un-reduced motifs. Its claims "BWTandem has the
lowest copy-number error" and "BWTandem leads the 101–2,000 bp stratum" do
not survive the corrected measurement and are retracted.

## Headline numbers (reciprocal 50%, one-to-one)

Region truth (sensitivity equals the many-to-one baseline by construction;
precision is what the 1:1 assignment adds):

| Tool | Matched | Sens. (%) | Prec. (%) |
|---|--:|--:|--:|
| ULTRA | 554,179 | 31.05 | 17.23 |
| BWTandem | 348,464 | 19.52 | **8.68** (lowest) |
| TRF | 197,087 | 11.04 | 20.47 |
| tantan | 113,649 | 6.37 | 3.42 |
| TRASH | 1,241 | 0.07 | 24.17 |

Annotation truth:

| Tool | Matched | Sens. (%) | Prec. (%) | Δstart med (bp) | Period exact (%) | Copy med rel err (%) | 1–6 | 7–20 | 21–100 | 101–2000 |
|---|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|
| ULTRA | 1,122,910 | 62.92 | 34.91 | 1 | 59.61 | — | 57.77 | 72.32 | 80.90 | 8.36 |
| BWTandem | 1,050,831 | 58.88 | 26.18 | 2 | 58.66 | 34.15 | 59.21 | 58.24 | 62.38 | 46.39 |
| tantan | 961,163 | 53.85 | 28.95 | 0 | 73.54 | — | 60.31 | 49.27 | 41.41 | 1.67 |
| TRF | 518,719 | 29.06 | 53.87 | 0 | — | — | 17.87 | 33.88 | 77.60 | 62.40 |
| TRASH | 892 | 0.05 | 17.37 | 494 | 33.41 | — | 0.01 | 0.02 | 0.21 | 0.86 |

What the strict metric shows, stated plainly: ULTRA leads one-to-one
sensitivity in both arms and precision among the three high-recall callers;
BWTandem is second on sensitivity in both arms and **has the lowest
region-arm precision of the five tools** — the same per-region call
granularity (several calls per catalog region) disclosed as fragmentation
in the manuscript's limitations, charged per extra call by the 1:1
assignment. In the 101–2,000 bp truth band **TRF leads at 62.40%** with
BWTandem second at 46.39%; ULTRA (8.36%) and tantan (1.67%) emit no call
above period 100, so the calls matching those truths are their shorter-
period calls. Each tool's matched-pair statistics are computed over its
*own* matched set, which differ in size and composition (e.g. the 21–100 bp
band is 11.4% of ULTRA's matched pairs but 6.8% of tantan's); the
period-exact ordering (tantan 73.54 > ULTRA 59.61 ≈ BWTandem 58.66) should
be read with that selection effect, and tantan's smaller matched set, in
mind. Copy error is measurable only for BWTandem (34.15% median relative
error against the rescaled annotation copy count).
