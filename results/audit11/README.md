# Issue #11 specificity audit (single reader)

Stratified blinded audit of the 809,886 regenerated full-range BWTandem calls
overlapping neither the adotto catalog nor any of TRF, ULTRA, tantan, TRASH
(Methods 2.2.3; sampler `../../scripts/scoring/sample_specificity_audit.py`,
seed 20260827, 100 calls per period stratum 1-6 / 7-20 / 21-100 / 101-2000 bp).

- `reviewer_sheet.tsv`: the source-blinded sheet given to the reader
  (coordinates, proposed period, sequence with flanks; no stratum totals).
- `verdicts_reviewer2.tsv`: the completed single-reader verdicts, 400/400
  (one `SUUPORTED` typo normalised to SUPPORTED in aggregation; the reader confirmed the SUPPORTED intent on 2026-08-31).
  The protocol's second reader was not completed before submission;
  `verdicts_reviewer1.tsv` therefore does not exist here and the manuscript
  reports the audit as single-reader.
- `answer_key.tsv`: the sealed key (sample_id to source BED line), opened only
  after all 400 verdicts were recorded.
- `aggregate_reviewer2_20260831.txt`: per-stratum aggregation with Wilson 95%
  CIs. Headline: 4 SUPPORTED / 346 UNSUPPORTED / 50 UNSURE, supported rate
  1.0% (CI 0.4-2.9% on definitive verdicts), 0 supported above period 20 bp.

The 400 dot-plot renderings (32 MB) are not deposited; they regenerate
deterministically from the sampler, the deposited sheet and the hg38 FASTA.
