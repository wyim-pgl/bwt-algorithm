# Kimi re-review after the day's edits — round 2

Run 2026-09-03 at commit bb64ad5, prose inlined, tool use forbidden. The lens was
deliberately not "find the old defects" but "find what the repairs did": a repair that
overshot, a repair that left its other half behind, and seams where a rewritten sentence
sits inside an unrewritten paragraph.

---

**Findings, most damaging first**

1. **Severity: HIGH — manuscript.md:113**
"tantan finished fastest" is contradicted by the paper's own Table 1a, where longdust ran in 0.47 h (line 125) against tantan's 0.9 h (line 123) — and by the manuscript's own opener, which credits longdust with "28 minutes" (line 103). mreps' partial run also ties tantan at 0.9 h (line 120). This is a leftover superlative from before the 2026 rows were added: the repair (adding the rows, and line 103's new sentence) left this sentence behind. Minimal fix: qualify it ("tantan was the fastest of the 2024 tools") or delete the superlative.

2. **Severity: HIGH — manuscript.md:140 vs lines 147–148**
Table 1c's caption states "TRASH's stratum row is not reported below" and "The 2026 tools likewise have no rows … both are reported at Table 1a's full-range rule only." The table directly below reports both: AniAnn's (205 calls, line 147) and TRASH (1,427 calls, line 148). The caption even argues AniAnn's numbers "describe whole arrays rather than calls comparable to this panel's" — while the panel includes them. The caption's first sentence ("the stratum BWTandem, TRF, TRASH and the re-run tantan entered") was updated for TRASH; the "not reported" sentence was not. Minimal fix: rewrite the caption to match the table, reusing the Table 1b formula ("reported here rather than withheld", line 129) or remove the rows.

3. **Severity: HIGH — manuscript.md:113**
Two citations point at values Table 1b does not contain: "the regenerated native matched-range value is 32.55%, Table 1b" and "within 1.9 at the native matched range (79.88 against 78.00, Table 1b)". Table 1b's BWTandem cells are 30.00% BP recall and 78.87% region recall (line 134) — the filtered arm — and lines 107 and 129 state explicitly that the native `--max-period 100` rerun is "reported as a sensitivity analysis rather than as the comparison", i.e. kept out of Table 1b. The value 32.55% appears nowhere else in the passage and is not checkable from it. A repair inserted native-arm caveats with "Table 1b" pointers as if the table had been switched to the native arm; it wasn't. Minimal fix: cite the sensitivity analysis (SLURM job 6143150 / Table 1d), not Table 1b.

4. **Severity: MEDIUM — manuscript.md:105**
"the 1.33 mean between-arm ratio" cannot be reproduced from the printed numbers: the mean of the three printed replicate ratios 1.30, 1.30, 1.41 is 1.3367, which rounds to 1.34; including the original pair's 1.40 gives 1.35. This is the overshoot case — a rounding fix that re-did the rounding wrongly. (The contrast with the 1.20 within-arm factor survives either way.) Minimal fix: print 1.34, or state the mean is computed from unrounded runtimes.

5. **Severity: MEDIUM — manuscript.md:140**
"the per-stratum recalls of Tables 1b and 1c sum to at least the full-range figure of Table 1a" is stated universally, but it fails for the newly added AniAnn's rows: 0.01 + 0.03 = 0.04 < 0.10 (45 of its 364 calls fall outside the two strata). The rule holds only for tools whose full searched range is 1–2,000 bp; the sentence was written when the panel had only such tools. Minimal fix: scope the claim ("for tools whose searched range is covered by the two strata") or footnote the AniAnn's exception.

6. **Severity: MEDIUM — manuscript.md:107, 109, 140, 150 (seam)**
The supersession of old-build values is stated three times in three registers: "The old-build value was removed rather than mixed with regenerated whole-genome results" (107), "The superseded operating-point values are not carried forward" (109), "the withdrawn earlier-build CSV, renderer, PNG and PDF are retained only under filenames containing `superseded`" (150) — and line 140's "BWTandem's Table 1b row is *now* that same post-hoc figure" narrates the edit itself. These read as editing-log entries left in Results prose; the fact is stated once per repair instead of once per manuscript. Minimal fix: state provenance once (caption or manifest reference) and delete the other two.

7. **Severity: LOW — manuscript.md:105 (seam, register)**
"That the seeding is done through the index turns out not to be what the flatness rests on." sits immediately after a sentence asserting the opposite design rationale ("Seeding candidate positions once from the index … is the design intended to decouple the two"), and the chatty "turns out not to be" clashes with the paragraph's hedged register. The same paragraph also says "unlogged" twice ("in an unlogged side measurement" / "retained as unlogged engineering observations"). Minimal fix: merge rationale and retraction into one sentence ("…was the intended mechanism, but the ablation below shows it is not the operative one"); drop one "unlogged".

8. **Severity: LOW — manuscript.md:105**
"differ in wall clock by 5.5%, 21 m 18 s against 22 m 29 s": the 71 s difference is 5.6% of the faster time and 5.3% of the slower; 5.5% matches neither natural denominator. Minimal fix: state the denominator ("the index arm is 5.3% slower").

9. **Severity: LOW — manuscript.md:109 vs 150 (terminology seam)**
Operating point P is "catch-all disabled" at line 109 but "P disables the optional periodicity pass" at line 150. Whether these name the same feature is not checkable from this passage; if they do, one term should be used in both places.

**Lens-by-lens accounting (explicit no-findings)**

- (1) Overshoot: finding 4 is the re-done rounding. No other repaired sentence was found that now claims something new or hedges into vacuity — the heavy qualifiers at 105 ("indicative rather than matched", "no material effect rather than a precise invariance") and 140 ("came from a shared node and is an upper bound") are checkable and appropriately scoped.
- (2) Left-behind halves: findings 1, 2, 3, 5. This was the productive lens; the 2026-row and native-arm repairs each updated exactly one of the two places that needed it.
- (3) Seams: findings 6, 7, 9. No tense shifts beyond those noted.
- Competitors on different terms: no new finding. The filtering-vs-rerunning asymmetry (107, 129), the tantan re-run provenance ("a parameter of ours rather than a property of the tool", 140), and the aborted ULTRA/TRF period-matched runs (105) are all disclosed; longdust and AniAnn's are described on the same full-range rule as BWTandem, with their structural absences explained in both text and captions (consistently, except at Table 1c — finding 2).
- Arithmetic recomputed and *confirmed*: 2.75-point gap (81.62−78.87); 1.01 points (79.88−78.87); 0.02-point H-gap; 5.22-point precision gap (53.66−48.44); 1.40 runtime ratio (443/317 min); 1.20 within-arm factor (397/330 min); 83- and 151-fold ULTRA figures (500/6, 34.7/0.23); 19× and 104× memory (28.08/1.45, 28.08/0.27); 4+346+50=400, 4/400=1.0%, 4/350=1.1% with Wilson CI ≈0.4–2.9% (the CI attaches correctly to the excluding-UNSURE rate); strata sums for BWTandem, TRF, TRASH and tantan all satisfy the ≥ rule.
- Not checkable from this passage: all Table 1d values (54.69/81.60/62.16/48.44), the Arabidopsis ablation numbers (161,187 calls, 3.39 vs 1.30 GiB), the cache timings (0.6 s / 9 s / 21 min), the 32.55% native BP recall (finding 3), and the commit/job identifiers.
