# Kimi re-review after the day's edits — round 3

Run 2026-09-03 at commit bb64ad5, prose inlined, tool use forbidden. The lens was
deliberately not "find the old defects" but "find what the repairs did": a repair that
overshot, a repair that left its other half behind, and seams where a rewritten sentence
sits inside an unrewritten paragraph.

---

Findings, ranked most damaging first. All recomputations done from the numbers in the passage itself.

1. **Severity: HIGH — Location: manuscript.md:311 (and 300)**
   The tantan banded loss is wrong relative to the table's own displayed values. Row: unfiltered 56.51, banded 55.59, "Loss under band −0.91". But 56.51 − 55.59 = 0.92, not 0.91, and the caption repeats the error ("Under the band the re-run loses 0.91 points"). This is exactly the "fixed by re-doing the rounding" pattern, and it is worse than a bare digit slip because the caption guarantees the rows are "directly differenceable" — a property the tantan pair now violates. Every other row recomputes exactly (BWTandem 58.55−40.71=17.84; TRF 58.50−57.63=0.87; ULTRA 57.73−55.84=1.89), so the table otherwise satisfies the guarantee. **Minimal fix:** recompute from the unrounded coverages and either correct both places to −0.92 or state explicitly that losses are computed before rounding; make the row and the caption agree.

2. **Severity: HIGH — Location: manuscript.md:328**
   A false comparative: "TRASH … reaches up to 120.64 GiB on maize. The template-guided NCRF is larger still at 80.96 GiB on Arabidopsis." 80.96 is not larger than 120.64. If the intended sense is footprint per megabase of genome (80.96 GiB / 132 Mb vs 120.64 GiB / 2.18 Gb), the sentence does not say so, and as written a reader is told the smaller number is the larger one. **Minimal fix:** "The template-guided NCRF reaches 80.96 GiB on the 132 Mb Arabidopsis genome — a smaller peak than TRASH's, though a far larger footprint per megabase."

3. **Severity: MEDIUM — Location: manuscript.md:326, 328, 330**
   The enumeration contradicts itself. Line 326 opens "Three limitations bound the results above. The first is period assignment…" and then, in the same paragraph, introduces "A fourth limitation is not about accuracy but about the implementation" — before "The second is memory" (328) and "The third concerns unmatched calls" (330) have appeared. Either there are four limitations and "Three" is stale, or the intermittent-regression-test passage is a separate implementation note that was inserted into the wrong paragraph and never renumbered. This is the clearest seam in the passage: a whole inserted block sitting ahead of the items it was numbered after. **Minimal fix:** move the intermittent-test text after the third limitation and open with "Four limitations…", or detach it as an unnumbered closing paragraph ("A separate issue concerns the implementation…").

4. **Severity: MEDIUM — Location: manuscript.md:315 vs 324**
   The job-records caveat was repaired in one place only. Line 324 now says "These inherited cost figures lack surviving job records", yet line 315 still builds the Section 4.1 design argument on the same inherited TRF timings — "TRF's turns out to be flat on maize (5.2 to 5.5 h; Tables 3A–3C) but not on human at wider ranges" — with no provenance flag, and uses them to conclude "the bound can sit at the top of the plausible range by default". Additionally, line 324 itself says "The inversion shows that genome size alone does not predict the observed runtime": "shows" is too strong a verb for figures the same sentence declares unverifiable. **Minimal fix:** attach the no-job-records caveat at first use in 4.1, and change "shows" to "is consistent with".

5. **Severity: MEDIUM — Location: manuscript.md:330 (revealed only at 332)**
   The audit paragraph lost its antecedents. "The regenerated full-range output has 893,480 calls absent from the other four tabulated tools" — the genome is never named in 4.4; only line 332 ("The unmatched human calls…") lets the reader infer it. And "the other four tabulated tools" are never identified; the benchmark names six competitors, so which four define the uniqueness denominator is not checkable from this passage. An audit's population definition is exactly the thing that must be unambiguous. **Minimal fix:** "The regenerated full-range human output has 893,480 calls absent from TRF, ULTRA, tantan and TRASH (the tools tabulated for human)…"

6. **Severity: MEDIUM — Location: manuscript.md:326**
   Citation and comparability gaps in the banded-loss summary. The sentence cites "(Tables 3B-b, 3C-b)" for losses on three classes — knob180, TR-1, CentC — but per the caption logic at line 300, knob180's banded table would be 3A-b, which is not cited. The TRASH figures quoted inline ("0.34 or less for TRASH … though 18.57 on TR-1") are not in those tables at all; per the caption they come from TRASH's "own period column through its own script", on intervals that are "fixed analysis windows rather than call boundaries". The caption carries that caveat; line 326 drops it while folding TRASH into the same comparison sentence as the call-boundary tools — competitor described on different terms, one paragraph later. **Minimal fix:** cite 3A-b and the caption source for the TRASH numbers, and carry "fixed-window intervals" into the sentence.

7. **Severity: LOW — Location: manuscript.md:330**
   Mixed denominators in the audit statistics. The headline "supported rate of 1.0% overall" is 4/400, while the parenthetical "Wilson 95% CI 0.4 to 2.9% on the definitive verdicts" is on 350 (4/350 = 1.14%; the CI itself recomputes correctly for 4/350). Both numbers are individually defensible and the basis is labelled, but pairing a point estimate on one denominator with an interval on another invites misreading, and the definitive-basis rate (1.1%) is never stated. **Minimal fix:** report the rate and CI on the same basis, e.g. "1.1% of definitive verdicts (Wilson 95% CI 0.4–2.9%)".

8. **Severity: LOW — Location: manuscript.md:328**
   "It is the largest de novo footprint on human" — TRASH ran on human (line 315: 107.6 h) and is "counted here as a de novo tool", but its human footprint is not given in this passage, so the superlative is not checkable from this passage. **Minimal fix:** state TRASH's human peak memory or scope the claim ("largest among the call-boundary de novo tools").

9. **Severity: LOW — Location: manuscript.md:300**
   "As in Table 3B-b, the 2026 tools have no banded rows." No reason is given here, and it is not obvious: AniAnn's has 17 calls and 81.42% unfiltered coverage (line 298), so it is bandable in principle unless it reports no per-call period. Whether the 3B-b caption supplies the reason is not checkable from this passage. **Minimal fix:** add the one-clause reason here too (e.g. "neither reports per-call periods").

10. **Severity: LOW — Location: manuscript.md:315**
    "(TRASH and SRF being exceptions, the former at 107.6 h on human, Section 1)" — SRF appears without introduction (it is not among the benchmarked tools), and a runtime datum is cited to Section 1, the Introduction; whether that number actually lives there is not checkable from this passage. **Minimal fix:** point the citation at the Results section/table carrying the 107.6 h figure, or drop the parenthetical.

Lens accounting:

- (1) Overshot repairs: findings 1 and 2; finding 7 is a milder case of a fix that changed the basis without changing the headline.
- (2) Caveats whose other half was left behind: findings 3, 4, 5, 6, 8.
- (3) Seams: finding 3 is the structural one. Beyond it, I see no further seam that rises to a scientific-integrity defect — the four separate "superseded" disclaimers (315, 319, 328, 336) are each correctly scoped to a different artefact, and the blunt insert at line 324 ("These inherited cost figures lack surviving job records.") is covered under finding 4 rather than being independently wrong.

Recomputations that passed (no finding): TRF/ULTRA/BWTandem banded losses; "finishes 14.88 points above BWTandem" (55.59−40.71); TRASH implied unfiltered CentC 58.68 vs BWTandem 58.55 = the "0.13 points more" at line 326; 4+346+50=400; the ~10⁴ extrapolation from 893,480 × ~1%; the Wilson interval for 4/350.
