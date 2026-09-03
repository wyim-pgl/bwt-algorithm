start 22:46:11
##### Q1 evidence-file consistency (22:46:11) #####
I'll read the README and audit it against the stated facts.Audit of `results/range_cost_attempts/README.md` (read only this file; cross-checks are internal arithmetic plus the FACTS supplied).

1. **LOW — ratio framing (b).** Offending: table cell "1.55× its 29.8 h 100 bp run, which covered the whole file" (and TRF cell "4.7× its 33.7 h 500 bp run"). The cell does hint ("which covered the whole file"), and the paragraph at lines 14–17 is explicit and correct ("divide a terminated partial run by a completed one… not like-for-like speed ratios"), so a full read cannot mislead. But a reader skimming only the table sees a bare "1.55×" that parses as a speedup. Minimal fix: add a "†" in both ratio cells pointing to the caveat paragraph.

2. **MEDIUM — partial-TSV path disclosure (c).** Offending: "stays on the cluster at the path recorded in the manifest" (line 56) and "likewise on the cluster path recorded in the manifest" (line 72). The README gives a literal absolute path for the published log (line 29) but only indirection via `results/manifest.tsv` for the TSV itself — the one artefact a reproducer most needs. The SHA-256 reference is good; the path should be equally concrete. Minimal fix: print the literal cluster path here, as is done for the log.

3. **MEDIUM — cross-method memory juxtaposition.** Offending: "batch-step MaxRSS 17,972,740 K (17.14 GB …; the published 100 bp run recorded 1.68 GB under GNU time)" (lines 64–66). A 10× gap between sacct MaxRSS of a local-binary run and GNU-time MaxRSS of a Singularity-wrapped run, with different binaries, is presented with no caveat that the two measurement methods (and builds) are not comparable; a hostile reviewer will read it as a memory regression claim or as noise. Minimal fix: one clause — "not comparable: different binary, environment, and measurement method".

4. **LOW — "GB" for a K/1024² value.** Offending: "17.14 GB by the manifest's K/1024² convention". Dividing KiB by 1024² yields GiB; sacct's "K" is KiB. The README half-acknowledges this ("the manifest's convention") but the unit label is wrong. Minimal fix: write "GiB" or "≈18.8 GB (KiB/1024²)".

5. **LOW — "version-matched" rests on self-report.** Offending: "The two are version- and input-matched" (line 26) vs "both self-reporting version 1.2.1" (line 42). Two different builds self-reporting the same version string are not verified to be the same code; the README's own evidence table undercuts the word "version-matched". Minimal fix: "self-reported-version-matched" or drop the claim to "invocation- and input-matched".

6. **LOW — unverifiable-from-artefact claims.** "beyond the 14-day partition limit" (line 68) and TRF "Threads 1", "379,077 lines", "33.7 h 500 bp run" (line 11) are not checkable from anything cited here. Not in the FACTS either; presumably supported elsewhere in the repo, but this file cites nothing for them. Minimal fix: cite the source (scontrol/partition doc, manifest row for the 500 bp run) or mark as cluster policy.

**Where there is no finding:** every number in the FACTS matches the file — job 6145581, `-t 2 -p 2000`, GCA_000001405.15, 1 d 22 h 15 m (= sacct 1-22:15:07), 138,425 calls / 138,426 lines incl. header, NC_000001.11 to 124,786,615 bp, ~4% of 3.25 Gb (124.79/3250 ≈ 3.8%), 5 h stall (10:05→15:20 on 2026-09-02), MaxRSS 17,972,740 K → 17.14 (KiB/1024²), 29:46:49 = 29.8 h, 1.55× (46.25/29.78 = 1.553), binary sizes 638,912 vs 573,672, TRF 6 d 13 h 57 m = 6.6 d = 4.7× 33.7 h (157.95/33.7 = 4.69). The "CANCELLED by the user… a decision, not a failure" framing, the extrapolation math (~3 Mb/h → >month), the header-line count explanation, and the honest correction of the earlier "same binary" error are all sound and match the FACTS. None.
[[elapsed 325977 ms]]
##### Q2 manuscript Methods (22:51:37) #####
[ASK_KIMI_RESPONSE_TIMEOUT] [ask.runtime] ask could not run because Kimi started and accepted the prompt but never returned a final response. Reduce the prompt scope (or for /kimi:ask and /kimi:rescue, retry with --background to detach). If the response still hangs after a fresh run, check local Kimi version and report upstream. Job 0e106b0f-697d-492f-a248-c8d4948e4802 was persisted as failed.
details: {"availability":"response_timeout","cause_code":"RESPONSE_TIMEOUT","cause_message":"ask.prompt timed out after 900000ms (race detected post-result)."}
[[elapsed 1064006 ms]]
##### Q3 manuscript Results+S1 (23:09:27) #####
Terminated
[[elapsed 1110428 ms]]
##### Q4 log swap provenance (23:28:15) #####
Findings (based solely on this file):

1. **HIGH — "every cell is identical" is asserted, not evidenced.** The claim that 6147179's scored cells match 6146742's (lines 44–47) is unverifiable: the crashed run's log is not deposited, no checksums of the two outputs are recorded, and the diff method isn't stated. The swap could have concealed a real difference (e.g. a nondeterministic scorer, or the `--adj ''` flag changing scoring, not just the corroborator). Minimal fix: deposit the truncated 6146742 log and its partial output alongside 6147179's, with sha256 of each and the verbatim diff (empty, ideally) of the cells feeding the tables.

2. **HIGH — the crash cause is suspiciously convenient and unproven.** The KeyError is attributed to "the wrapper had not passed `--adj ''`" (line 48), but nothing shows the traceback, and the same fix supposedly changes only a post-table section. If `--adj` affects any table-relevant block, the re-run's cells could differ legitimately. Minimal fix: quote the traceback and the scorer commit (`7b2113d`) location of the KeyError, and state explicitly that `--adj` is consumed only by the corroborator block, with a code reference.

3. **MEDIUM — repo-state provenance is asymmetric across the swap.** The crashed run is pinned (`7b2113d`, 0 dirty, converted-BED hashes in log), but the re-run sits at `89ada02` — an unspecified, later commit. The gap between them (which includes `bd32c8c`, the swap commit itself) is undescribed; scorer or wrapper changes in that interval would invalidate "same code" equivalence. Minimal fix: state what changed between `7b2113d` and `89ada02` (ideally: only the `--adj ''` wrapper fix), or show the scorers' hashes were unchanged.

4. **MEDIUM — the scratch deletion muddies the story.** "The re-run also removed `scripts/scoring/work/`" (lines 48–50) conflates the job's action with commit `bd32c8c`'s tree cleanup; a job does not remove a committed directory, a commit does. Also, `2b2beea` committing 104 MB "by mistake" is unverifiable here. Minimal fix: separate the two events ("job 6147179 regenerated outputs; commit bd32c8c deleted the mistakenly committed scratch dir from the tree") and note the scratch dir contained only intermediate artifacts, if so.

5. **LOW — re-run's own integrity fields are thinner.** 6147179 gets repo + dirty + date, but unlike 6146742 there's no mention of converted-BED hashes in its log, and the input FASTA / scorer versions for the re-run aren't re-pinned. Minimal fix: one line confirming the same hashes or that the log records them.

Bottom line: adequate for a trusting reader, not for a skeptical one. The minimum viable proof is the deposited crashed-run log + a byte-level or cell-level diff of the two human outputs + an explicit changelog of `7b2113d → 89ada02`. Without those, the swap rests on the authors' word exactly where the auditability matters most.
[[elapsed 470380 ms]]
end 23:36:09
DONE
