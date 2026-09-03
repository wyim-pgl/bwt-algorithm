# Codex adversarial scientific-integrity review — findings

Requested by team-lead. Repo tip at request time: `67a3f42` (branch
`perf/exp1-human-sensitivity`). This file is written and updated by the
`codex-review` teammate and completed in the main session. Deposited (tracked) on
2026-09-03, matching how the Kimi transcripts were deposited in `67f11b1`.

Status: **COMPLETE** (attempt 5, `gpt-5.6-sol`). Attempts 1-3 failed for
environment reasons and produced nothing; attempt 4 was aborted mid-run when the
requester chose a different model. The verbatim output below is attempt 5's and it
is final.

**Six of its twenty findings have since been reproduced against the real files by
the main session and are recorded in `quarantine.md` §6.1-6.6; the other fourteen are
NOT verified** (`todo.md` §B — 미검증). Do not act on an unverified finding.

---

## Run history

### Attempt 1 — FAILED (no usable output)

Command:
```
codex exec -m gpt-5.3-codex --config model_reasoning_effort="xhigh" \
  --sandbox read-only --skip-git-repo-check \
  -C /data/gpfs/assoc/pgl/devel/bwt-algorithm \
  "<adversarial review prompt, scoped to commits cb63bc9..ccf04dd>" \
  2>/dev/null | tee <scratchpad>/codex_review_output.txt
```

Exact error: the foreground Bash tool call was killed after its 10-minute
timeout (`Command timed out after 10m 0s`, exit code 143) before Codex
(running at `xhigh` reasoning effort, which explores extensively) produced
any output. The harness auto-backgrounded a *subsequent* diagnostic command
under task id `b327akvp6`; when that completed, its own output showed the
intended tee destination (`codex_review_output.txt`) did not exist on disk at
all — the pipeline was killed before `tee` durably wrote anything usable, and
the underlying `codex exec` process (confirmed via `ps`, running under an
internal `timeout 3000` wrapper) continued orphaned in the background,
disconnected from any captured output, and eventually exited with nothing
recoverable. **No Codex findings were obtained from attempt 1.** No review
content was substituted by the codex-review agent in its place.

### Attempt 2 — FAILED (unsupported model)

Relaunched fully detached (`nohup ... &; disown`) instead of relying on the
harness's foreground timeout/auto-backgrounding, writing stdout/stderr
directly to files so a repeat of the above failure mode is avoidable:

```
cd /data/gpfs/assoc/pgl/devel/bwt-algorithm
nohup codex exec -m gpt-5.3-codex --config model_reasoning_effort="xhigh" \
  --sandbox read-only --skip-git-repo-check \
  -C /data/gpfs/assoc/pgl/devel/bwt-algorithm \
  "<updated adversarial review prompt, scoped to cb63bc9..67a3f42,
     told which 6 defects Kimi already found/fixed (e8ac50c, 67a3f42) so it
     focuses on defects Kimi did not raise>" \
  > <scratchpad>/codex_review_output2.txt \
  2> <scratchpad>/codex_review_err2.txt < /dev/null &
disown
```

PID at launch: 12197.

Outcome: **failed within seconds, and the session then died before it was
noticed.** `codex_review_output2.txt` is 0 bytes; the stderr file ends with

```
ERROR: {"type":"error","status":400,"error":{"type":"invalid_request_error",
        "message":"The 'gpt-5.3-codex' model is not supported when using Codex
        with a ChatGPT account."}}
```

`~/.codex/config.toml` carries the migration notice `"gpt-5.3-codex" = "gpt-5.4"`,
so the model name in attempts 1 and 2 was retired. Probed 2026-09-03 10:56:
`gpt-5.4` and `gpt-5.6-sol` both answer normally on this account. **No Codex
findings were obtained from attempt 2 either.**

### Attempt 3 — FAILED (sandbox cannot start here)

Same detached pattern, model corrected to `gpt-5.4`, `--sandbox read-only`.
Codex started and stayed alive, but *every* command it tried to run failed:

```
exec /bin/bash -c pwd  exited 1 in 0ms:
bwrap: Can't bind mount /oldroot/ on /newroot/: Unable to mount source on destination: Invalid argument
```

`--sandbox read-only` forces Codex's bundled bubblewrap, which cannot construct
its mount namespace inside this Claude session (itself a SLURM job). Codex could
therefore not read a single repository file. Killed rather than left to produce a
review of nothing. This is an environment constraint, not a Codex failure — and
it is the reason the run history is worth keeping: a run that *looks* alive here
may be reading nothing at all.

### Attempt 4 — ABORTED (superseded by attempt 5)

`--sandbox danger-full-access` (no bubblewrap), executed inside a **throwaway
detached git worktree** at `67f11b1` under the session scratchpad, so full access
cannot touch the real working tree; the worktree shares the object store, so
`git log`/`git show` still work. Prompt scope updated to `cb63bc9..67f11b1` and
the prompt's repo paths rewritten to the worktree. Confirmed reading real content
(observed it diffing `results/range_cost_attempts/README.md`).

Aborted after ~5 minutes at the requester's instruction to use `gpt-5.6-sol`
instead of `gpt-5.4`. It had produced no final output yet (stdout still 0 bytes;
~454 kB of tool trace preserved as `codex_review_err4_gpt54_aborted.txt`), so
nothing was lost. Note for whoever repeats this: `pkill -f "codex exec -m gpt-5.4"`
also matches the issuing Bash tool's own command line and kills the caller — kill
by PID, or wait on `/proc/<pid>`.

### Attempt 5 — COMPLETE

Identical to attempt 4 (same prompt, same throwaway worktree at `67f11b1`,
`--sandbox danger-full-access`, `model_reasoning_effort=xhigh`) with
**`-m gpt-5.6-sol`**, the account's configured default model, chosen by the
requester. Confirmed reading real content within 75 s (observed it reading
`manuscript.md` Supplementary Tables S3/S4), and it went on to read the external
GNU-time logs under `/data/gpfs/assoc/pgl/filip/` that the read-only sandbox of
attempt 3 could not have reached. Produced 24,230 bytes / 20 findings. Verbatim
output below.

---

## Codex output (verbatim)

> Path note: Codex ran inside a throwaway detached worktree at `67f11b1` under the
> session scratchpad, so every `file:line` link below is prefixed with that
> worktree path. Strip the prefix through `review-wt/` to get the repository path;
> the line numbers are valid at `67f11b1`. The text is otherwise unedited.

```
## Verdict

The corrected record is still not publication-clean. The most serious new issue is that the Arabidopsis “TRASH (template)” accuracy row scores a three-output union but assigns it the cost of only one component run. I also found several residual Kimi issues—especially in `results/manifest.tsv`—plus provenance, unit, thread-count, and incomplete-output defects.

## Findings

1. **HIGH — The Table 2 TRASH “template” row combines three workloads but reports the cost of one.**

   **Location:** [manuscript.md:187](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:187>), [results/manifest.tsv:53](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:53>), [convert_to_bed.py:396](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/scoring/convert_to_bed.py:396>), [comparator_baselines.md:19](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:19>).

   The scored 591-row BED is the exact sorted union of TRASH’s CEN159-template, CEN178-template, and de-novo outputs. Yet Table 2 assigns it 6.31 h and 1.32 GiB, which reproduce only the CEN159 run. The other logs report 25:22:40/2.63 GiB for CEN178 and 5:47:30/1.29 GiB for de novo. Thus the accuracy numerator and cost denominator describe different workloads, and the statement that every competitor cost cell checks is false.

   **Minimal fix:** Report the three arms separately, or call the row a combined union and account for all constituent runs. Sequential total runtime is 37:28:39; template-only CEN159+CEN178 is 31:41:09. Peak memory among the three is 2.63 GiB.

2. **HIGH — The ULTRA manifest row still preserves several claims Kimi supposedly fixed.**

   **Location:** [results/manifest.tsv:67](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:67>).

   It still says “same binary 1.2.1,” “version- and input-matched,” describes termination at “1.55×” without the lower-bound/non-speed qualification, and attaches “chr1 centromeric satellite region” to the five-hour non-growth observation. The binaries are different files, the literal FASTAs differ, and the stall’s genomic location was not observed.

   **Minimal fix:** Replace the entire description with the corrected manuscript formulation: matching reported version string and sequence content, different binaries/environments, author-cancelled incomplete run, elapsed ratio only a lower bound, and only the last emitted position localized to the centromere.

3. **MEDIUM — The ULTRA output ends in a truncated record; count and endpoint currently come from different record sets.**

   **Location:** [range_cost_attempts/README.md:12](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:12>), [range_cost_attempts/README.md:66](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:66>), manifest row 67.

   The hashed external TSV contains a header, 138,425 complete nine-field records, then a six-field fragment without a final newline. The last complete record ends at 124,785,432; 124,786,615 comes from the truncated fragment. Therefore “138,425 calls … up to 124,786,615” quietly combines the complete-record count with a coordinate from an invalid tail record. The reported 138,426 “lines” is actually a newline count, not the number of logical lines.

   The TRF partial file is likewise unterminated and ends in a truncated record, so its “379,077 lines” wording also needs qualification.

   **Minimal fix:** State “138,425 complete records plus one truncated fragment; last complete endpoint 124,785,432; the fragment’s parsed End field is 124,786,615.” Distinguish newline counts from logical records for both tools.

4. **MEDIUM — “Same FASTA” is not literally true.**

   **Location:** [manuscript.md:78](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:78>), [range_cost_attempts/README.md:29](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:29>), [run_ultra_human_p2000.sbatch:9](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/run_ultra_human_p2000.sbatch:9>).

   The published run used GenBank-style `CM...` identifiers; the p2000 attempt used RefSeq-style `NC...` identifiers. All 455 sequence identifiers differ, and the files differ in size. I verified that the sequence lines themselves are byte-identical—455 sequences and 3,209,286,105 bases—so biological input content is matched, but the file artifact and headers are not.

   **Minimal fix:** Say “the same GRCh38 sequence content, with different RefSeq/GenBank accession headers,” not “same FASTA/file” or “changes ONLY the period.”

5. **MEDIUM — The 1.55× and 4.7× ratios are not marked as lower bounds everywhere.**

   **Location:** [manuscript.md:78](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:78>), [results/manifest.tsv:66](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:66>), [comparator_baselines.md:8](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:8>).

   The dagger treatment in `range_cost_attempts/README.md` is correct. It does not reach the manuscript, either manifest row, or `comparator_baselines.md`. “Incomplete” alone does not explain that these are elapsed-time lower bounds and not speed ratios between completed equivalent runs.

   **Minimal fix:** Add the same explicit caveat at every occurrence or replace the ratios with “at least 1.55×/4.7× the completed capped-run elapsed time before author cancellation.”

6. **MEDIUM — The GiB correction is incomplete across nearly every reporting surface.**

   **Location:** [manuscript.md:80](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:80>), [comparators2026/README.md:25](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:25>), [comparator_baselines.md:26](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:26>), [results/manifest.tsv:1](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:1>).

   A note now admits that “GB” values are numerically GiB, but manuscript prose and table headers, the 2026 README, comparator-baseline table, root README, and manifest column `peak_rss_gb_sacct` still say GB. That remains a unit error, not merely a disclosure issue.

   **Minimal fix:** Relabel every binary-converted value and column as GiB; use GB only for values divided by \(10^9\).

7. **MEDIUM — Two Arabidopsis memory cells are numerically inconsistent with the stated GiB conversion.**

   **Location:** [comparators2026/README.md:27](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:27>), [manuscript.md:197](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:197>), manifest rows 124–126.

   GNU time reports 66,020/66,340 KiB for the longdust arms, which round to 0.06 GiB, not 0.07. AniAnn’s reports 503,968 KiB, which rounds to 0.48 GiB, not 0.50. The reported numbers match decimal-GB conversion, contradicting the newly stated universal KiB/1024² rule.

   **Minimal fix:** Recompute these cells as 0.06 GiB and 0.48 GiB, or explicitly document a different conversion rule.

8. **MEDIUM — “Comparable cost cells” overstates peak-memory comparability.**

   **Location:** [manuscript.md:99](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:99>), [manuscript.md:115](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:115>), [comparators2026/README.md:14](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:14>).

   BWTandem memory is SLURM cgroup `MaxRSS`; 2026 competitor memory is GNU-time process accounting. The manuscript itself concedes these mechanisms are not interchangeable and may miss aggregate worker memory. Therefore runtime/input scope may be comparable, but the peak-memory cells are not method-matched.

   **Minimal fix:** Restrict “comparable” to input scope and wall clock, and explicitly mark peak-memory measurements as methodologically non-comparable.

9. **MEDIUM — AniAnn’s thread counts are wrong in the manifest.**

   **Location:** [results/manifest.tsv:123](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:123>) and rows 126, 129, 132.

   All AniAnn’s rows say `threads=1`, but the deposited documentation and surviving run scripts/logs show `-j 2` under a two-CPU allocation.

   **Minimal fix:** Set those thread fields to 2 and regenerate the manifest checksum.

10. **MEDIUM — “Every competitor run used one Singularity container” is directly false.**

    **Location:** [manuscript.md:76](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:76>), [manuscript.md:481](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:481>).

    The new longdust and AniAnn’s runs were Conda/local executions, and ULTRA p2000 used a local micromamba binary. Supplementary Table S1 itself records those exceptions.

    **Minimal fix:** Limit the statement to the original historical benchmark runs and enumerate the later exceptions.

11. **MEDIUM — Cost-provenance prose remains overbroad and asymmetrically favorable to BWTandem.**

    **Location:** [manuscript.md:80](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:80>), [results/README.md:52](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:52>), [comparator_baselines.md:14](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:14>).

    “With two exceptions, inherited” omits the new 2026 runs and reruns. “Every competitor cost cell” is false because of finding 1. `results/README.md` still says competitor costs “cannot be” checked, contradicting the new surviving-log verification. Finally, describing competitor costs as “approximate” but BWTandem costs as “documented” is advocacy: both are documented single-run point estimates, and BWTandem itself shows material wall-time variation.

    **Minimal fix:** Describe both sides neutrally as documented point estimates, list measurement-method differences, and distinguish “not in the SLURM manifest” from “not externally checkable.”

12. **MEDIUM — The human TRASH cost comes from a nonzero-exit run, without disclosure.**

    **Location:** [comparator_baselines.md:29](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:29>), [manuscript.md:124](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:124>).

    The surviving log records exit status 1 after a Circos plotting error. The error appears to occur after data export, so this does not prove the scored BED is invalid—but the repo currently presents the run as an ordinary successful completion.

    **Minimal fix:** Disclose the exit status and failure stage, and deposit/check evidence showing that the scored export completed before plotting failed.

13. **MEDIUM — The evidence trail is cluster-checkable here but not repository-self-contained.**

    **Location:** [range_cost_attempts/README.md:66](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:66>), [range_cost_attempts/README.md:93](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:93>), [comparator_baselines.md:15](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:15>).

    The large p2000 outputs are externally hashed, but the repo does not contain a compact field-count/tail validation—the missing check that exposed finding 3. Historical GNU-time logs are outside the repo and are not included in `external_evidence.sha256`. The TRF cancellation/accounting record is acknowledged as unavailable. Thus a reader without this filesystem must trust author transcriptions for several core claims.

    **Minimal fix:** Deposit the small GNU-time logs and machine-generated summaries of external outputs, including `NF`, final-newline status, complete-record counts, endpoints, exit status, and SHA-256. Mark the unrecoverable TRF accounting value explicitly as not independently verifiable.

14. **MEDIUM — The manifest contains live-looking PENDING rows superseded by completed rows.**

    **Location:** [results/manifest.tsv:8](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:8>) and rows 16–19, versus completed rows 89–93.

    These duplicate table identifiers remain labeled `PENDING`, without a superseded/failed-startup status, even though the manuscript says the initial jobs failed at startup and were resubmitted.

    **Minimal fix:** Remove the placeholders or give them distinct `superseded-failed-startup` identifiers so automated consumers cannot treat them as current results.

15. **MEDIUM — “About 4% processed/done” is stronger than the output proves.**

    **Location:** [manuscript.md:105](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:105>), [comparator_baselines.md:43](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:43>).

    Evidence establishes the last emitted output coordinate at approximately 3.84% of FASTA bytes, not necessarily the algorithm’s internal processing frontier. This distinction matters because the extrapolation treats that output frontier as progress.

    **Minimal fix:** Say “the last emitted record lay approximately 4% through the input coordinate order; internal processing progress is unknown.”

16. **LOW — The `--adj` code-path explanation is substantively right but literally inaccurate.**

    **Location:** [comparators2026/README.md:78](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:78>), [score_table1.py:403](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/scoring/score_table1.py:403>).

    The README says `adj_rules` is used at line 408 “and nowhere earlier.” Line 403 reads it while deciding whether to print an informational note. That branch is inert in this wrapper because `EXTRA=[]`, so the numeric-equivalence conclusion remains valid.

    **Minimal fix:** Say “the first substantive scoring use is the loop at line 408; the earlier reference only gates an inactive note.”

17. **MEDIUM — The mreps provenance row contradicts the manuscript and its own log.**

    **Location:** [comparator_baselines.md:30](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:30>), [manuscript.md:120](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:120>).

    It says the human mreps cost is “not reported,” whereas the manuscript reports 0.9 h and 6.38 GiB for the partial chr4 run. The surviving log reproduces 54:41.08 and 6,691,220 KiB.

    **Minimal fix:** Record “0.9 h, 6.38 GiB; partial chr4-only run” in the verification table.

18. **LOW — BWTandem is inaccurately described as producing “per-copy” records.**

    **Location:** [comparators2026/README.md:111](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:111>), versus [manuscript.md:34](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:34>).

    BWTandem emits per-array/per-call records containing an aggregate copy count, not a structured record for every repeat copy.

    **Minimal fix:** Use “per-array structured records with motif, period, copy count, purity, and nucleotide-coordinate boundaries.”

19. **LOW — The asserted reason for longdust’s non-detection is unobserved.**

    **Location:** [comparators2026/README.md:104](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:104>).

    Zero overlaps under the tested settings establishes non-detection, not that “these satellites are not low-complexity under its model.” Thresholding, preprocessing, or another code path could explain it.

    **Minimal fix:** Report only the observation: “Neither tested longdust arm called these curated arrays.”

20. **MEDIUM — Some passive/absolute competitor-failure framing survives Kimi’s correction.**

    **Location:** [comparator_baselines.md:42](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:42>), [README.md:176](</data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/README.md:176>).

    “Rerun infeasible” and “competitors’ caps could not be lifted” convert author cancellation decisions into absolute properties of the tools. The experiments demonstrate non-completion within the authors’ chosen budgets, not impossibility.

    **Minimal fix:** Say “the authors stopped the attempts incomplete after…” or “did not complete within the tested allocation.”

## Seven-question disposition

1. **Neutrality/overstatement:** No; findings 1, 2, 4, 8, 10–12, and 18–20 remain.
2. **“Well over a month”:** The arithmetic is defensible as a heuristic. Using roughly 41 active hours gives about 3.0 Mb/h; linear projection is about 45 days for the whole file. Both faster and slower possibilities are explicitly acknowledged, including easier euchromatin versus additional centromeres and the final zero-output period. There is no honest statistical confidence interval available from one censored attempt. The only defect is calling the output frontier “processed” progress, finding 15.
3. **Ratio caveats everywhere:** No. Only `range_cost_attempts/README.md` handles them adequately; see finding 5.
4. **Numerical consistency:** No; findings 3, 6, 7, 9, 14, and 17 are concrete inconsistencies. The 1.55× and 4.7× arithmetic itself is correct.
5. **Raw evidence:** Correctly path-linked and locally checkable in this environment, but not sufficiently self-contained for an external repository reviewer; see findings 3 and 13. All deposited files covered by `results/manifest.sha256` passed checksum verification.
6. **Historical competitor provenance:** The human GNU-time arithmetic checks exactly for ULTRA, TRF, tantan, TRASH, and mreps. The overall provenance claim still fails because of the composite TRASH Table 2 cell, undisclosed TRASH nonzero exit, and incorrect mreps description.
7. **Scoring-log swap:** **No substantive issue found.** The deposited logs differ only in timestamp/host provenance and the crashed run’s trailing 14-line traceback; no scored numeric cell differs. The hashes quoted in the README match. Only the literal “nowhere earlier” code-path sentence needs the low-severity correction in finding 16.

## NOVEL FINDINGS (beyond Kimi’s list)

- The Table 2 TRASH “template” row is a three-output union costed as one run.
- The ULTRA and TRF partial outputs end in truncated records; ULTRA’s count and endpoint mix complete and partial records.
- ULTRA inputs have identical sequence content but different FASTA headers and identifiers.
- Ratio lower-bound caveats do not reach every occurrence.
- Two Arabidopsis memory values use decimal-GB-like numbers despite the stated GiB rule.
- GNU-time and SLURM-cgroup peak-memory cells are presented as comparable.
- AniAnn’s manifest thread counts are wrong.
- Universal Singularity claims contradict the new and p2000 runs.
- Competitor-cost provenance language is overbroad and asymmetrically framed.
- The human TRASH cost comes from an undisclosed exit-status-1 run.
- External evidence is not repository-self-contained.
- Superseded PENDING manifest rows remain live-looking.
- “4% processed” confuses last-output position with internal progress.
- The `--adj` explanation is literally overbroad, though its numerical conclusion is correct.
- The mreps provenance row says “not reported” despite a reported, log-verifiable partial cost.
- BWTandem’s output is mislabeled “per-copy.”
- longdust’s zero overlap is turned into an unsupported mechanistic explanation.

No repository files were modified. The worktree still contains the pre-existing untracked `docs/2026-09-03-codex-review-findings.md`.
```

---

## What the main session verified

Reproduced against the real files, and therefore actionable — full entries in
`quarantine.md` §6.1-6.6:

1. **§6.1 (finding 1, HIGH).** `Col-CEN_v1.2_trash.bed` (591 lines) is the exact
   deduplicated union of three TRASH runs: CEN159 232 + CEN178 238 + denovo 234 =
   704 lines, `sort -u` -> 591, set difference zero in both directions. Both Table 2
   TRASH rows nevertheless print 6.31 h / 1.32 GiB, which is CEN159 alone
   (6:18:29; 1,389,316 KiB / 1024^2 = 1.3249). Sequential total of the three is
   37:28:39; peak RSS across them is 2,761,680 KiB = 2.63 GiB; the de novo run's own
   cost is 5:47:30 / 1.29 GiB. The caption's "inherited values that could not be
   re-verified" exemption lapsed when `a38201b` located the logs.
2. **§6.2 (finding 2).** `results/manifest.tsv:67` still carries all four claims the
   manuscript corrected in `24cd12a`/`e8ac50c`/`67a3f42`.
3. **§6.3 (finding 9).** AniAnn's manifest rows say `threads=1`; the sbatch and every
   run log say `-j 2` under `--cpus-per-task=2`.
4. **§6.4 (finding 7).** Two of four 2026-tool memory cells do not follow the stated
   KiB/1024^2 rule (AniAnn's Col-CEN 503,968 KiB -> 0.48, printed 0.50; longdust
   66,020/66,340 KiB -> 0.06, printed 0.07), while the other two do.
5. **§6.5 (finding 11, part).** "roughly 5% more sequence" against the manuscript's own
   3.25 Gb / 3.15 GB is +3.2%, and the two figures mix bases with file bytes.
6. **§6.6 (Kimi round 1, confirmed in code).** `manuscript.md:66` states stride
   `n / 40,000` saturating at 300 above "roughly 24 Mb"; `bwtandem/tier3.py:31,63`
   gives 12 Mb, and no preset multiplier (`fast 1.6 / balanced 1.0 / sensitive 0.4`,
   `tier3.py:23-25`) yields 24 Mb.

## Not verified

The remaining fourteen findings were not checked and must not be used as the basis
for an edit until they are. They are listed in the verbatim output above and
summarised in `todo.md` §B. Two deserve a note:

- Finding 3 (truncated final record in the ULTRA/TRF partial outputs) would change how
  "138,425 calls ... to 124,786,615" should be worded. It is a claim about external
  multi-GB files and needs its own check.
- Finding 7 of the seven-question disposition states the scoring-log swap has **no
  substantive issue** — that is a negative result confirming `e8ac50c`, and it is the
  one place where Codex agreed with work already done.

## Overlap with Kimi

Codex was told Kimi's six findings and asked to mark overlap. Its own "NOVEL
FINDINGS" section lists seventeen items it considered outside that set. Its
disposition treats (a), (c) and the log-swap question as correctly fixed, and
reports (b), (d), (e) and (f) as **incompletely propagated** — fixed in
`manuscript.md` but surviving in `results/manifest.tsv`, `comparator_baselines.md`,
`results/README.md` and the root `README.md`. That residue is `quarantine.md` §6.2.

(Kimi's six already-applied findings, for reference, so this
section can mark overlap: (a) "every cell identical" claim across the
comparators2026 log swap had no deposited evidence; (b) passive voice made
the authors' own cancellation read as a competitor failure; (c) Results 3.1
"the matched measurement" contradicted Methods 2.2's "not
environment-matched"; (d) "version-matched"/"same binary" overstated what two
binaries both self-reporting "1.2.1" prove; (e) memory figures labelled "GB"
are actually GiB; (f) "spent inside the chromosome 1 centromeric satellite
region" for the stalled 5 hours asserted an unobserved fact, when only the
position of the last emitted calls is known.)
