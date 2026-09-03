# Codex review of the 12-item work plan (repository-grounded)

Run 2026-09-03 against a detached worktree at `4aeba8e`, with the plan inlined and the
instruction to check it against the repository rather than judge the prose. Its verdict:
"The proposed order should not be executed."

---

The proposed order should not be executed. At `4aeba8e`, several premises are already false, Group 1 is materially incomplete, and its staging/rehashing sequence can create exactly the stale-hash commit it is meant to prevent.

1. **The rehash workflow is backwards and can ship a corrupt evidence commit.**

   The plan says `git add → rehash in background → tests → commit`. The hasher rewrites both checksum files after staging and truncates them before rebuilding them; tests then inspect the working tree, not necessarily the staged versions. A passing guard can therefore accompany a commit containing the old hashes. The script itself says to run it last because this already shipped once: [hash_deposited_beds.sh:18](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/benchmark/hash_deposited_beds.sh:18), while the ledger preserves the unsafe order at [todo.md:58](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/todo.md:58).

   Minimal correction: finish every `results/` mutation, run the hasher in the foreground—or capture and `wait` for its PID and check its status—then stage everything, including both checksum files. Inspect `git diff --cached`, require no unstaged results changes, run the guards, and only then commit. Do not perform two rehash cycles unless an intermediate results commit is genuinely needed.

2. **The twelve-premise inventory is not reliable at `4aeba8e`.**

   The repository supports only some of the descriptions:

   | Plan component | Repository verdict |
   |---|---|
   | Group 1.1 TRASH costs | Numbers verified, but the proposed row semantics are still wrong; see finding 5. |
   | Group 1.2 manifest row 67 | Row exists, but the premise that the manuscript is fully corrected is false: it still says “same FASTA.” |
   | Group 1.3 AniAnn threads | Supported: four manifest rows say 1; sbatch/logs say `-j 2`. |
   | Group 1.4 tantan | Thread value is wrong, but “widened commands unstated” is false in the manifest: rows 68–71 already include `-w 500/500/200/2000`. They are missing from Supplementary Table S1. |
   | Group 1.5 PENDING rows | Supported, but they were failed/OOM submissions, not neutral pending jobs. |
   | Codex round-1 “14” | Asserted by the transcript, but not reconstructable from its verification list. |
   | Codex R2/R3 “17” | Cannot be audited: the two cited transcripts are absent from Git. |
   | Kimi MEDIUM/LOW | No deduplicated inventory or count exists. |
   | Item 6 | False premise: the four values are not the same statistic. |
   | Item 7 | Necessary, but depends on C-1, final figures, archive URL, and actual App Note conversion—not merely a “closed defect list.” |
   | Item 8 | Mostly already implemented. The missing work is release/archive work, not creation of license/package/CI infrastructure. |
   | Item 9 | Six figure tasks exist, but all six implementations are stubs and no final figure is present; this is a submission blocker, not merely something to “track.” |

3. **The “31 unverified findings” count is unauditable and already stale.**

   The deposited round-1 transcript contains 20 numbered findings. Its verification section explicitly maps only five Codex findings—1, 2, 7, 9 and part of 11—plus one Kimi stride finding, yet then says 14 Codex findings remain: [review transcript:337](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/docs/2026-09-03-codex-review-findings.md:337). One additional Codex disposition may have been counted implicitly, but it is not identified.

   More seriously, [todo.md:71](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/todo.md:71) points to `msreview/codex_out_R2.txt` and `R3.txt`; neither file nor the `msreview/` directory is tracked. Thus “17” cannot be checked for duplicates, overlap with round 1/Kimi, or current disposition.

   At least round-1 finding 12 is now substantively resolved by `e4ae632`: the human TRASH nonzero exit is disclosed and the log deposited. Findings 7, 10, 13, 15 and 18/19 have received partial corrections, but live repository surfaces remain inconsistent. R2-8, R2-9, R3-10 and R3-12 were adjudicated today, although without the R2/R3 transcripts it is impossible to establish whether any was included in the alleged remaining 17.

   Minimal correction: deposit the missing transcripts and build one deduplicated disposition table with stable ID, source round, exact claim and base commit, affected path/line, evidence command/file/hash, verdict, rationale, resolving commit or blocker, reviewer and date. Use verdicts such as `confirmed`, `false-positive`, `duplicate`, `already-fixed`, `out-of-scope`, and `unverifiable`. “Not found by grep” is not evidence for `false-positive`; it remains `unverifiable` unless positive counterevidence exists. Preserve transcripts verbatim and append dispositions separately.

4. **C-1 is a hidden dependency, and the current jobs are not three clean performance replicates.**

   The running script does correctly pin `0363d8b` and runs each 100/2,000 pair sequentially. But it requests only four CPUs, lacks `--exclusive`, and uses `set -uo pipefail`, not `set -euo pipefail`: [run_rangerep0363.sbatch:3](/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363/run_rangerep0363.sbatch:3), [line 7](/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363/run_rangerep0363.sbatch:7). A failed arm can therefore be followed by the next arm and a final successful shell exit.

   The logs show replicates 2 and 3 started one second apart on the same node, `cpu-28`: [job 6147699](/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363/rangerep0363_6147699.log:2), [job 6147700](/data/gpfs/assoc/pgl/devel/exp1_human/wp0/rangerep0363/rangerep0363_6147700.log:2). They are co-scheduled technical repeats, not independent timing replicates; memory bandwidth/cache contention can affect both arms.

   C-1 changes the range-replicate manifest rows, `manuscript.md` §3.1, Figure 2A CSV, handoff, title and caption. Group 1 is therefore not independent in the sense relevant to a consistency commit.

   Minimal correction: treat the current runs as exploratory unless accounting establishes acceptable isolation. Prefer rerunning one pair per exclusive node or serializing the replicate jobs. Add `set -euo pipefail`, per-arm exit checks, required nonempty BED/time outputs, expected call-count/hash checks, and accounting validation. Correct the stale job IDs in [todo.md:79](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/todo.md:79).

5. **The proposed TRASH fix still assigns the wrong meaning to the row.**

   The three component BEDs contain 232, 238 and 234 rows. Their full deduplicated union is exactly the published 591-row file. But the union of the two template runs alone is only 397 rows, because they overlap by 73 rows. Thus the existing 591-row result is not “TRASH (template)”; it is “CEN159 template + CEN178 template + de novo combined.”

   The logs verify:

   - CEN159: 6:18:29, 1,389,316 KiB.
   - CEN178: 25:22:40, 2,761,680 KiB.
   - De novo: 5:47:30, 1,357,252 KiB.

   See [the deposited log summary](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/competitor_logs/README.md:25) and the current table at [manuscript.md:189](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:189).

   Minimal correction: choose one of two defensible rows:

   - Keep 591 rows, rename it explicitly as the combined three-run union, and report 37:28:39/2.63 GiB.
   - Produce the true 397-row template-only union, rescore every accuracy column, and report 31:41:09/2.63 GiB.

   Merely changing the 591-row “template” cost to the three-run total leaves a false row label and double-counts the de novo execution conceptually.

6. **Group 1’s blast radius is much larger than the named rows.**

   These are the other surfaces that must be synchronized:

   | Item | Other live surfaces |
   |---|---|
   | TRASH | [manifest rows 52–53](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:52), [competitor-log README](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/competitor_logs/README.md:25), `todo.md`, `quarantine.md`, and the historical Codex transcript. The transcript should remain immutable but receive a resolving disposition. |
   | ULTRA row 67 | [manuscript.md:78](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:78) still says “same FASTA”; [range-cost README:29](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/range_cost_attempts/README.md:29) still says version/input matched and over-localizes the stall; the submitted sbatch and its log preserve the false “same file” provenance; [comparator_baselines.md:8](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparator_baselines.md:8), root README, and Figure 2 handoff preserve absolute/infeasible framing. Submitted/log artifacts should remain verbatim and get an erratum note. |
   | AniAnn | Manifest rows 123/126/129/132 need `2`. S1 and [comparators2026 README:12](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:12) are already correct. But those same files still carry the superseded `0.07/0.50 GB` values instead of `0.06/0.48 GiB`; the plan omits them. |
   | tantan | Manifest rows 68–71 already state the four `-w` values, so that premise is false. The missing command disclosure is [Supplementary Table S1](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:483), while the Introduction still falsely says tantan was left at default: [manuscript.md:38](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:38). External sbatch/logs show both allocation and active threads were one, not “two allocated/one active.” |
   | PENDING rows | Completed replacements are manifest rows 89–93, while [sacct_provenance.txt:94](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/sacct_provenance.txt:94) shows the originals failed/OOM at startup. [results/figures/README.md:3](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/figures/README.md:3) and `make_curve_figure_superseded.py` still say Figure 1 is pending. Use `superseded-failed-at-startup`, not a generic superseded label. |

   Splitting `threads` into `allocated_cpus` and `tool_threads` is a schema migration, not a four-row edit. Update the header, all rows with explicit null semantics, [results/README.md:45](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:45), every consumer, and the external `wp0/make_manifest.py`, which otherwise regenerates the old schema.

7. **The partial-output finding can now be adjudicated; the current wording is wrong.**

   Direct inspection of the external ULTRA TSV found a header, 138,425 complete nine-field records, then one six-field fragment without a final newline. The last complete endpoint is 124,785,432; 124,786,615 comes from the fragment. Thus “138,425 calls … up to 124,786,615” mixes two record populations. The TRF file is likewise unterminated and ends in a truncated 14-field record.

   Minimal correction: state “138,425 complete records plus one truncated fragment; last complete endpoint 124,785,432; fragment endpoint field 124,786,615.” Describe `138426` as a newline count, not a complete logical-line count. Give analogous qualified wording for the TRF partial output. This closes round-1 finding 3 and must be reflected in the manifest and range-cost README, not left in Group 2 limbo.

8. **Item 6 confuses four different statistics.**

   The checkable answer is:

   - **771 bp**: rounded 770.65 bp, BWTandem TR-1 unfiltered/raw-call mean offset in Table 3B.
   - **7,973 bp**: rounded 7,973.15 bp, BWTandem TR-1 after the 100–500 bp called-period band.
   - **4,265 bp**: rounded 4,265.30 bp, but for **knob180 banded**, not TR-1.
   - **4.3 kb**: rounded 4,282 bp for TR-1 at coordinate-only “gap 0” post-merge. That scorer merges overlapping calls even when the allowed positive gap is zero, so it is not the raw-call Table 3B offset.

   The first three are explicit in [table3bc_replacement.md:5](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/regen/table3bc_replacement.md:5) and [lines 16–29](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/regen/table3bc_replacement.md:16). The 4.3 kb sweep is [manuscript.md:238](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:238), and the manifest explicitly warns that the post-merge scorer is not the Table 3B/3C convention: [manifest.tsv:84](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/manifest.tsv:84).

   Minimal correction: replace “the same quantity appears as…” with four named statistic definitions. No numeric reconciliation is needed.

9. **Item 8 is mostly already done; the actual Application Note blockers are missing from the plan.**

   | Check | At `4aeba8e` |
   |---|---|
   | LICENSE | Exists, MIT: [LICENSE](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/LICENSE:1). |
   | CITATION.cff | Exists, CFF 1.2, version 0.9.0: [CITATION.cff](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/CITATION.cff:1). It lists only Won, while the paper has Filip and Won; reconcile whether intentional. |
   | Packaging metadata | `pyproject.toml` and `setup.py` both exist, with dependencies, package data and console entry point: [pyproject.toml:5](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/pyproject.toml:5). |
   | Version string | `0.9.0` exists in pyproject and CFF. No Git tag exists; [CHANGELOG.md:10](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/CHANGELOG.md:10) calls it “pending first tagged release.” `bwtandem/__init__.py` is empty and no CLI `--version` was found. |
   | CI workflow | Exists and tests native/fallback code, hashes, CLI and wheel install: [ci.yml](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/.github/workflows/ci.yml:7). It covers only Ubuntu/Python 3.11, while documentation claims Linux/macOS and Python ≥3.11. |
   | Release/archive | Missing: no tag, no immutable software archive DOI, no archived submitted test-data snapshot, and only a mutable GitHub URL in the abstract. |
   | Reliability | The native memory-layout-dependent regression remains unresolved and is described at [results/README.md:98](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:98). |

   Bioinformatics requires the submitted software version and test data to be archived at a stable repository such as Zenodo/Figshare/Software Heritage and the URL placed in the abstract; it also expects installation documentation, test data, and reproducibility of presented results and figures. Those—not creation of files that already exist—are the real submission tasks. See the official [Bioinformatics author guidelines](https://academic.oup.com/BIOINFORMATICS/pages/author-guidelines).

10. **The figure and Application Note conversion work is much less complete than the plan implies.**

   All six `plot_*.py` files end in `TODO`; there are no tracked final PNG/PDF outputs. The handoff already contains stale C-1 ratios, `GB`, the old “~5%” scope estimate, and “neither competitor could be run” wording: [HANDOFF_FILIP.md:63](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/figures/paper_figs/HANDOFF_FILIP.md:63).

   The conversion document itself begins by calling the Application Note direction superseded: [conversion plan:3](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/docs/2026-09-01-application-note-conversion-plan.md:3). The current manuscript is about 24,541 words and still contains all tables inline. A conversion plan is not a converted submission.

   Minimal correction: make item 9 an owned, dated blocking deliverable with acceptance criteria: six implemented scripts, deterministic regenerated CSVs, final PNG/PDF files, captions synchronized with final data, visual review, and page-fit verification. Then actually execute the App Note restructuring and remove the obsolete “superseded” decision banner.

11. **Today’s changes left multiple referee-visible contradictions that the plan does not mention.**

   At minimum:

   - [results/README.md:25](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:25) says the raw 68,840 CEN180 hits were not retained, although `e4ae632` deposited them.
   - [results/README.md:52](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:52) says competitor cost evidence is unavailable, although 40 logs were deposited.
   - [comparators2026/README.md:14](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/comparators2026/README.md:14) still says “~5% larger”; the corrected value is 3.92%/~4%. Its memory table still says GB and retains 0.07/0.50.
   - [manuscript.md:76](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:76) still claims every competitor run used one Singularity container, contradicting S1 and the 2026/ULTRA attempts.
   - `manuscript.md:99` still says ~5%; `manuscript.md:105` says ~4% was “processed,” although only the last emitted output coordinate is known.
   - [one_to_one/README.md:72](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/one_to_one/README.md:72) says TRF period metrics are disabled, while the deposited JSON now reports 63.53%. Manifest rows 106–120 still identify two older scorer hashes; [todo.md:54](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/todo.md:54) explicitly defers those two hashes, but Group 1 omits them. Do not blindly replace historical execution hashes: add correct provenance for the regenerated JSON or clearly separate historical and replacement results.
   - [MANUAL.md:82](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/MANUAL.md:82) documents a merge gap of `max(10, motif length)`, while the code uses 10× the period for 20–99 bp and 100× for ≥100 bp: [finder.py:373](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/finder.py:373).
   - [CONTRIBUTING.md:40](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/CONTRIBUTING.md:40) says run logs are not in Git, now false.
   - The manuscript contains no identifiable Funding or Conflict-of-Interest section.
   - The manifest still depends on 75 absolute cluster paths. Hashes prove identity to someone on this filesystem; they do not make those artifacts publicly available or permit an external referee to regenerate the figures.

   Minimal correction: add a repository-wide consistency sweep to the plan, driven by facts rather than named files, and make public reproducibility/archive status an explicit release gate.

## REVISED ORDER

1. **Adjudicate first, while C-1 runs.** Deposit R2/R3, deduplicate all Codex/Kimi findings, and produce an auditable disposition matrix. This prevents another fix batch from discovering additional `results/` edits after hashing.

2. **Settle semantic contracts.** Decide combined-vs-template TRASH; define manifest thread/allocation columns; define historical-vs-regenerated scorer provenance; record the four boundary statistics correctly.

3. **Validate C-1 before using it.** Check both arm exits and artifacts. Because two replicates share `cpu-28`, either justify them as co-scheduled technical runs or rerun with isolation. Freeze the final ratios and call counts.

4. **Apply one coordinated consistency change.** Update manuscript, manifest, READMEs, manuals, ledgers, figure data/handoff, pending/failure labels, memory units, stale scorer evidence, and partial-output wording together.

5. **Implement and render all six figures**, regenerate their data from the final evidence, and perform caption/data/visual checks.

6. **Execute the Application Note conversion**, including page fit, supplement construction, funding/conflict statements, and final figure placement.

7. **Finish release readiness.** Resolve or bound the flaky native test, reconcile CFF authors/date, add runtime version reporting, create a real `0.9.0` tag, reserve an archive DOI, and prepare the exact software/test-data snapshot.

8. **Rewrite the structured abstract last**, after C-1, figures, page fit and the stable archive URL are final.

9. **Hash and commit last:** foreground rehash → stage all files → inspect staged diff and require no unstaged results changes → guards/CI → commit → tag/archive.

So yes: landing Group 1 before C-1 is wrong in its current form. The problem is not the extra 600 seconds; it is that Group 1, C-1, Group 2 dispositions, and the figures mutate the same evidence graph. Only the already-present portions of item 8 are unnecessary. Correctness fixes remain necessary even when their tables move to the supplement.

No files were modified; the detached worktree remained clean, and the current `results/manifest.sha256` verifies at `4aeba8e`.
