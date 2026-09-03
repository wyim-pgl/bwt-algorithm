# Codex adversarial review — manuscript round 3 (Discussion and Supplement vs code)

Verbatim output of `codex exec -m gpt-5.6-sol` at `xhigh`, run 2026-09-03 in a detached
worktree at `67f11b1`. Fourteen findings, six HIGH. Same path-prefix note as round 2.

Deposited 2026-09-03 for the same reason.

---

## Adversarial review findings

All five requested checks produced at least one finding. I verified the critical code paths against both HEAD `67f11b1` and the benchmark commit `0363d8b`; the discrepancies below were already present in the code that generated the reported runs.

1. **HIGH — The “no ground-truth overfitting” rule was violated, although the Methods partially disclose it.**

   - **Manuscript location:** Discussion 4.4, [manuscript.md:324](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:324): “Three limitations bound the results above.” Ground-truth-guided configuration selection is not one of them.
   - **Repository evidence:** The development plan explicitly calls for diagnosing missed adotto regions on chr21/22, sweeping parameters on those chromosomes, validating code changes there, and running the “best config” genome-wide: [docs/2026-06-20-exp1-human-sensitivity.md:62](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/docs/2026-06-20-exp1-human-sensitivity.md:62). That is supervised tuning against the benchmark truth, despite the same document saying “No GT overfitting” at line 42.
   - The manuscript itself concedes outside this round’s prose that “BWTandem’s configuration was selected empirically” on the same adotto catalog and competitors received no equivalent search: [manuscript.md:94](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:94). The held-out-chromosome check is valuable and shows little measured inflation, but it does not make the categorical “no GT overfitting” rule true or remove the asymmetric optimization.
   - **Minimal fix:** Add configuration selection as a fourth limitation, call it supervised benchmark-guided tuning, make held-out chromosomes the primary human comparison, and either tune competitors symmetrically or explicitly retain the resulting asymmetry.

2. **HIGH — The evidence underlying the central 400-call specificity audit cannot be reproduced as claimed.**

   - **Manuscript location:** Data Availability, [manuscript.md:334](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:334): “the dot-plot renderings regenerate deterministically from the deposited sampler and sheet.”
   - **Repository evidence:** The repository explicitly says the 400 renderings are not deposited: [results/audit11/README.md:21](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/audit11/README.md:21). The alleged sampler only samples calls and writes `reviewer_sheet.tsv` and `answer_key.tsv`; it contains no dot-plot, unit-shift-identity, or random-reference renderer: [sample_specificity_audit.py:127](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/scoring/sample_specificity_audit.py:127). Repository-wide search found no other renderer.
   - Consequently, the exact visual stimuli from which 346 calls were declared unsupported cannot be inspected or regenerated. A reader could invent a new plotter, but that would not recover the visualization settings shown to the reviewer.
   - **Minimal fix:** Deposit the 400 images or the exact rendering program, complete invocation, software versions, and visual parameters; hash them in the manifest.

3. **HIGH — “Scoring scripts” and “ground-truth intervals” are not all available at the repository as stated.**

   - **Manuscript location:** [manuscript.md:334](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:334): “Source code, scoring scripts, ground-truth intervals … are available” at the repository.
   - **Repository evidence:**
     - The repository expressly says the human adotto ground truth is not deposited: [results/README.md:33](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:33).
     - The manifest names `rescore_tables_3bc.py` as the scorer for the regenerated hard-validated maize cells, but that file is absent from the git tree. Its provenance points only to an author-local absolute path: [table3bc_provenance.json:3](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/regen/table3bc_provenance.json:3).
     - The deposited maize wrappers import the also-absent `score_exp3.py` from `/data/gpfs/...`: [score_maize_consistent.py:19](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/scoring/score_maize_consistent.py:19) and [score_trash_maize.py:24](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/scoring/score_trash_maize.py:24).
   - **Minimal fix:** Deposit the exact hashed `score_exp3.py` and `rescore_tables_3bc.py`; either deposit the processed human truth or qualify the Availability sentence and give its immutable download hash plus derivation command.

4. **HIGH — The purported machine-readable run provenance does not record the environment that controls detection.**

   - **Manuscript location:** S2 says “The environment blocks and command line above are the complete configuration” and launch wrappers “add scheduling directives only”: [manuscript.md:479](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:479).
   - **Repository evidence:**
     - `run_with_provenance.sh` promises “exact command -> environment -> output checksum” at [line 4](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/benchmark/run_with_provenance.sh:4), but the JSON it writes contains no environment or input hash: [lines 88–107](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/benchmark/run_with_provenance.sh:88).
     - The deposited P/B/F/H and identity-sweep provenance similarly records only the bare command, not `TIER1_*`, `TIER2_*`, `CATCHALL_*`, `SAT_*`, or `BWT_INDEX_CACHE`: [bwt_human_F_p100_0363.provenance.json:7](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/regen/bwt_human_F_p100_0363.provenance.json:7).
     - The main launcher is not scheduling-only: it exports the scientific configuration at [regen_345.sbatch:55](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/scripts/benchmark/regen_345.sbatch:55). Submission uses `--export=ALL` at line 22, so undeclared inherited knobs were possible unless explicitly cleared.
     - No native-period-100 or sweep launcher is deposited. Their exact environments, including the claim that `BWT_INDEX_CACHE` was unset, are therefore unverifiable from the repository.
     - The available “exact core freeze” pins Python 3.11.15/numpy 2.4.6: [environment.core.lock.yml:1](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/environment.core.lock.yml:1), whereas the main-run provenance records 3.11.14/2.3.1. No exact lock for that main-run environment is supplied.
   - I could verify the three main configurations against non-deposited cluster logs; they match S2. That does not repair repository-level provenance, and the native/sweep environments remain unverified.
   - **Minimal fix:** Serialize a whitelist of every behavior-affecting environment variable and each input SHA-256 into every provenance JSON, submit with `--export=NONE`, deposit all launchers/logs, and provide one exact environment lock per run family.

5. **HIGH — S1.3 reports the wrong Tier-3 jitter tolerance for every benchmark run.**

   - **Manuscript location:** [manuscript.md:438](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:438): “clamps to its 0.04 ceiling in every run, since the tier’s internal maximum period is 100 kb.”
   - **Code evidence:** `finder.py` passes the user’s requested maximum, capped at 100 kb, into Tier 3: [finder.py:130](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/finder.py:130). The formula is `0.02 + 0.02*(max_period/100000)`, then clamped: [tier3.py:38](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/tier3.py:38).
   - Therefore MAXP=2,000 runs used approximately `0.0204`, and MAXP=100 runs approximately `0.02002`; none reached the 0.04 ceiling.
   - **Minimal fix:** Replace the sentence with the realized tolerances for MAXP 100 and 2,000 and explain that the CLI maximum feeds adaptive seeding.

6. **HIGH — S1.4 gives neither the implemented autocorrelation definition nor the implemented valid-base threshold.**

   - **Manuscript location:** [manuscript.md:446](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:446) defines identity as all character equalities divided by `w-p`; line 448 requires “at least 70%” unambiguous bases.
   - **Code evidence:** The scalar satellite calculation excludes comparisons involving non-ACGT bases and divides matches by the number of valid comparisons, not `w-p`: [autocorr.py:60](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/autocorr.py:60). Its validity floor is 0.80: [autocorr.py:41](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/autocorr.py:41), and satellite segmentation enforces that 80% floor: [finder.py:621](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/finder.py:621).
   - Catch-all then uses yet another denominator: matches divided by the full window, with ambiguous positions treated as non-support: [finder.py:708](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/finder.py:708). Thus the manuscript’s single formula does not describe both passes.
   - **Minimal fix:** Give separate scalar-satellite and windowed-catch-all definitions, including ambiguity handling and the actual 0.80 validity threshold.

7. **MEDIUM — The advertised “per-figure provenance manifest” does not exist, and five figure programs are stubs.**

   - **Manuscript location:** Abstract, [manuscript.md:23](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:23): “per-figure provenance manifest.”
   - **Repository evidence:** `results/manifest.tsv` contains no figure mapping. The figure workspace says implementation was handed off and that some inputs are hard-coded from manuscript cells: [paper_figs/README.md:3](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/figures/paper_figs/README.md:3). Each Fig. 1–4/S1/S2 program terminates at `# TODO: implement panels`, e.g. [plot_fig1_accuracy_tradeoff.py:10](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/figures/paper_figs/plot_fig1_accuracy_tradeoff.py:10).
   - The active two-panel Figure 1 does exist, but its README still says it is pending and its CSV is a placeholder: [results/figures/README.md:3](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/figures/README.md:3).
   - **Minimal fix:** Remove the abstract promise until a real per-panel source/output/renderer/hash manifest exists; implement and deposit every referenced figure and replace hard-coded cells with machine-derived inputs.

8. **MEDIUM — A manuscript-plus-repository reader cannot regenerate the complete benchmark.**

   - **Manuscript evidence:** S2 admits that the comparator container has no image checksum and that full mreps and NCRF commands were not preserved: [manuscript.md:481](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:481). Availability admits competitor outputs are not deposited: [manuscript.md:334](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:334).
   - **Repository evidence:** Most comparator and native/sweep manifest paths are private `/data/gpfs/...` paths. Only the three main whole-genome BWTandem BEDs are deposited. The P/B/F/H and identity-sweep BEDs are absent, despite their scores and external paths being cited.
   - This limitation is candidly acknowledged in places, so it is not hidden benchmark gaming. It nevertheless means exact commands, inputs, images and scoring cannot be reconstructed from the claimed deposit.
   - **Minimal fix:** Scope the reproducibility claim explicitly to the three deposited BWTandem rows, or deposit comparator logs/outputs, native/sweep BEDs, complete commands, immutable input hashes and a rebuildable container.

9. **MEDIUM — Discussion 4.3 falsely says FM-index use is confined to Tiers 2 and 3.**

   - **Manuscript location:** [manuscript.md:321](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:321): “BWTandem itself still scans in Tier 1 …; the index-based step is the seeding of Tiers 2 and 3.”
   - **Code/run evidence:** Every paper configuration sets `TIER1_FMSCAN=1`: [manuscript.md:456](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:456). Tier 1 consequently runs an additive FM-index motif-enumeration pass on the residual sequence: [tier1.py:325](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/tier1.py:325). The surviving run logs also print thousands of “Tier 1 FM-index scan” calls per chromosome.
   - **Minimal fix:** State that Tier 1 consists of a sliding scanner followed by an active additive FM-index enumeration, while Tiers 2/3 use FM-index k-mer seeding.

10. **MEDIUM — Discussion 4.2 gives an incomplete genome-by-genome account of the supplementary passes.**

   - **Manuscript location:** [manuscript.md:317](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:317) attributes Arabidopsis to gap-fill and human to catch-all, but says nothing about maize.
   - **Code/run evidence:** Satellite gap-fill is unconditional after all three tiers: [finder.py:285](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/finder.py:285). S2 enables catch-all for both human and maize: [manuscript.md:469](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:469).
   - Direct BED inspection found 404,622 maize tier-5 catch-all calls and 7,819 tier-4 gap-fill calls. Human likewise contains both passes; Arabidopsis contains gap-fill but no catch-all.
   - **Minimal fix:** State the full matrix: satellite gap-fill ran for all three genomes; catch-all ran for human and maize only. Distinguish “was present” from “was proven necessary” unless an ablation supports dependence.

11. **MEDIUM — S1.2 misidentifies the 2% primitive-period test as autocorrelation.**

   - **Manuscript location:** [manuscript.md:415](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:415): approximate periodicity at 2% is “implemented via autocorrelation.”
   - **Code evidence:** `_reduce_to_primitive` calls `smallest_period_str_approx(..., 0.02)`: [tier2.py:42](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/tier2.py:42). That routine compares each character with the prefix template modulo candidate period and divides mismatches by motif length: [motif_utils.py:916](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/motif_utils.py:916). Autocorrelation appears only in a later, restricted fallback for unreduced short consensuses.
   - **Minimal fix:** Describe the 2% operation as cyclic-template mismatch testing; document the separate short-consensus autocorrelation fallback independently.

12. **MEDIUM — The Limitations section omits an unresolved native-path correctness flake.**

   - **Manuscript location:** Discussion 4.4 claims exactly “Three limitations”: [manuscript.md:324](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:324).
   - **Repository evidence:** A ground-truth sensitivity test fails in roughly 60% of full-suite runs depending on process memory layout: [results/README.md:98](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/README.md:98). The root cause remains described as “a layout-conditional read somewhere in the native path,” and three native libraries remain outside the disable-native parity harness: [CLAUDE.md:88](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/CLAUDE.md:88).
   - The repository presents evidence that one 18.8 Mb chromosome remained byte-identical and asserts no figure is affected. That reduces demonstrated impact but does not resolve the native correctness defect.
   - **Minimal fix:** Disclose the unresolved native-path flake, evidence limiting its observed scope, and incomplete accelerator parity coverage.

13. **MEDIUM — S1.3’s “values in force” table silently excludes processed Col-CEN organellar chromosomes.**

   - **Manuscript location:** [manuscript.md:427](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:427) says the table gives values across “the chromosomes of the three genomes”; it reports Arabidopsis k-mer 24–25, stride 300, max occurrences roughly 730–1,090.
   - **File/code evidence:** The recorded Col-CEN input includes ChrC and ChrM, and the deposited BED contains 156 ChrC plus 291 ChrM calls. Because parameters are computed per sequence, the implemented formulas give:
     - ChrC: k=12, stride=20, max occurrences=200.
     - ChrM: k=13, stride=20, max occurrences=200.
   - These are outside every Arabidopsis range printed in S1.3. See the per-sequence formulas at [tier3.py:30](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/tier3.py:30) and the deposited output description at [results/beds/README.md:8](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/results/beds/README.md:8).
   - **Minimal fix:** Say the printed ranges cover the five nuclear chromosomes only, and add the organellar micro-mode values or exclude those sequences explicitly from the reported output count.

14. **LOW — “Up to two million bases” is not the implemented Tier-1 sampling cap.**

   - **Manuscript location:** [manuscript.md:413](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:413).
   - **Code evidence:** `_base_freqs` uses `step = n // 2_000_000` and samples `range(0,n,step)`: [tier1.py:360](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/bwtandem/tier1.py:360). Floor division does not cap the sample count; for reported Col-CEN Chr4 it samples about 2.16 million positions, and inputs just below 4 million bases would sample nearly all 4 million.
   - **Minimal fix:** Say “approximately two million” or use `ceil(n/2_000_000)` in the implementation.

## Checks that did pass

The three promised main BWTandem BED archives are present, their line counts exactly match 161,330/4,014,108/2,406,800, and all three full SHA-256 entries validate. I found no discrepancy in the S1.1 main sliding-pass gate or S1.5 index-storage description. The non-deposited cluster logs also confirm the stated main human/maize/Col-CEN environment settings, but not the native-period-100 or identity-sweep environments.

I did not modify the worktree. GitHub issue pages #14/#26–33 could not be independently fetched in this environment; `gh` was unavailable and web retrieval failed, so issue-related conclusions above rely only on the readable repository documents and `resume.md`.
