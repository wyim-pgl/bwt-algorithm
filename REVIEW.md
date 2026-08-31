# Pre-Submission Review — bwtandem

**Date:** 2026-08-27
**Reviewer:** Codex (rescue pass, full-repository audit)
**Codex session ID:** `01a041b5-a458-7fc1-9e86-7435f4eeba42`
**Resume:** `codex resume 01a041b5-a458-7fc1-9e86-7435f4eeba42`

---

## 1. Summary verdict

No—the project is not ready for paper submission yet. The core detector is functional, the native accelerators build, all 34 tests pass, the synthetic comparison is reproducible, and a fresh Chr4 run completes successfully. However, manuscript-facing genome benchmarks are now stale relative to a confirmed satellite-scanner bug fix; ColCEN coverage used a coordinate conversion error; several conclusions exceed what the stored evidence establishes; validation is almost entirely synthetic or annotation-coverage based; and key data, logs, documentation, and the entire Oropetium validation remain untracked. These are submission blockers, but the codebase is reasonably close once the benchmark suite is rerun and the validation/methods package is tightened.

## 2. Confirmed facts about current state

- The maintained `bwtandem` environment uses Python 3.11.15. The committed [environment.yml](environment.yml) requested Python 3.13 and originally omitted test and plotting dependencies.

- The complete test suite passes after the fixes: **34 passed in 5.11 seconds**. Five warnings concern class-scoped pytest fixtures becoming unsupported in pytest 10.

- The deterministic stress test produced:

  - 107/108 true repeats detected
  - 1 false negative
  - 0 false positives
  - 99.1% sensitivity
  - 100% precision
  - 99.5% F1

  This differs from the former README table.

- The four-tool synthetic comparison script reproduces the 44-repeat table exactly:

  | Tool | Sensitivity | Precision | F1 |
  |---|---:|---:|---:|
  | bwtandem | 100.0% | 100.0% | 100.0% |
  | TRF | 97.7% | 100.0% | 98.9% |
  | mreps | 68.2% | 7.6% | 13.6% |
  | ULTRA | 72.7% | 71.1% | 71.9% |

- The standalone ctypes libraries in [src/c_extensions](src/c_extensions) build and load successfully. The Cython extension also rebuilds successfully when NumPy's include directory is supplied.

- Tier 1's implementation in [tier1.py](src/tier1.py) is a direct period-k sliding scan, followed by mismatch extension and alignment refinement. It does not use the size-dependent FM-index enumeration described previously in README and CLAUDE.md.

- Tier 2 in [tier2.py](src/tier2.py) uses LCP-derived candidates plus BWT k-mer seeding. The coordinator scans periods 10–50 in the general phase and periods ≥20 through the long-unit phase.

- Tier 3 in [tier3.py](src/tier3.py) uses adaptive BWT seeding and anchor scanning for large arrays. Before this review, its constructor always received the default 100–100,000 bp interval rather than the CLI interval.

- A fresh post-fix Chr4 run completed with **4,827 calls in 177.52 seconds**:

  - Tier 1: 1,314 final calls
  - Tier 2: 3,508
  - Tier 3: 1
  - Satellite: 4
  - No all-`N` satellite motif

  The existing untracked [bwtandem_Chr4_v3.bed](results/bwtandem_Chr4_v3.bed) contains 4,826 calls, including six satellite calls and an erroneous all-`N` call at `Chr4:0–1005`.

- The existing Chr4 Venn values are internally reproducible from the v3 BED:

  - bwtandem/TRF: 2,879 shared, 2,473 bwtandem-only, 990 TRF-only
  - all four tools: 1,891 shared bins

  However, [venn_compare.py](scripts/venn_compare.py) previously selected an older 6,987-call BED rather than the v3 input used for the documented table.

- Oropetium artifacts independently confirm:

  - 625 contigs and 243,174,629 bp
  - 56,854 calls
  - 362 canonical AAACCCT calls on exactly 18 contigs
  - 357 period-155 calls on 53 contigs
  - 13 period-310 calls
  - 3 period-465 calls
  - 689.56-second runtime

  These data and their run script are all under the untracked [oropetium](oropetium) directory.

- The ColCEN log records 32,509 calls and 1,582.84 seconds, not the 32,439 calls and approximately 12 minutes formerly stated in `TRASH.md`.

- The stored TRF ColCEN log shows TRF 4.10.0-rc.2 remaining on sequence 1 until scheduler cancellation. It does not establish why TRF failed to finish.

- The synthetic matching procedure permits compatible periods within ±20% or integer multiples, and one prediction may satisfy multiple adjacent truth records. This is useful for regression testing but permissive for a publication-grade accuracy estimate.

## 3. Issues found

- **CRITICAL — Published genome benchmarks are stale after a correctness fix.** The satellite scanner counted ambiguous bases as matches, allowing long `N` blocks to appear highly periodic. The stored Chr4 v3 BED contains a confirmed false all-`N` satellite call. Chr4, ColCEN, telomere, TRASH, Oropetium, and downstream Venn results should be regenerated from one frozen commit.

- **CRITICAL — ColCEN coverage used mixed coordinate systems.** [analyze_centromere.py](ath_cen/bench_results/analyze_centromere.py) read 1-based inclusive GFF3 coordinates as if they were 0-based half-open BED coordinates. The reported 98.2% coverage must be recomputed.

- **HIGH — Manuscript evidence is missing from version control.** `TRASH.md`, `TRASH/`, `docs/`, `oropetium/`, the Chr4 v3 BED, and its log are untracked. A fresh clone cannot reproduce or inspect central claims.

- **HIGH — Accelerator-free behavior is incomplete.** Several [accelerators.py](src/accelerators.py) fallbacks return `None` or empty results. Cython is therefore a correctness dependency for the reported three-tier pipeline, not merely an optional speed enhancement.

- **HIGH — Validation does not estimate real-genome precision.** Chr4 reports only call counts and inter-tool overlap; more calls than TRF are not evidence of higher sensitivity. ColCEN measures annotation coverage but not false-positive bases outside annotations. Oropetium confirms expected motif lengths and contig counts but lacks a full truth set.

- **HIGH — Satellite post-processing remains methodologically risky.** A few sampled windows can cause an entire uncovered block of up to 100 kb to be called repetitive. Although ambiguous-only blocks are now rejected, the approach can still bridge heterogeneous or nonperiodic sequence. Publication use needs negative controls and false-positive analysis outside annotated satellites.

- **HIGH — Benchmark inputs and software are insufficiently frozen.** Dependency versions are unpinned; external tool binaries, exact commands, hardware, thread counts, memory, checksums, and random-data provenance are not captured in a single benchmark manifest. Several run scripts contain user-specific absolute scratch paths.

- **HIGH — Missing standard research-software release material.** There is no root license, citation file, versioned package metadata, changelog/release tag, archival DOI, manuscript-ready methods document, or automated CI configuration.

- **MEDIUM — Accuracy matching is permissive.** Integer-multiple/±20% period matching and allowing a merged prediction to satisfy multiple truth records can inflate apparent sensitivity. Exact-boundary and one-to-one matching results should also be reported.

- **MEDIUM — Benchmark comparisons are not fully parameter-matched.** ULTRA's long-repeat sensitivity is discussed while noting its default maximum period is 100 bp, despite truth periods reaching 1,000 bp. Tool parameters should be harmonized where possible and sensitivity analyses reported where they cannot be.

- **MEDIUM — Some conclusions were unsupported.** The repository does not demonstrate that TRF's cancellation was caused by "combinatorial explosion," that other tools failed specifically because of divergence thresholds, or that all remaining uncovered CEN180 units have random-level autocorrelation.

- **MEDIUM — TRASH comparison documentation was inconsistent with logs.** Its total count and runtime were wrong, yielding an incorrect approximately 11× speed claim.

- **MEDIUM — Oropetium validation is only partially reproducible.** Its totals match the BED, but the reference publication comparison, acquisition provenance, and claims such as the significance of contig 003 are not implemented as a checked analysis script.

- **MEDIUM — Multiprocess failures formerly produced partial output.** Per-contig exceptions were printed and ignored, allowing a zero-exit partial genome result.

- **MEDIUM — TRF-format coordinates were incorrect.** `to_trf_dat()` emitted the internal zero-based start, although TRF `.dat` uses 1-based coordinates.

- **LOW — Test documentation and code had drifted.** The stress test uses 20 kb sequences, not 50 kb, and its reported metrics were stale.

- **LOW — Pytest compatibility warning.** Five class-scoped fixtures are instance methods and should be changed before pytest 10.

- **LOW — Repository hygiene.** A tracked compiled `__pycache__` file is modified; the working tree contains nested `.git` metadata, binaries, archives, large raw genomes, backups, and scratch artifacts. These need an explicit data-archiving policy rather than indiscriminate inclusion.

## 4. Fixes applied

- [finder.py](src/finder.py)

  - Wired the requested CLI period interval into Tier 3.
  - Avoided constructing Tier 3 when the requested interval ends below 100 bp.
  - Restricted satellite periods to the requested interval.
  - Excluded ambiguous bases from satellite autocorrelation and required at least 80% valid A/C/G/T comparisons.
  - Rejected invalid tier selections instead of silently running all tiers.

- [main.py](src/main.py)

  - Added positive and ordered period validation.
  - Made multiprocess chromosome failures abort the run before partial output is written.

- [models.py](src/models.py)

  - Corrected TRF `.dat` start coordinates to 1-based.

- [test_pipeline_constraints.py](tests/test_pipeline_constraints.py)

  - Added tests for Tier 3 range wiring, disabled out-of-range Tier 3, invalid tier names, ambiguous satellite rejection, and TRF coordinates.

- [environment.yml](environment.yml)

  - Added pytest, matplotlib, and matplotlib-venn for the documented testing and plotting workflows.

- [venn_compare.py](scripts/venn_compare.py)

  - Changed the input from the legacy `bwt_Chr4.bed` to the documented v3 output.

- [README.md](README.md) and [CLAUDE.md](CLAUDE.md)

  - Corrected Tier 1's architecture description.
  - Documented Cython as required for complete detection.
  - Corrected stress-test size and metrics.
  - Marked genome benchmarks as requiring rerun.
  - Removed unsupported causal explanations.

- [TRASH.md](TRASH.md)

  - Corrected 32,439 to 32,509 calls.
  - Corrected runtime from approximately 12 to approximately 26 minutes.
  - Corrected the speed comparison from approximately 11× to approximately 5×.

- [analyze_centromere.py](ath_cen/bench_results/analyze_centromere.py)

  - Converted GFF3 starts to 0-based half-open coordinates.

Verification completed:

- 34/34 tests pass.
- Cython and all four ctypes libraries build/load.
- Synthetic four-tool results reproduce.
- Deterministic stress test reproduces.
- Fresh full Chr4 run completed without the all-`N` satellite false call.
- `git diff --check` passes.

## 5. Remaining open items requiring the author's decision before submission, prioritized

1. **Freeze a release candidate and rerun every manuscript benchmark.** At minimum: synthetic suite, random stress test, Chr4, ColCEN, telomeres, TRASH comparison, Oropetium, and Venn plots. Preserve command, commit hash, environment lock, stdout/stderr, wall time, peak memory, CPU model, and output checksum.

2. **Decide how to validate specificity on real genomes.** Suitable author choices include manually curated loci, simulated chromosomes with realistic repeat/background composition, held-out annotations, or expert-reviewed stratified samples. Call counts and overlap diagrams alone should not be described as sensitivity or precision.

3. **Decide whether the satellite scanner is part of the core method or an optional specialized mode.** If retained by default, validate whole-block boundary accuracy, non-centromeric false positives, ambiguous/masked regions, transposon-rich gaps, divergence gradients, and behavior across several satellite families.

4. **Design a fair external-tool benchmark.** Freeze exact versions and harmonize period, copy-number, mismatch, and array-length limits. Report both each tool's recommended settings and matched-range settings.

5. **Strengthen evaluation criteria.** Add one-to-one matching, exact/near-exact boundary accuracy, motif-length accuracy, copy-number error, base-level precision/recall, and stratification by period, divergence, indels, copy count, array length, masking, and GC content.

6. **Place manuscript evidence under reproducible version control or archival storage.** The current untracked Oropetium, TRASH, docs, and v3 result trees cannot remain the sole evidence. Large inputs should use checksummed download manifests or an external DOI rather than an unreviewed bulk commit.

7. **Create publication-ready methods and software metadata.** Add a license, `CITATION.cff`, semantic version, release tag, dependency lock/container digest, installation package metadata, benchmark driver, data-accession table, and archival DOI.

8. **Resolve repository hygiene deliberately.** Decide which nested TRASH repository, archives, binaries, raw genomes, backups, caches, scratch logs, and figures belong in the software repository versus a data archive. Do not delete them until their provenance and licensing are settled.

9. **Correct the remaining pytest deprecations and add CI.** CI should test a native build, the full synthetic suite, CLI validation, and at least one accelerator-availability check.

10. **Have all biological interpretations reviewed against primary sources.** In particular, verify the Oropetium published counts and contig claims, ColCEN telomere interpretations, CEN180 divergence statements, and explanations of competing-tool behavior before transferring them into a manuscript.
