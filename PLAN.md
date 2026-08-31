# bwtandem remediation plan for a submittable release

## Purpose and current-state baseline

This plan converts the current working tree into a reproducible, validated, manuscript-submittable release. It is grounded in the repository as inspected on 2026-08-27.

Current facts that constrain the sequence of work:

- The checkout is on `main` at `7a80ed7` (`docs: add Oropetium thomaeum cross-species validation results`). Despite the surrounding project context describing that commit as unpushed, local Git currently shows `origin/main` at the same object; this does **not** make the target repository correct because `origin` is `git@github.com:framazan/bwtandem.git`, an unrelated third-party repository.
- The intended publication repository is `wyim-pgl/bwt-algorithm`. No later phase may assume that commits, tags, CI, or releases are visible there until repository ownership/remotes are deliberately corrected and verified.
- The working tree contains uncommitted audit fixes in `src/finder.py`, `src/main.py`, `src/models.py`, `scripts/venn_compare.py`, `environment.yml`, `README.md`, and `CLAUDE.md`, plus the new untracked `tests/test_pipeline_constraints.py`. A tracked `src/c_extensions/__pycache__/build.cpython-311.pyc` is modified.
- Important evidence remains untracked: `REVIEW.md`, `TRASH.md`, `TRASH/`, `docs/`, `oropetium/`, `results/bwt_chr4_v3.log`, and `results/bwtandem_Chr4_v3.bed`. `resume.md` is also untracked and needs classification.
- `README.md` documents a Python 3.11 micromamba workflow, while `environment.yml` is named `bwt` and requests Python 3.13. Dependencies and external tools are not pinned.
- `README.md` correctly warns that the current Chr4 and Col-CEN values predate the ambiguous-base satellite fix. The stored `results/bwtandem_Chr4_v3.bed` is legacy evidence, not an acceptable result for the release candidate.
- The audit statement “34/34 tests pass” is currently **CANNOT VERIFY**. The re-review environment exposed neither `pytest` nor `conda`; the claim must be independently reproduced and logged in the release environment.

### Governing dependency chain

The mandatory serial path is:

1. establish the correct repository/remote and preserve the current work;
2. decide the scientific contracts and benchmark protocol;
3. implement and test all detector, coordinate, validation, and tooling fixes;
4. commit those fixes and freeze a release-candidate commit plus environment;
5. rerun benchmarks only from that frozen commit;
6. review results and claims, then tag/archive the release and submit.

No full-genome benchmark should be launched before Phase 3 exits. Any earlier run is diagnostic only and must not be used in the manuscript.

## Phase 0 — Establish the authoritative repository and evidence boundary

**Tracker scope:** #9, #22

### Work

1. Preserve and classify the current dirty tree before changing repository configuration.
   - Capture `git status --short`, `git diff --binary`, `git log --oneline --decorate -20`, and `git remote -v` in a dated, external recovery bundle.
   - Confirm that the intended target is `wyim-pgl/bwt-algorithm` and verify its default branch and current tip.
   - Compare the local `7a80ed7` history with the target repository before choosing how to publish it.

2. **REQUIRES AUTHOR DECISION — remote topology.** Choose one:
   - replace `origin` with `wyim-pgl/bwt-algorithm` and optionally retain `framazan/bwtandem` as a read-only `upstream`; this is simpler for routine pushes but changes the established local name; or
   - keep the third-party remote under a renamed `upstream` and add the project repository as `origin`; this preserves provenance explicitly; or
   - create a fresh clone from `wyim-pgl/bwt-algorithm` and apply the reviewed changes there; this is safest if histories differ, but requires careful patch transfer.
   - Do not push to `framazan/bwtandem`. Record the selected topology in contributor documentation.

3. Separate source-controlled material from archived data and generated results.
   - Commit source/tests/small documentation that are part of the method: relevant content from `docs/`, `tests/test_pipeline_constraints.py`, benchmark scripts, provenance manifests, and reviewed portions of `TRASH.md`/Oropetium analysis.
   - Keep large/raw or third-party material such as `oropetium/Othomaeum_386_v1.0.fa`, `oropetium/download.20260401.094229.zip`, Col-CEN genomes, external-tool binaries, and large result sets outside ordinary Git unless licensing and repository policy explicitly allow them.
   - Create a checksummed data manifest, proposed path `benchmarks/data_manifest.tsv`, with source URL/accession, license, retrieval date, byte size, and SHA-256.
   - Create `benchmarks/README.md` describing which outputs are generated, archived, or versioned.
   - Add/update `.gitignore` for `__pycache__/`, compiled extension outputs, transient logs, scheduler files, scratch archives, and regenerated benchmark outputs.
   - Remove the tracked bytecode file `src/c_extensions/__pycache__/build.cpython-311.pyc` from version control without treating it as source.
   - Classify `REVIEW.md` as audit documentation, `resume.md` as local/session material, and `TRASH/` as either licensed benchmark evidence or external archive content.

4. Preserve the audit fixes as a reviewable commit series on the authoritative repository, but do not call it the release candidate yet.
   - First commit the already-reviewed functional fixes and `tests/test_pipeline_constraints.py` together or in narrowly scoped commits.
   - Keep legacy benchmark outputs explicitly labeled as pre-fix; do not overwrite them until Phase 4.

### Commands to execute when authorized

```bash
git remote -v
git status --short
git diff --check
git diff --binary > /safe/external/location/bwtandem-working-tree-2026-08-27.patch
git log --oneline --decorate -20
git ls-files | sort
sha256sum <benchmark-inputs> > benchmarks/input_checksums.sha256
```

Remote edits, commits, and pushes are deliberately planning steps only here; the author/repository owner must execute them after choosing the topology.

### Parallelism

- Data/provenance inventory, licensing review, and `.gitignore` design can proceed in parallel.
- Remote correction/history reconciliation is serial and must precede publishing any commit.
- Classification of untracked artifacts must finish before committing them or uploading them to an archive.

### Deliverables

- Correct, verified connection to `wyim-pgl/bwt-algorithm`, with the third-party remote incapable of receiving accidental project pushes.
- Recovery bundle for the starting dirty tree.
- Committed audit fixes and `tests/test_pipeline_constraints.py` on the authoritative repository.
- `benchmarks/data_manifest.tsv`, `benchmarks/README.md`, input checksum file, and repository hygiene policy.
- Clean Git status except for explicitly documented external/generated artifacts.

### Exit criterion

The authoritative repository contains the current code fixes and regression tests in reviewable commits; the remote target is verified; every untracked path is either committed, ignored, or represented in a checksummed archive manifest; and `git status --short` is clean before Phase 1 decisions are finalized.

## Phase 1 — Freeze scientific claims and evaluation contracts

**Tracker scope:** #11, #16, #17

### Work

1. Define what “correct” means on real genomes (#11).
   - **REQUIRES AUTHOR DECISION — real-genome specificity truth.** Select and justify one or more of: manually curated loci, held-out annotations, realistic simulation embedded in real background, or an expert-reviewed stratified sample. Manual review gives biological credibility but limited scale; simulation gives exact truth but may not capture genome complexity; held-out annotations scale but inherit annotation bias.
   - Specify call-level and base-level precision/recall, boundary tolerance, motif-period accuracy, copy-number error, stratification, and confidence intervals.
   - Explicitly prohibit interpreting raw call counts or Venn overlap as sensitivity/precision.
   - Define negative controls for centromeric and non-centromeric sequence, masked/ambiguous sequence, transposon-rich regions, and heterogeneous satellites.

2. Define fair external-tool protocols (#16), coordinating with #3, #4, and #5.
   - **REQUIRES AUTHOR DECISION — comparator policy.** Choose whether the main table uses each tool's recommended defaults, a harmonized detectable range, or both. Defaults reflect user experience but can be range-incomparable; harmonized settings improve methodological fairness but may be nonstandard for a tool. Preferred publication practice is to report both as separately named analyses if computationally feasible.
   - Record exact versions and commands for bwtandem, TRF, mreps, ULTRA, TRASH, and any tools used by #3/#4/#5.
   - Define matched min/max period, minimum copies, array length, mismatch/indel limits, masking, and output normalization.
   - Specify treatment of truth outside a tool's supported range rather than counting it as an unexplained false negative.

3. Audit every comparative, causal, and biological statement (#17).
   - Review `README.md`, `CLAUDE.md`, `code_update.md`, `software_comparison.md`, `TRASH.md`, and the manuscript draft.
   - Maintain a claim-evidence table, proposed `docs/claim_evidence_matrix.md`, linking each retained claim to a script, frozen result, or primary source.
   - Mark historical results as historical and remove unsupported explanations of competitor behavior, biological significance, or failure cause.

### Commands to execute

```bash
rg -n "sensitivity|precision|outperform|faster|slower|failed|divergence|combinatorial|significant|highest" \
  README.md CLAUDE.md code_update.md software_comparison.md TRASH.md docs/
python3 scripts/compare_all_tools.py --help  # after Phase 3 adds a real CLI, if adopted
```

### Parallelism

- Real-genome truth design, comparator parameter inventory, and prose/claim audit can run in parallel.
- Final metrics and claims are serial dependencies for Phase 3 evaluation-code design and Phase 4 benchmark execution.

### Deliverables

- `docs/validation_protocol.md` defining truth sets, negative controls, metrics, matching, stratification, and acceptance thresholds.
- `benchmarks/tool_parameters.tsv` with exact versions, settings, supported ranges, and rationale.
- `docs/claim_evidence_matrix.md` and corrected manuscript-facing prose proposals.
- Written decisions resolving the truth-set and comparator-policy choices.

### Exit criterion

The authors approve a versioned validation protocol and comparator matrix; every planned manuscript claim has an identified evidentiary test; and #3, #4, and #5 have explicit metric/parameter requirements ready for implementation.

## Phase 2 — Correct detector, coordinate, and test-control behavior

**Tracker scope:** #8, #10, #12, #19, #20, #21, #23

### Work

1. Correct Col-CEN coordinate handling before any Col-CEN rerun (#8; required by #5).
   - Normalize GFF3, BED, TRF, mreps, and ULTRA inputs to 0-based half-open intervals in `ath_cen/bench_results/analyze_centromere.py`.
   - Replace the residual inclusive arithmetic at current lines 41, 127, 144, and 148.
   - Normalize TRF/mreps starts rather than comparing them directly with BED coordinates.
   - Add a dedicated test module, proposed `tests/test_centromere_coordinates.py`, covering single-base intervals, adjacent/non-overlapping intervals, exact boundary overlap, chromosome start, and equivalent records in every tool format.

2. Establish a safe accelerator contract (#10).
   - **REQUIRES AUTHOR DECISION — native accelerator policy.** Choose either:
     - hard requirement: fail at startup with an actionable error when `_accelerators` is unavailable; simplest and least likely to produce scientifically incomplete output; or
     - complete pure-Python fallback: implement and test behaviorally equivalent fallbacks; more portable but potentially slow and substantially more implementation work.
   - Update `src/accelerators.py`, CLI startup in `src/main.py`, `README.md`, `CLAUDE.md`, and container/package tests to match the selected contract.
   - Remove debug-print-on-import behavior in favor of a controlled exception or structured warning consistent with the contract.

3. Constrain and test satellite filling/merging (#12; feeds #4 and #5).
   - **REQUIRES AUTHOR DECISION — satellite scanner product role.** Choose:
     - core/default method: must pass whole-block boundary, real-genome specificity, masked-base, transposon, heterogeneous-period, divergence-gradient, and cross-family validation; or
     - optional specialized mode: expose an explicit CLI flag/preset, keep general tandem-repeat results independent, and report satellite benchmarks separately. Optional mode reduces default false-positive risk but complicates headline method scope.
   - Apply the same A/C/G/T validity mask and coverage threshold to `_merge_adjacent_repeats()` that `_fill_satellite_gaps()` currently applies.
   - Replace “three windows imply the whole block” with interval-wide support or a segmentation method, or document and validate a conservative maximum expansion rule.
   - Add tests in `tests/test_pipeline_constraints.py` or a new `tests/test_satellite_postprocessing.py` for all-`N`, mixed `N`/valid, `N`-rich bridge, transposon-like gap, heterogeneous block, valid CEN180-like block, period bounds, and selected-mode behavior.

4. Reject any invalid tier token (#19).
   - Change tier parsing/normalization in `src/main.py` and `src/finder.py` so `tier1,typo` fails rather than silently running Tier 1.
   - Extend `tests/test_pipeline_constraints.py` with mixed-valid/invalid CLI and constructor cases.

5. Make stress failures fatal and thresholds scientifically explicit (#20).
   - Refactor `tests/test_random_stress.py` so any detector exception fails the run, all attempted cases remain in accounting, and temporary cleanup uses `try/finally` without suppression.
   - Add a regression test that injects a failure and proves the command exits nonzero.
   - **REQUIRES AUTHOR DECISION — stress acceptance thresholds.** Approve publication/release thresholds by tier and overall. The current 80% sensitivity/50% precision pass gate is a smoke-test threshold, not automatically a manuscript accuracy standard; distinguish CI regression gates from claimed benchmark performance.

6. Remove pytest compatibility warnings (#21).
   - Convert the five class-scoped instance fixtures in `tests/test_ground_truth.py` to supported module/class fixtures without `self`.
   - Test against the current pinned pytest and the newest supported pytest.

7. Correct Venn half-open binning (#23).
   - In `scripts/venn_compare.py`, derive the last occupied bin from `end - 1`, validate nonempty intervals, and add tests for exact 500 bp boundaries.
   - Do not regenerate the figure yet; Phase 4 regenerates it from the frozen release candidate.

### Commands to execute

```bash
micromamba run -n bwtandem python3 -m pytest \
  tests/test_pipeline_constraints.py \
  tests/test_centromere_coordinates.py \
  tests/test_satellite_postprocessing.py -v
micromamba run -n bwtandem python3 -m pytest tests/test_random_stress.py -v
micromamba run -n bwtandem python3 -m pytest tests/ -v -W error::pytest.PytestRemovedIn9Warning
python3 -m src.main tests/fixtures/synth_mixed.fa --tiers tier1,typo -o /tmp/must-not-exist
git diff --check
```

Use actual available pytest warning classes if the pinned version names them differently; the acceptance requirement is zero fixture deprecation warnings.

### Parallelism

- Coordinate/parser fixes, accelerator contract, tier validation, stress-test control flow, fixture modernization, and Venn binning can be implemented in parallel in separate branches.
- Satellite redesign depends on the Phase 1 role/validation decision but can otherwise proceed alongside those mechanical fixes.
- Integration and the full test run are serial after all branches merge.
- Col-CEN reruns are strictly blocked until #8 is merged and boundary-tested.

### Deliverables

- Corrected `ath_cen/bench_results/analyze_centromere.py`, `src/accelerators.py`, `src/main.py`, `src/finder.py`, `scripts/venn_compare.py`, `tests/test_random_stress.py`, and `tests/test_ground_truth.py`.
- New/extended coordinate, satellite, tier-selection, stress-failure, accelerator-availability, and bin-boundary tests.
- Updated `README.md`/`CLAUDE.md` describing the chosen accelerator and satellite contracts.
- Passing targeted and full test logs from the development environment.

### Exit criterion

All seven issues are merged into the authoritative repository; coordinate boundary tests pass for every Col-CEN input format; no analysis can silently run with incomplete accelerators, invalid tiers, or skipped stress cases; satellite behavior matches the approved contract; pytest fixture warnings are gone; and no benchmark rerun has yet been used as release evidence.

## Phase 3 — Build the reproducible evaluation harness and freeze the release candidate

**Tracker scope:** #13, #15, #18

### Work

1. Create a locked, documented execution environment (#13).
   - Reconcile `environment.yml` with the documented Python 3.11 environment, rename it consistently, and pin compatible versions.
   - Prefer an explicit lock file and/or immutable Singularity/Apptainer container digest. Ensure the container builds `_accelerators.pyx` and all standalone `src/c_extensions/*.c` libraries.
   - Add a benchmark manifest/driver, proposed paths `benchmarks/manifest.yaml` and `scripts/run_benchmarks.py`, that records commit, dirty status, commands, input/output checksums, tool versions, threads, CPU, memory, wall time, and seeds.
   - Refactor user-specific absolute paths in Col-CEN/Oropetium/TRASH scripts into arguments or manifest entries.

2. Replace permissive-only matching with publication metrics (#15).
   - Keep current `periods_compatible()`/many-to-one matching only as a clearly named regression metric.
   - Add one-to-one assignment, exact and tolerance-based boundary metrics, base-level precision/recall, motif-period error, and copy-number error in a reusable evaluation module, proposed `scripts/evaluate_repeats.py`.
   - Add unit tests with adjacent truths, one merged prediction, duplicate predictions, integer-multiple motifs, near-period mismatches, and no-overlap cases.
   - Update `scripts/compare_all_tools.py` to emit machine-readable TSV/JSON plus manuscript tables.

3. Make Oropetium validation checked and reproducible (#18).
   - Replace the hard-coded scheduler-only workflow in `oropetium/run_bwtandem.sh` with manifest-driven inputs/outputs.
   - Add `scripts/analyze_oropetium.py` to reproduce every count in `README.md`: assembly size/contigs, total calls, AAACCCT calls and contigs, and period-155/310/465 counts.
   - Record accession/publication provenance and checksums in `benchmarks/data_manifest.tsv`.
   - Distinguish mechanically reproduced counts from biological interpretations requiring a primary-source check.

4. Independently reproduce the test claim before freezing.
   - Create the environment from the committed lock/container, build native extensions, and run the complete suite.
   - Capture exact commands, `python --version`, `pip/conda list --explicit` or lock hash, compiler version, native library load checks, stdout/stderr, exit status, and duration.
   - Run the random stress command independently; it is not collected as an ordinary pytest test in the current README workflow.
   - The historical “34/34 tests pass” statement remains **CANNOT VERIFY** until this evidence exists. If the corrected suite has a different test count, report the new collected count rather than preserving 34 as a target.

5. Freeze the release-candidate commit only after all Phase 0–3 changes and tests are committed.
   - Require clean `git status`, reviewed commits on `wyim-pgl/bwt-algorithm`, green local/CI tests, and a recorded commit SHA.
   - Use that exact SHA and immutable environment identifier for every Phase 4 run. Any detector/evaluator change after freeze invalidates affected benchmark runs and requires a new candidate SHA.

### Commands to execute

```bash
micromamba create -f environment.yml
micromamba run -n bwtandem python --version
micromamba run -n bwtandem python -m pip freeze
micromamba run -n bwtandem python src/c_extensions/build.py
micromamba run -n bwtandem python -c "import src._accelerators"
micromamba run -n bwtandem python -m pytest tests/ -v --junitxml=artifacts/tests/junit.xml \
  2>&1 | tee artifacts/tests/pytest.log
micromamba run -n bwtandem python -m tests.test_random_stress \
  2>&1 | tee artifacts/tests/random_stress.log
git status --porcelain
git rev-parse HEAD
sha256sum environment.yml <lock-file> > artifacts/release_candidate_environment.sha256
```

Adjust the environment name only if the committed environment adopts a different approved name; commands and documentation must agree.

### Parallelism

- Environment locking, matching implementation, Oropetium analyzer development, and benchmark-driver development can proceed in parallel after Phase 2 interfaces stabilize.
- Full acceptance tests are serial after integration.
- Release-candidate commit/environment freeze is strictly serial after all Phase 0–3 work is committed and accepted.

### Deliverables

- Pinned `environment.yml`, lock file and/or immutable container recipe/digest.
- `benchmarks/manifest.yaml`, benchmark driver, data/tool manifests, and machine-readable provenance capture.
- `scripts/evaluate_repeats.py`, updated `scripts/compare_all_tools.py`, and matching metric tests.
- `scripts/analyze_oropetium.py` plus Oropetium provenance and analysis tests.
- `artifacts/tests/pytest.log`, JUnit output, stress-test log, environment inventory, native-build/load log, and checksums.
- A clean, immutable release-candidate commit SHA on `wyim-pgl/bwt-algorithm`.

### Exit criterion

A clean release-candidate commit and locked environment reproduce the complete test suite with no unexpected skips or deprecation warnings, native components build/load, stress cases cannot disappear, all evaluation scripts are tested, and the exact test count/result is independently documented. Only now may full-genome reruns begin.

## Phase 4 — Execute frozen benchmarks and validate manuscript results

**Tracker scope:** #7

**Legacy benchmark integration:** #3 (Human GRCh38), #4 (Maize Mo17 satellite/CentC), and #5 (Arabidopsis Col-CEN) execute here under the Phase 1 protocol and Phase 3 frozen release candidate. They are not part of the #7–#23 uniqueness count; they are required benchmark workstreams folded into this phase.

### Work

1. Launch long-wall-clock benchmarks as soon as the Phase 3 gate opens.
   - **LONG WALL CLOCK — #3 Human GRCh38:** run bwtandem and approved comparators; collect recall/precision, memory, runtime, stratified metrics, and artifacts required by #3.
   - **LONG WALL CLOCK — #4 Maize Mo17 satellite/CentC:** run the approved satellite mode and fair comparator configuration; replace boundary artifacts with the Phase 1 metric.
   - **LONG WALL CLOCK — #5 Arabidopsis Col-CEN:** run only after #8 is fixed and boundary-tested; use the corrected analyzer and v2/approved truth reconciliation from #5.
   - **LONG WALL CLOCK — Oropetium:** regenerate the BED from the release candidate, then run `scripts/analyze_oropetium.py`.
   - Run Chr4, telomere, TRASH, synthetic five-fixture comparison, and deterministic random stress benchmarks from the same commit/environment.

2. Run independent workstreams in parallel where compute and licenses permit.
   - GRCh38, Mo17, Col-CEN, Oropetium, and Chr4 are independent after the freeze and should be submitted concurrently to separate scheduler jobs.
   - External tools for a given dataset may run in parallel if each receives the frozen input and parameter manifest.
   - Analysis for a dataset waits for all required tool outputs and checksum validation.

3. Enforce provenance and failure rules.
   - Every job must record the release-candidate SHA, environment/container digest, command, version, hostname/CPU, allocated threads/memory, peak RSS, wall time, exit code, stdout/stderr, input checksum, and output checksum.
   - A timeout, cancellation, parser error, missing output, or partial chromosome is a failed/missing result, never a zero-call result or causal explanation.
   - Do not hand-edit generated tables. Generate TSV/JSON first, then render Markdown/figures.

4. Regenerate stale artifacts (#7).
   - Replace or version `results/bwtandem_Chr4_v3.bed`, `results/bwt_chr4_v3.log`, and `results/venn_comparison_chr4.png` with release-candidate outputs; preserve legacy files only under an explicitly historical location.
   - Regenerate Col-CEN outputs under `ath_cen/bench_results/` or archive them with a manifest rather than mixing them with pre-fix files.
   - Regenerate Oropetium and TRASH results, telomere analyses, synthetic comparison, random stress results, and README tables.
   - Confirm the all-`N` Chr4 satellite call is absent and run masked/ambiguous negative controls.

5. Review results against acceptance criteria.
   - Apply one-to-one, boundary, base-level, period, and copy-number metrics from Phase 3.
   - Report defaults and matched-range comparator results separately per the author decision.
   - Investigate regressions before changing acceptance thresholds. Any code change creates a new release candidate and requires rerunning affected datasets.

### Representative commands

Exact invocations must come from the committed manifest; the command pattern is:

```bash
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset synthetic
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset chr4
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset colcen
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset oropetium
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset mo17
python3 scripts/run_benchmarks.py --manifest benchmarks/manifest.yaml --dataset grch38
python3 scripts/venn_compare.py
python3 scripts/analyze_oropetium.py --manifest benchmarks/manifest.yaml
sha256sum <all-inputs-and-outputs> > artifacts/benchmarks/checksums.sha256
```

If a unified driver is not implemented, the manifest must still contain literal executable commands with equivalent provenance capture.

### Parallelism

- Launch GRCh38, Mo17, Col-CEN, Oropetium, and Chr4 concurrently immediately after Phase 3 exits; these are the principal long-wall-clock tasks.
- Synthetic/stress, telomere analysis, TRASH analysis, and documentation table-generation can run alongside full-genome jobs.
- Per-dataset analysis and manuscript updates are serial after that dataset's runs finish.
- Final cross-dataset conclusions wait for all required workstreams.

### Deliverables

- Complete, checksummed outputs/logs for synthetic, stress, Chr4, Col-CEN, telomere, TRASH, Oropetium, Mo17, and GRCh38 workstreams.
- Machine-readable metrics and regenerated figures/tables.
- Updated #3, #4, and #5 with commands, artifacts, and conclusions supported by their revised metrics.
- Updated `README.md`, `TRASH.md`, benchmark documentation, and claim-evidence matrix using release-candidate results only.
- Explicit record of failures/timeouts and excluded claims.

### Exit criterion

Every required benchmark completed or is transparently reported as unavailable; all reported values derive from the same frozen commit/environment; #3, #4, and #5 meet their accepted metric requirements; stale pre-fix values are removed from manuscript-facing material; and every retained claim maps to checksummed evidence.

## Phase 5 — Package, archive, and approve the submission release

**Tracker scope:** #14

### Work

1. Add research-software release material.
   - **REQUIRES AUTHOR DECISION — license.** Select an author/institution-approved license. A permissive license improves reuse; a copyleft license requires derivative distribution under matching terms. Confirm third-party code/data compatibility before selection.
   - Add root `LICENSE`, `CITATION.cff`, `CHANGELOG.md`, semantic version, and package metadata (`pyproject.toml` preferred).
   - Add a methods/reproducibility document, proposed `docs/METHODS.md`, and data-accession table.
   - Update `README.md`, `CLAUDE.md`, and installation instructions so accelerator/satellite/environment contracts are consistent.

2. Add automated CI and release checks.
   - Create `.github/workflows/tests.yml` to build native extensions, run the complete collected suite, test CLI validation, verify accelerator failure/fallback behavior, and run a small synthetic smoke benchmark.
   - Add packaging/install tests from a clean checkout and, if supported, the pinned Python version(s).
   - Add a release checklist ensuring the benchmark commit matches the proposed tag.

3. Archive the release and evidence.
   - Tag the exact benchmarked release-candidate commit; do not rebuild the release from a different SHA.
   - Upload large inputs/results/logs/manifests to the selected archival repository and obtain a DOI or stable accession.
   - Put archive identifiers and checksums in `CITATION.cff`, `README.md`, `docs/METHODS.md`, and `benchmarks/data_manifest.tsv`.
   - Ensure GitHub release artifacts and archival artifacts correspond to the same source SHA.

4. Conduct final author acceptance.
   - Authors review biological interpretations against primary sources, especially Oropetium contig claims, CEN180 divergence/telomere statements, CentC conclusions, and competitor behavior.
   - Re-run the acceptance commands from a fresh clone of `wyim-pgl/bwt-algorithm` at the proposed tag.

### Commands to execute

```bash
git clone https://github.com/wyim-pgl/bwt-algorithm.git /fresh/path/bwt-algorithm
git -C /fresh/path/bwt-algorithm checkout <proposed-tag>
micromamba create -f /fresh/path/bwt-algorithm/environment.yml
micromamba run -n bwtandem python -m pytest tests/ -v
micromamba run -n bwtandem python -m tests.test_random_stress
python -m build
git status --porcelain
git rev-parse HEAD
sha256sum -c artifacts/benchmarks/checksums.sha256
```

### Parallelism

- License/institutional approval, citation metadata, packaging, methods writing, and CI implementation can proceed in parallel while Phase 4 jobs run, but must use placeholders until final results/DOIs exist.
- Final tag, archive publication, fresh-clone verification, and submission approval are serial after Phase 4 acceptance.

### Deliverables

- `LICENSE`, `CITATION.cff`, `CHANGELOG.md`, `pyproject.toml`, `docs/METHODS.md`, data-accession table, and `.github/workflows/tests.yml`.
- Green CI and fresh-clone verification logs.
- Semantic version tag, GitHub release, archived benchmark/data package, checksums, and DOI/stable accession.
- Final claim-evidence matrix signed off by the authors.

### Exit criterion

The exact benchmarked commit is tagged and archived; licensing/citation/package/CI material is complete; a fresh clone reproduces installation, native build, the full test suite, stress test, and checksum validation; all manuscript claims are supported by archived evidence; and the authors approve the biological interpretation and final manuscript tables.

## Verification and manuscript-submission acceptance

The manuscript must not be submitted until all of the following evidence exists:

- [ ] The canonical source and release tag are accessible at `wyim-pgl/bwt-algorithm`; no release action depends on or targets `framazan/bwtandem`.
- [ ] The release tag resolves to the exact commit used for every manuscript benchmark, and the benchmark working tree was clean.
- [ ] The current historical claim “34/34 tests pass” has been independently reproduced or replaced with the actual corrected-suite count. Until a documented run exists, it remains **CANNOT VERIFY**.
- [ ] The test record includes exact commands, Python version, dependency/lock versions, compiler/native build details, environment or container digest, stdout/stderr log, exit status, duration, JUnit output, and unexpected-skip count.
- [ ] Native extensions and all standalone C libraries build and load in the documented environment; unavailable-accelerator behavior matches the approved contract.
- [ ] The deterministic random stress test runs all planned cases, reports failures in its denominator, exits nonzero on any case exception, and meets the author-approved release thresholds.
- [ ] Col-CEN coordinate unit tests cover all parsers and boundary cases, and the corrected analyzer—not a legacy output—produces the manuscript values.
- [ ] Chr4, Col-CEN, Oropetium, TRASH, telomere, Mo17, and GRCh38 results have been regenerated from the frozen commit where required; synthetic and stress results have also been regenerated.
- [ ] Long-wall-clock job logs include command, version, hardware, threads, peak memory, runtime, exit status, and checksums; cancellation or timeout is not reported as a biological/tool result.
- [ ] Real-genome validation reports specificity/false positives as well as recall/coverage, using the approved Phase 1 truth design.
- [ ] Comparative tables clearly distinguish recommended-default from matched-range runs and do not penalize tools for undisclosed unsupported ranges.
- [ ] Accuracy evaluation includes one-to-one matching, boundary quality, base-level metrics, motif-period accuracy, and copy-number error; permissive many-to-one matching is labeled as secondary.
- [ ] Satellite calls pass masked/ambiguous, heterogeneous, transposon-rich, non-centromeric, divergence-gradient, boundary, and cross-family validation appropriate to the chosen core/optional role.
- [ ] The known all-`N` Chr4 call is absent from corrected outputs, and no stale pre-fix table or figure remains manuscript-facing.
- [ ] Every input/output has a checksum and provenance record; every large artifact has a stable archive location; every required small script/test/document exists in the release repository.
- [ ] `LICENSE`, `CITATION.cff`, package metadata, changelog, methods, data accessions, CI, semantic tag, and archival identifier are complete.
- [ ] The claim-evidence matrix contains no unsupported causal, competitive, or biological statement, and authors have reviewed interpretations against primary sources.
- [ ] A final clean-room/fresh-clone acceptance run succeeds and its logs are archived with the release.

## Issue-accounting matrix

Each audit issue appears in exactly one phase:

| Issue | Phase |
|---|---|
| #7 | Phase 4 — frozen benchmark reruns |
| #8 | Phase 2 — coordinate correctness |
| #9 | Phase 0 — evidence/version-control boundary |
| #10 | Phase 2 — accelerator contract |
| #11 | Phase 1 — real-genome validation contract |
| #12 | Phase 2 — satellite filling/merging |
| #13 | Phase 3 — environment and benchmark freeze |
| #14 | Phase 5 — release metadata, CI, archive |
| #15 | Phase 3 — publication-grade matching |
| #16 | Phase 1 — fair comparator protocol |
| #17 | Phase 1 — supported claims |
| #18 | Phase 3 — reproducible Oropetium analysis |
| #19 | Phase 2 — tier validation |
| #20 | Phase 2 — stress-test failure handling |
| #21 | Phase 2 — pytest compatibility |
| #22 | Phase 0 — repository hygiene/archive policy |
| #23 | Phase 2 — Venn bin boundaries |

Legacy benchmark issues fold into Phase 4: #3 is the GRCh38 workstream, #4 is the Mo17 satellite/CentC workstream, and #5 is the Col-CEN workstream. Their metric and parameter requirements are designed in Phase 1, their code dependencies are implemented in Phases 2–3, and their long-running executions occur only after the release-candidate freeze.
