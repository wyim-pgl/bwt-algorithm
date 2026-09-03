# Codex adversarial review — manuscript round 2 (Results vs deposited artefacts)

Verbatim output of `codex exec -m gpt-5.6-sol` at `xhigh`, run 2026-09-03 in a detached
worktree at `67f11b1`. Nine findings. The `file:line` links carry that worktree's path
prefix; strip through `review-wt/` for the repository path.

Deposited 2026-09-03 because `todo.md` cited this file while it existed only in a session
scratchpad, which made its findings unauditable (Codex plan review, finding 3).

---

I found three serious integrity failures: post-preregistered selective reporting, ground-truth-driven configuration selection, and a scorer that discards TRF’s valid period field. All deposited checksums pass; I found no evidence of artifact tampering.

1. **Severity: HIGH — Preregistered AniAnn’s metrics, and other unfavorable Arabidopsis cells, were selectively omitted**

   - **Location:** [manuscript.md:129](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:129), [manuscript.md:139](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:139), [manuscript.md:182](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:182), [manuscript.md:254](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:254), [manuscript.md:298](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:298); `docs/2026-09-01-longdust-anianns-benchmark-protocol.md:54–60`; `results/comparators2026/score_2026_human.txt:25–33`; `results/comparators2026/score_2026_colcen.txt:12`.
   - **Defect:** The preregistered protocol, committed before scoring, says AniAnn’s is scoreable for “**period-banded rules using its periodicity column where a table’s rule reads a period**.” The manuscript reverses that decision after results exist, excluding it because its periods are “array-level” and stating that the 2026 tools have no banded rows.
   - Deposited missing human cells are:

     - ≤100 bp: 114 calls, 0.01% region recall, 85.09% precision, 0.92% bp recall, 54.50% bp precision.
     - 101–2,000 bp: 205 calls, 0.03%, 81.46%, 2.69%, 10.03%.

   - More damagingly, the deposited Col-CEN log gives AniAnn’s **99.15% recall in both 150–200 and 150–400 bands**, exceeding BWTandem’s reported 91.31% and 94.68%, and exceeding TRF’s 97.55%. Section 3.2 nevertheless reports banded recall only for BWTandem, ULTRA and TRF.
   - The same scoring rule applied to the manifest’s maize AniAnn’s BED gives omitted values of:

     - knob180 100–500: 21/25, 87.40% coverage, 22,148 bp offset.
     - TR-1 100–500: 14/17, 67.98%, 83,740 bp.
     - CentC 100–200: 17/17, 80.89%, 14,917 bp.

     These all lead the corresponding BWTandem banded coverages of 64.20%, 17.33%, and 40.71%. They were not deposited as score outputs.

   - The omission is broader in Section 3.2: rescoring the manifest inputs gives tantan 61.68%/98.97%, mreps 75.28%/78.45%, and both TRASH modes 99.39%/99.39% in the two bands. The manuscript reports none of these; several are unfavorable to BWTandem.
   - **Minimal fix:** Honor the preregistration. Add AniAnn’s to every applicable period-conditioned panel, report all scoreable Table 2 banded cells, and deposit the maize and active Col-CEN scoring outputs. If array-level periods are now deemed invalid, withdraw every period-based AniAnn’s cell—including its CEN180 Count—and explicitly document the post-registration change.

2. **Severity: HIGH — The Arabidopsis configuration was selected using the evaluation truth**

   - **Location:** `manuscript.md:470,479`; `/data/gpfs/assoc/pgl/devel/exp1_human/filip_repro/catchall_experiment_results.md:12–24`; `results/manifest.tsv:25`.
   - **Defect:** Supplementary S2 assigns Arabidopsis a different configuration—catch-all off—and explains that it was disabled “**where it costs base pair precision without recovering centromeric arrays**.” The surviving experiment record explicitly scores both configurations against the same `colcen_cen180.bed` truth and concludes: “**Keep catch-all OFF for Col-CEN**.” Its measured comparison was:

     - off: 99.67% recall, 65.54% bp precision;
     - on: 99.68% recall, 60.72% bp precision.

     Thus the submitted Table 2 configuration was chosen after inspecting the target evaluation metric, specifically to avoid an unfavorable precision result. No held-out Arabidopsis satellite set exists, and no catch-all-on current-build row is deposited.
   - **Minimal fix:** Either use a configuration fixed independently of Col-CEN, or designate Col-CEN as tuning data and evaluate the selected configuration on held-out satellite truth. At minimum, regenerate and report both configurations in Table 2 and label the choice as in-sample.

3. **Severity: HIGH — TRF’s valid period column is discarded, hiding a period score that beats BWTandem**

   - **Location:** [manuscript.md:518](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:518), [manuscript.md:537](/data/gpfs/assoc/pgl/tmp/claude-3092957/-data-gpfs-assoc-pgl-devel-bwt-algorithm/ab324600-278d-496b-944b-b0f7ec02f088/scratchpad/review-wt/manuscript.md:537); `scripts/scoring/convert_to_bed.py:72–94`; `scripts/scoring/score_one_to_one.py:20–27,66–88,219–220`; `results/manifest.tsv:109,114`.
   - **Defect:** The manuscript says “**TRF’s period column is empty because its BED’s motif column holds the full array sequence**.” The converter explicitly writes `motif` to column 4 and `period` to column 5. The actual TRF BED begins with periods 6, 29, 76, 61 and 117 in column 5.
   - The one-to-one scorer is invoked with `--pred-col5 period --pred-motif-is-sequence`, but `load()` never stores column 5 when its semantics are `period`; `period_of()` only returns `len(motif)`. Because `--pred-motif-is-sequence` nulls the motif, every TRF period becomes missing.
   - Reusing the same 518,719 matched pairs while reading column 5 correctly gives TRF: **63.53% exact**, 6.97% integer-multiple, 15.31% within 20%, 14.19% outside. Its exact score therefore exceeds BWTandem’s 58.66% and ULTRA’s 59.61%.
   - **Minimal fix:** Store an explicit `period` field in the scorer, rerun S4, deposit corrected JSON, and report TRF’s period result.

4. **Severity: HIGH — The claimed “fully covered” accuracy deposit is not self-contained**

   - **Location:** `results/README.md:39–43,52–60,123–127`; `results/manifest.tsv:25–57,85–86`; `results/regen/table3bc_provenance.json:3,13–100`; `manuscript.md:105,184,189–220,242–309`.
   - **Defect:** The repository claims: “**Every recall, precision, coverage, detection-count, offset and fragmentation number … [is reproducible by] re-running the named script on the named file.**” It later admits that competitor outputs are not deposited. Specific failures are:

     - The active Table 3B/3C scorer, `rescore_tables_3bc.py`, is absent. Its provenance points to an external cluster copy, as does its dependency `score_exp3.py`.
     - Most comparator BEDs are external absolute paths. The human score logs preserve their derived cells, but a clean checkout cannot reproduce them.
     - Table 2 has no deposited score log for the active regenerated/re-run table. The 2026 log contains stale baseline rows. I could reproduce the table only because the external cluster inputs remain accessible.
     - Table 3A’s NCRF/TRASH cells and the narrative BWTandem count 475,732 have no deposited scoring output. They do recompute correctly from the currently reachable external files.
     - Historical plant runtime/memory evidence remains only in external GNU-time logs.
     - The original narrow range-cost arm’s 5 h 17 m and 35.24 GB have no manifest row: row 45 merely names job 5981960, while the active row 1b is a different job.
     - The seeding-ablation memories 3.39/1.30 GB and the 0.6 s/9 s/21 min cache timings have no deposited raw records; the latter are candidly labelled unlogged.

   - **Trace census:** Human Tables 1a–1e, S2 and S3 agree with their deposited score logs; all numerical S4 cells agree with the JSONs, subject to findings 3 and 7. Table 3A’s BWTandem/TRF/ULTRA/tantan cells agree with `maize_extra_evidence.json`; Table 3B/C cells agree with `table3bc_replacement.md`; 2026 unfiltered cells agree with their logs. Table 2 and the missing Table 3A cells agree only after using non-deposited cluster inputs.
   - **Minimal fix:** Deposit all active scorers, converted comparator BEDs or durable public archives, active score logs, and the raw timing records. Change “fully covered” to describe the actual partial deposit until that is done.

5. **Severity: MEDIUM — Every Table 1 operating point is score-informed, while only four of 44 tried configurations are shown**

   - **Location:** `manuscript.md:109,150–157,451–475`; `results/manifest.tsv:89–98`; `/data/gpfs/assoc/pgl/devel/exp1_human/loop/ledger.tsv:2–45`.
   - **Defect:** The external optimization ledger contains 44 scored configurations on adotto chromosomes 21/22, with explicit `accept`/`reject` decisions after recall and precision were observed. The reported shared gate parameters were selected through that process; B corresponds to the tested 0.76/max-period-50 arm, H to the 0.72 arm, and F combines the 0.72 arm with the subsequently tested three-copy filter.
   - This is not wholly concealed: `manuscript.md:94`, outside the requested Results scope, acknowledges that BWTandem was selected on the same catalog and competitors were not equivalently searched. Nevertheless, Table 1d is a selected four-point presentation, not the complete prospective sweep, and the ledger establishing that fact is not deposited.
   - The surviving launch wrapper confirms that P/B/F/H used the S2 configurations; I found **no configuration-label mismatch** and no evidence that the human scorer was modified to favor a particular arm.
   - **Minimal fix:** Deposit the complete ledger, mark Table 1d “post-selection/in-sample,” state the selection objective, and distinguish chosen operating points from a prospectively defined sensitivity curve.

6. **Severity: MEDIUM — The 400-call audit’s verdicts match, but its visual evidence and reader provenance are missing**

   - **Location:** `manuscript.md:111`; `results/audit11/README.md:3–22`; `scripts/scoring/sample_specificity_audit.py:4–24`; `results/audit11/aggregate_reviewer2_20260831.txt:6–15`.
   - **Defect:** The deposited aggregate does match the manuscript: 4 SUPPORTED, 346 UNSUPPORTED, 50 UNSURE; 4/400 = 1.0%, 4/350 = 1.14% → 1.1%, Wilson CI 0.4–2.9%; zero supported above period 20; 38/50 = 76% of UNSURE in the shortest stratum.
   - Protocol, sheet, answer key and raw verdict file exist. However:

     - The 400 dot plots actually judged are explicitly “not deposited.”
     - No dot-plot/unit-shift renderer is present in the repository; the sampler generates only the sheet and key. The README’s claim that plots regenerate from the sampler is therefore false as a repository-level instruction.
     - The 809,886-call audit population BED, its full hash and its generation log are absent.
     - The reader is identified only as opaque `reviewer2`; no identity, signature or dated attestation is deposited.
     - The planned second reader was not completed, although this deviation is disclosed.
     - The raw file contains two lowercase `supported` values plus `SUUPORTED`; the README documents only the typo, not the case normalization.

   - **Minimal fix:** Deposit the renderer, exact plotting command/version, population BED or immutable hash/archive, plot archive, and reader attestation. Complete the second-reader protocol or formally amend it.

7. **Severity: MEDIUM — S4 falsely states that BWTandem has the lowest one-to-one precision**

   - **Location:** `manuscript.md:518,525–527`; `results/one_to_one/one_to_one_bwtandem_r50.json:8–10`; `results/one_to_one/one_to_one_tantan_r50.json:8–10`.
   - **Defect:** The prose says “**BWTandem’s is the lowest of the five tools**.” Its own table and deposited JSON give BWTandem 8.68% and tantan 3.42%. Tantan, not BWTandem, is lowest. This error is adverse to BWTandem rather than favorable, but it is still a direct evidence contradiction.
   - **Minimal fix:** Replace the statement with “tantan is lowest; BWTandem is second-lowest.”

8. **Severity: LOW — Environment and scorer-hash provenance contains direct mismatches**

   - **Location:** `manuscript.md:465`; `results/regen/bwt_human_{P,B,F,H}_p100_0363.provenance.json:14–15`; `results/manifest.tsv:2,13,25,46,83,89–93`.
   - **Defect:** S2 says the regenerated whole-genome runs used Python 3.11.14/numpy 2.3.1. The four whole-genome P/B/F/H provenance records say Python **3.11.15** and numpy **2.4.6**.
   - Manifest rows 2/13/46/83 record scorer prefix `b85820c724137549`, but the committed and current `score_table1_regen.py` hash is `b85820bebdd1f139`. Table 2 row 25 similarly records `585d55ba96cb380a`, versus actual `585d55558609d6e5`. Obsolete PENDING rows 8 and 16–19 also remain beside completed rows 89–93 without a superseded label.
   - The P/B/F/H JSON commands omit the controlling environment variables, so the attributed settings are verifiable only from a non-deposited launch wrapper.
   - **Minimal fix:** Correct the versions and hashes, label/remove stale rows, deposit the exact P/B/F/H wrapper, and record environment variables in machine-readable provenance.

9. **Severity: LOW — Several reported band-loss differences disagree by 0.01 point with their deposited source**

   - **Location:** `manuscript.md:254,298–309`; `results/regen/table3bc_replacement.md:23–25,37,62`.
   - **Defect:** The manuscript subtracts displayed two-decimal cells, while the artifact reports deltas calculated before rounding:

     - tantan knob180: 71.92 − 71.42 = **0.50**, artifact 0.51;
     - TRASH knob180: 81.23 − 80.88 = **0.35**, artifact 0.34;
     - TRASH TR-1: 53.17 − 34.61 = **18.56**, artifact 18.57;
     - tantan CentC: 56.51 − 55.59 = **0.92**, artifact 0.91.

   - All other checked ratios, percentages, percentage-point differences, coverage fractions and core-hour comparisons recompute at the stated precision.
   - **Minimal fix:** Use one rounding convention throughout—preferably calculate from unrounded values and expose enough digits to reproduce the reported delta.

**Disposition of the five requested checks:** numerical tracing produced findings 3, 4 and 8; derived arithmetic produced finding 9 but no material arithmetic error; operating-point attribution itself passed, while score-informed selection produced findings 2 and 5; selective reporting produced finding 1 and the TRF omission in finding 3; the audit summary passed numerically but its reproducibility/provenance failed as finding 6.
