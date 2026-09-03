# Finding dispositions — one table for every review finding

Built 2026-09-03 at commit `8525c01`, after the Codex plan review found that the
"31 unverified findings" figure could not be audited: the round-2 and round-3
transcripts were cited by `todo.md` while living only in a session scratchpad.
They are deposited now, and this table is the deduplicated disposition record.

## Verdict vocabulary

Three states, never two. **Absence of an artefact routes only to the third.**

| Verdict | Requires |
|---|---|
| `CONFIRMED` | positive evidence reproducing the defect, plus the check that reproduced it |
| `REJECTED` | **positive contrary evidence** — a failed `grep` is not evidence |
| `BLOCKED-ON-MISSING-ARTIFACT` | the artefact needed to judge does not exist |
| `RESOLVED` | confirmed earlier and fixed; the resolving commit is named |
| `PARTIAL` | fixed on some surfaces and live on others; both named |

The rule exists because on 2026-09-03 three findings were recorded as CONFIRMED
and later withdrawn (`quarantine.md` §6.6, §6.16, §6.27). All three failed the
same way: the evidence *location* was recorded, the *check* was not, and an
absence was promoted into an error. Every row below therefore carries a check a
stranger can re-run, for rejections and confirmations alike.

## Codex round 1 — 20 findings

| ID | Claim | Verdict | Check / evidence | Status |
|---|---|---|---|---|
| R1-1 | Table 2 TRASH template row costs one run for a three-run union | `CONFIRMED` | sort -u of the three component BEDs = 591 lines, set difference 0 both ways; CEN159 log 6:18:29 / 1,389,316 KiB equals the printed cell | OPEN — quarantine §6.1, A-2 |
| R1-2 | manifest.tsv:67 keeps four claims the manuscript corrected | `CONFIRMED` | sed -n '67p' results/manifest.tsv contains 'same binary 1.2.1', 'version- and input-matched', bare '1.55x', centromere attribution | OPEN — §6.2, A-2 |
| R1-3 | ULTRA partial output ends in a truncated record; count and endpoint from different populations | `CONFIRMED` | Codex plan review re-verified: header + 138,425 nine-field records + one six-field fragment, no final newline; last complete endpoint 124,785,432 | OPEN |
| R1-4 | 'Same FASTA' is not literally true (GenBank vs RefSeq identifiers) | `CONFIRMED` | grep -c 'same FASTA' manuscript.md = 1, still present at 4aeba8e; sequence content identical (3,209,286,105 bases) but identifiers and file differ | OPEN |
| R1-5 | 1.55x / 4.7x not marked as lower bounds everywhere | `CONFIRMED` | grep -c '1.55' = 3 in manuscript.md and 3 in results/manifest.tsv; only range_cost_attempts/README.md carries the dagger | OPEN |
| R1-6 | GiB relabel incomplete across reporting surfaces | `PARTIAL` | manuscript.md fully relabelled in 63ccd07 (37 occurrences, 1 GB left and it is a file size); results/ surfaces still say GB | PARTIAL — manuscript done, results/ open |
| R1-7 | Two Arabidopsis memory cells do not follow the stated GiB rule | `PARTIAL` | manuscript Table 2 cells fixed to 0.06/0.48 in 2901ff3; results/comparators2026/README.md:27-28 still 0.07/0.50 | PARTIAL |
| R1-8 | 'Comparable cost cells' overstates peak-memory comparability | `CONFIRMED` | BWTandem memory is sacct cgroup, 2026 competitors are GNU-time process accounting; manuscript concedes the mechanisms differ elsewhere | OPEN |
| R1-9 | AniAnn's manifest thread counts are 1, runs used -j 2 | `CONFIRMED` | run_anianns.sbatch:5,22,23 and every run log say -j 2 under --cpus-per-task=2; manifest rows 123/126/129/132 say 1 | OPEN — §6.3, A-2 |
| R1-10 | 'Every competitor run used one Singularity container' is false | `PARTIAL` | manuscript.md:481 (Table S1 caption) corrected in 2901ff3; manuscript.md:76 still makes the claim | PARTIAL |
| R1-11 | Cost-provenance prose overbroad and asymmetrically favourable | `PARTIAL` | 'every competitor cost cell is checkable' falsified by R1-1 and recorded as quarantine §3.6; results/README.md:52 still says competitor costs cannot be checked, now false after e4ae632 | PARTIAL |
| R1-12 | Human TRASH cost comes from a nonzero-exit run, undisclosed | `RESOLVED` | Deposited log shows 'Exit status: 1' after a Circos error; disclosed in the Table 1a dagger footnote | e4ae632 |
| R1-13 | Evidence trail cluster-checkable but not repository-self-contained | `PARTIAL` | 40 GNU-time logs and sacct for 34 jobs deposited in e4ae632; manifest still references 75 absolute cluster paths for the large BEDs | PARTIAL |
| R1-14 | Manifest keeps live-looking PENDING rows beside completed replacements | `CONFIRMED` | rows 8 and 16-19; sacct_provenance.txt shows the originals failed/OOM at startup, so the label should be superseded-failed-at-startup | OPEN — A-2 |
| R1-15 | 'About 4% processed' is stronger than the output proves | `CONFIRMED` | manuscript.md:105 still reads 'about 4% of the assembly processed'; only the last emitted coordinate is known | OPEN |
| R1-16 | --adj code-path explanation literally overbroad | `CONFIRMED` | score_table1.py:403 reads adj_rules to gate an informational note; the README says 'nowhere earlier'. Numeric conclusion unaffected | OPEN — LOW |
| R1-17 | mreps provenance row says 'not reported' while the manuscript reports the cost | `CONFIRMED` | comparator_baselines.md:30 says 'not reported'; manuscript Table 1a gives 0.9 h / 6.38 GiB; deposited log reproduces 54:41.08 / 6,691,220 KiB | OPEN |
| R1-18 | BWTandem output described as 'per-copy' | `CONFIRMED` | grep finds 'per-copy' in manuscript.md (2) and comparators2026/README.md (1); output is per-array with an aggregate copy count | OPEN — LOW |
| R1-19 | longdust non-detection given an unobserved mechanistic reason | `CONFIRMED` | comparators2026/README.md:104 attributes zero overlap to the tool's model; zero overlap establishes non-detection only | OPEN — LOW |
| R1-20 | Passive/absolute competitor-failure framing survives | `CONFIRMED` | comparator_baselines.md:42 'rerun infeasible' and README.md:176 'caps could not be lifted' convert author cancellations into tool properties | OPEN |

## Codex round 2 — 9 findings

| ID | Claim | Verdict | Check / evidence | Status |
|---|---|---|---|---|
| R2-1 | Preregistered AniAnn's banded cells omitted; deposited log has them and they lead | `RESOLVED` | score_2026_human.txt:25-33 and score_2026_colcen.txt:12 hold the cells; added to Tables 1b, 1c and §3.2 | 05a76a6 (C-9a) |
| R2-2 | Arabidopsis catch-all setting chosen by scoring both arms on the evaluation truth | `RESOLVED` | catchall_experiment_results.md scores both against colcen_cen180.bed: off 65.54% vs on 60.72% bp precision; S2 now states this and the absence of a held-out set | d580840 (C-11a) |
| R2-3 | Scorer discards TRF's period column, hiding a rate that beats BWTandem | `RESOLVED` | convert_to_bed.py:93 writes period to col5 and the BED holds 6/29/76; load() kept col5 only when it meant copies. Re-scored: only the five period fields move, TRF 63.53% | 903245c, 2901ff3 |
| R2-4 | 'Fully covered' accuracy deposit is not self-contained | `PARTIAL` | e4ae632 deposited the three missing scorers, 40 competitor logs and sacct for 34 jobs; most comparator BEDs remain absolute cluster paths | PARTIAL |
| R2-5 | Every Table 1 operating point is score-informed; only 4 of 44 shown | `RESOLVED` | ledger.tsv has 44 scored rows and best.json states the ULTRA-referenced objective; Methods 2.2.3 now says post-selection and in-sample, and the ledger is deposited | d580840, e4ae632 |
| R2-6 | 400-call audit verdicts match but the visual evidence and reader provenance are absent | `CONFIRMED` | results/audit11/ holds README, sheet, answer key, aggregate and verdicts; the 400 dot plots, the renderer, the population BED and a reader attestation are not there | OPEN |
| R2-7 | S4 falsely states BWTandem has the lowest one-to-one precision | `RESOLVED` | one_to_one_tantan_r50.json gives 3.42% against BWTandem's 8.68%; caption now says second-lowest | 2901ff3 |
| R2-8 | Environment and scorer-hash provenance mismatches | `PARTIAL` | The P/B/F/H records are the period-100 panel, not the whole-genome runs, so S2's 3.11.14/2.3.1 was right; the later panel's environment is now disclosed. The two manifest scorer-hash mismatches are open | PARTIAL — 2901ff3 |
| R2-9 | Band-loss deltas disagree with the deposited artifact by 0.01 | `RESOLVED` | table3bc_replacement.md computes each loss before rounding; four values corrected including one I had introduced by re-doing the subtraction | 4aeba8e |

## Codex round 3 — 14 findings

| ID | Claim | Verdict | Check / evidence | Status |
|---|---|---|---|---|
| R3-1 | The no-GT-overfitting rule was violated, partially disclosed | `RESOLVED` | Same evidence as R2-5; Methods now names the campaign, its size and its objective | d580840 |
| R3-2 | 400-call audit evidence not reproducible as claimed | `DUPLICATE` | Same defect as R2-6 | see R2-6 |
| R3-3 | Scoring scripts and ground-truth intervals not all present | `RESOLVED` | rescore_tables_3bc.py, score_exp3.py and score_overlap.py deposited with cluster paths made env-overridable; raw BLAST hit set deposited | e4ae632 |
| R3-4 | Machine-readable run provenance does not record the controlling environment | `CONFIRMED` | regen_*.provenance.json records the command but not the env vars that select the configuration; the launch wrapper is not deposited | OPEN |
| R3-5 | S1.3 reports the wrong Tier-3 jitter tolerance for every run | `CONFIRMED` | manuscript.md:440 asserts a 0.04 ceiling in every run; tier3.py gives 0.02 + 0.02*(max_period/1e5) = 0.0204 at --max-period 2000 | OPEN — tied to §6.17 and C-1 |
| R3-6 | S1.4 gives neither the implemented autocorrelation definition nor the valid-base threshold | `PARTIAL` | The 70%/80% half is fixed in 2901ff3; the autocorrelation definition is still absent | PARTIAL |
| R3-7 | Per-figure provenance manifest does not exist and the figure programs are stubs | `PARTIAL` | The Abstract promise was removed in 2901ff3; all six plot_*.py still end at a TODO and no final figure is tracked | PARTIAL — blocks submission |
| R3-8 | A manuscript-plus-repository reader cannot regenerate the benchmark | `PARTIAL` | Substantially improved by e4ae632; the 75 absolute cluster paths remain | PARTIAL |
| R3-9 | Discussion 4.3 confines FM-index use to Tiers 2 and 3 | `RESOLVED` | Every configuration sets TIER1_FMSCAN=1 (manuscript.md:456) and tier1.py:325 runs the enumeration; 4.3 rewritten | 2901ff3 |
| R3-10 | Discussion 4.2 gives an incomplete genome-by-pass account | `RESOLVED` | Gap-fill runs on all three genomes, catch-all on human and maize; the matrix is now stated with ablated and non-ablated distinguished | 2901ff3 |
| R3-11 | S1.2 misidentifies the 2% primitive-period test as autocorrelation | `CONFIRMED` | motif_utils.py:915-923 compares s[i] against template[i % p] where template = s[:p] — a template-match test. Autocorrelation at lag p compares s[i] against s[i+p]. Different statistic | OPEN |
| R3-12 | Limitations omits the unresolved native-path flake | `RESOLVED` | Added as a fourth limitation with the exhausted diagnostics and the reason it is currently unreproducible | 4aeba8e |
| R3-13 | S1.3's values-in-force table excludes the Col-CEN organelles | `RESOLVED` | ChrC and ChrM take the short-sequence branch (k=12/13, stride 20, max occ 200), outside every printed range; now stated | 2901ff3 |
| R3-14 | 'Up to two million bases' is not the implemented Tier-1 sampling cap | `REJECTED` | **Positive contrary evidence**: tier1.py:361-370 `_base_freqs` uses `step = max(1, n // 2_000_000)`, which samples at most about two million bases. The manuscript's wording matches the implementation | REJECTED 2026-09-03 |

## Tally

| Verdict | Count |
|---|--:|
| `CONFIRMED` | 18 |
| `RESOLVED` | 13 |
| `PARTIAL` | 10 |
| `DUPLICATE` | 1 |
| `REJECTED` | 1 |
| **total** | **43** |

Nothing in the three Codex rounds is "unverified" any more. 13 were fixed by
today's commits and name the commit. 18 are confirmed and open, each naming the
check that reproduced it. 10 are partial — fixed on one surface and live on another,
which is the failure mode behind most of today's rework, so both sides are named
rather than the item being called done.

**One is REJECTED**, and it is the only one of 43: R3-14, where the manuscript's
"up to two million bases" turned out to match `tier1.py:361-370` exactly. Four
adversarial passes produced one false positive. That ratio is not comfortable
reading, and it should not be softened: the three withdrawals in `quarantine.md`
(§6.6, §6.16, §6.27) were **my own confirmations**, not reviewer errors — the
reviewers raised the questions and I confirmed them wrongly.

The Kimi rounds are not tabulated here. Every Kimi finding that was acted on is
in `quarantine.md` §6.7–6.16 with its own evidence; the ones not acted on are the
MEDIUM/LOW residue, which has no deduplicated inventory. That residue is the one
part of this record still incomplete, and it is the next thing to close.
