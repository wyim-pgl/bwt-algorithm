# Academic Integrity Verification Report

**Target:** `/data/gpfs/assoc/pgl/devel/bwt-algorithm/manuscript.md` (468 lines, 17 references)
**Executed:** 2026-08-17
**Verification tools available:** WebSearch unavailable. External verification performed via
Crossref REST API, the DOI resolver (`doi.org`), DBLP, Europe PMC (search + fullTextXML),
DataCite, ENA/BioStudies, PMC, and the GitHub API — all genuine external sources. PubMed
E-utilities, Zenodo's API and bioRxiv were unreachable/rate-limited from this host and were
routed around (see Audit Trail).

## Verification Mode

**Final Verification (Mode 2)** — 100% of registered references, 100% of citation contexts,
100% of statistical/data surfaces, 100% of registered claims.

## Verdict

# FAIL

Two SERIOUS and one MEDIUM issue. Per the verdict criteria, any SERIOUS or MEDIUM issue blocks.

**Important framing:** the reference list itself is clean — all 17 entries verified field-by-field
against a single resolved record each, with zero fabrications. Both SERIOUS issues are
**citation-context defects**: real, correctly-cited papers are made to say things they do not say.

## Verification Summary

| Category | Total | Passed | Issues |
|----------|-------|--------|--------|
| Reference Existence (A1) | 17 | 17 | 0 |
| Bibliographic Accuracy (A2) | 17 | 17 | 0 |
| Ghost Citations (A3) | -- | -- | 0 orphan / 0 dangling |
| Citation Context Accuracy (B1) | 17 refs / 29 in-text uses (100%) | 15 | **2 SERIOUS**, 1 MINOR |
| Statistical Data Accuracy (C1) | 41 cross-checks | 40 | **1 MEDIUM** |
| Internal Consistency (C2) | 38 recomputations | Pass | 0 inconsistencies |
| Claim Verification (E) | 41 of 41 registered (ALL_REGISTERED) | 38 | 1 MAJOR_DISTORTION, 2 MINOR_DISTORTION |

Semantic extraction coverage: `not_machine_detectable` (E1.1 detector not present in this
repository; claim registry built by model-mediated extraction over the full manuscript).

---

## Phase A: Reference Verification — 17/17 VERIFIED, 0 NOT_FOUND, 0 MISMATCH

Every reference reached **VERIFIED**. No entry was left in a gray zone. Each was checked
field-by-field (authors + count, year, title, venue, volume, issue, pages/article number, DOI)
against **one** resolved record, so mashups (PH) and single-field distortions (SH) would surface.

| # | Reference | DOI / locator | Resolver | Field check |
|---|-----------|---------------|----------|-------------|
| 1 | Benson 1999, NAR 27(2):573–580 | 10.1093/nar/27.2.573 | 302 → OUP | all fields exact |
| 2 | Chen et al. 2023, Nat Genet 55:1221–1231 | 10.1038/s41588-023-01419-6 | 302 → Nature | exact; issue (7) omitted, house style |
| 3 | English et al., Nat Biotech 43(3):431–442 (2025) | 10.1038/s41587-024-02225-z | 302 → Nature | exact; **online 2024-04-26 confirmed**, matching the manuscript's own explanatory note |
| 4 | Ferragina & Manzini 2000, FOCS:390–398 | 10.1109/sfcs.2000.892127 | 302 → IEEE | exact; Crossref `issued` is null, **year 2000 confirmed via DBLP**; event name "41st Annual Symposium on Foundations of Computer Science" confirmed |
| 5 | Frith 2011, NAR 39(4):e23 | 10.1093/nar/gkq1212 | 302 → OUP | all fields exact |
| 6 | Gao et al. 2019, Bioinformatics 35(14):i200–i207 | 10.1093/bioinformatics/btz376 | 302 → article 35/14/i200 | all 4 authors exact |
| 7 | **Grossi, Gupta & Vitter 2003, SODA:841–850** | **no DOI** | n/a | **VERIFIED via DBLP**: SODA 2003, pages 841–850, authors Roberto Grossi; Ankur Gupta; Jeffrey Scott Vitter. DBLP carries only an ACM-DL id, no DOI — **the manuscript correctly asserts none**. Ordinal confirmed: DBLP proceedings record reads "Proceedings of the **Fourteenth** Annual ACM-SIAM Symposium on Discrete Algorithms … 2003" |
| 8 | Harris et al. 2019, Bioinformatics 35(22):4809–4811 | 10.1093/bioinformatics/btz484 | 302 → article 35/22/4809 | all fields exact |
| 9 | Kasai et al. 2001, CPM:181–192 | 10.1007/3-540-48194-x_17 | 302 → Springer | 5 authors exact and in order; pages exact. Ordinal confirmed: DBLP proceedings record reads "Combinatorial Pattern Matching, **12th** Annual Symposium, CPM 2001". Crossref container is "Lecture Notes in Computer Science" (LNCS 2089); citing the conference is legitimate proceedings style |
| 10 | Kolpakov, Bana & Kucherov 2003, NAR 31(13):3672–3678 | 10.1093/nar/gkg617 | 302 → OUP | **Crossref lists only 1 author (incomplete legacy OUP metadata) — not a defect in the manuscript.** DBLP confirms all three: Roman M. Kolpakov; Ghizlane Bana; Gregory Kucherov. Volume/issue/pages/year exact |
| 11 | Mirkin 2007, Nature 447:932–940 | 10.1038/nature05977 | 302 → Nature | exact; issue (7147) omitted, house style |
| 12 | **Mori 2008, libdivsufsort, GitHub** | **no DOI** | n/a | **VERIFIED via GitHub API**: repo `y-256/libdivsufsort` exists, owner `y-256` = **"Yuta Mori"**, description **"A lightweight suffix-sorting library"** — verbatim match to the manuscript's title. URL resolves. **Year 2008 confirmed** from the repo CHANGELOG: "## 2.0.0 - 2008-08-23". **No DOI is asserted — correct, none exists** |
| 13 | Naish et al. 2021, Science 374(6569):eabi7489 | 10.1126/science.abi7489 | 302 → Science | exact incl. article number; first 3 authors match "et al." |
| 14 | Nurk et al. 2022, Science 376:44–53 | 10.1126/science.abj6987 | 302 → Science | exact; issue (6588) omitted, house style |
| 15 | Olson & Wheeler 2024, Bioinf. Adv. 4(1):vbae149 | 10.1093/bioadv/vbae149 | 302 → OUP | all fields exact |
| 16 | Willard 1998, COGD 8:219–225 | 10.1016/s0959-437x(98)80144-5 | 302 → Elsevier | exact; issue (2) omitted, house style |
| 17 | Wlodzimierz et al. 2023, Bioinformatics 39(5):btad308 | 10.1093/bioinformatics/btad308 | 302 → OUP | all fields exact |

### Five-Type Hallucination Taxonomy — reference population

| Type | Code | Hits | Notes |
|------|------|-----:|-------|
| Total Fabrication | TF | **0** | every title+author combination resolved to a real record |
| Plausible Author/Conference | PAC | **0** | no author attributed to a paper they did not write; the mreps co-authors, absent from Crossref, were confirmed present via DBLP |
| Incomplete Hallucination | IH | **0 confirmed** (2 flagged, both cleared) | refs 7 and 12 carry no DOI; both were deep-checked externally and both legitimately have no DOI. Neither falsely asserts one |
| Partial Hallucination (mashup) | PH | **0** | every field of every entry was matched against ONE resolved record; no cross-source blending found |
| Subtle Hallucination | SH | **0 at reference level** | no wrong year, expanded initial, or swapped venue. **See Phase B — 2 SH-pattern defects occur at the citation-context level instead** |

**DOI Misdirection check:** all 15 asserted DOIs were resolved through `doi.org`. Every one returned
302 to a publisher landing page for the correct article (e.g. `btz376` → `bioinformatics/article/35/14/i200`,
`btz484` → `35/22/4809`, independently confirming volume/issue/page). Zero misdirections.

### A3. Ghost Citations — clean

All 17 reference-list entries are cited in the body. All 17 distinct in-text citation names resolve
to a reference-list entry. **0 orphans, 0 dangling.**

---

## Phase B: Citation Context Verification — 2 SERIOUS

100% of in-text citation contexts checked. Fifteen are accurate; three are noted below.

### ✅ Contexts verified against retrieved source text

- **ULTRA "HMM-based Viterbi decoding … report a period and consensus"** (§1) — **CONFIRMED**
  from the ULTRA full text (Europe PMC PMC11580682): §2.1 "ULTRA's hidden Markov model",
  §2.2 "Viterbi annotation", "a collection of repetitive states that each model different
  repetitive periodicities."
- **tantan "designed for masking"** and **"defaults to a maximum period of 100 bp"** (§1, §2.2) —
  **CONFIRMED** from the tantan full text (PMC3045581): "We present a new repeat-masking method,
  tantan"; "We set w to 50 for proteins and **100 for DNA**", with the tantan README defining
  `-w` as "maximum tandem repeat period to consider".
- **TRF "imposes a maximum detectable period of 2,000 bp"** (§1) — **CONFIRMED** from the official
  Benson-Genomics-Lab TRF README: "Repeats with pattern size in the range from **1 to 2000 bases**
  are detected."
- **NCRF "requires motifs to be specified a priori"** (§1) — **CONFIRMED** from the NCRF abstract:
  "uncover putative tandem repeats of **specified motifs**".
- **TideHunter "designed for noisy long reads rather than assembled genomes"** (§1, §4.5) —
  **CONFIRMED** from the TideHunter abstract (RCA long-read data).
- **libdivsufsort "O(n log n) in the worst case"** (§2.1.1) — **CONFIRMED** from the repo README:
  "The algorithm runs in O(n log n) worst-case time".
- **adotto catalog, "release v1.2.1", Zenodo 13987414** (§2.2.3) — **CONFIRMED via DataCite**:
  DOI 10.5281/zenodo.13987414 = "Project Adotto Tandem-Repeat Regions and Annotations",
  **version v1.2.1**, creator ADAM ENGLISH, issued 2024-10-24. The GIAB association is confirmed by
  the English et al. abstract ("Genome in a Bottle (GIAB) HG002") and by the catalog's home
  repository `ACEnglish/adotto` — "Working space for the GIAB TR benchmarking project".
- **Ferragina & Manzini for the FM-index; Kasai et al. for linear-time LCP; Grossi et al. for the
  wavelet tree rank index** — each cited work is the correct primary source for the structure named.
- **Benson 1999 "k-tuple matching followed by wraparound dynamic programming"** (§1, §4.3) — the
  canonical TRF algorithm; consistent with the Benson abstract. Full text is not in the Europe PMC
  OA subset, so this rests on the abstract plus the tool's own documentation. No distortion detected.
- **mreps "Hamming distance", "does not model indels"** (§1) — consistent with the mreps abstract
  ("a resolution parameter that allows the program to identify 'fuzzy' repeats"). No distortion detected.

### ❌ IL-SERIOUS-1 — Quantitative claim not supported by the cited source (SH pattern)

**Location:** §1, line 33.
**Text:** "satellites tolerate substantial inter-copy divergence, **20–30% within a CEN180 array**
(Naish et al., 2021), and span megabases."

**What Naish et al. 2021 actually reports** (verified against the full text, PMC10164409):

> "Each satellite was compared with the PPM to calculate a 'variant distance' by summation of
> disagreeing nucleotide probabilities. Substantial sequence variation was observed between
> satellites and the PPM, **with a mean variant distance of 20.2** (Fig. 2A)."

and, for divergence between higher-order repeats:

> "the surrounding CEN180 show higher divergence and copy number variation between the
> higher-order repeats (**94.3 to 97.3% identity**)"

**The defect.** "Variant distance" is a *unitless sum of disagreeing nucleotide probabilities*
across a ~178 bp monomer — **not a percentage**. A mean variant distance of 20.2 over 178 bp is
≈11.3% divergence from consensus, and the paper's only explicit percentage figure for CEN180
inter-copy divergence is 2.7–5.7% (100 − 94.3/97.3). **No 20–30% figure appears anywhere in the
source**; a regex for any "2x to 3x" range returned zero hits in the full text.

This is the classic SH signature: the real number **20.2** survives, but its unit is silently
converted to a percentage and widened into a range. The manuscript's own ground truth is in tension
with the claim — it applies an 80% identity floor, so monomers at 20–30% divergence would be largely
excluded from its own truth set.

**Verdict:** MAJOR_DISTORTION (SERIOUS).
**Caveat recorded:** Naish's supplementary figures (fig. S9, table S4) were not retrievable. The
main-text divergence statistics, which are the paper's headline characterisation of CEN180
variation, contradict the stated range.

### ❌ IL-SERIOUS-2 — A claim attributed to a source that never makes it (SH pattern)

**Locations:** §1 line 33, and §4.3 line 301.
**Text (§1):** "TRF … can become **prohibitively slow on satellite-rich sequence**
(Wlodzimierz et al., 2023)."
**Text (§4.3):** "…with dense CEN180 arrays the most plausible driver …, **a pattern already noted
for TRF on satellite rich sequence generally** (Wlodzimierz et al., 2023)."

**What the TRASH paper actually says about TRF.** The complete TRASH full text was retrieved
(Europe PMC PMC10199239 — Abstract → Introduction → Methods → Results → Summary → References) and
**all 8 occurrences of "TRF" were read in context**. Every limitation TRASH attributes to TRF is
about *annotation granularity*, never runtime:

1. "Although these tools are effective for de novo identification of regions of tandem repeats,
   **they do not precisely annotate individual repeat locations, or HORs**."
2. "TRF identified 78% of CEN178 repeats in a monomer state, and the remaining repeats as dimers."
3. "in both RepeatModeller and TRF outputs, CEN178 repeats are not individually mapped … which
   makes it challenging to extract individual repeats for downstream analysis."

Word-frequency check across the full text: **"slow" 0, "prohibit" 0, "efficien*" 0,
"time-consum*" 0, "cost" 0, "intensive" 0.** The only runtime figures in the TRASH paper are
TRASH's own (under 3 h on Col-CEN; 11 h 24 min on CHM13 across 88 cores). TRF's own README, by
contrast, advertises "The program is very fast."

**The defect.** §4.3's "a pattern **already noted** … (Wlodzimierz et al., 2023)" explicitly claims
prior-literature support for a finding that the cited work does not report. The underlying
observation is sound and is the manuscript's *own* measurement (131.9 h on Arabidopsis against
5.5 h on maize) — the defect is the false attribution of independent corroboration, which makes an
un-replicated single-run observation look like an established result.

**Verdict:** UNVERIFIABLE against the cited source (SERIOUS). Both locations require correction.

### ⚠️ IL-MINOR-1 — Loose characterisation of TRASH

**Location:** §1, line 35. "satellite-specialized annotators like TRASH (Wlodzimierz et al., 2023)
**achieve high accuracy for known repeat families**" groups TRASH with template-guided tools,
whereas the TRASH abstract stresses annotation "**without prior knowledge of repeat composition**".
Mitigating: the manuscript corrects this within the same section (§1 line 39: "TRASH in de novo mode
needs neither a period bound nor a motif") and again at §4.1. Recommend rewording for consistency.

### B2. Citation Format Consistency — MINOR, informational only

Issue numbers are supplied for OUP-family journals (27(2), 39(4), 35(14), 35(22), 31(13), 4(1), 39(5))
but omitted for the Nature/Science/Elsevier entries (Chen 55, Mirkin 447, Nurk 376, Willard 8), whose
issue numbers do exist in Crossref (7, 7147, 6588, 2 respectively). Consistent within each family;
harmonise if the target journal requires it. **Not an accuracy defect.**

---

## Phase C: Data Verification

### C1/C2 — 78 numeric surfaces checked; 1 MEDIUM, 1 MINOR

**Table ↔ manifest ↔ prose tracing is unusually strong.** Representative confirmations, all traced to
a `results/manifest.tsv` row:

| Manuscript figure | Manifest evidence | Result |
|---|---|---|
| Table 1a BWTandem: 4,009,261 calls / 12.1 h / 21.86 GB | job **5983793**, lines 4009261, elapsed **12:05:48** (=12.10 h), RSS **21.86**, node cpu-50 | ✅ exact |
| Range-cost: "7 h 23 m against 5 h 17 m … 1.40×" | job 5983772 **07:22:38** vs job 5981960 **05:16:55**; ratio **1.3966** | ✅ exact |
| "no additional memory (35.28 against 35.24 GB)" | RSS 35.28 (5983772) vs 35.24 (5981960) | ✅ exact |
| Replicate ratios "1.30, 1.30 and 1.41" | range-rep jobs 6076848–53 → **1.3020 / 1.3034 / 1.4086** | ✅ exact |
| Replicate call counts "3,991,867 and 4,009,297" | all three p100 rows = 3991867; all three p2000 rows = 4009297 | ✅ exact |
| Per-pair memory deltas "+5.3, +4.9, −1.6 GB" | 43.85−38.53, 43.00−38.07, 38.32−39.96 | ✅ exact |
| TRF 2,000 bp attempt "6.6 days, 4.7× its 33.7-hour run" | job 6076847 elapsed **6-13:57:48** (=6.582 d = 157.96 h); 157.96/33.7 = **4.688** | ✅ exact |
| Table 1d runtimes 4.5 / 5.3 / 5.3 / 5.2 h | 04:28:54 / 05:18:40 / 05:16:55 / 05:11:20 | ✅ exact |
| Table 1d call counts (all four) | jobs 5982255/5982256/5981960/5982257 | ✅ exact |
| Table 2 BWTandem 0.51 h / 1.95 GB / "31 minutes" | job 5983792, **00:30:49**, RSS 1.95 | ✅ exact |
| Table 2 ULTRA 10.12 h / 4.20 GB; mreps 0.12 h / 0.17 GB | jobs 5981977 (**10:07:11**) and 5981987 (**00:06:57**) | ✅ exact |
| Table 3A/3B/3C BWTandem 15.4 h / 22.41 GB | job 5983794, **15:25:17**, RSS 22.41 | ✅ exact |
| Sweep arm call counts (5 arms) | job 5994683, lines 2766054/3299760/3317408/4009297/4216900 | ✅ exact |
| Ablation arms 111,796 / 111,277 / 111,549 | job 5983739, three rows | ✅ exact |
| Supp. Table S2: 875,224 + 3,116,623 = 3,991,847 | job 5985651, lines 3991847 | ✅ exact, sums |
| `figure_curve_data.csv` (8 rows × 4 metrics) | — | ✅ **all 32 values identical to Table 1d** |
| Region counts ≤ raw manifest lines (caption rule) | every row of Tables 1a, 2 | ✅ holds; tantan's stated 149,706 out-of-chromosome calls = 3,469,229 − 3,319,523 exactly |
| Table 1b post-hoc 3,910,922 + Table 1c 98,339 | = **4,009,261** = Table 1a | ✅ exact closure |
| Ground truth: 66,683 monomers; 177 bp consensus; 80–100% identity; 25/17/17 maize arrays; 5 centromeres | `results/ground_truth/` | ✅ all exact |
| §3.3.3 curated bases: CentC "12.6 Mb", knob180 "29.7 Mb" | computed 12.56 Mb / 29.66 Mb | ✅ exact |

**38 further recomputations** (core-hours 24.2/59.6/33.7; memory ratios 15×/80×; ULTRA 83-fold/151-fold;
unique-region ratio 1.7×; adjusted-precision spread 2.7 pts; raw-precision spread 13.7 pts; sweep span
18.9 pts and largest step 9.0 pts; Arabidopsis ratios 259×/129×/20×/1.6×; scaling 30×/23.7×; all nine
band-loss and TRF-lead deltas; all three TAG-rule factors 15.4/13.9/3.8; 8/11 = 72.7%; the S1.1 copy-count
table and its p=6/p=9 binding argument; the S1.3 GC formula at 0.159/0.164/0.153; the n/30,000 occurrence
cap at ~730–1,090 for Arabidopsis and the 3× budget claim; the 24 Mb stride-clamp threshold) — **all
reproduce.** Internal consistency is clean.

### ❌ IL-MEDIUM-1 — Ground-truth range misdescribed (stated twice)

**Locations:** §2.2.3 line 89, and Table 2 caption line 180.
**Text:** "the 66,683 monomers … **span 175–179 bp** at 80–100% identity" / "the 66,683 resulting
monomers **span 175–179 bp**".

**Measured directly from the deposited file** `results/ground_truth/colcen_cen180.bed` (n = 66,683):

| Quantity | Stated | Actual |
|---|---|---|
| Monomer length range | 175–179 bp | **120–206 bp** |
| Identity range | 80–100% | 80.0–100.0% ✅ |
| Within 175–179 bp | (implied 100%) | **94.59%** (63,077) |
| Outside 175–179 bp | — | **5.41%** (3,606: 3,349 below, 257 above) |

Percentiles: p1 = 157, p5 = 174, p50 = 178, p99 = 179. The modal cluster is correctly described; the
word "span" asserts the full range, and the true range is wrong at both ends by a wide margin.

**Impact.** No reported metric changes — the same coordinate set is applied identically to every row
of Table 2, as the manuscript states. The defect is descriptive, but it is load-bearing for the
argument immediately following it ("Monomer recall is therefore measured over **this conserved
subset**"), which a reader would use to judge how conservative the truth set is.

**Verdict:** MEDIUM (data inconsistent with the deposited source). Correct to "**94.6% fall within
175–179 bp (full range 120–206 bp)**" or similar.

### ⚠️ IL-MINOR-2 — Replicate memory upper bound overstated

**Location:** §3.1, line 98. "Replicate memory spanned **38.1–44.0 GB** by the job-accounting peak."
Manifest `range-rep` RSS values are 38.53, 43.85, 38.07, 43.00, 39.96, 38.32 → true span
**38.07–43.85 GB**. The lower bound rounds correctly to 38.1; the upper bound **43.85 rounds to 43.9,
not 44.0**. Overstatement of 0.15 GB. Does not affect the "no consistent range effect" conclusion,
which the per-pair deltas independently support.

### ⚠️ IL-MINOR-3 — Two manifest rows lack the scorer hash their own contract promises

`results/README.md` states every accuracy figure carries "the scorer, **the scorer's hash**", and the
manuscript's Data Availability repeats "scoring script and hash". Two rows name a scorer but carry an
empty `scorer_sha`:

| Row | Scorer | Backs |
|---|---|---|
| `3.3.3-flanks / flank-periodicity` | `analyze_maize_flanks.py` | §3.3.3 and §4.5 flank-periodicity figures (0.430, 0.268, 0.282, 0.367/0.329/0.405, 72.7%) |
| `S2 / unique-calls` | `analyze_unique_regions.py` | all of Supplementary Table S2 |

Both scripts are present in `scripts/scoring/`, so the figures remain traceable to named deposited
code — they are simply not hash-pinned. (The nine rows with *no* scorer are all cost-only or
did-not-complete rows, where none is expected — consistent with the README.)

### Data surfaces that cannot be traced to a manifest value — all disclosed

These are reported for completeness. Each is explicitly declared in the manuscript, so none is a
concealed gap:

1. **Competitor runtime and memory** (Tables 1a, 1d, 2, 3A, 3B, 3C). A manifest row exists for each,
   but `producer_job` is the literal string `published` and `threads`/`elapsed`/`peak_rss_gb_sacct`/
   `node` are empty. The cells are traceable to a *row*, not to a *value*. Disclosed at §2.2.1 and in
   `results/README.md`, with the two ULTRA/mreps Arabidopsis re-run exceptions correctly carrying job
   records (5981977, 5981987) — verified.
2. **Thread-scaling pair** (§2.2.1, §4.4: 15 h 12 m/13.32 GB and 9 h 21 m/21.38 GB). Declared "side
   measurements, not in the manifest" in the Supplementary Methods S2 table. The 1-thread call-count
   leg of "verified at one, two and four threads" is therefore unverifiable here; the 2- and 4-thread
   legs both trace to manifest rows at 4,009,261.
3. **Index-cache timings** (§3.1: 0.6 s / 9 s / 21 min). Declared "an unlogged side measurement".
4. **Prose-only derived statistics** with no dedicated manifest row, each derivable from a deposited
   BED via a named scorer but not separately pinned: the §4.2 corroboration percentages
   (53.5 / 59.7 / 59.2 / ~71.5 / ~71.6%); the §3.1 distance-class decomposition (13.8 / 19.7 / 14.8 /
   20.1 / 19.0 / 15.0 / 51.7%, and 123.7 / 96.3 / ~13 / 40.1 Mb); the §2.2 above-period-100 base shares
   (33.4% and 65.7%); the Table 1b caption pairing statistics (54.2 / 8.5 / 27.8%, 4.08 of 120.0 Mb);
   and the unbanded substring TAG count of 474,410 in §3.3.1, which appears in no table.
5. **"68,840 raw hits"** (§2.2.3) — the raw hit set is stated as not retained, so the number is
   unverifiable by construction. Disclosed in both the manuscript and `results/README.md`.

---

## Phase D: Originality Verification — NOT RUN

**[D1-SKIPPED: WebSearch unavailable]** Paragraph-level originality screening requires phrase search
against publicly indexed literature. The APIs reachable from this host (Crossref, DBLP, Europe PMC,
DataCite, ENA, GitHub) do not provide full-text phrase matching over the open web, so no D1 sampling
was performed and **no originality claim is made in either direction**. D2 (self-plagiarism) likewise
not run. This is a coverage gap, not a clean result.

---

## Phase E: Claim Verification — 41 of 41 registered claims (ALL_REGISTERED)

| Verdict | Count | Proportion |
|---|--:|--:|
| VERIFIED | 38 | 92.7% |
| MINOR_DISTORTION | 2 | 4.9% |
| MAJOR_DISTORTION | 1 | 2.4% |
| UNVERIFIABLE | 0 | 0% |
| UNVERIFIABLE_ACCESS | 0 | 0% |

- **MAJOR_DISTORTION (1):** the Naish 20–30% CEN180 divergence claim (IL-SERIOUS-1).
- **MINOR_DISTORTION (2):** the replicate-memory upper bound (IL-MINOR-2) and the CEN180 length span
  (IL-MEDIUM-1, which is carried as MEDIUM in the issue list on the data-surface axis).
- The §4.3 TRF-runtime attribution (IL-SERIOUS-2) is recorded as a Phase B source-support failure
  rather than a Phase E numeric distortion — the number is the authors' own and reproduces.

### Claim-strength / scope observations (advisory, not counted in the gate)

The manuscript is unusually disciplined about not exceeding its evidence. Verified examples: it
labels the 0.20 pt Arabidopsis coverage spread, the 0.35 pt CentC lead and the 0.02 pt recall gap
as **ties rather than rankings** and confirms "no claim in the Discussion rests on one"; it states
"we have not isolated seeding's contribution by ablation"; it concedes the corroboration analysis
"does not run in BWTandem's favour"; it flags the flank-periodicity selection effect and the absent
matched control; and it declares "Every figure in this paper is a single-run point estimate, and we
report no confidence intervals or significance tests." **No novelty/primacy claim (E5) requiring
search-bounded wording was found.** No `ADV-E4` scope-broadening rows.

### C3 / C4 — not applicable

No Figure Package with a `figure_table_trace[]` block and no passport with
`experiment_intake_declaration` / `experiment_provenance[]` exists in this repository; this is a
direct-invocation integrity check on a manuscript, not an ARS-pipeline artifact bundle.
**[C3: LEGACY — no Figure Package; trace-unavailable note]** for Figure 1, whose underlying values
were instead verified cell-for-cell against `results/figures/figure_curve_data.csv` and Table 1d.
**[C4-NOT-APPLICABLE: no passport in scope]**. **[E6-SKIPPED: no revision evidence]**.

---

## Issue List (Sorted by Severity)

### SERIOUS (Must Fix)

| ID | # | Category | Location | Issue | Correct Information | Source |
|----|---|----------|----------|-------|--------------------|--------|
| IL-SERIOUS-1 | 1 | Citation context / claim | §1, line 33 | "20–30% inter-copy divergence within a CEN180 array" attributed to Naish et al. 2021; the source reports no such figure | Source reports a **mean variant distance of 20.2** (a unitless sum of disagreeing nucleotide probabilities over a ~178 bp monomer, ≈11%) and **94.3–97.3% identity** between higher-order repeats. Either re-derive and attribute the percentage to the authors' own analysis, or restate the source's figures | Europe PMC / PMC10164409 full text |
| IL-SERIOUS-2 | 2 | Citation context | §1 line 33; §4.3 line 301 | TRF being "prohibitively slow on satellite-rich sequence" / "a pattern **already noted**" attributed to Wlodzimierz et al. 2023, which makes no runtime claim about TRF | TRASH attributes to TRF only annotation-granularity limits ("do not precisely annotate individual repeat locations, or HORs"; "78% … in a monomer state, and the remaining … as dimers"). Present the observation as this paper's own measurement, uncited or self-cited | Europe PMC / PMC10199239 full text, all 8 TRF mentions |

### MEDIUM (Must Fix)

| ID | # | Category | Location | Issue | Correct Information | Source |
|----|---|----------|----------|-------|--------------------|--------|
| IL-MEDIUM-1 | 1 | Data surface | §2.2.3 line 89; Table 2 caption line 180 | Deposited CEN180 monomers stated to "span 175–179 bp"; they span 120–206 bp, with 5.41% outside the stated band | "94.6% fall within 175–179 bp; full range 120–206 bp at 80–100% identity" | `results/ground_truth/colcen_cen180.bed`, n = 66,683 |

### MINOR (Recommended Fix)

| ID | # | Category | Location | Issue | Suggestion |
|----|---|----------|----------|-------|-----------|
| IL-MINOR-1 | 1 | Citation context | §1, line 35 | TRASH grouped with tools that "achieve high accuracy for known repeat families"; TRASH is de novo | Reword to match §1 line 39 and §4.1, which already state this correctly |
| IL-MINOR-2 | 2 | Data surface | §3.1, line 98 | Replicate memory "38.1–44.0 GB"; manifest max is 43.85 | Change upper bound to 43.9 GB |
| IL-MINOR-3 | 3 | Provenance | `results/manifest.tsv` | `3.3.3-flanks` and `S2` rows name a scorer but carry no `scorer_sha`, against the stated contract | Add the two hashes, or soften the README/Data-Availability wording |
| IL-MINOR-4 | 4 | Format | References | Issue numbers present for OUP entries, omitted for Nature/Science/Elsevier entries that have them | Harmonise per target-journal style |

---

## Tool Limitation Disclaimer

> **Phase D (originality) was not executed** — WebSearch is unavailable on this host and no reachable
> API offers open-web full-text phrase matching. No plagiarism screening of any kind was performed,
> and none should be inferred. Professional detection (Turnitin / iThenticate) is required before
> submission.
>
> Phase A/B verification used Crossref, the DOI resolver, DBLP, Europe PMC, DataCite, ENA/BioStudies,
> PMC and the GitHub API in place of WebSearch. These are authoritative for bibliographic metadata and
> for open-access full text, but two full texts (Benson 1999, Kolpakov et al. 2003) are outside the
> Europe PMC OA subset, so their citation contexts were assessed from abstracts plus official tool
> documentation rather than from the article body. Naish et al.'s supplementary figures were not
> retrievable; the finding against it rests on the main text.

## Verification Audit Trail

| Ref / claim | Query | Source reached | Determination |
|---|---|---|---|
| All 15 DOI-bearing refs | `api.crossref.org/works/<DOI>?mailto=…` | HTTP 200 ×15 | full field extraction |
| All 15 DOIs | `curl -sI https://doi.org/<DOI>` | 302 ×15 → publisher | no misdirection |
| Ref 7 (SODA 2003) | DBLP `publ/api?q=High-order+entropy-compressed+text+indexes` | SODA 2003, pp. 841–850, 3 authors, no DOI | VERIFIED |
| SODA ordinal | DBLP proceedings query | "Fourteenth Annual ACM-SIAM Symposium … 2003" | "14th" VERIFIED |
| Ref 9 ordinal | DBLP proceedings query | "12th Annual Symposium, CPM 2001" | "12th" VERIFIED |
| Ref 4 year (Crossref null) | DBLP | FOCS, year 2000, pp. 390–398 | VERIFIED |
| Ref 10 co-authors (Crossref shows 1) | DBLP | Kolpakov; Bana; Kucherov | VERIFIED |
| Ref 12 (libdivsufsort) | GitHub API `repos/y-256/libdivsufsort`, `users/y-256`, CHANGELOG.md | "A lightweight suffix-sorting library"; owner "Yuta Mori"; "2.0.0 - 2008-08-23" | VERIFIED, no DOI asserted |
| ULTRA HMM/Viterbi | Europe PMC PMC11580682 fullTextXML | §2.1 HMM, §2.2 Viterbi | context VERIFIED |
| tantan masking + w=100 | Europe PMC PMC3045581 fullTextXML + tantan README | "w … 100 for DNA" | context VERIFIED |
| TRF max period 2000 | Benson-Genomics-Lab/TRF README | "range from 1 to 2000 bases" | context VERIFIED |
| adotto v1.2.1 | DataCite `10.5281/zenodo.13987414` (Zenodo API unreachable) | "Project Adotto…", version v1.2.1 | context VERIFIED |
| Data accessions | ENA portal + BioStudies + GitHub | PRJEB46164 = the Naish study; E-MTAB-10272 = ONT Col-0 centromeres; schatzlab/Col-CEN | all VERIFIED |
| Assembly accessions | ENA browser XML | GCA_022117705.1 = Zm-Mo17-…-T2T; GCA_000001405.15 = GRCh38 | VERIFIED |
| **TRASH ↔ TRF runtime** | Europe PMC PMC10199239 fullTextXML; all 8 "TRF" contexts; word counts for slow/prohibit/efficien/cost | no runtime claim about TRF | **NOT SUPPORTED → IL-SERIOUS-2** |
| **Naish ↔ 20–30% divergence** | PMC10164409 (via Semantic Scholar OA pointer; EPMC XML 404, bioRxiv 429) | "mean variant distance of 20.2"; "94.3 to 97.3% identity"; zero "2x–3x" matches | **NOT SUPPORTED → IL-SERIOUS-1** |
| All table/manifest cells | local recomputation over `results/manifest.tsv`, `figure_curve_data.csv`, `ground_truth/` | 78 surfaces | 1 MEDIUM, 1 MINOR; remainder exact |

**Unreachable from this host** (routed around, not silently skipped): PubMed E-utilities (empty
response), Zenodo API and web (HTTP 000), bioRxiv (HTTP 429 on four attempts), Europe PMC fullTextXML
for Benson/mreps/Naish (404 — outside the OA subset).
