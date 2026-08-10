# Manuscript revision work plan (2026-08-03, rev. 3)

**Revision 2 supersedes the morning version.** It incorporates the WP0 audit
completed later the same day — full details in `resume.md` (repo root) and
`/data/gpfs/assoc/pgl/devel/exp1_human/wp0/WP0_AUDIT_ALL_TABLES.md`.
The audit resolved some v1 items, invalidated others, and found three larger
defect classes the v1 plan did not know about.

**Revision 3** adds WP3.4 (adotto blind-spot mapping), WP6 (evaluation
methodology strengthening: Truvari rescoring, genome-scale synthetic
benchmark, period-range cost curve, direct CHM13 run, ablation table) and
packaging (bioconda/CI/container) under WP5.

**Hard rule: honest science only — no benchmark gaming, no GT overfitting, no
hiding unfavorable results.** Several audit findings are unfavorable to
BWTandem or correct errors that flattered it. They stay in.

**Target journal: Bioinformatics (OUP) Original Paper.** Every direct
competitor (TideHunter, NCRF, TRASH; ULTRA in the sister journal) published
there. Genome Biology is not viable (no new biology). NAR only if WP3 lands a
real biological result. The manuscript header (L20) still names Genome
Biology / Briefings / NAR — fix during the writing phase.

| # | Work package | Status after 2026-08-03 audit |
|---|---|---|
| WP0 | Benchmark integrity | audit done; fixes + rescoring OPEN |
| WP1 | Precision operating-point curve | open; blocked on WP0-B scoring rule |
| WP2 | Memory / boundary / period assignment | reframed by audit |
| WP3 | Biology (unique regions, CHM13 corroboration) | open, analysis-only |
| WP4 | Novelty framing + prior art + journal retarget | open, writing |
| WP5 | Submission mechanics (Abstract, figures, LICENSE/DOI, packaging) | open, writing |
| WP6 | Evaluation methodology strengthening | added 2026-08-03 (rev. 3) |

Compute constraints: proxy (chr21+22) ~29 min; full human ~7 h; maize ~11 h;
Col-CEN ~7 min; <=2 jobs/node. All runs must use build `d52a4ff` or later
(`706fb76` uninitialised-read fix). `sbatch`/`squeue` are not on PATH — use
`/cm/shared/apps/slurm/current/bin/`. Filip's
`/data/gpfs/assoc/pgl/filip/bwtandem_results/` is read-only; the patched
converter lives in `exp1_human/wp0/convert_to_bed.py`.

---

## WP0 — Benchmark integrity (BLOCKER, do first)

The audit found **three defect classes** plus three smaller errors. Nothing
else is worth doing until these are settled; every downstream number depends
on them.

### WP0.1 RESOLVED — ULTRA maize empty BEDs were a converter bug

~~v1 plan (re-run ULTRA on maize) is obsolete.~~ All three ULTRA maize runs
exited 0 and annotated 2.16 / 8.26 / 8.01 % of the genome. The maize runs used
`--bed` (4 columns, no header); `convert_ultra_tsv()` required `len(row) >= 6`
and dropped every row silently. Also `next(reader)` ate the first data row.
**Fixed** in `wp0/convert_to_bed.py`; recovered 1,114,556 / 1,891,671 /
1,895,076 regions (`wp0/beds/ultra_*.bed`). TSV path regression-checked.

Remaining: re-score the maize ULTRA row from the recovered BEDs (part of
WP0-A below), and fix the copied runtime/memory cells — Tables 3B/3C list
0.2 h / 2.23 GB, which is the exp3A run. Real: **3B = 34 h 43 m / 10.15 GB,
3C = 12 h 51 m / 3.54 GB**. ULTRA is the second-slowest tool in 3B, not the
fastest.

### WP0-A — Defect A: the maize tables score BWTandem under a stricter rule

`compute_metrics.py` scores each tool's per-experiment BED as given,
unfiltered — that produced the published TRF / tantan / TRASH rows. The
BWTandem row came from the band-filtered path (`score_exp3.py`). One table,
two rules. Confirmed by reproducing every published cell under each rule.

Under one consistent rule (unfiltered), BWTandem is **not** an outlier:

| class | tool | cov % | offset (bp) | frag |
|---|---|---|---|---|
| knob180 | TRASH | **81.23** | 9845 | 0.05 |
| | **BWTandem** | 80.32 | **274.5** | 0.00 |
| | TRF | 80.01 | 651.9 | 0.00 |
| | ULTRA | 78.71 | 203.4 | 0.00 |
| TR-1 | TRASH | **53.17** | 10561 | 0.06 |
| | **BWTandem** | 50.49 | 924.5 | 0.00 |
| | TRF | 45.77 | 884.9 | 0.00 |
| | ULTRA | 36.36 | 475.9 | 0.00 |
| CentC | **BWTandem** | **59.02** | 211.5 | 0.00 |
| | TRASH | 58.68 | 5319 | 0.05 |
| | TRF | 58.50 | **69.6** | 0.01 |
| | ULTRA | 57.73 | 58.1 | 0.00 |

Actions:
1. Regenerate every Tables 3A–3C cell under **both** rules (scripts already in
   `wp0/`: `score_maize_consistent.py`, `score_maize_all.py`, etc.). Report
   both; the gap between them is itself a result.
2. Rewrite the surviving real finding: band filtering costs BWTandem far more
   than anyone else (knob180 −15.6 pp vs TRF −0.6; TR-1 **−33.0 pp**; CentC
   −17.8 pp). BWTandem covers arrays best but often assigns a period outside
   the nominal monomer band — a **period-assignment** problem, not a boundary
   problem. This replaces Discussion 4.4's "Tier 3 tolerates degenerate
   flanks" passage (L453), which currently explains an artifact.
3. Manuscript cells that must change: Table 3C text (L388 — "41.54 %" is the
   band-filter artifact; consistent-rule figure is 59.02 %, best in class) and
   Table 3B offsets (BWTandem knob180 274.5 / TR-1 924.5, not 4,273 / 8,110).

### WP0-B — Defect B: tools were not run over comparable period ranges ★ largest

Invisible in the BEDs; only in `benchmarking_results/*/logs/*_run.log`:

| tool | Col-CEN | human GRCh38 | maize Mo17 |
|---|---|---|---|
| **BWTandem** | **1-2000** | **1-2000** | **1-2000** |
| ULTRA | 1-**100** | 1-**100** | 6 / 500 / 200 per exp |
| TRF | 2-**200** | 2-**500** | 6 / 500 / 200 per exp |
| mreps | maxperiod 6, plus a second run 150-400 | — | — |
| tantan | `-f4`, no period parameter | | |

Consequences already established:

- **Table 2 (Arabidopsis) is invalidated for the ULTRA row.** CEN180 monomers
  are 178 bp; ULTRA ran with max period 100, then the table filters to period
  150-200 — a window ULTRA's output cannot intersect by construction. "CEN180
  Count = 0", "Monomer Recall 0.62 %", and L176 ("mreps, ULTRA, and tantan all
  failed to meaningfully detect CEN180") are unsupported as written.
- **Table 1 (human) bp metrics are contaminated.** 35.36 % of BWTandem's bp
  come from period >100, a range ULTRA and tantan could never enter (their
  calls: 0 % >100; TRF: 67.63 % of bp >100). Region metrics barely move (2.47 %
  of regions), but bp recall 43.42 vs 38.14 is not a fair comparison, and
  adjusted precision is biased upward because neither corroborator can
  corroborate a call above period 100.

**Agreed plan (was awaiting go-ahead): do NOT re-run everything at 1-2000.**
Extrapolating ULTRA's own scaling (maize: 13 m 48 s at `-p 6`, 12 h 51 m at
`-p 200`, 34 h 43 m at `-p 500`) gives ~8 days for human at p=2000, and TRF
already took 131.9 h on Col-CEN at p=200. Fix per table instead:

| task | compute | effect |
|---|---|---|
| **Table 1 matched-range rescore** (filter every tool to p<=100, rescore) | **0** | removes bp-metric contamination |
| **Adjusted precision recomputed at matched range** | **0** | removes upward bias |
| Report p=101-2000 as a separate BWTandem-vs-TRF stratum | 0 | turns hidden advantage into explicit one |
| **Re-run ULTRA on Col-CEN at `-p 500`** (134 Mb; est. ~2 h) | ~2 h | makes the Table 2 ULTRA row valid |
| **mreps Col-CEN provenance check** — determine which of the two runs (`-maxperiod 6` vs `-minperiod 150 -maxperiod 400`) was scored into Table 2; if the 6 bp one, the mreps row has the same defect | 0 | — |
| **Table 2 consistent rescore** — Arabidopsis has not been rescored the way maize now has | 0 | — |
| Methods: period-range table; replace "All tools were run with default parameters" (L65 — false) | 0 | discloses Defect B |
| Fix Intro L25: TRF "imposes a maximum detectable period of 2,000 bp" while our own TRF runs used 200/500/6 — reframe as the steep cost of wide ranges, citing the measured ULTRA scaling above | 0 | — |
| Table 3 (maize): caps 6/200/500 match each experiment band — disclose only | 0 | — |

**Defensible framing** (not "competitors detected nothing"): competitors must
be told the period range in advance and pay steeply for a wide one, while
BWTandem covered 1-2000 in a single 6.6 h run. That is the honest version of
the advantage.

**Immediate next step: Table 1 matched-range rescore** (zero compute), with
the ULTRA Col-CEN `-p 500` re-run submitted to SLURM in parallel.

### WP0-C — Defect C: the TAG column measures substring presence

`compute_tag_metrics()` counts a call if `'TAG' in motif or 'CTA' in motif`,
on a field holding a 3-mer for BWTandem/ULTRA/tantan and a multi-kb consensus
for TRF/TRASH — it largely measures motif length. Under a fair definition
(period 3, leading 3-mer in {TAG,AGT,GTA,CTA,TAC,ACT}):

| tool | TAG arrays | TAG bp | longest | published count |
|---|---|---|---|---|
| BWTandem | 35,017 | 2,774,893 | 235.43 kb | 39,299 |
| tantan | 11,752 | 1,579,013 | 235.40 kb | **181,205** |
| ULTRA | 7,031 | 1,873,561 | 235.40 kb | **0** |
| TRF | 1,379 | 2,854,569 | 562.79 kb | 5,263 |

Real picture: TRF and BWTandem recover comparable TAG sequence (2.85 vs
2.77 Mb) — TRF in 1,379 long calls vs BWTandem's 35,017 short ones.
Fragmentation, not sensitivity. L262 ("tantan reported far more arrays than
any other tool") is an artifact and must be rewritten.

### WP0-D — Remaining single-cell fixes

1. **TRASH human template mode** — Table 1 prints the same BED twice (de novo
   and template rows byte-identical; only one human TRASH BED exists on disk).
   Run template mode on human, or collapse Table 1 to a single TRASH row with
   a note.
2. **tantan maize per-experiment filtering** — its 3A/3B/3C BEDs are one
   unfiltered file under three names (identical 68,279,539 bytes). Under a
   consistent band-1-6 filter tantan drops from 54.9 Mb to **13.7 Mb**, below
   BWTandem's 41.4 Mb, so L260's "highest total for tantan" does not survive.
   Fix in the WP0-A regeneration and rewrite L260.
3. **Pin every table cell** to (source BED path, git commit, SLURM job id,
   thread count, env vars) in one `results/` manifest; re-derive every cell
   from it. The manuscript's BWTandem row (80.46/50.65/79.68/43.42/31.53)
   matches neither the pre-fix nor post-fix row in
   `docs/2026-06-25-catch-all-benchmark-for-filip.md` (post-fix `d52a4ff`:
   80.54/50.52/79.5/43.42/31.51); memory differs most (12.99 GB vs 36.3 GB —
   thread count is the likely cause; feeds WP2.1).

**Checked and dismissed (do not re-investigate):** TRASH's fixed analysis
windows and `none_identified` / period `N/A` rows looked like they'd inflate
TRASH; measured effect is nil (knob180 81.23 -> 81.14 %, CentC unchanged) —
they sit outside the GT arrays.

**WP0 exit criterion:** every cell in Tables 1, 2, 3A–3C traces to a named
file produced by a named commit, under a scoring rule the Methods describes,
with competitor period ranges disclosed in a Methods table.

---

## WP1 — Precision (the headline-number problem)

### The actual problem (unchanged, but restated after the audit)

Table 1 reports BWTandem at (80.46 % region recall, 50.65 % region precision)
and ULTRA at (81.62 %, 53.65 %) — a dominated point. The v2.1 operating point
(**59.04 % recall / 58.98 % precision**, above tantan's 57.66 % and ULTRA's
53.65 %, at nearly double TRF's recall) exists but is not reported.

New dependency from the audit: **the WP1.1 runs must happen after WP0-B sets
the final scoring rule** (matched p<=100 range + separate p=101-2000 stratum),
or the curve will be drawn on numbers the paper no longer uses.

### WP1.1 Operating-point curve instead of a single point (highest leverage)

Re-measure 4-5 configurations end to end on the fixed build and report
BWTandem as a curve, competitors as points.

| Label | Config | Expected (pre-fix figures) |
|---|---|---|
| BWTandem-P | v2.1 base, catch-all OFF | 59.0 % recall / 59.0 % prec |
| BWTandem-B | + `CATCHALL_SCAN=1 CATCHALL_MIN_IDENTITY=0.76 CATCHALL_MAX_P=50` | 72.9 / 52.7 |
| BWTandem-F | + `CATCHALL_MIN_IDENTITY=0.72 CATCHALL_MIN_COPIES=3` (current paper) | 80.5 / 50.5 |
| BWTandem-H | + `CATCHALL_MIN_IDENTITY=0.72` | 82.2 / 48.4 |
| (optional) | `CATCHALL_MIN_IDENTITY=0.68` | 84.4 / 42.6 |

Cost: 3-4 full-genome human runs on `d52a4ff` (~7 h each, <=2 concurrent).
Col-CEN is 7 min so add it free; maize only if a maize claim depends on it.
Remember to export `TIER1_FMSCAN=1` at submit time — `run_fullgenome.sbatch`
does not set it, and a June session lost a full-genome run to exactly this.

Deliverables:
- **Figure (recall-precision):** region-level and bp-level panels, BWTandem as
  a connected curve, competitors as points.
- **Table 1 restructured:** default row = BWTandem-F, plus BWTandem-P.
- **Text:** state plainly which BWTandem points are and are not
  Pareto-dominated (F and H **are** dominated by ULTRA in region space — say
  so; carry the argument on bp recall, adjusted precision, unique regions,
  runtime, and the p=101-2000 stratum from WP0-B).

### WP1.2 Decompose bp precision (analysis-only, no runs)

31.53 % bp precision has two very different possible causes: (a) extra bp on
the flanks of otherwise-correct regions (boundary calibration), or (b) whole
spurious regions (false discovery). From existing BEDs
(`exp1_human/diag/bwt_merged.bed`, `data/adotto_primary.bed`): partition
BWTandem bp into inside-adotto / within 500 bp / >500 bp away, as a stacked
bar per period class, with ULTRA and tantan as controls. If (a) dominates,
this reframes weaknesses #1 and #2 as one calibration issue and feeds WP2.2.
If (b) dominates, drop the "low precision is benign" framing and lead with
BWTandem-P as default.

### WP1.3 Do NOT re-attempt intrinsic discriminators

Failed and documented: `CATCHALL_MIN_EXCESS` (2026-07-05), phase/coherence
family (2026-07-08, AUC ~0.5), `CATCHALL_MIN_ENTROPY` (~0 gain). The env-lever
frontier is at its ceiling.

---

## WP2 — Memory, boundary offset, period assignment (reframed by the audit)

### WP2.1 Settle the memory number, then characterise it

12.99 GB (manuscript) vs 36.3 GB (benchmark doc) must be resolved — probable
cause is thread count, since BWTandem holds one FM-index per in-flight
chromosome. Measure peak RSS at 1, 2, and 4 threads on human; plot peak RSS vs
longest chromosome length across all three genomes. Deliverable: supplementary
figure plus one Methods sentence — "peak memory scales with threads x longest
chromosome, so it is tunable; at one thread BWTandem needs ~X GB". The
wavelet-tree rewrite stays in Future directions.

### WP2.2 Are the wide maize boundaries overcall or real peripheral monomers?

**Audit note:** under the consistent (unfiltered) rule, BWTandem's maize
offsets are 274.5 / 924.5 / 211.5 bp — no longer an outlier (ULTRA 203-476,
TRF 70-885, TRASH 5319-10561). The dramatic offsets in the published tables
were the band filter. So this item shrinks to: verify with `src/autocorr.py`
whether the residual flank extension beyond curated coordinates is genuine
degenerate peripheral monomers (windowed identity at 180/358/156 bp vs array
interior and random intergenic controls, background ~0.18). If yes, it is a
result (curated coordinates are conservative) and feeds WP3. If no, implement
a boundary-trim step and re-score the 59 arrays.

### WP2.3 REFRAMED — from "why 41.54 % CentC" to the period-assignment problem

The v1 premise is dead: 41.54 % was the band-filter artifact; under one rule
BWTandem's CentC coverage (59.02 %) is best in class. The real limitation the
audit exposed is the flip side: band filtering costs BWTandem −15.6 to
−33.0 pp because it frequently calls the right array at a period **outside**
the nominal monomer band (e.g. a knob180 array called at 90 bp or 360 bp
primitive period).

Actions:
1. Quantify it: for the 59 curated arrays, extract BWTandem's called period(s)
   vs the nominal band; report how often the primitive period is a harmonic or
   sub-harmonic of the family period (autocorr can confirm the true period).
2. If harmonic confusion dominates, evaluate a post-processing
   period-normalisation step (report both the called primitive period and the
   array-level dominant period) — analysis first, code only if the data says so.
3. This becomes a Discussion subsection on period assignment, replacing the
   "degenerate flanks" passage — it is an honest limitation with a measured
   cost, which is exactly what a reviewer wants to see handled openly.

---

## WP3 — Biology (decides Bioinformatics vs NAR)

Everything needed is already on disk; no new genome runs.

### WP3.1 Characterise the 863,654 human unique regions

Inputs: `exp1_human/diag/bwt_merged.bed`, `beds/{ultra,tantan,trf}.bed`,
`data/adotto_primary.bed`. Define the unique set (overlaps no other tool, no
adotto region); profile period distribution, motif composition, GC, entropy,
array length; test enrichment for gene/exon/intron/intergenic, segmental
duplications, subtelomeric/pericentromeric bands against length- and
GC-matched shuffled background; cross-reference pathogenic STR loci.
**Note:** recompute the unique set after the WP0-B matched-range rescore —
the count will change. If the unique calls turn out to be uniformly
low-complexity with no enrichment, that is also a real finding, and it means
BWTandem-P must become the paper's default configuration.

### WP3.2 Independent-assembly corroboration (the validation that matters)

Cross-assembly recurrence does not share the "related heuristics corroborate
each other's errors" failure mode the Discussion already concedes (L453).
Lift the unique regions from GRCh38 to CHM13 (`exp1_human/chm13/` holds the
annotation BEDs) and measure recurrence; compare against adotto-supported
BWTandem calls (positive control) and shuffled intervals (negative control).
Deliverable: "N % of BWTandem-unique regions recur in an independent T2T
assembly, vs M % for shuffled controls" — the strongest single sentence
available for the precision defence.

### WP3.3 Optional: Arabidopsis monomer-level structure

Per-monomer divergence profile along each Col-CEN centromere from the 1,380
CEN180 calls. Novelty is methodological (Naish et al. 2021 already described
the structure). Include only if WP3.1/3.2 come in thin.

### WP3.4 Map adotto's structural blind spots onto the unique calls

The adotto construction is fully documented
(github.com/ACEnglish/adotto, `regions/` + `DataDescription.md`), so the
"catalog gap vs genuine false discovery" question is partly answerable by
enumeration, not just by recurrence:

1. **Region sources are an STR/VNTR-centric union:** UCSC simpleRepeat
   (i.e. TRF), Illumina/ExpansionHunter catalog, ensembleTR, adVNTR; v0.3
   added TRGT, pbsv, Vamos.
2. **A 58 % monomorphic filter:** v0.2 removed every region with no observed
   non-SNP variant across 172 haplotype-resolved long-read assemblies. Real
   but population-invariant arrays are absent from the catalog *by design*.
3. **Motif annotation is TRF at maxperiod 500** — period >500 bp arrays cannot
   be annotated (connects to Defect B's p>100 problem).
4. **Span filters:** <10 bp and >50 kbp regions removed; chr1-22/X/Y only.

Analysis (compute ~0, existing BEDs):
- Quantify adotto's coverage envelope: period distribution of `annos`
  (confirm the 500 cap), span distribution (10 bp-50 kb).
- Split the BWTandem-unique set (from WP3.1) into **outside the envelope**
  (period >500, span >50 kb, centromeric/satellite classes absent from the
  contributing tools — plausible catalog gaps) vs **inside the envelope**
  (contributing tools could have called these but none did — the hard-to-
  defend set).
- Layer WP3.2's CHM13 recurrence on the split: "N % of outside-envelope
  unique calls recur in CHM13" turns the GT-bias argument from assertion
  into measurement. Inside-envelope non-recurring calls are reported
  honestly as likely false discoveries.
- Optional: adotto's pVCF of 86 haplotypes is on Zenodo (variants v0.1) —
  test whether outside-envelope unique regions are variable.
- Writing: replace the generic "catalog is incomplete" sentence (L448) with
  the three concrete mechanisms (58 % monomorphic filter, TRF maxperiod 500,
  50 kb span cap), cited to English et al. 2024.

---

## WP4 — Novelty framing, prior art, journal retarget (writing)

### WP4.0 Retarget the manuscript to Bioinformatics

Fix L20 (target journal line). Check Bioinformatics formatting limits
(Original Paper: ~7 pages main text-equivalent; Abstract <=200 words) during
the writing phase; current main text is roughly in range but figures will add
length.

### WP4.1 Literature pass (verify every citation exists before adding)

Missing and near-certain to be raised:
- **srf** (Heng Li group) — de novo satellite discovery without a catalog;
  occupies exactly the niche Discussion 4.1 claims. Must cite + differentiate.
- **RepeatModeler2 / RepeatMasker** — the de facto de novo repeat-family pipeline.
- **StringDecomposer / HiCAT / alpha-CENTAURI** — monomer/HOR decomposition.
- Suffix-tree/suffix-array tandem repeat algorithms (Stoye & Gusfield line) —
  the algorithmic ancestry mreps already builds on.
- Any prior FM-index/BWT-based tandem repeat work — search explicitly; if it
  exists and is uncited the novelty claim collapses.

### WP4.2 Differentiation table (new Table 1 or Supplementary Table S1)

Axes: input (assembly / reads / k-mer counts) x prior knowledge (none / motif
/ library) x output (mask / consensus motif / monomer decomposition) x period
range x de novo. BWTandem's defensible niche: **assembly input, no prior
catalog, periods 1-2,000 bp in a single pass, emits consensus motifs.** srf is
satellite-focused, RepeatModeler is TE-family-focused,
StringDecomposer/HiCAT need a monomer supplied, tantan emits masks not motifs,
NCRF needs the motif a priori. Show it, don't assert it.

### WP4.3 Promote the algorithmic argument into the Introduction

Discussion 4.3 (L450) is the novelty claim — a k-mer sampled inside an array
recurs at every copy, so occurrence positions form an arithmetic progression
whose common difference *is* the period; the period is readable from the
occurrence list before any alignment, whereas alignment-based methods must
propose a period and test it by DP. It is buried on page 12. Move it to the
Introduction and state the contribution as "reading the period off the
occurrence structure", not "using an FM-index".

### WP4.4 Supporting evidence already in hand

The Arabidopsis-vs-maize runtime inversion (TRF: 131.9 h on 134 Mb vs 5.2 h on
2.18 Gb; BWTandem: 0.31 h vs 6.6 h) — make it a figure (runtime vs repeat
density, all tools), not a sentence. Keep the existing hedge (the effect was
not isolated to k-mer seeding). **Audit caveat:** this figure must use the
corrected competitor numbers — ULTRA's real maize 3B runtime is 34 h 43 m, not
the 0.2 h currently printed, and competitor period caps belong in the caption
or a companion table.

---

## WP5 — Submission mechanics (all currently missing, all blocking)

1. **Abstract** — section exists, body is empty.
2. **Figures** — zero figures exist (`figures/` holds only an mp4/gif).
   Planned set: (F1) recall-precision operating-point curve (WP1.1);
   (F2) runtime vs repeat density (WP4.4); (F3) maize consistent-vs-band
   scoring, both rules (WP0-A); (F4) unique-region characterisation (WP3.1);
   (F5) CHM13 corroboration (WP3.2/WP6.4); (F6) period-assignment / harmonic
   analysis (WP2.3); (F7) period-range cost curve (WP6.3); (F8) synthetic
   precision/recall vs divergence (WP6.2). Supplements: memory scaling
   (WP2.1), bp-precision decomposition (WP1.2).
3. **Availability** — no LICENSE file, no release tag, no Zenodo DOI. Choose
   license, tag the submission commit, mint the DOI.
4. **Methods completeness** — full env-var listing per experiment (the
   `CATCHALL_*` / `SAT_FILL_*` / `TIER1_FMSCAN` settings actually used), the
   competitor period-range table (WP0-B), and the corrected run parameters.
5. **Author block** — placeholder names/degrees/sign-offs still blank.
6. **Packaging** — no bioconda recipe, no CI, no container. Bioinformatics
   de facto expects one-command installation: write a bioconda recipe (or at
   minimum a working `pip install .`), wire GitHub Actions to run the
   existing pytest suite on push, and publish a Docker/Singularity image
   (extend `benchmarking/Dockerfile`). Add a `CITATION.cff`. This is cheap
   relative to its effect on the "software availability" review criterion.

---

## WP6 — Evaluation methodology strengthening (added rev. 3)

Four additions beyond the audit-driven plan, chosen 2026-08-03.

### WP6.1 Rescore with the Truvari protocol

adotto's own documentation tells users to compare callers with **Truvari**.
If we score with home-grown overlap scripts while the benchmark authors ship
an official protocol, the first review question is "why not Truvari?".
Convert BWTandem and competitor outputs to Truvari-compatible input and
reproduce Table 1's headline metrics under the official protocol; report
both scorings if they differ materially. Effort: conversion scripts +
validation against the existing numbers; no new genome runs.

### WP6.2 Genome-scale synthetic benchmark (divergence-controlled)

The only evaluation that fully escapes the adotto circularity: plant tandem
arrays of known period (spanning 1-2,000 bp) and controlled inter-copy
divergence (0-30 %) into background sequence, so ground truth is exact by
construction. Measure precision/recall as a function of divergence for
BWTandem and the fast competitors (tantan, ULTRA, TRF — all cheap on a
synthetic contig set). This extends the existing 44-repeat synthetic test
(README) to genome scale and answers "is the 50.65 % precision real?" with
zero GT ambiguity. Deliverable: F8 (precision/recall vs divergence curves) +
one Results paragraph.

### WP6.3 Period-range cost curve (turns Defect B into a selling point)

Already-owned numbers: ULTRA on maize took 13 m 48 s at `-p 6`, 12 h 51 m at
`-p 200`, 34 h 43 m at `-p 500` — a steeply superlinear cost of widening the
period range — while BWTandem covers 1-2,000 in one flat-cost pass. Add at
most 1-2 ULTRA points (e.g. `-p 100`) and optionally TRF's Col-CEN point
(131.9 h at 200), then plot runtime vs period cap with BWTandem as a flat
line. Compute ~0. Deliverable: F7 — the empirical version of the WP4.3
argument, and the honest frame for Defect B.

### WP6.4 Direct BWTandem run on CHM13 (strengthens WP3.2)

WP3.2 as planned measures recurrence via liftover of GRCh38 calls to CHM13
coordinates. A direct BWTandem run on the CHM13 assembly (~7 h, one job)
removes the coordinate-conversion critique entirely: call recurrence at
orthologous loci measured from two independent de novo annotations. Submit
to SLURM alongside the WP1.1 runs (same build, same config as the human
paper run). Deliverable: the WP3.2 sentence upgraded to "independently
detected in both assemblies", not "lifted and overlapped".

### WP6.5 Ablation table

A methods paper without an ablation gets asked for one. Most numbers already
exist and only need compiling into one table: three-tier pipeline
contributions (tier-by-tier recall on human), `TIER1_FMSCAN` on/off
(57.6 % vs 44.4 % region recall), satellite gap-fill on/off (92.5 -> 99.8 %
CHM13 centromere bp recall), catch-all on/off (the WP1.1 P-vs-F rows).
Fill gaps with proxy runs (chr21+22, ~29 min each) rather than full genomes
where possible. Deliverable: one Table or Supplementary Table + 3-4
sentences in Results.

---

## Sequencing

**Phase 1 (week 1) — unblock.** WP0-B matched-range rescore (compute 0, start
now) + ULTRA Col-CEN `-p 500` to SLURM in parallel. WP0-A maize regeneration
from existing `wp0/` scripts. WP0-D single-cell fixes. WP4.1 literature pass
in parallel (no compute). Start WP1.1 runs the moment the scoring rule is
pinned — they are the long pole; submit the WP6.4 CHM13 run in the same
batch (second slot, <=2 jobs/node).

**Phase 2 (week 2) — analysis on existing data.** WP1.2, WP2.2, WP2.3.1,
WP3.1, WP3.2, WP3.4, WP6.1 (Truvari conversion), WP6.2 (synthetic design +
competitor runs on synthetic contigs), WP6.3 (cost curve from existing
logs) — all analysis-only and mutually independent.

**Phase 3 (week 3) — decide and fix.** WP2.3 period-normalisation if the
harmonic analysis says so; WP2.2 boundary trim only if flanks are background;
WP2.1 memory characterisation; WP6.5 ablation compilation (+ proxy runs);
re-score anything touched; build the `results/` manifest (WP0-D.3) last,
over final numbers.

**Phase 4 (week 4) — write.** Abstract; eight main figures; WP4.0/4.2/4.3
rewrite of Intro and Discussion; Methods disclosure tables; WP5 availability
+ packaging (bioconda recipe, CI, container).

## Do not do

- New intrinsic precision discriminators (WP1.3 — failed twice, documented).
- Wavelet-tree memory rewrite (Future directions only).
- Additional env-var sweeps hoping to beat the frontier (documented ceiling).
- More genomes as *benchmark claims*. Three is enough; a fourth adds runtime,
  not argument. (Exception: WP6.4's CHM13 run — same species as GRCh38, used
  for cross-assembly corroboration, not as a new benchmark row.)
- Re-running every competitor at p=2000 to "match" BWTandem (~8 days for
  ULTRA on human; TRF may not finish). The matched-range-down rescore plus
  the p=101-2000 stratum is the honest fix.

## Journal decision gate (updated)

- WP3.2 gives a strong corroboration number **and** the WP0 rescoring keeps
  BWTandem competitive at matched range -> NAR is a reasonable shot,
  Bioinformatics as fallback.
- WP3.1/3.2 come in thin but WP0/WP1/WP2 are clean -> **Bioinformatics
  Original Paper** (the default target), BWTandem-P reported as default if
  WP1.2 shows false discovery dominates.
- The WP0-B matched-range rescore erases BWTandem's bp-recall edge entirely
  -> reconsider before writing: the honest story then rests on the operating
  point curve (WP1) and single-pass 1-2000 coverage (WP0-B framing); submit
  to Bioinformatics Advances or BMC Bioinformatics rather than delay further.
