# Venue and positioning analysis (2026-09-01, literature-scout, verified)

Question: where can a capability-framed (non-accuracy-leading) BWTandem
paper land, and does the positioning survive the 2026 literature?
Verification level is stated per claim; the scout self-corrected three of
its own initial statements during verification, which are recorded here.

## The finding that outranks the venue table: two 2026 papers sit on our claimed territory

**longdust** — Li H, Li B, "Finding low-complexity DNA sequences with
longdust", *Bioinformatics* 42(3):btag112 (2026), DOI
10.1093/bioinformatics/btag112. Bibliography verified via Europe PMC;
numbers verified against the arXiv:2509.07357 full text (the preprint of
the same paper, not the OUP typeset version):
- CHM13 in **0.94 CPU-h / 0.47 GB** (Table 1) — versus our human
  12.1 h / 21.86 GB. The cost axis is no longer our strength.
- Verbatim: *"TRF would not finish in days with the default option -l.
  With -l30, TRF could run in 11 hours using 17.82GB"* and *"TRF did not
  finish in a week with option '-l20'"* (gorilla) — our TRF-cost
  observation is independently reproduced in print; citing it strengthens
  us.
- Honest-framing precedent, verbatim: *"We do not intend to brand
  longdust as the best tool for LC finding."*
- **Our defense line, in the authors' own words**: *"In comparison to
  tandem repeat finding tools such as TRF and ULTRA, longdust is unable
  to report the repeat units in case of tandem repeats."*

**AniAnn's** — Sweeten A, Schatz MC, Phillippy AM, "AniAnn's:
alignment-free annotation of tandem repeat arrays using fast average
nucleotide identity estimates", *Bioinformatics* btag581 (2026), DOI
10.1093/bioinformatics/btag581. Cite WITHOUT volume/issue (advance
access; a previously reported "42(8)" did not verify). Verified against
the PMC full text + the marbl/anianns README:
- The "+25% F1 over TRASH2" is real but is **βSat class-classification**
  F1, not detection/boundary F1 (HSat2/3: "up to 5%").
- Its Table 1 is **HG002 chr21 diploid**, not a whole genome — its
  2.11 GB is not comparable to our whole-genome numbers.
- Output (README-verified): array intervals + periodicity score +
  monomer length + class label. **Not reported** (absence inference from
  full text + README + CLI options, no explicit limitation sentence):
  consensus motif sequence, copy count, purity/mismatch rate.
- Explicit limitations usable verbatim: classification *"must be done
  through the use of a satellite k-mer database"* (not fully de novo);
  window sizes 1000/2000/4000 bp with *"This case can occur when
  detecting very small arrays (i.e., substantially smaller than the
  window size)"* (short periods out of principle).
- Unfavorable symmetry to disclose: AniAnn's criticizes TRF as *"limited
  to finding repeats less than 2,000 base pairs"* — our benchmarked
  reporting cap is also 2,000 bp (the search window reaches 100 kb but
  nothing above 2 kb is benchmarked). AniAnn's own targets are ≥47 kb
  arrays, so the territorial overlap is smaller than it first looks.

**Positioning consequence** (supersedes "low cost, wide range" as the
lead): BWTandem's defensible unique claim is **single-pass structured TR
records — motif sequence, period, copy count, purity — across 1 bp–2 kb**;
longdust yields intervals without repeat units (their words), AniAnn's
yields boundaries/periodicity/class without per-copy structure (output
spec), and the measured no-competitor band narrows from ">500 bp" to
**500–2,000 bp** — and a same-day follow-up (job 6146343) found the adotto
catalog holds **zero truth annotations above primitive period 500 bp**, so
that band is not adjudicable against adotto at all; the >500 bp evidence
base is the satellite experiments (Col-CEN, maize, CHM13 HOR arrays), and
above 2 kb is AniAnn's/longdust territory with no accuracy evidence of
ours. The release-build range cost is sublinear rather than near-flat:
the 20× widening costs 1.77–1.82× runtime (mean 1.79). It stays as a
supporting scalability point — the field review (Liu & Li 2026,
*Brief. Bioinform.* 27(1):bbag031) names scalability, not long-period
coverage, among open problems.

## Venue ranking (capability framing)

| # | Venue | Fit | Notes |
|---|---|---|---|
| 1 | **GigaScience** (full Research) | Best | Review criteria are "reproducibility, usability and utility", impact explicitly excluded — our deposited BEDs/manifests/audit are the reviewed asset, not garnish. No practical length cap. Precedent: AlcoR (2023, giad101), no accuracy-lead claim. APC ~$2.6k (3rd-party figure). |
| 2 | **Bioinformatics** (Original Paper) | Good, conditional | longdust proves honest framing passes (2026). **Precondition: cite and benchmark against longdust + AniAnn's** — the Li and Phillippy labs are at the center of the reviewer pool. ≤7 pages/~5,000 words; ≥2-year free availability required. |
| 3 | Bioinformatics Advances | Fallback only | Wheeler (ULTRA author) is an Associate Editor there — structural reviewer-network convergence; otherwise fine (ULTRA itself was capability-framed). |
| 4 | BMC Bioinformatics | Good | pytrf (2025, 26:151) is a near-isomorphic precedent (runtime-lead only). |
| 5 | PeerJ / GigaByte | Viable | Soundness-only review / cheap short-format. |
| — | NAR Genom. Bioinform. | Structurally excluded | Methods track requires "statistically significant levels of improvement" over SOTA — contradicts our premise. |
| — | Genome Research | Excluded | Requires an accompanying biological finding (SRF had one; we do not). |
| — | JOSS | Excluded | Does not review benchmarks/claims — our core asset converts to nothing. |

**Format: full paper at either top venue, not an Application Note.** The
contribution's center of mass is the provable benchmark; ~2,600 App Note
words cannot hold the longdust/AniAnn's/ULTRA/TRASH/TRF/tantan
comparison, and cutting it leaves only "no advantage over existing
tools". (The 2026-09-01 Application Note conversion plan is retained as a
fallback design only.)

## Review-risk notes

- Wheeler-lab exposure is on balance favorable (they framed ULTRA the
  same honest way) **provided our ULTRA invocations are beyond
  reproach**: command lines are already in Supplementary Table S1, and
  the running `-p 2000` measurement (job 6145581) closes the
  "unfavorable settings" attack surface. Do not file reviewer-exclusion
  requests — reads defensive.
- Submitting to *Bioinformatics* without engaging longdust/AniAnn's is
  the highest-probability desk-level failure mode; both are 2026 papers
  the editors just handled.

## Action items this creates (author-gated)

1. Add longdust + AniAnn's to Related Work with the verbatim quotes
   above; decide whether to add them to the benchmark (longdust emits
   intervals — scoreable on region metrics; AniAnn's targets satellites —
   scoreable on the Col-CEN/maize experiments).
2. Reframe abstract lead: structured single-pass TR records first,
   range-cost second; narrow the unique band to 500–2,000 bp measured.
3. Venue decision: GigaScience (default) vs Bioinformatics (if willing
   to run the two extra comparisons).
