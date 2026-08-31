# TRASH vs bwtandem Comparison

## Overview

**TRASH** (Tandem Repeat Annotation and Structural Hierarchy; Wlodzimierz et al., 2023, *Bioinformatics*) is an R-based tool for centromere satellite repeat annotation. It uses k-mer frequency analysis to decompose tandem repeat arrays into individual monomers and classify them into repeat families (e.g., CEN178, CEN159, 5S rDNA). It also supports Higher-Order Repeat (HOR) analysis.

**bwtandem** uses BWT/FM-index to detect tandem repeat arrays across all period lengths (1 bp to 100 kbp) without requiring pre-defined repeat family templates.

## Tool Characteristics

| Feature | TRASH | bwtandem |
|---------|-------|----------|
| Algorithm | k-mer frequency + sliding window | BWT/FM-index + 3-tier pipeline |
| Input requirement | FASTA + optional repeat templates | FASTA only |
| Output granularity | Individual monomers | Continuous repeat arrays |
| Family classification | Yes (CEN178, CEN159, 5S rDNA, etc.) | No (reports motif sequence) |
| HOR analysis | Yes | No |
| Period range | Satellite repeats (>50 bp typical) | All (1 bp – 100 kbp) |
| Telomere detection | No (short repeats not targeted) | Yes (Tier 1: 1–9 bp) |
| Template-free | Partial (classifies without templates, but templates improve accuracy) | Yes |
| Language | R | Python + required Cython for complete detection |

## Comparison on Arabidopsis Col-CEN (132 Mb)

### Data Sources

- **TRASH (figshare)**: Pre-computed from [figshare](https://doi.org/10.6084/m9.figshare.22250185) (manuscript default settings with CEN178, CEN159, 5S rDNA templates)
- **TRASH (our run)**: Run locally with default settings, no templates, 8 cores
- **bwtandem**: v3 run with all tiers, default parameters, 8 threads

### Detection Summary

| Metric | TRASH (figshare) | TRASH (our run) | bwtandem |
|--------|-----------------|-----------------|----------|
| Total repeat calls | 96,793 monomers | 79,130 monomers | 32,509 arrays |
| CEN178 monomers | 64,791 | — (no templates) | — (detects as arrays) |
| CEN159 monomers | 1,479 | — | — |
| 5S rDNA | 1,262 | — | — |
| Unclassified | 29,261 | 79,130 (all) | — |
| Short repeats (1–9 bp) | — | — | 9,214 |
| Medium repeats (10–99 bp) | — | — | 19,762 |
| Long repeats (≥100 bp) | — | — | 3,533 |
| **Runtime** | not reported | **~2h 10min** | **~26 min** |

Note: TRASH counts individual repeat monomers (~178 bp each), while bwtandem counts continuous repeat arrays (which may contain many monomers). TRASH without templates found 18% fewer monomers than with templates, as template-guided classification recovers additional divergent copies.

### Per-base Coverage Comparison (All Repeats)

Three-way per-base coverage across all 5 nuclear chromosomes:

| Chr | TRASH (figshare) | TRASH (our run) | bwtandem | Fig∩Our | Fig∩bwt | Our∩bwt |
|-----|-----------------|-----------------|----------|---------|---------|---------|
| Chr1 | 2,717,291 | 2,580,675 | 3,261,891 | 2,531,400 | 2,588,790 | 2,510,155 |
| Chr2 | 2,375,420 | 2,317,318 | 2,809,075 | 2,294,249 | 2,334,346 | 2,299,938 |
| Chr3 | 2,511,689 | 2,454,798 | 2,970,718 | 2,425,479 | 2,446,308 | 2,408,349 |
| Chr4 | 2,921,037 | 2,864,382 | 3,372,094 | 2,832,259 | 2,858,341 | 2,821,840 |
| Chr5 | 2,609,074 | 2,544,753 | 3,129,527 | 2,516,962 | 2,528,266 | 2,488,916 |

### Centromere Satellite Coverage (TRASH CEN178+CEN159 vs bwtandem)

Per-base comparison of TRASH figshare CEN178/CEN159 annotations vs bwtandem (period ≥ 100 bp):

| Chr | TRASH (bp) | bwtandem (bp) | Shared (bp) | TRASH-only (bp) | bwt-only (bp) | Jaccard |
|-----|-----------|--------------|-------------|-----------------|---------------|---------|
| Chr1 | 2,477,280 | 2,768,790 | 2,462,930 | 14,350 | 305,860 | 88.5% |
| Chr2 | 2,176,321 | 2,408,593 | 2,170,338 | 5,983 | 238,255 | 89.9% |
| Chr3 | 2,094,416 | 2,602,610 | 2,081,518 | 12,898 | 521,092 | 79.6% |
| Chr4 | 2,748,442 | 3,051,011 | 2,724,223 | 24,219 | 326,788 | 88.6% |
| Chr5 | 2,175,603 | 2,713,845 | 2,147,504 | 28,099 | 566,341 | 78.3% |
| **Total** | **11,672,062** | **13,544,849** | **11,586,513** | **85,549** | **1,958,336** | **85.0%** |

### Key Findings

1. **bwtandem covers 99.3% of TRASH's centromere bases** (11,586,513 / 11,672,062). Only 85,549 bp (0.7%) of TRASH's CEN178/CEN159 annotations are missed by bwtandem.

2. **bwtandem detects 1.96 Mb additional centromere-like sequence** not in TRASH's annotations. This likely includes:
   - Degenerate centromere copies below TRASH's classification threshold
   - Inter-monomer gap regions within continuous arrays
   - Variant-length monomers not matching CEN178/CEN159 templates

3. **85% Jaccard overlap** across all 5 chromosomes indicates strong agreement between two fundamentally different approaches (k-mer frequency vs. BWT/FM-index).

4. **bwtandem was ~5x faster in these runs** (~26 min vs ~2h 10min on 8 cores) and requires no template sequences.

5. **TRASH without templates finds 18% fewer monomers** (79,130 vs 96,793), showing that template-guided classification is important for TRASH's completeness. bwtandem achieves high coverage without any templates.

6. **Complementary strengths**:
   - TRASH excels at monomer decomposition and family classification — essential for studying centromere evolution and HOR structure
   - bwtandem excels at comprehensive detection without prior knowledge — essential for novel genomes where repeat families are unknown

## When to Use Each Tool

| Use Case | Recommended Tool |
|----------|-----------------|
| Centromere satellite monomer annotation | TRASH |
| Higher-order repeat (HOR) analysis | TRASH |
| Repeat family classification with templates | TRASH |
| De novo tandem repeat discovery (no templates) | bwtandem |
| Full-spectrum repeat detection (telomere to satellite) | bwtandem |
| Novel/non-model organism genomes | bwtandem |
| Quick whole-genome scan | bwtandem |
| Detailed centromere structural analysis | TRASH + bwtandem |

## Reproduction

### TRASH

```bash
# Install
git clone https://github.com/vlothec/TRASH
cd TRASH && bash TRASH_install.sh --def

# Run on Col-CEN
bash TRASH_run.sh ColCEN.fasta --o /output/path --par 8 --def

# With CEN178 template classification
bash TRASH_run.sh ColCEN.fasta --o /output/path --par 8 --seqt template.csv --def
```

### bwtandem

```bash
python3 -m src.main ColCEN.fasta --tiers tier1,tier2,tier3 --format bed -t 8 -v
```

## References

- Wlodzimierz P, Hong M, Henderson IR (2023). TRASH: Tandem Repeat Annotation and Structural Hierarchy. *Bioinformatics*, 39(5):btad308. [DOI:10.1093/bioinformatics/btad308](https://doi.org/10.1093/bioinformatics/btad308)
- TRASH GitHub: https://github.com/vlothec/TRASH
- TRASH figshare data: [10.6084/m9.figshare.22250185](https://doi.org/10.6084/m9.figshare.22250185)
