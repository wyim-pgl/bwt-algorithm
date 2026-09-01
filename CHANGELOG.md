# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.9.0] — pending first tagged release

First versioned release candidate. Everything below describes the state of the
codebase at this version; there is no earlier tagged release to diff against.

### Added
- 3-tier BWT/FM-index tandem repeat detection pipeline (`src/`): Tier 1
  (short STRs, motifs 1–9 bp, FM-index enumeration or sliding-window scan),
  Tier 2 (medium/imperfect repeats ≥10 bp, LCP + k-mer seeding), Tier 3
  (long repeats, periods 100 bp–100 kbp, sparse k-mer seeding).
- Output formats: BED (default), VCF, TRF `.dat`, STRfinder CSV.
- Cython accelerator extension (`src/_accelerators.pyx`) plus four C
  extensions, each with a faithful pure-Python fallback pinned to
  byte-identical output by parity tests (`BWT_DISABLE_NATIVE=1`).
- Satellite interior gap-fill backstop for divergent higher-order-repeat
  arrays, and an opt-in catch-all periodicity pass for diverged STRs
  (`CATCHALL_SCAN=1`).
- FM-index cache keyed by sequence hash (`BWT_INDEX_CACHE`).
- Environment-variable overrides for all detection thresholds, pinned to the
  documentation by `tests/test_env_var_docs.py`; unset always reproduces
  baseline behaviour.
- Package metadata (`pyproject.toml`, wheel ships the compiled accelerator),
  `bwtandem` console entry point, MIT `LICENSE`, `CITATION.cff`.
- CI: native build + full test suite, pure-Python parity leg, CLI smoke test,
  and a clean-environment wheel-install job asserting the compiled extension
  loads from `site-packages`.
- Benchmark provenance under `results/`: per-table manifest with hashes,
  deposited whole-genome BEDs, ground-truth coordinate sets, and audit
  evidence.

### Fixed
- Uninitialised traceback pointer table in the C alignment path that made
  C and Python alignment disagree on 31% of regions (2026-07-09).
- Four output-correctness bugs: STRfinder CSV field count, C align path
  dropping per-copy detail, satellite gap-fill running past the sentinel and
  emitting non-ACGT motifs, non-integer BED tier column (2026-08-06).
- Half-open interval binning at exact Venn bin boundaries in
  `scripts/venn_compare.py` (#23).
- Wheel packaging: packages declared explicitly so the built wheel ships
  `src` with the native extension instead of a `py3-none-any` wheel (#14).
