# Superseded documentation

Kept for the record, superseded 2026-09-01 by the rewritten top-level
[README](../../README.md) and [MANUAL](../../MANUAL.md):

- `setup.md` — early conda setup + ULTRA build notes. Predates packaging:
  installation is now `pip install .` inside the env; the cythonize command
  it shows is built into `setup.py`.
- `SETUP_LINUX.md` — native Linux/HPC recipe from before the packaging work.
  ⚠️ Its claim that the pure-Python fallbacks "return empty results and the
  tool reports 0 repeats" describes a pre-2026-07 defect; the fallbacks have
  been faithful (parity-tested byte-identical) since the accelerator-parity
  campaign.
- `code_update.md` — an early development log of the Tier 1 sliding-window
  conversion; project history, not current documentation.
