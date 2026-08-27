"""Regression tests for CLI/pipeline constraints and output-format coordinates.

These cover fixes that were made during a pre-submission audit but reached this
branch without tests: TRF .dat coordinates, tier-selection validation, Tier 3
period wiring, argument validation, and the fail-closed multiprocess path.
"""
import os
import pathlib
import subprocess
import sys
import textwrap

import pytest

from src.finder import TandemRepeatFinder
from src.models import TandemRepeat

REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


# --- Tier 3 period wiring -------------------------------------------------

def test_tier3_honors_requested_period_range():
    finder = TandemRepeatFinder(
        "ACGT" * 300,
        min_period=250,
        max_period=400,
        enabled_tiers={"tier3"},
    )
    try:
        assert finder.tier3 is not None
        assert finder.tier3.min_length == 250
        assert finder.tier3.max_length == 400
    finally:
        finder.cleanup()


def test_tier3_floor_stays_at_100_when_min_period_is_lower():
    finder = TandemRepeatFinder(
        "ACGT" * 300,
        min_period=1,
        max_period=2000,
        enabled_tiers={"tier3"},
    )
    try:
        assert finder.tier3 is not None
        assert finder.tier3.min_length == 100
        assert finder.tier3.max_length == 2000
    finally:
        finder.cleanup()


def test_tier3_is_not_built_when_range_ends_below_100():
    finder = TandemRepeatFinder(
        "ACGT" * 100,
        min_period=1,
        max_period=9,
        enabled_tiers={"tier3"},
    )
    try:
        assert finder.tier3 is None
    finally:
        finder.cleanup()


# --- Tier selection validation -------------------------------------------

def test_all_invalid_tier_selection_is_rejected():
    # This used to fall through to "enable every tier".
    with pytest.raises(ValueError, match="Unknown tier"):
        TandemRepeatFinder("ACGT" * 100, enabled_tiers={"typo"})


def test_mixed_valid_and_invalid_tier_selection_is_rejected():
    # This used to silently run tier1 alone.
    with pytest.raises(ValueError, match="Unknown tier"):
        TandemRepeatFinder("ACGT" * 100, enabled_tiers={"tier1", "typo"})


def test_empty_token_from_trailing_comma_is_tolerated():
    finder = TandemRepeatFinder("ACGT" * 100, enabled_tiers={"tier1", ""})
    try:
        assert finder.enabled_tiers == {"tier1"}
    finally:
        finder.cleanup()


def test_all_expands_to_every_tier():
    finder = TandemRepeatFinder("ACGT" * 100, enabled_tiers={"all"})
    try:
        assert finder.enabled_tiers == {"tier1", "tier2", "tier3"}
    finally:
        finder.cleanup()


# --- Output format coordinates -------------------------------------------

def test_trf_dat_uses_one_based_start_and_inclusive_end():
    # Internal intervals are 0-based half-open; TRF .dat is 1-based inclusive,
    # so [9, 29) must be written as 10..29.
    repeat = TandemRepeat(
        chrom="chr1", start=9, end=29, motif="AC", copies=10.0,
        length=20, tier=1,
    )
    fields = repeat.to_trf_dat().split()
    assert fields[:3] == ["10", "29", "2"]


def test_bed_start_stays_zero_based():
    repeat = TandemRepeat(
        chrom="chr1", start=9, end=29, motif="AC", copies=10.0,
        length=20, tier=1,
    )
    assert repeat.to_bed().split("\t")[1] == "9"


# --- CLI argument validation ---------------------------------------------

def _run_cli(tmp_path, *args):
    fasta = tmp_path / "in.fa"
    fasta.write_text(">c\n" + "ACGT" * 100 + "\n")
    return subprocess.run(
        [sys.executable, "-m", "src.main", str(fasta),
         "-o", str(tmp_path / "out"), *args],
        capture_output=True, text=True, cwd=REPO_ROOT,
    )


@pytest.mark.parametrize("args", [
    ("--min-period", "0"),
    ("--max-period", "0"),
    ("--min-period", "-1"),
    ("--min-period", "50", "--max-period", "10"),
])
def test_cli_rejects_invalid_period_intervals(tmp_path, args):
    proc = _run_cli(tmp_path, *args)
    assert proc.returncode != 0
    assert not (tmp_path / "out.bed").exists()


@pytest.mark.parametrize("tiers", ["typo", "tier1,typo"])
def test_cli_rejects_unknown_tier_names(tmp_path, tiers):
    proc = _run_cli(tmp_path, "--tiers", tiers)
    assert proc.returncode != 0
    assert "unknown tier" in proc.stderr.lower()
    assert not (tmp_path / "out.bed").exists()


def test_cli_accepts_a_valid_run(tmp_path):
    proc = _run_cli(tmp_path, "--tiers", "tier1")
    assert proc.returncode == 0
    assert (tmp_path / "out.bed").exists()


# --- Fail-closed multiprocess path ---------------------------------------

def test_multiprocess_failure_writes_no_partial_output(tmp_path):
    """A worker exception must abort the run before any output is written."""
    fasta = tmp_path / "two.fa"
    fasta.write_text(">a\n" + "ACGT" * 200 + "\n>b\n" + "ACGT" * 200 + "\n")

    # Break _process_chromosome for one sequence only, so the other succeeds
    # and would previously have been written out on its own.
    driver = tmp_path / "driver.py"
    driver.write_text(textwrap.dedent(f"""
        import sys
        sys.path.insert(0, {str(REPO_ROOT)!r})
        sys.argv = ["src.main", {str(fasta)!r}, "-t", "2",
                    "-o", {str(tmp_path / 'out')!r}]
        import src.main as M
        _real = M._process_chromosome
        def boom(chrom, *a, **k):
            if chrom == "b":
                raise RuntimeError("injected worker fault")
            return _real(chrom, *a, **k)
        M._process_chromosome = boom
        M.main()
    """))

    proc = subprocess.run([sys.executable, str(driver)],
                          capture_output=True, text=True, cwd=REPO_ROOT)
    # The run must abort because of the injected fault, not because the driver
    # failed to import -- otherwise this test passes with the fix reverted.
    assert "injected worker fault" in proc.stderr, proc.stderr
    assert proc.returncode != 0
    assert not (tmp_path / "out.bed").exists()
    # And the successful sequence must not have been written on its own.
    assert not list(tmp_path.glob("out.*"))
