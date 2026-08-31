import numpy as np
import pytest

from src.finder import TandemRepeatFinder
from src.models import TandemRepeat


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


def test_invalid_tier_selection_is_rejected():
    with pytest.raises(ValueError, match="No valid tiers selected"):
        TandemRepeatFinder("ACGT" * 100, enabled_tiers={"typo"})


def test_satellite_scan_rejects_ambiguous_only_sequence():
    # Seed a nearby long-period call so the proximity gate is active; an N-only
    # gap must still not be accepted as a perfectly autocorrelated satellite.
    finder = TandemRepeatFinder(
        "A" * 400 + "N" * 1000 + "C" * 400,
        min_period=100,
        max_period=300,
        enabled_tiers={"tier3"},
    )
    try:
        seed = TandemRepeat(
            chrom="chr1", start=0, end=400, motif="A" * 100,
            copies=4.0, length=400, tier=3,
        )
        filled = finder._fill_satellite_gaps([seed])
        assert len(filled) == 1
    finally:
        finder.cleanup()


def test_trf_dat_uses_one_based_start_and_inclusive_end():
    repeat = TandemRepeat(
        chrom="chr1", start=9, end=29, motif="AC", copies=10.0,
        length=20, tier=1,
    )
    fields = repeat.to_trf_dat().split()
    assert fields[:3] == ["10", "29", "2"]
