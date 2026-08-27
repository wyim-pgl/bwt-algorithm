"""Bin-boundary behavior for scripts/venn_compare.py.

Regions are half-open [start, end). Deriving the last bin from `end` instead of
`end - 1` credited an extra bin to every region ending exactly on a bin
boundary, inflating that tool's area in the Venn diagram.
"""
import importlib.util
import pathlib
import re

import pytest

SCRIPT = pathlib.Path(__file__).resolve().parent.parent / "scripts" / "venn_compare.py"


def _bin_regions():
    """Load bin_regions without importing matplotlib at module scope."""
    src = SCRIPT.read_text()
    m = re.search(r"def bin_regions.*?\n    return bins", src, re.S)
    assert m, "bin_regions not found in venn_compare.py"
    ns = {}
    exec(compile(m.group(0), str(SCRIPT), "exec"), ns)
    return ns["bin_regions"]


BIN = 500


@pytest.fixture(scope="module")
def bin_regions():
    return _bin_regions()


def _bins(fn, start, end, size=BIN):
    return sorted(p for _, p in fn([("c", start, end)], size))


def test_region_ending_on_a_bin_boundary_does_not_claim_the_next_bin(bin_regions):
    assert _bins(bin_regions, 0, BIN) == [0]
    assert _bins(bin_regions, 0, 2 * BIN) == [0, 1]
    assert _bins(bin_regions, BIN, 2 * BIN) == [1]


def test_region_crossing_a_boundary_claims_both(bin_regions):
    assert _bins(bin_regions, BIN - 1, BIN + 1) == [0, 1]


def test_single_base_regions(bin_regions):
    assert _bins(bin_regions, 0, 1) == [0]
    assert _bins(bin_regions, BIN - 1, BIN) == [0]
    assert _bins(bin_regions, BIN, BIN + 1) == [1]


def test_empty_and_inverted_regions_claim_nothing(bin_regions):
    assert _bins(bin_regions, 0, 0) == []
    assert _bins(bin_regions, 100, 100) == []
    assert _bins(bin_regions, 200, 100) == []


def test_bins_match_the_bases_actually_covered(bin_regions):
    """Brute force: a region must claim exactly the bins its bases fall in."""
    for start, end in ((0, 1), (0, 500), (0, 501), (499, 1000),
                       (500, 1500), (123, 4567), (999, 1000)):
        expected = sorted({pos // BIN for pos in range(start, end)})
        assert _bins(bin_regions, start, end) == expected, (start, end)


def test_regions_are_deduplicated_across_inputs(bin_regions):
    fn = bin_regions
    assert len(fn([("c", 0, 1000), ("c", 0, 1000)], BIN)) == 2
    assert len(fn([("c", 0, 1000), ("d", 0, 1000)], BIN)) == 4
