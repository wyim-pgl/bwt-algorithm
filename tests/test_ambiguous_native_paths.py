"""An ambiguous base must not be evidence of periodicity in the native paths.

`bwtandem/autocorr.py` was corrected first, but the detector's hot paths compare
bases directly and counted `N == N` as agreement, so a seed could extend
straight through an assembly gap or a masked block. Two of these are worse than
raw equality:

  * the 2-bit packing in `_accelerators.pyx` has no code for an ambiguous base
    -- `_base_map` leaves N, `$` and every IUPAC code at 0, which is also A's
    code -- so the packed comparator scored `N == A` as a match and disagreed
    with the byte comparator on the same input;
  * `tier1_scan.c` checked only the first k bases of the motif for N while the
    run itself extended on raw equality, so a run starting on real sequence
    could end deep inside a gap.

The rule pinned here is that an ambiguous base matches nothing, not even
itself. The .pyx and the Python fallback must apply it identically or the
parity tests compare two different definitions.
"""
import numpy as np
import pytest

from bwtandem import accelerators as A

_native = pytest.importorskip("bwtandem._accelerators")


def _arr(s):
    return np.frombuffer(s.encode(), dtype=np.uint8)


# --- hamming ---------------------------------------------------------------

@pytest.mark.parametrize("s1,s2,expected", [
    ("ACGT", "ACGT", 0),
    ("AAAAAAAA", "AAAAAAAA", 0),      # 64-bit SWAR path
    ("ACGTACGTACGTACGT", "ACGTACGTACGTACGT", 0),
    ("NNNN", "NNNN", 4),              # ambiguous matches nothing, not itself
    ("NNNNNNNN", "NNNNNNNN", 8),      # ... including on the SWAR path
    ("AAAA", "NNNN", 4),
    ("ACNT", "ACNT", 1),              # only the ambiguous column
    ("ACGT", "ACGA", 1),
    ("AC$T", "AC$T", 1),              # the sentinel is ambiguous too
])
def test_hamming_charges_ambiguous_positions(s1, s2, expected):
    assert _native.hamming_distance(_arr(s1), _arr(s2)) == expected


def test_hamming_is_order_independent_of_pack_sequence():
    """The rule must not depend on some other call having initialised a table.

    An earlier attempt used a lookup table initialised inside `pack_sequence`;
    calling `hamming_distance` first then saw every base as ambiguous and
    scored ACGT vs ACGT as 4 mismatches.
    """
    assert _native.hamming_distance(_arr("ACGT"), _arr("ACGT")) == 0


# --- packed vs byte comparator agreement -----------------------------------

@pytest.mark.parametrize("seq,unit", [
    ("A" * 40 + "N" * 40, 20),
    ("ACGTACGTAC" * 3 + "N" * 30 + "ACGTACGTAC" * 3, 10),
    ("AC" * 60, 4),
    ("ACGT" * 20 + "$" + "ACGT" * 20, 4),
])
def test_packed_and_byte_scan_agree(seq, unit):
    arr = _arr(seq)
    packed = _native.pack_sequence(arr)
    byte_path = _native.scan_unit_repeats(arr, len(seq), unit, 2, 2, None)
    packed_path = _native.scan_unit_repeats(arr, len(seq), unit, 2, 2, packed)
    assert byte_path == packed_path


def test_a_run_does_not_extend_into_a_gap():
    """A*40 + N*40 at unit 20 was returned as one call spanning both halves."""
    seq = "A" * 40 + "N" * 40
    arr = _arr(seq)
    for packed in (None, _native.pack_sequence(arr)):
        calls = _native.scan_unit_repeats(arr, len(seq), 20, 2, 2, packed)
        for start, end in calls:
            assert end <= 40, f"call {start}-{end} runs into the N block"


# --- extension and anchor scanning -----------------------------------------

def _extend(fn, seq, seed, period):
    return fn(_arr(seq), seed, period, len(seq), 0.2)


def test_extension_stops_at_a_gap_in_both_implementations():
    """Extension must not run the length of a gap.

    It is mismatch-tolerant, so the global budget still absorbs one all-N copy
    before the run breaks -- 60 bp of real array plus one 10 bp period. What it
    must not do is treat the whole gap as agreement and consume all 60 bp of it.
    """
    seq = "ACGTACGTAC" * 6 + "N" * 60          # 60 bp of period-10, then a gap
    native = _extend(_native.extend_with_mismatches, seq, 0, 10)
    fallback = _extend(A._extend_with_mismatches_py, seq, 0, 10)
    assert native == fallback
    assert native[1] <= 70, f"extended to {native[1]}, {native[1] - 60} bp into the gap"
    assert native[1] < len(seq)


def test_anchor_scan_stops_at_a_gap_in_both_implementations():
    seq = "ACGTACGTAC" * 6 + "N" * 60
    native = _native.anchor_scan_boundaries(_arr(seq), 0, 10, len(seq), 0.7, 5, 20)
    fallback = A._anchor_scan_boundaries_py(_arr(seq), 0, 10, len(seq), 0.7, 5, 20)
    assert native == fallback
    assert native[1] <= 60, f"anchor scan reached {native[1]}, into the gap at 60"


def test_an_all_ambiguous_region_supports_nothing():
    seq = "N" * 200
    arr = _arr(seq)
    native = _native.extend_with_mismatches(arr, 0, 10, len(seq), 0.2)
    fallback = A._extend_with_mismatches_py(arr, 0, 10, len(seq), 0.2)
    assert native == fallback
    # One period is the seed itself; nothing beyond it may be claimed.
    assert native[2] == 1, f"{native[2]} copies claimed inside a pure N run"
