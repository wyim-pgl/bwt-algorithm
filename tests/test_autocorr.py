"""Unit tests for the shared autocorrelation primitives (src/autocorr.py).

These helpers are extracted from three previously-duplicated inline blocks:
  - finder._fill_satellite_gaps  (scalar identity per window/period)
  - finder._catchall_periodicity_fill  (cumsum windowed counts + run extraction)
  - tier2._autocorr_seed  (identical cumsum windowed counts + run extraction)

The refactor must be BEHAVIOR-PRESERVING, so the key tests assert the helpers
reproduce the exact numpy computation the inline code performed (same int64
cumsum, same int8-view run detection).
"""
import numpy as np

from src.autocorr import (
    DEFAULT_MIN_VALID_FRAC,
    autocorr_identity,
    valid_base_mask,
    windowed_match_counts,
    windowed_valid_counts,
    contiguous_true_runs,
)

# Bases 65..68 are A, B, C, D -- B and D are not nucleotides, so an array drawn
# from that range exercises only the raw-arithmetic path. Real sequence is ACGT.
ACGT = np.array([65, 67, 71, 84], dtype=np.uint8)


def _arr(s):
    return np.frombuffer(s.encode(), dtype=np.uint8)


# --- autocorr_identity (scalar) ---

def test_autocorr_identity_perfect_period():
    # "ATATAT" is perfectly periodic at period 2 -> identity 1.0
    a = _arr("ATATATATAT")
    assert autocorr_identity(a, 2) == 1.0


def test_autocorr_identity_partial():
    # period 1 over "AAAAC": compares a[:4]=AAAA vs a[1:]=AAAC -> 3/4 match
    a = _arr("AAAAC")
    assert autocorr_identity(a, 1) == 0.75


def test_autocorr_identity_period_too_long_is_zero():
    a = _arr("ACGT")
    assert autocorr_identity(a, 4) == 0.0
    assert autocorr_identity(a, 10) == 0.0


def test_autocorr_identity_matches_bruteforce():
    rng = np.random.RandomState(0)
    a = ACGT[rng.randint(0, 4, 500)]
    for p in (1, 3, 7, 13):
        total = a.size - p
        expected = float(np.sum(a[:total] == a[p:]) / total)
        assert autocorr_identity(a, p) == expected


def test_autocorr_identity_raw_mode_is_plain_arithmetic():
    """ignore_ambiguous=False keeps the original mean-of-equality contract."""
    rng = np.random.RandomState(0)
    a = rng.randint(65, 69, 500).astype(np.uint8)  # A, B, C, D
    for p in (1, 3, 7, 13):
        total = a.size - p
        expected = float(np.sum(a[:total] == a[p:]) / total)
        assert autocorr_identity(a, p, ignore_ambiguous=False) == expected


# --- ambiguous bases must not manufacture periodicity ---

def test_all_n_run_is_not_periodic():
    """N == N used to score 1.0 at every period, which emitted solid-N motifs."""
    a = _arr("N" * 400)
    for p in (1, 10, 100, 178):
        assert autocorr_identity(a, p) == 0.0
        assert autocorr_identity(a, p, ignore_ambiguous=False) == 1.0


def test_n_rich_window_with_a_real_seed_is_refused():
    # A genuine period-4 stretch followed by a long gap: the motif is ACGT but
    # the window as a whole cannot be judged.
    a = _arr("ACGT" * 10 + "N" * 160)
    assert autocorr_identity(a, 4) == 0.0


def test_identity_is_scored_over_valid_positions_only():
    # Period 1 over "AAAAA N CCCCC" gives 10 comparisons, of which (A,N) and
    # (N,C) are unusable. Raw arithmetic scores those two as mismatches (8/10);
    # masking drops them from both numerator and denominator (8/8), which is
    # the honest answer -- the gap is no evidence either way.
    a = _arr("AAAAA" + "N" + "CCCCC")
    assert autocorr_identity(a, 1, ignore_ambiguous=False) == 0.8
    assert autocorr_identity(a, 1, min_valid_frac=0.0) == 1.0


def test_min_valid_frac_gate_is_configurable():
    a = _arr("ACGT" * 10 + "N" * 20)
    assert autocorr_identity(a, 4) == 0.0            # default 0.8 gate refuses
    assert autocorr_identity(a, 4, min_valid_frac=0.0) > 0.9  # ungated scores


def test_sentinel_is_not_a_valid_base():
    assert not valid_base_mask(_arr("$"))[0]
    assert valid_base_mask(_arr("ACGT")).all()
    assert not valid_base_mask(_arr("NRYWacgt")).any()


def test_windowed_counts_exclude_ambiguous_matches():
    a = _arr("N" * 400)
    got = windowed_match_counts(a, 10, 50)
    assert got is not None and got.max() == 0
    raw = windowed_match_counts(a, 10, 50, ignore_ambiguous=False)
    assert raw is not None and raw.max() == 50


def test_windowed_valid_counts_aligns_with_match_counts():
    a = _arr("ACGT" * 100 + "N" * 200)
    m = windowed_match_counts(a, 4, 40)
    v = windowed_valid_counts(a, 4, 40)
    assert m is not None and v is not None
    assert m.shape == v.shape
    assert v[0] == 40          # fully real sequence at the start
    assert v[-1] == 0          # deep inside the gap
    assert DEFAULT_MIN_VALID_FRAC == 0.8


# --- windowed_match_counts (vector) ---

def test_windowed_match_counts_reproduces_inline_block():
    """Must equal the exact inline computation from finder/tier2 on ACGT input."""
    rng = np.random.RandomState(42)
    a = ACGT[rng.randint(0, 4, 1000)]
    n = a.size
    for period, window in ((1, 12), (5, 20), (13, 52), (20, 40)):
        s = a[:n]
        eq = (s[:n - period] == s[period:n])
        cs = np.empty(eq.size + 1, dtype=np.int64)
        cs[0] = 0
        np.cumsum(eq, out=cs[1:])
        m = eq.size - window
        expected = cs[window:window + m] - cs[:m] if m > 0 else None
        got = windowed_match_counts(a, period, window)
        if expected is None:
            assert got is None
        else:
            assert got is not None
            np.testing.assert_array_equal(got, expected)


def test_windowed_match_counts_window_semantics():
    # perfectly periodic: every window of any valid size is a full match
    a = _arr("AC" * 100)  # period-2 perfect
    got = windowed_match_counts(a, 2, 10)
    assert got is not None
    assert np.all(got == 10)


def test_windowed_match_counts_too_short_returns_none():
    a = _arr("ACGTACGT")  # 8 bp
    # period 2, window 20: eq has size 6, m = 6-20 < 0 -> None
    assert windowed_match_counts(a, 2, 20) is None


# --- contiguous_true_runs ---

def test_contiguous_true_runs_basic():
    mask = np.array([0, 1, 1, 0, 0, 1, 0], dtype=bool)
    rs, re = contiguous_true_runs(mask)
    np.testing.assert_array_equal(rs, [1, 5])
    np.testing.assert_array_equal(re, [3, 6])  # exclusive ends


def test_contiguous_true_runs_boundaries():
    # True at both ends must be captured
    mask = np.array([1, 1, 0, 1], dtype=bool)
    rs, re = contiguous_true_runs(mask)
    np.testing.assert_array_equal(rs, [0, 3])
    np.testing.assert_array_equal(re, [2, 4])


def test_contiguous_true_runs_all_false():
    mask = np.zeros(5, dtype=bool)
    rs, re = contiguous_true_runs(mask)
    assert len(rs) == 0
    assert len(re) == 0


def test_contiguous_true_runs_all_true():
    mask = np.ones(4, dtype=bool)
    rs, re = contiguous_true_runs(mask)
    np.testing.assert_array_equal(rs, [0])
    np.testing.assert_array_equal(re, [4])


def test_contiguous_true_runs_reproduces_inline_block():
    """Must equal the exact np.diff(int8-view) inline computation."""
    rng = np.random.RandomState(7)
    for _ in range(20):
        hit = rng.rand(200) > 0.5
        d = np.diff(hit.view(np.int8))
        run_s = np.nonzero(d == 1)[0] + 1
        run_e = np.nonzero(d == -1)[0] + 1
        if hit[0]:
            run_s = np.concatenate(([0], run_s))
        if hit[-1]:
            run_e = np.concatenate((run_e, [hit.size]))
        gs, ge = contiguous_true_runs(hit)
        np.testing.assert_array_equal(gs, run_s)
        np.testing.assert_array_equal(ge, run_e)
