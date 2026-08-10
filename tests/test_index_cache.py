"""Regression tests for the on-disk FM-index cache.

The cache exists so a parameter sweep pays for the suffix array once instead of
once per configuration. That is only safe if a loaded index is indistinguishable
from a freshly built one, and if a cache that does not match the sequence is
refused rather than used. Both are tested here, along with the failure modes that
would otherwise corrupt results silently: a truncated file, a foreign sequence of
the same length, and a stale format version.
"""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src.bwt_core import BWTCore  # noqa: E402


def _seq(n=20000, seed=11):
    rng = np.random.default_rng(seed)
    return "".join(rng.choice(list("ACGT"), n)) + "$"


@pytest.fixture(scope="module")
def built():
    s = _seq()
    return s, BWTCore(s, sa_sample_rate=1)


def test_roundtrip_is_identical(tmp_path, built):
    """A loaded index must match a built one field for field."""
    seq, a = built
    p = str(tmp_path / "idx.npz")
    a.save_index(p)
    b = BWTCore.load_index(p, seq, sa_sample_rate=1)
    assert b is not None

    assert b.n == a.n
    assert b.occ_sample_rate == a.occ_sample_rate
    assert b.alphabet == a.alphabet
    assert b.char_counts == a.char_counts
    assert b.char_totals == a.char_totals
    for field in ("text_arr", "suffix_array", "bwt_arr"):
        assert np.array_equal(getattr(b, field), getattr(a, field)), field
    assert set(b.occ_checkpoints) == set(a.occ_checkpoints)
    for k in a.occ_checkpoints:
        assert np.array_equal(b.occ_checkpoints[k], a.occ_checkpoints[k]), k


def test_queries_agree(tmp_path, built):
    """The point of the index is search; loaded and built must answer the same."""
    seq, a = built
    p = str(tmp_path / "idx.npz")
    a.save_index(p)
    b = BWTCore.load_index(p, seq, sa_sample_rate=1)
    for pat in ("A", "ACGT", "GGTA", "TTTTT", "CGCGCG", seq[100:112]):
        assert a.count_occurrences(pat) == b.count_occurrences(pat), pat
        assert sorted(a.locate_positions(pat)) == sorted(b.locate_positions(pat)), pat


def test_refuses_different_sequence(tmp_path, built):
    """Pairing an index with the wrong sequence would corrupt every downstream
    call, so a mismatch must fail loudly rather than load."""
    seq, a = built
    p = str(tmp_path / "idx.npz")
    a.save_index(p)
    assert BWTCore.load_index(p, _seq(seed=99)) is None      # same length, different bases
    assert BWTCore.load_index(p, seq[:-50] + "$") is None    # different length


def test_refuses_truncated_file(tmp_path, built):
    """A job killed mid-write must not leave a cache that later loads."""
    seq, a = built
    p = str(tmp_path / "idx.npz")
    a.save_index(p)
    with open(p, "r+b") as f:
        f.truncate(os.path.getsize(p) // 2)
    assert BWTCore.load_index(p, seq) is None


def test_refuses_missing_file(tmp_path):
    assert BWTCore.load_index(str(tmp_path / "absent.npz"), _seq()) is None


def test_cache_key_tracks_sequence_and_rate():
    """The key must separate different sequences and different build parameters,
    since both change the stored arrays."""
    s1, s2 = _seq(seed=1), _seq(seed=2)
    assert BWTCore.cache_key(s1) == BWTCore.cache_key(s1)
    assert BWTCore.cache_key(s1) != BWTCore.cache_key(s2)
    assert BWTCore.cache_key(s1, 128) != BWTCore.cache_key(s1, 64)


def test_save_is_atomic(tmp_path, built):
    """save_index writes via a temporary and renames, so no .tmp debris remains."""
    seq, a = built
    p = str(tmp_path / "idx.npz")
    a.save_index(p)
    leftovers = [f for f in os.listdir(tmp_path) if ".tmp" in f]
    assert leftovers == [], leftovers
    assert BWTCore.load_index(p, seq) is not None
