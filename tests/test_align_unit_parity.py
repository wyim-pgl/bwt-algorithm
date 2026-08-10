"""`align_unit_to_window`: the Cython loop and the Python banded DP must agree.

Why this exists (2026-08-06): of the five accelerator symbols the tiers consume, this
was the only one with no C-vs-Python value comparison — `test_accel_parity.py:151`
checks that the binding exists, not that it returns the same thing. That gap mattered
more than the others because this symbol's Python counterpart is a separate,
hand-written banded DP (`motif_utils.py:273-424`), the same code shape as the
`align_accel.c` loop that disagreed with its twin on 31% of random regions until
2026-07-09 (see `docs/2026-07-09-nondeterminism-uninitialised-ptr-table.md`).

The two implementations agree on all 2000 generated cases today, so this is a
regression guard rather than a bug report. Cases are drawn to hit the paths that
matter: clean repeats, single and multi-base deletions, insertions, and substitution
clusters, each followed by a random tail so the window is longer than the motif.
"""
import numpy as np
import pytest

import src.accelerators as acc
import src.motif_utils as mu
from src.motif_utils import MotifUtils

pytestmark = pytest.mark.skipif(
    not acc.NATIVE_AVAILABLE,
    reason="compiled src/_accelerators is absent; only one implementation to compare",
)

BASES = "ACGT"


def _native(motif, window, max_indel, tol):
    res = acc.align_unit_to_window(motif.encode(), window.encode(), max_indel, tol)
    if res is None:
        return None
    consumed, unit_seq, mm, ins, dele, ops, obs, ed = res
    return (consumed, unit_seq, mm, ins, dele,
            tuple(map(tuple, ops)), tuple(map(tuple, obs)), ed)


def _python(motif, window, max_indel, tol):
    """Run MotifUtils' own banded DP.

    `MotifUtils._align_unit_to_window` tries the accelerator first, so the
    module-level binding is neutralised for the duration of the call to force the
    Python body. Restored in a finally block — leaking this would silently disable
    the accelerator for every later test in the session.
    """
    saved = mu.align_unit_to_window
    mu.align_unit_to_window = lambda *a, **k: None
    try:
        r = MotifUtils._align_unit_to_window(motif, window, max_indel, tol)
    finally:
        mu.align_unit_to_window = saved
    if r is None:
        return None
    return (r.consumed, r.unit_sequence, r.mismatch_count, r.insertion_length,
            r.deletion_length, tuple(map(tuple, r.operations)),
            tuple(map(tuple, r.observed_bases)), r.edit_distance)


def _random_case(rng):
    m = int(rng.integers(2, 30))
    motif = "".join(BASES[i] for i in rng.integers(0, 4, size=m))
    w = list(motif)
    mode = rng.random()
    if mode < 0.30 and m > 3:                          # deletion
        k = int(rng.integers(1, min(4, m)))
        p = int(rng.integers(0, m - k))
        del w[p:p + k]
    elif mode < 0.55:                                  # insertion
        k = int(rng.integers(1, 4))
        p = int(rng.integers(0, len(w) + 1))
        w[p:p] = [BASES[i] for i in rng.integers(0, 4, size=k)]
    elif mode < 0.85:                                  # substitutions
        for _ in range(int(rng.integers(1, 4))):
            p = int(rng.integers(0, len(w)))
            w[p] = BASES[int(rng.integers(0, 4))]
    tail = "".join(BASES[i] for i in rng.integers(0, 4, size=int(rng.integers(0, 5))))
    window = "".join(w) + tail
    max_indel = max(1, min(10, m // 2 if m >= 4 else 1))   # mirrors align_repeat_region
    tol = max(1, int(m * 0.2))
    return motif, window, max_indel, tol


def test_align_unit_to_window_c_and_python_agree():
    rng = np.random.default_rng(20260805)
    bad = []
    for _ in range(2000):
        motif, window, max_indel, tol = _random_case(rng)
        native = _native(motif, window, max_indel, tol)
        python = _python(motif, window, max_indel, tol)
        if native != python:
            bad.append((motif, window, max_indel, tol, native, python))
    assert not bad, (
        f"{len(bad)}/2000 windows disagree between the Cython align_unit_to_window "
        f"and MotifUtils._align_unit_to_window.\nfirst: {bad[0]}"
    )
