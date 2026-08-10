"""The C alignment accelerator and the pure-Python loop must agree.

`MotifUtils.align_repeat_region` runs `libalign_accel`'s C loop when it loads and
its own Python loop otherwise, so the two must return the same summary or the
tool's answer depends on whether a build artefact happens to be present.

They did not. Measured 2026-07-09 on 1500 random repeat regions, they disagreed
on 31% of them, from three causes:

  * `align_accel.c` read uninitialised heap in its traceback (column 0 of the
    pointer table was written to a scratch row, not to the table the traceback
    reads), which also made the whole caller nondeterministic;
  * its extension bound was `end + 3*motif_len + 4*max_indel` where the Python
    loop used `max(end, start + motif_len*min_copies) + max(3*motif_len, 4*max_indel)`;
  * `_consensus_from_counts` broke ties by insertion order while the C — and
    `build_consensus_motif_array` — take the lexicographically smallest base.
    A different consensus cascades into every later copy's alignment.

With all three fixed the disagreement is zero. This test keeps it there.
"""
import numpy as np
import pytest

import src.motif_utils as mu
from src.motif_utils import MotifUtils

pytestmark = pytest.mark.skipif(
    mu._c_align_lib is None,
    reason="libalign_accel is not loaded; there is only one implementation to compare",
)

BASES = "ACGT"
# align_repeat_region caches its C text pointer by id(sequence); hold every
# sequence so CPython cannot recycle an id into a stale cache entry.
_KEEP = []


def summary_of(r):
    if r is None:
        return None
    return (r.consensus, r.motif_len, r.copies, r.consumed_length,
            round(r.mismatch_rate, 9), r.total_insertions, r.total_deletions)


def full_summary_of(r):
    """Every field of the summary, including the per-copy detail.

    `variations` reaches the user as the STRfinder Variations column;
    `copy_sequences` and `error_counts` are the per-copy record behind it. The C
    path returned all three empty while the Python path filled them, so what the
    tool reported depended on whether a build artefact was present.
    """
    if r is None:
        return None
    return summary_of(r) + (r.max_errors_per_copy, r.total_mismatches,
                            tuple(r.variations), tuple(r.copy_sequences),
                            tuple(r.error_counts))


def align(seq, start, end, motif, mism, max_indel, min_copies, use_c, summarize=summary_of):
    saved = mu._c_align_lib
    if not use_c:
        mu._c_align_lib = None
    try:
        return summarize(MotifUtils.align_repeat_region(
            seq, start, end, motif, mismatch_fraction=mism,
            max_indel=max_indel, min_copies=min_copies))
    finally:
        mu._c_align_lib = saved


def random_region(rng):
    """A tandem array carrying substitutions, insertions and deletions."""
    period = int(rng.integers(2, 30))
    n_copies = int(rng.integers(2, 12))
    unit = "".join(BASES[i] for i in rng.integers(0, 4, size=period))

    copies = []
    for _ in range(n_copies):
        c = list(unit)
        mode = rng.random()
        if mode < 0.30 and period > 3:
            k = int(rng.integers(1, min(4, period)))
            pos = int(rng.integers(0, period - k))
            del c[pos:pos + k]
        elif mode < 0.55:
            k = int(rng.integers(1, 4))
            pos = int(rng.integers(0, len(c) + 1))
            c[pos:pos] = [BASES[i] for i in rng.integers(0, 4, size=k)]
        elif mode < 0.85:
            for _ in range(int(rng.integers(1, 4))):
                p = int(rng.integers(0, len(c)))
                c[p] = BASES[int(rng.integers(0, 4))]
        copies.append("".join(c))

    left = "".join(BASES[i] for i in rng.integers(0, 4, size=int(rng.integers(0, 8))))
    right = "".join(BASES[i] for i in rng.integers(0, 4, size=int(rng.integers(0, 8))))
    seq = left + "".join(copies) + right
    _KEEP.append(seq)

    max_indel = rng.choice([None, 0, 1, 2, 3, 5])
    return (seq, len(left), len(seq) - len(right), unit,
            float(rng.choice([0.0, 0.05, 0.1, 0.2, 0.35])),
            None if max_indel is None else int(max_indel),
            int(rng.integers(1, 4)))


def test_c_and_python_align_loops_agree():
    rng = np.random.default_rng(20260709)
    mismatches = []
    for trial in range(600):
        seq, start, end, motif, mism, max_indel, min_copies = random_region(rng)
        c = align(seq, start, end, motif, mism, max_indel, min_copies, use_c=True)
        p = align(seq, start, end, motif, mism, max_indel, min_copies, use_c=False)
        if c != p:
            mismatches.append((trial, motif, start, end, mism, max_indel, min_copies, c, p))

    assert not mismatches, (
        f"{len(mismatches)}/600 regions align differently with and without "
        f"libalign_accel; first: {mismatches[0]}"
    )


def test_c_and_python_per_copy_detail_agrees():
    """The same comparison over every field, per-copy detail included."""
    rng = np.random.default_rng(20260805)
    mismatches = []
    for trial in range(600):
        seq, start, end, motif, mism, max_indel, min_copies = random_region(rng)
        c = align(seq, start, end, motif, mism, max_indel, min_copies,
                  use_c=True, summarize=full_summary_of)
        p = align(seq, start, end, motif, mism, max_indel, min_copies,
                  use_c=False, summarize=full_summary_of)
        if c != p:
            mismatches.append((trial, motif, start, end, mism, max_indel, min_copies, c, p))

    assert not mismatches, (
        f"{len(mismatches)}/600 regions report different per-copy detail with and "
        f"without libalign_accel; first: {mismatches[0]}"
    )


def test_the_c_path_actually_reports_detail():
    """Guard the original defect: parity alone passes if both paths go empty."""
    rng = np.random.default_rng(4242)
    with_variations = 0
    for _ in range(200):
        seq, start, end, motif, mism, max_indel, min_copies = random_region(rng)
        r = MotifUtils.align_repeat_region(seq, start, end, motif,
                                           mismatch_fraction=mism,
                                           max_indel=max_indel,
                                           min_copies=min_copies)
        if r is None:
            continue
        assert len(r.copy_sequences) == r.copies
        assert len(r.error_counts) == r.copies
        if r.variations:
            with_variations += 1

    assert with_variations > 0, (
        "the C path reported no variations at all over 200 mutated repeat "
        "regions, which is what the empty-list bug looked like"
    )


def resolved_indel(motif, max_indel):
    """The indel bound align_repeat_region hands the C loop."""
    if max_indel is not None:
        return max(0, max_indel)
    m = len(motif)
    return max(1, min(10, m // 2 if m >= 4 else 1))


def test_the_c_detail_buffers_are_never_too_small():
    """The C loop cannot grow its own detail buffers, so the caller sizes them.

    An undersized bound corrupts nothing — the C reports the overflow and the
    caller re-runs the Python loop, which is the reference implementation — but
    it silently throws the acceleration away on exactly the regions with the most
    detail to report. So assert the C path never has to bail: wherever Python
    produces a summary, so must C.
    """
    rng = np.random.default_rng(99)
    bailed = []
    for trial in range(600):
        seq, start, end, motif, mism, max_indel, min_copies = random_region(rng)
        c = MotifUtils._align_repeat_region_c(
            seq, start, end, motif, mism, resolved_indel(motif, max_indel), min_copies)
        p = align(seq, start, end, motif, mism, max_indel, min_copies, use_c=False)
        if p is not None and c is None:
            bailed.append((trial, motif, start, end, mism, max_indel, min_copies))

    assert not bailed, (
        f"{len(bailed)}/600 regions fell back to the Python loop; the detail "
        f"buffer bound in _c_detail_capacities is too small. First: {bailed[0]}"
    )


def test_consensus_tie_breaks_to_the_smallest_base():
    """All three consensus builders must agree on ties, or alignments diverge."""
    from collections import Counter

    counts = [Counter({"T": 2, "A": 2}), Counter({"G": 1, "C": 1}), Counter()]
    assert MotifUtils._consensus_from_counts(counts, "NNN") == "ACN"

    # the array builder, on a 1-mer with a 1-1 tie between A and T
    arr = np.frombuffer(b"AT", dtype=np.uint8)
    consensus, _, _ = MotifUtils.build_consensus_motif_array(arr, 0, 1, 2)
    assert consensus.tobytes() == b"A"


def test_deletion_count_is_arithmetically_consistent():
    """copies * motif_len - consumed == deletions - insertions.

    The uninitialised traceback used to report zero deletions for alignments that
    demonstrably had them (17 copies of a 3-mer consuming 44 bases needs 7).
    """
    seq = "TCATGT" + "ACG" * 10 + "TGCATGCATGCATGC"
    _KEEP.append(seq)
    r = MotifUtils.align_repeat_region(seq, 6, 36, "ACG",
                                       mismatch_fraction=0.1, min_copies=2)
    assert r is not None
    expected_gap = r.copies * r.motif_len - r.consumed_length
    assert expected_gap == r.total_deletions - r.total_insertions, (
        f"copies={r.copies} motif_len={r.motif_len} consumed={r.consumed_length} "
        f"del={r.total_deletions} ins={r.total_insertions}"
    )
