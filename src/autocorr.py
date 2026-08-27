"""Shared autocorrelation primitives for periodicity detection.

Three detection paths independently re-derived the same "how self-similar is the
sequence at offset ``period``" computation:

  - ``finder._fill_satellite_gaps``   — a scalar identity per (window, period)
  - ``finder._catchall_periodicity_fill`` — a cumsum-based sliding-window count
  - ``tier2._autocorr_seed``          — the identical cumsum sliding-window count

This module factors those out so the period-detection math lives in one place.
Every function reproduces the exact numpy operations the inline code used (same
``int64`` cumsum, same ``int8``-view run detection), so substituting them is
behavior-preserving.

All functions take a 1-D ``uint8`` view of the sequence (e.g. ``BWTCore.text_arr``).

Ambiguous bases
---------------
Sequences reach these functions uppercased, with assembly gaps and (under
``--mask soft``/``both``) masked regions written as ``N``; the BWT sentinel
``$`` can also be in view. Raw equality counts ``N == N`` as a match, so a run
of ``N`` is trivially "periodic" at *every* period and scores identity 1.0.
That is how solid-``N`` motifs were emitted at telomeres, and a motif-only
guard does not stop it: an A/C/G/T motif at the start of a window is enough to
authorize a call that then spans an ``N``-rich block.

So comparisons where either side is not A/C/G/T are excluded by default
(``ignore_ambiguous=True``) rather than counted as agreement, and
``autocorr_identity`` additionally refuses to score a window that does not have
enough valid comparisons to judge. Pass ``ignore_ambiguous=False`` for the
plain arithmetic.
"""
from typing import Optional, Tuple

import numpy as np

# uint8 lookup: True exactly for b'A', b'C', b'G', b'T'.
_VALID_BASE = np.zeros(256, dtype=bool)
_VALID_BASE[[65, 67, 71, 84]] = True

#: Default minimum fraction of comparisons that must be A/C/G/T vs A/C/G/T for
#: :func:`autocorr_identity` to return a score at all.
DEFAULT_MIN_VALID_FRAC = 0.8


def valid_base_mask(arr: np.ndarray) -> np.ndarray:
    """Boolean mask marking positions holding an unambiguous A/C/G/T base."""
    return _VALID_BASE[arr]


def autocorr_identity(arr: np.ndarray, period: int,
                      ignore_ambiguous: bool = True,
                      min_valid_frac: float = DEFAULT_MIN_VALID_FRAC) -> float:
    """Fraction of bases that match their copy ``period`` positions ahead.

    Returns ``0.0`` when ``period`` is not shorter than ``arr`` (no overlap).
    This is the scalar primitive used by the satellite gap-fill's per-window
    best-period search and by the satellite merge check.

    With ``ignore_ambiguous`` (the default), only positions where *both* bases
    are A/C/G/T are compared, and the result is the identity over those
    positions alone. If fewer than ``min_valid_frac`` of the comparisons are
    valid the window cannot be judged and ``0.0`` is returned -- an N-rich gap
    yields no evidence of periodicity rather than perfect evidence.
    """
    total = arr.size - period
    if total <= 0:
        return 0.0

    if not ignore_ambiguous:
        return float(np.sum(arr[:total] == arr[period:]) / total)

    vm = valid_base_mask(arr)
    valid = vm[:total] & vm[period:]
    valid_count = int(np.count_nonzero(valid))
    if valid_count <= 0 or valid_count < min_valid_frac * total:
        return 0.0
    matches = int(np.count_nonzero((arr[:total] == arr[period:]) & valid))
    return matches / valid_count


def windowed_match_counts(
    arr: np.ndarray, period: int, window: int,
    ignore_ambiguous: bool = True,
) -> Optional[np.ndarray]:
    """Sliding-window count of period-``period`` matches, in O(n) via a cumsum.

    Returns an array ``winsum`` where ``winsum[i]`` is the number of matching
    bases between ``arr[i:i+window]`` and ``arr[i+period:i+period+window]`` — i.e.
    the local autocorrelation strength at every window start. Its length is
    ``(arr.size - period) - window``.

    Returns ``None`` when the array is too short to contain a full window
    (``arr.size - period <= window``). Used by the catch-all and Tier-2
    autocorrelation seeders to find locally periodic stretches.

    With ``ignore_ambiguous`` (the default), a comparison involving a non-ACGT
    base contributes 0 rather than 1, so an all-``N`` window scores 0 instead of
    a full-match ``window``. Callers that turn a window count into a claimed
    interval should also require enough valid comparisons -- see
    :func:`windowed_valid_counts`.
    """
    counts = _windowed_counts(arr, period, window, ignore_ambiguous)
    return counts


def windowed_valid_counts(
    arr: np.ndarray, period: int, window: int
) -> Optional[np.ndarray]:
    """Sliding-window count of comparisons where both bases are A/C/G/T.

    Same shape and alignment as :func:`windowed_match_counts`, so
    ``valid[i] / window`` is the fraction of window ``i`` that carries usable
    sequence. Lets a caller reject a window whose apparent periodicity rests on
    too few real bases.
    """
    n = arr.size
    if n - period <= 0:
        return None
    vm = valid_base_mask(arr)
    return _cumsum_windows(vm[:n - period] & vm[period:n], window)


def _windowed_counts(arr: np.ndarray, period: int, window: int,
                     ignore_ambiguous: bool) -> Optional[np.ndarray]:
    n = arr.size
    if n - period <= 0:
        return None
    eq = (arr[:n - period] == arr[period:n])
    if ignore_ambiguous:
        vm = valid_base_mask(arr)
        eq &= vm[:n - period] & vm[period:n]
    return _cumsum_windows(eq, window)


def _cumsum_windows(eq: np.ndarray, window: int) -> Optional[np.ndarray]:
    cs = np.empty(eq.size + 1, dtype=np.int64)
    cs[0] = 0
    np.cumsum(eq, out=cs[1:])
    m = eq.size - window
    if m <= 0:
        return None
    return cs[window:window + m] - cs[:m]


def contiguous_true_runs(mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Start/end indices of maximal contiguous ``True`` runs in a boolean mask.

    Returns ``(run_starts, run_ends)`` with exclusive ends, so run ``k`` spans
    ``mask[run_starts[k]:run_ends[k]]``. Vectorized (no per-index Python) — the
    same ``np.diff`` over an ``int8`` view both seeders used inline. Empty arrays
    when there are no ``True`` positions.
    """
    d = np.diff(mask.view(np.int8))
    run_s = np.nonzero(d == 1)[0] + 1
    run_e = np.nonzero(d == -1)[0] + 1
    if mask[0]:
        run_s = np.concatenate(([0], run_s))
    if mask[-1]:
        run_e = np.concatenate((run_e, [mask.size]))
    return run_s, run_e
