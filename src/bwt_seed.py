"""Shared BWT k-mer seeding core for tandem repeat detection.

Used by both Tier 2 and Tier 3.  The algorithm is:
  1. Sample k-mers from the text at a configurable stride.
  2. Use FM-index backward_search / locate_positions to find all occurrences.
  3. Detect periodic runs in the position arrays (arithmetic progressions).
  4. Extend seed positions with mismatch tolerance.
  5. Return raw candidate regions for tier-specific post-processing.
"""
from __future__ import annotations  # Allow string-form type hints for Python 3.9 and below

import os                           # Environment-variable overrides
from dataclasses import dataclass   # Decorator for defining immutable data containers
from typing import List, Optional, Set, Tuple  # Generic types for type hints

import numpy as np  # NumPy for array operations and numerical computation

from .accelerators import extend_with_mismatches, find_periodic_runs
from .bwt_core import BWTCore, effective_length  # FM-index core module (BWT, backward_search, locate_positions, etc.)


class KTupleTable:
    """Sorted k-tuple table: the non-index way to answer "where does this k-mer occur".

    Exists for the seeding ablation (BWT_SEED_KTUPLE), which asks what the
    FM-index buys the seeding step. It answers the same question as
    backward_search plus a suffix-array slice and must return the same position
    set, so the only thing that differs between the two arms is cost.

    Every k-mer of the text is 2-bit encoded into one uint64 (k <= 32), the codes
    are sorted once, and a lookup is a pair of binary searches over that sorted
    array. Positions containing a base other than A, C, G or T are excluded, as
    the FM-index path excludes them at its own k-mer filter.
    """

    _CODE = np.full(256, 255, dtype=np.uint8)
    for _b, _v in ((b"A", 0), (b"C", 1), (b"G", 2), (b"T", 3)):
        _CODE[_b[0]] = _v
        _CODE[_b.lower()[0]] = _v

    def __init__(self, text_arr: np.ndarray, k: int, n: int):
        if k > 32:
            raise ValueError(f"KTupleTable supports k <= 32, got {k}")
        self.k = k
        codes = self._CODE[text_arr[:n]]
        valid = codes != 255
        # A k-mer is usable only if all k of its bases are valid; a rolling
        # minimum over the validity flags is the cheapest way to say that.
        ok = valid[: n - k + 1].copy()
        for off in range(1, k):
            ok &= valid[off : n - k + 1 + off]
        packed = np.zeros(n - k + 1, dtype=np.uint64)
        safe = np.where(valid, codes, 0).astype(np.uint64)
        for off in range(k):
            packed <<= np.uint64(2)
            packed |= safe[off : n - k + 1 + off]
        self._pos = np.flatnonzero(ok).astype(np.int64)
        self._key = packed[self._pos]
        order = np.argsort(self._key, kind="stable")
        self._key = self._key[order]
        self._pos = self._pos[order]

    def encode(self, kmer: str) -> Optional[int]:
        code = 0
        for ch in kmer:
            v = self._CODE[ord(ch)]
            if v == 255:
                return None
            code = (code << 2) | int(v)
        return code

    def positions(self, kmer: str) -> Optional[np.ndarray]:
        """Every start position of `kmer`, or None if it contains a non-ACGT base."""
        code = self.encode(kmer)
        if code is None:
            return None
        lo = int(np.searchsorted(self._key, np.uint64(code), side="left"))
        hi = int(np.searchsorted(self._key, np.uint64(code), side="right"))
        return self._pos[lo:hi]


_SEED_KTUPLE = os.environ.get("BWT_SEED_KTUPLE", "0") == "1"


def _ktuple_table(bwt: BWTCore, k: int, n: int) -> KTupleTable:
    """Build the ablation's table once per (index, k) and hang it off the index."""
    cache = getattr(bwt, "_ktuple_cache", None)
    if cache is None:
        cache = {}
        setattr(bwt, "_ktuple_cache", cache)
    if k not in cache:
        cache[k] = KTupleTable(bwt.text_arr, k, n)
    return cache[k]


@dataclass
class SeedCandidate:
    """A raw repeat candidate from BWT seeding."""
    start: int      # Start position of the repeat array after extension (0-based)
    end: int        # End position of the repeat array after extension (0-based)
    period: int     # Detected repeat unit length (bp)
    copies: int     # Detected number of copies
    motif: str      # Representative motif string (extracted from full_start position)
    seed_pos: int   # Original seed position that generated this candidate (value of i)


def bwt_kmer_seed_scan(
    bwt: BWTCore,                        # FM-index object
    min_period: int,                     # Minimum repeat unit length to detect
    max_period: int,                     # Maximum repeat unit length to detect
    kmer_size: int = 16,                 # Length of k-mers to sample (default 16 bp)
    stride: int = 50,                    # K-mer sampling interval (default 50 bp)
    min_copies: int = 2,                 # Minimum number of copies to qualify as a repeat
    allowed_mismatch_rate: float = 0.20, # Allowed mismatch rate during extension (0.0-0.5)
    tolerance_ratio: float = 0.03,       # Tolerance ratio for period jitter in periodic run detection (default 3%)
    max_occurrences: int = 5000,         # Maximum k-mer occurrence count (skip if exceeded, likely low-complexity)
    covered_mask: Optional[np.ndarray] = None,  # Boolean mask of positions already found by previous tiers
    show_progress: bool = False,         # Whether to print progress
    label: str = "",                     # Label string for progress messages
) -> List[SeedCandidate]:
    """BWT k-mer seeding scan for tandem repeat candidates.

    Samples k-mers from the text at `stride` intervals, uses FM-index to
    find all occurrences, detects periodic runs in the occurrence positions,
    and extends seed hits with mismatch tolerance.

    Parameters
    ----------
    bwt : BWTCore
        The BWT/FM-index built on the sequence.
    min_period, max_period : int
        Period range to search for.
    kmer_size : int
        Length of k-mers to sample (default 16).  Should be shorter than
        min_period so that each repeat copy contains at least one full k-mer.
    stride : int
        Distance between consecutive k-mer samples (default 50).
    min_copies : int
        Minimum number of tandem copies required.
    allowed_mismatch_rate : float
        Mismatch tolerance for extension (0.0-0.5).
    tolerance_ratio : float
        Tolerance for period jitter in periodic run detection (default 3%).
    max_occurrences : int
        Skip k-mers with more occurrences than this (likely low-complexity).
    covered_mask : np.ndarray or None
        Boolean mask of positions already found by a previous tier.
        Positions where covered_mask[i] is True are skipped as seed origins.
    show_progress : bool
        Whether to print progress.
    label : str
        Label prefix for progress messages.

    Returns
    -------
    list[SeedCandidate]
        Raw repeat candidates.  Callers perform tier-specific post-processing
        (primitive period reduction, HOR detection, DP refinement, etc.).
    """
    text_arr = bwt.text_arr  # Original sequence used for BWT (numpy uint8 array)
    n = effective_length(text_arr)

    # Return immediately if sequence is too short to contain repeats meeting minimum copy count
    if n < min_period * min_copies:
        return []

    # Clamp k-mer size if larger than min_period (k-mer must fit within one copy)
    effective_kmer = min(kmer_size, min_period)  # Limit to not exceed repeat unit length
    if effective_kmer < 6:
        effective_kmer = min(6, min_period)  # Minimum 6 bp: shorter k-mers cause explosive FM-index hits

    candidates: List[SeedCandidate] = []               # List of detected repeat candidates
    seen_regions: Set[Tuple[int, int]] = set()         # Set of (start//bucket, period) keys for duplicate removal
    seen_kmers: Set[str] = set()                       # Set of already-queried k-mers (avoid re-querying)
    bwt_queries = 0  # FM-index query counter (for progress reporting)

    i = 0  # Current sampling position (incremented by stride)
    while i <= n - effective_kmer:  # Iterate until the final full-length k-mer start
        # Skip positions already found by previous tiers
        if covered_mask is not None and covered_mask[i]:
            i += stride  # Position already covered; move to next sample position
            continue

        # Extract k-mer at current position
        kmer_arr = text_arr[i:i + effective_kmer]  # Sequence slice for the k-mer
        kmer_str = kmer_arr.tobytes().decode('ascii', errors='replace').upper()  # Normalize input case

        # Skip if k-mer was already queried or contains non-DNA bases
        if kmer_str in seen_kmers:  # Check cache to avoid duplicate queries
            i += stride
            continue
        if not all(c in 'ACGT' for c in kmer_str):  # Skip if contains N, lowercase, or other non-DNA characters
            i += stride
            continue
        seen_kmers.add(kmer_str)  # Register queried k-mer in cache

        # --- Query all occurrence positions of the k-mer ---
        # Default: FM-index backward search. Under the REV-2 seeding ablation
        # (BWT_SEED_KTUPLE=1) the same occurrence set comes from a sorted k-tuple
        # table instead, so the two arms differ in how seeds are located and in
        # nothing else. Both paths must yield the same positions; the ablation is
        # a cost measurement, not an accuracy one.
        bwt_queries += 1  # Increment occurrence-query counter
        if _SEED_KTUPLE:
            table = _ktuple_table(bwt, effective_kmer, n)
            hits = table.positions(kmer_str)
            if hits is None:  # k-mer holds a non-ACGT base; the FM path skips these too
                i += stride
                continue
            occ_count = int(hits.size)
        else:
            sp, ep = bwt.backward_search(kmer_str)  # Returns suffix array range [sp, ep]
            if sp == -1:  # Not found (should not happen in theory, but defensive check)
                i += stride
                continue
            occ_count = ep - sp + 1  # Calculate total occurrence count of the k-mer

        if occ_count < min_copies or occ_count > max_occurrences:
            # Skip if too few occurrences (no repeat) or too many (low-complexity sequence)
            i += stride
            continue

        if _SEED_KTUPLE:
            positions = hits.tolist()
        else:
            positions = bwt.suffix_array[sp:ep + 1].tolist()  # Reuse the SA interval from backward_search
        if len(positions) < min_copies:  # Skip if actual position count is below minimum copies
            i += stride
            continue

        # --- Detect arithmetic progressions (periodic runs) in occurrence position array ---
        pos_arr = np.array(sorted(positions), dtype=np.int64)
        # (Sorting is required for accurate arithmetic progression detection)

        # find_periodic_runs: detect evenly-spaced groups within tolerance_ratio in position array
        patterns = find_periodic_runs(
            pos_arr, min_period, max_period, min_copies, tolerance_ratio
        )  # Returns: [(run_start, run_end, period), ...] -- list of periodic k-mer occurrence runs

        for run_start, run_end, period in patterns:  # Process each periodic run
            run_start = int(run_start)  # Convert numpy type to Python int
            run_end = int(run_end)      # Convert numpy type to Python int
            period = int(period)        # Convert numpy type to Python int

            # Remove duplicate candidates for the same region: generate key from start position and period
            region_key = (run_start // max(period, 1), period)  # Key to identify the repeat region
            if region_key in seen_regions:  # Skip if this region was already processed
                continue

            # Skip if overlap with already covered region exceeds 50% (prevent redundant detection)
            if covered_mask is not None:
                covered_count = int(np.sum(
                    covered_mask[run_start:min(run_end + period, n)]
                ))  # Count already covered positions within this span
                span = run_end + period - run_start  # Total length of the run (including last k-mer)
                if span > 0 and covered_count > span * 0.5:
                    continue  # Skip if more than 50% is already covered

            # --- Mismatch-tolerant extension: extend the actual boundaries of the repeat array in both directions ---
            res = extend_with_mismatches(
                text_arr, run_start, period, n, allowed_mismatch_rate
            )  # Returns: (arr_start, arr_end, copies, full_start, full_end) or None

            if res is not None:
                arr_start, arr_end, copies, full_start, full_end = res
                # arr_start/arr_end: alignment-based boundaries, full_start/full_end: actual extended boundaries
            else:
                # Fallback when extension fails: use raw boundaries of the k-mer run
                full_start = run_start
                # run_end is the start position of the last k-mer, so add period to get end position
                full_end = run_end + period
                copies = max(1, (full_end - full_start) // period)  # Estimate copy count

            if copies < min_copies:  # Do not register as candidate if below minimum copy count
                continue

            # Extract motif string from the confirmed repeat region
            motif_start = max(0, full_start)  # Clamp to prevent negative index
            motif_arr = text_arr[motif_start:motif_start + period]  # Sequence slice for one copy
            motif_str = motif_arr.tobytes().decode('ascii', errors='replace')  # Convert byte array to string

            # Create SeedCandidate object and add to candidate list
            candidates.append(SeedCandidate(
                start=full_start,   # Extended repeat array start position
                end=full_end,       # Extended repeat array end position
                period=period,      # Detected repeat unit length
                copies=copies,      # Detected copy count
                motif=motif_str,    # Representative motif string
                seed_pos=i,         # Original seed position that generated this candidate
            ))
            seen_regions.add(region_key)  # Mark this region as processed

            # Mark discovered region in the mask: prevent re-seeding in the same span later
            if covered_mask is not None:
                covered_mask[full_start:min(full_end, n)] = True  # Mark this span as covered

        i += stride  # Move to next sample position

    # Print progress (only when show_progress is True)
    if show_progress:
        print(
            f"  [{label}] BWT k-mer seeding: {len(candidates)} candidates from "
            f"{bwt_queries} FM-index queries (kmer={effective_kmer}, stride={stride})",
            flush=True,
        )  # Print candidate count, FM-index query count, effective k-mer size, and sampling stride

    return candidates  # Return list of all detected raw repeat candidates
