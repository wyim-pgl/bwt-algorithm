"""What `_fill_satellite_gaps` emits must be a valid interval over real bases.

Two defects, both from the pass treating `bwt.text_arr` as the sequence when it
actually carries a trailing `$` sentinel and is one byte longer:

  * a gap running to the end of the sequence was reported with
    `end == len(sequence) + 1` — an interval past the end of the chromosome,
    which bedtools rejects. On the published human BED, 19 of the 11,803
    satellite intervals end 1 bp past their chromosome;
  * the motif is a raw window off the text, and nothing checked it was DNA, so
    `N` (an assembly gap is trivially periodic, N == N) or `$` could be reported
    as the repeat unit. 198 satellite rows on that same BED have non-ACGT motifs.

`_catchall_periodicity_fill` already had both guards; this pins them on the
satellite pass, which unlike the catch-all is on by default.
"""
import numpy as np

from src.bwt_core import effective_length
from src.finder import TandemRepeatFinder
from src.models import TandemRepeat

RNG = np.random.RandomState(770509)
BASES = np.frombuffer(b"ACGT", dtype=np.uint8)


def _rand_seq(n):
    return BASES[RNG.randint(0, 4, n)].tobytes().decode()


def _diverged_array(monomer, n_copies, p_mut):
    m = np.frombuffer(monomer.encode(), dtype=np.uint8).copy()
    out = []
    for _ in range(n_copies):
        c = m.copy()
        mask = RNG.rand(len(c)) < p_mut
        c[mask] = BASES[RNG.randint(0, 4, int(mask.sum()))]
        out.append(c)
    return np.concatenate(out).tobytes().decode()


def _finder(seq):
    return TandemRepeatFinder(seq, chromosome="syn", min_period=1,
                              max_period=8000, enabled_tiers={"tier1"})


def _anchor(seq, start, end):
    """A large-HOR-motif call, which is what marks a block as near-satellite."""
    motif = seq[start:start + 600]
    return TandemRepeat(chrom="syn", start=start, end=end, motif=motif,
                        copies=(end - start) / 600, length=end - start, tier=2,
                        consensus_motif=motif, mismatch_rate=0.02)


def test_text_arr_really_is_one_longer_than_the_sequence():
    """The premise: if this stops holding, the two tests below prove nothing."""
    seq = _rand_seq(2000)
    finder = _finder(seq)
    assert len(finder.bwt.text_arr) == len(seq) + 1
    assert effective_length(finder.bwt.text_arr) == len(seq)


def test_terminal_gap_does_not_run_past_the_end_of_the_sequence():
    monomer = _rand_seq(171)
    covered_head = _diverged_array(monomer, 30, 0.02)
    uncovered_tail = _diverged_array(monomer, 40, 0.10)
    seq = covered_head + uncovered_tail
    n = len(seq)

    finder = _finder(seq)
    anchors = [_anchor(seq, 0, len(covered_head))]

    filled = finder._fill_satellite_gaps(anchors)
    new = [r for r in filled if r not in anchors]
    assert new, "no gap was filled, so the bound is untested"
    assert any(r.end > len(covered_head) + 3000 for r in new), \
        "the fill did not reach the tail, so the terminal bound is untested"

    for r in new:
        assert r.end <= n, f"interval {r.start}-{r.end} ends past the sequence ({n} bp)"
        assert r.start >= 0


def test_assembly_gap_is_not_reported_as_a_repeat_motif():
    monomer = _rand_seq(171)
    left = _diverged_array(monomer, 30, 0.02)
    n_run = "N" * 4000                      # trivially periodic: N == N at every offset
    right = _diverged_array(monomer, 30, 0.02)
    seq = left + n_run + right
    lo, hi = len(left), len(left) + len(n_run)

    finder = _finder(seq)
    anchors = [_anchor(seq, 0, lo), _anchor(seq, hi, len(seq))]

    filled = finder._fill_satellite_gaps(anchors)
    new = [r for r in filled if r not in anchors]

    for r in new:
        motif = r.consensus_motif or r.motif
        assert motif and all(c in "ACGT" for c in motif), \
            f"non-DNA motif emitted: {motif[:40]!r} at {r.start}-{r.end}"


def _overlap(r, lo, hi):
    return max(0, min(r.end, hi) - max(r.start, lo))


def test_assembly_gap_is_not_claimed_as_satellite():
    """The motif guard above is necessary but not sufficient.

    A window at the edge of the gap carries real satellite, so the *motif* comes
    out ACGT and passes the check above -- while the emitted interval still
    swallows the N run behind it. Nothing may claim the gap itself.
    """
    monomer = _rand_seq(171)
    left = _diverged_array(monomer, 30, 0.02)
    n_run = "N" * 4000
    right = _diverged_array(monomer, 30, 0.02)
    seq = left + n_run + right
    lo, hi = len(left), len(left) + len(n_run)

    finder = _finder(seq)
    anchors = [_anchor(seq, 0, lo), _anchor(seq, hi, len(seq))]
    try:
        filled = finder._fill_satellite_gaps(anchors)
    finally:
        finder.cleanup()

    for r in [x for x in filled if x not in anchors]:
        covered = _overlap(r, lo, hi)
        assert covered == 0, (
            f"call {r.start}-{r.end} claims {covered} bp of the {hi - lo} bp "
            "assembly gap"
        )


def test_mostly_n_block_with_a_real_edge_is_refused():
    """A short real seed must not authorize a call over a mostly-N block.

    This is the shape seen on Arabidopsis Chr4:3,955,856-3,957,243, which was
    emitted as a 102 bp satellite with 13.6 copies while being 72% N.
    """
    monomer = _rand_seq(171)
    flank = _diverged_array(monomer, 30, 0.02)
    seq = flank + _rand_seq(200) + "N" * 3000 + _rand_seq(200) + flank
    lo = len(flank) + 200
    hi = lo + 3000

    finder = _finder(seq)
    anchors = [_anchor(seq, 0, len(flank)),
               _anchor(seq, len(seq) - len(flank), len(seq))]
    try:
        filled = finder._fill_satellite_gaps(anchors)
    finally:
        finder.cleanup()

    for r in [x for x in filled if x not in anchors]:
        assert _overlap(r, lo, hi) == 0, \
            f"call {r.start}-{r.end} rests on the {hi - lo} bp N block"


def test_two_arrays_are_not_merged_across_an_assembly_gap():
    """The merge path scores the gap by autocorrelation; N == N used to pass."""
    monomer = _rand_seq(171)
    left = _diverged_array(monomer, 20, 0.02)
    n_run = "N" * 2000
    right = _diverged_array(monomer, 20, 0.02)
    seq = left + n_run + right
    lo, hi = len(left), len(left) + len(n_run)

    finder = _finder(seq)
    motif = seq[:171]
    calls = [
        TandemRepeat(chrom="syn", start=0, end=lo, motif=motif,
                     copies=20.0, length=lo, tier=2,
                     consensus_motif=motif, mismatch_rate=0.02),
        TandemRepeat(chrom="syn", start=hi, end=len(seq), motif=motif,
                     copies=20.0, length=len(seq) - hi, tier=2,
                     consensus_motif=motif, mismatch_rate=0.02),
    ]
    try:
        merged = finder._merge_adjacent_repeats(calls)
    finally:
        finder.cleanup()

    for r in merged:
        assert not (r.start < lo and r.end > hi), \
            f"call {r.start}-{r.end} was merged across the {hi - lo} bp gap"
