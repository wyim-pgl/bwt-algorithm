"""BED column 6 must be an integer for every pass that can emit a record.

`TandemRepeat.tier` is declared `int` and columns 1-3 hold 1/2/3, but the two
post-tier passes used to pass `tier="satellite"` / `tier="catchall"`. The
satellite pass runs by default, so `int(fields[5])` blew up on ordinary
centromeric output; the same value also goes out as VCF `TIER=`.
"""
import numpy as np
import pytest

from src.finder import TandemRepeatFinder
from src.models import TandemRepeat, TIER_SATELLITE, TIER_CATCHALL

RNG = np.random.RandomState(20260805)
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


def bed_tier_field(repeat: TandemRepeat) -> str:
    return repeat.to_bed().split("\t")[5]


def test_satellite_pass_emits_an_integer_tier():
    monomer = _rand_seq(171)
    clean = _diverged_array(monomer, 30, 0.02)
    gap = _diverged_array(monomer, 80, 0.42)
    seq = clean + gap + clean
    L1, n = len(clean), len(clean) * 2 + len(gap)
    L2 = L1 + len(gap)

    finder = _finder(seq)
    anchors = [
        TandemRepeat(chrom="syn", start=0, end=L1, motif=seq[:600],
                     copies=L1 / 600, length=L1, tier=2,
                     consensus_motif=seq[:600], mismatch_rate=0.02),
        TandemRepeat(chrom="syn", start=L2, end=n, motif=seq[L2:L2 + 600],
                     copies=(n - L2) / 600, length=n - L2, tier=2,
                     consensus_motif=seq[L2:L2 + 600], mismatch_rate=0.02),
    ]
    filled = finder._fill_satellite_gaps(anchors)
    new = [r for r in filled if r not in anchors]
    assert new, "satellite pass emitted nothing, so this proves nothing"

    for r in new:
        assert r.tier == TIER_SATELLITE
        assert int(bed_tier_field(r)) == TIER_SATELLITE


def test_catchall_pass_emits_an_integer_tier(monkeypatch):
    seq = _diverged_array("ACGTT", 60, 0.10) + _rand_seq(200)
    finder = _finder(seq)
    monkeypatch.setenv("CATCHALL_MIN_IDENTITY", "0.72")

    filled = finder._catchall_periodicity_fill([])
    assert filled, "catch-all pass emitted nothing, so this proves nothing"

    for r in filled:
        assert r.tier == TIER_CATCHALL
        assert int(bed_tier_field(r)) == TIER_CATCHALL


def test_every_pass_id_is_distinct():
    """A shared id would make column 6 ambiguous rather than merely non-numeric."""
    ids = [1, 2, 3, TIER_SATELLITE, TIER_CATCHALL]
    assert len(set(ids)) == len(ids)
    assert all(isinstance(i, int) for i in ids)


@pytest.mark.parametrize("tier", [1, 2, 3, TIER_SATELLITE, TIER_CATCHALL])
def test_bed_and_vcf_tier_fields_parse_as_int(tier):
    r = TandemRepeat(chrom="chr1", start=0, end=10, motif="AT", copies=5.0,
                     length=10, tier=tier)
    assert int(r.to_bed().split("\t")[5]) == tier
    info = dict(kv.split("=", 1) for kv in r.to_vcf_info().split(";"))
    assert int(info["TIER"]) == tier
