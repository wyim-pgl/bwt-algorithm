#!/usr/bin/env python3
"""REV-27 -- CEN180 monomer recall stratified by monomer-to-consensus identity.

Section 2.2.3 states that monomers diverged below 80% identity from the 177 bp
consensus are outside the truth set, so the headline 99.7% is a recall over the
conserved subset. The reviewer asks whether that restriction inflates it.

The truth set itself cannot answer the question for monomers below 80%: the raw
blastn hit set was not retained, so there is nothing on disk at lower identity to
score against. What it can answer is whether recall degrades *within* the range
it covers. If recall is flat from 80% to 100% the restriction is not doing much
work; if it falls steeply toward the lower edge, the headline is optimistic and
the shape says by roughly how much.

Column 5 of colcen_cen180.bed carries each monomer's percent identity to the
consensus, so the stratification is a partition of the existing truth set and
introduces no new alignment.

Recall is computed exactly as score_colcen.py computes its unfiltered
`recall(pub)` column: a monomer counts as recovered if any call of any reported
period overlaps it by at least one base. Per-stratum monomer counts are printed
so a thin stratum cannot be read as a trend.

Usage: score_cen180_identity_strata.py [LABEL:PATH ...]
"""
import os
import sys
from collections import defaultdict

GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
HERE = os.path.dirname(os.path.abspath(__file__))
CHROMS = {"Chr1", "Chr2", "Chr3", "Chr4", "Chr5"}

# Left-closed, right-open except the last, which takes 100.0 itself.
STRATA = [(80, 85), (85, 90), (90, 95), (95, 99), (99, 100.01)]


def load_calls(path):
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] not in CHROMS:
                continue
            try:
                out[p[0]].append((int(p[1]), int(p[2])))
            except ValueError:
                continue
    return {c: merge(v) for c, v in out.items()}


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def load_monomers():
    """(chrom, start, end, identity) for every truth-set monomer."""
    rows = []
    with open(f"{GT}/colcen_cen180.bed") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[0] not in CHROMS:
                continue
            try:
                rows.append((p[0], int(p[1]), int(p[2]), float(p[4])))
            except ValueError:
                continue
    return rows


def recall(calls, monomers):
    """Fraction of the given monomers overlapped by at least one call."""
    if not monomers:
        return None
    by_chrom = defaultdict(list)
    for chrom, ms, me, _ in monomers:
        by_chrom[chrom].append((ms, me))
    hit = 0
    for chrom, mons in by_chrom.items():
        iv = calls.get(chrom, [])
        if not iv:
            continue
        mons.sort()
        j = 0
        for ms, me in mons:
            while j < len(iv) and iv[j][1] <= ms:
                j += 1
            if j < len(iv) and iv[j][0] < me:
                hit += 1
    return 100.0 * hit / len(monomers)


def main():
    monomers = load_monomers()
    ids = [m[3] for m in monomers]
    print(f"truth set: {len(monomers):,} monomers, identity "
          f"{min(ids):.1f}-{max(ids):.1f}%\n")

    cases = [
        ("BWTandem", f"{B}/bwtandem/bwt_colcen.bed"),
        ("TRF", f"{B}/trf/Col-CEN_v1.2_output.bed"),
        ("ULTRA (re-run)", f"{HERE}/../../../exp1_human/wp0/beds/ultra_colcen_p500.bed"),
        ("tantan (re-run)",
         "/data/gpfs/assoc/pgl/devel/exp1_human/wp0/fixcampaign/tantan/"
         "tantan_colcen_w500.bed"),
    ]
    for spec in sys.argv[1:]:
        label, _, path = spec.partition(":")
        if path:
            cases.append((label, path))

    buckets = [(lo, hi, [m for m in monomers if lo <= m[3] < hi])
               for lo, hi in STRATA]

    head = f"{'tool':18s}" + "".join(
        f"{f'{lo}-{hi if hi<=100 else 100}':>12s}" for lo, hi, _ in buckets) \
        + f"{'all':>12s}"
    print(head)
    print(f"{'n monomers':18s}"
          + "".join(f"{len(b):>12,d}" for _, _, b in buckets)
          + f"{len(monomers):>12,d}")
    print("-" * len(head))
    for label, path in cases:
        if not os.path.exists(path):
            print(f"{label:18s} {'(missing)':>12s}")
            continue
        calls = load_calls(path)
        cells = []
        for _, _, b in buckets:
            r = recall(calls, b)
            cells.append("n/a" if r is None else f"{r:.2f}%")
        allr = recall(calls, monomers)
        print(f"{label:18s}" + "".join(f"{c:>12s}" for c in cells)
              + f"{allr:>11.2f}%")


if __name__ == "__main__":
    main()
