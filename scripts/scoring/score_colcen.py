#!/usr/bin/env python3
"""WP0-B — rescore Table 2 (Arabidopsis Col-CEN) with the mreps re-run included.

Reimplements compute_metrics.py's three CEN180 metrics without bedtools:

  centromere coverage : fraction of each centromere interval covered by ANY tool
                        call, averaged over the 5 centromeres (no period filter)
  CEN180 count        : tool calls inside a centromere with period 150-200
  monomer recall      : fraction of ground-truth CEN180 monomers overlapped by a
                        tool call with period 150-200

The 150-200 window is the published definition and is period-assignment
sensitive: a tool that calls CEN180 as its ~356 bp dimer scores zero on it while
covering the array perfectly. ULTRA did exactly that on a synthetic 178 bp array,
so a second window (150-400, accepting the dimer) is reported alongside.

Validation first: the implementation must reproduce the published BWTandem row
(84.59 % coverage, 1,380 CEN180 calls, 99.73 % monomer recall) before any new
number from it is trusted.
"""
import os
import sys
from collections import defaultdict

GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
HERE = os.path.dirname(os.path.abspath(__file__))

CHROMS = {"Chr1", "Chr2", "Chr3", "Chr4", "Chr5"}


def load_bed(path, chroms=None):
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            if chroms and p[0] not in chroms:
                continue
            try:
                out[p[0]].append((int(p[1]), int(p[2]), p))
            except ValueError:
                continue
    return out


def period_of(parts):
    """Period from column 5, falling back to motif length.

    `int()` and not `int(float())`, which matters: filip's converted BEDs are 6
    columns with an integer period in column 5, but a BED written by BWTandem
    itself is 8 columns whose column 5 is a copy count like "46.3" and which has
    no period column at all — there the period is len(motif). int(float("46.3"))
    silently returns 46 and every period-banded metric then counts copy numbers
    instead of periods. That produced a CEN180 count of 21 for a BED whose
    coordinates are identical to the published one's, which is what exposed it.
    """
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            pass
    return len(parts[3]) if len(parts) >= 4 else None


def merge(iv):
    out = []
    for s, e in sorted(iv):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return out


def centromere_coverage(calls, cens):
    """Mean per-centromere covered fraction, using all calls (no period filter)."""
    fracs = []
    merged = {c: merge([(s, e) for s, e, _ in v]) for c, v in calls.items()}
    for chrom, cs, ce in cens:
        cov = 0
        for s, e in merged.get(chrom, []):
            if e <= cs:
                continue
            if s >= ce:
                break
            cov += min(e, ce) - max(s, cs)
        fracs.append(cov / (ce - cs) if ce > cs else 0.0)
    return 100.0 * sum(fracs) / len(fracs) if fracs else 0.0


def band_calls(calls, lo, hi):
    out = defaultdict(list)
    for chrom, v in calls.items():
        for s, e, parts in v:
            p = period_of(parts)
            if p is not None and lo <= p <= hi:
                out[chrom].append((s, e))
    return out


def count_in_centromeres(band, cens):
    n = 0
    for chrom, cs, ce in cens:
        for s, e in band.get(chrom, []):
            if s < ce and e > cs:
                n += 1
    return n


def cen180_bp_precision(calls, monomers):
    """intersect_bp(CEN180 GT, merged tool) / merged tool bp.

    Matches compute_bp_coverage()'s precision at compute_metrics.py:428: the tool
    side is merged, the GT side is NOT (each monomer contributes its own overlap),
    and no period filter is applied — the whole call set is used.
    """
    merged = {c: merge([(s, e) for s, e, _ in v]) for c, v in calls.items()}
    tool_bp = sum(e - s for v in merged.values() for s, e in v)
    if not tool_bp:
        return 0.0
    by_chrom = defaultdict(list)
    for chrom, ms, me, _ in monomers:
        by_chrom[chrom].append((ms, me))
    inter = 0
    for chrom, mons in by_chrom.items():
        iv = merged.get(chrom, [])
        if not iv:
            continue
        mons.sort()
        j = 0
        for ms, me in mons:
            while j < len(iv) and iv[j][1] <= ms:
                j += 1
            k = j
            while k < len(iv) and iv[k][0] < me:
                inter += min(iv[k][1], me) - max(iv[k][0], ms)
                k += 1
    return 100.0 * inter / tool_bp


def monomer_recall(band, monomers):
    """Fraction of GT monomers overlapped by at least one banded call."""
    hit = 0
    idx = {c: merge(v) for c, v in band.items()}
    for chrom, ms, me, _ in monomers:
        for s, e in idx.get(chrom, []):
            if e <= ms:
                continue
            if s >= me:
                break
            hit += 1
            break
    return 100.0 * hit / len(monomers) if monomers else 0.0


def load_flat(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3 and p[0] in CHROMS:
                try:
                    rows.append((p[0], int(p[1]), int(p[2]), p))
                except ValueError:
                    pass
    return rows


def main():
    cens = [(c, s, e) for c, s, e, _ in load_flat(f"{GT}/colcen_centromeres.bed")]
    monomers = load_flat(f"{GT}/colcen_cen180.bed")
    print(f"GT: {len(cens)} centromeres, {len(monomers):,} CEN180 monomers\n")

    cases = [
        ("BWTandem", f"{B}/bwtandem/bwt_colcen.bed", "84.59 / 1380 / 57.17 / 99.73"),
        ("TRF", f"{B}/trf/Col-CEN_v1.2_output.bed", "84.39 /  191 / 82.79 / 97.55"),
        ("tantan", f"{B}/tantan/Col-CEN_v1.2_output.bed", " 1.04 /    0 /  0.23 /  0.35"),
        ("mreps (published, maxp 6)", f"{B}/mreps/Col-CEN_v1.2_output.bed", " 1.58 /    0 /  0.56 /  5.76"),
        ("mreps (RE-RUN 150-400)", f"{HERE}/beds/mreps_colcen_150_400.bed", "  -- new --"),
        ("ULTRA (published, -p 100)", f"{HERE}/beds/ultra_colcen.bed", " 2.25 /    0 /  0.21 /  0.62"),
    ]
    ultra_new = f"{HERE}/beds/ultra_colcen_p500.bed"
    if os.path.exists(ultra_new):
        cases.append(("ULTRA (RE-RUN -p 500)", ultra_new, "  -- new --"))

    # Extra BEDs from the command line, as LABEL:PATH. Used for the gap-fill
    # ablation, whose whole point is to be scored by exactly the code that
    # produced the published Col-CEN row rather than by a second implementation.
    for spec in sys.argv[1:]:
        label, _, path = spec.partition(":")
        if path:
            cases.append((label, path, "  -- ablation --"))

    # The published row is internally inconsistent, and reproducing it requires
    # matching that: compute_cen180_count() filters to period 150-200, but
    # compute_cen180_monomer_recall() is called with the FULL tool bed at
    # compute_metrics.py:429 — extract_cen180() (the 150-200 filter) is defined
    # but never used at that call site. So "count" and "recall" in the same row
    # are computed over different call sets.
    print(f"{'tool':28s} {'regions':>9s} {'cov%':>7s} {'n150-200':>9s} {'bpPrec%':>8s} | "
          f"{'recall(pub)':>11s} {'rec150-200':>11s} {'rec150-400':>11s}   "
          f"published (cov/count/bpPrec/recall)")
    for label, path, pub in cases:
        if not os.path.exists(path):
            print(f"{label:28s} {'(missing)':>9s}")
            continue
        calls = load_bed(path, CHROMS)
        n_reg = sum(len(v) for v in calls.values())
        cov = centromere_coverage(calls, cens)
        b1 = band_calls(calls, 150, 200)
        b2 = band_calls(calls, 150, 400)
        allc = {c: [(s, e) for s, e, _ in v] for c, v in calls.items()}
        print(f"{label:28s} {n_reg:>9,d} {cov:>6.2f}% "
              f"{count_in_centromeres(b1, cens):>9,d} "
              f"{cen180_bp_precision(calls, monomers):>7.2f}% | "
              f"{monomer_recall(allc, monomers):>10.2f}% "
              f"{monomer_recall(b1, monomers):>10.2f}% "
              f"{monomer_recall(b2, monomers):>10.2f}%   {pub}")


if __name__ == "__main__":
    main()
