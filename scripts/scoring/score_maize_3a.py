#!/usr/bin/env python3
"""WP0-A + WP0-C — regenerate maize Table 3A under one consistent rule.

Two defects meet in this table.

**WP0-A (scoring rule).** `compute_metrics.py:457` reads:

    eval_bed = bed
    if tool == 'bwtandem':
        os.system(f"awk '$5 <= 6' {bed} > {eval_bed}")

so the period band is applied **to BWTandem alone**; every other tool is scored on
its BED as given. tantan's exp3A file is not even microsatellite-specific — its
3A/3B/3C BEDs are one unfiltered file under three names (identical 68,279,539 bytes),
so its "SSR bp" is every repeat it found at any period.

**WP0-C (TAG definition).** `compute_tag_metrics()` counts a call when
`'TAG' in motif or 'CTA' in motif`. For BWTandem/ULTRA/tantan the motif column holds
one primitive period, but for TRF/TRASH/NCRF it holds a consensus up to tens of kb,
where those 3-mers appear by chance. The column therefore measures motif length more
than TAG content.

This prints, for every tool:
  * the published rule, to verify the table reproduces
  * band 1-6 applied to everyone (the rule the table's title claims)
  * TAG under both the substring definition and a fair one (period 3 and the
    leading 3-mer a rotation of TAG)

Usage: score_maize_3a.py
"""
import os

B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
HERE = os.path.dirname(os.path.abspath(__file__))
MZ = "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic"

# (label, path, applies-band-under-published-rule, published "SSR bp / regions / TAG / longest")
TOOLS = [
    ("TRF",       f"{B}/trf/{MZ}_trf_exp3A_microsatellite.bed",      False,
     "4,793,515 / 29,559 / 5,263 / 562.79"),
    ("ULTRA",     f"{HERE}/beds/ultra_maize_exp3A.bed",              False,
     "0 / 0 / 0 / 0.00  (recovered from a 0-byte BED, WP0.1)"),
    ("BWTandem",  os.environ.get("BWT_MAIZE_BED", f"{B}/bwtandem/bwt_maize.bed"), True,
     "41,007,584 / 1,498,365 / 39,299 / 175.60"),
    ("tantan",    f"{B}/tantan/{MZ}_tantan_exp3A_microsatellite.bed", False,
     "54,933,424 / 1,303,967 / 181,205 / 210.09"),
    ("NCRF",      f"{B}/ncrf/{MZ}_ncrf_TAG_exp3A_microsatellite.bed", False,
     "2,604,968 / 214 / 214 / 235.40"),
    ("TRASH-dn",  f"{B}/trash/{MZ}_trash_denovo_exp3A_microsatellite.bed", False,
     "47,519,576 / 3,034 / 1,893 / 490.00"),
    ("TRASH-tpl", f"{B}/trash/{MZ}_trash_templates_TAG_exp3A_microsatellite.bed", False,
     "47,525,434 / 3,016 / 1,884 / 490.00"),
]

# The six rotations/complements of TAG. A period-3 array of TAG can be written
# starting at any of its three phases, on either strand.
TAG_ROT = {"TAG", "AGT", "GTA", "CTA", "TAC", "ACT"}


def period_of(parts):
    """Period from column 5, falling back to motif length.

    `int()` and not `int(float())`: filip's converted BEDs are 6 columns with an
    integer period in column 5, while a BED written by BWTandem itself is 8
    columns whose column 5 is a copy count like "46.3" and which carries no
    period column — there the period is len(motif). int(float("46.3")) silently
    returns 46, and every period-banded metric then counts copy numbers instead
    of periods.
    """
    if len(parts) >= 5:
        try:
            return int(parts[4])
        except ValueError:
            pass
    return len(parts[3]) if len(parts) >= 4 else None


def merged_bp(by_chrom):
    """Total bp after merging overlapping intervals, per chromosome."""
    total = 0
    for iv in by_chrom.values():
        iv.sort()
        cs, ce = iv[0]
        for s, e in iv[1:]:
            if s <= ce:
                ce = max(ce, e)
            else:
                total += ce - cs
                cs, ce = s, e
        total += ce - cs
    return total


def scan(path, band=None):
    """SSR bp/regions under the optional band, plus both TAG counts.

    `regions` is the raw line count but `bp` is computed AFTER merging overlaps —
    that is what compute_metrics.py:total_regions_and_bp does (it counts lines,
    then pipes the file through `bedtools sort | bedtools merge` for the bp). A
    non-merged sum overstates every tool that emits nested or overlapping calls,
    by 0.7 % for TRF up to 1.7 % for tantan on this genome.
    """
    from collections import defaultdict
    n = 0
    by_chrom = defaultdict(list)
    tag_pub = tag_fair = 0
    long_pub = long_fair = 0
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            span = end - start
            per = period_of(p)
            if band is not None:
                if per is None or not (band[0] <= per <= band[1]):
                    continue
            n += 1
            by_chrom[p[0]].append((start, end))
            motif = (p[3] if len(p) > 3 else "").upper()
            if "TAG" in motif or "CTA" in motif:
                tag_pub += 1
                long_pub = max(long_pub, span)
            if per == 3 and motif[:3] in TAG_ROT:
                tag_fair += 1
                long_fair = max(long_fair, span)
    return n, merged_bp(by_chrom), tag_pub, long_pub, tag_fair, long_fair


def main():
    print("########## PUBLISHED RULE (band on BWTandem only) ##########")
    print(f"{'tool':10s} {'SSR bp':>14s} {'regions':>12s} {'TAG(sub)':>10s} "
          f"{'longest kb':>11s}   published (bp / regions / TAG / longest)")
    for label, path, banded, pub in TOOLS:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"{label:10s} {'(missing/empty)':>14s}   {pub}")
            continue
        n, bp, tp, lp, _, _ = scan(path, (1, 6) if banded else None)
        print(f"{label:10s} {bp:>14,d} {n:>12,d} {tp:>10,d} {lp/1000:>11.2f}   {pub}")

    print("\n########## CONSISTENT RULE — band 1-6 for every tool ##########")
    print(f"{'tool':10s} {'SSR bp':>14s} {'regions':>12s} {'TAG(sub)':>10s} "
          f"{'TAG(fair)':>10s} {'longest kb':>11s}")
    for label, path, _, _ in TOOLS:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"{label:10s} {'(missing/empty)':>14s}")
            continue
        n, bp, tp, lp, tf, lf = scan(path, (1, 6))
        print(f"{label:10s} {bp:>14,d} {n:>12,d} {tp:>10,d} {tf:>10,d} "
              f"{lf/1000:>11.2f}")

    print("\n########## TAG ARRAYS — substring vs fair definition (no band) ##########")
    print("substring = compute_tag_metrics(): 'TAG' in motif or 'CTA' in motif")
    print("fair      = period 3 and the leading 3-mer in {TAG,AGT,GTA,CTA,TAC,ACT}\n")
    print(f"{'tool':10s} {'TAG(sub)':>12s} {'TAG(fair)':>12s} {'ratio':>8s} "
          f"{'longest(fair) kb':>18s}")
    for label, path, _, _ in TOOLS:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            print(f"{label:10s} {'(missing/empty)':>12s}")
            continue
        _, _, tp, _, tf, lf = scan(path, None)
        ratio = f"{tp/tf:.1f}x" if tf else "n/a"
        print(f"{label:10s} {tp:>12,d} {tf:>12,d} {ratio:>8s} {lf/1000:>18.2f}")


if __name__ == "__main__":
    main()
