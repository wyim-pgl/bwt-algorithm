#!/usr/bin/env python3
"""REV-36 -- catalog-blind family recovery on maize.

Section 4.1 claims BWTandem's output contains repeat families beyond the
annotated set and concedes that "how large that set is relative to the annotated
one is something we assert rather than measure". This measures the half that can
be measured: given only the tool's own output and no template, do the curated
families come back out, and how many comparably large clusters are not curated?

Nothing here consults the curated coordinates until after clustering. The
pipeline is: take the reported motifs, reduce each to a strand- and
rotation-invariant canonical form, group by that form, rank groups by the total
sequence they span, and only then ask which curated family each of the top
groups sits on.

Canonicalisation is the same idea as the pipeline's own: the lexicographically
smallest rotation of the motif or its reverse complement. Two calls of the same
family with a substitution still land in different groups, so this
under-clusters, and the family totals below are lower bounds. That is the
conservative direction for the question being asked.

A group is called curated if at least half its spanned bases fall inside the
curated arrays of one family; otherwise it is uncurated. "Spurious" is not
claimed for those: the maize annotation is not exhaustive, so an uncurated
cluster may be a real family, a false positive, or an artefact of the
under-clustering above, and the three are not separated here.

Usage: catalog_blind_families.py [BED] [TOPN]
"""
import os
import sys
from collections import defaultdict

GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
B = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"

FAMILIES = {
    "knob180": f"{GT}/mo17_knob180_arrays.bed",
    "TR-1": f"{GT}/mo17_tr1_arrays.bed",
    "CentC": f"{GT}/mo17_centc_arrays.bed",
}

COMP = str.maketrans("ACGTacgt", "TGCATGCA")


def revcomp(s):
    return s.translate(COMP)[::-1]


def least_rotation(s):
    """Booth's algorithm: index of the lexicographically minimal rotation.

    Linear rather than the obvious quadratic scan over rotations, which matters
    here: the maize call set is millions of motifs and reaches 2,000 bp, where
    building and comparing every rotation does not finish.
    """
    s2 = s + s
    n2 = len(s2)
    f = [-1] * n2
    k = 0
    for j in range(1, n2):
        sj = s2[j]
        i = f[j - k - 1]
        while i != -1 and sj != s2[k + i + 1]:
            if sj < s2[k + i + 1]:
                k = j - i - 1
            i = f[i]
        if sj != s2[k + i + 1]:
            if sj < s2[k]:
                k = j
            f[j - k] = -1
        else:
            f[j - k] = i + 1
    return k


def canonical(motif):
    """Lexicographically smallest rotation of the motif or its complement."""
    m = motif.upper()
    if not m or any(c not in "ACGT" for c in m):
        return None
    best = None
    for variant in (m, revcomp(m)):
        n = len(variant)
        k = least_rotation(variant)
        rot = (variant + variant)[k:k + n]
        if best is None or rot < best:
            best = rot
    return best


def load_family_intervals():
    fam = {}
    for name, path in FAMILIES.items():
        by_chrom = defaultdict(list)
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                try:
                    by_chrom[p[0]].append((int(p[1]), int(p[2])))
                except ValueError:
                    continue
        for c in by_chrom:
            by_chrom[c].sort()
        fam[name] = by_chrom
    return fam


def overlap_bp(by_chrom, intervals_by_chrom):
    """Total bp of `intervals_by_chrom` falling inside `by_chrom`."""
    total = 0
    for chrom, ivs in intervals_by_chrom.items():
        ref = by_chrom.get(chrom, [])
        if not ref:
            continue
        j = 0
        for s, e in sorted(ivs):
            while j < len(ref) and ref[j][1] <= s:
                j += 1
            k = j
            while k < len(ref) and ref[k][0] < e:
                total += min(ref[k][1], e) - max(ref[k][0], s)
                k += 1
    return total


def main():
    bed = sys.argv[1] if len(sys.argv) > 1 else f"{B}/bwtandem/bwt_maize.bed"
    topn = int(sys.argv[2]) if len(sys.argv) > 2 else 25
    min_unit = int(sys.argv[3]) if len(sys.argv) > 3 else 0

    groups = defaultdict(lambda: {"bp": 0, "n": 0, "iv": defaultdict(list)})
    skipped = 0
    filtered = 0
    with open(bed) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                continue
            key = canonical(p[3])
            if key is None:
                skipped += 1
                continue
            # Unit-length floor, read off the tool's own motif. This stays
            # catalog-blind: it uses no curated coordinate, only the length of
            # the repeat unit the tool itself reported. Without it the ranking is
            # entirely microsatellite, since a 3 bp unit recurs across the genome
            # and outweighs any satellite family in total spanned sequence.
            if len(key) < min_unit:
                filtered += 1
                continue
            g = groups[key]
            g["bp"] += e - s
            g["n"] += 1
            g["iv"][p[0]].append((s, e))

    print(f"source : {bed}")
    print(f"groups : {len(groups):,} canonical motifs, "
          f"{skipped:,} calls skipped for non-ACGT motifs"
          + (f", {filtered:,} calls below the {min_unit} bp unit floor" if min_unit else "")
          + "\n")

    fam = load_family_intervals()
    ranked = sorted(groups.items(), key=lambda kv: -kv[1]["bp"])[:topn]

    print(f"{'rank':>4s} {'period':>7s} {'calls':>10s} {'span Mb':>9s} "
          f"{'motif (truncated)':24s} {'assignment':16s} {'in-family %':>11s}")
    hits = defaultdict(list)
    for i, (key, g) in enumerate(ranked, 1):
        best, best_frac = None, 0.0
        for name, iv in fam.items():
            ov = overlap_bp(iv, g["iv"])
            frac = ov / g["bp"] if g["bp"] else 0.0
            if frac > best_frac:
                best, best_frac = name, frac
        assign = best if best_frac >= 0.5 else "uncurated"
        if best_frac >= 0.5:
            hits[best].append(i)
        motif = key if len(key) <= 22 else key[:19] + "..."
        print(f"{i:>4d} {len(key):>7d} {g['n']:>10,d} {g['bp']/1e6:>9.3f} "
              f"{motif:24s} {assign:16s} {100*best_frac:>10.1f}%")

    print()
    for name in FAMILIES:
        r = hits.get(name)
        print(f"  {name:9s} recovered at rank(s) {r}" if r
              else f"  {name:9s} NOT recovered in the top {topn}")
    n_unc = topn - sum(len(v) for v in hits.values())
    print(f"\n  {n_unc} of the top {topn} clusters are uncurated; see the module "
          f"docstring for why that is not a false-positive count.")

    # How badly does exact-canonical grouping split a family? For each curated
    # family, the tool bases that fall inside it, and the share of those bases
    # held by its single largest cluster. A low share means the family is
    # present in the output but scattered across many canonical groups, which
    # caps how high any one of them can rank above.
    print(f"\n  {'family':9s} {'curated Mb':>11s} {'tool bp inside Mb':>18s} "
          f"{'clusters':>9s} {'largest cluster share':>22s}")
    for name, iv in fam.items():
        cur_bp = sum(e - s for v in iv.values() for s, e in v)
        shares = []
        for key, g in groups.items():
            ov = overlap_bp(iv, g["iv"])
            if ov:
                shares.append(ov)
        tot = sum(shares)
        top = max(shares) if shares else 0
        print(f"  {name:9s} {cur_bp/1e6:>11.2f} {tot/1e6:>18.2f} "
              f"{len(shares):>9,d} {(100.0*top/tot if tot else 0):>21.1f}%")


if __name__ == "__main__":
    main()
