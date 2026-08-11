#!/usr/bin/env python3
"""WP3.1 — what are the regions only BWTandem reports?

Table 1 credits BWTandem with 863,654 regions no other tool calls, nearly twice the
count unique to tantan or ULTRA, and §4.4 already notes the count is ambiguous: a
shared tendency toward false positives explains it as readily as novel detection.
The BEDs can settle which.

For every call BWTandem makes that no other tool overlaps, this reports the period
distribution, the length distribution, motif entropy, and whether the region falls
inside the adotto catalog. The comparison that matters is against BWTandem's
*non*-unique calls: if the unique set is systematically shorter, lower-entropy and
outside the catalog, it is over-calling in low-complexity sequence. If it looks like
the rest of the output but sits in catalog gaps, it is detection the other tools
missed.

Usage: analyze_unique_regions.py [--tool BED] [--label NAME]
"""
import argparse
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
EXP1 = os.path.dirname(HERE)
BEDTOOLS = os.environ.get(
    "BEDTOOLS", "/data/gpfs/assoc/pgl/bin/bedtools2/bin/bedtools")
GT = os.path.join(EXP1, "data/adotto_primary.bed")
FILIP = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds"
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

OTHERS = [
    f"{HERE}/beds/ultra_human_GCA.bed",
    f"{FILIP}/tantan/GCA_000001405.15_GRCh38_genomic_output.bed",
    f"{FILIP}/trf/GCA_000001405.15_GRCh38_genomic_output.bed",
    f"{FILIP}/trash/GCA_000001405.15_GRCh38_genomic_trash.bed",
]


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


def entropy(seq):
    if not seq:
        return 0.0
    c = Counter(seq)
    n = len(seq)
    return -sum((v / n) * math.log2(v / n) for v in c.values())


def prep(src, dst):
    n = 0
    tmp = dst + ".u"
    with open(src) as fin, open(tmp, "w") as fout:
        for line in fin:
            p = line.rstrip("\n").split("\t")
            if len(p) < 3 or p[0] not in CHROMS:
                continue
            fout.write(line)
            n += 1
    with open(dst, "w") as out:
        subprocess.run([BEDTOOLS, "sort", "-i", tmp], stdout=out, check=True)
    os.remove(tmp)
    return n


def describe(label, rows):
    if not rows:
        print(f"  {label}: (none)")
        return
    lens = sorted(r[0] for r in rows)
    pers = [r[1] for r in rows if r[1] is not None]
    ents = [r[2] for r in rows]
    n = len(rows)

    def q(v, f):
        return v[min(len(v) - 1, int(len(v) * f))]

    print(f"  {label}: n={n:,d}  bp={sum(lens):,d}")
    print(f"    length  median {q(lens,0.5)}  q25 {q(lens,0.25)}  q75 {q(lens,0.75)}"
          f"  max {lens[-1]:,d}")
    pc = Counter()
    for p in pers:
        pc["1-2" if p <= 2 else "3-6" if p <= 6 else "7-20" if p <= 20
           else "21-100" if p <= 100 else ">100"] += 1
    dist = "  ".join(f"{k}:{100*v/max(1,len(pers)):.1f}%"
                     for k, v in sorted(pc.items(), key=lambda x: -x[1]))
    print(f"    period  {dist}")
    print(f"    motif entropy  mean {sum(ents)/n:.3f}  "
          f"<1.0 bit: {100*sum(1 for e in ents if e < 1.0)/n:.1f}%")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tool", default=f"{HERE}/out/bwt_hg38_catchF_p100.bed")
    ap.add_argument("--label", default="BWTandem-F")
    ap.add_argument("--workdir", default=f"{HERE}/work_wp31")
    a = ap.parse_args()

    os.makedirs(a.workdir, exist_ok=True)
    tool = os.path.join(a.workdir, "tool.bed")
    n_tool = prep(a.tool, tool)
    others = []
    for i, o in enumerate(OTHERS):
        if not os.path.exists(o):
            continue
        d = os.path.join(a.workdir, f"other{i}.bed")
        prep(o, d)
        others.append(d)
    gt = os.path.join(a.workdir, "gt.bed")
    prep(GT, gt)

    print(f"tool  : {a.label}  ({a.tool})")
    print(f"calls : {n_tool:,d}")
    print(f"others: {len(others)} tool BEDs\n")
    sys.stdout.flush()

    # unique = overlaps none of the other tools
    uniq = os.path.join(a.workdir, "unique.bed")
    with open(uniq, "w") as out:
        subprocess.run([BEDTOOLS, "intersect", "-a", tool, "-b", *others, "-v"],
                       stdout=out, check=True)
    # and the complement, for comparison
    shared = os.path.join(a.workdir, "shared.bed")
    with open(shared, "w") as out:
        subprocess.run([BEDTOOLS, "intersect", "-a", tool, "-b", *others, "-u"],
                       stdout=out, check=True)

    def rows_of(path, in_gt=None):
        rows = []
        with open(path) as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if len(p) < 3:
                    continue
                rows.append((int(p[2]) - int(p[1]), period_of(p),
                             entropy(p[3] if len(p) > 3 else "")))
        return rows

    uniq_rows = rows_of(uniq)
    shared_rows = rows_of(shared)
    print(f"unique to {a.label}: {len(uniq_rows):,d} "
          f"({100*len(uniq_rows)/max(1,n_tool):.1f}% of calls)")
    print(f"shared with >=1 tool: {len(shared_rows):,d}\n")

    describe("UNIQUE", uniq_rows)
    print()
    describe("SHARED", shared_rows)

    # how many unique calls are inside the catalog?
    for name, src in (("unique", uniq), ("shared", shared)):
        r = subprocess.run([BEDTOOLS, "intersect", "-a", src, "-b", gt, "-u"],
                           capture_output=True, text=True, check=True)
        hit = sum(1 for x in r.stdout.split("\n") if x.strip())
        tot = len(uniq_rows) if name == "unique" else len(shared_rows)
        print(f"\n  {name} calls inside adotto: {hit:,d} / {tot:,d} "
              f"({100*hit/max(1,tot):.1f}%)")

    print("\nReading: if UNIQUE is much shorter, lower-entropy and largely outside "
          "the catalog\nrelative to SHARED, the unique count reflects "
          "low-complexity over-calling rather than\nnovel detection.")


if __name__ == "__main__":
    main()
