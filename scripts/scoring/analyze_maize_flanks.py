#!/usr/bin/env python3
"""WP2.2 — are BWTandem's maize flanks real peripheral monomers or overcall?

BWTandem's calls extend beyond the curated coordinates of the 59 maize satellite
arrays. The draft has, at various points, guessed both ways: that the extension is
Tier 3 tolerating degenerate flanks (an excuse), and that curated coordinates are
conservative (a result). Neither was measured.

This measures it. For each curated array, take the sequence BWTandem calls outside
the curated interval and compute windowed autocorrelation identity at the family
period — 180 bp for knob180, 156 bp for TR-1, 156 bp for CentC — against two
references computed the same way:

  interior  the curated array itself, i.e. what a real array scores
  random    intergenic sequence away from any curated array, i.e. background

A flank scoring near interior is periodic sequence the curated coordinates omit. A
flank scoring near background is overcall, and the fix is a boundary trim.

Family periods: knob180 monomers are 180 bp and TR-1 is a 156 bp repeat
(Ananiev et al. 1998); CentC is 156 bp (Ananiev et al. 1998). Both 180 and 156 are
tested on every flank so a flank periodic at the *other* family's period is not
scored as background.

Usage: analyze_maize_flanks.py [--min-flank 500] [--max-flank 20000]
"""
import argparse
import os
import random
import sys
from collections import defaultdict

import numpy as np

sys.path.insert(0, "/data/gpfs/assoc/pgl/devel/bwt-algorithm")
from src.autocorr import autocorr_identity  # noqa: E402

GT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/ground_truth"
BWT = "/data/gpfs/assoc/pgl/filip/bwtandem_results/beds/bwtandem/bwt_maize.bed"
FASTA = ("/data/gpfs/assoc/pgl/filip/old_bwt/bwtandem/Zm/"
         "GCA_022117705.1_Zm-Mo17-REFERENCE-CAU-T2T-assembly_genomic.fna")

FAMILIES = [
    ("knob180", f"{GT}/mo17_knob180_arrays.bed", 180),
    ("TR-1",    f"{GT}/mo17_tr1_arrays.bed",     358),
    ("CentC",   f"{GT}/mo17_centc_arrays.bed",   156),
]
# knob180 monomers are 180 bp, TR-1 is ~358 bp, CentC is 156 bp. All three are
# tested on every window so a flank periodic at a neighbouring family's period is
# not scored as background. A first run used 156 for TR-1 as well, which put TR-1
# array interiors at 0.284 against a 0.268 background — i.e. the reference itself
# failed to separate from noise, which is the signature of a wrong period rather
# than of a non-periodic array.
TEST_PERIODS = (180, 358, 156)


def load_bed(path):
    out = defaultdict(list)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 3:
                try:
                    out[p[0]].append((int(p[1]), int(p[2])))
                except ValueError:
                    pass
    return out


def merge(iv):
    if not iv:
        return []
    iv = sorted(iv)
    out = [list(iv[0])]
    for s, e in iv[1:]:
        if s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]


def read_fasta(path, wanted):
    """Load only the sequences in `wanted` into memory."""
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name in wanted:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
                if len(seqs) == len(wanted):
                    break
            elif name in wanted:
                buf.append(line.strip())
        if name in wanted and name not in seqs:
            seqs[name] = "".join(buf).upper()
    return seqs


def best_identity(seq, periods=TEST_PERIODS):
    """Max autocorrelation identity over the candidate family periods.

    autocorr_identity takes a numpy byte array, not a string, and compares
    elementwise — so the sequence is encoded once here rather than per period.
    Windows that are mostly N are rejected: N==N counts as a match and would
    score a gap as perfectly periodic.
    """
    if len(seq) < 2 * max(periods):
        return None
    if seq.count("N") > len(seq) * 0.1:
        return None
    arr = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    vals = []
    for p in periods:
        try:
            vals.append(autocorr_identity(arr, p))
        except Exception:
            pass
    return max(vals) if vals else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-flank", type=int, default=500)
    ap.add_argument("--max-flank", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20260804)
    a = ap.parse_args()

    calls = load_bed(BWT)
    gts = {name: load_bed(path) for name, path, _ in FAMILIES}
    wanted = {c for g in gts.values() for c in g}
    print(f"loading {len(wanted)} sequences from the assembly ...", flush=True)
    seqs = read_fasta(FASTA, wanted)
    print(f"loaded {len(seqs)}: "
          + ", ".join(f"{k} {len(v)/1e6:.0f} Mb" for k, v in sorted(seqs.items()))
          + "\n", flush=True)

    rng = random.Random(a.seed)
    rows = []
    for fam, path, _ in FAMILIES:
        for chrom, arrays in gts[fam].items():
            seq = seqs.get(chrom)
            if not seq:
                continue
            merged = merge(calls.get(chrom, []))
            for gs, ge in arrays:
                # calls overlapping this array, merged
                ov = merge([(s, e) for s, e in merged if s < ge and e > gs])
                if not ov:
                    continue
                left = max(0, min(s for s, _ in ov))
                right = max(e for _, e in ov)
                for side, (fs, fe) in (("left", (left, gs)),
                                       ("right", (ge, right))):
                    span = fe - fs
                    if span < a.min_flank:
                        continue
                    fs2, fe2 = fs, min(fe, fs + a.max_flank)
                    ident = best_identity(seq[fs2:fe2])
                    if ident is None:
                        continue
                    rows.append((fam, chrom, gs, ge, side, fe2 - fs2, ident))

    # references
    interior, background = defaultdict(list), []
    for fam, path, _ in FAMILIES:
        for chrom, arrays in gts[fam].items():
            seq = seqs.get(chrom)
            if not seq:
                continue
            for gs, ge in arrays:
                sub = seq[gs:min(ge, gs + a.max_flank)]
                v = best_identity(sub)
                if v is not None:
                    interior[fam].append(v)
    allgt = defaultdict(list)
    for g in gts.values():
        for c, iv in g.items():
            allgt[c].extend(iv)
    for chrom, seq in seqs.items():
        occupied = merge(allgt[chrom])
        for _ in range(40):
            for _try in range(20):
                s = rng.randrange(0, max(1, len(seq) - a.max_flank))
                e = s + 10000
                if any(s < oe and e > os_ for os_, oe in occupied):
                    continue
                sub = seq[s:e]
                if sub.count("N") > len(sub) * 0.1:
                    continue
                v = best_identity(sub)
                if v is not None:
                    background.append(v)
                break

    def summ(v):
        if not v:
            return "n=0"
        v = sorted(v)
        return (f"n={len(v):3d}  median {v[len(v)//2]:.3f}  "
                f"q25 {v[len(v)//4]:.3f}  q75 {v[3*len(v)//4]:.3f}  "
                f"max {v[-1]:.3f}")

    print("=== reference identities (max over periods "
          + " / ".join(str(p) for p in TEST_PERIODS) + ") ===")
    for fam, _, _ in FAMILIES:
        print(f"  {fam+' interior':<20} {summ(interior[fam])}")
    print(f"  {'random intergenic':<20} {summ(background)}")

    bg = sorted(background)
    bg_p95 = bg[int(len(bg) * 0.95)] if bg else 1.0
    print(f"\n  background 95th percentile = {bg_p95:.3f} "
          f"(a flank above this is periodic, not background)\n")

    print("=== flanks called by BWTandem outside the curated coordinates ===")
    print(f"{'family':<9} {'n':>4} {'median':>8} {'q25':>7} {'q75':>7} "
          f"{'>bg95':>7} {'>=interior median':>18}")
    for fam, _, _ in FAMILIES:
        v = sorted(r[6] for r in rows if r[0] == fam)
        if not v:
            print(f"{fam:<9} {'n=0':>4}")
            continue
        im = sorted(interior[fam])
        imed = im[len(im) // 2] if im else 1.0
        above_bg = sum(1 for x in v if x > bg_p95)
        above_int = sum(1 for x in v if x >= imed)
        print(f"{fam:<9} {len(v):>4} {v[len(v)//2]:>8.3f} {v[len(v)//4]:>7.3f} "
              f"{v[3*len(v)//4]:>7.3f} {100*above_bg/len(v):>6.1f}% "
              f"{100*above_int/len(v):>17.1f}%")

    tot = sorted(r[6] for r in rows)
    if tot:
        print(f"\nall flanks: {summ(tot)}")
        print(f"  above background 95th pct: "
              f"{100*sum(1 for x in tot if x > bg_p95)/len(tot):.1f}%")
    print(f"\nflank spans measured: {len(rows)} "
          f"(min {a.min_flank} bp, truncated at {a.max_flank} bp)")


if __name__ == "__main__":
    main()
