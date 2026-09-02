#!/usr/bin/env python3
"""Convert longdust and AniAnn's raw outputs to the scorer BED contract (#16).

Contract (period_of in score_table1.py / score_colcen.py): column 5 is the
integer period, falling back to len(column 4).

- **longdust** emits 3-column intervals and, by its own paper, cannot report
  repeat units — the converted BED stays 3-column ON PURPOSE so every
  period-conditional metric is structurally empty for it rather than faked.
- **AniAnn's** writes one 9-column BED per sequence (header line, score
  column = array periodicity/monomer length). Converted:
  `chrom  start  end  .  <period>  AniAnns` — coordinates and period carry
  over; there is no motif sequence to carry (per-copy structure is not part
  of its output).

Usage:
  convert_2026_tools.py longdust IN.bed OUT.bed
  convert_2026_tools.py anianns  IN_DIR OUT.bed
"""
import glob
import os
import sys


def convert_longdust(src, dst):
    n = 0
    with open(src) as fin, open(dst, "w") as fout:
        for i, line in enumerate(fin, 1):
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                s, e = int(p[1]), int(p[2])
            except ValueError:
                sys.exit(f"FATAL {src}:{i}: bad coordinates")
            if e <= s:
                sys.exit(f"FATAL {src}:{i}: empty/inverted interval")
            fout.write(f"{p[0]}\t{s}\t{e}\n")
            n += 1
    return n


def convert_anianns(src_dir, dst):
    beds = sorted(glob.glob(os.path.join(src_dir, "*.bed")))
    if not beds:
        sys.exit(f"FATAL: no per-sequence BEDs under {src_dir}")
    n = 0
    with open(dst, "w") as fout:
        for bed in beds:
            with open(bed) as fin:
                for i, line in enumerate(fin, 1):
                    if line.startswith("#"):
                        continue
                    p = line.rstrip("\n").split("\t")
                    if len(p) < 5:
                        continue
                    try:
                        s, e = int(p[1]), int(p[2])
                        period = int(float(p[4]))
                    except ValueError:
                        sys.exit(f"FATAL {bed}:{i}: bad row")
                    if e <= s:
                        sys.exit(f"FATAL {bed}:{i}: empty/inverted interval")
                    fout.write(f"{p[0]}\t{s}\t{e}\t.\t{period}\tAniAnns\n")
                    n += 1
    if n == 0:
        sys.exit(f"FATAL: {src_dir} produced zero rows")
    return n


def main():
    if len(sys.argv) != 4 or sys.argv[1] not in ("longdust", "anianns"):
        sys.exit(__doc__)
    tool, src, dst = sys.argv[1:4]
    n = convert_longdust(src, dst) if tool == "longdust" else convert_anianns(src, dst)
    print(f"{tool}: {n} rows -> {dst}")


if __name__ == "__main__":
    main()
