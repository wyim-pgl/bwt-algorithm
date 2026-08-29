#!/usr/bin/env python3
"""Stratified blinded sample for the human specificity audit (issue #11).

Protocol (adopted 2026-08-27): from the regenerated human BED, take the calls
that overlap neither the adotto catalog nor any comparator tool ("BWTandem-
only"), stratify by reported period into 1-6 / 7-20 / 21-100 / 101-2000 bp,
and draw 100 per stratum with a fixed seed. Two review files are written:
the reviewer sheet is SOURCE-blinded (tier, score and originating record are
withheld) but NOT period-blinded -- the proposed period is on the sheet
because the reviewer needs it to test the periodicity; describe the audit as
source-blinded, not fully blinded. The sealed key maps sample IDs back to the
full call records.

A call is judged SUPPORTED only when the proposed period produces >= 3
visible tandem copies in a dot plot / unit-shift alignment and stronger
within-call periodicity than both flanks; two reviewers, disagreements
adjudicated blind to stratum totals; report per-stratum confirmation with
Wilson 95% CIs.

Usage:
  sample_specificity_audit.py BWT_ONLY_BED FASTA OUTDIR [--per-stratum 100] [--seed 20260827]

BWT_ONLY_BED is produced first, e.g. with bedtools:
  bedtools intersect -v -a regen_human.bed -b adotto.bed trf.bed ultra.bed tantan.bed trash.bed
"""
import argparse
import os
import random
import sys
from collections import defaultdict

STRATA = ((1, 6), (7, 20), (21, 100), (101, 2000))


def read_fasta_index(fasta):
    """chrom -> (offset, length, linebases, linewidth) from a .fai (build if absent)."""
    fai = fasta + ".fai"
    if not os.path.exists(fai):
        sys.exit(f"FATAL: {fai} missing. Create it first: samtools faidx {fasta}")
    idx = {}
    with open(fai) as f:
        for line in f:
            p = line.split("\t")
            idx[p[0]] = (int(p[2]), int(p[1]), int(p[3]), int(p[4]))
    return idx


def fetch(fh, idx, chrom, start, end):
    """Binary-mode fetch. .fai offsets/line widths are BYTE offsets, so the
    file must be opened 'rb': text mode's newline translation makes seek()
    non-portable, and CRLF files would leave \r behind. The exact-length
    assertion also rejects a truncated FASTA instead of silently returning a
    short context."""
    if chrom not in idx:
        sys.exit(f"FATAL: chromosome {chrom!r} not in the .fai index")
    off, length, lb, lw = idx[chrom]
    start, end = max(0, start), min(length, end)
    fh.seek(off + start // lb * lw + start % lb)
    need = end - start
    out = []
    while need > 0:
        chunk = fh.read(need + need // lb + 2)
        if not chunk:
            break
        chunk = chunk.replace(b"\n", b"").replace(b"\r", b"")
        out.append(chunk[:need])
        need -= len(chunk[:need])
    seq = b"".join(out).decode("ascii")
    if len(seq) != end - start:
        sys.exit(f"FATAL: truncated FASTA -- wanted {end - start} bp at "
                 f"{chrom}:{start}-{end}, got {len(seq)}")
    return seq.upper()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("bwt_only_bed")
    ap.add_argument("fasta")
    ap.add_argument("outdir")
    ap.add_argument("--per-stratum", type=int, default=100)
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--allow-underfilled", action="store_true",
                    help="proceed when a stratum has fewer records than "
                         "--per-stratum (changes the preregistered 4x100 "
                         "design; must be disclosed in the audit methods)")
    args = ap.parse_args()

    strata = defaultdict(list)
    with open(args.bwt_only_bed) as f:
        for i, line in enumerate(f, 1):
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            try:
                s0, e0 = int(p[1]), int(p[2])
            except ValueError:
                sys.exit(f"FATAL: {args.bwt_only_bed}:{i}: non-numeric coordinates")
            if not 0 <= s0 < e0:
                sys.exit(f"FATAL: {args.bwt_only_bed}:{i}: invalid interval {s0}-{e0}")
            period = len(p[3])
            for lo, hi in STRATA:
                if lo <= period <= hi:
                    strata[f"{lo}-{hi}"].append((i, p))
                    break
    for k in (f"{lo}-{hi}" for lo, hi in STRATA):
        if not strata.get(k):
            sys.exit(f"FATAL: stratum {k} is empty -- wrong input?")

    short = {k: len(v) for k, v in strata.items() if len(v) < args.per_stratum}
    if short and not args.allow_underfilled:
        sys.exit("FATAL: underfilled strata " + str(short) +
                 f" (< {args.per_stratum}). This changes the preregistered "
                 "design; rerun with --allow-underfilled to accept and "
                 "disclose it.")

    rng = random.Random(args.seed)
    # The output directory must not exist: rerunning over a distributed
    # reviewer sheet or an already-sealed key would corrupt the audit.
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        sys.exit(f"FATAL: {args.outdir} already exists -- refusing to "
                 "overwrite a possibly distributed sheet/sealed key.")
    idx = read_fasta_index(args.fasta)
    fh = open(args.fasta, "rb")

    sheet = open(os.path.join(args.outdir, "reviewer_sheet.tsv"), "w")
    key = open(os.path.join(args.outdir, "answer_key.tsv"), "w")
    sheet.write("sample_id\tchrom\tstart\tend\tperiod_bp\t"
                "left_flank_bp\tright_flank_bp\tmasked_frac\t"
                "sequence_with_flanks\n")
    key.write("sample_id\tsource_line\tfull_record\n")

    sampled = []
    for name, pool in sorted(strata.items()):
        take = rng.sample(pool, min(args.per_stratum, len(pool)))
        if len(take) < args.per_stratum:
            print(f"WARNING: stratum {name} has only {len(pool)} calls "
                  f"(< {args.per_stratum}); taking all", file=sys.stderr)
        sampled.extend((name, rec) for rec in take)
    rng.shuffle(sampled)  # blind the reviewer to stratum ordering

    for n, (name, (line_no, p)) in enumerate(sampled, 1):
        sid = f"S{n:03d}"
        chrom, s, e = p[0], int(p[1]), int(p[2])
        period = len(p[3])
        flank = period  # one repeat unit of flank each side, per protocol
        lf = min(flank, s)                       # clipped at the contig start
        rf = min(flank, idx[chrom][1] - e)       # clipped at the contig end
        seq = fetch(fh, idx, chrom, s - lf, e + rf)
        masked = sum(1 for c in seq if c not in "ACGT")
        sheet.write(f"{sid}\t{chrom}\t{s}\t{e}\t{period}\t{lf}\t{rf}\t"
                    f"{masked / max(1, len(seq)):.3f}\t{seq}\n")
        key.write(f"{sid}\t{line_no}\t" + "\t".join(p) + "\n")

    sheet.close(); key.close()
    os.chmod(os.path.join(args.outdir, "answer_key.tsv"), 0o600)
    print(f"{len(sampled)} samples -> {args.outdir}/reviewer_sheet.tsv "
          f"(blinded) + answer_key.tsv (keep sealed until adjudication)")
    print(f"strata: " + ", ".join(f"{k}:{len(v)}" for k, v in sorted(strata.items())))
    print(f"seed {args.seed} -- record it in the audit methods.")


if __name__ == "__main__":
    main()
