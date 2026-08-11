#!/usr/bin/env python3
"""WP0 — rebuild the human mreps BED from the raw mreps output.

The published Table 1 mreps row is scored from a BED whose chromosome column is
broken: of 18,071,983 calls carrying a `chr*` name, **every single one says
chr4**, and their coordinates stop at 48.6 Mb (chr4 is 190.2 Mb). Since the
ground truth spans all 24 chromosomes, a tool whose calls all land on one
chromosome cannot score above that chromosome's share of the catalog — which is
what the published 1.75 % region recall actually measures.

The raw mreps output is intact: 455 `Processing sequence '<id> ...'` records,
with real GenBank accessions (CM000664.2 = chr2, etc.) and unplaced contigs
(KI270518.1, ...). So this is a conversion defect, not a run defect, and it is
fixed by re-converting rather than by re-running mreps (which would cost hours).

Two things the original converter gets wrong on this file:

1. mreps truncates the FASTA header to 80 characters, so the sequence name line
   reads `'CM000664.2 Homo sapiens chromosome 2, GRCh38 reference primary assemb'`
   — the closing quote is missing. `line.split("'")[1]` still yields the right
   text here, but the missing terminator makes the parse fragile.
2. The run concatenated every split file into ONE output stream. Any sequence
   whose header line failed to parse silently inherits the previous sequence's
   name, so one bad parse mislabels everything after it. That is consistent with
   a single chromosome name absorbing millions of calls.

This parser therefore refuses to inherit: a repeat record emitted before any
sequence header, or after a header it could not parse, is counted and dropped
rather than attributed to the wrong chromosome.

Usage:
  reconvert_mreps_human.py IN.txt OUT.bed [--primary-only]
"""
import argparse
import os
import sys
import time

# GRCh38 GenBank (GCA_000001405.15) primary chromosome accessions.
GENBANK = {
    'CM000663.2': 'chr1',  'CM000664.2': 'chr2',  'CM000665.2': 'chr3',
    'CM000666.2': 'chr4',  'CM000667.2': 'chr5',  'CM000668.2': 'chr6',
    'CM000669.2': 'chr7',  'CM000670.2': 'chr8',  'CM000671.2': 'chr9',
    'CM000672.2': 'chr10', 'CM000673.2': 'chr11', 'CM000674.2': 'chr12',
    'CM000675.2': 'chr13', 'CM000676.2': 'chr14', 'CM000677.2': 'chr15',
    'CM000678.2': 'chr16', 'CM000679.2': 'chr17', 'CM000680.2': 'chr18',
    'CM000681.2': 'chr19', 'CM000682.2': 'chr20', 'CM000683.2': 'chr21',
    'CM000684.2': 'chr22', 'CM000685.2': 'chrX',  'CM000686.2': 'chrY',
}
# RefSeq equivalents, in case the same script is pointed at the GCF output.
REFSEQ = {f'NC_0000{i:02d}.{v}': f'chr{i}' for i, v in [
    (1, 11), (2, 12), (3, 12), (4, 12), (5, 10), (6, 12), (7, 14), (8, 11),
    (9, 12), (10, 11), (11, 10), (12, 12), (13, 11), (14, 9), (15, 10),
    (16, 10), (17, 11), (18, 10), (19, 10), (20, 11), (21, 9), (22, 11)]}
REFSEQ['NC_000023.11'] = 'chrX'
REFSEQ['NC_000024.10'] = 'chrY'

PRIMARY = set(GENBANK.values())


def clean_chrom(raw):
    tok = raw.split()[0] if raw.split() else ''
    return REFSEQ.get(tok, GENBANK.get(tok, tok))


def parse(in_path, out_path, primary_only):
    chrom = None                 # None = "no valid sequence header seen yet"
    n_seq = n_out = n_orphan = n_bad = n_skip = 0
    seen = {}
    size = os.path.getsize(in_path)
    t0 = time.time()
    last = 0

    with open(in_path, 'r', errors='replace') as fin, \
            open(out_path, 'w') as fout:
        for line in fin:
            if line.startswith('Processing sequence'):
                # mreps truncates this at 80 chars, so the closing quote may be
                # absent — take everything after the opening quote.
                q = line.find("'")
                if q < 0:
                    chrom = None
                    n_bad += 1
                    continue
                chrom = clean_chrom(line[q + 1:])
                n_seq += 1
                seen[chrom] = 0
                continue

            if '->' not in line or ':' not in line or '<' not in line:
                continue

            if chrom is None:
                # A repeat record with no sequence context. Never inherit — that
                # is precisely how the published BED ended up all-chr4.
                n_orphan += 1
                continue

            if primary_only and chrom not in PRIMARY:
                n_skip += 1
                continue

            try:
                left, right = line.split(':', 1)
                s_str, e_str = left.split('->')
                start = int(s_str.strip()) - 1      # mreps is 1-based inclusive
                end = int(e_str.strip())
                period = right.split('<')[1].split('>')[0].strip()
                motif = 'N/A'
                if ']' in right:
                    toks = right.split(']')[1].split()
                    if len(toks) >= 2:
                        motif = ''.join(toks[1:])
            except (ValueError, IndexError):
                n_bad += 1
                continue

            if end <= start:
                n_bad += 1
                continue

            fout.write(f"{chrom}\t{start}\t{end}\t{motif}\t{period}\tmreps\n")
            n_out += 1
            seen[chrom] = seen.get(chrom, 0) + 1

            if n_out - last >= 20_000_000:
                last = n_out
                el = time.time() - t0
                print(f"  ... {n_out:,d} written, {el/60:.1f} min elapsed",
                      flush=True)

    return n_seq, n_out, n_orphan, n_bad, n_skip, seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input')
    ap.add_argument('output')
    ap.add_argument('--primary-only', action='store_true',
                    help='keep only chr1-22,X,Y (what Table 1 scores)')
    a = ap.parse_args()

    print(f"input      : {a.input} ({os.path.getsize(a.input):,d} bytes)")
    print(f"output     : {a.output}")
    print(f"primary    : {a.primary_only}")
    sys.stdout.flush()

    t0 = time.time()
    n_seq, n_out, n_orphan, n_bad, n_skip, seen = parse(
        a.input, a.output, a.primary_only)
    el = time.time() - t0

    print(f"\nsequences parsed   : {n_seq}")
    print(f"records written    : {n_out:,d}")
    print(f"orphan (no header) : {n_orphan:,d}   <- dropped, not inherited")
    print(f"unparseable        : {n_bad:,d}")
    print(f"non-primary skipped: {n_skip:,d}")
    print(f"elapsed            : {el/60:.1f} min")
    prim = {k: v for k, v in seen.items() if k in PRIMARY}
    print(f"\nprimary chromosomes with calls: {len(prim)} / 24")
    for k in sorted(prim, key=lambda c: (len(c), c)):
        print(f"  {k:<6} {prim[k]:,d}")
    if len(prim) < 24:
        print("\nWARNING: not all 24 primary chromosomes carry calls.")


if __name__ == '__main__':
    main()
