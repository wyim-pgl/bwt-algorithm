#!/bin/bash
#SBATCH --job-name=ultra_h_p2000
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=13-23:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/regen/ultra_p2000/ultra_h_p2000_%j.log
# Range-cost measurement: the published human ULTRA run was
# `ultra -t 2 -o <tsv> <fasta>` (default -p 100, 29.8 h, GNU-time 1.68 GB).
# This run changes ONLY the maximum period to 2000, matching BWTandem's
# search ceiling, on the same input FASTA the competitor baselines consumed.
# If it does not finish inside the partition limit, the hourly progress log
# below still documents how far it got (cf. the TRF 2000 bp attempt,
# terminated incomplete at 6.6 days).
set -uo pipefail
ULTRA=/data/gpfs/home/wyim/micromamba/bin/ultra
FASTA=/data/gpfs/assoc/pgl/filip/bwtandem_results/human.fna
OUT=/data/gpfs/assoc/pgl/devel/exp1_human/regen/ultra_p2000/ultra_human_p2000.tsv

echo "=== provenance ==="
echo "host       : $(hostname)"
echo "jobid      : ${SLURM_JOB_ID:-none}"
echo "binary     : $ULTRA ($(stat -c %y "$ULTRA"))"
echo "version    : $("$ULTRA" --help 2>&1 | grep -i version | tr -s ' ')"
echo "fasta      : $FASTA ($(stat -c %s "$FASTA") bytes; RefSeq-accession flavour, 455 sequences incl. scaffolds — the same file the published competitor baselines consumed)"
echo "command    : $ULTRA -t 2 -p 2000 -o $OUT $FASTA"
echo "started    : $(date -Is)"
echo "=================="

( while true; do sleep 3600; s=$(stat -c %s "$OUT" 2>/dev/null || echo 0); echo "PROGRESS $(date -Is) out_bytes=$s"; done ) &
MON=$!
/usr/bin/time -v "$ULTRA" -t 2 -p 2000 -o "$OUT" "$FASTA"
st=$?
kill $MON 2>/dev/null
echo "finished   : $(date -Is) exit=$st"
exit $st
