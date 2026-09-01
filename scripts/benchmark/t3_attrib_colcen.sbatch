#!/bin/bash
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/t3attrib_%x_%j.out
#SBATCH --error=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/t3attrib_%x_%j.err

# Tier 3 attribution control (DIAGNOSTIC, not manuscript): two Tier-3-only
# arms on the SAME commit, so N handling is held fixed and only the adaptive
# parameterization differs -- max-period 2000 gives tolerance_ratio 0.0204,
# max-period 100000 restores the deposited runs' 0.0400. The wide arm is
# post-filtered to periods <= 2000 for comparison. Separates "Tier 3 moved
# because of tolerance" from "Tier 3 moved because of N handling".
#
# Submit: sbatch --job-name=t3a_narrow --export=ALL,ARM=narrow scripts/benchmark/t3_attrib_colcen.sbatch
#         sbatch --job-name=t3a_wide   --export=ALL,ARM=wide   scripts/benchmark/t3_attrib_colcen.sbatch
set -euo pipefail
: "${ARM:?set ARM=narrow|wide}"
REPO=/data/gpfs/assoc/pgl/devel/bwt-algorithm
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
OUT=/data/gpfs/assoc/pgl/devel/exp1_human/regen/t3attrib_${ARM}
FASTA=/data/gpfs/assoc/pgl/filip/bwt/colcen/colcen.fa
cd "$REPO"
echo "commit: $(git rev-parse HEAD)  arm: $ARM  started: $(date -Is)"
case "$ARM" in
  narrow) MAXP=2000 ;;
  wide)   MAXP=100000 ;;
  *) echo "FATAL: bad ARM"; exit 1 ;;
esac
$PY -m bwtandem.main "$FASTA" --tiers tier3 --min-period 1 --max-period $MAXP \
    --threads 2 --format bed -o "$OUT" -v
echo "BED lines: $(wc -l < "$OUT.bed")"
if [[ "$ARM" == wide ]]; then
    awk -F'\t' 'length($4) <= 2000' "$OUT.bed" > "$OUT.le2000.bed"
    echo "post-filtered to p<=2000: $(wc -l < "$OUT.le2000.bed") lines"
fi
