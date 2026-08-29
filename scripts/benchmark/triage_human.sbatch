#!/bin/bash
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --job-name=triage_human
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/triage_%x_%j.out
#SBATCH --error=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/triage_%x_%j.err

# Pre-rescoring triage of the regenerated human BED (workflow: gate -> triage
# -> score; the gate passed 2026-08-28, scoring already ran -- this closes the
# skipped middle step). UNTRACKED during the campaign, same as the other
# benchmark scripts.
set -euo pipefail
REPO=/data/gpfs/assoc/pgl/devel/bwt-algorithm
R=/data/gpfs/assoc/pgl/devel/exp1_human/regen
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
cd "$REPO"
echo "commit: $(git rev-parse HEAD)  started: $(date -Is)"
$PY scripts/benchmark/triage_regen.py \
    "$R/regen_human.bed" \
    "$REPO/results/beds/bwtandem_human.bed.gz" \
    /data/gpfs/assoc/pgl/devel/exp1_human/data/hg38_primary.fa \
    --provenance "$R/regen_human.provenance.json" \
    | tee "$R/analysis/triage_human.txt"
echo "triage done: $R/analysis/triage_human.txt"
