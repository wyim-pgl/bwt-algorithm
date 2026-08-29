#!/bin/bash
#SBATCH --account=cpu-s1-pgl-0
#SBATCH --partition=cpu-s1-pgl-0
#SBATCH --job-name=triage_maize
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/triage_%x_%j.out
#SBATCH --error=/data/gpfs/assoc/pgl/devel/exp1_human/regen/logs/triage_%x_%j.err

# Pre-rescoring triage of the regenerated maize BED (gate passed 2026-08-29).
# UNTRACKED during the campaign, same as the other benchmark scripts.
set -euo pipefail
REPO=/data/gpfs/assoc/pgl/devel/bwt-algorithm
R=/data/gpfs/assoc/pgl/devel/exp1_human/regen
PY=/data/gpfs/assoc/pgl/bin/conda/conda_envs/bwtandem/bin/python
cd "$REPO"
echo "commit: $(git rev-parse HEAD)  started: $(date -Is)"
$PY scripts/benchmark/triage_regen.py \
    "$R/regen_maize.bed" \
    "$REPO/results/beds/bwtandem_maize.bed.gz" \
    /data/gpfs/assoc/pgl/filip/bwt/zmays/zmays.fna \
    --provenance "$R/regen_maize.provenance.json" \
    | tee "$R/analysis/triage_maize.txt"
echo "triage done: $R/analysis/triage_maize.txt"
