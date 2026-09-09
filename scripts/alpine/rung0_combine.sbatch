#!/bin/bash
#SBATCH --job-name=rung0-combine
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Rung 0, stage 3 of 3: everything that needs the whole frame at once. Reads every cached slice
# (refusing if one is missing), scores both reliabilities, draws the nulls, writes the
# per-condition and dose-strata tables, pools the noise sums into the decomposition, and draws
# the figures into the task folder. No scan of the table happens here; the cache is the input.
#
# Submitted with --dependency=afterok:<array job id> so it waits for every slice.
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"
source scripts/alpine/rung0_env.sh

python scripts/delta_reproducibility.py \
    --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
    --stage combine \
    --slices "$RUNG0_SLICES" \
    --cache-dir "$RUNG0_CACHE" \
    --n-perm 500 \
    --out-dir "$OUT_DIR"
