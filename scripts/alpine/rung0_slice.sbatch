#!/bin/bash
#SBATCH --job-name=rung0-slice
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G                   # twice the engine's limit: what sits outside it (parquet
                                     # readers, the returned Arrow table, pandas) killed the
                                     # assign job at 64G/48G (job 32341467)
#SBATCH --time=04:00:00
#SBATCH --array=0-15
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.out
#
# Rung 0, stage 2 of 3: one job-array task per slice of the genes. Each task scans the whole
# table once, keeps the genes whose hash lands in its slice, aggregates them by (line, drug,
# dose, gene) in ONE group-by that serves both the reliabilities (the two half means) and the
# noise decomposition (the variance across plates and the mean squared standard error), and
# caches its slice's frame, per-condition noise sums and a bounded sample. Sixteen tasks on
# sixteen nodes make the wall clock one scan instead of sixteen in a row, and each task holds a
# sixteenth of the group table. A task that finds its completion record already written skips,
# so a failed index is resubmitted alone with --array=<index>.
#
# --slices here MUST equal the array size, and the combine job must pass the same number.
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"
source scripts/alpine/rung0_env.sh

python scripts/delta_reproducibility.py \
    --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
    --stage slice \
    --slices "$RUNG0_SLICES" \
    --slice-index "$SLURM_ARRAY_TASK_ID" \
    --cache-dir "$RUNG0_CACHE" \
    --duckdb-memory 60GB \
    --duckdb-threads 4 \
    --noise-sample-rows 2000000 \
    --out-dir "$OUT_DIR"
