#!/bin/bash
#SBATCH --job-name=rung0-slice
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=72G                    # Alpine charges memory per core at 3,840 MB, so memory IS
                                     # the core request: 120G was billed as 32 cores and sat in
                                     # the queue. Peak RSS runs 35-40G ABOVE the engine's own
                                     # limit whatever the slice size (98G at 60G, job
                                     # 32341611_0; 75G at 40G, job 32345592_0, killed), so the
                                     # engine limit below is what keeps a task inside this
#SBATCH --time=01:30:00              # short enough to backfill into the scheduler's gaps; a
                                     # sixteenth slice took 51 min at 4 threads
#SBATCH --array=0-31
#SBATCH --output=logs/%x-%A_%a.out
#SBATCH --error=logs/%x-%A_%a.out
#
# Rung 0, stage 2 of 3: one job-array task per slice of the genes. Each task scans the whole
# table once, keeps the genes whose hash lands in its slice, aggregates them by (line, drug,
# dose, gene) in ONE group-by that serves both the reliabilities (the two half means) and the
# noise decomposition (the variance across plates and the mean squared standard error), and
# caches its slice's frame, per-condition noise sums and a bounded sample. Thirty-two tasks
# make the wall clock one scan instead of thirty-two in a row, and each task holds a
# thirty-second of the group table. Smaller tasks slot into the queue sooner: the first
# attempt at sixteen slices of 120G each was billed as sixteen 32-core jobs and trickled in
# one at a time. A task that finds its completion record already written skips, so a failed
# index is resubmitted alone with --array=<index>.
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
    --duckdb-memory 28GB \
    --duckdb-threads 4 \
    --noise-sample-rows 2000000 \
    --out-dir "$OUT_DIR"
