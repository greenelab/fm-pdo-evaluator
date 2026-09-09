#!/bin/bash
#SBATCH --job-name=rung0-assign
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G                   # the first attempt at 64G with the engine at 48G was
                                     # OOM-killed by the cgroup (job 32341467): the parquet
                                     # readers and the group keys sit outside the engine's
                                     # limit, so the job needs twice the engine's allowance
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Rung 0, stage 1 of 3: the split assignment. One scan over the key columns alone (line, drug,
# dose, plate, and baseMean for the untestable share) -- about a quarter of an hour -- and the
# plate-to-half assignment becomes an explicit table every slice joins against: plates sorted by
# id within each (line, drug, dose) triple, assigned alternately. It also writes the pool
# description (plates per triple, per half) into the task folder.
#
# Chain: this job -> rung0_slice.sbatch (array, afterok) -> rung0_combine.sbatch (afterok).
# scripts/alpine/submit_rung0_chain.sh submits all three with the dependencies set.
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"
source scripts/alpine/rung0_env.sh

python scripts/delta_reproducibility.py \
    --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
    --stage assign \
    --cache-dir "$RUNG0_CACHE" \
    --duckdb-memory 60GB \
    --duckdb-threads 4 \
    --out-dir "$OUT_DIR"
