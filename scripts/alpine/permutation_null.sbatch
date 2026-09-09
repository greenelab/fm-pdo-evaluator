#!/bin/bash
#SBATCH --job-name=permutation-null
#SBATCH --partition=acpu      # was amilan; CURC retired that partition (2026-08-24)
#SBATCH --qos=cpu-normal      # acpu accepts only cpu-normal / cpu-long, not normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G                    # DuckDB memory_limit 36G + the aggregated result frame
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Rung-0 permutation-based exact permutation null (docs/tasks/rung0-assay-reliability/design.md,
# decision history). Reuses the same DuckDB split-half build as delta_reproducibility.sbatch,
# then runs THREE permutation-based nulls over the half-profile pivots, carrying the
# half-profile-sharing dependence the bootstrapped p-values had to assume away by construction
# in each case: any-pair (a permutation with no fixed points, validating the pooled aggregate),
# same-drug (a per-drug permutation, carrying the same-drug stratum's ~32-drug clustering that
# `p_vs_same_drug` is read against), and diff-drug (a permutation constrained to change both
# line and drug, matching the `p_vs_null` stratum). CPU only; reads the Tahoe DE parquet
# already on scratch. Same entry point and CLI as before -- no flag changes for the extension.
#
# Authenticate first (drug_metadata is fetched from HF):  hf auth login
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"

# Run once per gene set. Rung 0 promotes two reliabilities, and the permutation check measures
# the dependence the bootstrap ignores for each of them -- a design effect established on the
# all-gene statistic says nothing about the responder one, whose draws reuse a different and
# smaller set of genes per condition. No drug list and no panel, for the reason given in
# delta_reproducibility.sbatch.
source scripts/alpine/rung0_env.sh
for GENE_SET in all responder; do
    echo "=== permutation null, gene set: $GENE_SET ==="
    python scripts/permutation_null.py \
        --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
        --cache-dir "$RUNG0_CACHE" \
        --slices "$RUNG0_SLICES" \
        --gene-set "$GENE_SET" \
        --out-dir "$OUT_DIR"
done
