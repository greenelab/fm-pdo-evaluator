#!/bin/bash
#SBATCH --job-name=permutation-null
#SBATCH --partition=acpu      # was amilan; CURC retired that partition (2026-08-24)
#SBATCH --qos=cpu-normal      # acpu accepts only cpu-normal / cpu-long, not normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G                    # DuckDB memory_limit 36G + the aggregated result frame
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Rung-0 permutation-based exact permutation null (docs/tasks/rung0-replicate-ceiling/design.md,
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

# Activate the `stack` env by PATH instead of the module system. `module load anaconda`
# fails inside batch jobs (Lmod: module unknown) for reasons NOT established -- see
# ralpine's submit comment; do not assert a cause here. This pattern is the one that is
# PROVEN to work: 10_sarcoma_organoids_2024_context and 12_sarcoma_organoids_2024_score ran with it on 2026-08-22.
# `conda activate` is deliberately not used: this install ships no etc/profile.d/conda.sh.
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export CONDA_PREFIX="/projects/$USER/software/anaconda/envs/stack"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1   # else stdout block-buffers to the log file and a multi-hour
                           # job shows 0 bytes until it exits, making "working" and
                           # "wedged" indistinguishable (job 31655278, 2026-08-25).
export HF_HUB_CACHE="/scratch/alpine/$USER/hf"   # cache off /projects
# do NOT set HF_HOME to scratch -- that hides the token `hf auth login` saved under the HF cache:
export HF_TOKEN="${HF_TOKEN:-$(cat "${XDG_CACHE_HOME:-$HOME/.cache}/huggingface/token" 2>/dev/null || true)}"

# Pinned to the pyproject floor with an upper bound so a fallback install can't drift the
# persistent stack env to an untested major version. Lock-based provisioning (containers/uv on
# Alpine) is a declared follow-up, not done here.
python -c "import duckdb" 2>/dev/null || pip install -q "duckdb>=1.0,<2"

# Run once per gene set. Rung 0 promotes two reliabilities, and the permutation check measures
# the dependence the bootstrap ignores for each of them -- a design effect established on the
# all-gene statistic says nothing about the responder one, whose draws reuse a different and
# smaller set of genes per condition. No drug list and no panel, for the reason given in
# delta_reproducibility.sbatch.
OUT_DIR="${OUT_DIR:-docs/tasks/rung0-assay-reliability}"
for GENE_SET in all responder; do
    echo "=== permutation null, gene set: $GENE_SET ==="
    python scripts/permutation_null.py \
        --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
        --gene-set "$GENE_SET" \
        --out-dir "$OUT_DIR"
done
