#!/bin/bash
#SBATCH --job-name=register-tranche
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# One-time: ingest the Tahoe pseudobulk DE pool as tranche tahoe100m-pseudobulk-de.v1.
# Re-hashes all 1,026 shards (~83G, IO-bound) and cross-checks each against the sha256
# etag recorded at download time, so corruption since 2026-07-24 fails loudly here.
# The record is pulled back and committed by hand -- see the task plan.
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"

# PATH-based env activation; `module load anaconda` fails in batch jobs (cause not
# established -- see ralpine's submit comment; do not assert one here).
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export CONDA_PREFIX="/projects/$USER/software/anaconda/envs/stack"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1

# Pinned to the pyproject floor with an upper bound so a fallback install can't drift the
# persistent stack env to an untested major version. Lock-based provisioning (containers/uv on
# Alpine) is a declared follow-up, not done here.
python -c "import pydantic" 2>/dev/null || pip install -q "pydantic>=2.7,<3"

python scripts/register_tranche.py \
    --data-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
    --tranche-id tahoe100m-pseudobulk-de.v1 \
    --source "tahoebio/Tahoe-100M:pseudobulk_differential_expression" \
    --ingestion-date 2026-07-24 \
    --sample-count 50 --drug-count 32 \
    --description "Tahoe-100M pseudobulk DE shards on scratch; sample_count = cell lines (no patients in a cell-line corpus); drug_count = the declared drug panel. See docs/DATA.md." \
    --out data/tranches/tahoe100m-pseudobulk-de.v1.json
