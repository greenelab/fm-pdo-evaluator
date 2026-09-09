#!/bin/bash
#SBATCH --job-name=tahoe-pull
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G                    # snapshot_download bookkeeping only (no in-engine aggregation); modest
#SBATCH --time=12:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# One authenticated bulk pull of Tahoe's pseudobulk DE config to scratch, as a SINGLE process --
# deliberately NOT an array. The table is a flat 1026-file shard set with no drug partition, so
# there is no work to split across shards; concurrent readers would only trip the HF rate limit.
#
# The archived lineage's variant also built the rung-1 delta bundle; that half arrives with rung 1.
#
# Authenticate first so the download isn't rate-limited:  hf auth login
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"

module load anaconda
conda activate stack                 # data env: datasets + pyarrow + huggingface_hub + duckdb
export PYTHONPATH="$REPO/src"
export HF_HUB_CACHE="/scratch/alpine/$USER/hf"   # cache large downloads on scratch, not /projects
# do NOT set HF_HOME to scratch -- that hides the token `hf auth login` saved under the HF cache:
export HF_TOKEN="${HF_TOKEN:-$(cat "${XDG_CACHE_HOME:-$HOME/.cache}/huggingface/token" 2>/dev/null || true)}"

python scripts/download_tahoe_pseudobulk_de.py \
    --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de"
