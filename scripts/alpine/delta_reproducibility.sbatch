#!/bin/bash
#SBATCH --job-name=delta-repro
#SBATCH --partition=acpu      # was amilan; CURC retired that partition (2026-08-24)
#SBATCH --qos=cpu-normal      # acpu accepts only cpu-normal / cpu-long, not normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G                   # scoring every drug, not the superseded rung's 32:
                                     # the returned frame is tens of times larger, and the
                                     # nodes carry 245G, so headroom is free here and an
                                     # out-of-memory kill costs the whole 40-minute scan
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Rung-0 delta-reproducibility ceiling (docs/SPEC.md). Splits each (line, drug)'s replicate
# plates into two halves, aggregates each half's delta in DuckDB, and correlates the halves
# over the declared gene panel -- the delta's own test-retest reliability, the most any
# predictor at rung 1 could score against the real delta. CPU only; reads the Tahoe DE
# parquet already on scratch.
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

# Pinned to the pyproject floors with upper bounds so a fallback install can't drift the
# persistent stack env to an untested major version. Lock-based provisioning (containers/uv on
# Alpine) is a declared follow-up, not done here.
python -c "import duckdb, matplotlib" 2>/dev/null || pip install -q "duckdb>=1.0,<2" "matplotlib>=3.8,<4"

# --frame-cache: the DuckDB scan over every DE shard dominates this job (~40 min); the rest
# takes about a minute. Caching the built frame on scratch, keyed by a hash of its inputs, means
# a rerun that only adds an output or changes a figure skips the scan. Delete the cache dir to
# force a rebuild from the shards.
#
# NO --drugs-cid-file and NO --panel-file. Rung 0 measures the reliability of the assay at its
# full extent: every gene the table carries, every drug with plates enough to split, doses
# pooled. The superseded rung passed a 32-compound list and a 14,121-gene panel derived from
# inputs that are on no landed branch and cannot be rebuilt; passing either here would make this
# ceiling a measurement of a subset nobody at rung 0 can justify. A later rung that scores a
# subset computes its own restriction and supplies its own files.
OUT_DIR="${OUT_DIR:-docs/tasks/rung0-assay-reliability}"
python scripts/delta_reproducibility.py \
    --local-dir "/scratch/alpine/$USER/tahoe_pseudobulk_de" \
    --frame-cache "/scratch/alpine/$USER/rung0_frame_cache" \
    --duckdb-memory 110GB \
    --out-dir "$OUT_DIR"
