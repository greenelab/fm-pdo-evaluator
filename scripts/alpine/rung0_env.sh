# Sourced by every rung 0 job script: the environment the three stages share.
#
# Activate the `stack` env by PATH instead of the module system. `module load anaconda` fails
# inside batch jobs (Lmod: module unknown) for reasons NOT established -- see ralpine's submit
# comment; do not assert a cause here. This pattern is the one that is PROVEN to work.
# `conda activate` is deliberately not used: this install ships no etc/profile.d/conda.sh.
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export CONDA_PREFIX="/projects/$USER/software/anaconda/envs/stack"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1   # else stdout block-buffers to the log file and a multi-hour job
                           # shows 0 bytes until it exits (job 31655278, 2026-08-25).
export HF_HUB_CACHE="/scratch/alpine/$USER/hf"
# do NOT set HF_HOME to scratch -- that hides the token `hf auth login` saved under the HF cache:
export HF_TOKEN="${HF_TOKEN:-$(cat "${XDG_CACHE_HOME:-$HOME/.cache}/huggingface/token" 2>/dev/null || true)}"

# Pinned to the pyproject floors with upper bounds so a fallback install can't drift the
# persistent stack env to an untested major version.
python -c "import duckdb, matplotlib" 2>/dev/null || pip install -q "duckdb>=1.0,<2" "matplotlib>=3.8,<4"

# The three stages share one cache and one slice count. The cache lives on scratch: it is a
# speed optimisation, not the reproduction chain -- the pinned tranche plus the committed code
# is that -- and scratch is purged. The slice count is the job array's size.
export RUNG0_CACHE="${RUNG0_CACHE:-/scratch/alpine/$USER/rung0_cache_v3}"
export RUNG0_SLICES="${RUNG0_SLICES:-32}"
export OUT_DIR="${OUT_DIR:-docs/tasks/rung0-assay-reliability}"
mkdir -p logs
