#!/bin/bash
#SBATCH --job-name=tahoe-cids
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Prep: resolve which of the 229 GDSC2 CIDs actually exist as drugs in Tahoe, and
# write ONLY those to data/static/tahoe_target_cids.txt. This makes every downstream
# array shard non-empty by construction, so 01's afterok dependency is meaningful
# (a failed shard is a REAL failure, not "this shard happened to have no Tahoe drug").
set -euo pipefail

REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"

# ---- environment (data/scoring env: fmharness + datasets + pyarrow; NO torch) -------
module load anaconda
conda activate stack
export PYTHONPATH="$REPO/src"
export HF_HUB_CACHE="/scratch/alpine/$USER/hf"    # cache large downloads on scratch, not /projects
# do NOT set HF_HOME to scratch -- that hides the token `hf auth login` saved under the HF cache:
export HF_TOKEN="${HF_TOKEN:-$(cat "${XDG_CACHE_HOME:-$HOME/.cache}/huggingface/token" 2>/dev/null || true)}"

python - <<'PY'
from pathlib import Path
from datasets import load_dataset

repo = Path.cwd()
targets = {t for t in (repo / "data/static/gdsc2_auc_pubchem_cids.txt").read_text().split() if t}
dm = load_dataset("tahoebio/Tahoe-100M", "drug_metadata", split="train").to_pandas()
dm = dm[dm["pubchem_cid"].notna()].copy()
dm["cid"] = dm["pubchem_cid"].map(lambda c: str(int(c)))
mapped = sorted({c for c in dm["cid"] if c in targets}, key=int)
if not mapped:
    raise RuntimeError("no GDSC2 target CIDs matched Tahoe drug metadata")
out = repo / "data/static/tahoe_target_cids.txt"
out.write_text("\n".join(mapped) + "\n")
print(f"{len(mapped)} of {len(targets)} GDSC2 CIDs are Tahoe drugs -> {out}")
PY
