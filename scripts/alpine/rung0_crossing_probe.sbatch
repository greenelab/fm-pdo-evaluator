#!/bin/bash
#SBATCH --job-name=rung0-cross
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
set -euo pipefail
cd "${REPO:-$SLURM_SUBMIT_DIR}"
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export PYTHONUNBUFFERED=1
python - <<'PY'
import duckdb, glob, os
paths = sorted(glob.glob(f"/scratch/alpine/{os.environ['USER']}/tahoe_pseudobulk_de/**/*.parquet", recursive=True))
con = duckdb.connect()
con.execute("SET temp_directory='/scratch/alpine/%s/duckdb_tmp'" % os.environ["USER"])
con.execute("SET memory_limit='48GB'"); con.execute("SET threads=8")
con.execute("SET preserve_insertion_order=false")
# Is dose crossed with plate, or confounded with it? One small GROUP BY over the KEYS only --
# no gene dimension, so this is cheap.
q = """
WITH k AS (SELECT DISTINCT Cell_ID_DepMap AS line, drug, concentration AS dose, plate
           FROM read_parquet(?) WHERE plate IS NOT NULL)
SELECT n_plates, count(*) AS dose_conditions FROM (
  SELECT line, drug, dose, count(DISTINCT plate) AS n_plates FROM k GROUP BY 1,2,3)
GROUP BY 1 ORDER BY 1"""
print("plates per (line, drug, DOSE):", flush=True)
tot = splittable = 0
for n, c in con.execute(q, [paths]).fetchall():
    print(f"   plates={n:<4} dose-conditions={c:,}", flush=True)
    tot += c
    if n >= 2: splittable += c
print(f"total dose-conditions={tot:,}  with >=2 plates={splittable:,} "
      f"({100*splittable/max(tot,1):.1f}%)", flush=True)
PY
