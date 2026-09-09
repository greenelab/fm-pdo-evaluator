#!/bin/bash
#SBATCH --job-name=rung0-dosebal
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=00:40:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Dose is confounded with plate: 86.6% of (line, drug, dose) combinations sit on ONE plate
# (job 31996238). The reliabilities split plates by hash(plate) % 2 and pool dose, so if a
# condition's doses live on different plates the two halves can carry DIFFERENT DOSES -- and
# the split-half correlation would then be measuring reproducibility across a dose difference,
# not across replicates. That would depress it for a reason the design does not intend.
#
# Keys only, no gene dimension: one cheap scan.
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

con.execute("""CREATE TEMP TABLE k AS
  SELECT DISTINCT Cell_ID_DepMap AS line, drug, concentration AS dose, plate,
         hash(plate) % 2 AS half
  FROM read_parquet(?) WHERE plate IS NOT NULL""", [paths])
print("distinct (line, drug, dose, plate) keys:",
      con.execute("SELECT count(*) FROM k").fetchone()[0], flush=True)

rows = con.execute("""
  WITH per AS (
    SELECT line, drug,
           count(DISTINCT plate) AS n_plates,
           count(DISTINCT CASE WHEN half=0 THEN dose END) AS n_dose_h0,
           count(DISTINCT CASE WHEN half=1 THEN dose END) AS n_dose_h1,
           list_sort(list(DISTINCT CASE WHEN half=0 THEN dose END)) AS d0,
           list_sort(list(DISTINCT CASE WHEN half=1 THEN dose END)) AS d1
    FROM k GROUP BY line, drug)
  SELECT count(*) AS conditions,
         count(*) FILTER (WHERE n_plates >= 2) AS splittable,
         count(*) FILTER (WHERE n_plates >= 2 AND d0 = d1) AS same_dose_set,
         count(*) FILTER (WHERE n_plates >= 2 AND d0 <> d1) AS different_dose_set,
         count(*) FILTER (WHERE n_plates >= 2 AND n_dose_h0 <> n_dose_h1) AS different_dose_count
  FROM per""").fetchone()
names = ["conditions","splittable","same_dose_set","different_dose_set","different_dose_count"]
for n, v in zip(names, rows):
    print(f"{n:24s} {v:,}", flush=True)
sp = rows[1] or 1
print(f"\nof splittable conditions: {100*rows[2]/sp:.1f}% have the SAME dose set in both halves,",
      f"{100*rows[3]/sp:.1f}% differ", flush=True)
PY
