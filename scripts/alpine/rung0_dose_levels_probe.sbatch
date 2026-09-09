#!/bin/bash
#SBATCH --job-name=rung0-doselv
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=180G
#SBATCH --time=00:40:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Where does the within-dose replication actually sit? 7,641 (line, drug, dose) triples have at
# least two plates, but that is a union across dose levels -- it does not say whether one dose
# was systematically repeated or whether the replication is scattered. A within-dose reliability
# scored on a scattered base is a mixture over dose levels; scored on one level it is a single
# well-defined quantity. This decides which.
set -euo pipefail
cd "${REPO:-$SLURM_SUBMIT_DIR}"
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export PYTHONUNBUFFERED=1
python - <<'PY'
import duckdb, glob, os
paths = sorted(glob.glob(f"/scratch/alpine/{os.environ['USER']}/tahoe_pseudobulk_de/**/*.parquet", recursive=True))
con = duckdb.connect()
con.execute("SET temp_directory='/scratch/alpine/%s/duckdb_tmp'" % os.environ["USER"])
con.execute("SET memory_limit='150GB'"); con.execute("SET threads=8")
con.execute("SET preserve_insertion_order=false")
# One narrow projection straight into the per-key aggregate. The previous version built a
# DISTINCT temp table over 4.09 billion rows first and was OOM-killed at 64 GB; the keys are
# all that is needed, and DuckDB can aggregate them without materialising the distinct set.
con.execute("""CREATE TEMP TABLE k AS
  SELECT Cell_ID_DepMap AS line, drug, concentration AS dose, count(DISTINCT plate) AS n_plates
  FROM read_parquet(?) WHERE plate IS NOT NULL
  GROUP BY 1, 2, 3""", [paths])
print("dose-conditions:", con.execute("SELECT count(*) FROM k").fetchone()[0], flush=True)

print("\nreplicated (>=2 plates) dose-conditions, BY DOSE LEVEL:", flush=True)
rows = con.execute("""
  SELECT dose, count(*) AS dose_conditions,
         count(*) FILTER (WHERE n_plates >= 2) AS replicated,
         sum(n_plates) FILTER (WHERE n_plates >= 2) AS plates_in_replicated
  FROM k GROUP BY dose ORDER BY dose""").fetchall()
tot_rep = 0
for dose, dc, rep, plates in rows:
    tot_rep += rep or 0
    print(f"   dose={dose:<12} dose-conditions={dc:,}  replicated={rep:,} "
          f"({100*(rep or 0)/max(dc,1):.1f}%)  plates in those={plates or 0:,}", flush=True)
print(f"\ntotal replicated dose-conditions = {tot_rep:,}", flush=True)

print("\ndistinct DRUGS and LINES inside the replicated subset:", flush=True)
r = con.execute("""
  SELECT count(DISTINCT drug), count(DISTINCT line), count(*)
  FROM k WHERE n_plates >= 2""").fetchone()
print(f"   drugs={r[0]:,}  lines={r[1]:,}  (line,drug,dose) triples={r[2]:,}", flush=True)
PY
