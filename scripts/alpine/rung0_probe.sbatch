#!/bin/bash
#SBATCH --job-name=rung0-probe
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=00:45:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Scale probe: how big are the frames rung 0 now asks DuckDB to return, with no drug filter?
# The superseded run scored 32 of ~1,100 drugs; this task scores every splittable drug, so the
# returned (line, drug, gene) frame is tens of times larger and the noise decomposition adds a
# second aggregation grouped by dose. Both come back to pandas, so a frame that does not fit
# kills the job after the scan that dominates its runtime.
#
# ONE scan, approximate distinct counts. An exact COUNT(DISTINCT) over 4.09 billion rows is a
# full scan plus a hash set per query, and the first version of this probe spent five minutes
# per column and would have hit its wall clock before answering. HyperLogLog is within about
# two percent, which is far inside the precision this decision needs -- the question is whether
# a frame is 20 million rows or 500 million.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1

python - <<'PY'
import duckdb, glob, os, time
paths = sorted(glob.glob(f"/scratch/alpine/{os.environ['USER']}/tahoe_pseudobulk_de/**/*.parquet", recursive=True))
print("shards:", len(paths), flush=True)
con = duckdb.connect()
con.execute("SET temp_directory='/scratch/alpine/%s/duckdb_tmp'" % os.environ["USER"])
con.execute("SET memory_limit='72GB'")
con.execute("SET preserve_insertion_order=false")
con.execute("SET threads=16")

t = time.time()
row = con.execute("""
  SELECT approx_count_distinct(gene_name)                                   AS genes,
         approx_count_distinct(Cell_ID_DepMap)                              AS lines,
         approx_count_distinct(drug)                                        AS drugs,
         approx_count_distinct(plate)                                       AS plates,
         approx_count_distinct(concentration)                               AS doses,
         approx_count_distinct(Cell_ID_DepMap || '|' || drug)               AS conditions,
         approx_count_distinct(Cell_ID_DepMap || '|' || drug || '|' || gene_name) AS split_half_rows,
         approx_count_distinct(Cell_ID_DepMap || '|' || drug || '|' || CAST(concentration AS VARCHAR) || '|' || gene_name) AS noise_rows,
         count(*) FILTER (WHERE plate IS NULL)                              AS rows_no_plate,
         count(*) FILTER (WHERE padj IS NULL)                               AS rows_null_padj,
         count(*) FILTER (WHERE baseMean = 0)                               AS rows_zero_basemean
  FROM read_parquet(?)""", [paths]).fetchone()
names = ["genes","lines","drugs","plates","doses","conditions","split_half_rows","noise_rows",
         "rows_no_plate","rows_null_padj","rows_zero_basemean"]
for n, v in zip(names, row):
    print(f"{n:22s} {v:,}", flush=True)
print(f"one-scan elapsed: {time.time()-t:.0f}s", flush=True)

# Plates per condition: a small GROUP BY, cheap relative to the scan above.
t = time.time()
hist = con.execute("""
  SELECT n_plates, count(*) AS conditions FROM (
    SELECT Cell_ID_DepMap, drug, count(DISTINCT plate) AS n_plates
    FROM read_parquet(?) WHERE plate IS NOT NULL GROUP BY 1,2)
  GROUP BY n_plates ORDER BY n_plates""", [paths]).fetchall()
print("plates per condition:", flush=True)
total = splittable = even = 0
for n_plates, conditions in hist:
    print(f"   plates={n_plates:<4} conditions={conditions:,}", flush=True)
    total += conditions
    if n_plates >= 2: splittable += conditions
    if n_plates >= 2 and n_plates % 2 == 0: even += conditions
print(f"conditions total={total:,} splittable={splittable:,} equal-half={even:,}", flush=True)
print(f"histogram elapsed: {time.time()-t:.0f}s", flush=True)
PY
