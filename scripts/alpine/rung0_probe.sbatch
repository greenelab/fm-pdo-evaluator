#!/bin/bash
#SBATCH --job-name=rung0-probe
#SBATCH --partition=acpu
#SBATCH --qos=cpu-normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=00:40:00
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.out
#
# Scale probe: how big are the frames rung 0 now asks DuckDB to return, with no drug filter?
# The superseded run scored 32 of ~1,100 drugs; this task scores every splittable drug, which
# multiplies the returned (line, drug, gene) frame. Measured before the real job is submitted,
# because a job that dies on memory after 40 minutes of scanning costs the same 40 minutes.
set -euo pipefail
REPO="${REPO:-$SLURM_SUBMIT_DIR}"
cd "$REPO"
export PATH="/projects/$USER/software/anaconda/envs/stack/bin:$PATH"
export PYTHONPATH="$REPO/src"
export PYTHONUNBUFFERED=1

python - <<'PY'
import duckdb, glob, os
paths = sorted(glob.glob(f"/scratch/alpine/{os.environ['USER']}/tahoe_pseudobulk_de/**/*.parquet", recursive=True))
print("shards:", len(paths))
con = duckdb.connect()
con.execute("SET temp_directory='/scratch/alpine/%s/duckdb_tmp'" % os.environ["USER"])
con.execute("SET memory_limit='36GB'")
con.execute("SET preserve_insertion_order=false")
q = lambda s: con.execute(s, [paths]).fetchall()
print("cols:", [r[0] for r in con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [paths]).fetchall()])
print("total rows:", q("SELECT count(*) FROM read_parquet(?)")[0][0])
print("distinct genes:", q("SELECT count(DISTINCT gene_name) FROM read_parquet(?)")[0][0])
print("distinct lines:", q("SELECT count(DISTINCT Cell_ID_DepMap) FROM read_parquet(?)")[0][0])
print("distinct drugs:", q("SELECT count(DISTINCT drug) FROM read_parquet(?)")[0][0])
print("distinct plates:", q("SELECT count(DISTINCT plate) FROM read_parquet(?)")[0][0])
print("distinct doses:", q("SELECT count(DISTINCT concentration) FROM read_parquet(?)")[0][0])
print("(line,drug) conditions:", q("SELECT count(*) FROM (SELECT DISTINCT Cell_ID_DepMap, drug FROM read_parquet(?))")[0][0])
print("(line,drug) with >=2 plates:", q("SELECT count(*) FROM (SELECT Cell_ID_DepMap, drug FROM read_parquet(?) WHERE plate IS NOT NULL GROUP BY 1,2 HAVING count(DISTINCT plate)>=2)")[0][0])
print("split-half frame rows (line,drug,gene):", q("SELECT count(*) FROM (SELECT DISTINCT Cell_ID_DepMap, drug, gene_name FROM read_parquet(?))")[0][0])
print("noise frame rows (line,drug,dose,gene) with >=2 plates:", q("SELECT count(*) FROM (SELECT Cell_ID_DepMap, drug, concentration, gene_name FROM read_parquet(?) WHERE plate IS NOT NULL GROUP BY 1,2,3,4 HAVING count(DISTINCT plate)>=2)")[0][0])
print("plates per (line,drug) histogram:")
for row in q("SELECT n, count(*) FROM (SELECT Cell_ID_DepMap, drug, count(DISTINCT plate) AS n FROM read_parquet(?) WHERE plate IS NOT NULL GROUP BY 1,2) GROUP BY n ORDER BY n"):
    print("   plates=%s conditions=%s" % row)
PY
