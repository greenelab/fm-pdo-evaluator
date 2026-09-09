"""Delta reproducibility ceiling: how noisy is the real Tahoe delta itself?

This is rung 0 of the ladder in docs/SPEC.md: a reference every other rung reads its score
against. Split each (line, drug, dose) triple's replicate plates into two halves, aggregate
each half's delta, and correlate the two halves over every gene. That split-half
correlation is the delta's own test-retest reliability -- the most any predictor at rung 1
could score against the real delta. Spearman-Brown then lifts the half-data number to the
full-data reliability rung 1 is read against.

Reuses the DuckDB-over-local-parquet path (the Tahoe DE table is already local or on
scratch), so it needs no GPU.

  python scripts/delta_reproducibility.py --local-dir /scratch/alpine/$USER/tahoe_pseudobulk_de \\
      --drug-names-file <a file of Tahoe drug names, one per line> \\
      --panel-file results/rung1_panel/common_panel.txt --out-dir rung0_outputs
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

#: The columns that identify one scoreable unit. Dose is part of it: this screen ran 86.6% of
#: (line, drug, dose) combinations on a single plate, so splitting a condition's plates while
#: pooling dose puts different doses in the two halves for 99.7% of conditions -- a
#: dose-to-dose correlation, not a test-retest reliability. See decisions.md, 2026-09-01.
CONDITION_KEYS = ("patient", "drug", "dose")

#: Bumped whenever the built frame's schema, grouping or split rule changes, so a cached frame
#: from an earlier definition resolves to a different key instead of being silently reused.
#: v2 held dose fixed but split plates by ``hash(plate) % 2``, which puts a two-plate triple on
#: one side of the split about half the time; v3 alternates over the sorted plate ids, so every
#: replicated triple has a plate on each side. Same inputs, different meaning, different key.
FRAME_SCHEMA = "v3-dose-fixed-alternating"

#: How plates are assigned to halves, in one sentence, recorded beside every artifact that
#: depends on it. Within each (line, drug, dose) triple the distinct plate ids are sorted and
#: assigned alternately -- first to half 0, second to half 1, third to half 0, and so on. This
#: is deterministic, needs no random seed, gives every replicated triple a plate in each half,
#: and makes "equal halves" exactly "an even plate count".
SPLIT_RULE = "plates sorted by id within each (line, drug, dose) triple, assigned alternately"

DOSE_CANDIDATES = ("dose", "Dose", "drug_dose", "concentration", "dose_uM")

TAHOE = "tahoebio/Tahoe-100M"
DE = "pseudobulk_differential_expression"
REPL_CANDIDATES = ("plate", "Plate", "plate_barcode", "plate_id", "replicate", "batch")


def _target_names(repo: Path, cid_file: Path) -> list[str]:
    """Tahoe drug NAMES for the target PubChem CIDs (drug_metadata), matching the shortcut."""
    from datasets import load_dataset  # type: ignore  # Alpine-only

    cids = {t for t in cid_file.read_text().split() if t}
    dm = load_dataset(TAHOE, "drug_metadata", split="train").to_pandas()
    dm = dm[dm["pubchem_cid"].notna()].copy()
    dm["cid"] = dm["pubchem_cid"].map(lambda c: str(int(c)))
    return sorted(dm[dm["cid"].isin(cids)]["drug"].astype(str).unique())


def resolve_drug_names(repo: Path, args: argparse.Namespace) -> list[str] | None:
    """The drug list a run scores, or ``None`` meaning every drug in the pool.

    Rung 0 measures at the assay's full extent, so no drug file is the expected case and
    ``None`` is the expected answer. A file is still accepted -- a later rung asking for a
    restriction needs one -- but a missing default file is not an error here, because the
    superseded rung's compound list is not on any branch and cannot be rebuilt.
    """
    if getattr(args, "drug_names_file", None):
        return sorted(
            {ln.strip() for ln in Path(args.drug_names_file).read_text().splitlines() if ln.strip()}
        )
    # Checked as a string before it becomes a Path: Path("") is PosixPath("."), whose str() is
    # "." and is truthy, so testing the Path would send an empty default down the lookup path.
    raw = str(getattr(args, "drugs_cid_file", "") or "").strip()
    if not raw:
        return None
    cid_file = Path(raw)
    cid_file = cid_file if cid_file.is_absolute() else repo / cid_file
    if not cid_file.exists():
        return None
    return _target_names(repo, cid_file)


def _connect(tmp: Path, memory_limit: str = "36GB", threads: int | None = None):
    """A DuckDB connection configured to spill to ``tmp`` rather than exhaust memory.

    ``threads`` is worth setting low for the wide group-bys. DuckDB builds partial hash tables
    per thread, so sixteen threads over a billion groups multiplies the peak by roughly that
    factor before any of it can spill -- and the failure arrives as an out-of-memory error
    rather than as slow spilling.
    """
    import duckdb  # type: ignore  # Alpine-only

    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET temp_directory='{tmp}'")
    con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET preserve_insertion_order=false")
    if threads:
        con.execute(f"SET threads={int(threads)}")
    return con


def _gene_partition(n_parts: int, part: int, alias: str = "") -> str:
    """A SQL predicate selecting one slice of the genes, or nothing when there is one slice.

    Partitioning by ``hash(gene_name)`` rather than by file or by row is what makes the slices
    safe to aggregate independently: ``gene_name`` is part of every group key the noise
    decomposition uses, so a group's rows all land in the same slice and per-slice sums combine
    into the whole exactly. Splitting by shard would not have that property -- a gene's rows are
    spread across shards, and per-shard variances would each be computed on a fragment.
    """
    return "" if n_parts <= 1 else f" AND hash({alias}gene_name) % {int(n_parts)} = {int(part)}"


KEY_COLUMNS = ("patient", "drug", "gene_name", "dose")


def _compact_df(result) -> pd.DataFrame:
    """Fetch a DuckDB result as pandas with the key columns dictionary-encoded.

    Rung 0 scores every drug in the screen, not the superseded rung's 32, so the aggregated
    frame is tens of times larger and its three key columns repeat a few tens of thousands of
    distinct strings across hundreds of millions of rows. As Python objects those strings, not
    the fold changes, are what exhausts the job's memory. Arrow's dictionary encoding stores
    each distinct value once and an integer code per row, and pandas reads a dictionary array
    back as a Categorical, so the frame that reaches the scorer carries the same values at a
    fraction of the footprint. Nothing downstream changes: a Categorical indexes, groups and
    pivots exactly as the object column did.
    """
    import pyarrow as pa  # type: ignore  # Alpine-only

    # `.arrow()` returns a RecordBatchReader on some DuckDB versions and a Table on others;
    # `fetch_arrow_table()` is the one that is a Table everywhere, with the reader read out as
    # the fallback. Read off the installed version rather than assumed -- the first version of
    # this helper called `.arrow()` and failed on `column_names`.
    for attr in ("to_arrow_table", "fetch_arrow_table", "arrow"):
        if hasattr(result, attr):
            tbl = getattr(result, attr)()
            break
    else:  # pragma: no cover - every supported DuckDB exposes one of the three
        raise RuntimeError("this DuckDB build exposes no Arrow accessor")
    if hasattr(tbl, "read_all"):  # a RecordBatchReader on some versions
        tbl = tbl.read_all()
    cols = [
        tbl.column(i).dictionary_encode() if name in KEY_COLUMNS else tbl.column(i)
        for i, name in enumerate(tbl.column_names)
    ]
    return pa.Table.from_arrays(cols, names=tbl.column_names).to_pandas()


def _drug_predicate(target_names: list[str] | None, alias: str = "") -> tuple[str, list[object]]:
    """The drug filter, and the parameters it needs — empty when every drug is admitted.

    Rung 0 measures at the assay's full extent: every drug with plates enough to split. The
    superseded rung passed a 32-compound list derived from inputs no longer on any branch, so
    "no drug file" has to mean *no predicate*, not a predicate matching nothing.
    """
    if not target_names:
        return "", []
    return f" AND {alias}drug IN (SELECT unnest(?))", [list(target_names)]


def _pool_columns(con, paths: list[str], repl_col: str | None) -> tuple[list[str], str, str]:
    """The pool's column names, and the replicate and dose columns resolved against them."""
    schema = con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [paths]).df()
    cols = list(schema["column_name"])
    candidates = ([repl_col] if repl_col else []) + list(REPL_CANDIDATES)
    repl = next((c for c in candidates if c in cols), None)
    if repl is None:
        raise SystemExit(f"no replicate column in {cols}; pass --replicate-col")
    dose = next((c for c in DOSE_CANDIDATES if c in cols), None)
    if dose is None:
        raise SystemExit(
            f"no dose column in {cols}. Dose is part of the condition key: on this screen it is "
            "confounded with plate, so pooling it would compare different doses across the two "
            "halves. Pass a pool that carries one, or change the design."
        )
    return cols, repl, dose


def split_assignment(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    threads: int | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, str]:
    """Which half every plate of every (line, drug, dose) triple belongs to, and the pool's shape.

    One scan over the key columns alone -- about a quarter of an hour on the real screen -- and
    the split is then an explicit, committed table rather than a hash function evaluated inside
    the aggregations. ``SPLIT_RULE`` says how: plates sorted by id within a triple, assigned
    alternately. Under that rule a triple with two plates always has one on each side, which
    the hash split gave only when the two plate ids happened to hash to different parities.

    Returns ``(assignment, pool, replicate_col)``. ``assignment`` has one row per (triple, plate)
    with its ``half``; ``pool`` has one row per triple with the plate counts per half, the
    equal-halves flag (``n_plates_half0 == n_plates_half1``, which is the exact condition
    Spearman-Brown's correction assumes), and the share of rows DESeq2 could not test.
    """
    con = _connect(tmp, memory_limit, threads)
    cols, repl, dose = _pool_columns(con, paths, repl_col)
    untestable = (
        "avg(CASE WHEN baseMean = 0 THEN 1.0 ELSE 0.0 END)" if "baseMean" in cols else "NULL"
    )
    where, drug_params = _drug_predicate(target_names)
    con.execute(
        f"""CREATE OR REPLACE TEMP TABLE plate_keys AS
            SELECT Cell_ID_DepMap AS patient, drug, {dose} AS dose, {repl} AS plate,
                   count(*) AS n_rows, {untestable} AS frac_untestable
            FROM read_parquet(?)
            WHERE {repl} IS NOT NULL{where}
            GROUP BY 1, 2, 3, 4""",
        [paths, *drug_params],
    )
    assignment = con.execute(
        """SELECT patient, drug, dose, plate,
                  (row_number() OVER (PARTITION BY patient, drug, dose ORDER BY plate) - 1) % 2
                    AS half,
                  n_rows, frac_untestable
           FROM plate_keys ORDER BY patient, drug, dose, plate"""
    ).df()
    pool = con.execute(
        """SELECT patient, drug, dose,
                  sum(n_rows) AS n_rows,
                  count(*) AS n_plates,
                  count(*) FILTER (WHERE half = 0) AS n_plates_half0,
                  count(*) FILTER (WHERE half = 1) AS n_plates_half1,
                  count(*) FILTER (WHERE half = 0) = count(*) FILTER (WHERE half = 1)
                    AS n_plates_even,
                  sum(frac_untestable * n_rows) / sum(n_rows) AS frac_untestable
           FROM (SELECT *,
                        (row_number() OVER (PARTITION BY patient, drug, dose ORDER BY plate) - 1)
                          % 2 AS half
                 FROM plate_keys)
           GROUP BY patient, drug, dose ORDER BY patient, drug, dose"""
    ).df()
    assignment["half"] = assignment["half"].astype(int)
    pool["n_rows"] = pool["n_rows"].astype(int)
    for c in ("n_plates", "n_plates_half0", "n_plates_half1"):
        pool[c] = pool[c].astype(int)
    pool["n_plates_even"] = pool["n_plates_even"].astype(bool)
    return assignment, pool, repl


ASSIGNMENT_FILE = "assignment.parquet"
POOL_FILE = "pool.parquet"
ASSIGNMENT_SIDECAR = "assignment.json"


def write_assignment(
    cache_dir: Path, assignment: pd.DataFrame, pool: pd.DataFrame, repl: str
) -> Path:
    """Cache the split so every slice joins against the same table."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    assignment.to_parquet(cache_dir / ASSIGNMENT_FILE, index=False)
    pool.to_parquet(cache_dir / POOL_FILE, index=False)
    (cache_dir / ASSIGNMENT_SIDECAR).write_text(
        json.dumps({"replicate_col": repl, "split_rule": SPLIT_RULE, "schema": FRAME_SCHEMA}) + "\n"
    )
    return cache_dir / ASSIGNMENT_FILE


def read_assignment(cache_dir: Path) -> tuple[Path, pd.DataFrame, str]:
    """The cached split: its parquet path for joins, the pool table, and the replicate column."""
    side = cache_dir / ASSIGNMENT_SIDECAR
    if not (cache_dir / ASSIGNMENT_FILE).exists() or not side.exists():
        raise SystemExit(
            f"no split assignment in {cache_dir}: run the assign stage first "
            "(--stage assign), so every slice joins the same split"
        )
    meta = json.loads(side.read_text())
    if meta.get("schema") != FRAME_SCHEMA:
        raise SystemExit(
            f"the split in {cache_dir} was made under schema {meta.get('schema')!r}, this code "
            f"is {FRAME_SCHEMA!r}; rerun the assign stage rather than mixing the two"
        )
    return (
        cache_dir / ASSIGNMENT_FILE,
        pd.read_parquet(cache_dir / POOL_FILE),
        str(meta["replicate_col"]),
    )


def _ensure_assignment(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str,
    threads: int | None,
    assignment: Path | None,
) -> Path:
    """The assignment parquet to join against, building it under ``tmp`` when none is given."""
    if assignment is not None:
        return assignment
    a, pool, repl = split_assignment(paths, target_names, repl_col, tmp, memory_limit, threads)
    return write_assignment(tmp / "split", a, pool, repl)


#: What one slice's aggregate carries, per (line, drug, dose, gene). The two consumers read
#: different subsets: the reliabilities take the half means and the first group's p-values,
#: the decomposition takes the variance across plates and the mean squared standard error.
SLICE_FRAME_COLUMNS = ("patient", "drug", "dose", "gene_name", "lfc0", "lfc1", "padj0", "padj1")
SLICE_NOISE_COLUMNS = (
    "patient",
    "drug",
    "dose",
    "gene_name",
    "var_lfc",
    "mean_se2",
    "n_plates",
    "base_mean",
    "mean_lfc",
    "padj0",
)


def slice_aggregate(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    *,
    assignment: Path,
    n_parts: int = 1,
    part: int = 0,
    threads: int | None = None,
) -> tuple[pd.DataFrame, str]:
    """One slice of the genes, aggregated once for both the reliabilities and the decomposition.

    Reading the table is the entire cost -- every slice scans all four billion rows and keeps
    the genes whose hash lands in it -- so the split-half means and the across-plate variance
    come out of ONE group-by rather than two scans. Both group by exactly (line, drug, dose,
    gene), so nothing is lost by sharing the pass.

    Plates reach their halves through the cached assignment (``split_assignment``), joined on
    the full plate key. The join is against a table of a few tens of thousands of rows and costs
    nothing against the scan.

    ``padj0`` is the MINIMUM adjusted p-value over the FIRST half's rows -- the selection rule's
    own aggregate, one-sided so that selecting on the half a correlation is scored against
    cannot inflate it. ``padj1`` exists for the overlap diagnostic alone.

    Groups are returned when they have a fold change in both halves (a scoreable gene) or at
    least two distinct plates (a decomposable gene); the callers take their own subset with
    ``frame_from_slice`` and ``noise_from_slice``.
    """
    con = _connect(tmp, memory_limit, threads)
    cols, repl, dose = _pool_columns(con, paths, repl_col)
    padj0 = "min(t.padj) FILTER (WHERE h.half = 0)" if "padj" in cols else "CAST(NULL AS DOUBLE)"
    padj1 = "min(t.padj) FILTER (WHERE h.half = 1)" if "padj" in cols else "CAST(NULL AS DOUBLE)"
    se2 = "avg(t.lfcSE * t.lfcSE)" if "lfcSE" in cols else "CAST(NULL AS DOUBLE)"
    base = "avg(t.baseMean)" if "baseMean" in cols else "CAST(NULL AS DOUBLE)"
    if "padj" not in cols:
        print("no padj column in the pool: responder selection is unavailable for this run")
    if "lfcSE" not in cols:
        print("no lfcSE column in the pool: the noise decomposition is unavailable for this run")
    where, drug_params = _drug_predicate(target_names, alias="t.")
    agg = con.execute(
        f"""SELECT t.Cell_ID_DepMap AS patient, t.drug, t.{dose} AS dose, t.gene_name,
                   avg(t.log2FoldChange) FILTER (WHERE h.half = 0) AS lfc0,
                   avg(t.log2FoldChange) FILTER (WHERE h.half = 1) AS lfc1,
                   count(t.log2FoldChange) FILTER (WHERE h.half = 0) AS n0,
                   count(t.log2FoldChange) FILTER (WHERE h.half = 1) AS n1,
                   {padj0} AS padj0,
                   {padj1} AS padj1,
                   var_samp(t.log2FoldChange) AS var_lfc,
                   {se2} AS mean_se2,
                   count(DISTINCT t.{repl}) AS n_plates,
                   {base} AS base_mean,
                   avg(t.log2FoldChange) AS mean_lfc
            FROM read_parquet(?) t
            JOIN read_parquet(?) h
              ON t.Cell_ID_DepMap = h.patient AND t.drug = h.drug
             AND t.{dose} = h.dose AND t.{repl} = h.plate
            WHERE t.{repl} IS NOT NULL{where}{_gene_partition(n_parts, part, alias="t.")}
            GROUP BY 1, 2, 3, 4
            HAVING (count(t.log2FoldChange) FILTER (WHERE h.half = 0) > 0
                    AND count(t.log2FoldChange) FILTER (WHERE h.half = 1) > 0)
                OR count(DISTINCT t.{repl}) >= 2""",
        [paths, str(assignment), *drug_params],
    )
    return _compact_df(agg), repl


def frame_from_slice(agg: pd.DataFrame) -> pd.DataFrame:
    """The split-half frame: gene-conditions with a fold change in BOTH halves."""
    both = (agg["n0"].to_numpy() > 0) & (agg["n1"].to_numpy() > 0)
    return agg.loc[both, list(SLICE_FRAME_COLUMNS)].reset_index(drop=True)


def noise_from_slice(agg: pd.DataFrame) -> pd.DataFrame:
    """The decomposition's rows: gene-conditions with at least two plates, decomposed."""
    two = agg["n_plates"].to_numpy() >= 2
    return decompose_noise(agg.loc[two, list(SLICE_NOISE_COLUMNS)].reset_index(drop=True))


def build_split_half_frame(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    *,
    assignment: Path | None = None,
    n_parts: int = 1,
    threads: int | None = None,
) -> tuple[pd.DataFrame, str]:
    """Per (line, drug, dose, gene), the mean log2FoldChange in each plate half, via DuckDB.

    The in-process form: every slice is built here, one after another, and concatenated. The
    cluster runs the slices as a job array and concatenates them in the combine stage; the
    arithmetic is the same, which ``test_partitioned_frame_equals_one_pass`` pins.

    Only groups with a fold change in BOTH halves are returned. That is what every caller's
    ``dropna(subset=["lfc0", "lfc1"])`` does next, moved into the engine because at full extent
    it is not a small filter: DESeq2 could not test 59 percent of the screen's rows.
    """
    assignment = _ensure_assignment(
        paths, target_names, repl_col, tmp, memory_limit, threads, assignment
    )
    frames: list[pd.DataFrame] = []
    repl = ""
    for part in range(max(1, n_parts)):
        agg, repl = slice_aggregate(
            paths,
            target_names,
            repl_col,
            tmp,
            memory_limit,
            assignment=assignment,
            n_parts=n_parts,
            part=part,
            threads=threads,
        )
        frames.append(frame_from_slice(agg))
    de = pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]
    return _normalise_keys(de), repl


def build_noise_frame(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    sample_rows: int | None = None,
    n_parts: int = 1,
    threads: int | None = None,
    *,
    assignment: Path | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """Per (line, drug, dose, gene) with at least two plates: the decomposition's rows.

    ``lfcSE`` is the standard error of ONE plate's treated-versus-control contrast -- cell
    sampling at that row's cell counts. It cannot see plate-to-plate variation, which is the
    noise a model trained on other material actually meets and the noise the split-half
    measures. Under plate offsets plus independent sampling error, the sample variance of
    ``log2FoldChange`` across a condition's plates has expectation ``sigma2_plate +
    mean(lfcSE^2)`` exactly, for any set of per-plate standard errors.

    Dose is a grouping key, never pooled: a dose effect pooled into the across-plate variance
    would be reported as plate noise. With ``sample_rows`` the result is a bounded sample, an
    equal share from each gene slice.
    """
    assignment = _ensure_assignment(
        paths, target_names, repl_col, tmp, memory_limit, threads, assignment
    )
    parts = max(1, n_parts)
    per_part = None if sample_rows is None else max(1, int(sample_rows) // parts)
    frames: list[pd.DataFrame] = []
    for part in range(parts):
        agg, _ = slice_aggregate(
            paths,
            target_names,
            repl_col,
            tmp,
            memory_limit,
            assignment=assignment,
            n_parts=parts,
            part=part,
            threads=threads,
        )
        noise = noise_from_slice(agg)
        if per_part is not None and len(noise) > per_part:
            noise = noise.sample(per_part, random_state=seed).reset_index(drop=True)
        frames.append(noise)
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]


def decompose_noise(noise: pd.DataFrame) -> pd.DataFrame:
    """Each gene-condition's signed estimate of its between-plate variance.

        sigma2_plate_signed = var_samp(log2FoldChange across plates) - mean(lfcSE^2)

    Signed, not floored. With two plates the variance across them has one degree of freedom,
    so this per-gene difference is a very noisy estimate that lands either side of zero even
    with no plate effect at all. Flooring each gene at zero before averaging turns that noise
    into a positive bias: at two plates and no plate effect the floored mean is 0.48 of the
    within-plate variance, and a third of genes report a plate component that is not there.
    The estimate the run reports is therefore POOLED first and floored once
    (``pooled_plate_variance``); the per-gene column is kept signed so the figures can show
    the distribution honestly, straddling zero where there is nothing to find.

    ``between_plate_fraction_signed`` is the same difference as a share of the variance across
    plates, null where that variance is zero or missing.
    """
    out = noise.copy()
    var = out["var_lfc"].to_numpy(dtype=float)
    se2 = out["mean_se2"].to_numpy(dtype=float)
    diff = var - se2
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(var > 0, diff / var, np.nan)
    out["sigma2_plate_signed"] = diff
    out["between_plate_fraction_signed"] = frac
    return out


def pooled_plate_variance(var_lfc: np.ndarray, mean_se2: np.ndarray) -> tuple[float, float, int]:
    """The between-plate variance and its share, pooled over gene-conditions, floored ONCE.

        sigma2_plate = max(mean(var_lfc) - mean(lfcSE^2), 0)
        fraction     = sigma2_plate / mean(var_lfc)

    The mean of the signed per-gene differences is an unbiased estimate of the mean plate
    variance whatever the plate count, because expectation is linear and the floor has not
    been applied yet. Flooring the pooled value is then a constraint on one well-estimated
    number rather than on tens of millions of one-degree-of-freedom ones. Returns
    ``(sigma2_plate, fraction, n)`` over the rows where both inputs are finite and the
    variance across plates is positive; ``(nan, nan, 0)`` when there are none.
    """
    var = np.asarray(var_lfc, dtype=float)
    se2 = np.asarray(mean_se2, dtype=float)
    ok = np.isfinite(var) & (var > 0) & np.isfinite(se2)
    if not ok.any():
        return float("nan"), float("nan"), 0
    mean_var = float(np.mean(var[ok]))
    sigma2 = max(mean_var - float(np.mean(se2[ok])), 0.0)
    return sigma2, sigma2 / mean_var, int(ok.sum())


#: The per-condition sums one slice contributes. Sums, because sums add across slices; the
#: means are taken once, at the end, in ``combine_noise_partials``. The responder and
#: non-responder sums carry the squared standard errors too, so the pooled estimator can be
#: formed for each gene set and not only for all genes together.
NOISE_SUM_COLUMNS = (
    "n_gene_doses",
    "s_var",
    "s_se2",
    "n_responder_gene_doses",
    "s_var_resp",
    "s_se2_resp",
    "n_nonresponder_gene_doses",
    "s_var_non",
    "s_se2_non",
)


def noise_partials(noise: pd.DataFrame, alpha: float) -> pd.DataFrame:
    """Everything one slice contributes to the decomposition: per-condition sums.

    The slice key is part of the group key, so the slices partition the gene-conditions
    exactly -- no overlap, no omission -- and adding them gives the same answer a single pass
    would. A gene-condition contributes when its variance across plates is positive and finite
    and its squared standard error is finite; a responder is a gene the FIRST half called
    differentially expressed, the selection rule's own reading.
    """
    var = noise["var_lfc"].to_numpy(dtype=float)
    se2 = noise["mean_se2"].to_numpy(dtype=float)
    ok = np.isfinite(var) & (var > 0) & np.isfinite(se2)
    d = noise.loc[ok]
    var, se2 = var[ok], se2[ok]
    padj0 = d["padj0"].to_numpy(dtype=float)
    is_resp = np.isfinite(padj0) & (padj0 < alpha)
    per = (
        pd.DataFrame(
            {
                "patient": d["patient"].to_numpy(),
                "drug": d["drug"].to_numpy(),
                "dose": d["dose"].to_numpy(),
                "n_gene_doses": 1,
                "s_var": var,
                "s_se2": se2,
                "n_responder_gene_doses": is_resp.astype(int),
                "s_var_resp": np.where(is_resp, var, 0.0),
                "s_se2_resp": np.where(is_resp, se2, 0.0),
                "n_nonresponder_gene_doses": (~is_resp).astype(int),
                "s_var_non": np.where(is_resp, 0.0, var),
                "s_se2_non": np.where(is_resp, 0.0, se2),
            }
        )
        .groupby(list(CONDITION_KEYS), observed=True, sort=True, dropna=False)
        .sum(min_count=0)
    )
    return per.reset_index()


def _pooled_from_sums(n: np.ndarray, s_var: np.ndarray, s_se2: np.ndarray):
    """Pooled variance, plate component and share from sums, elementwise; nan where n is 0."""
    with np.errstate(invalid="ignore", divide="ignore"):
        var_mean = np.where(n > 0, s_var / n, np.nan)
        se2_mean = np.where(n > 0, s_se2 / n, np.nan)
        sigma2 = np.maximum(var_mean - se2_mean, 0.0)
        frac = np.where(var_mean > 0, sigma2 / var_mean, np.nan)
    return var_mean, se2_mean, sigma2, frac


def combine_noise_partials(per_condition: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Turn accumulated sums into the reported decomposition, once, at the end.

    Two weightings are reported, as the design asks: pooled over every gene-condition (a
    condition weighs as many of its genes as DESeq2 could test), and the mean over conditions
    of each condition's own pooled share (every condition weighs one). Where the two differ, the
    gap is uneven gene coverage rather than anything about plates.
    """
    c = (
        per_condition.groupby(list(CONDITION_KEYS), observed=True, sort=True, dropna=False)
        .sum()
        .reset_index()
    )
    by_condition = pd.DataFrame({k: c[k] for k in CONDITION_KEYS})
    by_condition["n_gene_doses"] = c["n_gene_doses"].astype(int)
    for suffix, n_col, v_col, s_col in (
        ("", "n_gene_doses", "s_var", "s_se2"),
        ("_responders", "n_responder_gene_doses", "s_var_resp", "s_se2_resp"),
        ("_nonresponders", "n_nonresponder_gene_doses", "s_var_non", "s_se2_non"),
    ):
        n = c[n_col].to_numpy(dtype=float)
        var_mean, se2_mean, sigma2, frac = _pooled_from_sums(
            n, c[v_col].to_numpy(dtype=float), c[s_col].to_numpy(dtype=float)
        )
        if suffix:
            by_condition[n_col] = c[n_col].astype(int)
        by_condition[f"var_lfc_mean{suffix}"] = var_mean
        by_condition[f"mean_se2_mean{suffix}"] = se2_mean
        by_condition[f"sigma2_plate_pooled{suffix}"] = sigma2
        by_condition[f"between_plate_fraction_pooled{suffix}"] = frac

    totals = {k: float(c[k].sum()) for k in NOISE_SUM_COLUMNS}
    overall: dict[str, float | int] = {"n_gene_conditions": int(totals["n_gene_doses"])}
    for suffix, n_key, v_key, s_key in (
        ("", "n_gene_doses", "s_var", "s_se2"),
        ("_responders", "n_responder_gene_doses", "s_var_resp", "s_se2_resp"),
        ("_nonresponders", "n_nonresponder_gene_doses", "s_var_non", "s_se2_non"),
    ):
        n = np.array([totals[n_key]])
        var_mean, se2_mean, sigma2, frac = _pooled_from_sums(
            n, np.array([totals[v_key]]), np.array([totals[s_key]])
        )
        if suffix:
            overall[f"n_gene_conditions{suffix}"] = int(totals[n_key])
        overall[f"var_lfc_mean{suffix}"] = float(var_mean[0])
        overall[f"mean_se2_mean{suffix}"] = float(se2_mean[0])
        overall[f"sigma2_plate_pooled{suffix}"] = float(sigma2[0])
        overall[f"between_plate_fraction_pooled{suffix}"] = float(frac[0])
    per_cond = by_condition["between_plate_fraction_pooled"].to_numpy(dtype=float)
    per_cond = per_cond[np.isfinite(per_cond)]
    overall["n_conditions_decomposed"] = int(per_cond.size)
    overall["between_plate_fraction_pooled_over_conditions"] = (
        float(np.mean(per_cond)) if per_cond.size else float("nan")
    )
    overall["frac_conditions_plate_dominated"] = (
        float(np.mean(per_cond > 0.5)) if per_cond.size else float("nan")
    )
    return pd.DataFrame([overall]), by_condition


NOISE_STRATA_COLUMNS = (
    "expression_quartile",
    "response_quartile",
    "n",
    "var_lfc_mean",
    "mean_se2_mean",
    "between_plate_fraction_pooled",
    "base_mean_mean",
    "abs_mean_lfc_mean",
)


def noise_strata_from_sample(noise: pd.DataFrame) -> pd.DataFrame:
    """The between-plate share by expression and response-size quartile, from the sample.

    Each stratum's share is the POOLED estimate over its rows -- the same estimator the run
    reports overall, applied within the stratum -- not a mean of per-gene shares, which would
    carry the per-gene truncation bias into every cell. Computed on the bounded sample rather
    than on every gene-condition: assigning quartiles over 175 million rows means sorting them
    twice, and a pooled share from a two-million-row sample is precise well past the third
    decimal anyone reads it to.
    """
    if noise.empty:
        return pd.DataFrame(columns=list(NOISE_STRATA_COLUMNS))
    var = noise["var_lfc"].to_numpy(dtype=float)
    se2 = noise["mean_se2"].to_numpy(dtype=float)
    ok = np.isfinite(var) & (var > 0) & np.isfinite(se2)
    d = noise.loc[ok].copy()
    d["abs_lfc"] = d["mean_lfc"].abs()
    d["expression_quartile"] = pd.qcut(
        d["base_mean"].rank(method="first"), 4, labels=[1, 2, 3, 4]
    ).astype(int)
    d["response_quartile"] = pd.qcut(
        d["abs_lfc"].rank(method="first"), 4, labels=[1, 2, 3, 4]
    ).astype(int)
    g = d.groupby(["expression_quartile", "response_quartile"], observed=True)
    out = g.agg(
        n=("var_lfc", "size"),
        var_lfc_mean=("var_lfc", "mean"),
        mean_se2_mean=("mean_se2", "mean"),
        base_mean_mean=("base_mean", "mean"),
        abs_mean_lfc_mean=("abs_lfc", "mean"),
    ).reset_index()
    sigma2 = np.maximum(
        out["var_lfc_mean"].to_numpy(dtype=float) - out["mean_se2_mean"].to_numpy(dtype=float),
        0.0,
    )
    out["between_plate_fraction_pooled"] = sigma2 / out["var_lfc_mean"].to_numpy(dtype=float)
    return out[list(NOISE_STRATA_COLUMNS)]


def masked_rowwise_pearson(
    a: np.ndarray, b: np.ndarray, min_genes: int, *, select: np.ndarray | None = None
) -> np.ndarray:
    """Pearson r per row between ``a`` and ``b``, over entries finite in both.

    Vectorized across rows; rows with fewer than ``min_genes`` shared finite entries or
    zero variance come back NaN.

    ``select`` is an optional boolean array of the same shape restricting each row to its own
    subset of columns -- rung 0's responder gene set. A false entry is treated exactly as a
    non-finite one, so every moment below (the count, both means, the covariance and both
    variances) is taken over the selected genes alone. Centring after masking is the part that
    matters: subtracting a mean computed over all genes would leave the selected columns
    off-centre and the correlation would not be the correlation of what was scored.
    """
    ok = np.isfinite(a) & np.isfinite(b)
    if select is not None:
        ok &= select
    n = ok.sum(axis=1)
    a0 = np.where(ok, a, 0.0)
    b0 = np.where(ok, b, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ma = a0.sum(axis=1) / n
        mb = b0.sum(axis=1) / n
        ac = np.where(ok, a - ma[:, None], 0.0)
        bc = np.where(ok, b - mb[:, None], 0.0)
        cov = (ac * bc).sum(axis=1)
        va = (ac**2).sum(axis=1)
        vb = (bc**2).sum(axis=1)
        r = cov / np.sqrt(va * vb)
    r[(n < min_genes) | (va <= 0) | (vb <= 0)] = np.nan
    return r


def dense_pivots(
    de: pd.DataFrame, panel: set[str], value_cols: tuple[str, ...]
) -> tuple[pd.MultiIndex, pd.Index, dict[str, np.ndarray]]:
    """Condition-by-gene matrices, built by scatter instead of ``pivot_table``.

    ``pivot_table`` groups and unstacks, and on this screen's 451 million rows its intermediates
    take more memory than the matrices they produce -- the first full-extent run reached 137 GB
    of a 200 GB allocation in that call alone. The result is a dense array either way, so this
    allocates it once and writes each row into place: linear in the rows, and the only large
    objects are the matrices themselves.

    Equivalent to ``pivot_table`` here, not merely similar. That function averages duplicate
    (row, column) pairs; the frame this reads was produced by a ``GROUP BY`` on exactly
    ``(patient, drug, gene_name)``, so no duplicates exist and a scatter and a mean agree. The
    index and columns are sorted the same way, so downstream alignment is unchanged, and a
    known-answer test asserts the two paths return identical arrays.
    """
    d = de[de["gene_name"].isin(panel)]

    def codes(col: str) -> tuple[np.ndarray, np.ndarray]:
        """Integer codes into the column's distinct values, plus those values."""
        series = d[col]
        if not isinstance(series.dtype, pd.CategoricalDtype):
            series = series.astype("category")
        cat = series.cat
        used, idx = np.unique(cat.codes.to_numpy(), return_inverse=True)
        return idx, np.asarray(cat.categories, dtype=object)[used]

    p_idx, p_names = codes("patient")
    g_idx, g_names = codes("gene_name")
    d_idx, d_names = codes("drug")
    s_idx, s_names = codes("dose")

    # One code per (patient, drug, dose), then sorted the way pivot_table sorts: by each key in
    # turn. The dictionary encoding the build returns orders categories by first appearance, not
    # lexicographically, so the sort is applied explicitly rather than assumed from the codes.
    # Dose is part of the key because this screen confounds it with plate -- see CONDITION_KEYS.
    triple = (
        p_idx.astype(np.int64) * d_names.size + d_idx.astype(np.int64)
    ) * s_names.size + s_idx.astype(np.int64)
    used, cond_idx = np.unique(triple, return_inverse=True)
    cond_dose = s_names[used % s_names.size]
    pair_used = used // s_names.size
    cond_patient = p_names[pair_used // d_names.size]
    cond_drug = d_names[pair_used % d_names.size]
    cond_order = np.lexsort((cond_dose, cond_drug, cond_patient))
    cond_rank = np.empty_like(cond_order)
    cond_rank[cond_order] = np.arange(cond_order.size)
    rows = cond_rank[cond_idx]

    # Columns keep the CATEGORY order, not a sorted one. Read off pandas rather than assumed:
    # pivot_table sorts its index but leaves a categorical column axis in category order, and
    # the dictionary encoding the build returns orders categories by first appearance. The gene
    # axis is therefore arbitrary but consistent -- which is all the scorer needs, since it
    # intersects the two halves' columns before correlating anything.
    index = pd.MultiIndex.from_arrays(
        [cond_patient[cond_order], cond_drug[cond_order], cond_dose[cond_order]],
        names=list(CONDITION_KEYS),
    )
    columns = pd.Index(g_names, name="gene_name")

    out: dict[str, np.ndarray] = {}
    for col in value_cols:
        mat = np.full((cond_order.size, g_names.size), np.nan, dtype=np.float64)
        mat[rows, g_idx] = d[col].to_numpy(dtype=float)
        out[col] = mat
    return index, columns, out


def score_split_half(
    de: pd.DataFrame,
    panel: set[str],
    min_genes: int = 50,
    *,
    select: np.ndarray | None = None,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """Per-(line, drug) split-half Pearson over the panel genes, plus the half pivots.

    ``select`` restricts each condition to its own responder genes; it must be aligned to the
    pivots this function returns, which is what ``padj_pivot`` exists to guarantee. Passing
    ``None`` reproduces the unselected statistic exactly.
    """
    # Both halves are scattered into ONE pass over the same rows, so they share an index and a
    # column axis by construction. The pivot_table path needed an index and column intersection
    # afterwards because it dropped all-NaN columns per half independently, which could leave
    # the two halves carrying different gene sets whenever their column counts happened to
    # match. That failure mode cannot arise here.
    index, columns, mats = dense_pivots(de, panel, ("lfc0", "lfc1"))
    piv0 = pd.DataFrame(mats["lfc0"], index=index, columns=columns)
    piv1 = pd.DataFrame(mats["lfc1"], index=index, columns=columns)
    r = masked_rowwise_pearson(
        piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float), min_genes, select=select
    )
    return r, piv0, piv1


def padj_pivot(de: pd.DataFrame, panel: set[str]) -> pd.DataFrame:
    """The first group's adjusted p-values, pivoted to the shape ``score_split_half`` returns.

    A separate function rather than a fourth return value so the three-tuple signature every
    existing caller unpacks -- including ``scripts/permutation_null.py`` -- keeps working. The
    caller reindexes it onto the scored pivots' rows and columns, which is the alignment step
    that guarantees a mask entry refers to the gene it appears to refer to.
    """
    index, columns, mats = dense_pivots(de, panel, ("padj0",))
    return pd.DataFrame(mats["padj0"], index=index, columns=columns)


def responder_overlap_table(de: pd.DataFrame, panel: set[str], alpha: float = 0.05) -> pd.DataFrame:
    """Per condition, how far the two plate groups agree on which genes responded.

    A DIAGNOSTIC and never an input to selection. It answers a question a reader will
    reasonably ask -- if the first group's responder call is noisy, what does the second group
    call? -- and it is exactly the quantity that must not steer the gene set, because keeping
    the genes both groups called is the pooled selection whose winner's curse the leakage
    control measures. Reported so the reader can see the agreement rate and judge the one-sided
    rule, not so the rule can be relaxed.

    Columns: ``patient``, ``drug``, ``n_first``, ``n_second``, ``n_both``, ``jaccard``.
    """
    # Counted by grouping the long frame, not by pivoting it. The counts are per condition, so
    # the two condition-by-gene matrices a pivot would build -- 7 GB each at this screen's size,
    # on top of everything else alive at this point -- are pure overhead for a per-row sum. The
    # first full-extent run was killed for memory a few steps after this call.
    d = de[de["gene_name"].isin(panel)]
    p0 = d["padj0"].to_numpy(dtype=float)
    p1 = d["padj1"].to_numpy(dtype=float)
    first = np.isfinite(p0) & (p0 < alpha)
    second = np.isfinite(p1) & (p1 < alpha)
    counts = (
        pd.DataFrame(
            {
                "patient": d["patient"].to_numpy(),
                "drug": d["drug"].to_numpy(),
                "dose": d["dose"].to_numpy(),
                "n_first": first,
                "n_second": second,
                "n_both": first & second,
                "n_union": first | second,
            }
        )
        .groupby(list(CONDITION_KEYS), observed=True, sort=True)
        .sum()
        .reset_index()
    )
    union = counts["n_union"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        jaccard = np.where(union > 0, counts["n_both"].to_numpy(dtype=float) / union, np.nan)
    counts["jaccard"] = np.round(jaccard, 4)
    return counts.drop(columns=["n_union"])


def responder_mask(piv_padj0: pd.DataFrame, alpha: float = 0.05) -> np.ndarray:
    """True where the FIRST plate group called that gene differentially expressed.

    One-sided by construction: the only input is ``padj0``, which the build aggregates over the
    first group's rows alone. There is deliberately no parameter here that could admit the
    second group -- a gene chosen using the half it is then scored against is chosen partly for
    noise that agreed by chance, and the correlation reports that agreement as reproducibility.
    A null adjusted p-value (a gene DESeq2 could not test) is not a responder.
    """
    v = piv_padj0.to_numpy(dtype=float)
    return np.isfinite(v) & (v < alpha)


def per_pair_table(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    r: np.ndarray,
    *,
    r_responder: np.ndarray | None = None,
    select: np.ndarray | None = None,
) -> pd.DataFrame:
    """The result's own data: one row per candidate (line, drug, dose) condition.

    The summary row summarizes these values; committing them makes the summary re-derivable
    anywhere -- mean, median, quartiles, positive fraction, the dose strata and the effect-size
    terciles all recompute from this table without cluster access. Columns: the split-half r
    whose mean is the declared statistic, each HALF's own mean absolute delta over the genes
    finite in both halves, and the shared finite-gene count. Rows scored NaN (fewer than
    ``min_genes`` shared genes) are kept, honestly NaN.

    The two half magnitudes are kept separate on purpose. The effect-size control ranks
    conditions by one half's magnitude and reads the correlation of the pair, then swaps
    halves; ranking by the magnitude of the two halves' SUM would select conditions whose
    halves happened to agree and make pure noise rise across the terciles.
    """
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    n = ok.sum(axis=1)
    denom = np.maximum(n, 1)
    mean_abs0 = np.where(ok, np.abs(a), 0.0).sum(axis=1) / denom
    mean_abs1 = np.where(ok, np.abs(b), 0.0).sum(axis=1) / denom
    out = pd.DataFrame(
        {
            "patient": piv0.index.get_level_values(0),
            "drug": piv0.index.get_level_values(1),
            "dose": piv0.index.get_level_values(2),
            "n_genes_scored": n,
            "mean_abs_half0": np.round(mean_abs0, 4),
            "mean_abs_half1": np.round(mean_abs1, 4),
            "r": np.round(r, 4),
        }
    )
    if r_responder is not None:
        out["r_responder"] = np.round(r_responder, 4)
    if select is not None:
        # The responders SCORED, not the responders selected: a gene the first group called but
        # whose second-group value is missing contributes to neither, so reporting the selection
        # count would overstate what the responder correlation was computed on.
        out["n_responders"] = (ok & select).sum(axis=1)
    return out


def dose_strata_table(per_pair: pd.DataFrame) -> pd.DataFrame:
    """Every candidate ceiling that can be formed from the per-condition table, side by side.

    Holding dose fixed made the scoreable unit a (line, drug, dose) triple, and the screen did
    not replicate its doses evenly, so what "the" reliability weights is a declared choice rather
    than an arithmetic fact. This table carries the candidates so the choice is read off
    committed numbers: each dose level on its own; all triples with equal weight (the summary
    row's statistic); and each (line, drug) pair weighted once, its triples averaged first.
    Within one dose level the last two coincide, so the per-pair weighting appears only for the
    pooled row. Both gene sets, each with its count, mean, median and Spearman-Brown value.
    """
    rows: list[dict[str, object]] = []

    def _row(dose: str, weighting: str, gene_set: str, r: np.ndarray) -> None:
        r = r[np.isfinite(r)]
        mean = float(np.mean(r)) if r.size else float("nan")
        rows.append(
            {
                "dose": dose,
                "weighting": weighting,
                "gene_set": gene_set,
                "n_pairs": int(r.size),
                "splithalf_mean_r": round(mean, 4),
                "splithalf_median_r": round(float(np.median(r)), 4) if r.size else float("nan"),
                "spearman_brown_full": round(spearman_brown_or_nan(mean), 4),
            }
        )

    gene_sets = [("all", "r")] + (
        [("responder", "r_responder")] if "r_responder" in per_pair else []
    )
    for gene_set, col in gene_sets:
        r = per_pair[col].to_numpy(dtype=float)
        doses = per_pair["dose"].astype(str).to_numpy()
        for dose in sorted(set(doses), key=lambda d: float(d) if _is_number(d) else d):
            _row(dose, "per_triple", gene_set, r[doses == dose])
        _row("all", "per_triple", gene_set, r)
        finite = np.isfinite(r)
        pair_means = (
            pd.DataFrame({"patient": per_pair["patient"], "drug": per_pair["drug"], "r": r})
            .loc[finite]
            .groupby(["patient", "drug"], observed=True)["r"]
            .mean()
            .to_numpy(dtype=float)
        )
        _row("all", "per_line_drug", gene_set, pair_means)
    return pd.DataFrame(rows)


def _is_number(text: str) -> bool:
    try:
        float(text)
    except ValueError:
        return False
    return True


def stratified_null_draws(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    n_perm: int = 500,
    seed: int = 0,
    min_genes: int = 50,
    *,
    select: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Mismatched-pair null correlations per stratum.

    any_pair: two different conditions (continuity with the archived lineage's first run).
    diff_drug: different line AND drug -- the generic-structure floor the ceiling clears.
    same_drug: same drug AT THE SAME DOSE, different line -- the line-specificity floor. Dose
    is held fixed here for the same reason it is held fixed in the condition: a same-drug pair
    at different doses is partly a dose contrast, and the floor would then be lower than the
    line-specificity it is meant to measure.

    With ``select``, a draw pairing condition *i*'s first group against condition *j*'s second
    group scores over **row i's** selected genes -- the row whose first group is used, which is
    the row the selection rule would have read. Using row *j*'s mask, or the union of the two,
    would apply a different rule to the null than to the observed value and the comparison would
    stop being like for like. The finiteness rule then intersects with row *j*'s second group,
    exactly as it does for a matched pair.
    """
    lines = piv0.index.get_level_values(0).to_numpy(dtype=str)
    drugs = piv0.index.get_level_values(1).to_numpy(dtype=str)
    doses = piv0.index.get_level_values(2).to_numpy(dtype=str)
    n = len(piv0)
    ii, jj = np.divmod(np.arange(n * n), n)
    off = ii != jj
    ii, jj = ii[off], jj[off]
    same_drug = drugs[ii] == drugs[jj]
    same_dose = doses[ii] == doses[jj]
    same_line = lines[ii] == lines[jj]
    strata = {
        "any_pair": np.ones(ii.size, dtype=bool),
        "diff_drug": ~same_drug & ~same_line,
        "same_drug": same_drug & same_dose & ~same_line,
    }
    rng = np.random.default_rng(seed)
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    out: dict[str, np.ndarray] = {}
    for name, mask in strata.items():
        avail = np.flatnonzero(mask)
        if avail.size == 0:
            out[name] = np.array([])
            continue
        pick = rng.choice(avail, size=min(n_perm, avail.size), replace=False)
        sel = select[ii[pick]] if select is not None else None
        r = masked_rowwise_pearson(a[ii[pick]], b[jj[pick]], min_genes, select=sel)
        out[name] = r[np.isfinite(r)]
    return out


def null_draw_table(nulls: dict[str, np.ndarray]) -> pd.DataFrame:
    """Every individual mismatched-pair correlation, long-format (stratum, r).

    The summary row reports only each stratum's mean, so the chance floors could be quoted
    but not seen. These are the draws those means average, committed so the floor
    distributions can be drawn and their means re-derived off-cluster.
    """
    return pd.DataFrame(
        {
            "stratum": np.repeat(list(nulls), [len(v) for v in nulls.values()]),
            "r": np.round(np.concatenate([np.asarray(v, dtype=float) for v in nulls.values()]), 4),
        }
    )


def example_pair_profiles(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    r: np.ndarray,
    quantiles: tuple[float, ...] = (0.05, 0.25, 0.5, 0.95),
    max_genes: int | None = None,
    seed: int = 0,
    select: np.ndarray | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Gene-level half-profiles for a few example conditions, so an r can be seen as a scatter.

    Every correlation in this analysis is one point in a distribution; nothing committed showed
    what the underlying agreement looks like gene by gene. This exports both halves' per-gene
    deltas for matched conditions spanning the reliability range (at ``quantiles`` of the sorted
    finite r, nearest-rank) plus the two mismatched comparisons the chance floors are built from,
    each formed by pairing the median condition's first half with another condition's second
    half: same drug and different line, then different drug and line. Deterministic -- selection
    is by sorted rank and by first-in-index-order among candidates.

    Every shared gene is exported by default and the file is written gzipped, which keeps it
    around 700 kilobytes -- a quarter of the plain-text size -- without any sampling error. An
    earlier version subsampled to 2,000 genes for size; on the real pool that put +/-0.04 of
    sampling error on correlations of 0.03-0.11, enough to reorder the examples (measured:
    r_full 0.028/0.071/0.109/0.354 came out as 0.062/0.054/0.069/0.357), so a panel captioned
    as spanning the reliability range would have shown correlations that do not rise.
    ``max_genes`` keeps that path available for pools where size forces it; with it set, the
    subsample is seeded and a rerun reproduces the file exactly.

    The index carries both ``r_full`` over every shared gene (the reported quantity, matching
    ``rung0_per_pair_r.csv`` for matched examples) and ``r_shown`` over exactly the exported
    points, so the committed file verifies against itself; with no subsampling the two are
    equal, which is itself the check that nothing was dropped.

    Returns (profiles, index): the long per-gene frame keyed by ``example_id``, and one row per
    example carrying its kind, the two conditions' labels, both gene counts, and both r values.
    """
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    lines = piv0.index.get_level_values(0).to_numpy(dtype=str)
    drugs = piv0.index.get_level_values(1).to_numpy(dtype=str)
    doses = piv0.index.get_level_values(2).to_numpy(dtype=str)
    order = np.flatnonzero(np.isfinite(r))[np.argsort(r[np.isfinite(r)])]
    if order.size == 0:
        index_cols = ["example_id", "kind", "patient0", "drug0", "dose0"]
        index_cols += ["patient1", "drug1", "dose1"]
        index_cols += ["n_genes_full", "r_full", "n_genes_shown", "r_shown"]
        index_cols += ["n_responders_shown", "r_responder_full"]
        return pd.DataFrame(
            columns=["example_id", "gene", "lfc0", "lfc1", "is_responder"]
        ), pd.DataFrame(columns=index_cols)

    def _at(q: float) -> int:
        return int(order[round(q * (order.size - 1))])

    picks: list[tuple[str, str, int, int]] = [
        (f"matched_q{round(q * 100):02d}", "matched", _at(q), _at(q)) for q in quantiles
    ]
    anchor = _at(0.5)
    for kind, mask in (
        (
            "same_drug_mismatch",
            (drugs == drugs[anchor]) & (doses == doses[anchor]) & (lines != lines[anchor]),
        ),
        ("diff_drug_mismatch", (drugs != drugs[anchor]) & (lines != lines[anchor])),
    ):
        candidates = np.flatnonzero(mask)
        if candidates.size:
            picks.append((kind, kind, anchor, int(candidates[0])))

    frames: list[pd.DataFrame] = []
    rows: list[dict[str, object]] = []
    genes = piv0.columns.to_numpy()
    rng = np.random.default_rng(seed)
    for example_id, kind, i, j in picks:
        shared = np.flatnonzero(np.isfinite(a[i]) & np.isfinite(b[j]))
        r_full = float(masked_rowwise_pearson(a[i][None, :], b[j][None, :], min_genes=1)[0])
        shown = (
            np.sort(rng.choice(shared, size=max_genes, replace=False))
            if max_genes is not None and shared.size > max_genes
            else shared
        )
        x, y = a[i][shown], b[j][shown]
        r_shown = float(masked_rowwise_pearson(x[None, :], y[None, :], min_genes=1)[0])
        # Which of the exported points are the FIRST condition's responders -- the same row the
        # selection rule reads, including for a mismatched example, so the marking means the same
        # thing everywhere. Exported so the design's second scatter can be drawn from the
        # committed table rather than recomputed from data the figure does not have.
        resp = select[i][shown] if select is not None else np.zeros(shown.size, dtype=bool)
        r_resp = (
            float(masked_rowwise_pearson(a[i][None, :], b[j][None, :], 1, select=select[[i]])[0])
            if select is not None
            else float("nan")
        )
        frames.append(
            pd.DataFrame(
                {
                    "example_id": example_id,
                    "gene": genes[shown],
                    "lfc0": np.round(x, 4),
                    "lfc1": np.round(y, 4),
                    "is_responder": resp,
                }
            )
        )
        rows.append(
            {
                "example_id": example_id,
                "kind": kind,
                "patient0": lines[i],
                "drug0": drugs[i],
                "dose0": doses[i],
                "patient1": lines[j],
                "drug1": drugs[j],
                "dose1": doses[j],
                "n_genes_full": int(shared.size),
                "r_full": round(r_full, 4),
                "n_genes_shown": int(shown.size),
                "r_shown": round(r_shown, 4),
                "n_responders_shown": int(resp.sum()),
                "r_responder_full": round(r_resp, 4) if np.isfinite(r_resp) else float("nan"),
            }
        )
    return pd.concat(frames, ignore_index=True), pd.DataFrame(rows)


def _tercile_masks(magnitude: np.ndarray, r: np.ndarray) -> list[np.ndarray]:
    """Three boolean masks cutting the finite conditions into thirds of ``magnitude``."""
    finite = np.isfinite(r) & np.isfinite(magnitude)
    if not finite.any():
        return [finite, finite, finite]
    edges = np.quantile(magnitude[finite], [1 / 3, 2 / 3])
    masks = []
    for t in (1, 2, 3):
        lo = -np.inf if t == 1 else edges[t - 2]
        hi = np.inf if t == 3 else edges[t - 1]
        masks.append(finite & (magnitude > lo) & (magnitude <= hi))
    return masks


def effect_size_terciles(
    per_pair: pd.DataFrame, ranked_by: str = "half0", r_col: str = "r"
) -> dict[str, float]:
    """Split-half mean r within terciles of ONE half's effect size (mean |delta|).

    The empirical positive control: an assay that cannot find more reproducibility where there
    is more signal is broken. Tercile 1 = smallest effects. Ranked by ``mean_abs_{ranked_by}``,
    one half's magnitude alone. Under no signal the second half is independent of the first, so
    conditioning on the first half's magnitude leaves the expected correlation at zero in every
    tercile; ranking by the magnitude of the two halves' SUM does not have that property, since
    a large sum is reached most easily when the halves agree, and pure noise would rise.
    """
    masks = _tercile_masks(
        per_pair[f"mean_abs_{ranked_by}"].to_numpy(dtype=float),
        per_pair[r_col].to_numpy(dtype=float),
    )
    r = per_pair[r_col].to_numpy(dtype=float)
    return {
        f"splithalf_mean_r_tercile{t}": round(float(np.mean(r[m])), 3) if m.any() else float("nan")
        for t, m in zip((1, 2, 3), masks, strict=True)
    }


def per_gene_reliability(
    piv0: pd.DataFrame, piv1: pd.DataFrame, min_pairs: int = 20
) -> pd.DataFrame:
    """The transpose diagnostic: each gene's delta correlated across pairs between halves.

    Unpromoted (see design.md): says which panel genes carry reproducible perturbation
    signal, as the evidence base for any future panel restriction.
    """
    a, b = piv0.to_numpy(dtype=float).T, piv1.to_numpy(dtype=float).T
    r = masked_rowwise_pearson(a, b, min_pairs)
    n = (np.isfinite(a) & np.isfinite(b)).sum(axis=1)
    return pd.DataFrame(
        {"gene": piv0.columns.to_numpy(), "n_pairs": n, "r": np.round(r, 4)}
    ).sort_values("r", ascending=False)


def spearman_brown_or_nan(r: float) -> float:
    """``spearman_brown`` with the guard its docstring asks callers to supply.

    The lift is undefined at r = -1 (a zero denominator) and meaningless below it. One guarded
    entry point means the two gene sets, the even-plate subset and any later rung all apply the
    same correction with the same guard, rather than each re-deciding what to do at the edge.
    """
    from fmharness.statistics import spearman_brown

    return spearman_brown(r) if r > -1 else float("nan")


def summarize(
    r: np.ndarray,
    nulls: dict[str, np.ndarray],
    seed: int = 0,
    *,
    label: str = "",
    even_mask: np.ndarray | None = None,
) -> dict:
    """One reliability's row: mean-over-conditions Pearson (the declared statistic), its nulls,
    p-values from the bootstrapped null aggregate, and the MDEs (SPEC rule 4).

    ``label`` prefixes every key, so the same function serves the all-gene and responder
    statistics and both land in ONE summary row. Two files would let one number be quoted
    without the other; one row with two prefixed families cannot.

    The Spearman-Brown lift is applied to the MEAN over conditions, not per condition and then
    averaged. ``2r/(1+r)`` is not linear, so the two differ, and the design's declared statistic
    is the mean.

    ``even_mask`` marks the conditions whose plate count splits into two equal groups, where the
    correction's equal-halves assumption holds exactly. Its corrected value is reported beside
    the full one, and the gap between them is the size of that assumption rather than an
    argument about it. The mask is over ``r`` BEFORE non-finite entries are dropped, so it is
    the caller's per-condition mask and needs no realignment here.
    """
    from fmharness.statistics import bootstrap_aggregate_pvalue, minimum_detectable_aggregate

    finite = np.isfinite(r)
    r_even = r[finite & even_mask] if even_mask is not None else np.array([])
    r = r[finite]
    mean = float(np.mean(r))
    nl = nulls["diff_drug"] if nulls["diff_drug"].size else nulls["any_pair"]
    p_boot, ci_lo, ci_hi = bootstrap_aggregate_pvalue(mean, nl, r.size, seed=seed)
    p_same = bootstrap_aggregate_pvalue(mean, nulls["same_drug"], r.size, seed=seed)[0]
    mean_even = float(np.mean(r_even)) if r_even.size else float("nan")
    out = {
        "n_pairs": int(r.size),
        "splithalf_mean_r": round(mean, 3),
        "splithalf_median_r": round(float(np.median(r)), 3),
        "splithalf_q1_r": round(float(np.quantile(r, 0.25)), 3),
        "splithalf_q3_r": round(float(np.quantile(r, 0.75)), 3),
        "spearman_brown_full": round(spearman_brown_or_nan(mean), 3),
        "n_pairs_even": int(r_even.size),
        "splithalf_mean_r_even_plates": round(mean_even, 3),
        "spearman_brown_full_even_plates": round(spearman_brown_or_nan(mean_even), 3)
        if r_even.size
        else float("nan"),
        "frac_pos": round(float(np.mean(r > 0)), 3),
        "null_any_pair_mean_r": round(float(np.mean(nulls["any_pair"])), 3),
        "null_diff_drug_mean_r": round(float(np.mean(nulls["diff_drug"])), 3),
        "null_same_drug_mean_r": round(float(np.mean(nulls["same_drug"])), 3),
        "null_n_draws": int(nl.size),
        "p_vs_null": round(p_boot, 4),
        "p_vs_same_drug": round(p_same, 4),
        "null_mean_ci_lo": round(ci_lo, 3),
        "null_mean_ci_hi": round(ci_hi, 3),
        "mde_80_vs_diff_drug": round(minimum_detectable_aggregate(r, nl, r.size, seed=seed), 4),
        "mde_80_vs_same_drug": round(
            minimum_detectable_aggregate(r, nulls["same_drug"], r.size, seed=seed), 4
        ),
    }
    return {f"{label}_{k}" if label else k: v for k, v in out.items()}


def write_per_gene_figure(per_gene: pd.DataFrame, out_png: Path) -> None:
    """Histogram of the per-gene split-half diagnostic (design.md, 'per-gene reliability').

    Unpromoted, same as the CSV it reads: says which panel genes carry reproducible
    perturbation signal, not the pair-level ceiling `write_figure` reports.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    r = per_gene["r"].to_numpy(dtype=float)
    r = r[np.isfinite(r)]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(r, bins=60, alpha=0.75, color="tab:blue")
    ax.axvline(0.0, color="k", lw=1.0, linestyle="--", label="zero")
    ax.axvline(float(np.median(r)), color="tab:red", lw=1.5, label=f"median ({np.median(r):.3f})")
    ax.set_xlabel("split-half r per gene, across (line, drug) conditions")
    ax.set_ylabel("count")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def _write_params_sidecar(result_path, args_ns, extra=None) -> None:
    """Record the git sha and every resolved argument beside the result.

    A ceiling used as a denominator has to be checkable against a rerun; a bare CSV is a number
    with no way back to the code and parameters that produced it. This script had none, which is
    how its value lived in doc prose for weeks with nothing behind it.
    """
    import json as _json
    import subprocess as _sp
    from pathlib import Path as _P

    try:
        sha = _sp.run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        ).stdout.strip()
    except Exception:
        sha = "unknown"
    import os as _os

    side = _P(str(result_path)).with_suffix(".params.json")
    side.write_text(
        _json.dumps(
            {
                "result": _P(str(result_path)).name,
                "git_sha": sha,
                "slurm_job_id": _os.environ.get("SLURM_JOB_ID", "local"),
                "args": {k: str(v) for k, v in vars(args_ns).items()},
                **(extra or {}),
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {side}")


def frame_cache_key(paths: list[str], names: list[str] | None, replicate_col: str | None) -> str:
    """Content key for a built split-half frame: the inputs that determine it, hashed.

    Naming the cache after its inputs is what makes it safe -- a different pool, drug set,
    replicate column or split rule resolves to a different file rather than silently reusing
    the wrong frame. ``FRAME_SCHEMA`` is part of the payload for exactly that reason.
    """
    drug_key = ["<all drugs>"] if not names else sorted(names)
    payload = "\n".join(
        [*sorted(paths), "--", *drug_key, "--", str(replicate_col), "--", FRAME_SCHEMA]
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def slice_paths(cache_dir: Path, key: str, n_parts: int, part: int) -> dict[str, Path]:
    """Where one slice's three products and its completion record live in the cache."""
    tag = f"{key}_{n_parts}_{part}"
    return {
        "frame": cache_dir / f"frame_{tag}.parquet",
        "noise": cache_dir / f"noise_{tag}.parquet",
        "sample": cache_dir / f"noise_sample_{tag}.parquet",
        "done": cache_dir / f"slice_{tag}.json",
    }


def _normalise_keys(de: pd.DataFrame) -> pd.DataFrame:
    """Give the key columns the dictionary-encoded dtype the builder produces.

    Parquet stores a dictionary-encoded float back as a plain float, so a frame read from cache
    had ``dose`` as float64 where a freshly built one had it as a Categorical. Everything
    downstream coped with either, which is precisely the kind of difference that goes unnoticed
    until something does not -- so the two paths are made to agree here instead.
    """
    for col in KEY_COLUMNS:
        if col in de.columns and not isinstance(de[col].dtype, pd.CategoricalDtype):
            de[col] = de[col].astype("category")
    return de


def _cache_dir_for(args: argparse.Namespace) -> Path:
    """The slice cache: the one given, else ``cache/`` under the output directory."""
    given = getattr(args, "cache_dir", None)
    return Path(given) if given else Path(args.out_dir) / "cache"


def _build_or_load_frame(
    paths: list[str], names: list[str] | None, args: argparse.Namespace, local: Path
) -> tuple[pd.DataFrame, str]:
    """The split-half frame, from the slice cache when every slice is there, else built.

    Building it scans every DE shard through DuckDB and dominates the run; everything after it
    takes minutes. The cluster builds the slices as a job array (``--stage slice``) and this
    reads them back; a local run with no slices cached builds them here, one after another, and
    caches each, so a second local run -- adding an output, changing a figure -- skips the scan.
    """
    cache_dir = _cache_dir_for(args)
    n_parts = max(1, int(getattr(args, "slices", 1)))
    key = frame_cache_key(paths, names, getattr(args, "replicate_col", None))
    missing = [
        p for p in range(n_parts) if not slice_paths(cache_dir, key, n_parts, p)["done"].exists()
    ]
    if not missing and (cache_dir / ASSIGNMENT_SIDECAR).exists():
        de, _, _, repl = load_slices(cache_dir, key, n_parts)
        print(f"loaded the split-half frame from {n_parts} cached slice(s) in {cache_dir}")
        return de, repl
    assignment = _ensure_cached_assignment(paths, names, args, local, cache_dir)
    for part in range(n_parts):
        run_slice(paths, names, args, local, cache_dir, assignment, part)
    de, _, _, repl = load_slices(cache_dir, key, n_parts)
    return de, repl


def _ensure_cached_assignment(
    paths: list[str],
    names: list[str] | None,
    args: argparse.Namespace,
    local: Path,
    cache_dir: Path,
) -> Path:
    """The split assignment in the cache, computed if it is not there yet."""
    if (cache_dir / ASSIGNMENT_SIDECAR).exists():
        path, _, _ = read_assignment(cache_dir)
        return path
    assignment, pool, repl = split_assignment(
        paths,
        names,
        getattr(args, "replicate_col", None),
        local.parent / "duckdb_tmp",
        getattr(args, "duckdb_memory", "36GB"),
        getattr(args, "duckdb_threads", None),
    )
    print(
        f"split assignment: {len(pool):,} (line, drug, dose) triples, "
        f"{int((pool['n_plates'] >= 2).sum()):,} with two or more plates; {SPLIT_RULE}"
    )
    return write_assignment(cache_dir, assignment, pool, repl)


def effect_size_tercile_table(
    per_pair: pd.DataFrame, n_boot: int = 2000, seed: int = 0, r_col: str = "r"
) -> pd.DataFrame:
    """The empirical control as a table, cross-fit, with an interval on each tercile mean.

    Two rankings, one per half: conditions in thirds by the first half's magnitude, scored by
    the correlation of the pair; then by the second half's magnitude. Both are reported and both
    must rise. Reading either half's magnitude alone keeps the control clean under no signal
    (see ``effect_size_terciles``); reporting both makes the asymmetry of the halves visible
    rather than a choice. Recomputable from the committed per-condition table alone, which is
    what lets the verification battery re-derive it.
    """
    rng = np.random.default_rng(seed)
    r = per_pair[r_col].to_numpy(dtype=float)
    rows = []
    for ranked_by in ("half0", "half1"):
        masks = _tercile_masks(per_pair[f"mean_abs_{ranked_by}"].to_numpy(dtype=float), r)
        for t, sel in zip((1, 2, 3), masks, strict=True):
            vals = r[sel]
            if vals.size:
                boot = np.array(
                    [np.mean(rng.choice(vals, size=vals.size, replace=True)) for _ in range(n_boot)]
                )
                lo, hi = float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))
                mean = float(np.mean(vals))
            else:
                lo = hi = mean = float("nan")
            rows.append(
                {
                    "ranked_by": ranked_by,
                    "tercile": t,
                    "n": int(vals.size),
                    "mean_r": round(mean, 4),
                    "ci_lo": round(lo, 4),
                    "ci_hi": round(hi, 4),
                }
            )
    return pd.DataFrame(rows)


def mde_curve_table(
    r_all: np.ndarray,
    r_resp: np.ndarray,
    nulls_all: dict[str, np.ndarray],
    nulls_resp: dict[str, np.ndarray],
    seed: int = 0,
) -> pd.DataFrame:
    """Minimum detectable effect against condition count, for both gene sets.

    A single MDE says whether this screen was powered; the curve says how much of that power is
    the screen's size rather than the effect's, which is the question the organoid rung will ask
    with a tenth of the conditions.
    """
    from fmharness.statistics import minimum_detectable_aggregate

    rows = []
    for gene_set, r, nulls in (("all", r_all, nulls_all), ("responder", r_resp, nulls_resp)):
        r = r[np.isfinite(r)]
        if r.size < 2:
            continue
        nl = nulls["diff_drug"] if nulls["diff_drug"].size else nulls["any_pair"]
        grid = sorted({*np.unique(np.geomspace(10, max(r.size, 11), 12).astype(int)), r.size})
        for n in grid:
            rows.append(
                {
                    "gene_set": gene_set,
                    "n_pairs": int(n),
                    "mde": round(minimum_detectable_aggregate(r, nl, int(n), seed=seed), 4),
                    "observed": bool(n == r.size),
                }
            )
    return pd.DataFrame(rows)


def leakage_table(min_genes: int, seed: int = 0) -> pd.DataFrame:
    """The one-sided rule beside the pooled one, on a pool with no signal at all.

    The design forbids selecting on the pooled data because it inflates the correlation by
    winner's curse: writing the halves as a and b, their sum and difference are independent, so
    selecting on a large |a + b| inflates var(a + b) alone and cov(a, b) = (var(a+b) -
    var(a-b))/4 goes positive with nothing generating it. Computed at run time on a signal-free
    pool so the figure shows the size of the effect being avoided rather than asserting it.
    """
    from fmharness.synthetic import planted_split_half_frame

    pool = planted_split_half_frame(
        n_lines=8, n_drugs=4, n_genes=2000, signal_sd=0.0, noise_sd=1.0, seed=seed
    )
    panel = set(pool["gene_name"].unique())
    _, piv0, piv1 = score_split_half(pool, panel, min_genes=min_genes)
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    one = responder_mask(padj_pivot(pool, panel).reindex(columns=piv0.columns).loc[piv0.index])
    k = max(int(one.sum(axis=1).mean()), min_genes)
    order = np.argsort(-np.abs(a + b), axis=1)
    pooled = np.zeros_like(one)
    np.put_along_axis(pooled, order[:, :k], True, axis=1)
    return pd.DataFrame(
        [
            {
                "rule": "one-sided",
                "mean_r": round(
                    float(np.nanmean(masked_rowwise_pearson(a, b, min_genes, select=one))), 4
                ),
                "genes_per_condition": int(one.sum(axis=1).mean()),
            },
            {
                "rule": "pooled",
                "mean_r": round(
                    float(np.nanmean(masked_rowwise_pearson(a, b, min_genes, select=pooled))), 4
                ),
                "genes_per_condition": int(pooled.sum(axis=1).mean()),
            },
        ]
    )


AUDIT_SUMS = "audit_checksums.json"


def write_audit_checksums(out_dir: Path) -> Path:
    """The sha256 of every table this run wrote, for the audit to cite and promotion to check.

    The audit reads these artifacts in the working tree, before they are committed (PROCESS
    section 1). Recording each one's checksum is what closes the window between what was audited
    and what gets committed: promotion refuses when a checksum has moved since.
    """
    import hashlib as _h

    sums = {
        p.name: _h.sha256(p.read_bytes()).hexdigest()
        for p in sorted(out_dir.rglob("*"))
        if p.is_file() and p.suffix in {".csv", ".gz", ".png", ".json"} and p.name != AUDIT_SUMS
    }
    path = out_dir / AUDIT_SUMS
    path.write_text(json.dumps(sums, indent=2, sort_keys=True) + "\n")
    print(f"wrote {path} ({len(sums)} artifacts)")
    return path


def run_slice(
    paths: list[str],
    names: list[str] | None,
    args: argparse.Namespace,
    local: Path,
    cache_dir: Path,
    assignment: Path,
    part: int,
) -> None:
    """One slice of the genes: aggregate once, cache the frame, the noise sums and a sample.

    The cluster runs this as one task of a job array, sixteen at a time on sixteen nodes, so
    the wall clock is one scan rather than sixteen in a row and each task holds a sixteenth of
    the group table. A finished slice is recorded by its completion file, written last, so a
    task killed mid-write leaves no half-slice that a later combine could mistake for a whole
    one; a rerun skips the slices that finished.
    """
    n_parts = max(1, int(getattr(args, "slices", 1)))
    key = frame_cache_key(paths, names, getattr(args, "replicate_col", None))
    files = slice_paths(cache_dir, key, n_parts, part)
    if files["done"].exists():
        print(f"  slice {part + 1}/{n_parts}: already cached, skipping")
        return
    agg, repl = slice_aggregate(
        paths,
        names,
        getattr(args, "replicate_col", None),
        local.parent / "duckdb_tmp",
        getattr(args, "duckdb_memory", "36GB"),
        assignment=assignment,
        n_parts=n_parts,
        part=part,
        threads=getattr(args, "duckdb_threads", None),
    )
    frame = frame_from_slice(agg)
    noise = noise_from_slice(agg)
    del agg
    per = noise_partials(noise, float(getattr(args, "padj_threshold", 0.05)))
    per_slice_sample = max(1, int(getattr(args, "noise_sample_rows", 2_000_000)) // n_parts)
    seed = int(getattr(args, "seed", 0))
    sample = noise.sample(min(per_slice_sample, len(noise)), random_state=seed).reset_index(
        drop=True
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(files["frame"], index=False)
    per.to_parquet(files["noise"], index=False)
    sample.to_parquet(files["sample"], index=False)
    record = {
        "part": part,
        "n_parts": n_parts,
        "replicate_col": repl,
        "schema": FRAME_SCHEMA,
        "n_frame_rows": len(frame),
        "n_noise_rows": len(noise),
        "n_sample_rows": len(sample),
    }
    files["done"].write_text(json.dumps(record) + "\n")
    print(
        f"  slice {part + 1}/{n_parts}: {len(frame):,} scoreable gene-conditions, "
        f"{len(noise):,} decomposable, cached in {cache_dir}"
    )


def load_slices(
    cache_dir: Path, key: str, n_parts: int
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    """Every slice's products, concatenated: frame, noise sums, sample, and the replicate column."""
    missing = [
        p for p in range(n_parts) if not slice_paths(cache_dir, key, n_parts, p)["done"].exists()
    ]
    if missing:
        raise SystemExit(
            f"{len(missing)} of {n_parts} slices are not in {cache_dir}: parts {missing}. "
            "Every slice must finish before the combine stage; rerun those array indices"
        )
    frames, noises, samples, repl = [], [], [], ""
    for part in range(n_parts):
        files = slice_paths(cache_dir, key, n_parts, part)
        repl = str(json.loads(files["done"].read_text())["replicate_col"])
        frames.append(pd.read_parquet(files["frame"]))
        noises.append(pd.read_parquet(files["noise"]))
        samples.append(pd.read_parquet(files["sample"]))
    de = _normalise_keys(pd.concat(frames, ignore_index=True))
    return (
        de,
        pd.concat(noises, ignore_index=True),
        pd.concat(samples, ignore_index=True),
        repl,
    )


def write_noise_outputs(
    per_condition: pd.DataFrame, sample: pd.DataFrame, args: argparse.Namespace, out_dir: Path
) -> pd.DataFrame:
    """The decomposition's tables, from the slices' sums, and the figure sample it returns.

    Every promoted number here is pooled over ALL gene-conditions from the per-condition sums
    the slices wrote; the sample is what the figures and the row-wise identity check are drawn
    from, and no promoted claim is read from it.
    """
    overall, by_condition = combine_noise_partials(per_condition)
    by_condition.round(6).to_csv(out_dir / "rung0_noise_by_condition.csv", index=False)
    sample.to_csv(out_dir / "rung0_noise_per_gene.csv.gz", index=False)
    noise_strata_from_sample(sample).round(5).to_csv(
        out_dir / "rung0_noise_strata.csv", index=False
    )
    noise_summary = {k: (v.item() if hasattr(v, "item") else v) for k, v in overall.iloc[0].items()}
    noise_summary["n_sample_rows"] = len(sample)
    noise_path = out_dir / "rung0_noise_decomposition.csv"
    pd.DataFrame([noise_summary]).round(6).to_csv(noise_path, index=False)
    _write_params_sidecar(
        noise_path,
        args,
        extra={
            **noise_summary,
            "estimator": "pooled: max(mean(var_lfc) - mean(lfcSE^2), 0), floored once after "
            "averaging over gene-conditions, never per gene",
            "split_rule": SPLIT_RULE,
        },
    )
    print(
        f"per-condition noise: {len(by_condition):,} dose-conditions; between-plate share "
        f"{noise_summary['between_plate_fraction_pooled']:.4f} pooled over gene-conditions, "
        f"{noise_summary['between_plate_fraction_pooled_over_conditions']:.4f} over conditions, "
        f"responders {noise_summary['between_plate_fraction_pooled_responders']:.4f}"
    )
    return sample


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local-dir", required=True, help="dir with the Tahoe DE parquet (on scratch)")
    ap.add_argument(
        "--drugs-cid-file",
        default="",
        help="optional PubChem CID list. Rung 0 measures at the assay's full extent, so the "
        "default is empty and every splittable drug is admitted.",
    )
    ap.add_argument(
        "--drug-names-file",
        default=None,
        help="one Tahoe drug name per line; bypasses the HuggingFace name lookup so fixtures "
        "and offline runs need no `datasets` import.",
    )
    ap.add_argument("--replicate-col", default=None, help="plate/replicate column (auto-detected)")
    ap.add_argument("--min-genes", type=int, default=50, help="min shared genes to score a pair")
    ap.add_argument("--padj-threshold", type=float, default=0.05, help="responder selection alpha")
    ap.add_argument(
        "--panel-file",
        default=None,
        help="one gene per line. Rung 0 scores every gene the table carries and passes none; a "
        "later rung computing a restriction of this ceiling supplies its own.",
    )
    ap.add_argument("--n-perm", type=int, default=500, help="mismatched-pair null draws")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default="rung0_outputs")
    ap.add_argument(
        "--stage",
        choices=("all", "assign", "slice", "combine"),
        default="all",
        help="which part of the run to do. The scan over the table is the whole cost, so on "
        "the cluster it is split three ways: `assign` scans the key columns once and writes the "
        "plate-to-half assignment; `slice` (one job-array task per --slice-index) scans the "
        "table once for its slice of the genes and caches that slice's frame and noise sums; "
        "`combine` reads every cached slice and does everything that needs the whole frame at "
        "once. `all` does the three in one process, which is how the fixture runs go.",
    )
    ap.add_argument(
        "--slices",
        type=int,
        default=1,
        help="how many slices of the GENES the scan is split into. The combination is exact: "
        "gene is part of every group key, so each gene-condition lands in exactly one slice "
        "and per-slice frames concatenate and per-slice sums add.",
    )
    ap.add_argument("--slice-index", type=int, default=None, help="which slice, for --stage slice")
    ap.add_argument(
        "--cache-dir",
        default=None,
        help="directory holding the split assignment and every slice's cached products, keyed "
        "by a hash of the inputs. Default: cache/ under --out-dir.",
    )
    ap.add_argument(
        "--duckdb-memory",
        default="36GB",
        help="DuckDB's memory_limit. The engine spills to --local-dir's parent when it runs "
        "out, so a value below the job's allocation trades speed for safety; set it well under "
        "the SBATCH --mem so the returned pandas frames still have room.",
    )
    ap.add_argument(
        "--duckdb-threads",
        type=int,
        default=None,
        help="DuckDB thread count. Fewer threads means fewer partial hash tables held at once, "
        "which is the difference between spilling and an out-of-memory error on the wide "
        "group-bys.",
    )
    ap.add_argument(
        "--noise-sample-rows",
        type=int,
        default=2_000_000,
        help="rows of the per-gene noise table to keep for the figures and the row-wise "
        "identity check, spread equally over the slices. Every reported number is pooled over "
        "all rows from the slices' sums; this only bounds what is committed.",
    )
    args = ap.parse_args()
    repo = Path(__file__).resolve().parent.parent

    local = Path(args.local_dir)
    local = local if local.is_absolute() else repo / local
    paths = sorted(str(p) for p in local.rglob("*.parquet") if DE in str(p))
    if not paths:
        raise SystemExit(f"no {DE} parquet under {local}")
    names = resolve_drug_names(repo, args)
    scope = "all drugs (no drug list given)" if names is None else f"{len(names)} target drugs"
    print(f"{scope}; reading {len(paths)} DE parquet files ...")

    out_dir = Path(args.out_dir) if Path(args.out_dir).is_absolute() else repo / args.out_dir
    args.out_dir = str(out_dir)
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = _cache_dir_for(args)
    n_parts = max(1, int(args.slices))
    key = frame_cache_key(paths, names, args.replicate_col)

    # --- assign: the split, as a table ------------------------------------------------------
    if args.stage in ("assign", "all"):
        _ensure_cached_assignment(paths, names, args, local, cache_dir)
        _write_split_tables(cache_dir, out_dir)
        if args.stage == "assign":
            return

    # --- slice: one job-array task ------------------------------------------------------------
    if args.stage == "slice":
        if args.slice_index is None:
            raise SystemExit("--stage slice needs --slice-index (the job array's task id)")
        assignment, _, _ = read_assignment(cache_dir)
        run_slice(paths, names, args, local, cache_dir, assignment, int(args.slice_index))
        return
    if args.stage == "all":
        assignment, _, _ = read_assignment(cache_dir)
        for part in range(n_parts):
            run_slice(paths, names, args, local, cache_dir, assignment, part)

    # --- combine: everything that needs the whole frame ---------------------------------------
    _, pool, _ = read_assignment(cache_dir)
    _write_split_tables(cache_dir, out_dir)
    de, per_condition, noise_sample, repl = load_slices(cache_dir, key, n_parts)
    if de.empty:
        raise SystemExit("no (line, drug, dose, gene) had both plate halves -- too few plates?")
    print(f"{len(de):,} gene-conditions have both plate halves")

    # Rung 0 scores every gene the table carries. --panel-file exists for a later rung asking
    # for a restriction of this ceiling; it is not passed here, and there is deliberately no
    # top-variance fallback -- silently scoring 2,000 genes when the design says "every gene"
    # is the class of error this task exists to undo.
    if args.panel_file:
        declared = {
            ln.strip() for ln in Path(args.panel_file).read_text().splitlines() if ln.strip()
        }
        panel = declared & set(de["gene_name"].unique())
        print(f"restricted to a supplied panel: {len(panel)} of {len(declared)} genes present")
    else:
        panel = set(de["gene_name"].unique())
        print(f"scoring every gene the table carries: {len(panel)}")

    # --- score, both gene sets -------------------------------------------------------------
    r_all, piv0, piv1 = score_split_half(de, panel, min_genes=args.min_genes)
    if not np.any(np.isfinite(r_all)):
        raise SystemExit("no (line, drug, dose) condition had enough shared genes to score")
    padj = padj_pivot(de, panel).reindex(columns=piv0.columns).loc[piv0.index]
    select = responder_mask(padj, alpha=args.padj_threshold)

    # Everything that still needs the long frame is done HERE, before the null draws, and then
    # the frame is released. At this screen's size it is hundreds of millions of rows and the
    # null step allocates arrays quadratic in the condition count on top of it.
    overlap = responder_overlap_table(de, panel, alpha=args.padj_threshold)
    overlap.to_csv(out_dir / "rung0_responder_overlap.csv", index=False)
    stride = max(1, len(de) // 200_000)
    delta_real = pd.DataFrame({"log2FoldChange": de["lfc0"].to_numpy(dtype=float)[::stride]})
    # Written out, not just held: every figure in this run is drawn from a table a reader can
    # open, and the build figure's fold-change panel was the one exception -- it took an
    # in-memory array that reached no file, so the distribution it showed was uncheckable.
    delta_real.to_csv(out_dir / "rung0_delta_sample.csv.gz", index=False)
    padj_sample = pd.DataFrame({"padj0": padj.to_numpy(dtype=float).ravel()}).dropna()
    padj_sample = padj_sample.sample(min(200_000, len(padj_sample)), random_state=args.seed)
    padj_sample.to_csv(out_dir / "rung0_padj_sample.csv.gz", index=False)
    del de, padj
    import gc

    gc.collect()
    r_resp = masked_rowwise_pearson(
        piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float), args.min_genes, select=select
    )
    print(
        f"all-gene: {int(np.isfinite(r_all).sum())} conditions scored; "
        f"responder: {int(np.isfinite(r_resp).sum())} conditions, "
        f"{float(select.sum(axis=1).mean()):.0f} genes each on average"
    )

    # --- pool composition, and the equal-halves subset ---------------------------------------
    # Joined on the full condition key, dose included. The pool table is keyed the same way, so
    # a stale two-part key here would silently mark the wrong conditions as equal-half.
    even_by_key = {
        (str(p), str(d), str(x)): bool(v)
        for p, d, x, v in zip(
            pool["patient"], pool["drug"], pool["dose"], pool["n_plates_even"], strict=True
        )
    }
    even_mask = np.array(
        [even_by_key.get((str(p), str(d), str(x)), False) for p, d, x in piv0.index], dtype=bool
    )
    print(f"{int(even_mask.sum())} of {even_mask.size} conditions split into equal halves")

    # --- nulls, both gene sets ---------------------------------------------------------------
    nulls_all = stratified_null_draws(
        piv0, piv1, n_perm=args.n_perm, seed=args.seed, min_genes=args.min_genes
    )
    nulls_resp = stratified_null_draws(
        piv0, piv1, n_perm=args.n_perm, seed=args.seed, min_genes=args.min_genes, select=select
    )
    for label, nulls in (("all", nulls_all), ("responder", nulls_resp)):
        for k, v in nulls.items():
            med = float(np.median(v)) if v.size else float("nan")
            print(f"null[{label:<9} {k:<10}] median r = {med:+.3f} over {v.size} draws")

    summary = {
        "replicate_col": repl,
        "n_genes": len(panel),
        "padj_threshold": args.padj_threshold,
        **summarize(r_all, nulls_all, args.seed, label="all", even_mask=even_mask),
        **summarize(r_resp, nulls_resp, args.seed, label="responder", even_mask=even_mask),
    }
    summary_path = out_dir / "rung0_reliability.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    _write_params_sidecar(
        summary_path,
        args,
        extra={
            "all_n_pairs": summary["all_n_pairs"],
            "responder_n_pairs": summary["responder_n_pairs"],
            "selection_rule": "padj < threshold in at least one of the first plate group's rows",
            "gene_inclusion": "every gene the table carries",
            "drug_inclusion": "every drug with a (line, drug, dose) triple on two or more plates",
            "dose_handling": "held fixed: a condition is a (line, drug, dose) triple and the "
            "split is between that triple's plates",
            "split_rule": SPLIT_RULE,
            "weighting": "every scored triple weighs one; the dose strata and the per-pair "
            "weighting are in rung0_dose_strata.csv",
        },
    )

    # --- evidence tables ---------------------------------------------------------------------
    per_pair = per_pair_table(piv0, piv1, r_all, r_responder=r_resp, select=select)
    per_pair["n_plates_even"] = even_mask
    per_pair.to_csv(out_dir / "rung0_per_pair_r.csv", index=False)
    dose_strata = dose_strata_table(per_pair)
    dose_strata.to_csv(out_dir / "rung0_dose_strata.csv", index=False)

    null_rows = [
        null_draw_table(nulls_all).assign(gene_set="all"),
        null_draw_table(nulls_resp).assign(gene_set="responder"),
    ]
    pd.concat(null_rows, ignore_index=True).to_csv(out_dir / "rung0_null_draws.csv", index=False)

    profiles, profile_index = example_pair_profiles(piv0, piv1, r_all, select=select)
    profiles.to_csv(out_dir / "rung0_example_pair_profiles.csv.gz", index=False)
    profile_index.to_csv(out_dir / "rung0_example_pair_index.csv", index=False)

    terciles = effect_size_tercile_table(per_pair, seed=args.seed)
    terciles.to_csv(out_dir / "rung0_effect_terciles.csv", index=False)

    mde_curve = mde_curve_table(r_all, r_resp, nulls_all, nulls_resp, seed=args.seed)
    mde_curve.to_csv(out_dir / "rung0_mde_curve.csv", index=False)

    leakage = leakage_table(args.min_genes, seed=args.seed)
    leakage.to_csv(out_dir / "rung0_leakage_control.csv", index=False)
    print(f"leakage control: {leakage.to_dict(orient='records')}")

    per_gene = per_gene_reliability(piv0, piv1)
    per_gene.to_csv(out_dir / "rung0_per_gene_reliability.csv", index=False)

    # --- the noise decomposition, from the slices' sums --------------------------------------
    noise = write_noise_outputs(per_condition, noise_sample, args, out_dir)

    # --- figures -------------------------------------------------------------------------------
    from fmharness import figures as fg
    from fmharness.synthetic import (
        noise_sd_for_reliability,
        planted_noise_frame,
        planted_split_half_frame,
    )

    pos = planted_split_half_frame(
        n_genes=2000, noise_sd=noise_sd_for_reliability(0.5, 4), n_responders=400, seed=args.seed
    )
    neg = planted_split_half_frame(n_genes=2000, signal_sd=0.0, noise_sd=1.0, seed=args.seed + 1)
    ctrl_rows = []
    for label, frame in (("positive (planted r = 0.5)", pos), ("negative (no signal)", neg)):
        cpanel = set(frame["gene_name"].unique())
        cr, cp0, cp1 = score_split_half(frame, cpanel, min_genes=args.min_genes)
        ctrl_rows.append(per_pair_table(cp0, cp1, cr).assign(control=label))
    control_per_pair = pd.concat(ctrl_rows, ignore_index=True)
    control_per_pair.to_csv(out_dir / "rung0_control_per_pair.csv", index=False)

    delta_syn = pd.DataFrame({"log2FoldChange": pos["lfc0"].to_numpy(dtype=float)})

    fg.fig_build(pool, delta_real, delta_syn, fig_dir / "01_build.png")
    fg.fig_split(pool, per_pair, fig_dir / "02_split.png")
    fg.fig_select(per_pair, padj_sample, leakage, fig_dir / "03_select.png", overlap=overlap)
    fg.fig_score(
        profiles, profile_index, per_pair, control_per_pair, summary, fig_dir / "04_score.png"
    )
    if not noise.empty:
        # The control pool plants a plate component of variance 0.25 on sampling variance
        # 0.25 at two plates -- the real screen's regime -- so the pooled estimator has a
        # known answer of 0.5 to recover, and the per-gene panel shows what the one-degree-of-
        # freedom spread looks like when the truth is known.
        control_noise = decompose_noise(planted_noise_frame(n_plates=2, seed=args.seed))
        control_noise.to_csv(out_dir / "rung0_control_noise.csv.gz", index=False)
        strata = pd.read_csv(out_dir / "rung0_noise_strata.csv")
        fg.fig_decompose(noise, control_noise, strata, fig_dir / "05_decompose.png")
    fg.fig_null(per_pair, pd.concat(null_rows, ignore_index=True), fig_dir / "06_null.png")
    fg.fig_terciles(terciles, fig_dir / "07_terciles.png")
    fg.fig_power(mde_curve, fig_dir / "08_power.png")
    write_per_gene_figure(per_gene, fig_dir / "09_per_gene_reliability.png")
    fg.fig_dose(per_pair, dose_strata, fig_dir / "10_dose.png")

    write_audit_checksums(out_dir)

    print("\n=== rung 0: the reliability of the assay ===")
    for k, v in summary.items():
        print(f"  {k:42s} {v}")
    print(
        f"\nall-gene    r = {summary['all_splithalf_mean_r']:.3f} "
        f"(Spearman-Brown {summary['all_spearman_brown_full']:.3f}) "
        f"over {summary['all_n_pairs']} conditions"
        f"\nresponder   r = {summary['responder_splithalf_mean_r']:.3f} "
        f"(Spearman-Brown {summary['responder_spearman_brown_full']:.3f}) "
        f"over {summary['responder_n_pairs']} conditions"
        f"\nwrote {summary_path}"
    )


def _write_split_tables(cache_dir: Path, out_dir: Path) -> None:
    """The split assignment and the pool description, as the committed CSVs."""
    _, pool, _ = read_assignment(cache_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pd.read_parquet(cache_dir / ASSIGNMENT_FILE).to_csv(
        out_dir / "rung0_split_assignment.csv", index=False
    )
    pool.to_csv(out_dir / "rung0_pool_description.csv", index=False)


if __name__ == "__main__":
    main()
