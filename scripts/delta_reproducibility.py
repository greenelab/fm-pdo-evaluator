"""Delta reproducibility ceiling: how noisy is the real Tahoe delta itself?

This is rung 0 of the ladder in docs/SPEC.md: a reference every other rung reads its score
against. Split each (line, drug) pair's replicate plates into two halves, aggregate each
half's delta, and correlate the two halves over the declared gene panel. That split-half
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

#: Bumped whenever the built frame's schema or grouping changes, so a cached frame from an
#: earlier definition resolves to a different key instead of being silently reused. The
#: dose-pooled frame and the dose-fixed frame have the same inputs and different meanings.
FRAME_SCHEMA = "v2-dose-fixed"

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


def _gene_partition(n_parts: int, part: int) -> str:
    """A SQL predicate selecting one slice of the genes, or nothing when there is one slice.

    Partitioning by ``hash(gene_name)`` rather than by file or by row is what makes the slices
    safe to aggregate independently: ``gene_name`` is part of every group key the noise
    decomposition uses, so a group's rows all land in the same slice and per-slice sums combine
    into the whole exactly. Splitting by shard would not have that property -- a gene's rows are
    spread across shards, and per-shard variances would each be computed on a fragment.
    """
    return "" if n_parts <= 1 else f" AND hash(gene_name) % {int(n_parts)} = {int(part)}"


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


def _drug_predicate(target_names: list[str] | None) -> tuple[str, list[object]]:
    """The drug filter, and the parameters it needs — empty when every drug is admitted.

    Rung 0 measures at the assay's full extent: every drug with plates enough to split. The
    superseded rung passed a 32-compound list derived from inputs no longer on any branch, so
    "no drug file" has to mean *no predicate*, not a predicate matching nothing.
    """
    if not target_names:
        return "", []
    return " AND drug IN (SELECT unnest(?))", [list(target_names)]


def build_split_half_frame(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
) -> tuple[pd.DataFrame, str]:
    """Per (line, drug, gene), the mean log2FoldChange in each of two plate halves, via DuckDB.

    Splits plates by ``hash(repl_col) % 2`` (deterministic, no RNG) and aggregates each half
    IN-ENGINE, so raw rows are never materialized. Returns the long frame plus the chosen
    replicate column.

    ``target_names`` of ``None`` or empty admits every drug in the pool.

    Only (line, drug, gene) groups with a fold change in BOTH halves are returned. That is
    exactly what every caller's ``dropna(subset=["lfc0", "lfc1"])`` does next, moved into the
    engine because it is not a small filter here: DESeq2 could not test 59 percent of this
    screen's rows (``baseMean`` zero, so a null fold change), and doing the drop in pandas means
    materialising 1.42 billion rows to keep a fraction of them. The semantics are unchanged --
    a group kept by one is kept by the other.

    The frame also carries ``padj0``: the MINIMUM Benjamini-Hochberg adjusted p-value over the
    FIRST group's (plate, dose) rows. The minimum is what the selection rule asks for -- a gene
    is a responder when the first group called it differentially expressed in at least one of
    its rows -- and it is deliberately one-sided: nothing here aggregates the second group's
    adjusted p-values, because selecting on the half a correlation is scored against inflates
    that correlation by winner's curse. ``min`` skips nulls, so a gene DESeq2 could not test
    comes back null and falls out by the same finiteness rule that governs the fold changes.
    """
    con = _connect(tmp, memory_limit)
    schema = con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [paths]).df()
    cols = list(schema["column_name"])
    print(f"DE columns: {cols}")
    candidates = ([repl_col] if repl_col else []) + list(REPL_CANDIDATES)
    chosen = next((c for c in candidates if c in cols), None)
    if chosen is None:
        raise SystemExit(f"no replicate column in {cols}; pass --replicate-col")
    print(f"splitting plates by hash({chosen}) % 2")
    where, drug_params = _drug_predicate(target_names)
    # FILTER attaches only to an aggregate call, so the no-padj fallback has to replace the whole
    # expression rather than just the function -- `CAST(NULL AS DOUBLE) FILTER (...)` is a parse
    # error, not a null column.
    if "padj" in cols:
        padj = f"min(padj) FILTER (WHERE hash({chosen}) % 2 = 0)"
        # padj1 is the SECOND group's minimum, and it exists for exactly one purpose: the
        # overlap diagnostic the design's select step declares, which shows how much the two
        # groups' responder sets agree. It is never an input to selection -- `responder_mask`
        # takes only padj0 and has no parameter that could admit this column, and a control
        # asserts that. Carrying it makes the diagnostic possible; using it for selection would
        # be the winner's curse the one-sided rule exists to avoid.
        padj1 = f"min(padj) FILTER (WHERE hash({chosen}) % 2 = 1)"
    else:
        padj = "CAST(NULL AS DOUBLE)"
        padj1 = "CAST(NULL AS DOUBLE)"
        print("no padj column in the pool: responder selection is unavailable for this run")
    dose = next((c for c in DOSE_CANDIDATES if c in cols), None)
    if dose is None:
        raise SystemExit(
            f"no dose column in {cols}. Dose is part of the condition key: on this screen it is "
            "confounded with plate, so pooling it would compare different doses across the two "
            "halves. Pass a pool that carries one, or change the design."
        )
    de = con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug, {dose} AS dose, gene_name,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 0) AS lfc0,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 1) AS lfc1,
                   {padj} AS padj0,
                   {padj1} AS padj1
            FROM read_parquet(?)
            WHERE {chosen} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug, {dose}, gene_name
            HAVING count(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 0) > 0
               AND count(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 1) > 0""",
        [paths, *drug_params],
    )
    return _compact_df(de), chosen


def _noise_select(repl: str, dose: str | None, base: str) -> tuple[str, str]:
    """The per (line, drug, dose, gene) select list and GROUP BY the decomposition rests on."""
    dose_sel = f"{dose} AS dose," if dose else "CAST(NULL AS DOUBLE) AS dose,"
    dose_grp = f", {dose}" if dose else ""
    sel = f"""SELECT Cell_ID_DepMap AS patient, drug, {dose_sel} gene_name,
                   var_samp(log2FoldChange) AS var_lfc,
                   avg(lfcSE * lfcSE) AS mean_se2,
                   count(DISTINCT {repl}) AS n_plates,
                   {base} AS base_mean,
                   avg(log2FoldChange) AS mean_lfc,
                   min(padj) FILTER (WHERE hash({repl}) % 2 = 0) AS padj0"""
    grp = (
        f"GROUP BY Cell_ID_DepMap, drug{dose_grp}, gene_name\n"
        f"            HAVING count(DISTINCT {repl}) >= 2"
    )
    return sel, grp


def _noise_setup(
    paths: list[str], repl_col: str | None, tmp: Path, memory_limit: str, threads: int | None
) -> tuple[object, str, str, str]:
    """The connection and the SQL pieces every noise query shares."""
    con = _connect(tmp, memory_limit, threads)
    cols = list(
        con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [paths]).df()["column_name"]
    )
    candidates = ([repl_col] if repl_col else []) + list(REPL_CANDIDATES)
    repl = next((c for c in candidates if c in cols), None)
    if repl is None:
        raise SystemExit(f"no replicate column in {cols}; pass --replicate-col")
    if "lfcSE" not in cols:
        raise SystemExit("the pool carries no lfcSE column; the noise decomposition needs it")
    dose = next((c for c in DOSE_CANDIDATES if c in cols), None)
    if dose is None:
        print(
            "WARNING: no dose column in the pool. Grouping by (line, drug, gene) only, so a dose "
            "effect would be charged to plate noise -- the decomposition is not interpretable."
        )
    base = "avg(baseMean)" if "baseMean" in cols else "CAST(NULL AS DOUBLE)"
    sel, grp = _noise_select(repl, dose, base)
    return con, repl, sel, grp


def noise_slice(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    n_parts: int = 1,
    part: int = 0,
    threads: int | None = None,
) -> pd.DataFrame:
    """One slice of the per (line, drug, dose, gene) noise table, decomposed.

    Reading the table is what costs: a scan of the 83 GB pool takes about twenty minutes, and
    the earlier design ran three separate aggregations over it -- overall figures, per-condition
    figures, and a sample -- each re-scanning every slice. That is 48 scans, sixteen hours, and
    it timed out at six. Every one of those aggregations is a summary of the SAME per
    gene-condition rows, so this returns the rows for one slice and the caller computes all
    three from them in memory. One scan per slice instead of three.

    The slice is small even though the table is not: about 11 million rows at sixteen slices,
    because the two-plate requirement drops 87% of the groups. It is the pre-filter group table,
    not the result, that does not fit -- which is why the slicing is needed at all.
    """
    con, repl, sel, grp = _noise_setup(paths, repl_col, tmp, memory_limit, threads)
    where, drug_params = _drug_predicate(target_names)
    frame = _compact_df(
        con.execute(
            f"""{sel}
            FROM read_parquet(?)
            WHERE {repl} IS NOT NULL{where}{_gene_partition(n_parts, part)}
            {grp}""",
            [paths, *drug_params],
        )
    )
    return decompose_noise(frame)


def noise_partials(noise: pd.DataFrame, alpha: float) -> tuple[dict[str, float], pd.DataFrame]:
    """Everything one slice contributes: overall sums, and per-condition sums.

    Sums rather than means, because sums are what combine across slices. The slice key is part
    of the group key, so the slices partition the gene-conditions exactly -- no overlap, no
    omission -- and adding them gives the same answer a single pass would.
    """
    frac = noise["between_plate_fraction"].to_numpy(dtype=float)
    ok = np.isfinite(frac)
    d = noise.loc[ok].copy()
    fr = frac[ok]
    overall = {
        "n": float(ok.sum()),
        "frac": float(fr.sum()),
        "sigma2": float(d["sigma2_plate"].to_numpy(dtype=float).sum()),
        "se2": float(d["mean_se2"].to_numpy(dtype=float).sum()),
        "dominated": float((fr > 0.5).sum()),
    }
    is_resp = np.isfinite(d["padj0"].to_numpy(dtype=float)) & (
        d["padj0"].to_numpy(dtype=float) < alpha
    )
    var = d["var_lfc"].to_numpy(dtype=float)
    per = (
        pd.DataFrame(
            {
                "patient": d["patient"].to_numpy(),
                "drug": d["drug"].to_numpy(),
                "dose": d["dose"].to_numpy(),
                "n_gene_doses": 1,
                "s_frac": fr,
                "s_var": var,
                "s_sigma2": d["sigma2_plate"].to_numpy(dtype=float),
                "s_se2": d["mean_se2"].to_numpy(dtype=float),
                "s_var_resp": np.where(is_resp, var, np.nan),
                "s_var_non": np.where(is_resp, np.nan, var),
                "s_frac_resp": np.where(is_resp, fr, np.nan),
                "n_responder_gene_doses": is_resp.astype(int),
                "n_nonresponder_gene_doses": (~is_resp).astype(int),
            }
        )
        .groupby(list(CONDITION_KEYS), observed=True, sort=True, dropna=False)
        .sum(min_count=0)
    )
    return overall, per.reset_index()


def combine_noise_partials(
    totals: dict[str, float], per_condition: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Turn accumulated sums into the reported means, once, at the end."""
    n = max(totals["n"], 1.0)
    overall = pd.DataFrame(
        [
            {
                "n_gene_conditions": int(totals["n"]),
                "between_plate_fraction_mean": totals["frac"] / n,
                "sigma2_plate_mean": totals["sigma2"] / n,
                "mean_se2_mean": totals["se2"] / n,
                "frac_plate_dominated": totals["dominated"] / n,
            }
        ]
    )
    c = (
        per_condition.groupby(list(CONDITION_KEYS), observed=True, sort=True, dropna=False)
        .sum()
        .reset_index()
    )
    n_c = c["n_gene_doses"].to_numpy(dtype=float)
    n_r = c["n_responder_gene_doses"].to_numpy(dtype=float)
    n_n = c["n_nonresponder_gene_doses"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        by_condition = pd.DataFrame(
            {
                "patient": c["patient"],
                "drug": c["drug"],
                "dose": c["dose"],
                "n_gene_doses": c["n_gene_doses"].astype(int),
                "between_plate_fraction_mean": c["s_frac"].to_numpy(dtype=float) / n_c,
                "var_lfc_mean": c["s_var"].to_numpy(dtype=float) / n_c,
                "sigma2_plate_mean": c["s_sigma2"].to_numpy(dtype=float) / n_c,
                "mean_se2_mean": c["s_se2"].to_numpy(dtype=float) / n_c,
                "var_lfc_mean_responders": np.where(
                    n_r > 0, c["s_var_resp"].to_numpy(dtype=float) / n_r, np.nan
                ),
                "var_lfc_mean_nonresponders": np.where(
                    n_n > 0, c["s_var_non"].to_numpy(dtype=float) / n_n, np.nan
                ),
                "between_plate_fraction_mean_responders": np.where(
                    n_r > 0, c["s_frac_resp"].to_numpy(dtype=float) / n_r, np.nan
                ),
                "n_responder_gene_doses": c["n_responder_gene_doses"].astype(int),
            }
        )
    return overall, by_condition


def noise_strata_from_sample(noise: pd.DataFrame) -> pd.DataFrame:
    """The between-plate share by expression and response-size quartile, from the sample.

    Computed on the bounded sample rather than on every gene-condition, and that is a deliberate
    trade rather than a shortcut. Assigning quartiles over 1.4 billion rows means sorting them
    twice, which is exactly what exhausted the engine's memory when this ran in SQL. A quartile
    MEAN from a two-million-row sample has a standard error in the fourth decimal place -- far
    inside the precision anyone reads a variance share to -- while the headline figures the
    design promotes are still computed over every row.
    """
    if noise.empty:
        return pd.DataFrame(
            columns=[
                "expression_quartile",
                "response_quartile",
                "n",
                "between_plate_fraction_mean",
                "base_mean_mean",
                "abs_mean_lfc_mean",
            ]
        )
    d = noise.dropna(subset=["between_plate_fraction"]).copy()
    d["abs_lfc"] = d["mean_lfc"].abs()
    d["expression_quartile"] = pd.qcut(
        d["base_mean"].rank(method="first"), 4, labels=[1, 2, 3, 4]
    ).astype(int)
    d["response_quartile"] = pd.qcut(
        d["abs_lfc"].rank(method="first"), 4, labels=[1, 2, 3, 4]
    ).astype(int)
    out = (
        d.groupby(["expression_quartile", "response_quartile"], observed=True)
        .agg(
            n=("between_plate_fraction", "size"),
            between_plate_fraction_mean=("between_plate_fraction", "mean"),
            base_mean_mean=("base_mean", "mean"),
            abs_mean_lfc_mean=("abs_lfc", "mean"),
        )
        .reset_index()
    )
    return out


def noise_by_condition(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    alpha: float = 0.05,
    n_parts: int = 1,
    threads: int | None = None,
) -> pd.DataFrame:
    """Per (line, drug), the noise in that condition -- and the same split by responder status.

    Two of the design's four hypotheses are about noise per condition rather than noise overall:
    that noise is higher where the correlations are lower, and that it is higher in the responder
    genes than across all genes. Neither can be answered from a single screen-wide mean, so this
    returns the per-condition aggregate a reader can join to the per-condition correlations, with
    the responder split computed from the same first-half adjusted p-value the selection rule
    uses.

    Partitioned by gene like ``noise_aggregate``, and combined the same way: each slice returns
    per-condition SUMS and counts, which add across slices to the whole. Means are taken once at
    the end, so the result is identical to a single pass over everything.
    """
    con, repl, sel, grp = _noise_setup(paths, repl_col, tmp, memory_limit, threads)
    where, drug_params = _drug_predicate(target_names)
    partials: list[pd.DataFrame] = []
    for part in range(max(1, n_parts)):
        partials.append(
            con.execute(
                f"""
                WITH g AS ({sel}
                    FROM read_parquet(?)
                    WHERE {repl} IS NOT NULL{where}{_gene_partition(n_parts, part)}
                    {grp}),
                d AS (SELECT patient, drug, dose, var_lfc, mean_se2,
                             greatest(var_lfc - mean_se2, 0.0) AS sigma2_plate,
                             CASE WHEN var_lfc > 0
                                  THEN greatest(var_lfc - mean_se2, 0.0) / var_lfc END
                               AS between_plate_fraction,
                             (padj0 IS NOT NULL AND padj0 < {alpha}) AS is_responder
                      FROM g)
                SELECT patient, drug, dose,
                       count(*) AS n_gene_doses,
                       sum(between_plate_fraction) AS s_frac,
                       sum(var_lfc) AS s_var,
                       sum(sigma2_plate) AS s_sigma2,
                       sum(mean_se2) AS s_se2,
                       sum(CASE WHEN is_responder THEN var_lfc END) AS s_var_resp,
                       sum(CASE WHEN NOT is_responder THEN var_lfc END) AS s_var_non,
                       sum(CASE WHEN is_responder THEN between_plate_fraction END)
                         AS s_frac_resp,
                       count(*) FILTER (WHERE is_responder) AS n_responder_gene_doses,
                       count(*) FILTER (WHERE NOT is_responder) AS n_nonresponder_gene_doses
                FROM d WHERE between_plate_fraction IS NOT NULL
                GROUP BY patient, drug, dose""",
                [paths, *drug_params],
            ).df()
        )
    combined = (
        pd.concat(partials, ignore_index=True)
        .groupby(list(CONDITION_KEYS), observed=True, sort=True)
        .sum()
        .reset_index()
    )
    n = combined["n_gene_doses"].to_numpy(dtype=float)
    n_resp = combined["n_responder_gene_doses"].to_numpy(dtype=float)
    n_non = combined["n_nonresponder_gene_doses"].to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        out = pd.DataFrame(
            {
                "patient": combined["patient"],
                "drug": combined["drug"],
                "dose": combined["dose"],
                "n_gene_doses": combined["n_gene_doses"].astype(int),
                "between_plate_fraction_mean": combined["s_frac"].to_numpy(dtype=float) / n,
                "var_lfc_mean": combined["s_var"].to_numpy(dtype=float) / n,
                "sigma2_plate_mean": combined["s_sigma2"].to_numpy(dtype=float) / n,
                "mean_se2_mean": combined["s_se2"].to_numpy(dtype=float) / n,
                "var_lfc_mean_responders": np.where(
                    n_resp > 0, combined["s_var_resp"].to_numpy(dtype=float) / n_resp, np.nan
                ),
                "var_lfc_mean_nonresponders": np.where(
                    n_non > 0, combined["s_var_non"].to_numpy(dtype=float) / n_non, np.nan
                ),
                "between_plate_fraction_mean_responders": np.where(
                    n_resp > 0, combined["s_frac_resp"].to_numpy(dtype=float) / n_resp, np.nan
                ),
                "n_responder_gene_doses": combined["n_responder_gene_doses"].astype(int),
            }
        )
    return out


def build_noise_frame(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
    sample_rows: int | None = None,
    n_parts: int = 1,
    threads: int | None = None,
) -> pd.DataFrame:
    """Per (line, drug, dose, gene) with at least two plates: the ingredients of the noise split.

    ``lfcSE`` is the standard error of ONE plate's treated-versus-control contrast -- cell
    sampling at that row's cell counts. It cannot see plate-to-plate variation, which is the
    noise a model trained on other material actually meets and the noise the split-half
    measures. Under plate offsets plus independent sampling error, the sample variance of
    ``log2FoldChange`` across a condition's plates has expectation ``sigma2_plate +
    mean(lfcSE^2)`` -- exactly, for any set of per-plate standard errors -- so those two
    quantities are what this returns and ``decompose_noise`` subtracts.

    **Dose is a grouping key here, not pooled.** The reliabilities pool over dose because a
    condition means "this drug at this screen's dose design"; this aggregation cannot, because a
    dose effect pooled into the across-plate variance would be reported as plate noise.

    With ``sample_rows`` the result is a bounded sample: an equal share drawn from each gene
    slice, which makes it a stratified sample by gene rather than a plain reservoir. That is the
    right shape for what it feeds -- figures and per-gene checks -- and every promoted number is
    computed over all rows by ``noise_aggregate`` regardless.
    """
    con, repl, sel, grp = _noise_setup(paths, repl_col, tmp, memory_limit, threads)
    where, drug_params = _drug_predicate(target_names)
    parts = max(1, n_parts)
    per_part = None if sample_rows is None else max(1, int(sample_rows) // parts)
    frames: list[pd.DataFrame] = []
    for part in range(parts):
        limit = (
            f"\n        USING SAMPLE reservoir({per_part} ROWS) REPEATABLE (0)" if per_part else ""
        )
        frames.append(
            _compact_df(
                con.execute(
                    f"""{sel}
                    FROM read_parquet(?)
                    WHERE {repl} IS NOT NULL{where}{_gene_partition(n_parts, part)}
                    {grp}{limit}""",
                    [paths, *drug_params],
                )
            )
        )
    return pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]


def decompose_noise(noise: pd.DataFrame) -> pd.DataFrame:
    """Split each delta's variance into its between-plate and within-plate parts.

        sigma2_plate = var_samp(log2FoldChange across plates) - mean(lfcSE^2), floored at zero

    The floor is what makes the result a variance rather than a difference: the estimator is
    unbiased, so on data with no plate effect it lands either side of zero and half its values
    would be negative unaided.

    ``between_plate_fraction`` is ``sigma2_plate / var_lfc`` -- the share of a delta's
    replicate variance that plate effects account for. Null where ``var_lfc`` is zero or
    missing, rather than an arbitrary zero or one, since a delta with no variance at all has no
    share to report.
    """
    out = noise.copy()
    var = out["var_lfc"].to_numpy(dtype=float)
    se2 = out["mean_se2"].to_numpy(dtype=float)
    sigma2 = np.maximum(var - se2, 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(var > 0, sigma2 / var, np.nan)
    out["sigma2_plate"] = sigma2
    out["between_plate_fraction"] = frac
    return out


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
    """The result's own data: one row per candidate (line, drug) condition.

    The headline row summarizes these values; committing them makes the summary re-derivable
    anywhere -- mean, median, quartiles, positive fraction, and the effect-size terciles all
    recompute from this table without cluster access. Columns: the split-half r whose mean is
    the declared statistic, the pair's effect size (mean |delta| over genes finite in both
    halves -- the same quantity ``effect_size_terciles`` stratifies on), and the shared
    finite-gene count. Rows scored NaN (fewer than ``min_genes`` shared genes) are kept,
    honestly NaN.
    """
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    n = ok.sum(axis=1)
    mean_abs = np.where(ok, np.abs(a + b) / 2.0, 0.0).sum(axis=1) / np.maximum(n, 1)
    out = pd.DataFrame(
        {
            "patient": piv0.index.get_level_values(0),
            "drug": piv0.index.get_level_values(1),
            "dose": piv0.index.get_level_values(2),
            "n_genes_scored": n,
            "mean_abs_delta": np.round(mean_abs, 4),
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

    any_pair: two different pairs (continuity with the archived lineage's first run).
    diff_drug: different line AND drug -- the generic-structure floor the ceiling clears.
    same_drug: same drug, different line -- the line-specificity floor.

    With ``select``, a draw pairing condition *i*'s first group against condition *j*'s second
    group scores over **row i's** selected genes -- the row whose first group is used, which is
    the row the selection rule would have read. Using row *j*'s mask, or the union of the two,
    would apply a different rule to the null than to the observed value and the comparison would
    stop being like for like. The finiteness rule then intersects with row *j*'s second group,
    exactly as it does for a matched pair.
    """
    lines = piv0.index.get_level_values(0).to_numpy(dtype=str)
    drugs = piv0.index.get_level_values(1).to_numpy(dtype=str)
    n = len(piv0)
    ii, jj = np.divmod(np.arange(n * n), n)
    off = ii != jj
    ii, jj = ii[off], jj[off]
    same_drug = drugs[ii] == drugs[jj]
    same_line = lines[ii] == lines[jj]
    strata = {
        "any_pair": np.ones(ii.size, dtype=bool),
        "diff_drug": ~same_drug & ~same_line,
        "same_drug": same_drug & ~same_line,
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
        ("same_drug_mismatch", (drugs == drugs[anchor]) & (lines != lines[anchor])),
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


def effect_size_terciles(piv0: pd.DataFrame, piv1: pd.DataFrame, r: np.ndarray) -> dict[str, float]:
    """Split-half mean r within terciles of per-pair effect size (mean |delta|).

    The empirical positive control: an assay that cannot find more reproducibility where
    there is more signal is broken. Tercile 1 = smallest effects.
    """
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    mean_abs = np.where(ok, np.abs(a + b) / 2.0, 0.0).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)
    finite = np.isfinite(r)
    edges = np.quantile(mean_abs[finite], [1 / 3, 2 / 3])
    out: dict[str, float] = {}
    for t in (1, 2, 3):
        lo = -np.inf if t == 1 else edges[t - 2]
        hi = np.inf if t == 3 else edges[t - 1]
        sel = finite & (mean_abs > lo) & (mean_abs <= hi)
        out[f"splithalf_mean_r_tercile{t}"] = round(float(np.mean(r[sel])), 3)
    return out


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


DOSE_CANDIDATES = ("dose", "Dose", "drug_dose", "concentration", "dose_uM")


def pool_description(
    paths: list[str], target_names: list[str] | None, repl: str, tmp: Path
) -> pd.DataFrame:
    """Measured composition of the consumed pool (design: 'measured not asserted'):
    per (line, drug) the replicate-row count, distinct plates per half, and dose levels
    when a dose column exists.

    ``n_plates_even`` marks the conditions whose plate count splits into two equal groups.
    Spearman-Brown assumes equal halves, and three quarters of this screen's conditions split
    one plate against two; the corrected value is reported again over these conditions, where
    the correction is exact, and the gap between the two is the size of the assumption.
    """
    con = _connect(tmp)
    cols = list(
        con.execute("DESCRIBE SELECT * FROM read_parquet(?) LIMIT 0", [paths]).df()["column_name"]
    )
    dose = next((c for c in DOSE_CANDIDATES if c in cols), None)
    # The share of this condition's rows DESeq2 could not test (baseMean zero, so a null fold
    # change). Measured per condition because it is large and uneven here -- 59 percent of the
    # screen's rows overall -- and a reader comparing gene counts across conditions needs to
    # see it rather than infer it.
    untestable = (
        "avg(CASE WHEN baseMean = 0 THEN 1.0 ELSE 0.0 END)" if "baseMean" in cols else "NULL"
    )
    where, drug_params = _drug_predicate(target_names)
    if dose is None:
        raise SystemExit("no dose column in the pool; dose is part of the condition key")
    return con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug, {dose} AS dose,
                   count(*) AS n_rows,
                   count(DISTINCT {repl}) AS n_plates,
                   count(DISTINCT {repl}) % 2 = 0 AS n_plates_even,
                   count(DISTINCT {repl}) FILTER (WHERE hash({repl}) % 2 = 0) AS n_plates_half0,
                   count(DISTINCT {repl}) FILTER (WHERE hash({repl}) % 2 = 1) AS n_plates_half1,
                   {untestable} AS frac_untestable
            FROM read_parquet(?)
            WHERE {repl} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug, {dose} ORDER BY patient, drug, dose""",
        [paths, *drug_params],
    ).df()


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

    Naming the cache after its inputs is what makes it safe -- a different pool, drug set or
    replicate column resolves to a different file rather than silently reusing the wrong frame.
    """
    drug_key = ["<all drugs>"] if not names else sorted(names)
    payload = "\n".join(
        [*sorted(paths), "--", *drug_key, "--", str(replicate_col), "--", FRAME_SCHEMA]
    )
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


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


def _build_or_load_frame(
    paths: list[str], names: list[str] | None, args: argparse.Namespace, local: Path
) -> tuple[pd.DataFrame, str]:
    """The split-half frame, from cache when one matches these inputs.

    Building it scans every DE shard through DuckDB and dominates the run (~40 minutes on the
    real pool); everything after it takes about a minute. Adding an output or changing a figure
    therefore does not need the scan repeated, so ``--frame-cache`` keeps the built frame beside
    the data and reuses it when the inputs hash the same.
    """
    cache_dir = Path(args.frame_cache) if args.frame_cache else None
    cache = None
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        key = frame_cache_key(paths, names, args.replicate_col)
        cache = cache_dir / f"split_half_{key}.parquet"
        sidecar = cache.with_suffix(".json")
        if cache.exists() and sidecar.exists():
            de = _normalise_keys(pd.read_parquet(cache))
            repl = str(json.loads(sidecar.read_text())["replicate_col"])
            print(f"loaded the split-half frame from {cache} ({len(de):,} rows, replicate {repl})")
            return de, repl

    # getattr, not attribute access: the frame-cache control constructs a minimal Namespace to
    # exercise the cache key, and a helper that demands every CLI flag makes that test carry the
    # whole argument list rather than the two fields it is about.
    de, repl = build_split_half_frame(
        paths,
        names,
        args.replicate_col,
        local.parent / "duckdb_tmp",
        getattr(args, "duckdb_memory", "36GB"),
    )
    if cache is not None:
        de.to_parquet(cache, index=False)
        cache.with_suffix(".json").write_text(json.dumps({"replicate_col": repl}) + "\n")
        print(f"cached the split-half frame at {cache}")
    return de, repl


def effect_size_tercile_table(
    piv0: pd.DataFrame, piv1: pd.DataFrame, r: np.ndarray, n_boot: int = 2000, seed: int = 0
) -> pd.DataFrame:
    """The empirical control as a table, with an interval on each tercile mean.

    ``effect_size_terciles`` returns three bare numbers, which cannot say whether a rise across
    them is real or within noise. This carries the same three means with a bootstrap interval
    and a count, so the figure can show the rise with its uncertainty and a reader can tell the
    two cases apart.
    """
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    ok = np.isfinite(a) & np.isfinite(b)
    mean_abs = np.where(ok, np.abs(a + b) / 2.0, 0.0).sum(axis=1) / np.maximum(ok.sum(axis=1), 1)
    finite = np.isfinite(r)
    edges = np.quantile(mean_abs[finite], [1 / 3, 2 / 3])
    rng = np.random.default_rng(seed)
    rows = []
    for t in (1, 2, 3):
        lo = -np.inf if t == 1 else edges[t - 2]
        hi = np.inf if t == 3 else edges[t - 1]
        sel = finite & (mean_abs > lo) & (mean_abs <= hi)
        vals = r[sel]
        boot = np.array(
            [np.mean(rng.choice(vals, size=vals.size, replace=True)) for _ in range(n_boot)]
        )
        rows.append(
            {
                "tercile": t,
                "n": int(vals.size),
                "mean_r": round(float(np.mean(vals)), 4),
                "ci_lo": round(float(np.quantile(boot, 0.025)), 4),
                "ci_hi": round(float(np.quantile(boot, 0.975)), 4),
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


def run_noise(
    paths: list[str], names: list[str] | None, args: argparse.Namespace, local: Path, out_dir: Path
) -> pd.DataFrame:
    """The whole noise decomposition, written to ``out_dir``, returning the figure sample.

    One scan of the pool per slice, not three. Reading the table is the entire cost -- about
    twenty minutes for 83 GB -- and the three things reported (overall figures, per-condition
    figures, a sample for the figures) are three summaries of the same per gene-condition rows.
    An earlier version ran them as three separate aggregations, which meant 48 scans and a
    six-hour timeout after finishing only the first sixteen.

    Each slice's contribution is cached as it completes. A slice is twenty minutes of work whose
    result is a few megabytes; losing all of them because the job ran out of wall clock is a bad
    trade, and with the cache a rerun continues instead of restarting.
    """
    tmp = local.parent / "duckdb_tmp"
    parts = max(1, args.noise_partitions)
    cache_dir = Path(args.noise_cache) if args.noise_cache else out_dir / "noise_slices"
    cache_dir.mkdir(parents=True, exist_ok=True)
    per_slice_sample = max(1, args.noise_sample_rows // parts)

    totals = {"n": 0.0, "frac": 0.0, "sigma2": 0.0, "se2": 0.0, "dominated": 0.0}
    per_frames: list[pd.DataFrame] = []
    samples: list[pd.DataFrame] = []
    for part in range(parts):
        tag = f"{FRAME_SCHEMA}_{parts}_{part}"
        o_path = cache_dir / f"overall_{tag}.json"
        c_path = cache_dir / f"per_condition_{tag}.parquet"
        s_path = cache_dir / f"sample_{tag}.parquet"
        if o_path.exists() and c_path.exists() and s_path.exists():
            overall_part = json.loads(o_path.read_text())
            per_frames.append(pd.read_parquet(c_path))
            samples.append(pd.read_parquet(s_path))
            print(f"  noise slice {part + 1}/{parts}: reused cached partial")
        else:
            noise = noise_slice(
                paths,
                names,
                args.replicate_col,
                tmp,
                args.duckdb_memory,
                n_parts=parts,
                part=part,
                threads=args.duckdb_threads,
            )
            overall_part, per = noise_partials(noise, args.padj_threshold)
            sample = noise.sample(
                min(per_slice_sample, len(noise)), random_state=args.seed
            ).reset_index(drop=True)
            o_path.write_text(json.dumps(overall_part) + "\n")
            per.to_parquet(c_path, index=False)
            sample.to_parquet(s_path, index=False)
            per_frames.append(per)
            samples.append(sample)
            print(f"  noise slice {part + 1}/{parts}: {int(overall_part['n']):,} gene-conditions")
        for key in totals:
            totals[key] += float(overall_part[key])

    overall, by_condition = combine_noise_partials(totals, pd.concat(per_frames, ignore_index=True))
    noise_summary = {k: float(v) for k, v in overall.iloc[0].items()}
    by_condition.round(6).to_csv(out_dir / "rung0_noise_by_condition.csv", index=False)
    print(f"per-condition noise: {len(by_condition):,} dose-conditions")

    # The design says the between-plate share is aggregated "over genes within a condition and
    # over conditions". The figure above is a flat mean over gene-conditions, which weights a
    # condition by how many of its genes were testable -- and testability varies a lot here,
    # since DESeq2 could not test most rows. Both are reported; where they disagree, the gap is
    # uneven gene coverage rather than anything about plates.
    per_cond = by_condition["between_plate_fraction_mean"].to_numpy(dtype=float)
    per_cond = per_cond[np.isfinite(per_cond)]
    noise_summary["between_plate_fraction_mean_over_conditions"] = round(
        float(np.mean(per_cond)), 5
    )
    noise_summary["n_conditions_decomposed"] = int(per_cond.size)

    noise = pd.concat(samples, ignore_index=True)
    noise.to_csv(out_dir / "rung0_noise_per_gene.csv.gz", index=False)
    noise_strata_from_sample(noise).round(5).to_csv(out_dir / "rung0_noise_strata.csv", index=False)

    # The median comes from the stratified sample. An exact median over 175 million values means
    # sorting them; a median from a sample stratified by gene is precise well past the third
    # decimal it is reported to. The mean beside it is exact over every row.
    frac_sample = noise["between_plate_fraction"].to_numpy(dtype=float)
    frac_sample = frac_sample[np.isfinite(frac_sample)]
    noise_summary["between_plate_fraction_median_sampled"] = (
        round(float(np.median(frac_sample)), 5) if frac_sample.size else float("nan")
    )

    noise_path = out_dir / "rung0_noise_decomposition.csv"
    pd.DataFrame([noise_summary]).round(5).to_csv(noise_path, index=False)
    _write_params_sidecar(noise_path, args, extra=noise_summary)
    print(
        "between-plate share: "
        f"{noise_summary['between_plate_fraction_mean']:.4f} over gene-conditions, "
        f"{noise_summary['between_plate_fraction_mean_over_conditions']:.4f} over conditions"
    )
    print(f"noise sample for figures: {len(noise):,} gene-conditions")
    return noise


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
        "--frame-cache",
        default=None,
        help="directory holding the built split-half frame, keyed by a hash of the inputs. "
        "The build scans every DE shard and dominates the run; with a cache, adding an output "
        "or changing a figure reruns in about a minute instead of repeating the scan.",
    )
    ap.add_argument(
        "--duckdb-memory",
        default="36GB",
        help="DuckDB's memory_limit. The engine spills to --local-dir's parent when it runs "
        "out, so a value below the job's allocation trades speed for safety; set it well under "
        "the SBATCH --mem so the returned pandas frames still have room.",
    )
    ap.add_argument(
        "--noise-sample-rows",
        type=int,
        default=2_000_000,
        help="rows of the per-gene noise table to keep for the figures and the row-wise "
        "identity check. Every reported number is aggregated over all rows in the engine; this "
        "only bounds what comes back.",
    )
    ap.add_argument(
        "--skip-noise",
        action="store_true",
        help="do not compute the noise decomposition; read what an earlier --only-noise pass "
        "wrote into --out-dir, so the figures still show it",
    )
    ap.add_argument(
        "--noise-cache",
        default=None,
        help="directory holding each noise slice's partial result. A slice is about twenty "
        "minutes of work and a few megabytes of output; caching them means a rerun continues "
        "rather than restarting, and a job that runs out of wall clock loses one slice.",
    )
    ap.add_argument(
        "--noise-partitions",
        type=int,
        default=1,
        help="split the noise decomposition into this many slices of the GENES and combine the "
        "slices. The combination is exact -- the slice key is part of the group key, so every "
        "gene-condition lands in exactly one slice and per-slice sums add. Needed at full "
        "extent, where the single-slice group table did not fit at 140 GB.",
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
        "--only-noise",
        action="store_true",
        help="compute ONLY the noise decomposition and exit. It groups four billion rows into "
        "roughly one and a half billion and needs most of the machine to do it; sharing an "
        "allocation with the reliability pass exhausted memory twice, once on each side.",
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
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    if args.only_noise:
        run_noise(paths, names, args, local, out_dir)
        return

    de, repl = _build_or_load_frame(paths, names, args, local)
    n_rows_built = len(de)
    de = de.dropna(subset=["lfc0", "lfc1"])
    if de.empty:
        raise SystemExit("no (line, drug, gene) had both plate halves -- too few plates per pair?")
    print(f"built {n_rows_built:,} rows; {len(de):,} have both plate halves")

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
        raise SystemExit("no (line, drug) pair had enough shared genes to score")
    padj = padj_pivot(de, panel).reindex(columns=piv0.columns).loc[piv0.index]
    select = responder_mask(padj, alpha=args.padj_threshold)

    # Everything that still needs the long frame is done HERE, before the null draws, and then
    # the frame is released. At this screen's size it is 451 million rows and the null step
    # allocates arrays quadratic in the condition count on top of it -- holding both put the
    # first attempt at 178 GB of a 200 GB allocation. Nothing about the numbers changes; the
    # frame is simply not alive during the step that needs the room.
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
    pool = pool_description(paths, names, repl, local.parent / "duckdb_tmp")
    # Joined on the full condition key, dose included. The pool description is keyed the same
    # way, so a stale two-part key here would silently mark the wrong conditions as equal-half.
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
            "drug_inclusion": "every drug with at least two distinct plates",
            "dose_handling": "pooled for the reliabilities, held fixed for the decomposition",
        },
    )

    # --- evidence tables ---------------------------------------------------------------------
    per_pair = per_pair_table(piv0, piv1, r_all, r_responder=r_resp, select=select)
    per_pair["n_plates_even"] = even_mask
    per_pair.to_csv(out_dir / "rung0_per_pair_r.csv", index=False)

    null_rows = [
        null_draw_table(nulls_all).assign(gene_set="all"),
        null_draw_table(nulls_resp).assign(gene_set="responder"),
    ]
    pd.concat(null_rows, ignore_index=True).to_csv(out_dir / "rung0_null_draws.csv", index=False)

    profiles, profile_index = example_pair_profiles(piv0, piv1, r_all, select=select)
    profiles.to_csv(out_dir / "rung0_example_pair_profiles.csv.gz", index=False)
    profile_index.to_csv(out_dir / "rung0_example_pair_index.csv", index=False)

    terciles = effect_size_tercile_table(piv0, piv1, r_all, seed=args.seed)
    terciles.to_csv(out_dir / "rung0_effect_terciles.csv", index=False)

    mde_curve = mde_curve_table(r_all, r_resp, nulls_all, nulls_resp, seed=args.seed)
    mde_curve.to_csv(out_dir / "rung0_mde_curve.csv", index=False)

    leakage = leakage_table(args.min_genes, seed=args.seed)
    leakage.to_csv(out_dir / "rung0_leakage_control.csv", index=False)
    print(f"leakage control: {leakage.to_dict(orient='records')}")

    per_gene = per_gene_reliability(piv0, piv1)
    per_gene.to_csv(out_dir / "rung0_per_gene_reliability.csv", index=False)
    pool.to_csv(out_dir / "rung0_pool_description.csv", index=False)

    # --- the noise decomposition -------------------------------------------------------------
    sample_path = out_dir / "rung0_noise_per_gene.csv.gz"
    if args.skip_noise:
        # The noise phase runs as its own process; read what it wrote, so the figures below are
        # drawn from the same tables a reader opens rather than from a second computation.
        noise = pd.read_csv(sample_path) if sample_path.exists() else pd.DataFrame()
        if noise.empty:
            print("no noise tables present: the decompose figure will be skipped")
        else:
            print(f"read the noise sample written by the --only-noise pass: {len(noise):,} rows")
    else:
        noise = run_noise(paths, names, args, local, out_dir)

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
        control_noise = decompose_noise(planted_noise_frame(seed=args.seed))
        control_noise.to_csv(out_dir / "rung0_control_noise.csv.gz", index=False)
        fg.fig_decompose(noise, control_noise, fig_dir / "05_decompose.png")
    fg.fig_null(per_pair, pd.concat(null_rows, ignore_index=True), fig_dir / "06_null.png")
    fg.fig_terciles(terciles, fig_dir / "07_terciles.png")
    fg.fig_power(mde_curve, fig_dir / "08_power.png")
    write_per_gene_figure(per_gene, fig_dir / "09_per_gene_reliability.png")

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


if __name__ == "__main__":
    main()
