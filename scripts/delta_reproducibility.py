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


def _connect(tmp: Path, memory_limit: str = "36GB"):
    """A DuckDB connection configured to spill to ``tmp`` rather than exhaust memory."""
    import duckdb  # type: ignore  # Alpine-only

    tmp.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect()
    con.execute(f"SET temp_directory='{tmp}'")
    con.execute(f"SET memory_limit='{memory_limit}'")
    con.execute("SET preserve_insertion_order=false")
    return con


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
    else:
        padj = "CAST(NULL AS DOUBLE)"
        print("no padj column in the pool: responder selection is unavailable for this run")
    de = con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug, gene_name,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 0) AS lfc0,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 1) AS lfc1,
                   {padj} AS padj0
            FROM read_parquet(?)
            WHERE {chosen} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug, gene_name""",
        [paths, *drug_params],
    ).df()
    return de, chosen


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
    d = de[de["gene_name"].isin(panel)]
    piv0 = d.pivot_table(index=["patient", "drug"], columns="gene_name", values="lfc0")
    piv1 = d.pivot_table(index=["patient", "drug"], columns="gene_name", values="lfc1")
    common = piv0.index.intersection(piv1.index)
    piv0, piv1 = piv0.loc[common], piv1.loc[common]
    # pivot_table drops all-NaN columns per half, so the two halves can carry different gene
    # sets even after the index intersection above -- callers that skip main()'s dropna could
    # otherwise silently correlate misaligned genes whenever the column COUNTS happen to match.
    cols = piv0.columns.intersection(piv1.columns)
    piv0, piv1 = piv0[cols], piv1[cols]
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
    d = de[de["gene_name"].isin(panel)]
    return d.pivot_table(index=["patient", "drug"], columns="gene_name", values="padj0")


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


def per_pair_table(piv0: pd.DataFrame, piv1: pd.DataFrame, r: np.ndarray) -> pd.DataFrame:
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
    return pd.DataFrame(
        {
            "patient": piv0.index.get_level_values(0),
            "drug": piv0.index.get_level_values(1),
            "n_genes_scored": n,
            "mean_abs_delta": np.round(mean_abs, 4),
            "r": np.round(r, 4),
        }
    )


def stratified_null_draws(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    n_perm: int = 500,
    seed: int = 0,
    min_genes: int = 50,
) -> dict[str, np.ndarray]:
    """Mismatched-pair null correlations per stratum.

    any_pair: two different pairs (continuity with the archived lineage's first run).
    diff_drug: different line AND drug -- the generic-structure floor the ceiling clears.
    same_drug: same drug, different line -- the line-specificity floor.
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
        r = masked_rowwise_pearson(a[ii[pick]], b[jj[pick]], min_genes)
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
    order = np.flatnonzero(np.isfinite(r))[np.argsort(r[np.isfinite(r)])]
    if order.size == 0:
        index_cols = ["example_id", "kind", "patient0", "drug0", "patient1", "drug1"]
        index_cols += ["n_genes_full", "r_full", "n_genes_shown", "r_shown"]
        return pd.DataFrame(columns=["example_id", "gene", "lfc0", "lfc1"]), pd.DataFrame(
            columns=index_cols
        )

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
        frames.append(
            pd.DataFrame(
                {
                    "example_id": example_id,
                    "gene": genes[shown],
                    "lfc0": np.round(x, 4),
                    "lfc1": np.round(y, 4),
                }
            )
        )
        rows.append(
            {
                "example_id": example_id,
                "kind": kind,
                "patient0": lines[i],
                "drug0": drugs[i],
                "patient1": lines[j],
                "drug1": drugs[j],
                "n_genes_full": int(shared.size),
                "r_full": round(r_full, 4),
                "n_genes_shown": int(shown.size),
                "r_shown": round(r_shown, 4),
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


def summarize(r: np.ndarray, nulls: dict[str, np.ndarray], seed: int = 0) -> dict:
    """The headline row: mean-over-pairs Pearson (the declared statistic), its nulls,
    p-values from the bootstrapped null aggregate, and the MDEs (SPEC rule 4)."""
    from fmharness.statistics import (
        bootstrap_aggregate_pvalue,
        minimum_detectable_aggregate,
        spearman_brown,
    )

    r = r[np.isfinite(r)]
    mean = float(np.mean(r))
    nl = nulls["diff_drug"] if nulls["diff_drug"].size else nulls["any_pair"]
    p_boot, ci_lo, ci_hi = bootstrap_aggregate_pvalue(mean, nl, r.size, seed=seed)
    p_same = bootstrap_aggregate_pvalue(mean, nulls["same_drug"], r.size, seed=seed)[0]
    sb = spearman_brown(mean) if mean > -1 else float("nan")
    return {
        "n_pairs": int(r.size),
        "splithalf_mean_r": round(mean, 3),
        "splithalf_median_r": round(float(np.median(r)), 3),
        "splithalf_q1_r": round(float(np.quantile(r, 0.25)), 3),
        "splithalf_q3_r": round(float(np.quantile(r, 0.75)), 3),
        "spearman_brown_full": round(sb, 3),
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
    dose_expr = f"count(DISTINCT {dose})" if dose else "NULL"
    where, drug_params = _drug_predicate(target_names)
    return con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug,
                   count(*) AS n_rows,
                   count(DISTINCT {repl}) AS n_plates,
                   count(DISTINCT {repl}) % 2 = 0 AS n_plates_even,
                   count(DISTINCT {repl}) FILTER (WHERE hash({repl}) % 2 = 0) AS n_plates_half0,
                   count(DISTINCT {repl}) FILTER (WHERE hash({repl}) % 2 = 1) AS n_plates_half1,
                   {dose_expr} AS n_dose_levels
            FROM read_parquet(?)
            WHERE {repl} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug ORDER BY patient, drug""",
        [paths, *drug_params],
    ).df()


def write_figure(r: np.ndarray, nulls: dict[str, np.ndarray], out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(7, 4.5))
    bins = np.linspace(-0.3, 0.8, 56)
    ax.hist(r[np.isfinite(r)], bins=bins, density=True, alpha=0.65, label="matched pairs")
    ax.hist(nulls["diff_drug"], bins=bins, density=True, alpha=0.45, label="diff-drug null")
    ax.hist(nulls["same_drug"], bins=bins, density=True, alpha=0.45, label="same-drug null")
    ax.axvline(float(np.nanmean(r)), color="k", lw=1.5, label="mean (headline)")
    ax.set_xlabel("split-half Pearson r per (line, drug) pair")
    ax.set_ylabel("density")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


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


def frame_cache_key(paths: list[str], names: list[str], replicate_col: str | None) -> str:
    """Content key for a built split-half frame: the inputs that determine it, hashed.

    Naming the cache after its inputs is what makes it safe -- a different pool, drug set or
    replicate column resolves to a different file rather than silently reusing the wrong frame.
    """
    payload = "\n".join([*sorted(paths), "--", *sorted(names), "--", str(replicate_col)])
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def _build_or_load_frame(
    paths: list[str], names: list[str], args: argparse.Namespace, local: Path
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
            de = pd.read_parquet(cache)
            repl = str(json.loads(sidecar.read_text())["replicate_col"])
            print(f"loaded the split-half frame from {cache} ({len(de):,} rows, replicate {repl})")
            return de, repl

    de, repl = build_split_half_frame(paths, names, args.replicate_col, local.parent / "duckdb_tmp")
    if cache is not None:
        de.to_parquet(cache, index=False)
        cache.with_suffix(".json").write_text(json.dumps({"replicate_col": repl}) + "\n")
        print(f"cached the split-half frame at {cache}")
    return de, repl


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local-dir", required=True, help="dir with the Tahoe DE parquet (on scratch)")
    ap.add_argument("--drugs-cid-file", default="data/static/tahoe_target_cids.txt")
    ap.add_argument(
        "--drug-names-file",
        default=None,
        help="one Tahoe drug name per line; bypasses the HuggingFace name lookup so fixtures "
        "and offline runs need no `datasets` import.",
    )
    ap.add_argument("--replicate-col", default=None, help="plate/replicate column (auto-detected)")
    ap.add_argument(
        "--n-hvg",
        type=int,
        default=2000,
        help="top HVGs by variance, used only without --panel-file",
    )
    ap.add_argument("--min-genes", type=int, default=50, help="min shared genes to score a pair")
    ap.add_argument(
        "--panel-file",
        default=None,
        help="one gene per line; pins the ceiling to the SAME panel it will be a denominator "
        "for. Without it the ceiling is top-HVG and not comparable to a panel-scored rung.",
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
    args = ap.parse_args()
    repo = Path(__file__).resolve().parent.parent

    local = Path(args.local_dir)
    local = local if local.is_absolute() else repo / local
    paths = sorted(str(p) for p in local.rglob("*.parquet") if DE in str(p))
    if not paths:
        raise SystemExit(f"no {DE} parquet under {local}")
    if args.drug_names_file:
        names = sorted(
            {ln.strip() for ln in Path(args.drug_names_file).read_text().splitlines() if ln.strip()}
        )
    else:
        cid_file = Path(args.drugs_cid_file)
        cid_file = cid_file if cid_file.is_absolute() else repo / cid_file
        names = _target_names(repo, cid_file)
    print(f"{len(names)} target drugs; reading {len(paths)} DE parquet files ...")

    out_dir = Path(args.out_dir) if Path(args.out_dir).is_absolute() else repo / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    de, repl = _build_or_load_frame(paths, names, args, local)
    de = de.dropna(subset=["lfc0", "lfc1"])
    if de.empty:
        raise SystemExit("no (line, drug, gene) had both plate halves -- too few plates per pair?")
    de["mean"] = (de["lfc0"].to_numpy() + de["lfc1"].to_numpy()) / 2.0

    # The gene set has to match whatever this ceiling will be a DENOMINATOR for. Scoring the
    # ceiling on top-2000 HVG while rung 1 scores on the 14,121-gene common panel makes
    # "fraction of achievable" a ratio between two different measurements -- the same class of
    # error as the panel bug itself. --panel-file pins it to the panel actually used.
    if args.panel_file:
        panel = {ln.strip() for ln in Path(args.panel_file).read_text().splitlines() if ln.strip()}
        hvg = panel & set(de["gene_name"].unique())
        print(
            f"scoring the ceiling on the supplied panel: {len(hvg)} of {len(panel)} genes present"
        )
    else:
        gene_var = de.groupby("gene_name")["mean"].var()
        hvg = set(gene_var.sort_values(ascending=False).index[: args.n_hvg])
        print(f"scoring the ceiling on top-{args.n_hvg} HVG (no --panel-file given)")

    # per (line, drug): split-half Pearson r between the two plate halves over the HVG genes.
    # `r` stays aligned with piv0/piv1's rows (NaNs and all) -- Task 4's tercile/per-gene
    # diagnostics need that alignment; `r_fin` is the finite-only view for aggregate stats.
    r, piv0, piv1 = score_split_half(de, hvg, min_genes=args.min_genes)
    r_fin = r[np.isfinite(r)]
    if r_fin.size == 0:
        raise SystemExit("no (line, drug) pair had enough shared HVG genes to score")

    # NEGATIVE CONTROL, STRATIFIED. A split-half correlation has a nonzero floor because genes
    # share structure whether or not two halves come from the same perturbation. But a single
    # "mismatched pair" null conflates two very different floors, and the first run showed why:
    # it drew two random pairs, so whenever they happened to share a DRUG the correlation was
    # high -- drug effects dominate the delta -- and the null came back at median 0.139 with 23%
    # of draws exceeding the observed 0.299. That null is inflated by same-drug matches and
    # cannot be read as a floor for reproducibility.
    #
    # Three strata, the same distinction design.md's "Why three null strata" draws:
    #   any_pair      two random pairs. Mixed; reported only for continuity with the first run.
    #   diff_drug     different line AND different drug -- the floor from generic gene structure.
    #                 This is the one the CEILING must clear to be a ceiling at all.
    #   same_drug     same drug, different line -- the floor for LINE specificity. A split-half
    #                 above this says the pair's delta is reproducible beyond its drug's effect.
    nulls = stratified_null_draws(
        piv0, piv1, n_perm=args.n_perm, seed=args.seed, min_genes=args.min_genes
    )
    for k, v in nulls.items():
        med_k = float(np.median(v)) if v.size else float("nan")
        print(f"null[{k:<10}] median r = {med_k:+.3f} over {v.size} draws")

    summary = {
        "replicate_col": repl,
        "n_genes": len(hvg),
        **summarize(r, nulls, args.seed),
        **effect_size_terciles(piv0, piv1, r),
    }
    summary_path = out_dir / "rung0_delta_reproducibility.csv"
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    _write_params_sidecar(summary_path, args, extra={"n_pairs": summary["n_pairs"]})

    per_pair_table(piv0, piv1, r).to_csv(out_dir / "rung0_per_pair_r.csv", index=False)
    null_draw_table(nulls).to_csv(out_dir / "rung0_null_draws.csv", index=False)
    profiles, profile_index = example_pair_profiles(piv0, piv1, r)
    # gzipped: every shared gene, ~700 KB instead of ~2.4 MB; pandas reads it by extension
    profiles.to_csv(out_dir / "rung0_example_pair_profiles.csv.gz", index=False)
    profile_index.to_csv(out_dir / "rung0_example_pair_index.csv", index=False)

    per_gene = per_gene_reliability(piv0, piv1)
    per_gene.to_csv(out_dir / "rung0_per_gene_reliability.csv", index=False)
    write_per_gene_figure(per_gene, out_dir / "rung0_per_gene_reliability.png")
    pool_description(paths, names, repl, local.parent / "duckdb_tmp").to_csv(
        out_dir / "rung0_pool_description.csv", index=False
    )
    write_figure(r, nulls, out_dir / "rung0_ceiling.png")

    print("\n=== delta reproducibility ceiling (real Tahoe delta, plate split-half) ===")
    for k, v in summary.items():
        print(f"  {k:22s} {v}")
    print(
        f"\nrung-0 ceiling: split-half mean r = {summary['splithalf_mean_r']:.3f}, "
        f"Spearman-Brown full-data = {summary['spearman_brown_full']:.3f}.\nwrote {summary_path}"
    )


if __name__ == "__main__":
    main()
