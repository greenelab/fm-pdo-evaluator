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
    de = con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug, gene_name,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 0) AS lfc0,
                   avg(log2FoldChange) FILTER (WHERE hash({chosen}) % 2 = 1) AS lfc1,
                   {padj} AS padj0,
                   {padj1} AS padj1
            FROM read_parquet(?)
            WHERE {chosen} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug, gene_name""",
        [paths, *drug_params],
    ).df()
    return de, chosen


def build_noise_frame(
    paths: list[str],
    target_names: list[str] | None,
    repl_col: str | None,
    tmp: Path,
    memory_limit: str = "36GB",
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
    dose effect pooled into the across-plate variance would be reported as plate noise, and the
    claim the decomposition exists to make would be wrong.
    """
    con = _connect(tmp, memory_limit)
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
    dose_sel = f"{dose} AS dose," if dose else "CAST(NULL AS DOUBLE) AS dose,"
    dose_grp = f", {dose}" if dose else ""
    base = "avg(baseMean)" if "baseMean" in cols else "CAST(NULL AS DOUBLE)"
    where, drug_params = _drug_predicate(target_names)
    return con.execute(
        f"""SELECT Cell_ID_DepMap AS patient, drug, {dose_sel} gene_name,
                   var_samp(log2FoldChange) AS var_lfc,
                   avg(lfcSE * lfcSE) AS mean_se2,
                   count(DISTINCT {repl}) AS n_plates,
                   {base} AS base_mean,
                   avg(log2FoldChange) AS mean_lfc
            FROM read_parquet(?)
            WHERE {repl} IS NOT NULL{where}
            GROUP BY Cell_ID_DepMap, drug{dose_grp}, gene_name
            HAVING count(DISTINCT {repl}) >= 2""",
        [paths, *drug_params],
    ).df()


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
    p0 = padj_pivot(de, panel)
    p1 = de[de["gene_name"].isin(panel)].pivot_table(
        index=["patient", "drug"], columns="gene_name", values="padj1"
    )
    cols = p0.columns.intersection(p1.columns)
    rows = p0.index.intersection(p1.index)
    a = p0.loc[rows, cols].to_numpy(dtype=float)
    b = p1.loc[rows, cols].to_numpy(dtype=float)
    first = np.isfinite(a) & (a < alpha)
    second = np.isfinite(b) & (b < alpha)
    both = first & second
    union = (first | second).sum(axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        jaccard = np.where(union > 0, both.sum(axis=1) / union, np.nan)
    return pd.DataFrame(
        {
            "patient": rows.get_level_values(0),
            "drug": rows.get_level_values(1),
            "n_first": first.sum(axis=1),
            "n_second": second.sum(axis=1),
            "n_both": both.sum(axis=1),
            "jaccard": np.round(jaccard, 4),
        }
    )


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


def frame_cache_key(paths: list[str], names: list[str] | None, replicate_col: str | None) -> str:
    """Content key for a built split-half frame: the inputs that determine it, hashed.

    Naming the cache after its inputs is what makes it safe -- a different pool, drug set or
    replicate column resolves to a different file rather than silently reusing the wrong frame.
    """
    drug_key = ["<all drugs>"] if not names else sorted(names)
    payload = "\n".join([*sorted(paths), "--", *drug_key, "--", str(replicate_col)])
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


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
        "--skip-noise",
        action="store_true",
        help="skip the noise decomposition's second full scan (local smoke runs only)",
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
    even_by_key = dict(
        zip(zip(pool["patient"], pool["drug"], strict=True), pool["n_plates_even"], strict=True)
    )
    even_mask = np.array([bool(even_by_key.get((p, d), False)) for p, d in piv0.index], dtype=bool)
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

    profiles, profile_index = example_pair_profiles(piv0, piv1, r_all)
    profiles.to_csv(out_dir / "rung0_example_pair_profiles.csv.gz", index=False)
    profile_index.to_csv(out_dir / "rung0_example_pair_index.csv", index=False)

    terciles = effect_size_tercile_table(piv0, piv1, r_all, seed=args.seed)
    terciles.to_csv(out_dir / "rung0_effect_terciles.csv", index=False)

    mde_curve = mde_curve_table(r_all, r_resp, nulls_all, nulls_resp, seed=args.seed)
    mde_curve.to_csv(out_dir / "rung0_mde_curve.csv", index=False)

    overlap = responder_overlap_table(de, panel, alpha=args.padj_threshold)
    overlap.to_csv(out_dir / "rung0_responder_overlap.csv", index=False)

    leakage = leakage_table(args.min_genes, seed=args.seed)
    leakage.to_csv(out_dir / "rung0_leakage_control.csv", index=False)
    print(f"leakage control: {leakage.to_dict(orient='records')}")

    per_gene = per_gene_reliability(piv0, piv1)
    per_gene.to_csv(out_dir / "rung0_per_gene_reliability.csv", index=False)
    pool.to_csv(out_dir / "rung0_pool_description.csv", index=False)

    padj_sample = pd.DataFrame({"padj0": padj.to_numpy(dtype=float).ravel()}).dropna()
    padj_sample = padj_sample.sample(min(200_000, len(padj_sample)), random_state=args.seed)
    padj_sample.to_csv(out_dir / "rung0_padj_sample.csv.gz", index=False)

    # --- the noise decomposition -------------------------------------------------------------
    noise_summary: dict[str, float] = {}
    if args.skip_noise:
        print("skipping the noise decomposition (--skip-noise)")
        noise = pd.DataFrame()
    else:
        noise = decompose_noise(
            build_noise_frame(paths, names, args.replicate_col, local.parent / "duckdb_tmp")
        )
        frac = noise["between_plate_fraction"].to_numpy(dtype=float)
        frac = frac[np.isfinite(frac)]
        noise_summary = {
            "n_gene_conditions": len(noise),
            "between_plate_fraction_mean": round(float(np.mean(frac)), 4),
            "between_plate_fraction_median": round(float(np.median(frac)), 4),
            "sigma2_plate_mean": round(
                float(np.mean(noise["sigma2_plate"].to_numpy(dtype=float))), 5
            ),
            "mean_se2_mean": round(float(np.mean(noise["mean_se2"].to_numpy(dtype=float))), 5),
            "frac_plate_dominated": round(float(np.mean(frac > 0.5)), 4),
        }
        noise_path = out_dir / "rung0_noise_decomposition.csv"
        pd.DataFrame([noise_summary]).to_csv(noise_path, index=False)
        _write_params_sidecar(noise_path, args, extra=noise_summary)
        # The per-gene rows are large; committed compressed, and the summary above is what a
        # promoted claim is read from.
        noise.to_csv(out_dir / "rung0_noise_per_gene.csv.gz", index=False)
        print(f"noise: {noise_summary}")

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

    delta_real = pd.DataFrame(
        {"log2FoldChange": de["lfc0"].to_numpy(dtype=float)[:: max(1, len(de) // 200_000)]}
    )
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
