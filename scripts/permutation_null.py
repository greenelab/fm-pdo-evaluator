"""Exact permutation null for rung 0's split-half reliabilities.

Rung 0's promoted p-values (`scripts/delta_reproducibility.py`, `summarize`) bootstrap the
observed mean split-half correlation against a mismatched-pair null pool, but that pool's
draws are not exchangeable: `stratified_null_draws` reuses each (line, drug) pair's
half-profiles across many mismatched-pair comparisons, so treating them as an i.i.d. pool
understates their true dependence. `verification.md`'s "Write-up caveat" bounded that
exposure theoretically (roughly 100 bootstrap standard errors of headroom, so the dependence
would need to inflate the null's variance ~3,000-fold to change the conclusion) but did not
measure it directly. This script carries the dependence by construction instead of assuming
it away: it permutes the pairing between the two half-profile pivots -- a permutation with no
fixed points, so every mismatched draw uses a real profile but never the correct partner --
and builds the null distribution of the MEAN mismatched correlation directly from n_perm such
permutations. The ratio of that permutation-null variance to the i.i.d.-pool bootstrap's
assumed variance is the design effect the exchangeable-pool treatment ignores.

That any-pair permutation validates the POOLED aggregate, but the promoted p-values are
per-stratum (`p_vs_null` vs a diff-drug mismatched-pair pool, `p_vs_same_drug` vs a same-drug
pool that clusters over ~32 drugs) and an any-pair permutation mixes both mismatch types
freely, carrying neither stratum's dependence specifically. `stratified_permutation_null`
builds two more permutation nulls that do: `sample_within_drug_permutation` (mismatched pairs
stay inside one drug -- the `same_drug` stratum) and `sample_cross_permutation` (mismatched
pairs change both line and drug -- the `diff_drug` stratum).

Reuses the measurement core in `scripts/delta_reproducibility.py` (`build_split_half_frame`,
`score_split_half`, `stratified_null_draws`, `masked_rowwise_pearson`) rather than
reimplementing it, loaded the same way `tests/test_rung0_controls.py` loads it.

Deliberate divergence from `delta_reproducibility.py`: without `--panel-file`, this script
scores on all genes present in the data rather than falling back to a top-HVG subset (there is
no `--n-hvg` here). A verification tool should score on exactly the panel the ceiling it is
checking used, never a different gene set of its own choosing -- the sbatch job always passes
`--panel-file`, so the fallback exists only for ad hoc local runs.

  python scripts/permutation_null.py --local-dir /scratch/alpine/$USER/tahoe_pseudobulk_de \\
      --drug-names-file <a file of Tahoe drug names, one per line> \\
      --panel-file results/rung1_panel/common_panel.txt --out-dir rung0_outputs
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from itertools import pairwise
from pathlib import Path

import numpy as np
import pandas as pd

_SPEC = importlib.util.spec_from_file_location(
    "delta_reproducibility",
    Path(__file__).resolve().parents[1] / "scripts" / "delta_reproducibility.py",
)
assert _SPEC is not None and _SPEC.loader is not None
dr = importlib.util.module_from_spec(_SPEC)
sys.modules["delta_reproducibility"] = dr
_SPEC.loader.exec_module(dr)


def sample_permutation(rng: np.random.Generator, n: int, max_tries: int = 1000) -> np.ndarray:
    """A permutation of ``range(n)`` with no fixed points, by rejection sampling.

    A uniform random permutation is a permutation with probability ~1/e (inclusion-exclusion),
    so ``max_tries`` is generous headroom for the expected ~e attempts, not a real limit --
    except at n=1, where no permutation exists and every attempt fails by construction.
    """
    identity = np.arange(n)
    for _ in range(max_tries):
        perm = rng.permutation(n)
        if not np.any(perm == identity):
            return perm
    raise RuntimeError(f"failed to sample a permutation of size {n} within {max_tries} tries")


def sample_within_drug_permutation(rng: np.random.Generator, drugs: np.ndarray) -> np.ndarray:
    """A permutation sigma with ``sigma(i) != i`` and ``drugs[sigma(i)] == drugs[i]`` for every
    row in a drug group of size >= 2 -- a per-drug permutation, built by rejection-sampling
    `sample_permutation` independently within each drug's row indices.

    A drug group of size 1 has no permutation (the same degeneracy `sample_permutation` hits at
    n=1): those rows are left mapped to themselves rather than raising, and are excluded from
    any aggregate over the result -- the mask of included rows is directly computable from
    ``drugs``' own per-row group counts (``>= 2``), so it is not returned separately here.
    """
    n = len(drugs)
    sigma = np.arange(n)
    order = np.argsort(drugs, kind="stable")
    sorted_drugs = drugs[order]
    boundaries = np.flatnonzero(np.r_[True, sorted_drugs[1:] != sorted_drugs[:-1], True])
    for start, end in pairwise(boundaries):
        idx = order[start:end]
        if idx.size < 2:
            continue  # singleton drug group: no permutation exists, row stays a fixed point
        sigma[idx] = idx[sample_permutation(rng, idx.size)]
    return sigma


def sample_cross_permutation(
    rng: np.random.Generator, drugs: np.ndarray, lines: np.ndarray, max_sweeps: int = 200
) -> np.ndarray:
    """A permutation sigma with ``sigma(i) != i``, ``drugs[sigma(i)] != drugs[i]``, AND
    ``lines[sigma(i)] != lines[i]`` for every row -- exactly `stratified_null_draws`'s
    ``diff_drug`` mask (different line and different drug).

    Starts from a uniform random permutation, then repeatedly sweeps the violating positions
    and swaps each with a uniformly random other position (redrawing on a self-swap, which can
    never resolve the violation it was meant to fix), until no violations remain. With the real
    pool (~33 drugs, ~50 lines, ~1,600 rows) violations are sparse -- most random partners
    already differ on both drug and line -- so this converges in a handful of sweeps;
    ``max_sweeps`` bounds a pathological composition (e.g. one drug shared by most rows) so a
    non-converging repair fails loudly instead of looping forever.

    Caveat, MEASURED not assumed: the repair converges to SOME permutation satisfying every
    constraint, but the swap process is a local repair, not `sample_permutation`'s plain
    rejection sampling (which is exactly uniform over permutations by construction) -- and it
    is NOT uniform over the constraint-satisfying set. An empirical probe (brute-force
    enumerating all 448 valid permutations on a small 3-drug x 3-line, 9-row fixture, then
    50,000 draws from this sampler) found a clear departure: chi-square goodness-of-fit against
    uniform gave statistic 2857 on 447 degrees of freedom (expectation ~447 under uniform),
    counts ranging 51-209 against a mean of 112, and a coefficient of variation of counts
    (0.24) roughly 2.5x the ~0.09 a uniform multinomial would produce at this sample size --
    not a borderline result. The exact p-value this feeds is therefore exact only under
    exchangeability of the *draws this sampler actually generates*, not under a formal
    uniform-over-the-constraint-set guarantee; every p_exact this script has reported (real run
    and fixtures alike) sits far from any decision boundary (<=0.007), where a biased sampling
    distribution over permutations sharing the SAME symmetric constraint structure as the
    observed matched pairing is unlikely to flip significance, but this is a real limitation
    of the "exact" framing, not a cosmetic one -- a provably uniform cross-permutation sampler
    (e.g. Metropolis-Hastings with a detailed-balance-respecting proposal) is future work, not
    done here.
    """
    n = len(drugs)
    if n < 2:
        raise ValueError(f"sample_cross_permutation needs at least 2 rows, got {n}")
    sigma = rng.permutation(n)
    idx = np.arange(n)

    def _violations(s: np.ndarray) -> np.ndarray:
        return (s == idx) | (drugs[s] == drugs) | (lines[s] == lines)

    for _ in range(max_sweeps):
        bad_idx = np.flatnonzero(_violations(sigma))
        if bad_idx.size == 0:
            return sigma
        for i in bad_idx.tolist():
            j = i
            while j == i:  # a self-swap is a no-op and can never fix row i's violation
                j = int(rng.integers(0, n))
            sigma[i], sigma[j] = sigma[j], sigma[i]
    # The check above runs at the START of each sweep; the LAST sweep's swaps may have repaired
    # every violation with no further iteration to notice. Check once more before giving up, so
    # a last-pass repair is honored rather than discarded.
    if not _violations(sigma).any():
        return sigma
    raise RuntimeError(
        f"sample_cross_permutation failed to repair every violation within {max_sweeps} sweeps "
        f"(n={n}) -- the diff-drug constraint may be unsatisfiable for this drug/line "
        "composition (e.g. too few distinct drugs or lines)"
    )


def permutation_null(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    r: np.ndarray,
    min_genes: int,
    n_perm: int,
    seed: int,
    *,
    pools: dict[str, np.ndarray] | None = None,
    select: np.ndarray | None = None,
) -> tuple[dict[str, float], np.ndarray]:
    """The exact-permutation null for the split-half mean, carrying the half-profile-sharing
    dependence by construction instead of treating the mismatched-pair pool as exchangeable.

    ``piv0``/``piv1``/``r`` are exactly what `score_split_half` returns: ``r`` may carry NaNs
    (rows below ``min_genes`` shared finite entries, or zero variance), which are dropped here
    together with the corresponding pivot rows before deranging. For each of ``n_perm``
    permutations sigma of the row order, every row ``i`` is paired with half-1's row
    ``sigma(i)`` -- never its own match, by construction -- and the mean Pearson r over that
    mismatched pairing is one null draw. ``design_effect`` is the ratio of the actual sampling
    variance of the mean under permutation resampling to the variance an i.i.d. pool of the
    same size would have (``stratified_null_draws``'s ``any_pair`` stratum, the pool a
    permutation's composition matches): the number the write-up caveat in
    the superseded branch's verification write-up could only bound theoretically.

    ``pools``, if given, must be `stratified_null_draws`'s return value computed on the same
    (``piv0_f``, ``piv1_f``, ``n_perm``, ``seed``, ``min_genes``) this function would use --
    the caller (`main`, sharing one call with `stratified_permutation_null`) passes it to avoid
    recomputing an identical pool twice; omitted, it is computed here exactly as before, so
    every existing call site (and the standalone-callability known-answer tests) is unaffected.

    Returns the summary dict (one row's worth of columns) and the raw array of ``n_perm``
    permutation means.
    """
    finite = np.isfinite(r)
    piv0_f, piv1_f = piv0.loc[finite], piv1.loc[finite]
    a = piv0_f.to_numpy(dtype=float)
    b = piv1_f.to_numpy(dtype=float)
    # The mask is over the caller's full row set, so it is filtered to the finite rows with the
    # pivots -- and it is indexed by the row whose FIRST half is used, never by the permuted
    # partner, so a permuted pairing scores the same genes the selection rule would have read.
    sel = select[finite] if select is not None else None
    n = a.shape[0]
    if n < 2:
        raise ValueError(
            f"permutation_null needs at least 2 finite (line, drug) pairs to permute, got {n}"
        )
    if n_perm < 2:
        raise ValueError(f"permutation_null needs n_perm >= 2, got {n_perm}")
    observed_mean = float(np.mean(r[finite]))

    rng = np.random.default_rng(seed)
    perm_means = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        sigma = sample_permutation(rng, n)
        perm_means[k] = float(
            np.nanmean(dr.masked_rowwise_pearson(a, b[sigma], min_genes, select=sel))
        )

    # An all-NaN row-correlation vector for a single permutation (every mismatched pairing in
    # that draw falling below min_genes shared finite entries) would put a NaN into perm_means.
    # np.mean/np.var/np.std then propagate that NaN into perm_mean_mean, perm_mean_sd, and
    # design_effect, while p_exact's `>=` comparison silently treats NaN as False -- so a
    # broken draw would still produce an ordinary-looking, wrong number instead of an error.
    assert not np.isnan(perm_means).any(), (
        "a permutation produced zero scoreable pairs -- investigate before trusting "
        "design_effect/p_exact"
    )

    if pools is None:
        pools = dr.stratified_null_draws(
            piv0_f, piv1_f, n_perm=n_perm, seed=seed, min_genes=min_genes, select=sel
        )
    pool = pools["any_pair"]
    var_pool = float(np.var(pool, ddof=1))
    se_iid = float(np.sqrt(var_pool / n))
    design_effect = float(np.var(perm_means, ddof=1) / (var_pool / n))
    p_exact = float((1 + np.sum(perm_means >= observed_mean)) / (1 + n_perm))
    perm_mean_mean = float(np.mean(perm_means))
    perm_mean_sd = float(np.std(perm_means, ddof=1))
    z_permutation = (
        (observed_mean - perm_mean_mean) / perm_mean_sd if perm_mean_sd > 0 else float("nan")
    )

    summary = {
        "n_pairs": int(n),
        "n_perm": int(n_perm),
        "observed_mean": round(observed_mean, 4),
        "perm_mean_mean": round(perm_mean_mean, 4),
        "perm_mean_sd": round(perm_mean_sd, 4),
        "p_exact": round(p_exact, 4),
        "se_iid_pool": round(se_iid, 5),
        "design_effect": round(design_effect, 3),
        "z_permutation": round(z_permutation, 2),
    }
    return summary, perm_means


#: Mirrors `fmharness.statistics.bootstrap_aggregate_pvalue`'s `min_null_draws` default: a
#: pool this small cannot support a variance estimate `design_effect`'s denominator needs, so
#: below this many finite draws the design effect is reported as nan rather than a numerically
#: fragile (or silently nan/inf) ratio.
MIN_NULL_DRAWS_FOR_DESIGN_EFFECT = 10


def stratified_permutation_null(
    piv0: pd.DataFrame,
    piv1: pd.DataFrame,
    r: np.ndarray,
    min_genes: int,
    n_perm: int,
    seed: int,
    *,
    pools: dict[str, np.ndarray] | None = None,
    select: np.ndarray | None = None,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    """Stratum-preserving permutation nulls for the promoted PER-STRATUM p-values.

    `permutation_null`'s any-pair permutation validates the pooled aggregate, but the promoted
    p-values in `delta_reproducibility.summarize` (``p_vs_null`` vs the diff-drug pool,
    ``p_vs_same_drug`` vs the same-drug pool, which clusters over ~32 drugs) are per-stratum: an
    any-pair permutation mixes same- and diff-drug mismatches freely and carries neither
    stratum's dependence specifically. This builds two more permutation nulls that DO, by
    construction:

    - ``same_drug``: `sample_within_drug_permutation` -- every draw's mismatched pairing stays
      inside one drug, exactly `stratified_null_draws`'s ``same_drug`` mask (same drug,
      different line). Rows in a singleton drug group have no permutation and are excluded from
      both the null and the observed comparator, which is therefore computed over the SAME row
      subset (``observed_mean_same_drug_rows``, ``n_rows_same_drug``) rather than the full pool.
      ``same_drug_rows_equal_n`` records whether that subset is the FULL finite pool (every drug
      has >= 2 rows): only when it is does ``design_effect_same_drug`` transfer directly to the
      promoted ``p_vs_same_drug`` (computed over all ``n`` rows) -- when it does not, the design
      effect describes only the multi-row subset and that scope must be stated wherever it is
      cited. ``observed_mean_diff_drug_rows`` is recorded alongside it purely so the summary is
      self-describing per stratum; it equals the global ``observed_mean`` by construction, since
      the diff-drug stratum (unlike same-drug) excludes no rows.
    - ``diff_drug``: `sample_cross_permutation` -- every draw changes both line and drug,
      exactly `stratified_null_draws`'s ``diff_drug`` mask, over all finite rows.

    Each stratum's ``design_effect`` is the permutation-null variance of the mean over the
    matching `stratified_null_draws` pool's variance at that stratum's own row count -- the same
    construction as `permutation_null`'s design effect, applied per stratum. If a stratum's pool
    has fewer than `MIN_NULL_DRAWS_FOR_DESIGN_EFFECT` finite draws (mirroring
    `bootstrap_aggregate_pvalue`'s ``min_null_draws`` spirit), that design effect is reported as
    nan, with a printed warning, rather than a ratio built on too few draws to estimate a
    variance from.

    ``pools``, if given, must be `stratified_null_draws`'s return value computed on the same
    (``piv0_f``, ``piv1_f``, ``n_perm``, ``seed``, ``min_genes``) this function would use --
    the caller (`main`, sharing one call with `permutation_null`) passes it to avoid recomputing
    an identical pool twice; omitted, it is computed here exactly as before.

    Returns a summary dict with ``_same_drug``/``_diff_drug``-suffixed keys (plus
    ``n_rows_same_drug``, ``same_drug_rows_equal_n``, and the two per-stratum observed-mean
    keys) and a dict of the two raw perm-mean arrays, keyed ``"same_drug"``/``"diff_drug"``.
    """
    finite = np.isfinite(r)
    piv0_f, piv1_f = piv0.loc[finite], piv1.loc[finite]
    a = piv0_f.to_numpy(dtype=float)
    b = piv1_f.to_numpy(dtype=float)
    # Filtered with the pivots, and indexed by the row whose FIRST half is used -- the same
    # one-sided rule the observed statistic applies (see `permutation_null`).
    sel = select[finite] if select is not None else None
    r_f = r[finite]
    n = a.shape[0]
    if n < 2:
        raise ValueError(
            f"stratified_permutation_null needs at least 2 finite (line, drug) pairs, got {n}"
        )
    if n_perm < 2:
        raise ValueError(f"stratified_permutation_null needs n_perm >= 2, got {n_perm}")

    lines = piv0_f.index.get_level_values(0).to_numpy(dtype=str)
    drugs = piv0_f.index.get_level_values(1).to_numpy(dtype=str)
    observed_mean = float(np.mean(r_f))

    drug_counts = pd.Series(drugs).map(pd.Series(drugs).value_counts())
    multi_mask = (drug_counts >= 2).to_numpy()
    n_multi = int(multi_mask.sum())
    if n_multi < 2:
        raise ValueError(
            "stratified_permutation_null needs at least 2 rows in >=2-row drug groups, got "
            f"{n_multi}"
        )
    observed_mean_multi = float(np.mean(r_f[multi_mask]))

    rng = np.random.default_rng(seed)
    perm_means_same = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        sigma = sample_within_drug_permutation(rng, drugs)
        row_r = dr.masked_rowwise_pearson(
            a[multi_mask],
            b[sigma[multi_mask]],
            min_genes,
            select=None if sel is None else sel[multi_mask],
        )
        perm_means_same[k] = float(np.nanmean(row_r))
    # Same NaN-propagation hazard `permutation_null` guards against: a broken draw (zero
    # scoreable pairs among the multi-row rows) must fail loudly, not silently corrupt
    # perm_mean_mean/perm_mean_sd/design_effect while p_exact's `>=` treats NaN as False.
    assert not np.isnan(perm_means_same).any(), (
        "a within-drug permutation produced zero scoreable pairs among multi-row drug groups"
    )

    perm_means_diff = np.empty(n_perm, dtype=float)
    for k in range(n_perm):
        sigma = sample_cross_permutation(rng, drugs, lines)
        row_r = dr.masked_rowwise_pearson(a, b[sigma], min_genes, select=sel)
        perm_means_diff[k] = float(np.nanmean(row_r))
    assert not np.isnan(perm_means_diff).any(), "a cross permutation produced zero scoreable pairs"

    if pools is None:
        pools = dr.stratified_null_draws(
            piv0_f, piv1_f, n_perm=n_perm, seed=seed, min_genes=min_genes, select=sel
        )

    def _stratum(
        name: str, perm_means: np.ndarray, pool: np.ndarray, n_used: int, observed: float
    ) -> dict:
        pool_finite = pool[np.isfinite(pool)]
        if pool_finite.size < MIN_NULL_DRAWS_FOR_DESIGN_EFFECT:
            print(
                f"WARNING: {name} stratum's stratified_null_draws pool has only "
                f"{pool_finite.size} finite draws (< {MIN_NULL_DRAWS_FOR_DESIGN_EFFECT}); "
                "design_effect set to nan rather than computed from too few draws to estimate "
                "a variance from"
            )
            design_effect = float("nan")
        else:
            var_pool = float(np.var(pool_finite, ddof=1))
            design_effect = float(np.var(perm_means, ddof=1) / (var_pool / n_used))
        p_exact = float((1 + np.sum(perm_means >= observed)) / (1 + n_perm))
        return {
            "perm_mean_mean": round(float(np.mean(perm_means)), 4),
            "perm_mean_sd": round(float(np.std(perm_means, ddof=1)), 4),
            "p_exact": round(p_exact, 4),
            "design_effect": round(design_effect, 3)
            if np.isfinite(design_effect)
            else design_effect,
        }

    same = _stratum("same_drug", perm_means_same, pools["same_drug"], n_multi, observed_mean_multi)
    diff = _stratum("diff_drug", perm_means_diff, pools["diff_drug"], n, observed_mean)

    summary = {
        "n_rows_same_drug": n_multi,
        "same_drug_rows_equal_n": bool(n_multi == n),
        "observed_mean_same_drug_rows": round(observed_mean_multi, 4),
        "observed_mean_diff_drug_rows": round(observed_mean, 4),
        **{f"{k}_same_drug": v for k, v in same.items()},
        **{f"{k}_diff_drug": v for k, v in diff.items()},
    }
    return summary, {"same_drug": perm_means_same, "diff_drug": perm_means_diff}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--local-dir", required=True, help="dir with the Tahoe DE parquet (on scratch)")
    ap.add_argument(
        "--drugs-cid-file",
        default="",
        help="optional PubChem CID list. Empty by default: this check validates the reliability "
        "computed at the assay's FULL extent, so it must score the same pool. The superseded "
        "32-compound list still sits untracked in the Alpine checkout, and a default pointing "
        "at it would have silently validated 33 drugs' dependence against a number computed "
        "over every drug.",
    )
    ap.add_argument(
        "--drug-names-file",
        default=None,
        help="one Tahoe drug name per line; bypasses the HuggingFace name lookup so fixtures "
        "and offline runs need no `datasets` import.",
    )
    ap.add_argument(
        "--panel-file", default=None, help="one gene per line; same panel the ceiling was scored on"
    )
    ap.add_argument(
        "--cache-dir",
        default=None,
        help="directory holding the split assignment and the cached slices -- the same cache "
        "scripts/delta_reproducibility.py writes. With every slice present this job reads the "
        "frame back in minutes; without it, it repeats that job's full scan of every shard.",
    )
    ap.add_argument(
        "--slices",
        type=int,
        default=1,
        help="how many gene slices the cached frame was built in; must match the run that "
        "wrote the cache, or the slices will not be found.",
    )
    ap.add_argument("--duckdb-memory", default="36GB", help="DuckDB's memory_limit for a rebuild")
    ap.add_argument("--duckdb-threads", type=int, default=None, help="DuckDB thread count")
    ap.add_argument("--replicate-col", default=None, help="plate/replicate column (auto-detected)")
    ap.add_argument("--min-genes", type=int, default=50, help="min shared genes to score a pair")
    ap.add_argument("--n-perm", type=int, default=500, help="permutation null draws")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-dir", default="rung0_outputs")
    ap.add_argument(
        "--gene-set",
        choices=("all", "responder"),
        default="all",
        help="which of rung 0's two reliabilities to check. 'responder' scores each condition "
        "over the genes its FIRST plate group called differentially expressed, the same "
        "one-sided rule the observed statistic uses, and writes its outputs under a "
        "_responder suffix so the two runs cannot overwrite each other.",
    )
    args = ap.parse_args()
    repo = Path(__file__).resolve().parent.parent

    local = Path(args.local_dir)
    local = local if local.is_absolute() else repo / local
    paths = sorted(str(p) for p in local.rglob("*.parquet") if dr.DE in str(p))
    if not paths:
        raise SystemExit(f"no {dr.DE} parquet under {local}")
    names = dr.resolve_drug_names(repo, args)
    scope = "all drugs (no drug list given)" if names is None else f"{len(names)} target drugs"
    print(f"{scope}; reading {len(paths)} DE parquet files ...")

    out_dir = Path(args.out_dir) if Path(args.out_dir).is_absolute() else repo / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    args.out_dir = str(out_dir)

    de, repl = dr._build_or_load_frame(paths, names, args, local)
    de = de.dropna(subset=["lfc0", "lfc1"])
    if de.empty:
        raise SystemExit("no (line, drug, dose, gene) had both plate halves -- too few plates?")

    if args.panel_file:
        panel_all = {
            ln.strip() for ln in Path(args.panel_file).read_text().splitlines() if ln.strip()
        }
        panel = panel_all & set(de["gene_name"].unique())
        print(f"scoring on the supplied panel: {len(panel)} of {len(panel_all)} genes present")
    else:
        panel = set(de["gene_name"].unique())
        print(f"no --panel-file given; scoring on all {len(panel)} genes present in the data")

    r, piv0, piv1 = dr.score_split_half(de, panel, min_genes=args.min_genes)
    select = None
    if args.gene_set == "responder":
        padj = dr.padj_pivot(de, panel).reindex(columns=piv0.columns).loc[piv0.index]
        select = dr.responder_mask(padj)
        r, piv0, piv1 = dr.score_split_half(de, panel, min_genes=args.min_genes, select=select)
        print(
            f"responder gene set: {int(select.sum(axis=1).mean())} genes per condition on average"
        )
    if not np.any(np.isfinite(r)):
        raise SystemExit("no (line, drug, dose) condition had enough shared genes to score")

    # Computed once and shared with both nulls below: `permutation_null` and
    # `stratified_permutation_null` each independently call `stratified_null_draws` on the same
    # (finite piv0, finite piv1, n_perm, seed, min_genes) -- identical inputs and a deterministic
    # RNG make the two calls produce identical pools, so computing it once here is pure cost
    # removal, not a behavior change (each function still computes it itself if called without
    # `pools`, as the existing known-answer tests do).
    finite = np.isfinite(r)
    pools = dr.stratified_null_draws(
        piv0.loc[finite],
        piv1.loc[finite],
        n_perm=args.n_perm,
        seed=args.seed,
        min_genes=args.min_genes,
        select=None if select is None else select[finite],
    )
    summary, perm_means = permutation_null(
        piv0, piv1, r, args.min_genes, args.n_perm, args.seed, pools=pools, select=select
    )
    strat_summary, strat_perm_means = stratified_permutation_null(
        piv0, piv1, r, args.min_genes, args.n_perm, args.seed, pools=pools, select=select
    )
    summary_row = {
        "replicate_col": repl,
        "n_genes": len(panel),
        "gene_set": args.gene_set,
        **summary,
        **strat_summary,
    }

    suffix = "" if args.gene_set == "all" else "_responder"
    summary_path = out_dir / f"rung0_permutation_summary{suffix}.csv"
    pd.DataFrame([summary_row]).to_csv(summary_path, index=False)
    dr._write_params_sidecar(summary_path, args, extra={"n_pairs": summary["n_pairs"]})

    pd.DataFrame({"perm_mean": perm_means}).to_csv(
        out_dir / f"rung0_permutation_perm_means{suffix}.csv", index=False
    )
    pd.DataFrame({"perm_mean": strat_perm_means["same_drug"]}).to_csv(
        out_dir / f"rung0_permutation_perm_means_same_drug{suffix}.csv", index=False
    )
    pd.DataFrame({"perm_mean": strat_perm_means["diff_drug"]}).to_csv(
        out_dir / f"rung0_permutation_perm_means_diff_drug{suffix}.csv", index=False
    )

    # The design declares this figure so the design effect is SEEN rather than asserted -- the
    # permutation null's spread against the independent-pool spread the bootstrap assumed. It is
    # drawn here because this is the only job that has the permutation draws.
    draws_path = out_dir / "rung0_null_draws.csv"
    reliability_path = out_dir / "rung0_reliability.csv"
    if draws_path.exists() and reliability_path.exists():
        from fmharness import figures as fg

        fig_dir = out_dir / "figures"
        fig_dir.mkdir(parents=True, exist_ok=True)
        rel = pd.read_csv(reliability_path).iloc[0].to_dict()
        fg.fig_permutation_vs_bootstrap(
            pd.DataFrame({"perm_mean": perm_means}),
            pd.read_csv(draws_path),
            {**rel, **summary_row},
            fig_dir / f"10_permutation_vs_bootstrap{suffix}.png",
        )
        print(f"wrote {fig_dir}/10_permutation_vs_bootstrap{suffix}.png")
    else:
        print(
            "skipping the design-effect figure: it is drawn from rung0_null_draws.csv and "
            "rung0_reliability.csv, which this job did not produce -- run "
            "scripts/delta_reproducibility.py into the same --out-dir first"
        )

    print("\n=== permutation-based exact permutation null (rung 0, final verification step) ===")
    for k, v in summary_row.items():
        print(f"  {k:22s} {v}")
    print(
        f"\nany-pair: observed mean r = {summary['observed_mean']:.4f}, permutation null mean = "
        f"{summary['perm_mean_mean']:.4f} +/- {summary['perm_mean_sd']:.4f} (sd), "
        f"p_exact = {summary['p_exact']:.4f}, design effect = {summary['design_effect']:.3f} "
        f"(se_iid_pool = {summary['se_iid_pool']:.5f})."
    )
    print(
        f"same-drug ({strat_summary['n_rows_same_drug']} rows in >=2-row drug groups): "
        f"observed mean r = {strat_summary['observed_mean_same_drug_rows']:.4f}, "
        f"within-drug null mean = {strat_summary['perm_mean_mean_same_drug']:.4f} +/- "
        f"{strat_summary['perm_mean_sd_same_drug']:.4f} (sd), p_exact = "
        f"{strat_summary['p_exact_same_drug']:.4f}, design effect = "
        f"{strat_summary['design_effect_same_drug']:.3f}."
    )
    transfer_note = (
        "(transfers directly to the promoted p_vs_same_drug)"
        if strat_summary["same_drug_rows_equal_n"]
        else "(applies to the multi-line subset only -- state this wherever the number is cited)"
    )
    print(
        f"same-drug design effect measured over {strat_summary['n_rows_same_drug']} of "
        f"{summary['n_pairs']} rows {transfer_note}"
    )
    print(
        f"diff-drug: observed mean r = {summary['observed_mean']:.4f}, cross null mean = "
        f"{strat_summary['perm_mean_mean_diff_drug']:.4f} +/- "
        f"{strat_summary['perm_mean_sd_diff_drug']:.4f} (sd), p_exact = "
        f"{strat_summary['p_exact_diff_drug']:.4f}, design effect = "
        f"{strat_summary['design_effect_diff_drug']:.3f}.\nwrote {summary_path}"
    )


if __name__ == "__main__":
    main()
