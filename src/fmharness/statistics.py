"""Shared significance-testing helpers.

Every p-value in this repo that compares a REPORTED AGGREGATE (a mean or median over many
(line, drug) pairs) against a null must resample the null to that same aggregate, at the same
pair count -- not compare the aggregate to the spread of individual null draws. An aggregate's
sampling distribution is roughly sqrt(n) tighter than a single draw's, so the wrong comparison
inflates p by orders of magnitude and can turn a real signal into an apparent null.

This was first found and fixed for the rung-0 replicate ceiling on the archived lineage
(branch ``rung0-replicate-ceiling-old-lineage``), where it had been independently
reintroduced in several scripts. One shared, tested helper is how it stays fixed
everywhere at once.
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np


def bootstrap_aggregate_pvalue(
    observed_agg: float,
    null_draws: np.ndarray,
    n_obs: int,
    *,
    agg: Callable[[np.ndarray], float] = np.mean,
    n_boot: int = 2000,
    seed: int = 0,
    min_null_draws: int = 10,
) -> tuple[float, float, float]:
    """One-sided p (and a 95% interval) for ``observed_agg`` against a bootstrapped null.

    ``null_draws`` is a pool of INDIVIDUAL mismatched/shuffled-pair statistics (not aggregates).
    This resamples ``n_obs`` of them with replacement, applies ``agg`` (mean or median, matching
    however the observed statistic was aggregated), and repeats ``n_boot`` times to build the
    null's sampling distribution for an aggregate over ``n_obs`` items -- the same kind of
    quantity as ``observed_agg``, so the comparison is like-for-like.

    Returns ``(p, ci_lo, ci_hi)``, each ``nan`` if ``null_draws`` has fewer than
    ``min_null_draws`` entries (too few to bootstrap from).
    """
    null_draws = np.asarray(null_draws, dtype=np.float64)
    null_draws = null_draws[np.isfinite(null_draws)]
    if null_draws.size < min_null_draws or n_obs < 1:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    boot = np.array([agg(rng.choice(null_draws, size=n_obs, replace=True)) for _ in range(n_boot)])
    p = float((1 + np.sum(boot >= observed_agg)) / (1 + boot.size))
    lo, hi = float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))
    return p, lo, hi


def minimum_detectable_aggregate(
    observed_draws: np.ndarray,
    null_draws: np.ndarray,
    n_obs: int,
    *,
    agg: Callable[[np.ndarray], float] = np.mean,
    alpha: float = 0.05,
    power: float = 0.80,
    n_boot: int = 2000,
    seed: int = 0,
    min_null_draws: int = 10,
) -> float:
    """Smallest true aggregate detectable against ``null_draws`` at (``alpha``, ``power``).

    The same bootstrap as ``bootstrap_aggregate_pvalue``, read in the other direction.
    The detection threshold is the (1 - alpha) quantile of the null aggregate's sampling
    distribution at ``n_obs``; the MDE is the aggregate value whose own sampling
    distribution -- the observed draws' spread, recentred -- clears that threshold with
    probability ``power``. Valid for shift-equivariant aggregates (mean, median,
    quantiles), where shifting every draw by c shifts the aggregate by c.

    Reported beside every promoted comparison: a null result with no MDE cannot be told
    apart from an underpowered experiment.
    """
    null_draws = np.asarray(null_draws, dtype=np.float64)
    null_draws = null_draws[np.isfinite(null_draws)]
    observed_draws = np.asarray(observed_draws, dtype=np.float64)
    observed_draws = observed_draws[np.isfinite(observed_draws)]
    if null_draws.size < min_null_draws or observed_draws.size < 2 or n_obs < 1:
        return float("nan")
    rng = np.random.default_rng(seed)
    boot_null = np.array(
        [agg(rng.choice(null_draws, size=n_obs, replace=True)) for _ in range(n_boot)]
    )
    crit = float(np.quantile(boot_null, 1.0 - alpha))
    centred = observed_draws - float(np.mean(observed_draws))
    boot_centred = np.array(
        [agg(rng.choice(centred, size=n_obs, replace=True)) for _ in range(n_boot)]
    )
    # P(agg + theta >= crit) >= power  <=>  theta >= crit - Q_{boot_centred}(1 - power)
    return float(crit - np.quantile(boot_centred, 1.0 - power))


def spearman_brown(r: float) -> float:
    """Lift a split-half reliability to full-data reliability: ``2r / (1 + r)``.

    The one implementation of the lift the ladder states wherever it is used, so a rung
    cannot quietly compute a different correction. Undefined at r = -1; callers guard.
    """
    return 2.0 * r / (1.0 + r)
