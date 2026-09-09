"""Known-answer validation for the statistics rung 0 reports (SPEC rule 4).

Each test plants an answer and requires the real, shipped function to recover it, and
plants nothing and requires null. This file exists because, on the archived lineage, a
null test compared an aggregate against a distribution of single draws and reported that
a reproducible ceiling had failed; running against a known answer catches that class of
defect immediately. The recurring error is aggregate-vs-per-item, so units are pinned
explicitly throughout.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from fmharness.statistics import (
    bootstrap_aggregate_pvalue,
    minimum_detectable_aggregate,
    spearman_brown,
)

pytestmark = pytest.mark.known_answer


def test_bootstrap_aggregate_pvalue_returns_null_when_there_is_no_signal() -> None:
    rng = np.random.default_rng(1)
    observed = rng.normal(0.14, 0.13, 1600)
    null = rng.normal(0.14, 0.13, 500)
    p, _, _ = bootstrap_aggregate_pvalue(float(np.mean(observed)), null, observed.size)
    assert p > 0.05, f"no-signal case must not be significant, got p={p}"


def test_bootstrap_aggregate_pvalue_recovers_a_planted_difference() -> None:
    rng = np.random.default_rng(2)
    observed = rng.normal(0.30, 0.13, 1600)
    null = rng.normal(0.14, 0.13, 500)
    p, _, _ = bootstrap_aggregate_pvalue(float(np.mean(observed)), null, observed.size)
    assert p < 0.01, f"planted difference must be recovered, got p={p}"


def test_the_wrong_form_is_the_one_that_fails_to_recover_it() -> None:
    # Pins the aggregate-vs-per-item defect: with the SAME planted data, comparing the
    # aggregate to individual null draws looks non-significant where the real function
    # does not.
    rng = np.random.default_rng(2)
    observed = rng.normal(0.30, 0.13, 1600)
    null = rng.normal(0.14, 0.13, 500)
    wrong = float(np.mean(null >= np.mean(observed)))
    correct, _, _ = bootstrap_aggregate_pvalue(float(np.mean(observed)), null, observed.size)
    assert wrong > 0.05, "the defective form should look non-significant on planted signal"
    assert correct < 0.01
    assert wrong > correct * 10, "the defect inflates p by orders of magnitude"


def test_spearman_brown_lifts_half_data_reliability_as_documented() -> None:
    assert spearman_brown(0.0) == 0.0
    assert abs(spearman_brown(1.0) - 1.0) < 1e-12
    assert spearman_brown(0.3) > 0.3, "correcting half-data reliability must raise it"
    assert spearman_brown(0.5) > spearman_brown(0.3)


def test_minimum_detectable_aggregate_matches_the_normal_closed_form() -> None:
    # Normal null (mu0, sigma) and normal observed spread (sigma), agg = mean over n:
    #   MDE = mu0 + (z_{1-alpha} + z_{power}) * sigma / sqrt(n)
    rng = np.random.default_rng(0)
    mu0, sigma, n = 0.03, 0.10, 400
    null = rng.normal(mu0, sigma, 2000)
    observed = rng.normal(0.20, sigma, 2000)
    mde = minimum_detectable_aggregate(observed, null, n)
    z95, z80 = 1.6449, 0.8416
    expected = mu0 + (z95 + z80) * sigma / math.sqrt(n)
    assert abs(mde - expected) < 0.005, f"mde={mde}, closed form={expected}"


def test_an_aggregate_at_the_mde_clears_the_null() -> None:
    # Self-consistency with the p-value: the MDE sits above the null's critical value,
    # so a result at the MDE must come out significant at the same alpha.
    rng = np.random.default_rng(3)
    null = rng.normal(0.03, 0.10, 2000)
    observed = rng.normal(0.10, 0.10, 2000)
    n = 400
    mde = minimum_detectable_aggregate(observed, null, n)
    p, _, _ = bootstrap_aggregate_pvalue(mde, null, n)
    assert p < 0.05, f"an aggregate at the MDE must be significant, got p={p}"


def test_minimum_detectable_aggregate_returns_nan_on_too_few_null_draws() -> None:
    assert math.isnan(minimum_detectable_aggregate(np.ones(50), np.ones(3), 10))
