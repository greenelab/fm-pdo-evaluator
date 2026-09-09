"""Known-answer tests for the permutation-based exact permutation null (SPEC rule 4).

This is rung 0's final verification step (docs/tasks/rung0-replicate-ceiling/design.md,
decision history), added to carry -- by construction -- the dependence the bootstrapped
p-values in `scripts/delta_reproducibility.py` had to assume away: the stratified null draws
in that script reuse the same half-profiles across many mismatched-pair comparisons, so they
are not an exchangeable i.i.d. pool even though the bootstrap treats them as one
(docs/tasks/rung0-replicate-ceiling/verification.md, "Write-up caveat"). Permutations of the
pairing preserve every row exactly once while eliminating the possibility of a matched pair,
so the resulting null distribution of the mean carries the real dependence structure rather
than an assumed one. Every test here runs the REAL functions from
`scripts/permutation_null.py` (which itself reuses, not reimplements, the measurement core in
`scripts/delta_reproducibility.py`) on synthetic replicate pools with planted, known answers.
"""

# pandas ships no PEP-561 type stubs in this environment; under strict mode that turns every
# pandas call site into a cascade of reportUnknown* noise about *pandas'* types, not ours. Same
# suppression, same rationale as the rest of this project's pyright strict config where it
# touches scientific-Python packages -- the rules that check our own code stay on.
# pyright: reportMissingTypeStubs=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownParameterType=false, reportUnknownLambdaType=false

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

# `write_fixture_pool` lives in tests/test_rung0_controls.py; the leading underscore marks it
# as private to that module, not a boundary this sibling test file needs to respect -- both
# test the same rung 0 task, and re-declaring the fixture builder here would be the duplication
# PROCESS §3 tells tests not to reintroduce. The suppression covers importing a sibling test
# module's private helper for exactly that reason.
from tests.test_rung0_controls import write_fixture_pool  # pyright: ignore[reportPrivateUsage]

pytestmark = pytest.mark.known_answer

_SPEC = importlib.util.spec_from_file_location(
    "permutation_null",
    Path(__file__).resolve().parents[1] / "scripts" / "permutation_null.py",
)
assert _SPEC is not None and _SPEC.loader is not None
dn = importlib.util.module_from_spec(_SPEC)
sys.modules["permutation_null"] = dn
_SPEC.loader.exec_module(dn)


def test_sample_permutation_has_no_fixed_points() -> None:
    """A rejection-sampled permutation has zero fixed points, at several sizes and seeds."""
    for n in (2, 3, 5, 10, 50):
        for seed in range(5):
            rng = np.random.default_rng(seed * 1000 + n)
            sigma = dn.sample_permutation(rng, n)
            assert sigma.shape == (n,)
            assert sorted(sigma.tolist()) == list(range(n)), f"not a permutation: {sigma}"
            assert not np.any(sigma == np.arange(n)), f"n={n} seed={seed}: fixed point in {sigma}"


def test_planted_signal_pool_clears_the_permutation_null(tmp_path: Path) -> None:
    # signal_sd = noise_sd = 1, 8 plates -> expected split-half r ~ 0.8 (same planted shape as
    # test_score_positive_planted_reliability_is_recovered in tests/test_rung0_controls.py).
    # Mismatching the pairing via a permutation destroys that pair-specific signal entirely
    # (drug_sd = 0.0, the default -- no drug-shared component survives a mismatch either), so
    # the permutation null should sit near zero while the observed mean sits near 0.8.
    path = write_fixture_pool(tmp_path, signal_sd=1.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    summary, _ = dn.permutation_null(piv0, piv1, r, min_genes=50, n_perm=99, seed=0)
    assert summary["p_exact"] < 0.05, (
        f"planted signal must clear the null, got p={summary['p_exact']}"
    )
    assert summary["observed_mean"] > summary["perm_mean_mean"] + 0.1, (
        f"observed {summary['observed_mean']} not far above perm null mean "
        f"{summary['perm_mean_mean']}"
    )


def test_zero_signal_pool_is_not_significant(tmp_path: Path) -> None:
    # No planted signal at all: matched and mismatched pairs are drawn from the same
    # generative process, so the permutation null and the observed mean should sit at the
    # same (near-zero) place. Deterministic under the pinned seed; if this ever flips the
    # failure message below records the actual p rather than a bare assertion.
    path = write_fixture_pool(tmp_path, signal_sd=0.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    summary, _ = dn.permutation_null(piv0, piv1, r, min_genes=50, n_perm=200, seed=0)
    assert summary["p_exact"] > 0.01, (
        f"zero-signal pool must not clear the permutation null, got p={summary['p_exact']} "
        f"(observed_mean={summary['observed_mean']}, perm_mean_mean={summary['perm_mean_mean']})"
    )


def test_permutation_null_rejects_too_few_pairs_or_too_few_permutations() -> None:
    """A single finite pair has no permutation (rejection sampling can never satisfy "no
    fixed points" for n=1), and n_perm=1 gives a null distribution with no spread -- both
    are degenerate inputs the function should refuse before entering the permutation loop,
    rather than failing obscurely inside `sample_permutation` or downstream variance math."""
    import pandas as pd

    one_pair = pd.DataFrame([[1.0, 2.0, 3.0]], index=["a"])
    with pytest.raises(ValueError, match="at least 2"):
        dn.permutation_null(one_pair, one_pair, np.array([0.9]), min_genes=1, n_perm=99, seed=0)

    two_pairs = pd.DataFrame([[1.0, 2.0], [3.0, 4.0]], index=["a", "b"])
    with pytest.raises(ValueError, match="n_perm"):
        dn.permutation_null(
            two_pairs, two_pairs, np.array([0.9, 0.8]), min_genes=1, n_perm=1, seed=0
        )


def test_design_effect_sits_in_a_sane_band_on_near_independent_rows(tmp_path: Path) -> None:
    """The design effect the exchangeable-pool bootstrap ignores, measured on data with no
    shared drug or line component (drug_sd = 0.0, the default -- rows are as close to
    independent as this fixture can make them).

    This is a known answer for "no pathological inflation," not a precise theoretical
    prediction: the mean-r statistic still shares gene-panel structure across rows through the
    finite gene count, and only 12 (line, drug) pairs are available (4 lines x 3 drugs, the
    fixture default), so the permutation-null variance is itself estimated noisily at
    n_perm=200. [0.2, 5.0] is a loose band that would catch an order-of-magnitude inflation
    factor -- the failure mode the write-up caveat in verification.md was worried about --
    while tolerating ordinary Monte Carlo noise at this sample size.
    """
    path = write_fixture_pool(tmp_path, signal_sd=0.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    summary, _ = dn.permutation_null(piv0, piv1, r, min_genes=50, n_perm=200, seed=0)
    assert 0.2 <= summary["design_effect"] <= 5.0, (
        f"design effect {summary['design_effect']} outside the sane band [0.2, 5.0]"
    )


# --- Stratum-preserving permutation nulls -----------------------------------------------
#
# The any-pair permutation above validates the POOLED aggregate, but the promoted p-values
# (`delta_reproducibility.summarize`'s `p_vs_null` and `p_vs_same_drug`) are per-stratum:
# `p_vs_same_drug` is read against a same-drug mismatched-pair pool that clusters over ~32
# drugs, and `p_vs_null` against a diff-drug pool. An any-pair permutation mixes same- and
# diff-drug mismatches freely and carries neither stratum's dependence specifically. These
# tests cover the two stratum-constrained permutation samplers and the stratified null they
# build, on the same synthetic fixtures as the rest of this file and `test_rung0_controls.py`.


def test_sample_within_drug_permutation_excludes_singletons_and_permutes_the_rest() -> None:
    """Uneven drug-group sizes, including a singleton drug with no permutation: every row in a
    >=2-row drug group maps to a DIFFERENT row of the SAME drug; the singleton row has no valid
    target and must map to itself, i.e. sit outside any >=2-row aggregate."""
    drugs = np.array(["D0", "D0", "D0", "D1", "D1", "D2"])  # D2 is a singleton
    multi = np.array([True, True, True, True, True, False])
    rng = np.random.default_rng(0)
    sigma = dn.sample_within_drug_permutation(rng, drugs)
    assert sigma.shape == (6,)
    assert sorted(sigma.tolist()) == list(range(6)), f"not a permutation: {sigma}"
    for i in range(6):
        if multi[i]:
            assert sigma[i] != i, f"row {i} (drug {drugs[i]}) mapped to itself"
            assert drugs[sigma[i]] == drugs[i], (
                f"row {i} (drug {drugs[i]}) mapped to drug {drugs[sigma[i]]}"
            )
        else:
            assert sigma[i] == i, f"singleton row {i} (drug {drugs[i]}) must map to itself"


def test_sample_cross_permutation_satisfies_all_three_constraints() -> None:
    """`sample_cross_permutation` on a fixture shaped like a real (line, drug) pivot index (4
    lines x 3 drugs, one row per combination -- the diff-drug stratum's constraints are easily
    satisfiable here): every row's partner differs in row index, drug, AND line."""
    lines = np.array([f"L{i // 3}" for i in range(12)])
    drugs = np.array([f"D{i % 3}" for i in range(12)])
    rng = np.random.default_rng(0)
    sigma = dn.sample_cross_permutation(rng, drugs, lines)
    assert sigma.shape == (12,)
    assert sorted(sigma.tolist()) == list(range(12)), f"not a permutation: {sigma}"
    idx = np.arange(12)
    assert not np.any(sigma == idx), "some row mapped to itself"
    assert not np.any(drugs[sigma] == drugs), "some row mapped within its own drug"
    assert not np.any(lines[sigma] == lines), "some row mapped within its own line"


def test_stratified_null_drug_effects_survive_within_drug_die_across_drugs(
    tmp_path: Path,
) -> None:
    """Ordering known answer, the stratified analog of
    `test_null_positive_planted_components_recover_the_stratum_ordering` in
    `test_rung0_controls.py`: with a planted drug-shared component (drug_sd > 0) plus
    line-specific signal, a within-drug permutation's mismatched pairs still share the drug
    component (same drug, different line), while a cross permutation's mismatched pairs share
    nothing (different drug AND different line). The within-drug permutation-null MEAN must
    therefore sit above the cross permutation-null MEAN."""
    path = write_fixture_pool(
        tmp_path, n_lines=6, n_drugs=4, signal_sd=0.7, drug_sd=0.7, noise_sd=0.7
    )
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2", "D3"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    strat, _ = dn.stratified_permutation_null(piv0, piv1, r, min_genes=50, n_perm=150, seed=0)
    same_mean, diff_mean = strat["perm_mean_mean_same_drug"], strat["perm_mean_mean_diff_drug"]
    assert same_mean > diff_mean + 0.1, (
        f"within-drug null mean {same_mean:.3f} !> cross null mean {diff_mean:.3f} by 0.1"
    )


def test_stratified_zero_signal_pool_is_not_significant(tmp_path: Path) -> None:
    """Zero-signal negative control for both stratified nulls, and the same design-effect
    sanity band as `test_design_effect_sits_in_a_sane_band_on_near_independent_rows`, applied
    to each stratum."""
    path = write_fixture_pool(tmp_path, signal_sd=0.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    strat, _ = dn.stratified_permutation_null(piv0, piv1, r, min_genes=50, n_perm=200, seed=0)
    assert strat["p_exact_same_drug"] > 0.01, (
        f"zero-signal same-drug stratum must not clear its null, got p={strat['p_exact_same_drug']}"
    )
    assert strat["p_exact_diff_drug"] > 0.01, (
        f"zero-signal diff-drug stratum must not clear its null, got p={strat['p_exact_diff_drug']}"
    )
    assert 0.2 <= strat["design_effect_same_drug"] <= 5.0, (
        f"same-drug design effect {strat['design_effect_same_drug']} outside [0.2, 5.0]"
    )
    assert 0.2 <= strat["design_effect_diff_drug"] <= 5.0, (
        f"diff-drug design effect {strat['design_effect_diff_drug']} outside [0.2, 5.0]"
    )


def test_stratified_planted_signal_clears_both_nulls(tmp_path: Path) -> None:
    """Planted-signal positive control for both stratified nulls, at n_perm=99 -- the same
    fixture as `test_planted_signal_pool_clears_the_permutation_null`."""
    path = write_fixture_pool(tmp_path, signal_sd=1.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    strat, _ = dn.stratified_permutation_null(piv0, piv1, r, min_genes=50, n_perm=99, seed=0)
    assert strat["p_exact_same_drug"] < 0.05, (
        f"planted signal must clear the same-drug null, got p={strat['p_exact_same_drug']}"
    )
    assert strat["p_exact_diff_drug"] < 0.05, (
        f"planted signal must clear the diff-drug null, got p={strat['p_exact_diff_drug']}"
    )


def test_stratified_summary_reports_transfer_scope_and_diff_drug_observed_mean(
    tmp_path: Path,
) -> None:
    """`same_drug_rows_equal_n` records whether the same-drug design effect (measured over the
    >=2-row-drug-group subset) transfers directly to the promoted `p_vs_same_drug` (measured
    over every row): the fixture's default composition gives every drug all 4 lines, so no row
    is excluded and the two counts coincide. `observed_mean_diff_drug_rows` is recorded per
    stratum purely so the summary is self-describing -- by construction (the diff-drug stratum
    excludes no rows) it equals the global observed mean exactly."""
    path = write_fixture_pool(tmp_path, signal_sd=1.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    strat, _ = dn.stratified_permutation_null(piv0, piv1, r, min_genes=50, n_perm=99, seed=0)
    assert strat["n_rows_same_drug"] == 12, "fixture default: 4 lines x 3 drugs, no singletons"
    assert strat["same_drug_rows_equal_n"] is True
    r_fin = r[np.isfinite(r)]
    assert strat["observed_mean_diff_drug_rows"] == round(float(np.mean(r_fin)), 4)


def test_stratified_design_effect_is_nan_when_the_matching_pool_is_too_small(
    tmp_path: Path,
) -> None:
    """A stratum's `stratified_null_draws` pool below `MIN_NULL_DRAWS_FOR_DESIGN_EFFECT` (10,
    mirroring `bootstrap_aggregate_pvalue`'s `min_null_draws` spirit) must not feed
    `np.var(ddof=1)` silently -- that stratum's design_effect is nan instead of a numerically
    fragile ratio built on too few draws. n_perm=5 caps both the permutation draws AND
    `stratified_null_draws`' pool size at <=5 for every stratum, regardless of how many
    candidate mismatched pairs the fixture actually has."""
    path = write_fixture_pool(tmp_path, signal_sd=0.0, noise_sd=1.0)
    de, _ = dn.dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dn.dr.score_split_half(de, set(de["gene_name"].unique()))
    strat, _ = dn.stratified_permutation_null(piv0, piv1, r, min_genes=50, n_perm=5, seed=0)
    assert np.isnan(strat["design_effect_same_drug"]), (
        f"expected nan for a <10-draw pool, got {strat['design_effect_same_drug']}"
    )
    assert np.isnan(strat["design_effect_diff_drug"]), (
        f"expected nan for a <10-draw pool, got {strat['design_effect_diff_drug']}"
    )
