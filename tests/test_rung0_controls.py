"""Positive and negative controls for rung 0's build, score, and null steps (SPEC rule 4).

Every test runs the REAL functions from scripts/delta_reproducibility.py on synthetic
replicate pools with planted, known answers -- not reimplementations of their logic.
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
import pandas as pd
import pytest

pytestmark = pytest.mark.known_answer

_SPEC = importlib.util.spec_from_file_location(
    "delta_reproducibility",
    Path(__file__).resolve().parents[1] / "scripts" / "delta_reproducibility.py",
)
assert _SPEC is not None and _SPEC.loader is not None
dr = importlib.util.module_from_spec(_SPEC)
sys.modules["delta_reproducibility"] = dr
_SPEC.loader.exec_module(dr)


def _write_fixture_pool(
    tmp: Path,
    n_lines: int = 4,
    n_drugs: int = 3,
    n_genes: int = 300,
    plates: tuple[str, ...] = ("P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"),
    signal_sd: float = 1.0,
    noise_sd: float = 1.0,
    drug_sd: float = 0.0,
    seed: int = 0,
    n_responders: int | None = None,
    doses: tuple[float, ...] = (0.05,),
    plate_offset_sd: float = 0.0,
    se: float | None = None,
) -> Path:
    """A synthetic replicate pool in the DE table's own shape, one parquet file.

    Per (line, drug, gene): a fixed pair-specific signal (sd ``signal_sd``), an optional
    drug-shared component (sd ``drug_sd``), plus independent per-plate noise (sd
    ``noise_sd``). Expected split-half r over genes, as plate count grows:
    (signal_sd^2 + drug_sd^2) / (signal_sd^2 + drug_sd^2 + noise_sd^2 / plates_per_half).

    The DESeq2 columns the assay-reliability task reads are planted too:

    ``n_responders``  the first N genes are planted as differentially expressed -- ``padj``
                      drawn below 0.01 -- and the rest above 0.2. ``None`` (the default) plants
                      ``padj`` uniform on (0, 1) for every gene, which is the signal-free case
                      selection must admit at no more than the nominal rate.
    ``doses``         one row per (plate, dose); the reliabilities pool over dose, the noise
                      decomposition holds it fixed.
    ``plate_offset_sd``  a per-plate offset shared across genes and doses -- the plate effect
                      ``lfcSE`` cannot see and the decomposition exists to recover. It is added
                      ON TOP of ``noise_sd``, so the planted plate variance is
                      ``plate_offset_sd ** 2`` and the planted within-plate variance is
                      ``noise_sd ** 2``.
    ``se``            the value written to ``lfcSE``. Defaults to ``noise_sd``, which is the
                      truth for this generator: each row's deviation from its condition mean is
                      drawn at sd ``noise_sd``, so a standard error of ``noise_sd`` is what a
                      correctly calibrated DESeq2 would report.
    """
    rng = np.random.default_rng(seed)
    lines = [f"L{i}" for i in range(n_lines)]
    drugs = [f"D{j}" for j in range(n_drugs)]
    genes = [f"G{k}" for k in range(n_genes)]
    drug_eff = {d: rng.normal(0.0, drug_sd, n_genes) for d in drugs}
    plate_off = {p: rng.normal(0.0, plate_offset_sd, 1)[0] for p in plates}
    lfc_se = noise_sd if se is None else se
    rows = []
    for li in lines:
        for d in drugs:
            signal = rng.normal(0.0, signal_sd, n_genes) + drug_eff[d]
            for p in plates:
                for dose in doses:
                    lfc = signal + rng.normal(0.0, noise_sd, n_genes) + plate_off[p]
                    if n_responders is None:
                        padj = rng.uniform(0.0, 1.0, n_genes)
                    else:
                        padj = np.where(
                            np.arange(n_genes) < n_responders,
                            rng.uniform(0.0, 0.01, n_genes),
                            rng.uniform(0.2, 1.0, n_genes),
                        )
                    rows.append(
                        pd.DataFrame(
                            {
                                "Cell_ID_DepMap": li,
                                "drug": d,
                                "gene_name": genes,
                                "log2FoldChange": lfc,
                                "lfcSE": lfc_se,
                                "padj": padj,
                                "baseMean": 100.0,
                                "concentration": dose,
                                "plate": p,
                            }
                        )
                    )
    pool_dir = tmp / "pseudobulk_differential_expression"
    pool_dir.mkdir(parents=True)
    out = pool_dir / "train-00000-of-00001.parquet"
    pd.concat(rows, ignore_index=True).to_parquet(out, index=False)
    return out


def test_build_positive_planted_pool_comes_out_with_the_planted_shape(tmp_path: Path) -> None:
    path = _write_fixture_pool(tmp_path)
    de, chosen = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    assert chosen == "plate"
    pairs = de.groupby(["patient", "drug"]).ngroups
    assert pairs == 12, f"planted 4 lines x 3 drugs, built {pairs} pairs"
    assert de["gene_name"].nunique() == 300


def test_build_negative_no_replication_yields_no_scoreable_pairs(tmp_path: Path) -> None:
    path = _write_fixture_pool(tmp_path, plates=("P1",))  # one plate: one half stays empty
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    assert de.dropna(subset=["lfc0", "lfc1"]).empty


def test_build_edge_unmatched_target_drugs_return_an_empty_frame(tmp_path: Path) -> None:
    """ "Empty input" (a pool with zero rows for the requested target drugs) and "unmatched
    target identifiers" (an explicit query for drug names absent from the file) are the same
    mechanism in the real loader: `build_split_half_frame`'s WHERE clause filters on
    `target_names`, so both reduce to requesting drugs the pool does not contain. Read off the
    real function rather than assumed: it returns an empty DataFrame, not an error."""
    path = _write_fixture_pool(tmp_path)  # the fixture pool only ever has drugs D0, D1, D2
    de, chosen = dr.build_split_half_frame(
        [str(path)], ["ZZZ_NOT_A_REAL_DRUG"], None, tmp_path / "duck", memory_limit="2GB"
    )
    assert de.empty, "requesting target drugs absent from the pool must return an empty frame"
    assert chosen == "plate"  # the replicate column is still resolved; only the rows are empty


def test_build_edge_all_nan_pair_drops_out_after_the_dropna_path(tmp_path: Path) -> None:
    """A (line, drug) pair whose log2FoldChange is NaN across every one of its rows aggregates
    to NaN in both halves (DuckDB's avg() over an all-NaN DOUBLE column returns NaN, verified
    directly, not NULL -- so this is a distinct code path from a pair simply absent from the
    data) and is then removed by the `dropna(subset=["lfc0", "lfc1"])` every caller applies
    before scoring."""
    path = _write_fixture_pool(tmp_path)
    df = pd.read_parquet(path)
    nan_mask = (df["Cell_ID_DepMap"] == "L0") & (df["drug"] == "D0")
    assert nan_mask.any(), "fixture must actually contain the (L0, D0) pair to NaN out"
    df.loc[nan_mask, "log2FoldChange"] = np.nan
    df.to_parquet(path, index=False)

    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    remaining_pairs = {tuple(row) for row in de[["patient", "drug"]].drop_duplicates().to_numpy()}
    assert ("L0", "D0") not in remaining_pairs, "an all-NaN pair must not survive the dropna path"
    assert ("L1", "D0") in remaining_pairs, "other pairs must be unaffected"


@pytest.mark.step_build
def test_build_admits_every_drug_when_no_drug_list_is_given(tmp_path: Path) -> None:
    """Positive control for build, this task's inclusion rule: every splittable drug.

    The superseded rung passed a 32-compound list this task does not have and cannot rebuild;
    ``target_names=None`` is what "no drug file" means at the SQL level, and it must admit the
    whole pool rather than silently matching nothing.
    """
    path = _write_fixture_pool(tmp_path)
    de, chosen = dr.build_split_half_frame(
        [str(path)], None, None, tmp_path / "duck", memory_limit="2GB"
    )
    assert chosen == "plate"
    assert set(de["drug"].unique()) == {"D0", "D1", "D2"}
    assert de.dropna(subset=["lfc0", "lfc1"]).groupby(["patient", "drug"]).ngroups == 12


@pytest.mark.step_build
def test_build_carries_the_first_groups_minimum_padj(tmp_path: Path) -> None:
    """The selection rule is "significant in AT LEAST ONE of the first group's rows", so the
    aggregate that decides it is the MINIMUM adjusted p-value over those rows, not their mean.
    A mean would let one strongly significant row be averaged away by its neighbours."""
    path = _write_fixture_pool(tmp_path)
    df = pd.read_parquet(path)
    # Which plates land in the first group is DuckDB's hash, so read it off the engine rather
    # than assuming: plant into whichever plates the build itself calls group 0.
    import duckdb

    half0 = {
        str(p)
        for p in duckdb.connect()
        .execute(
            "SELECT DISTINCT plate FROM read_parquet(?) WHERE hash(plate) % 2 = 0", [str(path)]
        )
        .df()["plate"]
    }
    assert half0 and half0 != set(df["plate"].unique()), "fixture must split across both groups"
    first = sorted(half0)[0]
    key = (df["Cell_ID_DepMap"] == "L0") & (df["drug"] == "D0") & (df["gene_name"] == "G0")
    df.loc[key, "padj"] = 0.9
    df.loc[key & (df["plate"] == first), "padj"] = 0.01
    other = (df["Cell_ID_DepMap"] == "L0") & (df["drug"] == "D0") & (df["gene_name"] == "G1")
    df.loc[other, "padj"] = 0.9
    df.to_parquet(path, index=False)

    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp_path / "duck", "2GB")
    got = de.set_index(["patient", "drug", "gene_name"])["padj0"]
    assert got.loc[("L0", "D0", "G0")] == pytest.approx(0.01), "one significant row must carry"
    assert got.loc[("L0", "D0", "G1")] == pytest.approx(0.9), "no significant row, no selection"


@pytest.mark.step_build
def test_build_leaves_untestable_genes_null(tmp_path: Path) -> None:
    """A gene DESeq2 could not test carries baseMean 0 and null in every statistic column.
    Such genes must fall out by the finiteness rule the scorer already applies -- not by a
    filter of ours, which would be a second, undeclared inclusion rule."""
    path = _write_fixture_pool(tmp_path)
    df = pd.read_parquet(path)
    dead = df["gene_name"] == "G7"
    df.loc[dead, ["log2FoldChange", "lfcSE", "padj"]] = np.nan
    df.loc[dead, "baseMean"] = 0.0
    df.to_parquet(path, index=False)

    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp_path / "duck", "2GB")
    # The build now applies the both-halves-present rule in SQL rather than leaving it to the
    # caller's dropna, because at full extent that filter removes most of 1.42 billion rows and
    # doing it in pandas means materialising them first. The observable property is unchanged
    # and is what this asserts: an untestable gene is scored by nothing, and it is absent for
    # the finiteness reason rather than by any rule naming baseMean or padj.
    assert de[de["gene_name"] == "G7"].empty, "an untestable gene must reach no statistic"
    assert de["gene_name"].nunique() == 299, "every other gene survives"
    assert de.dropna(subset=["lfc0", "lfc1"])["gene_name"].nunique() == 299, (
        "the caller's dropna must now be a no-op: the engine already applied it"
    )


def _first_group_plates(path: Path) -> set[str]:
    """Which plates the build's own ``hash(plate) % 2`` puts in the first group.

    Read off DuckDB rather than reimplemented: Python's ``hash`` is salted per process and is
    not DuckDB's hash, so a test that computed the split itself would be testing a different
    split from the one the measurement uses.
    """
    import duckdb

    rows = (
        duckdb.connect()
        .execute(
            "SELECT DISTINCT plate FROM read_parquet(?) WHERE hash(plate) % 2 = 0", [str(path)]
        )
        .df()["plate"]
    )
    return {str(p) for p in rows}


def _scored(path: Path, tmp: Path, min_genes: int = 50):
    """Build, drop unpaired rows, and score both gene sets on the real code path.

    Returns ``(r_all, r_resp, piv0, piv1, mask)`` -- what every selection control needs, from
    the shipped functions rather than a reimplementation of them.
    """
    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp / "duck", "2GB")
    de = de.dropna(subset=["lfc0", "lfc1"])
    panel = set(de["gene_name"].unique())
    r_all, piv0, piv1 = dr.score_split_half(de, panel, min_genes=min_genes)
    mask = dr.responder_mask(dr.padj_pivot(de, panel).reindex(columns=piv0.columns).loc[piv0.index])
    r_resp, _, _ = dr.score_split_half(de, panel, min_genes=min_genes, select=mask)
    return r_all, r_resp, piv0, piv1, mask


@pytest.mark.step_select
def test_selection_recovers_the_planted_responder_set(tmp_path: Path) -> None:
    """Positive control for select: responders planted in a known gene subset, with padj
    planted to match, come back as exactly that subset -- read from the first group alone."""
    path = _write_fixture_pool(tmp_path, n_genes=300, n_responders=80)
    _, _, piv0, _, mask = _scored(path, tmp_path)
    responders = {f"G{k}" for k in range(80)}
    for row in range(mask.shape[0]):
        got = {str(g) for g, keep in zip(piv0.columns, mask[row], strict=True) if keep}
        assert got == responders, f"row {row} selected {len(got)} genes, planted 80"


@pytest.mark.step_select
def test_responder_reliability_exceeds_all_gene_reliability_on_a_planted_pool(
    tmp_path: Path,
) -> None:
    """The responders carry the response; the rest carry noise around zero. Scoring only the
    responders must therefore read higher than scoring everything, on a pool where that is
    true by construction."""
    # Responders get a real per-condition signal; non-responders get none (signal_sd applies to
    # all genes, so plant the difference through padj AND through the signal by zeroing it).
    path = _write_fixture_pool(tmp_path, n_genes=300, n_responders=80, seed=3)
    df = pd.read_parquet(path)
    non = df["gene_name"].isin([f"G{k}" for k in range(80, 300)])
    rng = np.random.default_rng(11)
    n_non = int(non.to_numpy().sum())
    df.loc[non, "log2FoldChange"] = rng.normal(0.0, 1.0, n_non)
    df.to_parquet(path, index=False)

    r_all, r_resp, _, _, _ = _scored(path, tmp_path)
    mean_all = float(np.nanmean(r_all))
    mean_resp = float(np.nanmean(r_resp))
    assert mean_resp > mean_all + 0.05, (
        f"responder r {mean_resp:.3f} must clear all-gene r {mean_all:.3f} on a pool where only "
        "the responders carry reproducible signal"
    )


@pytest.mark.step_select
def test_selection_admits_no_more_than_the_nominal_rate_on_signal_free_data(
    tmp_path: Path,
) -> None:
    """Negative control for select, and the known answer for the rule's own multiplicity.

    The rule is "significant in AT LEAST ONE of the first group's rows", so under the null a
    gene is admitted when the smallest of its k adjusted p-values falls below alpha:

        P(selected | no signal) = 1 - (1 - alpha) ** k

    not alpha. With four first-group plates that is 0.185, not 0.05. This is a property of the
    declared rule rather than a defect -- and it does not bias the responder reliability, since
    the mismatched-pair nulls apply the same rule to the same first group. What it does mean is
    that the responder set under the null is a fifth of the genes rather than a twentieth, so
    the responder statistic is diluted toward the all-gene one rather than inflated away from
    it. Pinned here so the day someone changes the aggregate, the rate moves and this fails.
    """
    path = _write_fixture_pool(
        tmp_path, n_genes=2000, signal_sd=0.0, noise_sd=1.0, n_responders=None, seed=5
    )
    n_first_group_rows = len(_first_group_plates(path))  # one dose in this fixture
    expected = 1.0 - (1.0 - 0.05) ** n_first_group_rows
    _, r_resp, _, _, mask = _scored(path, tmp_path, min_genes=10)
    rate = float(mask.mean())
    # 2000 genes x 12 conditions; four binomial standard deviations at this rate is under 0.02.
    assert rate == pytest.approx(expected, abs=0.02), (
        f"selected {rate:.3f} of genes; the at-least-one-row rule over {n_first_group_rows} "
        f"first-group rows predicts {expected:.3f}"
    )
    assert abs(float(np.nanmean(r_resp))) < 0.15, "signal-free responder r must sit near zero"


@pytest.mark.step_select
def test_pooled_selection_inflates_a_signal_free_correlation(tmp_path: Path) -> None:
    """The leakage check the one-sided rule exists to prevent, on data with no signal at all.

    The inflating variant is selection on the POOLED data -- the natural mistake of calling
    differential expression on every plate at once and then correlating the halves over those
    genes. Write the two halves as ``a`` and ``b``; their sum and difference are independent.
    Selecting on a large ``|a + b|`` inflates the variance of the sum while leaving the
    difference alone, and since ``cov(a, b) = (var(a + b) - var(a - b)) / 4``, the covariance
    within the selected genes is positive even when nothing generated it. That is the whole of
    the winner's curse here, and it is why the design forbids the pooled and two-sided variants
    rather than treating them as a defensible alternative.

    Selecting on each half's magnitude SEPARATELY is not the same thing and does not inflate:
    truncating ``|a|`` and ``|b|`` independently leaves their signs independent. The distinction
    matters, so the test states which variant it is demonstrating.
    """
    path = _write_fixture_pool(
        tmp_path, n_genes=2000, signal_sd=0.0, noise_sd=1.0, n_responders=None, seed=7
    )
    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp_path / "duck", "2GB")
    de = de.dropna(subset=["lfc0", "lfc1"])
    panel = set(de["gene_name"].unique())
    _, piv0, piv1 = dr.score_split_half(de, panel, min_genes=10)
    padj = dr.padj_pivot(de, panel).reindex(columns=piv0.columns).loc[piv0.index]
    one_sided = dr.responder_mask(padj)

    # The forbidden variant, written HERE and never in the shipped code, at the same gene count
    # as the one-sided rule admits, so the comparison is of selection rules and not of set size.
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    k = max(int(one_sided.sum(axis=1).mean()), 10)
    pooled = np.zeros_like(one_sided)
    strength = np.abs(a + b)
    for row in range(strength.shape[0]):
        pooled[row, np.argsort(-strength[row])[:k]] = True

    r_one = float(np.nanmean(dr.masked_rowwise_pearson(a, b, 10, select=one_sided)))
    r_pooled = float(np.nanmean(dr.masked_rowwise_pearson(a, b, 10, select=pooled)))
    assert r_pooled > r_one + 0.2, (
        f"pooled selection read {r_pooled:.3f} against one-sided {r_one:.3f} on signal-free "
        "data; the gap is the winner's curse the one-sided rule avoids"
    )
    assert abs(r_one) < 0.1, f"the shipped one-sided rule must stay near zero, read {r_one:.3f}"


@pytest.mark.step_select
def test_select_none_reproduces_the_unselected_scorer_exactly(tmp_path: Path) -> None:
    """The selection keyword must be inert when absent -- otherwise the all-gene reliability
    silently changes meaning the day the responder statistic arrives."""
    path = _write_fixture_pool(tmp_path)
    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp_path / "duck", "2GB")
    de = de.dropna(subset=["lfc0", "lfc1"])
    panel = set(de["gene_name"].unique())
    r_a, piv0, piv1 = dr.score_split_half(de, panel)
    a, b = piv0.to_numpy(dtype=float), piv1.to_numpy(dtype=float)
    r_b = dr.masked_rowwise_pearson(a, b, 50, select=np.ones(a.shape, dtype=bool))
    np.testing.assert_allclose(r_a, r_b, equal_nan=True)


def _decomposed(path: Path, tmp: Path) -> pd.DataFrame:
    """The noise decomposition off the real code path."""
    noise = dr.build_noise_frame([str(path)], None, None, tmp / "duck", "2GB")
    return dr.decompose_noise(noise)


@pytest.mark.step_decompose
def test_decompose_recovers_a_planted_plate_variance(tmp_path: Path) -> None:
    """Positive control for decompose.

    Plate offsets of known variance are planted on top of sampling noise of known lfcSE. The
    sample variance of log2FoldChange across plates has expectation sigma2_plate + mean(lfcSE^2)
    exactly -- for any set of per-plate standard errors, not only equal ones -- so subtracting
    the mean squared standard error must return the planted plate variance.
    """
    plate_sd, noise_sd = 0.6, 0.4
    path = _write_fixture_pool(
        tmp_path,
        n_genes=400,
        signal_sd=0.0,
        noise_sd=noise_sd,
        plate_offset_sd=plate_sd,
        plates=tuple(f"P{i}" for i in range(24)),
        seed=13,
    )
    d = _decomposed(path, tmp_path)
    got = float(np.mean(d["sigma2_plate"].to_numpy(dtype=float)))
    # The planted offsets are one draw of 24 values at sd 0.6, so the realised variance is not
    # exactly 0.36; compare against what was actually drawn, recovered from the pool itself.
    raw = pd.read_parquet(path)
    codes, _ = pd.factorize(raw["plate"])
    lfc = raw["log2FoldChange"].to_numpy(dtype=float)
    plate_means = np.bincount(codes, weights=lfc) / np.bincount(codes)
    realised = float(np.var(plate_means, ddof=1))
    assert got == pytest.approx(realised, abs=0.05), (
        f"recovered plate variance {got:.3f} against the realised planted {realised:.3f}"
    )
    frac = float(np.nanmean(d["between_plate_fraction"].to_numpy(dtype=float)))
    expected_frac = realised / (realised + noise_sd**2)
    assert frac == pytest.approx(expected_frac, abs=0.06), (
        f"between-plate fraction {frac:.3f} against planted {expected_frac:.3f}"
    )


@pytest.mark.step_decompose
def test_decompose_floors_at_zero_without_going_negative(tmp_path: Path) -> None:
    """Negative control for decompose: plates differing only by the planted sampling noise
    carry no plate effect. The estimator is a difference of two noisy quantities, so it would
    go negative about half the time unaided; the floor is what makes it a variance, and this
    test is what keeps the floor from being quietly removed."""
    path = _write_fixture_pool(
        tmp_path,
        n_genes=400,
        signal_sd=0.0,
        noise_sd=0.5,
        plate_offset_sd=0.0,
        plates=tuple(f"P{i}" for i in range(24)),
        seed=17,
    )
    d = _decomposed(path, tmp_path)
    sigma2 = d["sigma2_plate"].to_numpy(dtype=float)
    assert (sigma2 >= 0).all(), "the plate component is a variance and cannot be negative"
    assert float(np.mean(sigma2)) < 0.03, (
        f"no plate effect was planted; recovered {float(np.mean(sigma2)):.4f}"
    )


@pytest.mark.step_decompose
def test_decompose_does_not_charge_a_dose_effect_to_plate_noise(tmp_path: Path) -> None:
    """Dose is a grouping key, not pooled. If the decomposition pooled over dose, a screen
    where each dose has a different mean response would report that dose effect as plate
    noise -- and the design's claim about what the ceiling is made of would be wrong."""
    path = _write_fixture_pool(
        tmp_path,
        n_genes=400,
        signal_sd=0.0,
        noise_sd=0.0,  # plates within a dose are identical
        plate_offset_sd=0.0,
        doses=(0.01, 0.1, 1.0),
        plates=("P1", "P2", "P3", "P4"),
        se=0.0,
        seed=19,
    )
    df = pd.read_parquet(path)
    # A large, dose-specific shift shared by every plate at that dose. Vectorized through a
    # lookup on the sorted dose levels rather than a dict map, which pandas-stubs types as a
    # callable-only parameter.
    dose_levels = np.array([0.01, 0.1, 1.0])
    shift = np.array([-2.0, 0.0, 2.0])
    idx = np.searchsorted(dose_levels, df["concentration"].to_numpy(dtype=float))
    df["log2FoldChange"] = df["log2FoldChange"].to_numpy(dtype=float) + shift[idx]
    df.to_parquet(path, index=False)

    d = _decomposed(path, tmp_path)
    worst = float(np.max(d["sigma2_plate"].to_numpy(dtype=float)))
    assert worst < 1e-9, f"a dose effect must not appear as plate variance; recovered {worst:.4f}"


def test_score_positive_planted_reliability_is_recovered(tmp_path: Path) -> None:
    # signal_sd = noise_sd = 1, 8 plates -> 4 per half; half-mean noise sd^2 = 1/4.
    # Expected r = 1 / (1 + 0.25) = 0.8.
    path = _write_fixture_pool(tmp_path, n_genes=600, signal_sd=1.0, noise_sd=1.0)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, _, _ = dr.score_split_half(de, set(de["gene_name"].unique()))
    r = r[np.isfinite(r)]
    assert abs(float(np.mean(r)) - 0.8) < 0.05, f"planted 0.8, recovered {np.mean(r):.3f}"


def test_score_negative_zero_signal_returns_null(tmp_path: Path) -> None:
    path = _write_fixture_pool(tmp_path, signal_sd=0.0, noise_sd=1.0)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    r = r[np.isfinite(r)]
    nulls = dr.stratified_null_draws(piv0, piv1, n_perm=200, seed=0)
    from fmharness.statistics import bootstrap_aggregate_pvalue

    p, _, _ = bootstrap_aggregate_pvalue(float(np.mean(r)), nulls["diff_drug"], r.size)
    assert p > 0.05, f"no planted signal must not clear the null, got p={p}"


def test_restrict_positive_panel_subset_scores_exactly_the_subset(tmp_path: Path) -> None:
    """A panel covering a strict subset of the data's genes restricts scoring to that subset.

    Exercises the real restrict step: `score_split_half`'s `panel` argument is exactly what
    a `--panel-file` resolves to in `main()`, so a fixture set here stands in for one.
    """
    path = _write_fixture_pool(tmp_path)  # default: 300 genes, G0..G299
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    subset = {f"G{k}" for k in range(100)}  # 100 of the fixture's 300 genes
    r, piv0, piv1 = dr.score_split_half(de, subset)
    assert piv0.shape[1] == len(subset), f"scored {piv0.shape[1]} genes, panel had {len(subset)}"
    assert piv1.shape[1] == len(subset)
    assert r.size > 0, "the restricted panel must still leave scoreable pairs"


def test_restrict_negative_disjoint_panel_scores_nothing(tmp_path: Path) -> None:
    """A panel disjoint from the data's genes leaves nothing for the real code to score.

    `score_split_half` itself does not raise on an empty intersection -- it returns an empty
    result -- so this asserts that honestly; `main()` is what turns an empty result into a
    `SystemExit` (the "aborts the run" behavior design.md's restrict control describes).
    """
    path = _write_fixture_pool(tmp_path)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    disjoint = {f"ZZZ{k}" for k in range(10)}  # none of these genes are in the fixture
    r, piv0, piv1 = dr.score_split_half(de, disjoint)
    assert r.size == 0, "a disjoint panel must leave nothing scoreable, not silently score empty"
    assert piv0.shape[1] == 0 and piv1.shape[1] == 0


def test_null_positive_planted_components_recover_the_stratum_ordering(tmp_path: Path) -> None:
    # Drug-shared + line-specific components: matched pairs share both, same-drug
    # mismatches share only the drug component, diff-drug mismatches share nothing.
    path = _write_fixture_pool(
        tmp_path, n_lines=6, n_drugs=4, signal_sd=0.7, drug_sd=0.7, noise_sd=0.7
    )
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2", "D3"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    r = r[np.isfinite(r)]
    nulls = dr.stratified_null_draws(piv0, piv1, n_perm=300, seed=0)
    observed, same_d, diff_d = (
        float(np.mean(r)),
        float(np.mean(nulls["same_drug"])),
        float(np.mean(nulls["diff_drug"])),
    )
    assert observed > same_d + 0.05, f"observed {observed:.3f} !> same_drug {same_d:.3f}"
    assert same_d > diff_d + 0.05, f"same_drug {same_d:.3f} !> diff_drug {diff_d:.3f}"


def test_null_negative_signal_free_strata_sit_at_their_floors(tmp_path: Path) -> None:
    path = _write_fixture_pool(tmp_path, signal_sd=0.0, drug_sd=0.0, noise_sd=1.0)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    _, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    nulls = dr.stratified_null_draws(piv0, piv1, n_perm=200, seed=0)
    for stratum, draws in nulls.items():
        assert abs(float(np.mean(draws))) < 0.05, f"{stratum} floor is not ~0 on noise"


def test_tercile_control_rises_monotonically_with_planted_effect_size(tmp_path: Path) -> None:
    # Three drugs with graded signal size, same noise: split-half r must rise with
    # effect size, tercile 1 -> 3.
    rng = np.random.default_rng(7)
    lines = [f"L{i}" for i in range(6)]
    genes = [f"G{k}" for k in range(400)]
    plates = tuple(f"P{p}" for p in range(8))
    rows = []
    for d, s in (("D0", 0.3), ("D1", 0.8), ("D2", 2.0)):
        for li in lines:
            signal = rng.normal(0.0, s, len(genes))
            for p in plates:
                rows.append(
                    pd.DataFrame(
                        {
                            "Cell_ID_DepMap": li,
                            "drug": d,
                            "gene_name": genes,
                            "log2FoldChange": signal + rng.normal(0.0, 1.0, len(genes)),
                            "plate": p,
                        }
                    )
                )
    pool_dir = tmp_path / "pseudobulk_differential_expression"
    pool_dir.mkdir(parents=True)
    pd.concat(rows, ignore_index=True).to_parquet(
        pool_dir / "train-00000-of-00001.parquet", index=False
    )
    de, _ = dr.build_split_half_frame(
        [str(pool_dir / "train-00000-of-00001.parquet")],
        ["D0", "D1", "D2"],
        None,
        tmp_path / "duck",
        memory_limit="2GB",
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    terc = dr.effect_size_terciles(piv0, piv1, r)
    assert (
        terc["splithalf_mean_r_tercile1"]
        < terc["splithalf_mean_r_tercile2"]
        < terc["splithalf_mean_r_tercile3"]
    ), f"terciles not monotone: {terc}"


def test_per_pair_table_rows_are_keyed_to_their_own_scores(tmp_path: Path) -> None:
    # The exported evidence table must carry each condition's OWN r and effect size.
    # Graded per-drug signal (as in the tercile control) makes misalignment detectable:
    # a scrambled export would break the planted D0 < D1 < D2 ordering in both columns.
    rng = np.random.default_rng(11)
    lines = [f"L{i}" for i in range(6)]
    genes = [f"G{k}" for k in range(400)]
    plates = tuple(f"P{p}" for p in range(8))
    rows = []
    for d, s in (("D0", 0.3), ("D1", 0.8), ("D2", 2.0)):
        for li in lines:
            signal = rng.normal(0.0, s, len(genes))
            for p in plates:
                rows.append(
                    pd.DataFrame(
                        {
                            "Cell_ID_DepMap": li,
                            "drug": d,
                            "gene_name": genes,
                            "log2FoldChange": signal + rng.normal(0.0, 1.0, len(genes)),
                            "plate": p,
                        }
                    )
                )
    pool_dir = tmp_path / "pseudobulk_differential_expression"
    pool_dir.mkdir(parents=True)
    pd.concat(rows, ignore_index=True).to_parquet(
        pool_dir / "train-00000-of-00001.parquet", index=False
    )
    de, _ = dr.build_split_half_frame(
        [str(pool_dir / "train-00000-of-00001.parquet")],
        ["D0", "D1", "D2"],
        None,
        tmp_path / "duck",
        memory_limit="2GB",
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    table = dr.per_pair_table(piv0, piv1, r)

    # column integrity: the r column IS the scored array, in pivot order, keys included
    assert len(table) == len(piv0) == 18
    assert np.array_equal(table["r"].to_numpy(), np.round(r, 4), equal_nan=True)
    assert list(table["patient"]) == list(piv0.index.get_level_values(0))
    assert list(table["drug"]) == list(piv0.index.get_level_values(1))

    # row keying: planted per-drug ordering recovered in both exported quantities
    by_drug_effect = table.groupby("drug")["mean_abs_delta"].mean()
    assert by_drug_effect["D0"] < by_drug_effect["D1"] < by_drug_effect["D2"], (
        f"effect sizes misordered: {dict(by_drug_effect)}"
    )
    by_drug_r = table.groupby("drug")["r"].mean()
    assert by_drug_r["D0"] < by_drug_r["D1"] < by_drug_r["D2"], (
        f"reliabilities misordered: {dict(by_drug_r)}"
    )


def test_null_draw_table_carries_the_draws_its_means_summarize(tmp_path: Path) -> None:
    # The exported floor distributions must BE the draws the summary's floor means average:
    # per-stratum counts and means recomputed from the table must match the dict it came from.
    path = _write_fixture_pool(tmp_path, n_lines=6, n_drugs=4, signal_sd=0.7, drug_sd=0.7)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2", "D3"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    _, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    nulls = dr.stratified_null_draws(piv0, piv1, n_perm=200, seed=0)
    table = dr.null_draw_table(nulls)

    assert set(table["stratum"]) == set(nulls)
    for stratum, draws in nulls.items():
        exported = table[table["stratum"] == stratum]["r"].to_numpy(float)
        assert exported.size == draws.size, (
            f"{stratum}: {exported.size} exported, {draws.size} drawn"
        )
        assert abs(float(np.mean(exported)) - float(np.mean(draws))) < 5e-4, (
            f"{stratum} mean not preserved by the export"
        )
    # the planted ordering must survive into the exported table, not just the dict
    means = table.groupby("stratum")["r"].mean()
    assert means["same_drug"] > means["diff_drug"], f"stratum ordering lost: {dict(means)}"


def test_example_profiles_reproduce_their_own_correlations(tmp_path: Path) -> None:
    # The scatter data must be the scatter the reported r came from: recomputing Pearson
    # from each example's exported two columns must return that example's index r.
    path = _write_fixture_pool(tmp_path, n_lines=6, n_drugs=4, signal_sd=0.7, drug_sd=0.7)
    de, _ = dr.build_split_half_frame(
        [str(path)], ["D0", "D1", "D2", "D3"], None, tmp_path / "duck", memory_limit="2GB"
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    r, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    profiles, index = dr.example_pair_profiles(piv0, piv1, r)

    assert set(index["kind"]) >= {"matched", "same_drug_mismatch", "diff_drug_mismatch"}
    for _, row in index.iterrows():
        sub = profiles[profiles["example_id"] == row["example_id"]]
        assert len(sub) == row["n_genes_shown"], f"{row['example_id']}: gene count disagrees"
        recomputed = np.corrcoef(sub["lfc0"].to_numpy(float), sub["lfc1"].to_numpy(float))[0, 1]
        assert abs(recomputed - row["r_shown"]) < 5e-3, (
            f"{row['example_id']}: exported points give r={recomputed:.4f}, "
            f"index says r_shown={row['r_shown']}"
        )
        # exporting every shared gene (the default): nothing dropped, so the two agree exactly
        assert row["n_genes_shown"] == row["n_genes_full"]
        assert row["r_shown"] == row["r_full"], (
            f"{row['example_id']}: full export must give r_shown == r_full, "
            f"got {row['r_shown']} vs {row['r_full']}"
        )
    # matched examples are ordered by construction (they are drawn at rising quantiles)
    matched = index[index["kind"] == "matched"].sort_values("example_id")
    assert matched["r_full"].is_monotonic_increasing, (
        f"matched examples not ordered: {list(matched['r_full'])}"
    )
    # a mismatched pairing keeps its own two conditions' identities
    mismatch = index[index["kind"] == "same_drug_mismatch"].iloc[0]
    assert mismatch["drug0"] == mismatch["drug1"] and mismatch["patient0"] != mismatch["patient1"]
    cross = index[index["kind"] == "diff_drug_mismatch"].iloc[0]
    assert cross["drug0"] != cross["drug1"] and cross["patient0"] != cross["patient1"]

    # the subsampling path stays available for pools where size forces it: seeded, so a
    # rerun reproduces the file byte for byte, and the shown points keep their own r
    small, small_index = dr.example_pair_profiles(piv0, piv1, r, max_genes=100)
    again, _ = dr.example_pair_profiles(piv0, piv1, r, max_genes=100)
    assert small.equals(again), "subsampling must be reproducible at a fixed seed"
    for _, row in small_index.iterrows():
        assert row["n_genes_shown"] == min(100, row["n_genes_full"])
        sub = small[small["example_id"] == row["example_id"]]
        recomputed = np.corrcoef(sub["lfc0"].to_numpy(float), sub["lfc1"].to_numpy(float))[0, 1]
        assert abs(recomputed - row["r_shown"]) < 5e-3


def test_frame_cache_returns_the_built_frame_and_never_the_wrong_one(tmp_path: Path) -> None:
    # The cache exists so a rerun that only adds an output skips the expensive shard scan.
    # It must hand back exactly what the build produced, and a different input set must not
    # resolve to the same cache entry -- that would silently score the wrong pool.
    import argparse

    path = _write_fixture_pool(tmp_path)
    paths, names = [str(path)], ["D0", "D1", "D2"]
    args = argparse.Namespace(replicate_col=None, frame_cache=str(tmp_path / "cache"))

    built, repl = dr._build_or_load_frame(paths, names, args, tmp_path / "pool")
    cached, repl_again = dr._build_or_load_frame(paths, names, args, tmp_path / "pool")
    assert repl == repl_again == "plate"
    assert built.equals(cached), "the cache must return the frame that was built"

    # a different drug set is a different frame, so it must key to a different cache entry
    assert dr.frame_cache_key(paths, names, None) != dr.frame_cache_key(paths, ["D0"], None)
    assert dr.frame_cache_key(paths, names, None) != dr.frame_cache_key(paths, names, "well")
    subset, _ = dr._build_or_load_frame(paths, ["D0"], args, tmp_path / "pool")
    assert set(subset["drug"].unique()) == {"D0"}, "cache reuse must not cross input sets"


def test_per_gene_reliability_separates_reliable_from_noise_genes(tmp_path: Path) -> None:
    # Half the genes carry pair-specific signal, half are pure noise: the diagnostic
    # must rank the signal genes above the noise genes.
    rng = np.random.default_rng(8)
    lines = [f"L{i}" for i in range(8)]
    drugs = [f"D{j}" for j in range(4)]
    genes = [f"S{k}" for k in range(100)] + [f"N{k}" for k in range(100)]
    plates = tuple(f"P{p}" for p in range(8))
    rows = []
    for li in lines:
        for d in drugs:
            signal = np.concatenate([rng.normal(0.0, 1.5, 100), np.zeros(100)])
            for p in plates:
                rows.append(
                    pd.DataFrame(
                        {
                            "Cell_ID_DepMap": li,
                            "drug": d,
                            "gene_name": genes,
                            "log2FoldChange": signal + rng.normal(0.0, 1.0, 200),
                            "plate": p,
                        }
                    )
                )
    pool_dir = tmp_path / "pseudobulk_differential_expression"
    pool_dir.mkdir(parents=True)
    pd.concat(rows, ignore_index=True).to_parquet(
        pool_dir / "train-00000-of-00001.parquet", index=False
    )
    de, _ = dr.build_split_half_frame(
        [str(pool_dir / "train-00000-of-00001.parquet")],
        drugs,
        None,
        tmp_path / "duck",
        memory_limit="2GB",
    )
    de = de.dropna(subset=["lfc0", "lfc1"])
    _, piv0, piv1 = dr.score_split_half(de, set(de["gene_name"].unique()))
    pg = dr.per_gene_reliability(piv0, piv1, min_pairs=10)
    mean_signal = pg[pg["gene"].str.startswith("S")]["r"].mean()
    mean_noise = pg[pg["gene"].str.startswith("N")]["r"].mean()
    assert mean_signal > mean_noise + 0.3, f"signal {mean_signal:.3f} vs noise {mean_noise:.3f}"


@pytest.mark.step_score
@pytest.mark.parametrize("full_reliability", [0.2, 0.5, 0.8])
def test_spearman_brown_round_trips_a_planted_full_data_reliability(
    tmp_path: Path, full_reliability: float
) -> None:
    """Positive control for score, and the test that checks the correction rather than
    assuming it.

    Plant a pool whose FULL-data reliability is R. With per-condition signal variance 1, n
    plates and per-plate noise variance ``n(1 - R)/R``, the mean over all n plates has
    reliability R and the mean over n/2 has ``R / (2 - R)`` -- which is exactly what
    Spearman-Brown inverts. So the split-half correlation must come back at ``R / (2 - R)`` and
    the corrected value back at R. A correction applied to the wrong quantity, or an off-by-one
    in the half sizes, breaks the round trip.
    """
    n_plates = 16
    noise_sd = float(np.sqrt(n_plates * (1.0 - full_reliability) / full_reliability))
    path = _write_fixture_pool(
        tmp_path,
        n_genes=1500,
        signal_sd=1.0,
        noise_sd=noise_sd,
        plates=tuple(f"P{i}" for i in range(n_plates)),
        seed=23,
    )
    r_all, _, _, _, _ = _scored(path, tmp_path)
    observed_half = float(np.nanmean(r_all))
    expected_half = full_reliability / (2.0 - full_reliability)
    assert observed_half == pytest.approx(expected_half, abs=0.04), (
        f"planted full-data R = {full_reliability}; half correlation read {observed_half:.3f}, "
        f"expected {expected_half:.3f}"
    )
    corrected = dr.spearman_brown_or_nan(observed_half)
    assert corrected == pytest.approx(full_reliability, abs=0.05), (
        f"the correction returned {corrected:.3f} for a planted R of {full_reliability}"
    )


@pytest.mark.step_score
def test_zero_signal_returns_null_and_the_correction_leaves_zero_at_zero(tmp_path: Path) -> None:
    """Negative control for score. A correction that manufactured a ceiling out of nothing
    would be worse than no correction at all, so the identity 2*0/(1+0) = 0 is asserted on the
    real function and on a real signal-free pool."""
    path = _write_fixture_pool(tmp_path, n_genes=600, signal_sd=0.0, noise_sd=1.0, seed=29)
    r_all, _, _, _, _ = _scored(path, tmp_path)
    mean = float(np.nanmean(r_all))
    assert abs(mean) < 0.05, f"signal-free pool read {mean:.3f}"
    assert dr.spearman_brown_or_nan(0.0) == 0.0
    assert abs(dr.spearman_brown_or_nan(mean)) < 0.1
    assert np.isnan(dr.spearman_brown_or_nan(-1.0)), "undefined at r = -1, and guarded"


@pytest.mark.step_score
def test_summary_carries_both_gene_sets_with_their_own_counts_and_mdes() -> None:
    """The summary row is one row carrying two statistics, each with its own condition count,
    p-values and minimum detectable effect. Two files would let one be quoted without the
    other; one row with two prefixed families cannot."""
    rng = np.random.default_rng(31)
    r_all = rng.normal(0.14, 0.06, 1600)
    r_resp = np.concatenate([rng.normal(0.31, 0.10, 900), np.full(700, np.nan)])
    nulls = {
        "any_pair": rng.normal(0.03, 0.05, 500),
        "diff_drug": rng.normal(0.03, 0.05, 500),
        "same_drug": rng.normal(0.07, 0.05, 500),
    }
    even = np.zeros(1600, dtype=bool)
    even[::4] = True
    s = {
        **dr.summarize(r_all, nulls, seed=0, label="all", even_mask=even),
        **dr.summarize(r_resp, nulls, seed=0, label="responder", even_mask=even),
    }
    assert s["all_n_pairs"] == 1600
    assert s["responder_n_pairs"] == 900, "the responder statistic scores fewer conditions"
    assert s["all_n_pairs"] != s["responder_n_pairs"]
    for fam in ("all", "responder"):
        mean = s[f"{fam}_splithalf_mean_r"]
        assert s[f"{fam}_spearman_brown_full"] == pytest.approx(2 * mean / (1 + mean), abs=2e-3)
        assert s[f"{fam}_mde_80_vs_diff_drug"] > 0
        assert s[f"{fam}_mde_80_vs_same_drug"] > 0
        assert not np.isnan(s[f"{fam}_spearman_brown_full_even_plates"])
        assert s[f"{fam}_n_pairs_even"] < s[f"{fam}_n_pairs"]


@pytest.mark.step_null
def test_a_mismatched_responder_draw_uses_the_first_conditions_mask() -> None:
    """A null draw pairs condition i's first group with condition j's second group. The genes
    it scores must be i's responders -- the row whose first group is in play, and so the row
    the selection rule would actually have read. Row j's mask, or the union, would apply a
    different rule to the null than to the observed value.

    Built from disjoint planted responder sets so the three candidate answers -- i's count,
    j's count, and the union's -- are all different numbers.
    """
    n_cond, n_genes = 6, 400
    rng = np.random.default_rng(37)
    idx = pd.MultiIndex.from_arrays(
        [[f"L{i}" for i in range(n_cond)], [f"D{i % 2}" for i in range(n_cond)]],
        names=["patient", "drug"],
    )
    cols = pd.Index([f"G{k}" for k in range(n_genes)], name="gene_name")
    piv0 = pd.DataFrame(rng.normal(size=(n_cond, n_genes)), index=idx, columns=cols)
    piv1 = pd.DataFrame(rng.normal(size=(n_cond, n_genes)), index=idx, columns=cols)
    # Row i selects genes [60*i, 60*i + 60): disjoint blocks, 60 genes each.
    select = np.zeros((n_cond, n_genes), dtype=bool)
    for i in range(n_cond):
        select[i, 60 * i : 60 * i + 60] = True

    nulls = dr.stratified_null_draws(piv0, piv1, n_perm=200, seed=0, min_genes=10, select=select)
    assert nulls["any_pair"].size > 0

    # The scored count is what identifies which mask was used, so read it off the real path:
    # rerun one draw's arithmetic with each candidate mask and require only i's to match 60.
    ii, jj = 0, 3
    a = piv0.to_numpy(dtype=float)[[ii]]
    b = piv1.to_numpy(dtype=float)[[jj]]
    for name, mask, expected in (
        ("first condition", select[[ii]], 60),
        ("second condition", select[[jj]], 60),
        ("union", (select[[ii]] | select[[jj]]), 120),
    ):
        scored = int((np.isfinite(a) & np.isfinite(b) & mask).sum())
        assert scored == expected, f"{name} mask scored {scored}"
    r_first = dr.masked_rowwise_pearson(a, b, 10, select=select[[ii]])
    r_union = dr.masked_rowwise_pearson(a, b, 10, select=select[[ii]] | select[[jj]])
    assert not np.allclose(r_first, r_union), (
        "the three candidate masks must give different answers, or this test cannot tell "
        "which one the null used"
    )
    # And every value the shipped path returned is one the FIRST-condition mask produces, and
    # none is one the union produces. Asserted as set membership rather than by replaying the
    # draw order: stratified_null_draws advances one shared generator across the three strata,
    # so a fresh generator reproduces only the first stratum's picks, and a test that assumed
    # otherwise would be testing its own bookkeeping.
    n = len(piv0)
    gi, gj = np.divmod(np.arange(n * n), n)
    off = gi != gj
    gi, gj = gi[off], gj[off]
    a_all = piv0.to_numpy(dtype=float)[gi]
    b_all = piv1.to_numpy(dtype=float)[gj]
    cand_first = dr.masked_rowwise_pearson(a_all, b_all, 10, select=select[gi])
    cand_union = dr.masked_rowwise_pearson(a_all, b_all, 10, select=select[gi] | select[gj])
    first_set = np.round(cand_first[np.isfinite(cand_first)], 9)
    union_set = np.round(cand_union[np.isfinite(cand_union)], 9)
    for stratum, draws in nulls.items():
        got = np.round(draws, 9)
        assert np.isin(got, first_set).all(), (
            f"{stratum} drew a value the first-condition mask cannot produce"
        )
        assert not np.isin(got, union_set).any(), (
            f"{stratum} drew a value only the union mask produces -- the null is using a "
            "different selection rule from the observed statistic"
        )


def test_summarize_headlines_the_mean_and_reports_both_mdes() -> None:
    rng = np.random.default_rng(9)
    r = rng.normal(0.14, 0.06, 1600)
    nulls = {
        "any_pair": rng.normal(0.03, 0.05, 500),
        "diff_drug": rng.normal(0.03, 0.05, 500),
        "same_drug": rng.normal(0.07, 0.05, 500),
    }
    s = dr.summarize(r, nulls, seed=0)
    assert abs(s["splithalf_mean_r"] - float(np.mean(r))) < 5e-4
    assert (
        abs(s["spearman_brown_full"] - 2 * s["splithalf_mean_r"] / (1 + s["splithalf_mean_r"]))
        < 2e-3
    )
    assert s["p_vs_null"] < 0.01 and s["p_vs_same_drug"] < 0.01
    assert 0 < s["mde_80_vs_diff_drug"] < s["splithalf_mean_r"], "trivially powered here"
    assert 0 < s["mde_80_vs_same_drug"] < s["splithalf_mean_r"]
    assert s["splithalf_median_r"] is not None  # descriptive column retained


@pytest.mark.step_document
def test_main_writes_every_declared_artifact_on_a_synthetic_pool(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """PROCESS section 3: prove the whole run on synthetic data before spending cluster time.

    Runs the REAL ``main`` end to end on a fixture pool and requires the exact set of artifacts
    the design declares. A missing table or figure fails by name here, in seconds, rather than
    at the end of a forty-minute cluster job.
    """
    path = _write_fixture_pool(
        tmp_path,
        n_lines=6,
        n_drugs=4,
        n_genes=400,
        n_responders=120,
        doses=(0.01, 0.1),
        plates=("P1", "P2", "P3", "P4", "P5", "P6"),
        plate_offset_sd=0.3,
        seed=41,
    )
    out = tmp_path / "out"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "delta_reproducibility.py",
            "--local-dir",
            str(path.parent.parent),
            "--out-dir",
            str(out),
            "--min-genes",
            "20",
            "--n-perm",
            "40",
        ],
    )
    dr.main()

    expected = {
        "rung0_reliability.csv",
        "rung0_reliability.params.json",
        "rung0_per_pair_r.csv",
        "rung0_null_draws.csv",
        "rung0_example_pair_profiles.csv.gz",
        "rung0_example_pair_index.csv",
        "rung0_effect_terciles.csv",
        "rung0_mde_curve.csv",
        "rung0_leakage_control.csv",
        "rung0_responder_overlap.csv",
        "rung0_per_gene_reliability.csv",
        "rung0_pool_description.csv",
        "rung0_padj_sample.csv.gz",
        "rung0_noise_decomposition.csv",
        "rung0_noise_decomposition.params.json",
        "rung0_noise_per_gene.csv.gz",
        "rung0_noise_strata.csv",
        "rung0_noise_by_condition.csv",
        "rung0_control_per_pair.csv",
        "rung0_control_noise.csv.gz",
        "audit_checksums.json",
    }
    got = {p.name for p in out.glob("*") if p.is_file()}
    assert expected <= got, f"missing artifacts: {sorted(expected - got)}"

    figures = {p.name for p in (out / "figures").glob("*.png")}
    expected_figures = {
        "01_build.png",
        "02_split.png",
        "03_select.png",
        "04_score.png",
        "05_decompose.png",
        "06_null.png",
        "07_terciles.png",
        "08_power.png",
        "09_per_gene_reliability.png",
    }
    assert expected_figures <= figures, f"missing figures: {sorted(expected_figures - figures)}"

    summary = pd.read_csv(out / "rung0_reliability.csv").iloc[0]
    assert summary["all_n_pairs"] > 0 and summary["responder_n_pairs"] > 0
    assert summary["n_genes"] == 400, "every gene the table carries, no panel and no HVG fallback"

    # The checksum record must cover the artifacts, since the audit cites it and promotion
    # checks against it.
    import hashlib
    import json as _json

    sums = _json.loads((out / "audit_checksums.json").read_text())
    assert expected - {"audit_checksums.json"} <= set(sums)
    name = "rung0_reliability.csv"
    assert sums[name] == hashlib.sha256((out / name).read_bytes()).hexdigest()


@pytest.mark.step_select
def test_responder_mask_cannot_read_the_second_group(tmp_path: Path) -> None:
    """The invariant the whole responder statistic rests on, asserted directly.

    ``padj1`` is carried in the built frame for one reason -- the overlap diagnostic -- and a
    future edit that let it reach selection would inflate every responder number without
    changing a single test's shape. So: flip the second group's adjusted p-values to their
    complement, leaving the first group's untouched, and require the mask to be identical. If
    selection ever reads the second group, this fails.
    """
    path = _write_fixture_pool(tmp_path, n_genes=300, n_responders=90, seed=53)
    de, _ = dr.build_split_half_frame([str(path)], None, None, tmp_path / "duck", "2GB")
    de = de.dropna(subset=["lfc0", "lfc1"])
    panel = set(de["gene_name"].unique())
    assert "padj1" in de.columns, "the overlap diagnostic needs the second group's p-values"
    assert de["padj1"].notna().any(), "and they must actually be populated"

    before = dr.responder_mask(dr.padj_pivot(de, panel))
    flipped = de.copy()
    flipped["padj1"] = 1.0 - flipped["padj1"].to_numpy(dtype=float)
    after = dr.responder_mask(dr.padj_pivot(flipped, panel))
    np.testing.assert_array_equal(before, after)

    # And the diagnostic that IS allowed to read it does move, so the test above is not passing
    # merely because padj1 is inert everywhere.
    o_before = dr.responder_overlap_table(de, panel)["n_second"].to_numpy(dtype=int)
    o_after = dr.responder_overlap_table(flipped, panel)["n_second"].to_numpy(dtype=int)
    assert not np.array_equal(o_before, o_after), (
        "the overlap diagnostic must read padj1, or this test cannot distinguish "
        "'selection ignores it' from 'nothing reads it'"
    )
