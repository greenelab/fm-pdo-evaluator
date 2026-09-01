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
) -> Path:
    """A synthetic replicate pool in the DE table's own shape, one parquet file.

    Per (line, drug, gene): a fixed pair-specific signal (sd ``signal_sd``), an optional
    drug-shared component (sd ``drug_sd``), plus independent per-plate noise (sd
    ``noise_sd``). Expected split-half r over genes, as plate count grows:
    (signal_sd^2 + drug_sd^2) / (signal_sd^2 + drug_sd^2 + noise_sd^2 / plates_per_half).
    """
    rng = np.random.default_rng(seed)
    lines = [f"L{i}" for i in range(n_lines)]
    drugs = [f"D{j}" for j in range(n_drugs)]
    genes = [f"G{k}" for k in range(n_genes)]
    drug_eff = {d: rng.normal(0.0, drug_sd, n_genes) for d in drugs}
    rows = []
    for li in lines:
        for d in drugs:
            signal = rng.normal(0.0, signal_sd, n_genes) + drug_eff[d]
            for p in plates:
                lfc = signal + rng.normal(0.0, noise_sd, n_genes)
                rows.append(
                    pd.DataFrame(
                        {
                            "Cell_ID_DepMap": li,
                            "drug": d,
                            "gene_name": genes,
                            "log2FoldChange": lfc,
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
