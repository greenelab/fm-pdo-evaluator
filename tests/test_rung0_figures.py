"""Export controls for rung 0's figures (SPEC rule 4).

These are not appearance tests -- nothing here asks whether a figure is pretty. They ask the two
questions the task design makes of every figure: that the shipped function, given the committed
table it is declared to read, actually writes a real image; and that where a figure PRINTS a
number, that number can be recovered from what the figure exported.

The second question is the load-bearing one. ``fig_score`` prints each example condition's own
split-half correlation on its panel. A printed number no reader can recompute is an assertion,
not evidence, so the function writes the exact points it plotted to a companion CSV; the test
below recomputes Pearson from that CSV and requires it to reproduce the printed value. That is
what turns the design's "recomputable from the points plotted" into a checkable claim.

Every function is also run against a thin or empty table. The figure step runs at the very end of
a long cluster job, and a figure function that raises on one condition, one tercile, or a missing
optional column throws the whole run away.
"""

# pandas ships no PEP-561 type stubs in this environment; under strict mode that turns every
# pandas call site into a cascade of reportUnknown* noise about *pandas'* types, not ours. Same
# suppression, same rationale as the rest of this project's pyright strict config where it
# touches scientific-Python packages -- the rules that check our own code stay on.
# pyright: reportMissingTypeStubs=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownParameterType=false, reportUnknownLambdaType=false

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from fmharness import figures

pytestmark = pytest.mark.known_answer

PNG_MAGIC = b"\x89PNG"


def _assert_is_a_png(path: Path) -> None:
    """A figure that exists but is a stub or a truncated write is not evidence of anything."""
    assert path.exists(), f"no figure written at {path}"
    blob = path.read_bytes()
    assert blob[:4] == PNG_MAGIC, f"{path} is not a PNG (first bytes {blob[:8]!r})"
    assert len(blob) > 1000, f"{path} is only {len(blob)} bytes, too small to hold panels"


# --------------------------------------------------------------------------------------------
# fixtures: small tables with known answers, in the shape the run commits
# --------------------------------------------------------------------------------------------


def _pool(n: int = 12) -> pd.DataFrame:
    """A pool table with a deliberate odd-plate-count minority, so the imbalance is present."""
    rng = np.random.default_rng(0)
    n_plates = rng.choice([3, 4, 6], size=n)
    half0 = n_plates // 2
    return pd.DataFrame(
        {
            "patient": [f"line{i % 4}" for i in range(n)],
            "drug": [f"drug{i % 3}" for i in range(n)],
            "dose": [(0.05, 0.5, 5.0)[i % 3] for i in range(n)],
            "n_rows": rng.integers(500, 1500, size=n),
            "n_plates": n_plates,
            "n_plates_even": half0 == n_plates - half0,
            "n_plates_half0": half0,
            "n_plates_half1": n_plates - half0,
            "frac_untestable": rng.uniform(0.0, 0.4, size=n),
        }
    )


def _delta(n: int = 400, scale: float = 1.0, seed: int = 1) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    return pd.DataFrame({"log2FoldChange": rng.normal(0.0, scale, size=n)})


def _per_pair(n: int = 12) -> pd.DataFrame:
    rng = np.random.default_rng(2)
    return pd.DataFrame(
        {
            "patient": [f"line{i % 4}" for i in range(n)],
            "drug": [f"drug{i % 3}" for i in range(n)],
            "dose": [(0.05, 0.5, 5.0)[i % 3] for i in range(n)],
            "n_genes_scored": rng.integers(60, 900, size=n),
            "mean_abs_half0": rng.uniform(0.1, 1.2, size=n),
            "mean_abs_half1": rng.uniform(0.1, 1.2, size=n),
            "r": rng.uniform(0.05, 0.45, size=n),
            "r_responder": rng.uniform(0.25, 0.75, size=n),
            "n_responders": rng.integers(20, 400, size=n),
        }
    )


def _profiles(n_examples: int = 3, n_genes: int = 120) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Example profiles with a planted, known correlation per example.

    Each example's second half is its first half rotated toward independent noise by a planted
    weight, so the true correlation rises across the examples and the panels span the range the
    design asks the scatters to span.
    """
    rng = np.random.default_rng(3)
    frames: list[pd.DataFrame] = []
    index_rows: list[dict[str, object]] = []
    weights = np.linspace(0.2, 0.9, n_examples)
    for position in range(n_examples):
        signal = rng.normal(size=n_genes)
        noise = rng.normal(size=n_genes)
        weight = float(weights[position])
        lfc0 = signal + rng.normal(scale=0.5, size=n_genes)
        lfc1 = weight * signal + np.sqrt(1.0 - weight**2) * noise
        example_id = f"line{position}|drugA"
        frames.append(
            pd.DataFrame(
                {
                    "example_id": example_id,
                    "gene": [f"G{g:04d}" for g in range(n_genes)],
                    "lfc0": lfc0,
                    "lfc1": lfc1,
                }
            )
        )
        index_rows.append(
            {
                "example_id": example_id,
                "kind": "all_genes" if position % 2 == 0 else "responders",
                "n_genes_full": n_genes,
                "r_full": float(np.corrcoef(lfc0, lfc1)[0, 1]),
                "n_genes_shown": n_genes,
                "r_shown": float(np.corrcoef(lfc0, lfc1)[0, 1]),
            }
        )
    return pd.concat(frames, ignore_index=True), pd.DataFrame(index_rows)


def _summary() -> dict[str, float]:
    return {
        "all_splithalf_mean_r": 0.24,
        "all_spearman_brown_full": 0.387,
        "all_spearman_brown_full_even_plates": 0.401,
        "responder_splithalf_mean_r": 0.51,
        "responder_spearman_brown_full": 0.675,
        "responder_spearman_brown_full_even_plates": 0.69,
        "all_n_pairs": 12.0,
        "design_effect": 1.82,
    }


def _noise(n: int = 300) -> pd.DataFrame:
    rng = np.random.default_rng(4)
    sigma2_plate = rng.gamma(2.0, 0.05, size=n)
    mean_se2 = rng.gamma(2.0, 0.08, size=n)
    return pd.DataFrame(
        {
            "patient": [f"line{i % 4}" for i in range(n)],
            "drug": [f"drug{i % 3}" for i in range(n)],
            "dose": rng.choice([0.1, 1.0], size=n),
            "gene_name": [f"G{i:04d}" for i in range(n)],
            "var_lfc": sigma2_plate + mean_se2,
            "mean_se2": mean_se2,
            "n_plates": 2,
            "base_mean": rng.lognormal(3.0, 1.5, size=n),
            "mean_lfc": rng.normal(0.0, 0.8, size=n),
            "sigma2_plate_signed": sigma2_plate,
            "between_plate_fraction_signed": sigma2_plate / (sigma2_plate + mean_se2),
        }
    )


def _strata() -> pd.DataFrame:
    rows = []
    for e in (1, 2, 3, 4):
        for q in (1, 2, 3, 4):
            rows.append(
                {
                    "expression_quartile": e,
                    "response_quartile": q,
                    "n": 40,
                    "var_lfc_mean": 0.3 + 0.05 * e,
                    "mean_se2_mean": 0.2,
                    "between_plate_fraction_pooled": 1.0 - 0.2 / (0.3 + 0.05 * e),
                    "base_mean_mean": 10.0 * e,
                    "abs_mean_lfc_mean": 0.2 * q,
                }
            )
    return pd.DataFrame(rows)


def _null_draws(n: int = 200) -> pd.DataFrame:
    rng = np.random.default_rng(5)
    frames: list[pd.DataFrame] = []
    centres = {"any_pair": 0.01, "diff_drug": 0.02, "same_drug": 0.08}
    for gene_set in ("all", "responder"):
        for stratum, centre in centres.items():
            frames.append(
                pd.DataFrame(
                    {
                        "gene_set": gene_set,
                        "stratum": stratum,
                        "r": rng.normal(centre, 0.06, size=n),
                    }
                )
            )
    return pd.concat(frames, ignore_index=True)


# --------------------------------------------------------------------------------------------
# build
# --------------------------------------------------------------------------------------------


def test_fig_build_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_build(_pool(), _delta(), _delta(scale=0.4, seed=9), tmp_path / "build.png")
    assert out == tmp_path / "build.png"
    _assert_is_a_png(out)


def test_fig_build_skips_the_untestable_panel_when_the_column_is_absent(tmp_path: Path) -> None:
    pool = _pool().drop(columns=["frac_untestable"])
    out = figures.fig_build(pool, _delta(), _delta(seed=9), tmp_path / "build_no_untestable.png")
    _assert_is_a_png(out)


def test_fig_build_tolerates_one_condition_and_empty_delta_tables(tmp_path: Path) -> None:
    empty = pd.DataFrame({"log2FoldChange": pd.Series(dtype=float)})
    out = figures.fig_build(_pool(n=1), empty, empty, tmp_path / "build_thin.png")
    _assert_is_a_png(out)


def test_fig_build_tolerates_a_wholly_empty_pool(tmp_path: Path) -> None:
    out = figures.fig_build(
        _pool().iloc[0:0], _delta(), _delta(seed=9), tmp_path / "build_empty.png"
    )
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# split
# --------------------------------------------------------------------------------------------


def test_fig_split_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_split(_pool(), _per_pair(), tmp_path / "split.png")
    _assert_is_a_png(out)


def test_fig_split_tolerates_a_single_condition(tmp_path: Path) -> None:
    out = figures.fig_split(_pool(n=1), _per_pair(n=1), tmp_path / "split_one.png")
    _assert_is_a_png(out)


def test_fig_split_tolerates_empty_tables(tmp_path: Path) -> None:
    out = figures.fig_split(_pool().iloc[0:0], _per_pair().iloc[0:0], tmp_path / "split_empty.png")
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# select
# --------------------------------------------------------------------------------------------


def _leakage() -> pd.DataFrame:
    """The known answer the leakage panel exists to show: pooled selection inflates, one-sided
    does not, on a pool with no signal at all."""
    return pd.DataFrame({"rule": ["one-sided", "pooled"], "mean_r": [0.004, 0.31]})


def _padj_sample(n: int = 500) -> pd.DataFrame:
    rng = np.random.default_rng(6)
    return pd.DataFrame({"padj0": rng.beta(0.6, 3.0, size=n)})


def test_fig_select_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_select(_per_pair(), _padj_sample(), _leakage(), tmp_path / "select.png")
    _assert_is_a_png(out)


def test_fig_select_tolerates_a_missing_leakage_row(tmp_path: Path) -> None:
    partial = pd.DataFrame({"rule": ["one-sided"], "mean_r": [0.004]})
    out = figures.fig_select(_per_pair(), _padj_sample(), partial, tmp_path / "select_partial.png")
    _assert_is_a_png(out)


def test_fig_select_tolerates_empty_tables(tmp_path: Path) -> None:
    out = figures.fig_select(
        _per_pair().iloc[0:0],
        _padj_sample().iloc[0:0],
        _leakage().iloc[0:0],
        tmp_path / "select_empty.png",
    )
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# score -- including the recomputable printed correlation
# --------------------------------------------------------------------------------------------


def test_fig_score_writes_a_png_and_a_companion_values_csv(tmp_path: Path) -> None:
    profiles, index = _profiles()
    out = figures.fig_score(
        profiles, index, _per_pair(), _per_pair(n=8), _summary(), tmp_path / "score.png"
    )
    _assert_is_a_png(out)
    values_csv = tmp_path / "score.values.csv.gz"
    assert values_csv.exists(), "the printed correlation has no companion table to be checked from"
    exported = pd.read_csv(values_csv)
    assert list(exported.columns) == [
        "example_id",
        "gene_set",
        "gene",
        "lfc0",
        "lfc1",
        "r_printed",
    ]
    assert set(exported["example_id"]) == set(profiles["example_id"])


def test_fig_score_printed_r_recomputes_from_the_points_it_plotted(tmp_path: Path) -> None:
    """The design's claim, made checkable: the number on the panel comes from the points drawn."""
    profiles, index = _profiles()
    figures.fig_score(profiles, index, _per_pair(), None, _summary(), tmp_path / "score.png")
    exported = pd.read_csv(tmp_path / "score.values.csv.gz")
    for example_id, rows in exported.groupby("example_id"):
        recomputed = float(
            np.corrcoef(rows["lfc0"].to_numpy(dtype=float), rows["lfc1"].to_numpy(dtype=float))[
                0, 1
            ]
        )
        printed = float(rows["r_printed"].to_numpy(dtype=float)[0])
        assert rows["r_printed"].nunique() == 1, f"{example_id} printed more than one correlation"
        assert recomputed == pytest.approx(printed, abs=1e-4), (
            f"{example_id}: printed r {printed} does not recompute from the exported points "
            f"({recomputed})"
        )


def test_fig_score_printed_r_matches_the_planted_correlation(tmp_path: Path) -> None:
    """A planted correlation is recovered by the shipped function, not by a reimplementation."""
    genes = 500
    rng = np.random.default_rng(11)
    lfc0 = rng.normal(size=genes)
    lfc1 = 0.8 * lfc0 + np.sqrt(1.0 - 0.8**2) * rng.normal(size=genes)
    planted = float(np.corrcoef(lfc0, lfc1)[0, 1])
    profiles = pd.DataFrame(
        {
            "example_id": "lineX|drugY",
            "gene": [f"G{g:04d}" for g in range(genes)],
            "lfc0": lfc0,
            "lfc1": lfc1,
        }
    )
    index = pd.DataFrame(
        [
            {
                "example_id": "lineX|drugY",
                "kind": "all_genes",
                "n_genes_full": genes,
                "r_full": planted,
                "n_genes_shown": genes,
                "r_shown": planted,
            }
        ]
    )
    figures.fig_score(profiles, index, _per_pair(), None, _summary(), tmp_path / "planted.png")
    exported = pd.read_csv(tmp_path / "planted.values.csv.gz")
    assert float(exported["r_printed"].to_numpy(dtype=float)[0]) == pytest.approx(planted, abs=1e-4)


def test_fig_score_tolerates_no_control_pool_and_no_examples(tmp_path: Path) -> None:
    empty_profiles = pd.DataFrame(
        {column: pd.Series(dtype=object) for column in ("example_id", "gene")}
    ).assign(lfc0=pd.Series(dtype=float), lfc1=pd.Series(dtype=float))
    empty_index = pd.DataFrame({"example_id": pd.Series(dtype=object)})
    out = figures.fig_score(
        empty_profiles, empty_index, _per_pair(), None, {}, tmp_path / "score_empty.png"
    )
    _assert_is_a_png(out)
    exported = pd.read_csv(tmp_path / "score_empty.values.csv.gz")
    assert list(exported.columns) == [
        "example_id",
        "gene_set",
        "gene",
        "lfc0",
        "lfc1",
        "r_printed",
    ]
    assert len(exported) == 0


def test_fig_score_tolerates_a_single_example_and_no_scored_conditions(tmp_path: Path) -> None:
    profiles, index = _profiles(n_examples=1, n_genes=5)
    out = figures.fig_score(
        profiles, index, _per_pair().iloc[0:0], None, _summary(), tmp_path / "score_one.png"
    )
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# decompose
# --------------------------------------------------------------------------------------------


def test_fig_decompose_writes_a_png_with_its_control_panel(tmp_path: Path) -> None:
    out = figures.fig_decompose(_noise(), _noise(n=150), _strata(), tmp_path / "decompose.png")
    _assert_is_a_png(out)


def test_fig_decompose_tolerates_no_control_pool_and_no_strata(tmp_path: Path) -> None:
    out = figures.fig_decompose(_noise(), None, None, tmp_path / "decompose_no_control.png")
    _assert_is_a_png(out)


def test_fig_decompose_tolerates_too_few_rows_to_stratify(tmp_path: Path) -> None:
    out = figures.fig_decompose(_noise(n=2), None, _strata().iloc[0:0], tmp_path / "thin.png")
    _assert_is_a_png(out)


def test_fig_decompose_tolerates_a_negative_signed_component(tmp_path: Path) -> None:
    """The per-gene estimate is signed and lands below zero for many genes at two plates; the
    histogram clips it and the log-axis scatter drops non-positive rows rather than raising."""
    noise = _noise()
    noise["var_lfc"] = noise["mean_se2"] * 0.5
    noise["sigma2_plate_signed"] = noise["var_lfc"] - noise["mean_se2"]
    noise["between_plate_fraction_signed"] = noise["sigma2_plate_signed"] / noise["var_lfc"]
    out = figures.fig_decompose(noise, None, _strata(), tmp_path / "decompose_negative.png")
    _assert_is_a_png(out)


def test_fig_decompose_tolerates_an_empty_table(tmp_path: Path) -> None:
    out = figures.fig_decompose(_noise().iloc[0:0], None, None, tmp_path / "decompose_empty.png")
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# null
# --------------------------------------------------------------------------------------------


def test_fig_null_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_null(_per_pair(), _null_draws(), tmp_path / "null.png")
    _assert_is_a_png(out)


def test_fig_null_tolerates_a_single_gene_set(tmp_path: Path) -> None:
    draws = _null_draws()
    out = figures.fig_null(
        _per_pair(), draws.loc[draws["gene_set"] == "all"], tmp_path / "null_one.png"
    )
    _assert_is_a_png(out)


def test_fig_null_tolerates_empty_draws(tmp_path: Path) -> None:
    out = figures.fig_null(_per_pair(), _null_draws().iloc[0:0], tmp_path / "null_empty.png")
    _assert_is_a_png(out)


def test_fig_null_tolerates_no_matched_pairs_and_no_draws(tmp_path: Path) -> None:
    out = figures.fig_null(
        _per_pair().iloc[0:0], _null_draws().iloc[0:0], tmp_path / "null_nothing.png"
    )
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# permutation against bootstrap -- the design effect
# --------------------------------------------------------------------------------------------


def _perm_means(n: int = 500) -> pd.DataFrame:
    rng = np.random.default_rng(7)
    return pd.DataFrame({"perm_mean": rng.normal(0.02, 0.03, size=n)})


def test_fig_permutation_vs_bootstrap_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_permutation_vs_bootstrap(
        _perm_means(), _null_draws(), _summary(), tmp_path / "design_effect.png"
    )
    _assert_is_a_png(out)


def test_fig_permutation_vs_bootstrap_tolerates_a_summary_without_a_design_effect(
    tmp_path: Path,
) -> None:
    summary = _summary()
    del summary["design_effect"]
    out = figures.fig_permutation_vs_bootstrap(
        _perm_means(), _null_draws(), summary, tmp_path / "design_effect_bare.png"
    )
    _assert_is_a_png(out)


def test_fig_permutation_vs_bootstrap_tolerates_empty_inputs(tmp_path: Path) -> None:
    out = figures.fig_permutation_vs_bootstrap(
        pd.DataFrame({"perm_mean": pd.Series(dtype=float)}),
        _null_draws().iloc[0:0],
        {},
        tmp_path / "design_effect_empty.png",
    )
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# terciles and power
# --------------------------------------------------------------------------------------------


def _terciles() -> pd.DataFrame:
    """The known answer: the mean rises across the thirds under both rankings."""
    return pd.DataFrame(
        {
            "ranked_by": ["half0"] * 3 + ["half1"] * 3,
            "tercile": [1, 2, 3, 1, 2, 3],
            "mean_r": [0.11, 0.24, 0.41, 0.10, 0.22, 0.39],
            "ci_lo": [0.05, 0.18, 0.33, 0.04, 0.16, 0.31],
            "ci_hi": [0.17, 0.30, 0.49, 0.16, 0.28, 0.47],
            "n": [12, 12, 11, 12, 12, 11],
        }
    )


def _dose_strata() -> pd.DataFrame:
    rows = []
    for gene_set, base in (("all", 0.12), ("responder", 0.45)):
        for dose, shift in (("0.05", -0.02), ("0.5", 0.0), ("5.0", 0.03), ("all", 0.01)):
            rows.append(
                {
                    "dose": dose,
                    "weighting": "per_triple",
                    "gene_set": gene_set,
                    "n_pairs": 4,
                    "splithalf_mean_r": base + shift,
                    "splithalf_median_r": base + shift,
                    "spearman_brown_full": 2 * (base + shift) / (1 + base + shift),
                }
            )
        rows.append(
            {
                "dose": "all",
                "weighting": "per_line_drug",
                "gene_set": gene_set,
                "n_pairs": 4,
                "splithalf_mean_r": base,
                "splithalf_median_r": base,
                "spearman_brown_full": 2 * base / (1 + base),
            }
        )
    return pd.DataFrame(rows)


def test_fig_terciles_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_terciles(_terciles(), tmp_path / "terciles.png")
    _assert_is_a_png(out)


def test_fig_terciles_tolerates_a_single_tercile(tmp_path: Path) -> None:
    out = figures.fig_terciles(_terciles().iloc[0:1], tmp_path / "terciles_one.png")
    _assert_is_a_png(out)


def test_fig_terciles_tolerates_missing_intervals_and_an_empty_table(tmp_path: Path) -> None:
    bare = _terciles().assign(ci_lo=np.nan, ci_hi=np.nan)
    _assert_is_a_png(figures.fig_terciles(bare, tmp_path / "terciles_bare.png"))
    _assert_is_a_png(figures.fig_terciles(_terciles().iloc[0:0], tmp_path / "terciles_empty.png"))


def _mde_curve() -> pd.DataFrame:
    n_pairs = np.array([5, 10, 20, 40, 80])
    frames: list[pd.DataFrame] = []
    for gene_set, observed_n in (("all", 20), ("responder", 10)):
        frames.append(
            pd.DataFrame(
                {
                    "n_pairs": n_pairs,
                    "mde": 0.6 / np.sqrt(n_pairs),
                    "gene_set": gene_set,
                    "observed": n_pairs == observed_n,
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def test_fig_power_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_power(_mde_curve(), tmp_path / "power.png")
    _assert_is_a_png(out)


def test_fig_power_tolerates_a_curve_with_no_observed_row(tmp_path: Path) -> None:
    curve = _mde_curve().assign(observed=False)
    out = figures.fig_power(curve, tmp_path / "power_unmarked.png")
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# dose
# --------------------------------------------------------------------------------------------


def test_fig_dose_writes_a_png(tmp_path: Path) -> None:
    out = figures.fig_dose(_per_pair(), _dose_strata(), tmp_path / "dose.png")
    _assert_is_a_png(out)


def test_fig_dose_tolerates_one_dose_level_and_no_strata(tmp_path: Path) -> None:
    per_pair = _per_pair().assign(dose=5.0)
    out = figures.fig_dose(per_pair, pd.DataFrame(), tmp_path / "dose_one.png")
    _assert_is_a_png(out)


def test_fig_dose_tolerates_empty_tables_and_no_responder_column(tmp_path: Path) -> None:
    _assert_is_a_png(figures.fig_dose(_per_pair().iloc[0:0], _dose_strata(), tmp_path / "e.png"))
    bare = _per_pair().drop(columns=["r_responder"])
    _assert_is_a_png(figures.fig_dose(bare, _dose_strata(), tmp_path / "dose_bare.png"))


def test_fig_dose_uses_one_fixed_hue_per_dose_level_low_to_high() -> None:
    palette = figures.dose_palette(["5.0", "0.05", "0.5", "0.05"])
    assert list(palette) == ["0.05", "0.5", "5.0"], "sorted by dose, not by appearance"
    assert len(set(palette.values())) == 3


def test_fig_power_tolerates_an_empty_curve(tmp_path: Path) -> None:
    out = figures.fig_power(_mde_curve().iloc[0:0], tmp_path / "power_empty.png")
    _assert_is_a_png(out)


# --------------------------------------------------------------------------------------------
# the module's own contract
# --------------------------------------------------------------------------------------------


def test_every_figure_function_creates_the_output_directory_it_is_given(tmp_path: Path) -> None:
    """The run writes into ``<out-dir>/figures/``, which may not exist when the step starts."""
    nested = tmp_path / "figures" / "deeper"
    out = figures.fig_terciles(_terciles(), nested / "terciles.png")
    _assert_is_a_png(out)


def test_no_emojis_in_the_figures_module() -> None:
    """Project rule: no emojis anywhere, including in axis labels a reader will read."""
    source = (Path(figures.__file__)).read_text(encoding="utf-8")
    offending = [character for character in source if ord(character) > 0x2100]
    assert not offending, f"non-ascii pictographic characters in figures.py: {offending}"


def test_fig_score_draws_the_examples_over_both_gene_sets(tmp_path: Path) -> None:
    """The design declares the example scatters drawn twice -- all genes, then that condition's
    responders -- because the second is what a rung scoring responding genes is read against,
    and seeing them side by side is how a reader judges whether restricting to responders bought
    signal or only removed the easy agreement of shared zeros.

    The marking comes from the committed profile table. When the table carries no responder
    column the figure falls back to one row rather than inventing a panel, which is the state a
    run without responder selection leaves behind.
    """
    profiles, index = _profiles()
    rng = np.random.default_rng(3)
    profiles = profiles.assign(is_responder=rng.random(len(profiles)) < 0.4)
    figures.fig_score(profiles, index, _per_pair(), None, _summary(), tmp_path / "score.png")
    exported = pd.read_csv(tmp_path / "score.values.csv.gz")
    assert set(exported["gene_set"]) == {"all genes", "responders"}
    for example_id, part in exported.groupby("example_id"):
        n_all = int((part["gene_set"] == "all genes").sum())
        n_resp = int((part["gene_set"] == "responders").sum())
        assert n_resp < n_all, f"{example_id}: responders must be a subset, got {n_resp}/{n_all}"
        # And each panel's printed r must still recompute from the points that panel plotted.
        for gene_set, panel in part.groupby("gene_set"):
            if len(panel) > 2:
                got = float(np.corrcoef(panel["lfc0"], panel["lfc1"])[0, 1])
                assert got == pytest.approx(float(panel["r_printed"].iloc[0]), abs=1e-4), (
                    f"{example_id} / {gene_set}: printed r does not match its own points"
                )

    # No responder column: one row of scatters, and the companion table says so.
    figures.fig_score(_profiles()[0], index, _per_pair(), None, _summary(), tmp_path / "score2.png")
    plain = pd.read_csv(tmp_path / "score2.values.csv.gz")
    assert set(plain["gene_set"]) == {"all genes"}
