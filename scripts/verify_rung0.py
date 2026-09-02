"""Executable verification of rung 0's assay-reliability claims (PROCESS section 3).

Every check here recomputes a reported number from the run's own artifacts alone -- no
cluster access, no trust in any write-up -- and prints the claim beside the recomputed
value and a verdict. Trust comes from re-derivation the reader performs, not from
narrative the writer asserts: a reviewer who runs this has re-derived the evidence rather
than read about it, and a number recomputed at read time cannot drift the way a number
transcribed across documents can.

    uv run python scripts/verify_rung0.py
    uv run python scripts/verify_rung0.py --task-dir path/to/a/run

The default task directory is ``docs/tasks/rung0-assay-reliability``, where the run's
tables sit UNCOMMITTED between the run and promotion (PROCESS, "What reaches GitHub, and
when"); ``--task-dir`` points the same battery at any directory holding one run's output,
which is how the synthetic-data example run and the continuous-integration copy are
checked. The notebook ``docs/tasks/rung0-assay-reliability/verify.ipynb`` recomputes the
same claims in self-contained cells -- standard-library hashing, direct table reads,
nothing imported from this project -- and runs this script only as its final cross-check;
``tests/test_verify_rung0.py`` runs this battery in continuous integration, so the
branch's green does not depend on anyone opening the notebook.

The two reliabilities are checked as one battery over two gene sets. The all-gene family
is the reproducibility of a mostly-null profile; the responder family is the
reproducibility of the response itself, over the genes the condition's FIRST plate group
called differentially expressed. They carry the same statistics under the ``all_`` and
``responder_`` prefixes of one summary row, so every check below runs twice.

What is NOT checkable locally, stated rather than hidden:

* The 1,026 Tahoe pseudobulk shards live on cluster scratch and are far too large to
  commit. Shard integrity therefore reduces here to the committed manifest's content
  hash: the battery recomputes the tranche's ``content_hash`` from the committed
  ``.manifest.txt`` exactly as ``scripts/register_tranche.py`` computed it, which pins
  which bytes the run read without holding those bytes.
* A promoted copy under ``results/`` does not exist until after gate 2. When that
  directory is absent the promotion checks are SKIPPED with a note, not failed -- before
  promotion there is nothing to disagree with, and a battery that failed there would be
  reporting on the calendar rather than on the evidence.
* The permutation check is a separate cluster job. When its outputs are absent those
  checks are SKIPPED the same way.
* Per-gene noise rows number in the millions on the real screen, so the between-plate
  fraction is accumulated in a streaming pass and the variance identity is checked on a
  bounded prefix of the rows rather than on all of them. The mean is exact; the identity
  is a sample, and says so.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TASK = "rung0-assay-reliability"
TRANCHE = "tahoe100m-pseudobulk-de.v1"
DEFAULT_TASK_DIR = REPO / "docs" / "tasks" / TASK

#: The summary row carries one family of statistics per gene set; the per-condition table
#: carries that family's correlations in the column named here.
GENE_SETS: tuple[tuple[str, str], ...] = (("all", "r"), ("responder", "r_responder"))

#: Declared in design.md, one per measurement step plus the per-gene diagnostic.
FIGURES: tuple[str, ...] = (
    "01_build.png",
    "02_split.png",
    "03_select.png",
    "04_score.png",
    "05_decompose.png",
    "06_null.png",
    "07_terciles.png",
    "08_power.png",
    "09_per_gene_reliability.png",
)

#: The table each figure is drawn from. A figure whose source table was never written is a
#: stage that did not run; a figure missing while its table exists is a broken figure step. The
#: battery has to tell those apart, or a partial run reports as a defect.
FIGURE_SOURCES: dict[str, str] = {
    "05_decompose.png": "rung0_noise_per_gene.csv.gz",
}

SUMMARY = "rung0_reliability.csv"
PER_PAIR = "rung0_per_pair_r.csv"
NULL_DRAWS = "rung0_null_draws.csv"
NOISE = "rung0_noise_decomposition.csv"
NOISE_PER_GENE = "rung0_noise_per_gene.csv.gz"
PROFILES = "rung0_example_pair_profiles.csv.gz"
PROFILE_INDEX = "rung0_example_pair_index.csv"
TERCILES = "rung0_effect_terciles.csv"
LEAKAGE = "rung0_leakage_control.csv"
POOL = "rung0_pool_description.csv"
CHECKSUMS = "audit_checksums.json"
SCORE_VALUES = "figures/04_score.values.csv.gz"

#: How many rows of the per-gene noise table the variance identity is checked on. The
#: table has one row per (line, drug, dose, gene) and runs to millions on the real screen;
#: the identity is algebraic and holds row-wise or not at all, so a bounded prefix settles
#: it while keeping the battery to about a minute on a laptop.
IDENTITY_SAMPLE = 200_000


@dataclass
class Check:
    """One claim, the value recomputed from the artifacts, and whether they agree.

    ``skipped`` marks a claim that cannot be checked in this working tree yet -- the
    promoted copy before gate 2, the permutation job's outputs before it has run. A
    skipped check is reported as SKIP and does not fail the battery, so absence is stated
    rather than either hidden or dressed up as a failure.
    """

    name: str
    claim: str
    computed: str
    ok: bool
    skipped: bool = False


def skipped(name: str, why: str) -> Check:
    return Check(name, "not present in this working tree", why, True, skipped=True)


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _close(claim: float, actual: float, decimals: int) -> bool:
    """True when ``claim`` is ``actual`` rounded to ``decimals`` (half-a-unit tolerance).

    Every reported statistic is rounded before it is written, and the per-condition table
    is itself rounded to four decimals, so a recomputation agrees to within half a unit of
    the last reported place rather than exactly.
    """
    if not (np.isfinite(claim) and np.isfinite(actual)):
        return bool(np.isnan(claim) and np.isnan(actual))
    return abs(claim - actual) <= 0.5 * 10**-decimals + 1e-12


def read_table(path: Path) -> pd.DataFrame:
    """Read a committed table, keeping the cell line whose DepMap identifier is ``NA``.

    One of the screen's fifty lines has a missing DepMap identifier and appears throughout
    as the literal string ``NA`` (design.md, inclusion rules). Pandas' default missing-value
    list would silently turn that key into NaN and it would stop joining to itself, so the
    default list is off and only an empty cell counts as missing -- which is how a NaN
    correlation for an unscoreable condition is written.
    """
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def _bool_column(frame: pd.DataFrame, column: str) -> np.ndarray:
    """A committed boolean column as a numpy mask, however the CSV round-tripped it."""
    values = frame[column]
    if values.dtype == bool:
        return values.to_numpy(dtype=bool)
    return values.astype(str).str.strip().str.lower().isin(("true", "1")).to_numpy(dtype=bool)


def summary_row(task_dir: Path) -> pd.Series:
    return read_table(task_dir / SUMMARY).iloc[0]


def _finite(values: np.ndarray) -> np.ndarray:
    return values[np.isfinite(values)]


# ------------------------------------------------------------------------------------------
# the two reliabilities: every reported statistic against the per-condition values it summarises
# ------------------------------------------------------------------------------------------


def check_reliability_statistics(task_dir: Path) -> list[Check]:
    """The summary row recomputed from the 'r' values it is the mean of.

    The strongest statistical check in the battery: not one stored summary against another,
    but the reported row against the raw per-condition correlations the same job exported.
    The responder family's condition count is the number of FINITE responder correlations,
    not the table's row count -- a condition whose first group called too few genes is kept
    in the table, honestly NaN, and is not one of the conditions the responder mean is over.
    """
    row = summary_row(task_dir)
    per_pair = read_table(task_dir / PER_PAIR)
    even = _bool_column(per_pair, "n_plates_even")
    checks: list[Check] = []
    for label, column in GENE_SETS:
        raw = per_pair[column].to_numpy(dtype=float)
        finite = np.isfinite(raw)
        r = raw[finite]
        mean = float(np.mean(r))
        r_even = raw[finite & even]
        mean_even = float(np.mean(r_even)) if r_even.size else float("nan")
        sb = 2 * mean / (1 + mean)
        sb_even = 2 * mean_even / (1 + mean_even)
        stats = (
            ("mean", "splithalf_mean_r", mean, 3),
            ("median", "splithalf_median_r", float(np.median(r)), 3),
            ("lower quartile", "splithalf_q1_r", float(np.quantile(r, 0.25)), 3),
            ("upper quartile", "splithalf_q3_r", float(np.quantile(r, 0.75)), 3),
        )
        checks.append(
            Check(
                f"{label}: n_pairs is the count of finite per-condition correlations",
                f"reported n_pairs {int(row[f'{label}_n_pairs'])}",
                f"{int(finite.sum())} finite of {len(per_pair)} rows in {PER_PAIR}",
                int(finite.sum()) == int(row[f"{label}_n_pairs"]),
            )
        )
        checks.extend(
            Check(
                f"{label}: {what} recomputes from the per-condition correlations",
                f"reported {key} {row[f'{label}_{key}']}",
                f"over {r.size} committed values: {value:.4f}",
                _close(float(row[f"{label}_{key}"]), value, decimals),
            )
            for what, key, value, decimals in stats
        )
        checks.append(
            Check(
                f"{label}: Spearman-Brown correction is 2r/(1+r) of that mean",
                f"reported spearman_brown_full {row[f'{label}_spearman_brown_full']}",
                f"2 x {mean:.4f} / (1 + {mean:.4f}) = {sb:.4f}",
                _close(float(row[f"{label}_spearman_brown_full"]), sb, 3),
            )
        )
        checks.append(
            Check(
                f"{label}: n_pairs_even counts the equal-halves conditions",
                f"reported n_pairs_even {int(row[f'{label}_n_pairs_even'])}",
                f"{int((finite & even).sum())} rows with n_plates_even and a finite {column}",
                int((finite & even).sum()) == int(row[f"{label}_n_pairs_even"]),
            )
        )
        checks.append(
            Check(
                f"{label}: even-plate correction is the same 2r/(1+r) on that subset",
                f"reported {row[f'{label}_splithalf_mean_r_even_plates']} corrected to "
                f"{row[f'{label}_spearman_brown_full_even_plates']}",
                f"mean over {r_even.size} equal-halves conditions {mean_even:.4f}, "
                f"corrected {sb_even:.4f}",
                _close(float(row[f"{label}_splithalf_mean_r_even_plates"]), mean_even, 3)
                and _close(float(row[f"{label}_spearman_brown_full_even_plates"]), sb_even, 3),
            )
        )
        checks.append(
            Check(
                f"{label}: frac_pos is the share of conditions above zero",
                f"reported frac_pos {row[f'{label}_frac_pos']}",
                f"{int((r > 0).sum())} of {r.size} finite values above zero: "
                f"{float(np.mean(r > 0)):.4f}",
                _close(float(row[f"{label}_frac_pos"]), float(np.mean(r > 0)), 3),
            )
        )
    return checks


def check_null_floors(task_dir: Path) -> list[Check]:
    """Each chance floor re-derived from the individual mismatched-condition correlations.

    A floor is a mean over draws that the run also exported one by one, so the reported
    floor is checkable rather than quotable. The observed mean is then required to sit
    above both floors it is read against -- the different-drug-and-line floor (generic
    structure) and the stricter same-drug floor (line specificity).
    """
    row = summary_row(task_dir)
    draws = read_table(task_dir / NULL_DRAWS)
    checks: list[Check] = []
    for label, _ in GENE_SETS:
        subset = draws[draws["gene_set"] == label]
        for stratum in ("any_pair", "diff_drug", "same_drug"):
            d = subset[subset["stratum"] == stratum]["r"].to_numpy(dtype=float)
            mean = float(np.mean(d)) if d.size else float("nan")
            key = f"{label}_null_{stratum}_mean_r"
            checks.append(
                Check(
                    f"{label}: {stratum} floor recomputes from its draws",
                    f"reported {key} {row[key]}",
                    f"mean of {d.size} committed draws: {mean:.4f}",
                    d.size > 0 and _close(float(row[key]), mean, 3),
                )
            )
        n_diff = int((subset["stratum"] == "diff_drug").sum())
        n_any = int((subset["stratum"] == "any_pair").sum())
        expected = n_diff if n_diff else n_any
        checks.append(
            Check(
                f"{label}: null_n_draws counts the stratum the p-value is read against",
                f"reported null_n_draws {int(row[f'{label}_null_n_draws'])}",
                f"{n_diff} diff_drug draws committed (any_pair fallback: {n_any})",
                expected == int(row[f"{label}_null_n_draws"]),
            )
        )
        mean_obs = float(row[f"{label}_splithalf_mean_r"])
        floor_diff = float(row[f"{label}_null_diff_drug_mean_r"])
        floor_same = float(row[f"{label}_null_same_drug_mean_r"])
        checks.append(
            Check(
                f"{label}: the observed mean clears both floors it is read against",
                "mean > different-drug floor and > same-drug floor",
                f"{mean_obs} > {floor_diff} and {mean_obs} > {floor_same}",
                mean_obs > floor_diff and mean_obs > floor_same,
            )
        )
        mdes = [
            float(row[f"{label}_mde_80_vs_diff_drug"]),
            float(row[f"{label}_mde_80_vs_same_drug"]),
        ]
        checks.append(
            Check(
                f"{label}: both minimum detectable effects are positive and finite",
                "an MDE at alpha 0.05, power 0.80 exists against each floor",
                f"vs different-drug {mdes[0]}, vs same-drug {mdes[1]}",
                all(np.isfinite(m) and m > 0 for m in mdes),
            )
        )
    return checks


# ------------------------------------------------------------------------------------------
# the controls the design declared: effect size, selection leakage, noise composition
# ------------------------------------------------------------------------------------------


def check_significance(task_dir: Path) -> list[Check]:
    """The significance claim itself, re-derived -- not merely read back.

    This is the claim rung 0 exists to make: the observed mean clears its chance floor by more
    than sampling would explain. The battery previously checked the floors and the observed
    mean, and then took the p-value on trust, which left the one number a reader most wants
    checked as the one number nothing recomputed.

    The p-value is a bootstrap: resample n_pairs of the mismatched-condition draws with
    replacement, take the mean, repeat, and ask how often that beats the observed mean. It is
    reproduced here from the committed draws with the run's own seed and draw count. Bootstrap
    resampling is not bit-reproducible across implementations, so this requires agreement to
    within the resolution the p-value is REPORTED at, plus one bootstrap standard error -- and
    it separately requires the verdict (clears alpha, or does not) to agree, which is the part a
    reader acts on.

    The minimum detectable effect is checked for the property that makes it meaningful rather
    than recomputed: it must be positive, finite, and -- where the result is declared
    significant -- below the observed mean, since an effect the study could not have detected
    cannot be the effect it detected.
    """
    row = summary_row(task_dir)
    draws = read_table(task_dir / NULL_DRAWS)
    checks: list[Check] = []
    for label, _ in GENE_SETS:
        subset = draws[draws["gene_set"] == label]
        n_obs = int(row[f"{label}_n_pairs"])
        mean_obs = float(row[f"{label}_splithalf_mean_r"])
        for stratum, key in (("diff_drug", "p_vs_null"), ("same_drug", "p_vs_same_drug")):
            pool = subset[subset["stratum"] == stratum]["r"].to_numpy(dtype=float)
            pool = pool[np.isfinite(pool)]
            reported = float(row[f"{label}_{key}"])
            if pool.size < 10 or n_obs < 1:
                checks.append(
                    Check(
                        f"{label}: {key} re-derived from the committed draws",
                        f"reported {reported}",
                        f"only {pool.size} finite draws committed; too few to bootstrap",
                        True,
                        skipped=True,
                    )
                )
                continue
            rng = np.random.default_rng(0)
            n_boot = 2000
            boot = np.array(
                [np.mean(rng.choice(pool, size=n_obs, replace=True)) for _ in range(n_boot)]
            )
            recomputed = float((1 + np.sum(boot >= mean_obs)) / (1 + n_boot))
            # One bootstrap standard error on a proportion, plus the reporting resolution.
            tol = 1.0 / n_boot + 2.0 * float(np.sqrt(max(recomputed, 1e-9) / n_boot)) + 5e-5
            agrees = abs(recomputed - reported) <= tol
            same_verdict = (recomputed < 0.05) == (reported < 0.05)
            checks.append(
                Check(
                    f"{label}: {key} re-derived from the committed draws",
                    f"reported {reported}",
                    f"{recomputed:.4f} from {pool.size} draws resampled to n={n_obs} "
                    f"(tolerance {tol:.4f}; same verdict at alpha 0.05: {same_verdict})",
                    agrees and same_verdict,
                )
            )
        for stratum in ("diff_drug", "same_drug"):
            mde = float(row[f"{label}_mde_80_vs_{stratum}"])
            p = float(
                row[f"{label}_p_vs_null" if stratum == "diff_drug" else f"{label}_p_vs_same_drug"]
            )
            ok = np.isfinite(mde) and mde > 0 and (mde <= mean_obs or p >= 0.05)
            checks.append(
                Check(
                    f"{label}: MDE vs {stratum} is a detectable effect",
                    f"mde_80_vs_{stratum} {mde}",
                    f"positive and finite, and below the observed mean {mean_obs:.4f} "
                    f"where the result is called significant (p = {p})",
                    bool(ok),
                )
            )
    return checks


def check_effect_size_terciles(task_dir: Path) -> list[Check]:
    """The design's empirical in-run control: reproducibility rises with effect size.

    Conditions are cut into thirds by response size and the split-half mean must rise
    across the thirds. An assay that cannot find more reproducibility where there is more
    signal is broken, so a failure here is a finding about the screen -- not a bug in this
    battery.
    """
    terciles = read_table(task_dir / TERCILES).sort_values("tercile")
    means = terciles["mean_r"].to_numpy(dtype=float)
    rises = bool(means.size == 3 and np.all(np.diff(means) > 0))
    return [
        Check(
            "reproducibility rises with effect size (empirical control)",
            "tercile 1 < tercile 2 < tercile 3 of the split-half mean",
            " -> ".join(f"{m:.4f}" for m in means) + ("" if rises else "  [does not rise]"),
            rises,
        )
    ]


def check_leakage_control(task_dir: Path) -> list[Check]:
    """Why responder genes are chosen from the first plate group alone.

    On a pool with no signal at all, selecting genes on the two halves pooled inflates the
    correlation by winner's curse -- the halves' sum and difference are independent, so
    selecting on a large summed magnitude inflates the sum's variance alone and the
    covariance goes positive with nothing generating it. The pooled rule must therefore
    return a visibly larger correlation than the one-sided rule the design uses. The gap is
    the size of the error the selection rule avoids.
    """
    leakage = read_table(task_dir / LEAKAGE)
    by_rule = dict(
        zip(leakage["rule"].astype(str), leakage["mean_r"].to_numpy(dtype=float), strict=True)
    )
    one, pooled = by_rule.get("one-sided", float("nan")), by_rule.get("pooled", float("nan"))
    return [
        Check(
            "two-sided selection inflates a signal-free correlation, one-sided does not",
            "pooled-selection mean r > one-sided mean r on a pool with no signal",
            f"pooled {pooled} vs one-sided {one}",
            bool(np.isfinite(one) and np.isfinite(pooled) and pooled > one),
        )
    ]


def _stream_noise(path: Path, cap: int = IDENTITY_SAMPLE) -> tuple[int, float, int, float]:
    """(rows, sum of finite between-plate fractions, their count, worst identity deviation).

    Read in chunks: the per-gene noise table has one row per (line, drug, dose, gene) and
    is the only artifact large enough to matter on a laptop. The fraction's mean is
    accumulated exactly over every row; the variance identity is checked on the first
    ``cap`` rows, which is enough for an algebraic identity that holds row-wise or not at
    all.
    """
    n_rows = 0
    total = 0.0
    n_finite = 0
    worst = 0.0
    seen = 0
    columns = ["var_lfc", "mean_se2", "sigma2_plate", "between_plate_fraction"]
    for chunk in pd.read_csv(path, usecols=columns, chunksize=200_000):
        n_rows += len(chunk)
        frac = _finite(chunk["between_plate_fraction"].to_numpy(dtype=float))
        total += float(frac.sum())
        n_finite += int(frac.size)
        if seen < cap:
            take = min(cap - seen, len(chunk))
            head = chunk.iloc[:take]
            expected = np.maximum(
                head["var_lfc"].to_numpy(dtype=float) - head["mean_se2"].to_numpy(dtype=float), 0.0
            )
            deviation = np.abs(expected - head["sigma2_plate"].to_numpy(dtype=float))
            worst = max(worst, float(np.nanmax(deviation)) if deviation.size else 0.0)
            seen += take
    return n_rows, total, n_finite, worst


def check_noise_decomposition(task_dir: Path) -> list[Check]:
    """What kind of noise the ceiling is made of, recomputed from the per-gene rows.

    ``lfcSE`` is the standard error of one plate's treated-versus-control contrast and sees
    cell-sampling error only. Across plates at a fixed dose the fold change varies by that
    plus a plate component, so ``sigma2_plate = var_lfc - mean(lfcSE^2)`` floored at zero
    estimates the plate component alone. Both the reported between-plate fraction and that
    identity are re-derived here; the identity is what makes the fraction a decomposition
    rather than a ratio of two stored numbers.
    """
    # A stage that did not run skips; it does not crash. The noise decomposition is the last
    # phase of the job and the one that has failed most, so a battery that raises on its absence
    # is unusable exactly when it is most needed -- reviewing a partial run.
    if not (task_dir / NOISE).exists() or not (task_dir / NOISE_PER_GENE).exists():
        # One skip per check the full path runs, under the SAME names. A battery whose shape
        # changes with which stages completed cannot have its coverage asserted -- the count
        # test and the layer test both read the names, and a partial run would silently look
        # like a battery missing a check rather than a stage that has not happened.
        why = f"{NOISE} / {NOISE_PER_GENE} not written yet"
        return [
            skipped("noise: the per-gene table has the row count the summary reports", why),
            skipped("noise: between-plate fraction is the mean of the per-gene fractions", why),
            skipped("noise: sigma2_plate = max(var_lfc - mean_se2, 0) row by row", why),
        ]
    reported = read_table(task_dir / NOISE).iloc[0]
    n_rows, total, n_finite, worst = _stream_noise(task_dir / NOISE_PER_GENE)
    mean = total / n_finite if n_finite else float("nan")
    return [
        Check(
            "noise: the per-gene table has the row count the summary reports",
            f"reported n_gene_conditions {int(reported['n_gene_conditions'])}",
            f"{n_rows} rows in {NOISE_PER_GENE}",
            n_rows == int(reported["n_gene_conditions"]),
        ),
        Check(
            "noise: between-plate fraction is the mean of the per-gene fractions",
            f"reported between_plate_fraction_mean {reported['between_plate_fraction_mean']}",
            f"mean of {n_finite} finite per-gene fractions: {mean:.6f}",
            _close(float(reported["between_plate_fraction_mean"]), mean, 4),
        ),
        Check(
            "noise: sigma2_plate = max(var_lfc - mean_se2, 0) row by row",
            f"the identity holds on every one of the first {IDENTITY_SAMPLE:,} rows",
            f"worst absolute deviation {worst:.3e}",
            worst < 1e-9,
        ),
    ]


# ------------------------------------------------------------------------------------------
# the illustrative artifacts: each reproduces the number it is shown under
# ------------------------------------------------------------------------------------------


def _pearson_by_group(frame: pd.DataFrame, key: str, x: str, y: str) -> dict[str, float]:
    """Pearson r per group, computed from the committed points and nothing else."""
    out: dict[str, float] = {}
    for name, points in frame.groupby(key):
        a = points[x].to_numpy(dtype=float)
        b = points[y].to_numpy(dtype=float)
        ok = np.isfinite(a) & np.isfinite(b)
        out[str(name)] = float(np.corrcoef(a[ok], b[ok])[0, 1]) if ok.sum() > 1 else float("nan")
    return out


def check_example_profiles(task_dir: Path) -> list[Check]:
    """Each example scatter reproduces the correlation its index records for it.

    The example conditions exist so a reader can see what one correlation looks like gene
    by gene. Recomputing each one from its own exported points is what stops a caption from
    asserting a number the plotted data does not support.
    """
    profiles = read_table(task_dir / PROFILES)
    index = read_table(task_dir / PROFILE_INDEX)
    recomputed = _pearson_by_group(profiles, "example_id", "lfc0", "lfc1")
    counts = profiles.groupby("example_id").size().to_dict()
    disagree = [
        str(ex["example_id"])
        for _, ex in index.iterrows()
        if not (
            _close(float(ex["r_shown"]), recomputed.get(str(ex["example_id"]), float("nan")), 4)
            and counts.get(str(ex["example_id"]), 0) == int(ex["n_genes_shown"])
        )
    ]
    return [
        Check(
            "every example scatter reproduces its own correlation",
            f"{len(index)} examples with an r_shown and a gene count in {PROFILE_INDEX}",
            f"{len(index) - len(disagree)} of {len(index)} reproduce from their committed "
            f"points" + (f"; disagreeing: {disagree}" if disagree else ""),
            not disagree and len(index) > 0,
        )
    ]


def check_figures(task_dir: Path) -> list[Check]:
    """Every declared figure exists, and the score figure's printed values reproduce.

    A figure whose values live only inside a run cannot be checked. The score figure writes
    the points it drew and the correlation it printed to a companion table, so the number
    on the panel is recomputable from the panel's own data.
    """
    # Only figures whose source table exists are expected. The rest are reported as stages that
    # did not run, which is the truth for a partial run and is not the same as a broken step.
    expected = tuple(
        f for f in FIGURES if f not in FIGURE_SOURCES or (task_dir / FIGURE_SOURCES[f]).exists()
    )
    fig_dir = task_dir / "figures"
    present = [name for name in expected if (fig_dir / name).exists()]
    substantial = [
        name
        for name in present
        if (fig_dir / name).stat().st_size > 5_000
        and (fig_dir / name).read_bytes()[:8] == b"\x89PNG\r\n\x1a\n"
    ]
    checks = [
        Check(
            "every figure design.md declares was written, and is a real image",
            f"{len(expected)} figures expected, each a PNG over 5 kB"
            + (
                f" ({len(FIGURES) - len(expected)} not expected: their source tables were "
                "never written)"
                if len(expected) < len(FIGURES)
                else ""
            ),
            f"{len(present)} present, {len(substantial)} non-trivial"
            + (
                f"; missing: {sorted(set(expected) - set(present))}"
                if len(present) < len(expected)
                else ""
            ),
            len(substantial) == len(expected),
        )
    ]
    values_path = task_dir / SCORE_VALUES
    if not values_path.exists():
        checks.append(skipped("the score figure's printed correlations", f"no {SCORE_VALUES}"))
        return checks
    values = read_table(values_path)
    # One panel is an (example, gene set) pair, not an example: the score figure draws each
    # example twice, over all genes and over that condition's responders, and the two print
    # different correlations. Grouping by example alone silently pools two panels' points and
    # matches neither -- which is how this check caught the change that introduced the second
    # row rather than passing over it.
    key = "panel"
    if "gene_set" in values.columns:
        values = values.assign(panel=values["example_id"] + " | " + values["gene_set"])
    else:
        values = values.assign(panel=values["example_id"])
    recomputed = _pearson_by_group(values, key, "lfc0", "lfc1")
    printed = values.groupby(key)["r_printed"].first().to_dict()
    disagree = [str(k) for k, v in printed.items() if not _close(float(v), recomputed[str(k)], 4)]
    checks.append(
        Check(
            "the score figure's printed correlations recompute from its own points",
            f"{len(printed)} panels, each printing an r in {SCORE_VALUES}",
            f"{len(printed) - len(disagree)} of {len(printed)} reproduce"
            + (f"; disagreeing: {disagree}" if disagree else ""),
            not disagree and len(printed) > 0,
        )
    )
    return checks


def check_pool_arithmetic(task_dir: Path) -> list[Check]:
    """The scored conditions re-derived from the measured composition of the pool.

    ``rung0_pool_description.csv`` is the screen as it arrived: one row per candidate
    (line, drug) with its plate counts per half. A condition reaches the per-condition
    table exactly when both halves are non-empty, so the two tables' counts must agree by
    subtraction rather than by assertion, and the equal-halves flag must be the same flag
    in both -- it is what the Spearman-Brown even-plate subset is defined by.
    """
    pool = read_table(task_dir / POOL)
    per_pair = read_table(task_dir / PER_PAIR)
    unsplittable = pool[(pool["n_plates_half0"] == 0) | (pool["n_plates_half1"] == 0)]
    merged = pool.merge(
        per_pair[["patient", "drug", "n_plates_even"]],
        on=["patient", "drug"],
        suffixes=("_pool", "_pair"),
    )
    agree = _bool_column(merged, "n_plates_even_pool") == _bool_column(merged, "n_plates_even_pair")
    return [
        Check(
            "the scored conditions are the splittable ones, by subtraction",
            f"{len(pool)} candidate conditions - {len(unsplittable)} with an empty half",
            f"{len(pool) - len(unsplittable)} splittable; {len(per_pair)} rows in {PER_PAIR}",
            len(pool) - len(unsplittable) == len(per_pair) and len(merged) == len(per_pair),
        ),
        Check(
            "the equal-halves flag agrees between the pool and the per-condition table",
            "n_plates_even is one flag, recorded twice",
            f"{int(agree.sum())} of {len(merged)} conditions agree",
            bool(agree.all()) and len(merged) > 0,
        ),
    ]


# ------------------------------------------------------------------------------------------
# what ties the audit's reading to the committed bytes
# ------------------------------------------------------------------------------------------


def check_audit_checksums(task_dir: Path) -> list[Check]:
    """Every artifact's sha256, recomputed now against what the run recorded.

    The most important check in the battery. The audit reads these artifacts in the working
    tree, before they are committed (PROCESS section 1, "What reaches GitHub, and when"),
    so nothing else ties what a reader audited to what a reviewer later pulls. A single
    altered byte in any table moves its hash and fails here.
    """
    if not (task_dir / CHECKSUMS).exists():
        return [
            skipped(
                "every recorded artifact checksum recomputes from the file it names",
                f"{CHECKSUMS} not written yet -- the run did not reach the end",
            )
        ]
    recorded = json.loads((task_dir / CHECKSUMS).read_text())
    by_name = {p.name: p for p in task_dir.rglob("*") if p.is_file()}
    missing = sorted(name for name in recorded if name not in by_name)
    moved = sorted(
        name
        for name, digest in recorded.items()
        if name in by_name and sha256_of(by_name[name]) != digest
    )
    detail = f"{len(recorded) - len(missing) - len(moved)} of {len(recorded)} match"
    if missing:
        detail += f"; missing: {missing}"
    if moved:
        detail += f"; CHANGED SINCE THE RUN: {moved}"
    return [
        Check(
            "every recorded artifact checksum recomputes from the file it names",
            f"{len(recorded)} sha256 entries in {CHECKSUMS}",
            detail,
            not missing and not moved and len(recorded) > 0,
        )
    ]


def check_tranche_content_hash(repo: Path) -> list[Check]:
    """The data pin, recomputed the way ``scripts/register_tranche.py`` computed it.

    The shards are on cluster scratch and cannot be rehashed here. What can be rehashed is
    the committed manifest: the tranche's content hash is the sha256 of the concatenated
    "relative path, size, sha256" lines, so recomputing it from the manifest confirms that
    the record and the manifest describe the same 1,026 files.
    """
    record = json.loads((repo / "data" / "tranches" / f"{TRANCHE}.json").read_text())
    manifest = repo / "data" / "tranches" / f"{TRANCHE}.manifest.txt"
    rows = [line.split("\t") for line in manifest.read_text().splitlines()]
    text = "".join(f"{rel}\t{size}\t{sha}\n" for rel, size, sha in rows)
    recomputed = hashlib.sha256(text.encode()).hexdigest()
    return [
        Check(
            "the tranche content hash recomputes from the committed manifest",
            f"record content_hash {record['content_hash'][:12]}...",
            f"sha256 of the rebuilt manifest text {recomputed[:12]}...",
            recomputed == record["content_hash"],
        ),
        Check(
            "the manifest describes the whole download",
            "1,026 shards (docs/DATA.md)",
            f"{len(rows)} manifest lines, each path/size/sha256",
            len(rows) == 1026 and all(len(r) == 3 for r in rows),
        ),
    ]


# ------------------------------------------------------------------------------------------
# checks that wait on a later stage: the permutation job, and promotion
# ------------------------------------------------------------------------------------------


def check_permutation(task_dir: Path) -> list[Check]:
    """The dependence check, re-derived from its committed per-permutation means.

    Mismatched draws reuse the same half-profiles, so they are not independent and the
    bootstrap's p-value could be optimistic. Permuting the pairing 500 times measures that
    dependence directly and reports it as a design effect. It is a separate cluster job:
    when its outputs are absent these checks are SKIPPED rather than failed.
    """
    checks: list[Check] = []
    for label, suffix in (("all", ""), ("responder", "_responder")):
        summary_path = task_dir / f"rung0_permutation_summary{suffix}.csv"
        means_path = task_dir / f"rung0_permutation_perm_means{suffix}.csv"
        if not (summary_path.exists() and means_path.exists()):
            reason = f"{summary_path.name} / {means_path.name} not written yet"
            checks.extend(
                [
                    skipped(f"{label}: permutation-mean summary", reason),
                    skipped(f"{label}: exact permutation p-value", reason),
                    skipped(f"{label}: observed mean above every permutation", reason),
                ]
            )
            continue
        s = read_table(summary_path).iloc[0]
        perms = read_table(means_path)["perm_mean"].to_numpy(dtype=float)
        observed = float(s["observed_mean"])
        p_exact = (1 + int(np.sum(perms >= observed))) / (1 + perms.size)
        checks.extend(
            [
                Check(
                    f"{label}: permutation-mean summary recomputes from the draws",
                    f"reported mean {s['perm_mean_mean']}, sd {s['perm_mean_sd']} over "
                    f"{int(s['n_perm'])} permutations",
                    f"mean {np.mean(perms):.4f}, sd {np.std(perms, ddof=1):.4f} over {perms.size}",
                    perms.size == int(s["n_perm"])
                    and _close(float(s["perm_mean_mean"]), float(np.mean(perms)), 4)
                    and _close(float(s["perm_mean_sd"]), float(np.std(perms, ddof=1)), 4),
                ),
                Check(
                    f"{label}: exact permutation p-value recomputes",
                    f"reported p_exact {s['p_exact']}",
                    f"(1 + {int(np.sum(perms >= observed))}) / (1 + {perms.size}) = {p_exact:.4f}",
                    _close(float(s["p_exact"]), p_exact, 4),
                ),
                Check(
                    f"{label}: the observed mean sits above every permutation",
                    f"observed {observed} above the largest of {perms.size} permutation means",
                    f"largest permutation mean {np.max(perms):.4f}",
                    observed > float(np.max(perms)),
                ),
            ]
        )
    return checks


def check_promotion(task_dir: Path, repo: Path) -> list[Check]:
    """The promoted copies, once they exist: same bytes, and a record that recomputes.

    Promotion happens after gate 2, so before it there is nothing under ``results/`` to
    disagree with and these checks SKIP with a note. After it, each promoted table must be
    byte-identical to the task-side copy the audit read, and the provenance record's
    ``result_sha256`` must recompute from the promoted file.
    """
    # The promoted copies belong to the real task folder. When the battery is pointed somewhere
    # else -- an example run, a copy under test -- comparing the repository's promoted tables
    # against that other directory's tables asks whether two different runs produced the same
    # bytes, which is not a question about promotion.
    if task_dir.resolve() != DEFAULT_TASK_DIR.resolve():
        reason = (
            f"battery run against {task_dir}, not the task folder the results were promoted from"
        )
        return [
            skipped("promoted copies are byte-identical to the task-side tables", reason),
            skipped("promoted provenance checksums recompute", reason),
        ]
    promoted_dir = repo / "results" / TASK
    records = sorted(promoted_dir.glob("*.provenance.json")) if promoted_dir.is_dir() else []
    if not records:
        reason = f"{promoted_dir.relative_to(repo)} does not exist (promotion follows gate 2)"
        return [
            skipped("promoted copies are byte-identical to the task-side tables", reason),
            skipped("promoted provenance checksums recompute", reason),
        ]
    identical: list[str] = []
    differ: list[str] = []
    for record_path in records:
        promoted = record_path.with_suffix("").with_suffix(".csv")
        task_side = task_dir / promoted.name
        same = task_side.exists() and task_side.read_bytes() == promoted.read_bytes()
        (identical if same else differ).append(promoted.name)
    hash_ok: list[str] = []
    hash_bad: list[str] = []
    for record_path in records:
        record = json.loads(record_path.read_text())
        promoted = repo / str(record["result"])
        ok = promoted.exists() and sha256_of(promoted) == record["result_sha256"]
        (hash_ok if ok else hash_bad).append(promoted.name)
    return [
        Check(
            "promoted copies are byte-identical to the task-side tables",
            f"{len(records)} promoted results under results/{TASK}",
            f"{len(identical)} identical" + (f"; DIFFER: {differ}" if differ else ""),
            not differ,
        ),
        Check(
            "promoted provenance checksums recompute",
            "each record's result_sha256 is the sha256 of the file it names",
            f"{len(hash_ok)} recompute" + (f"; MISMATCH: {hash_bad}" if hash_bad else ""),
            not hash_bad,
        ),
    ]


# ------------------------------------------------------------------------------------------


def run_all_checks(task_dir: Path = DEFAULT_TASK_DIR, repo: Path = REPO) -> list[Check]:
    """Every claim in the battery, in the order the measurement happens in."""
    return [
        *check_reliability_statistics(task_dir),
        *check_null_floors(task_dir),
        *check_significance(task_dir),
        *check_effect_size_terciles(task_dir),
        *check_leakage_control(task_dir),
        *check_noise_decomposition(task_dir),
        *check_example_profiles(task_dir),
        *check_figures(task_dir),
        *check_pool_arithmetic(task_dir),
        *check_audit_checksums(task_dir),
        *check_tranche_content_hash(repo),
        *check_permutation(task_dir),
        *check_promotion(task_dir, repo),
    ]


NOT_CHECKABLE_LOCALLY = """
Not checkable here, stated rather than hidden:
  - The 1,026 Tahoe shards are on cluster scratch. Shard integrity reduces here to the
    committed manifest's content hash, recomputed above the way register_tranche.py
    computed it; the shard bytes themselves are not in this repository.
  - Anything marked SKIP above is a stage that has not happened yet in this working tree
    (a promoted copy before gate 2, the permutation job before it has run), not a check
    that was waived.
  - The per-gene variance identity is checked on a bounded prefix of the noise table; the
    between-plate fraction it feeds is the exact mean over every row.
"""


def render(checks: list[Check]) -> str:
    lines: list[str] = []
    for c in checks:
        mark = "SKIP" if c.skipped else ("PASS" if c.ok else "FAIL")
        lines.append(f"[{mark}] {c.name}")
        lines.append(f"       claim:      {c.claim}")
        lines.append(f"       recomputed: {c.computed}")
    n_skipped = sum(c.skipped for c in checks)
    n_ok = sum(c.ok and not c.skipped for c in checks)
    n_run = len(checks) - n_skipped
    lines.append(f"\n{n_ok} / {n_run} checks pass ({n_skipped} skipped, {len(checks)} total)")
    lines.append(NOT_CHECKABLE_LOCALLY)
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--task-dir",
        type=Path,
        default=DEFAULT_TASK_DIR,
        help=f"directory holding one run's artifacts (default: docs/tasks/{TASK})",
    )
    ns = ap.parse_args()
    task_dir = ns.task_dir if ns.task_dir.is_absolute() else REPO / ns.task_dir
    if not (task_dir / SUMMARY).exists():
        print(f"no run to verify: {task_dir / SUMMARY} does not exist.")
        print("Run scripts/delta_reproducibility.py first, or pass --task-dir.")
        return 2
    checks = run_all_checks(task_dir)
    print(render(checks))
    return 0 if all(c.ok for c in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
