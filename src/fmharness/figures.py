"""Figures for rung 0, the assay-reliability measurement.

Rung 0 asks how much of a drug's transcriptional response survives being measured twice. The
screen's replicate plates for one (patient line, drug) condition are split into two groups, each
group is differentially expressed against its own controls, and the two half-profiles are
correlated -- once over all genes, once over the genes the FIRST group alone called responders.
The figures here are how a reader checks that measurement rather than taking it on trust.

Two rules from the task design hold for every function in this module.

**A figure is drawn from a committed table.** Each function takes an already-computed table and
an output path, and reads nothing else -- no database, no raw counts, no statistic recomputed
from the screen. If a number is drawn on a figure it came from a table that ships beside the
figure, so a reader can recompute it. One case looks like an exception and is in fact the rule
enforced: the example scatters print each condition's own correlation, and because that number is
computed here it is written back out with the exact points it was computed from, so the printed
value stays checkable.

**A figure that has a control shows it beside the real data on shared axes.** Real data alone
says what the screen looks like; real data beside a pool with a planted answer says whether the
machinery reads it correctly. Every such panel names which side is which in its title or legend.

Every function returns the path it wrote, and every function tolerates a thin or empty table
without raising: the figure step runs at the end of a long cluster job, and a crash there throws
away the whole run.
"""

# pandas ships no PEP-561 type stubs in this environment; under strict mode that turns every
# pandas call site into a cascade of reportUnknown* noise about *pandas'* types, not ours. Same
# suppression, same rationale as the rest of this project's pyright strict config where it
# touches scientific-Python packages -- the rules that check our own code stay on.
# pyright: reportMissingTypeStubs=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownParameterType=false, reportUnknownLambdaType=false

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any, cast

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure

# The non-interactive backend is chosen BEFORE pyplot is imported, so a figure renders
# identically on a headless cluster node and on a laptop.
matplotlib.use("Agg")

import matplotlib.pyplot as plt

__all__ = [
    "SCORING_THRESHOLD",
    "fig_build",
    "fig_decompose",
    "fig_null",
    "fig_permutation_vs_bootstrap",
    "fig_power",
    "fig_score",
    "fig_select",
    "fig_split",
    "fig_terciles",
]

# The scoring threshold the design fixes: a condition is scored only when at least this many
# genes are finite in both halves, and a responder set is used only when it reaches the same
# size. A correlation over a handful of genes is noise with a number attached.
SCORING_THRESHOLD = 50

_REAL_COLOR = "tab:blue"
_RESPONDER_COLOR = "tab:red"
_CONTROL_COLORS = ("#8c564b", "#7f7f7f", "#bcbd22")
_CONTROL_COLOR = "tab:orange"
_LFC_UNITS = "log2 fold change (drug vs control)"


# --------------------------------------------------------------------------------------------
# Reading a table. Every column any figure draws arrives through one of these four accessors,
# as plain numpy, which keeps the drawing code vectorised and keeps pandas' untyped surface out
# of the type checker's way.
# --------------------------------------------------------------------------------------------


def _numeric(frame: pd.DataFrame, column: str) -> np.ndarray:
    """One column as floats, aligned row-for-row; all-``nan`` when the column is absent.

    Alignment matters wherever two columns are read together -- the stratified panels correlate
    a quantity with a stratifier row by row, so nothing may be dropped here.
    """
    if column not in frame.columns:
        return np.full(len(frame), np.nan, dtype=float)
    coerced = cast(Any, pd.to_numeric(cast(Any, frame)[column], errors="coerce"))
    return np.asarray(coerced.to_numpy(dtype=float), dtype=float)


def _finite(frame: pd.DataFrame, column: str) -> np.ndarray:
    """The finite values of one column, or an empty array when the column is absent.

    Absence is tolerated rather than raised on: some tables carry optional columns (the
    untestable-gene fraction, for one), and a missing optional column should cost a panel rather
    than the run.
    """
    values = _numeric(frame, column)
    return values[np.isfinite(values)]


def _labels(frame: pd.DataFrame, column: str) -> np.ndarray:
    """One column as strings, aligned row-for-row; empty strings when the column is absent."""
    if column not in frame.columns:
        return np.array([""] * len(frame), dtype=str)
    raw = cast(Any, frame)[column].tolist()
    return np.array([str(value) for value in raw], dtype=str)


def _flags(frame: pd.DataFrame, column: str) -> np.ndarray:
    """One column as booleans, aligned row-for-row; all false when the column is absent."""
    if column not in frame.columns:
        return np.zeros(len(frame), dtype=bool)
    filled = cast(Any, frame)[column].fillna(False)
    return np.asarray(filled.to_numpy(dtype=bool), dtype=bool)


def _ordered_unique(values: np.ndarray) -> list[str]:
    """Distinct labels in the order they first appear, so panels keep the table's own order."""
    return [str(value) for value in dict.fromkeys(values.tolist())]


# --------------------------------------------------------------------------------------------
# Small drawing helpers -- bins, empty panels, quartile strata, writing the file out.
# --------------------------------------------------------------------------------------------


def _shared_bins(arrays: Sequence[np.ndarray], n_bins: int = 30) -> list[float]:
    """Common bin edges spanning every array, so overlaid histograms are comparable.

    Two histograms drawn on separate binnings cannot be read against each other, which defeats
    the point of putting a control beside the real data.
    """
    filled = [array for array in arrays if array.size]
    if not filled:
        return _even_bins(0.0, 1.0, n_bins)
    stacked = np.concatenate(filled)
    low, high = float(np.min(stacked)), float(np.max(stacked))
    if not (np.isfinite(low) and np.isfinite(high)) or high <= low:
        pad = max(abs(low), 1.0) * 0.05
        return _even_bins(low - pad, low + pad, n_bins)
    return _even_bins(low, high, n_bins)


def _even_bins(low: float, high: float, n_bins: int) -> list[float]:
    """``n_bins`` equal-width bins between two edges."""
    return [float(edge) for edge in np.linspace(low, high, n_bins + 1)]


def _integer_bins(arrays: Sequence[np.ndarray]) -> list[float]:
    """Bin edges centred on the integers, for counts of plates, genes or conditions.

    Plate counts are small integers, and the whole point of the split figure is that one plate
    and two plates are visibly different bars, which continuous bins would smear together.
    """
    filled = [array for array in arrays if array.size]
    if not filled:
        return [-0.5, 0.5]
    stacked = np.concatenate(filled)
    low = int(np.floor(float(np.min(stacked))))
    high = int(np.ceil(float(np.max(stacked))))
    return [float(edge) for edge in np.arange(low - 0.5, high + 1.5, 1.0)]


def _pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation over the pairs finite in both arrays; ``nan`` where undefined."""
    keep = np.isfinite(x) & np.isfinite(y)
    x_kept, y_kept = x[keep], y[keep]
    if x_kept.size < 2:
        return float("nan")
    if float(np.std(x_kept)) == 0.0 or float(np.std(y_kept)) == 0.0:
        return float("nan")
    return float(np.corrcoef(x_kept, y_kept)[0, 1])


def _number(value: object) -> float:
    """A summary entry as a float, or ``nan`` when it is missing or is not a number."""
    if value is None or isinstance(value, bool):
        return float("nan")
    if isinstance(value, int | float | np.floating | np.integer):
        return float(value)
    return float("nan")


def _legend(ax: Axes, *, fontsize: int = 7, ncol: int = 1) -> None:
    """Draw a legend only when something is labelled, so an empty panel stays quiet."""
    handles, _ = ax.get_legend_handles_labels()
    if handles:
        ax.legend(frameon=False, fontsize=fontsize, ncol=ncol)


def _note_empty(ax: Axes, message: str) -> None:
    """Say plainly that a panel had nothing to draw, rather than showing bare axes."""
    ax.text(0.5, 0.5, message, ha="center", va="center", fontsize=8, transform=ax.transAxes)


def _category_bar(ax: Axes, frame: pd.DataFrame, column: str, xlabel: str) -> None:
    """How many conditions each patient line, or each drug, contributes, tallest first.

    Reads as the screen's balance. If a handful of lines or drugs carry most of the conditions,
    any reliability measured over the pool describes those few more than it describes the
    screen, and that limits how far the number generalises.
    """
    ax.set_ylabel("conditions (count)")
    if column not in frame.columns or len(frame) == 0:
        _note_empty(ax, f"no {column} values in the pool")
        ax.set_xlabel(xlabel)
        return
    unique, counts = np.unique(_labels(frame, column), return_counts=True)
    order = np.argsort(-counts, kind="stable")
    unique, counts = unique[order], counts[order]
    positions = np.arange(unique.size, dtype=float)
    ax.bar(positions, counts.astype(float), color=_REAL_COLOR, alpha=0.85)
    ax.set_xticks(positions)
    if unique.size <= 20:
        ax.set_xticklabels([str(label) for label in unique.tolist()], rotation=90, fontsize=6)
    else:
        ax.set_xticklabels([])
        xlabel = f"{xlabel} ({unique.size} levels, labels omitted)"
    ax.set_xlabel(xlabel)


def _quartile_profile(
    quantity: np.ndarray, stratifier: np.ndarray
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    """Mean, standard error and count of ``quantity`` within quartiles of ``stratifier``.

    Returns empty arrays when there are too few rows, or too few distinct stratifier values, to
    cut into quartiles at all -- a stratification of three rows is not a stratification.
    """
    keep = np.isfinite(quantity) & np.isfinite(stratifier)
    quantity, stratifier = quantity[keep], stratifier[keep]
    empty = np.empty(0, dtype=float)
    if quantity.size < 4 or np.unique(stratifier).size < 2:
        return [], empty, empty, empty
    frame = pd.DataFrame({"quantity": quantity, "stratifier": stratifier})
    quartile = cast(Any, pd.qcut(cast(Any, frame)["stratifier"], 4, duplicates="drop"))
    grouped = cast(Any, frame).groupby(quartile, observed=True)["quantity"]
    aggregated = grouped.agg(["mean", "sem", "size"])
    labels = [f"{edge.left:.3g} to {edge.right:.3g}" for edge in aggregated.index]
    means = np.asarray(aggregated["mean"].to_numpy(dtype=float), dtype=float)
    sems = np.asarray(aggregated["sem"].to_numpy(dtype=float), dtype=float)
    counts = np.asarray(aggregated["size"].to_numpy(dtype=float), dtype=float)
    return labels, means, np.nan_to_num(sems, nan=0.0), counts


def _quartile_panel(
    ax: Axes,
    labels: list[str],
    means: np.ndarray,
    sems: np.ndarray,
    counts: np.ndarray,
    xlabel: str,
    title: str,
) -> None:
    """Mean between-plate share within each quartile, with its standard error."""
    if means.size:
        positions = np.arange(means.size, dtype=float)
        ax.errorbar(
            positions, means, yerr=sems, fmt="o-", color=_REAL_COLOR, capsize=3, lw=1.2, ms=5
        )
        ax.set_xticks(positions)
        ticks = [
            f"{label}\nn={int(count)}"
            for label, count in zip(labels, counts.tolist(), strict=False)
        ]
        ax.set_xticklabels(ticks, fontsize=6)
    else:
        _note_empty(ax, "too few rows to cut into quartiles")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("mean between-plate share of variance")
    ax.set_title(title, fontsize=9)


def _finish(fig: Figure, out: Path) -> Path:
    """Write the figure and return its path, which is what every public function returns."""
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# --------------------------------------------------------------------------------------------
# build -- the screen as it arrives
# --------------------------------------------------------------------------------------------


def fig_build(
    pool: pd.DataFrame,
    delta_real: pd.DataFrame,
    delta_synthetic: pd.DataFrame,
    out: Path,
) -> Path:
    """What the screen looks like before anything is measured from it.

    A reliability is only as general as the pool it was measured on, so the reviewer meets the
    pool first: how many replicate plates each (patient line, drug) condition has, and how the
    conditions are spread over lines and over drugs. A pool where a few lines carry most of the
    conditions supports a narrower claim than an even one, and a pool of mostly two-plate
    conditions cannot support a split-half at all.

    The fold-change panel puts the real screen's differential expression beside a synthetic
    control pool with planted structure, on shared bins. That is this step's control: the
    synthetic pool went through the same builder, so if its distribution did not come out with
    the shape that was planted, nothing downstream is reading what it claims to read.

    The last panel, drawn only when the pool table carries ``frac_untestable``, is the fraction
    of genes DESeq2 could not test at all (zero base mean). Genes that were never testable are
    not evidence about reproducibility either way, and a condition where that fraction is large
    is a condition whose correlation rests on few genes.
    """
    has_untestable = "frac_untestable" in pool.columns
    fig = plt.figure(figsize=(15.0, 8.0), layout="constrained")
    grid = fig.add_gridspec(2, 3)

    ax_plates = fig.add_subplot(grid[0, 0])
    plates = _finite(pool, "n_plates")
    if plates.size:
        ax_plates.hist(plates, bins=_integer_bins([plates]), color=_REAL_COLOR, alpha=0.85)
    else:
        _note_empty(ax_plates, "no conditions in the pool")
    ax_plates.set_xlabel("replicate plates per condition (count)")
    ax_plates.set_ylabel("conditions (count)")
    ax_plates.set_title("(a) replicate depth per condition", fontsize=9)

    ax_patient = fig.add_subplot(grid[0, 1])
    _category_bar(ax_patient, pool, "patient", "patient line")
    ax_patient.set_title("(b) conditions per patient line", fontsize=9)

    ax_drug = fig.add_subplot(grid[0, 2])
    _category_bar(ax_drug, pool, "drug", "drug")
    ax_drug.set_title("(c) conditions per drug", fontsize=9)

    ax_delta = fig.add_subplot(grid[1, 0:2] if has_untestable else grid[1, 0:3])
    real = _finite(delta_real, "log2FoldChange")
    synthetic = _finite(delta_synthetic, "log2FoldChange")
    if real.size or synthetic.size:
        bins = _shared_bins([real, synthetic], n_bins=60)
        if real.size:
            ax_delta.hist(
                real, bins=bins, density=True, alpha=0.6, color=_REAL_COLOR, label="real screen"
            )
        if synthetic.size:
            ax_delta.hist(
                synthetic,
                bins=bins,
                density=True,
                alpha=0.5,
                color=_CONTROL_COLOR,
                label="synthetic control pool (planted structure)",
            )
    else:
        _note_empty(ax_delta, "no differential-expression rows")
    ax_delta.set_xlabel(_LFC_UNITS)
    ax_delta.set_ylabel("density (genes per unit log2 fold change)")
    ax_delta.set_title(
        "(d) fold changes: real screen against the synthetic control pool, shared bins",
        fontsize=9,
    )
    _legend(ax_delta, fontsize=8)

    if has_untestable:
        ax_untestable = fig.add_subplot(grid[1, 2])
        untestable = _finite(pool, "frac_untestable")
        if untestable.size:
            upper = max(1.0, float(np.max(untestable)))
            ax_untestable.hist(
                untestable, bins=_even_bins(0.0, upper, 30), color=_REAL_COLOR, alpha=0.85
            )
        else:
            _note_empty(ax_untestable, "no finite values")
        ax_untestable.set_xlabel("fraction of genes DESeq2 could not test (zero base mean)")
        ax_untestable.set_ylabel("conditions (count)")
        ax_untestable.set_title("(e) untestable-gene fraction per condition", fontsize=9)

    fig.suptitle("build: the replicate pool the reliability is measured on", fontsize=11)
    return _finish(fig, out)


# --------------------------------------------------------------------------------------------
# split -- one condition becomes two half-profiles
# --------------------------------------------------------------------------------------------


def fig_split(pool: pd.DataFrame, per_pair: pd.DataFrame, out: Path) -> Path:
    """Whether the two halves each condition is split into are comparable, and large enough.

    A split-half correlation assumes two halves of equal precision. With an odd number of
    replicate plates that assumption is false by construction: one half gets one plate and the
    other gets two, and the half built from a single plate is the noisier of the two. Panel (a)
    draws the two group sizes side by side so the imbalance is read off the figure rather than
    argued about, and it is why the score step reports a corrected mean over the even-plate-count
    conditions alongside the corrected mean over all of them.

    Panel (b) is the other precondition: how many genes are finite in BOTH halves, which is how
    many points each condition's correlation is actually computed from. The line at 50 genes is
    the threshold below which a condition is not scored at all.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4), layout="constrained")
    ax_sizes, ax_genes = axes[0], axes[1]

    half0 = _finite(pool, "n_plates_half0")
    half1 = _finite(pool, "n_plates_half1")
    if half0.size or half1.size:
        ax_sizes.hist(
            [half0, half1],
            bins=_integer_bins([half0, half1]),
            label=["group 1 (responders selected here)", "group 2 (held out)"],
            color=[_REAL_COLOR, _CONTROL_COLOR],
        )
    else:
        _note_empty(ax_sizes, "no conditions in the pool")
    ax_sizes.set_xlabel("replicate plates in the group (count)")
    ax_sizes.set_ylabel("conditions (count)")
    ax_sizes.set_title("(a) group sizes: the one-plate-against-two imbalance", fontsize=9)
    _legend(ax_sizes, fontsize=8)

    scored = _finite(per_pair, "n_genes_scored")
    if scored.size:
        ax_genes.hist(scored, bins=_shared_bins([scored], n_bins=40), color=_REAL_COLOR, alpha=0.85)
    else:
        _note_empty(ax_genes, "no scored conditions")
    ax_genes.axvline(
        float(SCORING_THRESHOLD),
        color="k",
        lw=1.5,
        linestyle="--",
        label=f"scoring threshold ({SCORING_THRESHOLD} genes)",
    )
    ax_genes.set_xlabel("genes finite in both halves (count)")
    ax_genes.set_ylabel("conditions (count)")
    ax_genes.set_title("(b) genes each condition's correlation is computed over", fontsize=9)
    _legend(ax_genes, fontsize=8)

    fig.suptitle("split: one condition becomes two half-profiles", fontsize=11)
    return _finish(fig, out)


# --------------------------------------------------------------------------------------------
# select -- the responder set, chosen from the first half alone
# --------------------------------------------------------------------------------------------


def fig_select(
    per_pair: pd.DataFrame,
    padj_sample: pd.DataFrame,
    leakage: pd.DataFrame,
    out: Path,
    overlap: pd.DataFrame | None = None,
) -> Path:
    """Why the responder genes are chosen from one half only, and what that choice is worth.

    The responder reliability asks how reproducible a drug's response is over the genes that
    responded -- but "responded" has to be decided from the FIRST half alone. Panel (a) is the
    adjusted p-values that decision is made on; panel (b) is how many genes clear it per
    condition, against the same 50-gene floor, since a responder set below it is too small to
    correlate over.

    Panel (c) is the reason for the one-sided rule, measured rather than asserted. Both bars come
    from the same signal-free pool, where the true reliability is zero. The one-sided bar selects
    on the first half and scores on the second, and sits at its null. The pooled bar selects on
    both halves at once and then scores those same halves; it is inflated, and that inflation is
    the winner's curse -- genes chosen partly for noise the second half happens to share. The gap
    between the bars is exactly the bias the one-sided rule avoids, which is why nothing
    downstream is permitted to select on both halves.

    Panel (d), when an ``overlap`` table is supplied, is how far the two halves agree on which
    genes responded -- the Jaccard overlap of their responder sets per condition. It is a
    DIAGNOSTIC and never an input to selection: keeping the genes both halves called IS the
    pooled rule panel (c) measures the bias of. It is drawn because a reader will reasonably ask
    how noisy a one-sided responder call is, and the honest answer is a number rather than a
    reassurance.
    """
    n_panels = 3 if overlap is None else 4
    fig, axes = plt.subplots(1, n_panels, figsize=(5.0 * n_panels, 4.6), layout="constrained")
    axes = np.atleast_1d(axes)
    ax_padj, ax_responders, ax_leakage = axes[0], axes[1], axes[2]

    padj = _finite(padj_sample, "padj0")
    if padj.size:
        ax_padj.hist(padj, bins=_even_bins(0.0, 1.0, 40), color=_REAL_COLOR, alpha=0.85)
    else:
        _note_empty(ax_padj, "no adjusted p-values sampled")
    ax_padj.set_xlabel("adjusted p-value in group 1 (Benjamini-Hochberg)")
    ax_padj.set_ylabel("genes (count)")
    ax_padj.set_title("(a) the statistic responders are selected on", fontsize=9)

    responders = _finite(per_pair, "n_responders")
    if responders.size:
        ax_responders.hist(
            responders, bins=_shared_bins([responders], n_bins=40), color=_REAL_COLOR, alpha=0.85
        )
    else:
        _note_empty(ax_responders, "no scored conditions")
    ax_responders.axvline(
        float(SCORING_THRESHOLD),
        color="k",
        lw=1.5,
        linestyle="--",
        label=f"threshold ({SCORING_THRESHOLD} genes)",
    )
    ax_responders.set_xlabel("responder genes selected in group 1 (count)")
    ax_responders.set_ylabel("conditions (count)")
    ax_responders.set_title("(b) responder set size per condition", fontsize=9)
    _legend(ax_responders, fontsize=8)

    rules = ("one-sided", "pooled")
    lookup = dict(
        zip(_labels(leakage, "rule").tolist(), _numeric(leakage, "mean_r").tolist(), strict=False)
    )
    heights = np.array([float(lookup.get(rule, np.nan)) for rule in rules], dtype=float)
    positions = np.arange(len(rules), dtype=float)
    drawable = np.isfinite(heights)
    if drawable.any():
        colors = [_REAL_COLOR, _CONTROL_COLOR]
        ax_leakage.bar(
            positions[drawable],
            heights[drawable],
            color=[colors[i] for i in np.flatnonzero(drawable).tolist()],
            alpha=0.85,
        )
    else:
        _note_empty(ax_leakage, "no leakage measurement")
    ax_leakage.axhline(0.0, color="k", lw=1.0)
    ax_leakage.set_xticks(positions)
    ax_leakage.set_xticklabels(
        ["one-sided\n(select on group 1)", "pooled\n(select on both halves)"], fontsize=8
    )
    ax_leakage.set_xlabel("selection rule, both on the same signal-free pool")
    ax_leakage.set_ylabel("mean responder split-half Pearson r")
    ax_leakage.set_title(
        "(c) leakage: the pooled bar is the winner's curse\nthe one-sided rule avoids",
        fontsize=9,
    )

    fig.suptitle("select: responders chosen from the first half alone", fontsize=11)
    return _finish(fig, out)


# --------------------------------------------------------------------------------------------
# score -- the two reliabilities
# --------------------------------------------------------------------------------------------

_SUMMARY_MARKS = (
    ("all_splithalf_mean_r", _REAL_COLOR, "-", "all genes, raw split-half mean r"),
    ("all_spearman_brown_full", _REAL_COLOR, "--", "all genes, Spearman-Brown corrected"),
    (
        "all_spearman_brown_full_even_plates",
        _REAL_COLOR,
        ":",
        "all genes, corrected, even-plate conditions",
    ),
    ("responder_splithalf_mean_r", _RESPONDER_COLOR, "-", "responders, raw split-half mean r"),
    (
        "responder_spearman_brown_full",
        _RESPONDER_COLOR,
        "--",
        "responders, Spearman-Brown corrected",
    ),
    (
        "responder_spearman_brown_full_even_plates",
        _RESPONDER_COLOR,
        ":",
        "responders, corrected, even-plate conditions",
    ),
)


def _mark_summary_lines(ax: Axes, summary: Mapping[str, object]) -> None:
    """Draw the six reported reliabilities as labelled vertical lines on the r histograms.

    Raw split-half mean, its Spearman-Brown correction to full length, and that correction over
    the even-plate-count conditions only -- for each of the two gene sets. The even-plate value
    sits beside the all-conditions value so a reader can see for themselves whether the
    equal-halves assumption behind the correction moved the answer.
    """
    for key, color, style, label in _SUMMARY_MARKS:
        value = _number(summary.get(key))
        if np.isfinite(value):
            ax.axvline(value, color=color, lw=1.4, linestyle=style, label=f"{label} = {value:.3f}")


def _example_ids(profiles: pd.DataFrame, profile_index: pd.DataFrame, limit: int = 4) -> list[str]:
    """Up to ``limit`` example conditions, in the order the index lists them."""
    available = _ordered_unique(_labels(profiles, "example_id")) if len(profiles) else []
    if not available:
        return []
    wanted = _ordered_unique(_labels(profile_index, "example_id")) if len(profile_index) else []
    ordered = [name for name in wanted if name in set(available)]
    ordered += [name for name in available if name not in set(ordered)]
    return ordered[:limit]


def fig_score(
    profiles: pd.DataFrame,
    profile_index: pd.DataFrame,
    per_pair: pd.DataFrame,
    control_per_pair: pd.DataFrame | None,
    summary: Mapping[str, object],
    out: Path,
) -> Path:
    """The reliabilities themselves: one condition at a time, then the whole screen.

    The top row is what a split-half correlation is, drawn: one condition's genes, the first
    half's fold change against the second half's, one point per gene. A tight diagonal cloud is a
    reproducible response; a round cloud is a condition whose measured response is mostly noise.
    Each panel prints its own correlation, computed from exactly the points drawn, and those
    points and that number are written to a companion CSV beside the image (``<name>.values.csv``)
    so anyone who doubts the printed value can recompute it.

    The bottom row is every condition at once: the distribution of the all-gene correlation and
    of the responder correlation. The design expects the responder distribution to sit to the
    right -- most genes do not respond to most drugs, so an all-gene correlation is diluted by
    genes that carry no signal to reproduce. The reported means are marked as vertical lines, raw
    and Spearman-Brown corrected, with the even-plate-count subset's corrected value beside them.

    Where a control pool is supplied, its correlations are drawn beneath on the same x-axis. A
    signal-free control sits at zero; a pool planted at a known reliability sits at the planted
    value. Either way the real distribution is read against something rather than against
    nothing.
    """
    example_ids = _example_ids(profiles, profile_index)
    n_cols = max(1, len(example_ids))
    control_r = _finite(control_per_pair, "r") if control_per_pair is not None else None
    n_rows = 2 if control_r is None else 3

    # The bottom row carries six labelled reference lines in its legend, so the figure stays
    # wide enough for them even when only one example scatter sits above.
    fig = plt.figure(figsize=(max(11.0, 4.0 * n_cols), 3.8 * n_rows), layout="constrained")
    grid = fig.add_gridspec(n_rows, n_cols)

    profile_ids = _labels(profiles, "example_id")
    profile_genes = _labels(profiles, "gene")
    profile_lfc0 = _numeric(profiles, "lfc0")
    profile_lfc1 = _numeric(profiles, "lfc1")
    kinds = dict(
        zip(
            _labels(profile_index, "example_id").tolist(),
            _labels(profile_index, "kind").tolist(),
            strict=False,
        )
    )

    exported: list[pd.DataFrame] = []
    for position, example_id in enumerate(example_ids):
        ax = fig.add_subplot(grid[0, position])
        selected = (profile_ids == example_id) & np.isfinite(profile_lfc0)
        selected &= np.isfinite(profile_lfc1)
        x, y = profile_lfc0[selected], profile_lfc1[selected]
        r_printed = _pearson(x, y)
        if x.size:
            ax.scatter(x, y, s=6, alpha=0.4, color=_REAL_COLOR, edgecolors="none")
            low = float(np.minimum(x.min(), y.min()))
            high = float(np.maximum(x.max(), y.max()))
            ax.plot([low, high], [low, high], color="k", lw=0.8, linestyle="--", label="identity")
        else:
            _note_empty(ax, "no genes finite in both halves")
        kind = kinds.get(example_id, "")
        ax.set_title(f"{example_id}{f' ({kind})' if kind else ''}", fontsize=9)
        ax.set_xlabel(f"group 1 {_LFC_UNITS}")
        ax.set_ylabel(f"group 2 {_LFC_UNITS}")
        ax.text(
            0.04,
            0.95,
            f"r = {r_printed:.4f}\nn = {x.size} genes",
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8,
        )
        _legend(ax, fontsize=7)
        exported.append(
            pd.DataFrame(
                {
                    "example_id": np.repeat(example_id, x.size),
                    "gene": profile_genes[selected],
                    "lfc0": x,
                    "lfc1": y,
                    "r_printed": np.repeat(r_printed, x.size),
                }
            )
        )
    if not example_ids:
        _note_empty(fig.add_subplot(grid[0, 0]), "no example profiles exported")

    ax_hist = fig.add_subplot(grid[1, :])
    r_all = _finite(per_pair, "r")
    r_responder = _finite(per_pair, "r_responder")
    bins = _shared_bins(
        [r_all, r_responder, control_r if control_r is not None else np.empty(0)], n_bins=40
    )
    if r_all.size:
        ax_hist.hist(
            r_all, bins=bins, density=True, alpha=0.55, color=_REAL_COLOR, label="all genes"
        )
    if r_responder.size:
        ax_hist.hist(
            r_responder,
            bins=bins,
            density=True,
            alpha=0.45,
            color=_RESPONDER_COLOR,
            label="responder genes (selected in group 1)",
        )
    if not r_all.size and not r_responder.size:
        _note_empty(ax_hist, "no scored conditions")
    _mark_summary_lines(ax_hist, summary)
    ax_hist.set_xlabel("split-half Pearson r between the two half-profiles, per condition")
    ax_hist.set_ylabel("density (conditions per unit r)")
    ax_hist.set_title("real screen: per-condition reliability, both gene sets", fontsize=9)
    _legend(ax_hist, fontsize=7, ncol=2)

    if control_r is not None:
        ax_control = fig.add_subplot(grid[2, :], sharex=ax_hist)
        # The control pools are drawn SEPARATELY. Merging them makes one bimodal blur under a
        # title that is wrong for half of it: the positive pool should sit at its planted value
        # and the negative pool at zero, and the whole point of the panel is that a reader can
        # see those two places are different. Split on the `control` column when the run
        # labelled them, and say so plainly when it did not.
        pools: list[tuple[str, np.ndarray]] = []
        if control_per_pair is not None and "control" in control_per_pair.columns:
            labels = control_per_pair["control"].astype(str).to_numpy()
            all_r = _numeric(control_per_pair, "r")
            for name in dict.fromkeys(labels):
                vals = all_r[(labels == name) & np.isfinite(all_r)]
                if vals.size:
                    pools.append((str(name), vals))
        elif control_r.size:
            pools.append(("control pool (unlabelled)", control_r))
        if pools:
            for i, (name, vals) in enumerate(pools):
                ax_control.hist(
                    vals,
                    bins=bins,
                    density=True,
                    alpha=0.55,
                    color=_CONTROL_COLORS[i % len(_CONTROL_COLORS)],
                    label=f"{name} (n={vals.size})",
                )
                ax_control.axvline(
                    float(np.mean(vals)),
                    color=_CONTROL_COLORS[i % len(_CONTROL_COLORS)],
                    lw=1.4,
                    linestyle="--",
                )
        else:
            _note_empty(ax_control, "control pool has no scored conditions")
        ax_control.set_xlabel("split-half Pearson r per condition (same axis as the panel above)")
        ax_control.set_ylabel("density (conditions per unit r)")
        ax_control.set_title(
            "control pools: a planted reliability and a pool with none, on the axis above",
            fontsize=9,
        )
        _legend(ax_control, fontsize=7)

    fig.suptitle("score: split-half reliability, per condition and over the screen", fontsize=11)
    written = _finish(fig, out)

    columns = ["example_id", "gene", "lfc0", "lfc1", "r_printed"]
    table = (
        pd.concat(exported, ignore_index=True)
        if exported
        else pd.DataFrame({column: [] for column in columns})
    )
    table.to_csv(written.parent / f"{written.stem}.values.csv", index=False)
    return written


# --------------------------------------------------------------------------------------------
# decompose -- what kind of noise the ceiling is made of
# --------------------------------------------------------------------------------------------


def _fraction_hist(ax: Axes, fraction: np.ndarray, color: str, label: str, title: str) -> None:
    """One between-plate-fraction histogram, on the fixed zero-to-one axis both panels share."""
    if fraction.size:
        ax.hist(fraction, bins=_even_bins(0.0, 1.0, 40), color=color, alpha=0.85, label=label)
    else:
        _note_empty(ax, "no decomposed gene-condition rows")
    ax.set_xlabel("between-plate share of the fold-change variance (0 = none, 1 = all)")
    ax.set_ylabel("gene-condition rows (count)")
    ax.set_title(title, fontsize=9)
    _legend(ax, fontsize=8)


def fig_decompose(
    noise: pd.DataFrame,
    control_noise: pd.DataFrame | None,
    out: Path,
) -> Path:
    """What the reliability ceiling is made of: plate effects, or cell sampling.

    A condition's two half-profiles disagree for two separable reasons. Each plate carries its
    own batch offset, so profiles built from different plates differ even with infinite cells --
    that is the between-plate variance. And each plate counts a finite number of cells, so its
    fold changes carry sampling error, which DESeq2 already reports as ``lfcSE``. Which of the
    two dominates decides what would actually raise the ceiling: more plates, or deeper plates.

    Panel (a) is the between-plate share of the total. Panel (b) is the same comparison one gene
    at a time, plate variance against mean squared standard error, with the identity line drawn:
    points above the line are rows whose disagreement is batch, points below are rows whose
    disagreement is counting. Reading which side of that line the mass sits on is the whole
    purpose of the panel, so both axes are log-scaled over the positive rows, where variances
    span orders of magnitude.

    Panels (c) and (d) ask whether that share is constant across the screen or depends on where
    in it you stand -- lowly expressed genes are sampling-limited, and strongly responding genes
    have more real signal for a batch offset to be measured against. A control pool with a
    planted split, when supplied, is drawn beside panel (a) on shared axes: a decomposition that
    cannot recover a known split is not evidence about the real one.
    """
    # Row-aligned for the stratified panels, filtered only for the histogram: panels (c) and (d)
    # pair each row's fraction with that same row's expression and response size, so dropping
    # rows column by column would silently mismatch them.
    fraction_by_row = _numeric(noise, "between_plate_fraction")
    fraction = fraction_by_row[np.isfinite(fraction_by_row)]
    control_fraction = (
        _finite(control_noise, "between_plate_fraction") if control_noise is not None else None
    )
    n_cols = 3 if control_fraction is not None else 2
    fig = plt.figure(figsize=(5.2 * n_cols, 8.6), layout="constrained")
    grid = fig.add_gridspec(2, n_cols)

    ax_fraction = fig.add_subplot(grid[0, 0])
    _fraction_hist(
        ax_fraction, fraction, _REAL_COLOR, "real screen", "(a) between-plate share, real screen"
    )

    if control_fraction is not None:
        ax_control = fig.add_subplot(grid[0, 1], sharex=ax_fraction, sharey=ax_fraction)
        _fraction_hist(
            ax_control,
            control_fraction,
            _CONTROL_COLOR,
            "control pool (planted split)",
            "(a, control) between-plate share, control pool",
        )
        ax_scatter = fig.add_subplot(grid[0, 2])
    else:
        ax_scatter = fig.add_subplot(grid[0, 1])

    plate_var = _numeric(noise, "sigma2_plate")
    sampling_var = _numeric(noise, "mean_se2")
    finite = np.isfinite(plate_var) & np.isfinite(sampling_var)
    positive = finite & (plate_var > 0.0) & (sampling_var > 0.0)
    drawn = positive if positive.any() else finite
    x, y = sampling_var[drawn], plate_var[drawn]
    if x.size:
        ax_scatter.scatter(x, y, s=6, alpha=0.35, color=_REAL_COLOR, edgecolors="none")
        low = float(np.minimum(x.min(), y.min()))
        high = float(np.maximum(x.max(), y.max()))
        ax_scatter.plot(
            [low, high],
            [low, high],
            color="k",
            lw=1.0,
            linestyle="--",
            label="identity: plate effects equal cell sampling",
        )
    else:
        _note_empty(ax_scatter, "no decomposed gene-condition rows")
    if positive.any():
        ax_scatter.set_xscale("log")
        ax_scatter.set_yscale("log")
        dropped = int(finite.sum() - positive.sum())
        if dropped:
            ax_scatter.text(
                0.02,
                0.02,
                f"{dropped} rows at or below zero are not on the log axes",
                transform=ax_scatter.transAxes,
                fontsize=7,
            )
    ax_scatter.set_xlabel(
        "within-plate sampling variance, mean(lfcSE^2)\n(squared log2 fold change)"
    )
    ax_scatter.set_ylabel("between-plate variance, sigma^2_plate\n(squared log2 fold change)")
    ax_scatter.set_title(
        "(b) above the line, plate effects dominate;\nbelow it, cell sampling does", fontsize=9
    )
    _legend(ax_scatter, fontsize=7)

    ax_expression = fig.add_subplot(grid[1, 0])
    labels, means, sems, counts = _quartile_profile(fraction_by_row, _numeric(noise, "base_mean"))
    _quartile_panel(
        ax_expression,
        labels,
        means,
        sems,
        counts,
        "expression quartile (DESeq2 base mean, normalised counts)",
        "(c) between-plate share by expression",
    )

    ax_response = fig.add_subplot(grid[1, 1])
    labels, means, sems, counts = _quartile_profile(
        fraction_by_row, np.abs(_numeric(noise, "mean_lfc"))
    )
    _quartile_panel(
        ax_response,
        labels,
        means,
        sems,
        counts,
        "response-size quartile (absolute mean log2 fold change)",
        "(d) between-plate share by response size",
    )

    fig.suptitle("decompose: plate effects against cell sampling", fontsize=11)
    return _finish(fig, out)


# --------------------------------------------------------------------------------------------
# null -- what the reliabilities are read against
# --------------------------------------------------------------------------------------------

_GENE_SET_COLUMN = {"all": "r", "responder": "r_responder"}
_STRATUM_STYLE = (
    ("any_pair", "tab:green", "any mismatched pair"),
    ("diff_drug", "tab:orange", "different drug and line (generic-structure floor)"),
    ("same_drug", "tab:purple", "same drug, different line (line-specificity floor)"),
)


def fig_null(per_pair: pd.DataFrame, null_draws: pd.DataFrame, out: Path) -> Path:
    """The significance claim in one picture: matched conditions against three mismatched nulls.

    A split-half correlation is evidence only if it beats what you get by correlating one
    condition's first half against a DIFFERENT condition's second half. Three such nulls are
    drawn, in increasing strictness. Any mismatched pair is the loosest. Different drug and
    different line is the generic-structure floor -- whatever any perturbation does to any line.
    Same drug, different line is the strictest, because two lines given the same drug already
    share that drug's generic response; clearing this one is what makes a reliability specific to
    the line rather than to the drug.

    One panel per gene set, with each stratum's mean marked. The matched distribution has to sit
    to the right of all three, and visibly so, or the reliability is measuring shared structure
    rather than a reproducible line-specific response.
    """
    draw_gene_set = _labels(null_draws, "gene_set")
    draw_stratum = _labels(null_draws, "stratum")
    draw_r = _numeric(null_draws, "r")
    available = set(_ordered_unique(draw_gene_set))
    ordered = [name for name in ("all", "responder") if name in available]
    ordered += sorted(available - set(ordered) - {""})
    gene_sets = ordered or ["all"]

    fig, axes = plt.subplots(
        1, len(gene_sets), figsize=(6.4 * len(gene_sets), 4.8), layout="constrained", squeeze=False
    )
    for position, gene_set in enumerate(gene_sets):
        ax = axes[0][position]
        matched = _finite(per_pair, _GENE_SET_COLUMN.get(gene_set, "r"))
        in_gene_set = draw_gene_set == gene_set
        strata = {
            name: draw_r[in_gene_set & (draw_stratum == name) & np.isfinite(draw_r)]
            for name, _, _ in _STRATUM_STYLE
        }
        bins = _shared_bins([matched, *strata.values()], n_bins=45)
        if matched.size:
            ax.hist(
                matched,
                bins=bins,
                density=True,
                alpha=0.6,
                color=_REAL_COLOR,
                label="matched condition (real)",
            )
            ax.axvline(
                float(np.mean(matched)),
                color=_REAL_COLOR,
                lw=1.6,
                label=f"matched mean = {float(np.mean(matched)):.3f}",
            )
        for name, color, label in _STRATUM_STYLE:
            draws = strata[name]
            if not draws.size:
                continue
            ax.hist(draws, bins=bins, density=True, alpha=0.35, color=color, label=label)
            ax.axvline(
                float(np.mean(draws)),
                color=color,
                lw=1.4,
                linestyle="--",
                label=f"{name} null mean = {float(np.mean(draws)):.3f}",
            )
        if not matched.size and not any(draws.size for draws in strata.values()):
            _note_empty(ax, "no matched conditions and no null draws")
        ax.set_xlabel("split-half Pearson r between two half-profiles")
        ax.set_ylabel("density (pairs per unit r)")
        ax.set_title(f"gene set: {gene_set}", fontsize=10)
        _legend(ax, fontsize=7)

    fig.suptitle("null: matched conditions against the mismatched-pair floors", fontsize=11)
    return _finish(fig, out)


def fig_permutation_vs_bootstrap(
    perm_means: pd.DataFrame,
    null_draws: pd.DataFrame,
    summary: Mapping[str, object],
    out: Path,
) -> Path:
    """The dependence the bootstrap ignores, drawn instead of asserted.

    The p-values are read off a bootstrap that resamples null draws as if each were independent.
    They are not: mismatched draws reuse the same half-profiles many times over, so each draw
    shares a whole half-profile with many of its neighbours. That dependence makes the true
    sampling distribution of the null mean WIDER than the independent bootstrap believes, and a
    null that is too narrow makes every p-value look better than it is.

    So the permutation check permutes the pairing itself, which preserves the reuse structure,
    and its distribution of the mean is drawn here against the normal curve the independent
    bootstrap assumes -- same centre, standard deviation of the draws over the square root of the
    condition count. The two widths are meant to be compared directly: their ratio is the design
    effect, and it is the factor by which the bootstrap's error bar is too small.

    The bootstrap curve is built from the different-drug stratum of the all-gene draws, which is
    the stratum the reported p-value is read against.
    """
    perm = _finite(perm_means, "perm_mean")
    stratum = _labels(null_draws, "stratum")
    gene_set = _labels(null_draws, "gene_set")
    values = _numeric(null_draws, "r")
    selected = (stratum == "diff_drug") & np.isfinite(values)
    if (selected & (gene_set == "all")).any():
        selected &= gene_set == "all"
    draws = values[selected]
    n_pairs = _number(summary.get("all_n_pairs"))

    centre = float(np.mean(draws)) if draws.size else float("nan")
    draw_sd = float(np.std(draws, ddof=1)) if draws.size > 1 else float("nan")
    bootstrap_sd = float("nan")
    if np.isfinite(draw_sd) and np.isfinite(n_pairs) and n_pairs >= 1.0:
        bootstrap_sd = draw_sd / float(np.sqrt(n_pairs))

    fig, ax = plt.subplots(1, 1, figsize=(8.4, 5.2), layout="constrained")
    if perm.size:
        ax.hist(
            perm,
            bins=_shared_bins([perm], n_bins=40),
            density=True,
            alpha=0.6,
            color=_REAL_COLOR,
            label="permutation null: mean over a permuted pairing",
        )
    else:
        _note_empty(ax, "no permutation means")
    if np.isfinite(centre) and np.isfinite(bootstrap_sd) and bootstrap_sd > 0.0:
        span = 4.0 * bootstrap_sd
        low = min(centre - span, float(np.min(perm)) if perm.size else centre - span)
        high = max(centre + span, float(np.max(perm)) if perm.size else centre + span)
        curve_x = np.linspace(low, high, 400)
        curve_y = np.exp(-0.5 * ((curve_x - centre) / bootstrap_sd) ** 2) / (
            bootstrap_sd * np.sqrt(2.0 * np.pi)
        )
        ax.plot(
            curve_x,
            curve_y,
            color=_CONTROL_COLOR,
            lw=2.0,
            label="i.i.d. bootstrap: assumed sampling distribution of the mean",
        )
        ax.axvline(centre, color="k", lw=1.0, linestyle=":", label=f"centre = {centre:.4f}")

    perm_sd = float(np.std(perm, ddof=1)) if perm.size > 1 else float("nan")
    lines = [
        f"permutation sd of the mean = {perm_sd:.4f}",
        f"i.i.d. bootstrap sd of the mean = {bootstrap_sd:.4f}",
    ]
    design_effect = _number(summary.get("design_effect"))
    lines.append(
        f"design effect = {design_effect:.3f}"
        if np.isfinite(design_effect)
        else "design effect not reported in the summary"
    )
    # Bottom left, clear of the legend above it and of both distributions' tails.
    ax.text(
        0.02, 0.02, "\n".join(lines), transform=ax.transAxes, ha="left", va="bottom", fontsize=8
    )
    ax.set_xlabel("mean split-half Pearson r over the conditions, one value per resample")
    ax.set_ylabel("density (resamples per unit mean r)")
    ax.set_title(
        "design effect: the gap between the two widths is the dependence\n"
        "the i.i.d. bootstrap ignores",
        fontsize=10,
    )
    _legend(ax, fontsize=8)
    return _finish(fig, out)


# --------------------------------------------------------------------------------------------
# cross-cutting: the empirical control, and power
# --------------------------------------------------------------------------------------------


def fig_terciles(terciles: pd.DataFrame, out: Path) -> Path:
    """The empirical in-run control: reproducibility has to rise where there is more signal.

    Conditions are sorted by how large a response the drug produced and cut into thirds. If the
    assay is working, the third with the largest responses reproduces best -- a strong, real
    transcriptional change is easier to measure twice than a weak one. A flat or falling line
    here would say the correlations are driven by something other than drug response, and would
    invalidate both reliabilities regardless of how they compare to their nulls.

    The confidence intervals are what say whether the rise is real or within noise; a rise
    smaller than the intervals is not yet evidence of anything.
    """
    fig, ax = plt.subplots(1, 1, figsize=(6.6, 4.6), layout="constrained")
    tercile = _numeric(terciles, "tercile")
    mean_r = _numeric(terciles, "mean_r")
    ci_lo = _numeric(terciles, "ci_lo")
    ci_hi = _numeric(terciles, "ci_hi")
    counts = _numeric(terciles, "n")
    keep = np.isfinite(tercile) & np.isfinite(mean_r)
    if keep.any():
        x, y = tercile[keep], mean_r[keep]
        lower = np.nan_to_num(y - ci_lo[keep], nan=0.0)
        upper = np.nan_to_num(ci_hi[keep] - y, nan=0.0)
        yerr = np.vstack([np.clip(lower, 0.0, None), np.clip(upper, 0.0, None)])
        ax.errorbar(x, y, yerr=yerr, fmt="o-", color=_REAL_COLOR, capsize=4, lw=1.4, ms=6)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [
                f"{int(value)}\nn={int(count)}" if np.isfinite(count) else f"{int(value)}"
                for value, count in zip(x.tolist(), counts[keep].tolist(), strict=False)
            ],
            fontsize=8,
        )
    else:
        _note_empty(ax, "no tercile summary")
    ax.set_xlabel("tercile of response size (1 = weakest response, 3 = strongest)")
    ax.set_ylabel("mean split-half Pearson r (95% confidence interval)")
    ax.set_title("empirical control: reproducibility must rise with effect size", fontsize=10)
    return _finish(fig, out)


def fig_power(mde_curve: pd.DataFrame, out: Path) -> Path:
    """How much of the answer is the screen's size rather than the effect's.

    The minimum detectable effect is the smallest true mean correlation this design would call
    significant at alpha = 0.05 with 80% power, read off the same null bootstrap the p-values
    come from. Plotted against the number of conditions, it falls as the screen grows. Each gene
    set's own observed condition count is marked, so a null result can be told apart from an
    underpowered one: if the observed count sits where the curve is still high, the honest
    statement is that the screen cannot yet resolve an effect of that size.
    """
    fig, ax = plt.subplots(1, 1, figsize=(7.2, 4.8), layout="constrained")
    palette = [_REAL_COLOR, _RESPONDER_COLOR, "tab:green", "tab:purple"]
    curve_gene_set = _labels(mde_curve, "gene_set")
    curve_n_pairs = _numeric(mde_curve, "n_pairs")
    curve_mde = _numeric(mde_curve, "mde")
    curve_observed = _flags(mde_curve, "observed")
    gene_sets = _ordered_unique(curve_gene_set) if len(mde_curve) else []

    for position, gene_set in enumerate(gene_sets):
        color = palette[position % len(palette)]
        rows = (curve_gene_set == gene_set) & np.isfinite(curve_n_pairs) & np.isfinite(curve_mde)
        if rows.any():
            n_pairs, mde = curve_n_pairs[rows], curve_mde[rows]
            order = np.argsort(n_pairs, kind="stable")
            ax.plot(n_pairs[order], mde[order], "o-", color=color, lw=1.4, ms=4, label=gene_set)
        marked = curve_n_pairs[
            (curve_gene_set == gene_set) & curve_observed & np.isfinite(curve_n_pairs)
        ]
        if marked.size:
            ax.axvline(
                float(marked[0]),
                color=color,
                lw=1.2,
                linestyle="--",
                label=f"{gene_set}: observed {int(marked[0])} conditions",
            )
    if not gene_sets:
        _note_empty(ax, "no minimum-detectable-effect curve")
    ax.set_xlabel("conditions in the comparison (count)")
    ax.set_ylabel("minimum detectable mean split-half r\n(alpha = 0.05, power = 0.80)")
    ax.set_title("power: minimum detectable effect against screen size", fontsize=10)
    _legend(ax, fontsize=8)
    return _finish(fig, out)
