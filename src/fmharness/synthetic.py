"""Synthetic replicate pools with planted, known answers.

Rung 0's figures put real data beside a planted answer on shared axes: a histogram of the real
screen shows what the assay looks like, and the same histogram of a pool whose reliability we
chose shows whether the machinery reads it correctly. This module is where that pool comes from,
so the run and its controls plant the same generative model rather than two similar ones.

The model, per (cell line, drug, gene):

    delta = signal + plate_offset + sampling_noise

``signal`` is fixed for the condition (variance ``signal_sd ** 2``), ``plate_offset`` is shared by
every gene on a plate, and ``sampling_noise`` is independent per plate and gene. A split-half
correlation over ``p`` plates per half therefore has expectation

    signal_sd ** 2 / (signal_sd ** 2 + noise_sd ** 2 / p)

which is the closed form the score controls check against. Inverting it: to plant a FULL-data
reliability ``R`` over ``n`` plates, set ``noise_sd = sqrt(n * (1 - R) / R)`` with unit signal --
then the half correlation is ``R / (2 - R)`` and Spearman-Brown lifts it back to ``R``.
"""

# pyright: reportMissingTypeStubs=false, reportUnknownMemberType=false, reportUnknownArgumentType=false, reportUnknownVariableType=false, reportUnknownParameterType=false, reportUnknownLambdaType=false

from __future__ import annotations

import numpy as np
import pandas as pd


def noise_sd_for_reliability(full_reliability: float, n_plates: int) -> float:
    """The per-plate noise that plants ``full_reliability`` over ``n_plates`` at unit signal."""
    if not 0.0 < full_reliability < 1.0:
        raise ValueError(f"full_reliability must lie strictly in (0, 1), got {full_reliability}")
    return float(np.sqrt(n_plates * (1.0 - full_reliability) / full_reliability))


def planted_split_half_frame(
    n_lines: int = 20,
    n_drugs: int = 10,
    n_genes: int = 2000,
    n_plates: int = 4,
    signal_sd: float = 1.0,
    noise_sd: float = 1.0,
    n_responders: int | None = None,
    seed: int = 0,
) -> pd.DataFrame:
    """A pool already aggregated to the two plate halves, in the built frame's own shape.

    Returns the columns ``build_split_half_frame`` returns -- ``patient``, ``drug``,
    ``gene_name``, ``lfc0``, ``lfc1``, ``padj0`` -- so a control pool flows through exactly the
    scoring, selection and null code the real data flows through. The plate aggregation is done
    here in closed form (the mean of ``n_plates / 2`` draws has noise variance
    ``noise_sd ** 2 / (n_plates / 2)``) rather than by writing plates and re-aggregating: this
    module feeds the figures, and the SQL that does the real aggregation is exercised separately
    by the build controls in ``tests/test_rung0_controls.py``.

    ``padj0`` is planted the way the first plate group's adjusted p-values behave: below 0.01 for
    the first ``n_responders`` genes and above 0.2 for the rest, or uniform on (0, 1) for every
    gene when ``n_responders`` is None, which is the signal-free case.
    """
    if n_plates < 2 or n_plates % 2:
        raise ValueError(f"n_plates must be an even number of at least 2, got {n_plates}")
    rng = np.random.default_rng(seed)
    per_half = n_plates / 2.0
    half_noise = noise_sd / np.sqrt(per_half)
    n_cond = n_lines * n_drugs

    signal = rng.normal(0.0, signal_sd, size=(n_cond, n_genes))
    lfc0 = signal + rng.normal(0.0, half_noise, size=(n_cond, n_genes))
    lfc1 = signal + rng.normal(0.0, half_noise, size=(n_cond, n_genes))
    if n_responders is None:
        padj0 = rng.uniform(0.0, 1.0, size=(n_cond, n_genes))
    else:
        is_resp = np.arange(n_genes) < n_responders
        padj0 = np.where(
            is_resp[None, :],
            rng.uniform(0.0, 0.01, size=(n_cond, n_genes)),
            rng.uniform(0.2, 1.0, size=(n_cond, n_genes)),
        )

    lines = np.repeat([f"SYNL{i}" for i in range(n_lines)], n_drugs * n_genes)
    drugs = np.tile(np.repeat([f"SYND{j}" for j in range(n_drugs)], n_genes), n_lines)
    genes = np.tile([f"SYNG{k}" for k in range(n_genes)], n_cond)
    return pd.DataFrame(
        {
            "patient": lines,
            "drug": drugs,
            "gene_name": genes,
            "lfc0": lfc0.ravel(),
            "lfc1": lfc1.ravel(),
            "padj0": padj0.ravel(),
        }
    )


def planted_noise_frame(
    n_cond: int = 40,
    n_genes: int = 2000,
    n_plates: int = 4,
    plate_sd: float = 0.5,
    within_sd: float = 0.5,
    seed: int = 0,
) -> pd.DataFrame:
    """A pool for the decomposition control, in ``build_noise_frame``'s own shape.

    Plants a plate component of variance ``plate_sd ** 2`` on top of sampling error of variance
    ``within_sd ** 2``, and reports the second as ``mean_se2`` -- which is what a correctly
    calibrated ``lfcSE`` would say. ``decompose_noise`` must then recover ``plate_sd ** 2``, and
    the between-plate fraction must come back at ``plate_sd**2 / (plate_sd**2 + within_sd**2)``.
    """
    rng = np.random.default_rng(seed)
    n = n_cond * n_genes
    # var_samp over n_plates draws of (plate offset + sampling error): its expectation is
    # plate_sd^2 + within_sd^2, and its own sampling spread is the chi-square factor below.
    total = plate_sd**2 + within_sd**2
    dof = n_plates - 1
    var_lfc = total * rng.chisquare(dof, size=n) / dof
    return pd.DataFrame(
        {
            "patient": np.repeat([f"SYNL{i}" for i in range(n_cond)], n_genes),
            "drug": "SYND0",
            "dose": 0.05,
            "gene_name": np.tile([f"SYNG{k}" for k in range(n_genes)], n_cond),
            "var_lfc": var_lfc,
            "mean_se2": np.full(n, within_sd**2),
            "n_plates": n_plates,
            "base_mean": rng.lognormal(3.0, 1.0, size=n),
            "mean_lfc": rng.normal(0.0, 0.5, size=n),
        }
    )
