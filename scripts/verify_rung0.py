"""Executable verification of rung 0's promoted claims (PROCESS §3).

Every check here recomputes a promoted or reported number from the committed artifacts
alone -- no cluster access, no trust in any write-up -- and compares it to the claim.
Trust comes from re-derivation the reader performs, not narrative the writer asserts:
a reviewer who runs this has re-derived the evidence, not read about it.

    uv run python scripts/verify_rung0.py

prints one row per check (claim, recomputed value, PASS/FAIL) and exits nonzero on any
failure. The notebook ``docs/tasks/rung0-replicate-ceiling/verify.ipynb`` recomputes the
same claims in self-contained cells (standard-library hashing, direct table reads --
nothing imported from this module, so the reviewer reads exactly what is computed) and
runs this script only as its final cross-check; ``tests/test_verify_rung0.py`` runs this
battery in continuous integration so the branch stays green independent of anyone
opening the notebook.

What is NOT checkable locally, stated rather than hidden: the gene and drug panel files
live on Alpine and are pinned by sha256 in the provenance record, so the declared panel
size (14,121) is a recorded input property here, not a recomputable one; and the 1,026
data shards themselves are on scratch, so shard integrity reduces locally to the
committed manifest's content hash matching the promoted record's ``data_commit``.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
TASK = "rung0-replicate-ceiling"
TRANCHE = "tahoe100m-pseudobulk-de.v1"


@dataclass
class Check:
    name: str
    claim: str
    computed: str
    ok: bool


def sha256_of(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def _close(claim: float, actual: float, decimals: int) -> bool:
    """True when ``claim`` is ``actual`` rounded to ``decimals`` (half-a-unit tolerance)."""
    return abs(claim - actual) <= 0.5 * 10**-decimals + 1e-12


def _numbers(line: str) -> list[float]:
    return [float(m.replace(",", "")) for m in re.findall(r"\d[\d,]*\.?\d*", line)]


def _table_row(text: str, label: str) -> list[float]:
    """Numbers from the value cell of the markdown table row whose label contains ``label``."""
    for line in text.splitlines():
        if line.startswith("|") and label in line:
            return _numbers(line.rsplit("|", 2)[-2])
    raise AssertionError(f"summary.md has no table row labeled {label!r}")


def _headline(repo: Path) -> pd.Series:
    frame = pd.read_csv(repo / "results" / TASK / "rung0_delta_reproducibility.csv")
    return frame.iloc[0]


def _permutation(repo: Path) -> pd.Series:
    frame = pd.read_csv(repo / "docs" / "tasks" / TASK / "rung0_permutation_summary.csv")
    return frame.iloc[0]


def check_promoted_hashes(repo: Path) -> list[Check]:
    """The provenance record's checksums, re-derived from the committed artifacts."""
    record = json.loads(
        (repo / "results" / TASK / "rung0_delta_reproducibility.provenance.json").read_text()
    )
    promoted = repo / str(record["result"])
    task_side = repo / "docs" / "tasks" / TASK / "rung0_delta_reproducibility.csv"
    result_sha = sha256_of(promoted)
    log_sha = sha256_of(repo / str(record["log"]))
    return [
        Check(
            "promoted CSV checksum",
            f"record.result_sha256 {record['result_sha256'][:12]}…",
            f"sha256 {result_sha[:12]}…",
            result_sha == record["result_sha256"],
        ),
        Check(
            "job log checksum",
            f"record.log_sha256 {record['log_sha256'][:12]}…",
            f"sha256 {log_sha[:12]}…",
            log_sha == record["log_sha256"],
        ),
        Check(
            "task-side copy byte-identical to promoted copy",
            "identical bytes",
            "identical" if task_side.read_bytes() == promoted.read_bytes() else "DIFFER",
            task_side.read_bytes() == promoted.read_bytes(),
        ),
    ]


def check_tranche_content_hash(repo: Path) -> list[Check]:
    """The data pin: sha256 of the committed shard manifest is the recorded data_commit."""
    record = json.loads(
        (repo / "results" / TASK / "rung0_delta_reproducibility.provenance.json").read_text()
    )
    tranche = json.loads((repo / "data" / "tranches" / f"{TRANCHE}.json").read_text())
    manifest = repo / "data" / "tranches" / f"{TRANCHE}.manifest.txt"
    recomputed = sha256_of(manifest)
    n_shards = len(manifest.read_text().splitlines())
    return [
        Check(
            "tranche content hash from committed manifest",
            f"tranche record {tranche['content_hash'][:12]}…",
            f"sha256(manifest) {recomputed[:12]}…",
            recomputed == tranche["content_hash"],
        ),
        Check(
            "promoted record pins the same data",
            f"environment.data_commit {record['environment']['data_commit'][:12]}…",
            f"tranche content_hash {tranche['content_hash'][:12]}…",
            record["environment"]["data_commit"] == tranche["content_hash"],
        ),
        Check("manifest shard count", "1,026 shards", f"{n_shards} lines", n_shards == 1026),
    ]


def check_headline_from_raw_values(repo: Path) -> list[Check]:
    """Every headline statistic recomputed from the 1,600 per-condition correlations.

    The strongest check in the battery: not one stored summary against another, but the
    reported row against the raw values it summarizes (`rung0_per_pair_r.csv`, exported by
    the same job that produced the row -- see verification.md's job 31955710 section).
    """
    row = _headline(repo)
    per_pair = pd.read_csv(repo / "docs" / "tasks" / TASK / "rung0_per_pair_r.csv")
    r = per_pair["r"].to_numpy(dtype=float)
    effect = per_pair["mean_abs_delta"].to_numpy(dtype=float)
    mean = float(np.mean(r))
    edges = np.quantile(effect, [1 / 3, 2 / 3])
    terciles = [
        round(float(np.mean(r[(effect > lo) & (effect <= hi)])), 3)
        for lo, hi in ((-np.inf, edges[0]), (edges[0], edges[1]), (edges[1], np.inf))
    ]
    stats = {
        "n_pairs": (float(len(per_pair)), float(row["n_pairs"])),
        "splithalf_mean_r": (round(mean, 3), float(row["splithalf_mean_r"])),
        "splithalf_median_r": (round(float(np.median(r)), 3), float(row["splithalf_median_r"])),
        "splithalf_q1_r": (round(float(np.quantile(r, 0.25)), 3), float(row["splithalf_q1_r"])),
        "splithalf_q3_r": (round(float(np.quantile(r, 0.75)), 3), float(row["splithalf_q3_r"])),
        "frac_pos": (round(float(np.mean(r > 0)), 3), float(row["frac_pos"])),
        "spearman_brown_full": (
            round(2 * mean / (1 + mean), 3),
            float(row["spearman_brown_full"]),
        ),
    }
    checks = [
        Check(
            f"{name} recomputes from the raw per-condition values",
            f"reported {reported}",
            f"from {len(per_pair)} raw values: {computed}",
            computed == reported,
        )
        for name, (computed, reported) in stats.items()
    ]
    checks.append(
        Check(
            "effect-size terciles recompute from the raw values",
            " / ".join(str(float(row[f"splithalf_mean_r_tercile{t}"])) for t in (1, 2, 3)),
            " / ".join(str(t) for t in terciles),
            terciles == [float(row[f"splithalf_mean_r_tercile{t}"]) for t in (1, 2, 3)],
        )
    )
    return checks


def check_null_draws(repo: Path) -> list[Check]:
    """The chance floors re-derived from the individual mismatched-pair correlations."""
    row = _headline(repo)
    draws = pd.read_csv(repo / "docs" / "tasks" / TASK / "rung0_null_draws.csv")
    checks: list[Check] = []
    for stratum, column in (
        ("diff_drug", "null_diff_drug_mean_r"),
        ("same_drug", "null_same_drug_mean_r"),
        ("any_pair", "null_any_pair_mean_r"),
    ):
        d = draws[draws["stratum"] == stratum]["r"].to_numpy(dtype=float)
        mean = round(float(np.mean(d)), 3)
        checks.append(
            Check(
                f"{stratum} floor recomputes from its draws",
                f"reported {row[column]} over {int(row['null_n_draws'])} draws",
                f"mean of {d.size} committed draws: {mean}",
                mean == float(row[column]) and d.size == int(row["null_n_draws"]),
            )
        )
    same = draws[draws["stratum"] == "same_drug"]["r"].to_numpy(dtype=float)
    diff = draws[draws["stratum"] == "diff_drug"]["r"].to_numpy(dtype=float)
    checks.append(
        Check(
            "same-drug floor sits above the different-drug floor",
            "shared drug response lifts the stricter floor",
            f"{np.mean(same):.4f} > {np.mean(diff):.4f}",
            float(np.mean(same)) > float(np.mean(diff)),
        )
    )
    return checks


def check_example_profiles(repo: Path) -> list[Check]:
    """The scatter data reproduces the correlations it is shown under."""
    task_dir = repo / "docs" / "tasks" / TASK
    profiles = pd.read_csv(task_dir / "rung0_example_pair_profiles.csv.gz")
    index = pd.read_csv(task_dir / "rung0_example_pair_index.csv")
    per_pair = pd.read_csv(task_dir / "rung0_per_pair_r.csv", keep_default_na=False)
    checks: list[Check] = []
    for _, ex in index.iterrows():
        points = profiles[profiles["example_id"] == ex["example_id"]]
        recomputed = round(
            float(
                np.corrcoef(
                    points["lfc0"].to_numpy(dtype=float), points["lfc1"].to_numpy(dtype=float)
                )[0, 1]
            ),
            4,
        )
        ok = (
            len(points) == int(ex["n_genes_shown"])
            and abs(recomputed - float(ex["r_shown"])) < 5e-3
            and float(ex["r_shown"]) == float(ex["r_full"])
        )
        checks.append(
            Check(
                f"example {ex['example_id']}: scatter reproduces its correlation",
                f"index r {ex['r_shown']} over {int(ex['n_genes_shown'])} genes",
                f"from the committed points: {recomputed}",
                ok,
            )
        )
    matched = index[index["kind"] == "matched"]
    agree = [
        float(
            per_pair[(per_pair["patient"] == ex["patient0"]) & (per_pair["drug"] == ex["drug0"])][
                "r"
            ].iloc[0]
        )
        == float(ex["r_full"])
        for _, ex in matched.iterrows()
    ]
    checks.append(
        Check(
            "matched examples agree with the per-condition table",
            "each example's r equals its row in rung0_per_pair_r.csv",
            f"{sum(agree)} of {len(agree)} agree",
            all(agree),
        )
    )
    checks.append(
        Check(
            "matched examples span the reliability range in order",
            "5th < 25th < 50th < 95th percentile",
            " < ".join(str(v) for v in matched.sort_values("example_id")["r_full"]),
            list(matched.sort_values("example_id")["r_full"]) == sorted(matched["r_full"]),
        )
    )
    return checks


def check_headline_consistency(repo: Path) -> list[Check]:
    """The promoted row's internal arithmetic: the lift, the quartiles, the terciles."""
    row = _headline(repo)
    mean = float(row["splithalf_mean_r"])
    sb = 2 * mean / (1 + mean)
    terciles = [float(row[f"splithalf_mean_r_tercile{t}"]) for t in (1, 2, 3)]
    return [
        Check(
            "Spearman-Brown lift equals 2r/(1+r) of the mean",
            f"spearman_brown_full {row['spearman_brown_full']}",
            f"2·{mean}/(1+{mean}) = {sb:.3f}",
            _close(float(row["spearman_brown_full"]), sb, 3),
        ),
        Check(
            "quartile ordering",
            "q1 <= median <= q3",
            f"{row['splithalf_q1_r']} <= {row['splithalf_median_r']} <= {row['splithalf_q3_r']}",
            float(row["splithalf_q1_r"])
            <= float(row["splithalf_median_r"])
            <= float(row["splithalf_q3_r"]),
        ),
        Check(
            "floor ordering",
            "diff-drug floor < same-drug floor < observed mean",
            f"{row['null_diff_drug_mean_r']} < {row['null_same_drug_mean_r']} < {mean}",
            float(row["null_diff_drug_mean_r"]) < float(row["null_same_drug_mean_r"]) < mean,
        ),
        Check(
            "tercile positive control is monotone",
            "tercile1 < tercile2 < tercile3",
            " < ".join(str(t) for t in terciles),
            terciles[0] < terciles[1] < terciles[2],
        ),
        Check(
            "both detection thresholds sit below the observed mean",
            "mde_80_vs_diff_drug and mde_80_vs_same_drug < mean",
            f"{row['mde_80_vs_diff_drug']}, {row['mde_80_vs_same_drug']} < {mean}",
            float(row["mde_80_vs_diff_drug"]) < mean and float(row["mde_80_vs_same_drug"]) < mean,
        ),
    ]


def check_pool_arithmetic(repo: Path) -> list[Check]:
    """The scored-pair count re-derived from the measured pool description."""
    # keep_default_na: one line key is literally "NA" (a missing DepMap id the grouping
    # carried through, documented in verification.md) and must count as a key, not as NaN
    pool = pd.read_csv(
        repo / "docs" / "tasks" / TASK / "rung0_pool_description.csv", keep_default_na=False
    )
    row = _headline(repo)
    unsplittable = pool[(pool["n_plates_half0"] == 0) | (pool["n_plates_half1"] == 0)]
    return [
        Check(
            "candidate conditions",
            "1,650 (line, drug) rows: 50 line keys x 33 drug names",
            f"{len(pool)} rows, {pool['patient'].nunique()} lines, {pool['drug'].nunique()} drugs",
            len(pool) == 1650 and pool["patient"].nunique() == 50 and pool["drug"].nunique() == 33,
        ),
        Check(
            "unscoreable conditions are exactly Ribociclib's",
            "50 rows with an empty half, all Ribociclib",
            f"{len(unsplittable)} rows, drugs {sorted(unsplittable['drug'].unique())}",
            len(unsplittable) == 50 and set(unsplittable["drug"]) == {"Ribociclib"},
        ),
        Check(
            "scored-pair arithmetic",
            "1,650 - 50 = n_pairs",
            f"{len(pool)} - {len(unsplittable)} = {len(pool) - len(unsplittable)}; "
            f"n_pairs {row['n_pairs']}",
            len(pool) - len(unsplittable) == int(row["n_pairs"]),
        ),
    ]


def check_permutation(repo: Path) -> list[Check]:
    """The permutation nulls, re-derived from the committed per-permutation means."""
    s = _permutation(repo)
    task_dir = repo / "docs" / "tasks" / TASK
    strata = {
        "any-pair": ("rung0_permutation_perm_means.csv", "observed_mean", ""),
        "same-drug": (
            "rung0_permutation_perm_means_same_drug.csv",
            "observed_mean_same_drug_rows",
            "_same_drug",
        ),
        "diff-drug": (
            "rung0_permutation_perm_means_diff_drug.csv",
            "observed_mean_diff_drug_rows",
            "_diff_drug",
        ),
    }
    checks: list[Check] = []
    for name, (filename, observed_col, suffix) in strata.items():
        perms = pd.read_csv(task_dir / filename)["perm_mean"].to_numpy(dtype=float)
        observed = float(s[observed_col])
        p_exact = (1 + int(np.sum(perms >= observed))) / (1 + perms.size)
        checks.extend(
            [
                Check(
                    f"{name}: permutation-mean summary",
                    f"mean {s['perm_mean_mean' + suffix]}, sd {s['perm_mean_sd' + suffix]} "
                    f"over {int(s['n_perm'])} permutations",
                    f"mean {np.mean(perms):.4f}, sd {np.std(perms, ddof=1):.4f} over {perms.size}",
                    perms.size == int(s["n_perm"])
                    and _close(float(s["perm_mean_mean" + suffix]), float(np.mean(perms)), 4)
                    and _close(float(s["perm_mean_sd" + suffix]), float(np.std(perms, ddof=1)), 4),
                ),
                Check(
                    f"{name}: exact permutation p",
                    f"p_exact{suffix} {s['p_exact' + suffix]}",
                    f"(1 + {int(np.sum(perms >= observed))}) / (1 + {perms.size}) = {p_exact:.3f}",
                    _close(float(s["p_exact" + suffix]), p_exact, 3),
                ),
                Check(
                    f"{name}: observed exceeds every permutation",
                    f"observed {observed} > max of {perms.size} permutation means",
                    f"max {np.max(perms):.4f}",
                    observed > float(np.max(perms)),
                ),
            ]
        )
    # The any-pair design effect and z-score are recomputable from committed values; the
    # per-stratum design effects need each stratum's pooled standard error, which only the
    # cluster run holds, so they are cross-checked in the summary-table check instead.
    perms = pd.read_csv(task_dir / "rung0_permutation_perm_means.csv")["perm_mean"].to_numpy(
        dtype=float
    )
    de = float(np.var(perms, ddof=1)) / float(s["se_iid_pool"]) ** 2
    z = (float(s["observed_mean"]) - float(np.mean(perms))) / float(np.std(perms, ddof=1))
    checks.append(
        Check(
            "any-pair: design effect and z from committed draws",
            f"design_effect {s['design_effect']}, z_permutation {s['z_permutation']}",
            f"var(perm)/se² = {de:.3f}, z = {z:.2f}",
            abs(de - float(s["design_effect"])) < 0.02 and abs(z - float(s["z_permutation"])) < 0.5,
        )
    )
    # Cross-mechanism floor agreement: permutation-null means vs the headline CSV's
    # bootstrapped floors -- two different sampling mechanisms landing on the same floors.
    row = _headline(repo)
    checks.append(
        Check(
            "stratum floors agree across mechanisms",
            f"bootstrap floors {row['null_same_drug_mean_r']} / {row['null_diff_drug_mean_r']}",
            f"permutation means {s['perm_mean_mean_same_drug']} / {s['perm_mean_mean_diff_drug']}",
            abs(float(s["perm_mean_mean_same_drug"]) - float(row["null_same_drug_mean_r"])) < 0.0015
            and abs(float(s["perm_mean_mean_diff_drug"]) - float(row["null_diff_drug_mean_r"]))
            < 0.0015,
        )
    )
    return checks


def check_per_gene_diagnostic(repo: Path) -> list[Check]:
    """verification.md's per-gene numbers, re-derived from the committed diagnostic CSV."""
    pg = pd.read_csv(repo / "docs" / "tasks" / TASK / "rung0_per_gene_reliability.csv")
    finite = pg[np.isfinite(pg["r"])]
    top5 = pg.nlargest(5, "r")
    stated_top = {"HSP90AA1": 0.79, "EGR1": 0.75, "HSPA1B": 0.72, "HSPH1": 0.71, "PLEC": 0.69}
    top_ok = list(top5["gene"]) == list(stated_top) and all(
        _close(stated_top[str(g)], float(r), 2)
        for g, r in zip(top5["gene"], top5["r"], strict=True)
    )
    return [
        Check(
            "per-gene coverage",
            "13,886 genes, 13,759 finite (127 with too few pairs)",
            f"{len(pg)} genes, {len(finite)} finite ({len(pg) - len(finite)} not)",
            len(pg) == 13886 and len(finite) == 13759,
        ),
        Check(
            "per-gene distribution",
            "97.0% positive; median 0.146, quartiles 0.089-0.230",
            f"{100 * float((finite['r'] > 0).mean()):.1f}% positive; "
            f"median {finite['r'].median():.3f}, "
            f"quartiles {finite['r'].quantile(0.25):.3f}-{finite['r'].quantile(0.75):.3f}",
            _close(0.970, float((finite["r"] > 0).mean()), 3)
            and _close(0.146, float(finite["r"].median()), 3)
            and _close(0.089, float(finite["r"].quantile(0.25)), 3)
            and _close(0.230, float(finite["r"].quantile(0.75)), 3),
        ),
        Check(
            "most reproducible genes are the stated stress-response set",
            "HSP90AA1 0.79, EGR1 0.75, HSPA1B 0.72, HSPH1 0.71, PLEC 0.69",
            ", ".join(f"{g} {r:.2f}" for g, r in zip(top5["gene"], top5["r"], strict=True)),
            top_ok,
        ),
    ]


def check_summary_table(repo: Path) -> list[Check]:
    """Every number in summary.md's evidence, matched to the artifact it comes from.

    This is the anti-drift check: a transcribed number that no longer matches its artifact
    fails here mechanically, instead of waiting for an audit pass to read for it.
    """
    text = (repo / "docs" / "tasks" / TASK / "summary.md").read_text()
    flat = re.sub(r"\s+", " ", text)
    row = _headline(repo)
    s = _permutation(repo)

    scored = _table_row(text, "Conditions scored")
    present = _table_row(text, "Panel genes present")
    reliability = _table_row(text, "Split-half reliability")
    ceiling = _table_row(text, "Spearman-Brown")
    positive = _table_row(text, "positive reliability")
    floor_diff = _table_row(text, "Mismatched-condition floor")
    floor_same = _table_row(text, "Same-drug floor")
    signif = _table_row(text, "Significance vs both floors")
    mde = _table_row(text, "Smallest detectable effect")[-2:]
    terciles = _numbers(re.search(r"rises with effect size: ([\d./ ]+) across", flat).group(1))  # type: ignore[union-attr]
    effects = re.search(
        r"\((\d\.\d+) for the original all-conditions check, (\d\.\d+) for the "
        r"across-drugs-and-lines version, and (\d\.\d+) for the within-drug version\)",
        flat,
    )
    assert effects is not None, "summary.md's per-stratum design-effect parenthetical not found"
    p_shuffle = re.search(r"p = (0\.\d+), the strongest claim", flat)
    assert p_shuffle is not None, "summary.md's shuffle-check p-value sentence not found"
    ratios = re.search(
        r"roughly (\d\.\d+) times the mismatched-condition detection threshold "
        r"\((0\.\d+)\) and roughly (\d\.\d+) times the stricter same-drug",
        flat,
    )
    assert ratios is not None, "summary.md's power-ratio sentence not found"

    return [
        Check(
            "summary: conditions scored",
            f"{scored[0]:.0f}",
            f"n_pairs {row['n_pairs']}",
            scored[0] == float(row["n_pairs"]),
        ),
        Check(
            "summary: panel genes present",
            f"{present[0]:.0f} (declared {present[1]:.0f} is pinned by the gene_panel "
            "hash, not recomputable locally)",
            f"n_genes {row['n_genes']}",
            present[0] == float(row["n_genes"]),
        ),
        Check(
            "summary: reliability mean, median, quartiles",
            f"{reliability[0]} ({reliability[1]}; {reliability[2]}-{reliability[3]})",
            f"{row['splithalf_mean_r']} ({row['splithalf_median_r']}; "
            f"{row['splithalf_q1_r']}-{row['splithalf_q3_r']})",
            reliability
            == [
                float(row["splithalf_mean_r"]),
                float(row["splithalf_median_r"]),
                float(row["splithalf_q1_r"]),
                float(row["splithalf_q3_r"]),
            ],
        ),
        Check(
            "summary: Spearman-Brown ceiling",
            f"{ceiling[0]}",
            f"{row['spearman_brown_full']}",
            ceiling[0] == float(row["spearman_brown_full"]),
        ),
        Check(
            "summary: conditions with positive reliability",
            f"{positive[0]}%",
            f"frac_pos {row['frac_pos']}",
            _close(positive[0] / 100, float(row["frac_pos"]), 3),
        ),
        Check(
            "summary: the two floors",
            f"{floor_diff[-1]} / {floor_same[-1]}",
            f"{row['null_diff_drug_mean_r']} / {row['null_same_drug_mean_r']}",
            floor_diff[-1] == float(row["null_diff_drug_mean_r"])
            and floor_same[-1] == float(row["null_same_drug_mean_r"]),
        ),
        Check(
            "summary: significance vs both floors",
            f"p = {signif[0]}",
            f"p_vs_null {row['p_vs_null']}, p_vs_same_drug {row['p_vs_same_drug']}",
            signif[0] == float(row["p_vs_null"]) and signif[0] == float(row["p_vs_same_drug"]),
        ),
        Check(
            "summary: detection thresholds",
            f"{mde[0]} / {mde[1]}",
            f"{row['mde_80_vs_diff_drug']} / {row['mde_80_vs_same_drug']}",
            _close(mde[0], float(row["mde_80_vs_diff_drug"]), 3)
            and _close(mde[1], float(row["mde_80_vs_same_drug"]), 3),
        ),
        Check(
            "summary: power ratios",
            f"{ratios.group(1)}x vs {ratios.group(2)}, {ratios.group(3)}x (stricter)",
            f"{float(row['splithalf_mean_r']) / float(row['mde_80_vs_diff_drug']):.1f}x, "
            f"{float(row['splithalf_mean_r']) / float(row['mde_80_vs_same_drug']):.1f}x",
            _close(
                float(ratios.group(1)),
                float(row["splithalf_mean_r"]) / float(row["mde_80_vs_diff_drug"]),
                1,
            )
            and _close(
                float(ratios.group(3)),
                float(row["splithalf_mean_r"]) / float(row["mde_80_vs_same_drug"]),
                1,
            ),
        ),
        Check(
            "summary: tercile control values",
            " / ".join(str(t) for t in terciles),
            f"{row['splithalf_mean_r_tercile1']} / {row['splithalf_mean_r_tercile2']} / "
            f"{row['splithalf_mean_r_tercile3']}",
            [float(t) for t in terciles]
            == [float(row[f"splithalf_mean_r_tercile{t}"]) for t in (1, 2, 3)],
        ),
        Check(
            "summary: per-stratum shuffle-check design effects",
            f"{effects.group(1)} / {effects.group(2)} / {effects.group(3)}",
            f"{s['design_effect']} / {s['design_effect_diff_drug']} / "
            f"{s['design_effect_same_drug']}",
            _close(float(effects.group(1)), float(s["design_effect"]), 2)
            and _close(float(effects.group(2)), float(s["design_effect_diff_drug"]), 2)
            and _close(float(effects.group(3)), float(s["design_effect_same_drug"]), 2),
        ),
        Check(
            "summary: shuffle-check p-value",
            f"p = {p_shuffle.group(1)}",
            f"p_exact {s['p_exact']} (recomputed independently above)",
            float(p_shuffle.group(1)) == float(s["p_exact"]),
        ),
    ]


def run_all_checks(repo: Path = REPO) -> list[Check]:
    return [
        *check_promoted_hashes(repo),
        *check_tranche_content_hash(repo),
        *check_headline_from_raw_values(repo),
        *check_null_draws(repo),
        *check_example_profiles(repo),
        *check_headline_consistency(repo),
        *check_pool_arithmetic(repo),
        *check_permutation(repo),
        *check_per_gene_diagnostic(repo),
        *check_summary_table(repo),
    ]


def render(checks: list[Check]) -> str:
    lines = []
    for c in checks:
        mark = "PASS" if c.ok else "FAIL"
        lines.append(f"[{mark}] {c.name}")
        lines.append(f"       claim:      {c.claim}")
        lines.append(f"       recomputed: {c.computed}")
    n_ok = sum(c.ok for c in checks)
    lines.append(f"\n{n_ok} / {len(checks)} checks pass")
    return "\n".join(lines)


def main() -> int:
    checks = run_all_checks()
    print(render(checks))
    return 0 if all(c.ok for c in checks) else 1


if __name__ == "__main__":
    sys.exit(main())
