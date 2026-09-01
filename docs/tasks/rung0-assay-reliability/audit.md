# Rung 0 — drift audit

**Date** 2026-09-01
**Commit audited** `c3e55b28069638c3020fcb87e7d4972fca73ffae` (`c3e55b2`, branch `rung0-assay-reliability`)
**Auditor** A fresh reader. I did not write this code or its design, and I worked from
`docs/audit.md`, `docs/tasks/rung0-assay-reliability/design.md`, `docs/SPEC.md`,
`docs/tasks/rung0-assay-reliability/decisions.md`, `docs/tasks/rung0-assay-reliability/plan.md`
and the landed tree — not from the conversation that produced them.
**Standard** [`docs/audit.md`](../../audit.md).

> **The cluster run has not finished.** No result artifact exists in the task folder
> (`ls docs/tasks/rung0-assay-reliability/` returns only `decisions.md`, `design.md`, `plan.md`,
> `summary.ipynb`, `verify.ipynb`). Every claim whose truth depends on a number the run has not
> produced is verdicted **PENDING RUN** in its own section below, never ALIGNED on faith and never
> DRIFT for being absent.

> **The tree moved during this audit.** I began reading at `8654f7c`; two commits (`8f09a7b`,
> `c3e55b2`) landed while I was reading. Every verdict below was re-checked against `c3e55b2`, and
> the two commits' diffs were read in full. `git status` shows `summary.ipynb` as modified, but
> `git diff -- docs/tasks/rung0-assay-reliability/summary.ipynb` is empty — an mtime-only change,
> not a content change.

---

## Method

### What was read

In full: `docs/audit.md`, `design.md`, `docs/SPEC.md` (rung 0 and the four project rules),
`decisions.md`, `plan.md`, `scripts/delta_reproducibility.py` (1,407 lines),
`scripts/permutation_null.py`, `src/fmharness/figures.py`, `src/fmharness/statistics.py`,
`src/fmharness/synthetic.py`, `scripts/alpine/delta_reproducibility.sbatch`,
`scripts/alpine/permutation_null.sbatch`, `docs/DATA.md`, `docs/STATE.md`,
`data/tranches/tahoe100m-pseudobulk-de.v1.json`, `.gitattributes`, `pyproject.toml`, and the
rule-4 section of `tests/test_project_rules.py`. Read by delegated readers and cross-checked
against my own greps: `scripts/verify_rung0.py`, `verify.ipynb`, `summary.ipynb`, and the five
rung-0 test files.

### What was recomputed, and the commands

**The whole suite passes.**

```
$ uv run pytest -q
..............ss.....s.................................................. [ 52%]
................................................................         [100%]
SKIPPED [1] tests/test_project_rules.py:73  no promoted results yet
SKIPPED [1] tests/test_project_rules.py:100 no provenance records yet
SKIPPED [1] tests/test_project_rules.py:355 no promoted results yet
```
133 passed, 3 skipped. The three skips are the promotion rules, which cannot bind before promotion.

**The tranche content hash recomputes from the committed manifest.**

```
$ python3 -c "import hashlib,pathlib; print(hashlib.sha256(pathlib.Path(
    'data/tranches/tahoe100m-pseudobulk-de.v1.manifest.txt').read_bytes()).hexdigest())"
9a8797a5698e2c56ec1b61bdd3d5f68d18a972e227e86b64ac341ef507f73dd6
$ wc -l data/tranches/tahoe100m-pseudobulk-de.v1.manifest.txt
    1026
```
Matches `content_hash` in `data/tranches/tahoe100m-pseudobulk-de.v1.json`, over exactly 1,026
shard lines, and matches the construction in `scripts/register_tranche.py:56-58`.

**Which declared figure functions have a caller.** Mechanical, by grep, not by reading:

```
$ grep -rn fig_permutation_vs_bootstrap scripts src tests
src/fmharness/figures.py:57      (in __all__)
src/fmharness/figures.py:976     (def)
tests/test_rung0_figures.py:439,440,446,451,457,458
```
No run script calls it. Same method found `write_figure` (`delta_reproducibility.py:864`) with no
caller at all, and `frac_untestable` produced by no code path.

**Notebook outputs and imports**, by script over the raw JSON:

```
summary.ipynb: 26 cells, 13 code cells, 0 outputs
   imports: pandas, IPython.display, pathlib
verify.ipynb : 38 cells, 19 code cells, 0 outputs
   imports: csv, gzip, hashlib, json, math, statistics, subprocess, sys, pathlib
```

**What the port actually changed.** The design claims select, decompose and the figures are the
only additions to the carried-over `delta_reproducibility.py`:

```
$ git diff 640a428..HEAD --stat -- scripts/delta_reproducibility.py
 scripts/delta_reproducibility.py | 970 ++++++++++++++++++++--- (900 insertions, 140 deletions)
$ diff <(git show 640a428:scripts/delta_reproducibility.py | grep '^def ') \
       <(grep '^def ' scripts/delta_reproducibility.py)
```
Sixteen functions added, of which `_compact_df`, `_drug_predicate`, `resolve_drug_names`,
`effect_size_tercile_table`, `mde_curve_table`, `spearman_brown_or_nan` and
`write_audit_checksums` belong to none of the three named additions. See D109.

**This audit's own arithmetic, recounted by script.** `docs/audit.md:86-89` records that a previous
pass was spent reconciling an audit's own counts, so the Counts section below was not tallied by
reading. A script parsed the verdict tables back out of this file:

```
verdict rows found: 138
D rows: 119 unique: 119     S rows: 19 unique: 19
duplicates: []              missing D: []      missing S: []
tally: {'ALIGNED': 106, 'DRIFT': 21, 'PENDING RUN': 10, 'DEVIATION-RECORDED': 1}   sum: 138
DRIFT ids: D6 D8 D22 D29 D30 D35 D39 D44 D56 D60 D63 D70 D74 D75 D78 D79 D83 D96 D97 D102 D109
PENDING ids: D17 D99 D100 D101 D104 D105 D106 D107 S11 S13
```
D1–D119 and S1–S19 each appear exactly once, with no gaps and no duplicates.

### Artifact checksums

`docs/audit.md` requires the checksum of every artifact the audit read. **No run artifact exists
to checksum.** In its place, the sha256 prefixes of every source file this audit verdicted
against, so a fix wave can tell whether a file moved under a verdict:

| sha256 (first 16) | File |
|---|---|
| `378d8c583a9d8aec` | `scripts/delta_reproducibility.py` |
| `27780a005e42724d` | `scripts/permutation_null.py` |
| `bc75b280ca31d7c9` | `scripts/verify_rung0.py` |
| `d26380ee63c585dd` | `scripts/promote_result.py` |
| `e1e578134cf31ac1` | `src/fmharness/figures.py` |
| `262aab6ed0483cb1` | `src/fmharness/synthetic.py` |
| `b7024779d526dfce` | `src/fmharness/statistics.py` |
| `0d85084410356f6a` | `scripts/alpine/delta_reproducibility.sbatch` |
| `48e9af48bc65f465` | `scripts/alpine/permutation_null.sbatch` |
| `5319051d72b3b62a` | `tests/test_rung0_controls.py` |
| `1726797c592dcaec` | `tests/test_rung0_figures.py` |
| `60b4d293efe5df1a` | `tests/test_verify_rung0.py` |
| `bb05664ebe800f61` | `tests/test_permutation_null.py` |
| `f0d59f0f5ba07bff` | `tests/test_statistics_known_answers.py` |
| `0ff1d59cdcc7cb05` | `tests/test_project_rules.py` |
| `2e3a4c5f5a4183a6` | `docs/tasks/rung0-assay-reliability/summary.ipynb` |
| `178e964df82c514f` | `docs/tasks/rung0-assay-reliability/verify.ipynb` |
| `ab4133ddffd5a00d` | `docs/tasks/rung0-assay-reliability/design.md` |
| `b7e9de3165960d20` | `docs/DATA.md` |
| `22094e07eec6c4f0` | `data/tranches/tahoe100m-pseudobulk-de.v1.json` |
| `9a8797a5698e2c56` | `data/tranches/tahoe100m-pseudobulk-de.v1.manifest.txt` |
| `c54e3df8ef24ab32` | `pyproject.toml` |
| `2f384d2d784a9b5c` | `.gitattributes` |

---

## Counts

**Canonical count: 138 claims.** One pass, one enumeration, no second count.

| Verdict | Count |
|---|---|
| ALIGNED | 106 |
| DEVIATION-RECORDED | 1 |
| DRIFT | 21 |
| PENDING RUN | 10 |
| **Total** | **138** |

By source: `design.md` D1–D119 (119 claims), `docs/SPEC.md` S1–S19 (19 claims).

PENDING RUN is not one of `docs/audit.md`'s three verdicts. It is used here only for claims whose
truth is a number the unfinished run has not produced; every such claim is listed again in its own
section and none is counted as passing.

---

## Clause verdicts — `design.md`

### Header and scope

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D1 | Task `rung0-assay-reliability`, status OPEN, branch `rung0-assay-reliability` | ALIGNED | `git rev-parse --abbrev-ref HEAD` → `rung0-assay-reliability`; `docs/SPEC.md:60` names it OPEN; `docs/STATE.md:20` ladder row "In progress" |
| D2 | Steps: build, split, select, score, decompose, null, document, promote | ALIGNED | `design.md:4`; the run implements each — `build_split_half_frame`, `hash(plate)%2` split, `responder_mask`, `score_split_half`, `decompose_noise`, `stratified_null_draws`, `verify_rung0.py`, `promote_result.py` |
| D3 | Supersedes the unmerged branch `rung0-replicate-ceiling` | ALIGNED | `docs/SPEC.md:60`; `decisions.md` first entry, dated 2026-09-01 |

### Data

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D4 | One tranche, registered `tahoe100m-pseudobulk-de.v1` | ALIGNED | `data/tranches/tahoe100m-pseudobulk-de.v1.json` `tranche_id` |
| D5 | Content hash over all 1,026 downloaded files | ALIGNED | Recomputed above: sha256 of the manifest = registered `content_hash`; `wc -l` = 1026 |
| D6 | Registered under `data/tranches/`, described in `docs/DATA.md` | **DRIFT** | Both files exist, but both describe the *superseded* rung, not this one. `data/tranches/…v1.json` sets `"drug_count": 32` with the description "drug_count = the declared drug panel" — this task declares no drug panel and the corpus carries ~1,100 perturbations (`docs/DATA.md:26`). `docs/DATA.md:69-71` states "rung 0's 1,600 scored pairs is 1,650 minus those 50 unscoreable Ribociclib rows" — that is the superseded branch's result, presented as rung 0's, in the registry document the design points a reader to. `docs/DATA.md` is dated "As of 2026-08-28", before the design |
| D7 | Per (cell line, drug, dose, plate) it carries each gene's log2 fold change against plate-matched solvent controls | ALIGNED | `docs/DATA.md:38-44`; DMSO plate-matched per `docs/DATA.md:28-29` |
| D8 | "reads **five** of its columns: `log2FoldChange`, `lfcSE`, `padj`, `baseMean`, and the keys (`gene_name`, `plate`, `concentration`, `Cell_ID_DepMap`, `drug`)" | **DRIFT** | Recount: the enumeration names **nine** columns, not five. The code reads all nine — `delta_reproducibility.py:178-185` (`Cell_ID_DepMap`, `drug`, `gene_name`, `log2FoldChange`, `padj`, and the replicate column), `:195-200` (`lfcSE`, `baseMean`, the dose column), `:846` (`concentration` via `DOSE_CANDIDATES`). The sentence reads as five items and is a nine-column list; this is exactly the wrong-count-that-reads-fine class `docs/audit.md:20-22` names |

### What is measured

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D9 | Replicate plates split into two groups per (line, drug) | ALIGNED | `delta_reproducibility.py:179-180`, `FILTER (WHERE hash(plate) % 2 = 0/1)` |
| D10 | Each group's per-gene log2FoldChange averaged | ALIGNED | `avg(log2FoldChange) FILTER (...)`, `:179-180` |
| D11 | The two averaged profiles correlated across genes (Pearson) | ALIGNED | `masked_rowwise_pearson`, `:375-406`; `score_split_half`, `:436-438` |
| D12 | Spearman-Brown corrected, `2r/(1+r)` | ALIGNED | `src/fmharness/statistics.py:103`; guarded entry `spearman_brown_or_nan`, `delta_reproducibility.py:751-760` |
| D13 | Read against stratified mismatched-condition nulls | ALIGNED | `stratified_null_draws`, `:548-595`; consumed in `summarize`, `:794-796` |
| D14 | All-gene reliability over every gene the table carries | ALIGNED | `main`, `:1219-1221`: with no `--panel-file`, `panel = set(de["gene_name"].unique())`; no top-variance fallback, and `:1209-1212` says so explicitly |
| D15 | Responder reliability over genes the **first** group calls DE: `padj < 0.05` in at least one of that group's (plate, dose) rows | ALIGNED | `min(padj) FILTER (WHERE hash(plate) % 2 = 0) AS padj0`, `:165`; `responder_mask` thresholds that minimum, `:494-504`. `min` over the group's rows is exactly "at least one". Pinned by `tests/test_rung0_controls.py:193` |
| D16 | Both correlations reported raw and Spearman-Brown corrected | ALIGNED | `summarize` emits `splithalf_mean_r` and `spearman_brown_full` under both `all_` and `responder_` prefixes, `:798-824`, called twice at `:1262-1263` |
| D17 | Three quarters of conditions split one plate against two | PENDING RUN | The number is a property of the full-extent pool, which does not exist yet. The 3/4 figure comes from the superseded 32-drug pool. `pool_description` emits `n_plates_half0`/`n_plates_half1` (`:854-855`) so the run can settle it |
| D18 | The corrected value reported again over even-plate-count conditions | ALIGNED | `n_plates_even` in `pool_description`, `:853`; joined to conditions in `main`, `:1240-1243`; `spearman_brown_full_even_plates` in `summarize`, `:807-809`; asserted by `tests/test_rung0_controls.py:927-928` |
| D19 | `lfcSE` is one plate's within-contrast standard error and cannot see plate-to-plate variation | ALIGNED | `docs/DATA.md:45-46`; restated at `delta_reproducibility.py:298-303` |
| D20 | `sigma^2_plate = var_across_plates(lfc) - mean(lfcSE^2)`, floored at zero, over (line, drug, dose, gene) with ≥2 plates | ALIGNED | SQL `var_samp(log2FoldChange)`, `avg(lfcSE*lfcSE)`, `HAVING count(DISTINCT plate) >= 2` at `:196-204`; `greatest(..., 0.0)` at `:253`; pandas twin `decompose_noise:367` |
| D21 | Dose held fixed, so a dose effect cannot masquerade as plate noise | ALIGNED | `_noise_select`, `:191-205`: `dose_grp = f", {dose}"` enters the `GROUP BY` at `:202`. `DOSE_CANDIDATES` (`:827`) resolves to `concentration`, which the table carries (`docs/DATA.md:44`). Pinned by `tests/test_rung0_controls.py:477` (a ±2.0 per-dose shift must return `max(sigma2_plate) < 1e-9`) |
| D22 | Reported as the between-plate fraction "aggregated over genes within a condition and over conditions" | **DRIFT** | `noise_aggregate:263` takes a flat `avg(between_plate_fraction)` over every (line, drug, dose, gene) row. That is a single mean over gene-conditions, not a mean over genes **within a condition** then over **conditions** — the two differ whenever conditions carry unequal gene counts, and only the second weights each condition equally. The 2026-09-01 decision entry (task 3/7, scale) records moving the aggregation into the engine and sampling the per-gene table; it does not record dropping the two-stage weighting |
| D23 | Stratified by expression (`baseMean`) and by response size | ALIGNED | `noise_aggregate:279-283`, `ntile(4) OVER (ORDER BY base_mean)` and `ntile(4) OVER (ORDER BY abs(mean_lfc))`, written to `rung0_noise_strata.csv` (`:1333`); figure panels `figures.py:861-885` |

### Inclusion rules

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D24 | All genes: every gene the table carries | ALIGNED | `main:1220`; the sbatch passes no `--panel-file` (`delta_reproducibility.sbatch:60-64`) |
| D25 | A gene contributes when finite in both groups | ALIGNED | `masked_rowwise_pearson:390`, `ok = np.isfinite(a) & np.isfinite(b)` |
| D26 | All-gene: a condition is scored at ≥50 qualifying genes | ALIGNED | `min_genes` default 50 (`:1148`), applied at `:405` |
| D27 | Responder genes: finite in both groups **and** called DE by the first group | ALIGNED | `ok &= select` at `:391-392`, ANDed into the finiteness mask before any moment |
| D28 | Responder: a condition is scored at ≥50 qualifying genes | ALIGNED | Same `n < min_genes` test at `:405`, with `n` counted after the mask |
| D29 | Drugs: every drug for which a condition has at least two distinct plates, "so a split exists" | **DRIFT** | No such filter is implemented, and the stated rule does not imply its consequence. Admission is by `dropna(subset=["lfc0","lfc1"])` (`:1204`), i.e. **at least one plate in each hash group**. A condition with two or three distinct plates that all hash to group 0 has "at least two distinct plates" and no split, and is silently dropped. The implemented rule is stricter than the declared one, so the design over-states which conditions are admitted |
| D30 | All fifty cell-line keys, including the literal string `NA` | **DRIFT** | The measurement path is clean: no filter on `Cell_ID_DepMap` anywhere in the SQL (`:178-185`, `:850-859`), Arrow keeps `'NA'` a string (`_compact_df:109-113`), and both verifiers guard the CSV round trip (`verify_rung0.py:152` `keep_default_na=False, na_values=[""]`; `verify.ipynb` uses `csv.DictReader`, no pandas). But `summary.ipynb`'s own reader does not: `table()` calls `pd.read_csv(p, **kw)` and **no call site ever passes `kw`**, so pandas' default `keep_default_na=True` turns the literal `NA` key into `NaN`. At `summary.ipynb` L180 the reviewer's notebook prints `pool['patient'].nunique()`, which excludes NaN — **49 cell lines, against the design's fifty**, in the one document written to convince a reader the fifty are there. The count "fifty" itself is PENDING RUN |
| D31 | Doses pooled: averaging within each plate group runs over the screen's three doses | ALIGNED | `GROUP BY Cell_ID_DepMap, drug, gene_name` at `:185` — no dose key, so `avg` spans every (plate, dose) row in the group |
| D32 | A dose-resolved reliability is out of scope | ALIGNED | No dose-keyed reliability anywhere; the only dose-keyed aggregation is the decomposition |
| D33 | Replicate unit is the plate; assignment by `hash(plate) % 2`, one fixed split per condition | ALIGNED | `:159` prints the rule; `:179-180` applies it. `hash` is deterministic and global, so the split is fixed and identical across conditions |

### Figure rules

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D34 | Every step declares a positive control, a negative control, and figures | ALIGNED | `design.md:96-167`; enforced for build/split/score/null by `tests/test_project_rules.py:289-341` (passes) |
| D35 | A figure is drawn from a committed table | **DRIFT** | Holds for every figure but one. `fig_build`'s fold-change panel (d) takes `delta_real` and `delta_synthetic`, built in memory at `main:1369-1372` from `de["lfc0"]` (strided) and `pos["lfc0"]`, and **written to no file** — `write_audit_checksums:1121-1125` only hashes files, so nothing carries them. A reader cannot recompute the distribution that panel displays, which is the one thing this rule exists to guarantee. Every other figure's inputs land as CSV (`rung0_pool_description.csv`, `rung0_per_pair_r.csv`, `rung0_padj_sample.csv.gz`, `rung0_leakage_control.csv`, `rung0_control_per_pair.csv`, `rung0_control_noise.csv.gz`, `rung0_noise_per_gene.csv.gz`, `rung0_null_draws.csv`, `rung0_effect_terciles.csv`, `rung0_mde_curve.csv`) |
| D36 | A figure with a control shows it in the same panel or the one beside it, on shared axes | ALIGNED | `figures.py:342` shared bins (build), `:717` `sharex=ax_hist` (score), `:808` `sharex/sharey` (decompose), `:937` shared bins (null) |
| D37 | Figures are produced by the run, not drawn by hand | ALIGNED | `main:1374-1387` |
| D38 | Figures land in `docs/tasks/rung0-assay-reliability/figures/` | ALIGNED | `main:1199` `fig_dir = out_dir / "figures"`; `delta_reproducibility.sbatch:59` `OUT_DIR` default is the task folder |
| D39 | The reviewer meets them in `summary.ipynb` in the order below | **DRIFT** | The order is right for what exists — `01_build` … `08_power` at `summary.ipynb` L184, 211, 256, 301, 343, 377, 434, 451, matching the design's step order. But the design's second **null** figure (D74) is never produced and therefore never met; and `09_per_gene_reliability.png`, which the run writes and both verifiers require, is displayed nowhere in the summary. The reviewer's path is neither complete over the declared figures nor closed over the produced ones |

### build

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D40 | positive: a synthetic pool with planted structure flows through the real DuckDB builder and comes out with the planted shape | ALIGNED | `tests/test_rung0_controls.py:119` `test_build_positive_planted_pool_comes_out_with_the_planted_shape` — replicate column resolves, 12 conditions survive, 300 genes; and `:176` `test_build_admits_every_drug_when_no_drug_list_is_given` |
| D41 | negative: a pool with no plate replication yields no scoreable conditions | ALIGNED | `tests/test_rung0_controls.py:131` `test_build_negative_no_replication_yields_no_scoreable_pairs` |
| D42 | figure: plates per condition, conditions per cell line, conditions per drug | ALIGNED | `figures.py:320-336`, panels (a), (b), (c) |
| D43 | figure: `log2FoldChange` for a handful of real conditions with the synthetic pool's histogram beside it | ALIGNED | `figures.py:338-364`, shared bins. *(Note: the run passes a strided sample over the whole pool, not "a handful of conditions"; the panel is honest about being the real screen. Recorded, not verdicted drift — the substance, real beside synthetic on shared bins, is there.)* |
| D44 | figure: histogram of the fraction of genes DESeq2 could not test (`baseMean` zero) per condition | **DRIFT** | The panel is conditional on a column nothing produces. `figures.py:316` `has_untestable = "frac_untestable" in pool.columns`; `grep -rn frac_untestable scripts src` returns hits only in `figures.py` and in the **test fixture** `tests/test_rung0_figures.py:68`. `pool_description` (`:849-861`) emits `n_rows`, `n_plates`, `n_plates_even`, `n_plates_half0`, `n_plates_half1`, `n_dose_levels` — no untestable fraction. On the real run this panel silently does not exist, and `tests/test_rung0_figures.py:197` asserts only that dropping the column still writes a PNG |

### split

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D45 | positive: a planted pool splits into two populated groups | ALIGNED | `tests/test_rung0_controls.py:119` asserts 12 conditions survive `dropna(subset=["lfc0","lfc1"])`, i.e. both groups populated. *(Marked `step_build`, not `step_split` — see R11.)* |
| D46 | negative: a single plate cannot split and yields no scoreable conditions | ALIGNED | `tests/test_rung0_controls.py:131`, same caveat |
| D47 | figure: histogram of the two group sizes, showing the one-against-two imbalance | ALIGNED | `figures.py:406-419`, panel (a), from `n_plates_half0`/`n_plates_half1` |
| D48 | figure: genes finite in both halves per condition, with the 50-gene threshold marked | ALIGNED | `figures.py:422-437`, panel (b); `SCORING_THRESHOLD = 50` at `:68`, drawn at `:427-433` |

### select

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D49 | positive: planted responders with matching `padj` are recovered by selection from the first group alone | ALIGNED | `tests/test_rung0_controls.py:280` `test_selection_recovers_the_planted_responder_set` — every mask row is exactly {G0…G79} |
| D50 | positive: responder reliability sits above the all-gene reliability on the same pool | ALIGNED | `tests/test_rung0_controls.py:292`, `mean_resp > mean_all + 0.05` |
| D51 | negative: on a signal-free pool, selection admits at no more than the nominal FDR | ALIGNED | `tests/test_rung0_controls.py:318`, `mask.mean() ≈ 1-(1-0.05)**k` within 0.02. The `1-(1-α)^k` form rather than α is recorded in `decisions.md` (2026-09-01, implementation task 2) as a property of the declared rule |
| D52 | negative: responder reliability sits at its null | ALIGNED | Same test, `abs(mean(r_resp)) < 0.15` |
| D53 | leakage check: selecting on **both halves** of that signal-free pool returns a visibly inflated correlation | DEVIATION-RECORDED | The implemented control demonstrates selection on the **pooled** `|a+b|`, not on each half separately: `leakage_table:1086-1088`; test `tests/test_rung0_controls.py:351` `test_pooled_selection_inflates_a_signal_free_correlation` asserts `r_pooled > r_one + 0.2`. `decisions.md` (2026-09-01, implementation task 2) records precisely this — that truncating `|a|` and `|b|` independently does not inflate Pearson, and that the design's "both halves, or on the pooled data" treats two different things as one |
| D54 | figure: histogram of group-1 `padj` | ALIGNED | `figures.py:483-489`, panel (a), from `rung0_padj_sample.csv.gz` (`main:1312-1314`) |
| D55 | figure: responders per condition, with the 50-gene threshold marked | ALIGNED | `figures.py:491-508`, panel (b), threshold at `:498-504` |
| D56 | figure: overlap between the first group's responders and the second group's, captioned as a diagnostic | **DRIFT** | The panel is allocated and never drawn. `figures.py:478` `n_panels = 3 if overlap is None else 4`, `:480` unpacks only `axes[0], axes[1], axes[2]`, and `overlap` appears nowhere else in the function body — `grep -n overlap src/fmharness/figures.py` returns only the signature (`:453`), the docstring (`:471-472`) and that count (`:478`); `grep -e jaccard -e n_first -e n_both` returns nothing. `main:1376` passes the table, so the run writes a four-panel figure whose fourth panel is blank, with no caption. The docstring at `:471-476` describes a panel that does not exist |
| D57 | figure: the two-sided leakage value beside the one-sided value on the same signal-free pool | ALIGNED | `figures.py:510-537`, panel (c), from `rung0_leakage_control.csv` |

### score

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D58 | positive: a pool planted at full-data reliability `R` returns a half correlation of `R/(2-R)` and corrects back to `R`, on both gene sets, raw and corrected | ALIGNED | `tests/test_rung0_controls.py:850`, parametrised over R ∈ {0.2, 0.5, 0.8}; closed form derived in `src/fmharness/synthetic.py:14-20` |
| D59 | negative: planted zero signal returns null, and the correction leaves zero at zero | ALIGNED | `tests/test_rung0_controls.py:887`; `spearman_brown_or_nan(0.0) == 0.0` |
| D60 | figure: scatter of first half against second half for example conditions spanning the reliability range, **drawn twice — all genes, then that condition's responders** | **DRIFT** | Each example is drawn **once**, over all genes shared between the halves. `figures.py:648-685` iterates `example_ids` and draws one scatter per example into `grid[0, position]`; there is no second row, no responder variant, and no responder mask reaches this function. `example_pair_profiles` (`delta_reproducibility.py:613-712`) exports `example_id`, `gene`, `lfc0`, `lfc1` — no responder flag — so the responder half of the declared figure cannot be drawn from the committed table either |
| D61 | figure: each panel's own correlation printed and recomputable from the points plotted | ALIGNED | Printed at `figures.py:665-673`; written back with its points to `<name>.values.csv` at `:737-743`; recomputation pinned by `tests/test_rung0_figures.py:291` and `:310`, and by `verify_rung0.py`'s check "the score figure's printed correlations recompute from its own points" |
| D62 | figure: histogram of the per-condition correlation for both gene sets on real data | ALIGNED | `figures.py:689-714`, `r` and `r_responder` on shared bins |
| D63 | figure: the same two histograms from the **positive-control** pool (mass at the planted value) **and** from the **negative-control** pool (mass at zero) beneath them on shared axes | **DRIFT** | Both control pools are merged into **one** histogram. `main:1357-1367` builds `pos` (planted r = 0.5) and `neg` (signal-free), scores each, and concatenates them into `control_per_pair` with a distinguishing `control` column. `fig_score` then reads `control_r = _finite(control_per_pair, "r")` (`figures.py:627`) and draws a single histogram labelled "control pool, all genes" (`:719-726`); the `control` column is never read. The declared reading — mass at the planted value beside mass at zero — is destroyed: a bimodal blur is drawn instead, and the panel title ("where a planted answer lands") is wrong for half the mass |
| D64 | figure: raw and Spearman-Brown corrected means marked, with the even-plate corrected mean alongside | ALIGNED | `_SUMMARY_MARKS`, `figures.py:547-569`; drawn by `_mark_summary_lines`, `:572-583` — six labelled lines, three per gene set |

### decompose

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D65 | positive: planted plate offsets on top of known `lfcSE` recover `sigma^2_plate` and the planted fraction within tolerance | ALIGNED | `tests/test_rung0_controls.py:416` `test_decompose_recovers_a_planted_plate_variance`, tolerances 0.05 and 0.06 stated in the assertion |
| D66 | negative: plates differing only by sampling noise return a plate component at zero, not negative | ALIGNED | `tests/test_rung0_controls.py:454`, `all(sigma2_plate >= 0)` and mean < 0.03 |
| D67 | figure: histogram of the between-plate fraction | ALIGNED | `figures.py:802-805`, panel (a) |
| D68 | figure: scatter of `sigma^2_plate` against `mean(lfcSE^2)` per gene, identity line drawn | ALIGNED | `figures.py:820-859`, panel (b), identity at `:830-837` |
| D69 | figure: the same fraction stratified by expression and by response size | ALIGNED | `figures.py:861-885`, panels (c) and (d), row-aligned per the comment at `:790-792` |
| D70 | figure: "**Each** with its control pool beside it" | **DRIFT** | Only panel (a) gets a control. `figures.py:807-815` adds one control panel, for the fraction histogram alone; panels (b), (c) and (d) are drawn from the real screen only. `control_noise` is computed once (`main:1381`) and reaches `fig_decompose` as a single frame whose only use is `_finite(control_noise, "between_plate_fraction")` (`figures.py:795-796`) — the control's `sigma2_plate`, `mean_se2`, `base_mean` and `mean_lfc` columns, which panels (b)–(d) would need, are never read |

### null

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D71 | positive: planted drug-shared and line-specific components recover matched > same-drug > different-drug | ALIGNED | `tests/test_rung0_controls.py:571`, chained margins `observed > same_drug + 0.05 > diff_drug + 0.10`; permutation side `tests/test_permutation_null.py:191` |
| D72 | negative: signal-free data sits at its floors and the observed clears neither | ALIGNED | `tests/test_rung0_controls.py:593`; `tests/test_permutation_null.py:83`, `:216` |
| D73 | figure: matched against all three null strata, one panel per gene set, each stratum's mean marked | ALIGNED | `figures.py:903-973`; `_STRATUM_STYLE:896-900` carries all three; means at `:947-964`; one panel per gene set at `:926-929` |
| D74 | figure: the permutation check's null distribution drawn on the same axes as the bootstrap's | **DRIFT** | The function exists and is tested; **no run script calls it**. `grep -rn fig_permutation_vs_bootstrap scripts src tests` returns `figures.py:57` (`__all__`), `figures.py:976` (def) and six lines in `tests/test_rung0_figures.py` — nothing in `scripts/`. `main` in `delta_reproducibility.py` draws figures 01–09 (`:1374-1387`) and has no permutation means to draw from; `scripts/permutation_null.py:511-532` writes `rung0_permutation_perm_means*.csv` and draws no figure at all. `verify_rung0.py:70-80` lists nine required figures and this is not among them, and `tests/test_rung0_controls.py:1087-1097` asserts the same nine — so nothing would notice its absence. The design effect is asserted in prose and never made visible, which is the exact thing this figure was declared to prevent |

### document

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D75 | positive: the verification battery recomputes **every** claim from the run's artifacts alone and reports pass | **DRIFT** | Two promoted claim families are not recomputed. (a) **The p-values.** `p_vs_null` and `p_vs_same_drug` appear nowhere in `scripts/verify_rung0.py` or in `verify.ipynb` (grep: zero hits in both). `check_null_floors` (`verify_rung0.py:260-319`) recomputes the three floors and then asserts only the ordering `mean > floor_diff and mean > floor_same` (`:299-306`) — the significance claim itself is never re-derived. (b) **The MDEs.** `:307-318` reads both keys and asserts only `isfinite(m) and m > 0`; `rung0_mde_curve.csv` is read by neither verifier. Two committed evidence tables — `rung0_responder_overlap.csv` and `rung0_mde_curve.csv` — are read only by `summary.ipynb` (L249, L444), so the numbers the summary prints from them are asserted and unverified. 51 checks pass; "every claim" is not among them |
| D76 | negative: a claim perturbed in the summary fails that mechanical check | ALIGNED | `tests/test_verify_rung0.py:127` `test_a_perturbed_claim_fails_the_battery` — rewrites `all_splithalf_mean_r` to `"0.123"` and requires failures named `"all: mean recomputes"` and `"checksum recomputes"` |
| D77 | figures: none of its own; every figure is placed in `summary.ipynb` beside the table it was drawn from | ALIGNED | `summary.ipynb` pairs each figure with the table above it (e.g. `rung0_pool_description.csv` L177 → `01_build.png` L184). Completeness of that set is D39's finding, not this one |

### promote

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D78 | positive: a promoted copy byte-identical to the task-side table passes, **and the provenance record's checksums recompute from the files and match the ones the audit recorded** | **DRIFT** | The first half is implemented (`tests/test_promote_result.py:71`, `:90`). The second half is not: `scripts/promote_result.py` never reads `audit_checksums.json` — `grep -n -e audit_checksum -e AUDIT_SUMS scripts/promote_result.py` returns nothing. `verify_rung0.py`'s `check_promotion` (`:692-735`) recomputes the record's own checksums but does not compare them to the audit's; and it SKIPs entirely until `results/rung0-assay-reliability/` exists |
| D79 | negative: promotion refuses when the two copies differ, when the record is incomplete, or **when a checksum has moved since the audit read it** | **DRIFT** | Two of the three refusals exist — `promote_result.py:99` (destination checksum differs) and `:105` (record already exists), with schema validation covering "incomplete". The third does not exist anywhere. This is the refusal `docs/audit.md:28-32` names as the thing that closes the window opened by having the audit read uncommitted artifacts; without it that window is unchecked, and the plan's Task 9 Step 2b, which specified it, was not carried out |

### The three spanning checks

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D80 | Exports: the per-condition table carries each condition's own correlations, responder count and effect size | ALIGNED | `per_pair_table:507-545` emits `r`, `r_responder`, `n_responders`, `mean_abs_delta`, `n_genes_scored`; graded plant pinned by `tests/test_rung0_controls.py:650` |
| D81 | Exports: the decomposition table's between-plate fractions recompute from the same rows | ALIGNED | `verify_rung0.py` check "noise: sigma2_plate = max(var_lfc - mean_se2, 0) row by row" (`:432`), over a 200,000-row prefix |
| D82 | Exports: the null-draw table preserves each stratum's count and mean | ALIGNED | `null_draw_table:598-610`; `tests/test_rung0_controls.py:707` |
| D83 | Exports: each example profile reproduces its own correlation from its exported points, **over both its full gene set and its marked responders** | **DRIFT** | The full-gene-set half holds: `r_full`/`r_shown` in the index (`:702-709`), pinned by `tests/test_rung0_controls.py:733`. The responder half does not exist: `example_pair_profiles` exports no responder marking (`:688-696` writes `example_id`, `gene`, `lfc0`, `lfc1` only), so there is no marked responder set for a reader to recompute over. Same root cause as D60 |
| D84 | Exports: the build cache returns the frame that was built and never one built from different inputs | ALIGNED | `frame_cache_key:945-953` keys on paths, drug set and replicate column; `tests/test_rung0_controls.py:782` asserts `built.equals(cached)` and that a different drug set or replicate column yields a different key |
| D85 | Empirical control: conditions in thirds by response size; the split-half mean must rise | ALIGNED | `effect_size_tercile_table:996-1030`; `tests/test_rung0_controls.py:605`; `verify_rung0.py:340` check "reproducibility rises with effect size" |
| D86 | Its figure is the tercile means with their confidence intervals | ALIGNED | `figures.py:1075-1113`; intervals from `ci_lo`/`ci_hi` in `rung0_effect_terciles.csv` |
| D87 | Power: every promoted comparison reports its MDE at α = 0.05, power = 0.80, from the same null bootstrap as its p-value | ALIGNED | `minimum_detectable_aggregate` defaults `alpha=0.05, power=0.80` (`statistics.py:61-62`) and reuses `bootstrap_aggregate_pvalue`'s construction (`:85-94`); emitted at `delta_reproducibility.py:819-822` from the same `nl` pool used for `p_vs_null` at `:794-795` |
| D88 | Each of the two reliabilities reports its own MDE at its own condition count | ALIGNED | `summarize` is called once per gene set with that set's own finite `r` (`:790-795`, `:819-822`, `:1262-1263`); pinned by `tests/test_rung0_controls.py:901` |
| D89 | Its figure is MDE against condition count with both observed counts marked | ALIGNED | `mde_curve_table:1033-1064` emits an `observed` flag; `figures.py:1116-1158` marks each gene set's observed count at `:1141-1151` |

### Statistics machinery

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D90 | `src/fmharness/statistics.py` carries `bootstrap_aggregate_pvalue`, `minimum_detectable_aggregate`, `spearman_brown` | ALIGNED | `statistics.py:22`, `:54`, `:97` |
| D91 | Each has known-answer tests | ALIGNED | `tests/test_statistics_known_answers.py` — 7 tests, module-level `pytest.mark.known_answer` at `:24`; `:27`/`:35`/`:43` cover the bootstrap, `:57` Spearman-Brown, `:64`/`:77`/`:89` the MDE |

### Null strata

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D92 | Three mismatched-condition nulls: any pair; different drug and line; same drug, different line | ALIGNED | `stratified_null_draws:578-582` |
| D93 | The reported p-values read each reliability against the second and third | ALIGNED | `summarize:794-796` — `p_vs_null` against `diff_drug`, `p_vs_same_drug` against `same_drug` |
| D94 | A mismatched responder draw uses the **first** condition's selected genes, intersected with the second condition's finite genes | ALIGNED | `stratified_null_draws:592-593` — `sel = select[ii[pick]]` indexes by the row whose first half is used; `masked_rowwise_pearson` then ANDs `isfinite(b[jj])`. Pinned by `tests/test_rung0_controls.py:932`, which asserts every draw is a value the first-condition mask can produce and none is union-only |
| D95 | An exact permutation check — 500 permutations of the pairing, once pooled and once within each stratum — reports the dependence as a design effect | ALIGNED | `permutation_null:163-258` (pooled, no fixed points via `sample_permutation:61`), `stratified_permutation_null:268-416` (`sample_within_drug_permutation`, `sample_cross_permutation`); `--n-perm` default 500 (`:434`), not overridden by the sbatch; `design_effect` at `:239` and `:394` |

### Run and promotion

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D96 | **One** cluster job (`scripts/alpine/delta_reproducibility.sbatch`) with the build cache enabled and no gene or drug file | **DRIFT** | The named job is correct in itself: `--frame-cache` at `delta_reproducibility.sbatch:62`, no `--panel-file`, no `--drugs-cid-file`, and the script's own `--drugs-cid-file` default is `""` (`delta_reproducibility.py:1137`). But the run is **two** jobs, not one: `scripts/alpine/permutation_null.sbatch` is a separate `#SBATCH` script that the plan's Task 8 submits with a dependency, and it has no build cache — `permutation_null.py` has no `--frame-cache` flag and calls `dr.build_split_half_frame` directly (`:461-463`), so it repeats the ~40-minute full scan, and its `for GENE_SET in all responder` loop (`permutation_null.sbatch:55`) repeats it **twice more** |
| D97 | then `scripts/permutation_null.py` for the permutation check | **DRIFT** | It runs, but not on the pool the design declares. `permutation_null.py:422` sets `--drugs-cid-file` default `"data/static/tahoe_target_cids.txt"` — **not** empty, unlike its sibling — and `permutation_null.sbatch:57-60` does not override it, while its own comment at `:52-53` claims "No drug list and no panel". `resolve_drug_names` returns `None` only if that path does not exist (`delta_reproducibility.py:63-64`); `docs/DATA.md:64-65` states the file is "on Alpine, in the repository checkout; not tracked on this branch". If it is present in the Alpine working tree, the permutation check silently scores the superseded rung's 33 drug names while the reliability it is checking scores every drug — the two would be computed on different pools, and nothing in the job would say so |
| D98 | Outputs land in the task folder | ALIGNED | `OUT_DIR` default `docs/tasks/rung0-assay-reliability` in both sbatch scripts (`:59` and `:54`) |
| D99 | Three results promoted with `scripts/promote_result.py` | PENDING RUN | Nothing promoted; `results/rung0-assay-reliability/` does not exist. The reversal making the decomposition the third promoted result is recorded in `decisions.md` (2026-09-01) |
| D100 | Each provenance record's inputs are the tranche content hash and nothing else | PENDING RUN | No provenance record exists |
| D101 | Its arguments record the inclusion choices and, for the responder row, the selection rule and its `padj` threshold | PENDING RUN | The material exists in the params sidecar — `_write_params_sidecar` extra at `main:1270-1277` carries `selection_rule`, `gene_inclusion`, `drug_inclusion`, `dose_handling`, and `args` carries `padj_threshold`. Whether it reaches the provenance record's `arguments` cannot be checked until promotion |
| D102 | The reviewer's path is `summary.ipynb`: the hypotheses, then each step's figures beside their table, then the conclusions | **DRIFT** | The spine is right (steps 1–6 at L162, 192, 219, 281, 309, 351) and the four hypotheses are stated up front (L96-103). But the hypotheses are paraphrased, not carried over: hypothesis 3 truncates `design.md:230` — "with lower correlations **over either all genes or the differentially expressed genes**" becomes "with lower correlations" — and, structurally, **hypothesis 3 is answered nowhere in the notebook**, while hypothesis 4's answer (L479) is the pointer "see the decomposition above" rather than a number. Two of four hypotheses have no answering cell, which is checkable now and independent of the run |
| D103 | `verify.ipynb` recomputes inline from committed artifacts and imports nothing from this project's own code | ALIGNED | Recomputed above: 19 code cells, imports are `csv`, `gzip`, `hashlib`, `json`, `math`, `statistics`, `subprocess`, `sys`, `pathlib` — nine, all standard library, no `fmharness`, no `scripts`, no numpy or pandas. The script runs only as the final cross-check cell, via `subprocess` |

### Expected result — the four hypotheses

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D104 | All-gene correlations will be low | PENDING RUN | No run |
| D105 | Responder correlations higher than all-gene | PENDING RUN | No run |
| D106 | Noise higher where correlations are lower | PENDING RUN | No run — and no cell computes it (D102) |
| D107 | Aggregate noise higher in responders than over all genes | PENDING RUN | No run — answered by pointer only (D102) |

### Ported apparatus

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D108 | Carried over by path; nothing re-typed; each file arrives with its tests | ALIGNED | Commit `640a428` "port: rung-0 apparatus from the superseded branch" adds all 24 listed paths in one change, code and tests together |
| D109 | `delta_reproducibility.py`: "the select and decompose steps and the figures are **the only additions** to a carried-over file" | **DRIFT** | `git diff 640a428..HEAD -- scripts/delta_reproducibility.py --stat` → 900 insertions, 140 deletions. Sixteen functions were added; seven belong to none of the three named additions: `_drug_predicate` and `resolve_drug_names` (the build's drug admission — a **build**-step change), `_compact_df` (Arrow dictionary encoding, recorded in `decisions.md` 2026-09-01 but not reflected in this sentence), `mde_curve_table` (power), `effect_size_tercile_table` (the empirical control with intervals), `spearman_brown_or_nan` (the score correction's guard), `write_audit_checksums` (the audit interface). The sentence is a completeness claim and it is false |
| D110 | `permutation_null.py` ported from `derangement_null.py`, renamed throughout — file, sbatch, test, output names | ALIGNED | `scripts/permutation_null.py`, `scripts/alpine/permutation_null.sbatch`, `tests/test_permutation_null.py`; outputs `rung0_permutation_*` at `permutation_null.py:520-532`. No `derangement` remains: `grep -rn derangement scripts src tests` returns nothing |
| D111 | `src/fmharness/statistics.py` — significance, power, Spearman-Brown | ALIGNED | Present; see D90 |
| D112 | `register_tranche.py`, `promote_result.py`, `src/fmharness/schema/` — provenance machinery | ALIGNED | All three present, all with tests |
| D113 | `scripts/alpine/ralpine`, `scripts/alpine/*.sbatch` — cluster boundary and jobs | ALIGNED | Present; boundary pinned by `tests/test_ralpine_boundary.py` |
| D114 | `verify_rung0.py`, `test_verify_rung0.py`, `verify.ipynb` — rebuilt around this task's outputs | ALIGNED | `8654f7c` rewrites the battery; `EXPECTED_CHECKS = 51` (`tests/test_verify_rung0.py:62`); `verify.ipynb` created there and its outputs stripped in `c3e55b2` |
| D115 | The seven listed test files exist | ALIGNED | All seven present; `uv run pytest -q` collects and passes them |
| D116 | `docs/DATA.md`, `data/tranches/` — dataset registry and pinned tranche | ALIGNED | Both present. Their *content* is D6's finding, not this one |

### Out of scope

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| D117 | A dose-resolved reliability is not built | ALIGNED | No dose key in `build_split_half_frame`'s `GROUP BY` (`:185`) and no dose-keyed reliability anywhere; the only dose-keyed aggregation is the decomposition, which the design puts in scope |
| D118 | Sensitivity of either reliability to the choice of split is not built | ALIGNED | One fixed `hash(plate) % 2` split; no sweep, no alternative-split parameter, no test of split sensitivity |
| D119 | An external deposit of the half-profile matrix is not built | ALIGNED | Only a handful of example conditions are exported (`example_pair_profiles`, `:613`), plus per-gene summary r (`per_gene_reliability`, `:735`). The full (condition × gene) pivots are never written |

---

## Clause verdicts — `docs/SPEC.md`

| # | Claim | Verdict | Evidence |
|---|---|---|---|
| S1 | Split each condition's replicate plates in two and correlate the two groups' per-gene deltas | ALIGNED | See D9–D11 |
| S2 | Per-pair Pearson, aggregated as the **mean over pairs** — the one statistic every delta rung reports | ALIGNED | `summarize:799-800` reports `splithalf_mean_r` as `float(np.mean(r))` over finite conditions |
| S3 | Spearman-Brown corrected to full-data reliability | ALIGNED | Applied to the mean, not per condition then averaged: `summarize:804` `spearman_brown_or_nan(mean)`; the rationale is stated at `:780-782` and pinned by `tests/test_rung0_controls.py:924` and `:1012-1015` (the identity `2m/(1+m)` against the reported mean, tolerance 2e-3, tighter than the ~0.005 Jensen gap a per-condition implementation would produce) |
| S4 | Against stratified mismatched-pair nulls | ALIGNED | See D92 |
| S5 | Alongside a planted-signal positive control | ALIGNED | `main:1357-1367` scores the planted pool in the run and commits `rung0_control_per_pair.csv`; `tests/test_rung0_controls.py:507` |
| S6 | At the assay's full extent, every drug with plates enough to split | ALIGNED | `_drug_predicate:116-125` returns an empty predicate when no drug list is given; sbatch passes none |
| S7 | Two gene sets: every gene the table carries, and each condition's first-group responders | ALIGNED | See D14, D15 |
| S8 | Selected from that group alone, never from the group the correlation is scored against | ALIGNED | Traced every path: `padj0` is the only selection input (`:165`, `:442-453`, `:494-504`); `padj1` exists solely for the overlap diagnostic (`:456-491`) and `responder_mask` has no parameter that could admit it. `main:1227-1228` and `permutation_null.py:481-482` both build the mask from `padj_pivot`. The forbidden pooled rule appears only inside `leakage_table` (`:1067-1106`) on a synthetic pool. Pinned by `tests/test_rung0_controls.py:1116`, which flips `padj1` to its complement, asserts the mask is unchanged, and — crucially — asserts the overlap diagnostic *does* move, so the test cannot pass vacuously |
| S9 | Promoted beside them: a decomposition of each delta's variance into between-plate and within-plate parts | ALIGNED (structure) | Computed and exported (`noise_aggregate`, `rung0_noise_decomposition.csv`); the reversal to promoting it is recorded in `decisions.md` 2026-09-01. The promotion itself is PENDING RUN under D99 |
| S10 | Rung 0 measures at full extent only; a later rung declares its own restriction | ALIGNED | `--panel-file` retained but unused by this task's jobs (`:1213-1221`, sbatch `:53-58`); the `restrict` step's removal from this task is recorded in `decisions.md` 2026-09-01 |
| S11 | Passing means both ceilings significantly above the mismatched-pair null | PENDING RUN | No numbers |
| S12 | The task is named in the spec tree | ALIGNED | `docs/SPEC.md:60` links `docs/tasks/rung0-assay-reliability/design.md`; `tests/test_project_rules.py:125` passes |
| S13 | Rule 1 — every promoted result carries its provenance record in `results/<task-slug>/` | PENDING RUN | `tests/test_project_rules.py:73` and `:100` skip: no promoted results |
| S14 | Rule 2 — reversals dated in the task's `decisions.md` | ALIGNED | 14 dated entries in `decisions.md`, including two explicit reversals (the decomposition's promotion; Spearman-Brown's restoration); `tests/test_project_rules.py` rule-2 tests pass |
| S15 | Rule 3 — the README stays in step | ALIGNED | `README.md` changed in this task's commits (`git diff 54a3fec...HEAD --stat` shows README.md +44/−…) and names the branch and design at `README.md:62-63`; the rule-3 tests pass |
| S16 | Rule 4 — every measurement step a task touches carries a positive and a negative control | ALIGNED | `tests/test_project_rules.py:289` passes over this task's design. *(Coverage note, not a departure: `MEASUREMENT_STEPS` at `:284` is SPEC's own list — `load, build, restrict, split, fit, score, null` — so the scan does not bind `select` or `decompose`, the two steps this task introduced. Both are nevertheless declared in the design and implemented as controls; see R12.)* |
| S17 | Rule 4 — every promoted comparison reports its MDE at the declared α and power beside its p-value | ALIGNED | `summarize:815-822` emits `p_vs_null`, `p_vs_same_drug`, `mde_80_vs_diff_drug`, `mde_80_vs_same_drug` in one row |
| S18 | Each step's marker is registered in `pyproject.toml` alongside the test that uses it | ALIGNED | `pyproject.toml:68-77`, nine markers. *(`step_split` is registered and used by no test — see R11.)* |
| S19 | Frame: a fraction-of-ceiling is valid only when numerator and denominator share a frame and a panel | ALIGNED | Not exercised by this rung, which computes no fraction; the tranche content hash is the frame key and is pinned |

---

## Reverse-direction findings

**Base used.** `git merge-base HEAD origin/main` → `0a9a12ade770d7399ebe4e8db6f8da48066d1ce4`. That diff is
**51 files**, but it is not this task's surface: the branch also carries the *unmerged project-documents
work* (commits `5ff88a4`…`54a3fec`, 15 files — `.coderabbit.yaml`, `codecov.yml`, `.github/workflows/ci.yml`,
`docs/PROCESS.md`, `docs/SPEC.md`, `docs/STATE.md`, `docs/adapter_contract.md`, `docs/environment.md`,
`src/fmharness/schema/*`, `uv.lock`, `README.md`, `pyproject.toml`, `.gitignore`, `tests/test_project_rules.py`),
which belongs to a different task. I therefore ran the attribution over **both**: the full
merge-base diff, and the task's own surface `git diff 54a3fec...HEAD --stat` — **44 files**, every one
of which I attributed below.

Files accounted for by a claim: `scripts/delta_reproducibility.py`, `scripts/permutation_null.py`,
`scripts/verify_rung0.py`, `scripts/promote_result.py`, `scripts/register_tranche.py`,
`scripts/download_tahoe_pseudobulk_de.py`, `src/fmharness/statistics.py`, `src/fmharness/figures.py`,
`src/fmharness/synthetic.py`, both sbatch jobs, `scripts/alpine/ralpine`,
`scripts/alpine/{00_target_cids,01_pseudobulk_shortcut,register_tranche}.sbatch`, the seven test files,
`tests/test_rung0_figures.py`, `data/tranches/*`, `docs/DATA.md`, `design.md`, `plan.md`,
`decisions.md`, `docs/audit.md`, `docs/decisions.md`, `docs/PROCESS.md`, `docs/SPEC.md`,
`docs/STATE.md`, `README.md`, `summary.ipynb`, `verify.ipynb` (D34–D116, S12–S18).

Changed under no claim, or contradicting one:

| # | Finding | Evidence |
|---|---|---|
| R1 | **`scripts/alpine/rung0_probe.sbatch` is undeclared.** A 77-line cluster job added by `d1250c2` and edited by `89deeec`, in neither the design's ported-apparatus table nor its run section, and with no decision entry | `git log --name-status`; `design.md:238-247` |
| R2 | **`write_figure` is dead code.** Defined at `delta_reproducibility.py:864`, called by nothing (`grep -rn write_figure scripts src tests` → the def, one docstring mention, and one plan line). Plan Task 6 said it would be replaced; it was left behind. It also draws a line labelled "mean (headline)" (`:875`), a term `decisions.md` retired on 2026-09-01 | grep |
| R3 | **An undeclared table and figure ship with the run.** `rung0_per_gene_reliability.csv` (`main:1308-1309`) and `09_per_gene_reliability.png` (`:1387`) appear in no claim in `design.md`. Both are nevertheless *required* by `verify_rung0.py:79` and by `tests/test_rung0_controls.py:1096`, so an undeclared artifact is enforced while a declared one (D74) is not. `summary.ipynb` never shows the figure | grep, `verify_rung0.py:70-80` |
| R4 | **Two stale cross-references in `scripts/permutation_null.py` contradict this task.** `:29-33` states "the sbatch job always passes `--panel-file`" — this task's sbatch passes none, by design. `:9`, `:186` cite `docs/tasks/rung0-replicate-ceiling/verification.md`, a path that does not exist in this tree; `:18-22`, `:112` describe "~32 drugs" and "~1,600 rows", the superseded pool | reading the file |
| R5 | **`scripts/alpine/permutation_null.sbatch:13` cites `docs/tasks/rung0-replicate-ceiling/design.md`**, a path not in this tree, and `:18-19` describes "the same-drug stratum's ~32-drug clustering" | reading the file |
| R6 | **Two docstrings in `delta_reproducibility.py` cite a `design.md` section that no longer exists.** `:741` "Unpromoted (see design.md)" and `:886-888` "(design.md, 'per-gene reliability')". The design has no per-gene reliability section | grep of `design.md` |
| R7 | **The module usage example contradicts the design.** `delta_reproducibility.py:13-15` shows the script invoked with `--drug-names-file` and `--panel-file`, the two things this rung must not pass | reading the file |
| R8 | **`src/fmharness/synthetic.py`'s stated model is not the one implemented.** The docstring at `:10-12` says `delta = signal + plate_offset + sampling_noise` with "`plate_offset` shared by every gene on a plate"; `planted_split_half_frame:69-71` plants signal and noise only, no plate offset. The closed forms the score controls check against are still correct, but the generative model a reader is told about is not the one in the code | reading the file |
| R9 | **`docs/STATE.md` points at a file that does not exist and states a rule the tree contradicts.** It says figures are "pointed at from `verification.md`" — no `verification.md` exists in the task folder (plan Task 9 declares one) — and that decisions are "dated and appended to the bottom of the task's own `design.md` and `plan.md`", where SPEC rule 2 and the tree put them in `decisions.md` | `ls docs/tasks/rung0-assay-reliability/`; `docs/STATE.md` "Where things live" |
| R10 | **Scope the decisions removed is still implemented.** `decisions.md` (2026-09-01) records that "the `restrict` step, its control and the restriction section leave this task", but `tests/test_rung0_controls.py:535` and `:553` still ship `test_restrict_positive_panel_subset_scores_exactly_the_subset` and `test_restrict_negative_disjoint_panel_scores_nothing`, and run on every suite | grep |
| R11 | **`step_split` is a registered marker no test uses.** `pyproject.toml:72` registers it; no test in the repository carries it, so `-m step_split` — the selector SPEC rule 4 tells a task to run for the steps it touched — selects nothing for the split step. `step_promote` is likewise carried only by `tests/test_project_rules.py:64` and `:88`, both of which skip before promotion | grep over `tests/` |
| R12 | **The rule-4 scan cannot see two of this task's steps.** `MEASUREMENT_STEPS` (`tests/test_project_rules.py:284`) is SPEC's list and omits `select` and `decompose`, the two measurement steps this task introduces and registers markers for. Not a departure from SPEC as written — recorded because the enforcement a reader would assume covers all eight declared steps covers four |
| R13 | **A large generated table is on course to be committed.** `--noise-sample-rows` defaults to 2,000,000 (`delta_reproducibility.py:1176`), written as `rung0_noise_per_gene.csv.gz`. `decisions.md` (2026-09-01) ruled out committing the ~100 MB half-profile matrix on exactly this ground; nothing bounds this file's size, and `.gitattributes` marks `docs/tasks/**/*.csv.gz` generated but does not keep it out of the repository |
| R14 | **`.gitattributes` and `uv.lock` changed under no claim.** `.gitattributes` (15 lines, `9015cbd`) implements PROCESS §1's review/evidence split, which no design claim names; `uv.lock` (+2,222 lines, `640a428`) is a dependency change with no entry in `decisions.md` |
| R15 | **`verify.ipynb` carries two small defects** a reader will meet: a tuple-wrapped `verdict()` argument at file line ~718-723 that prints `("… match",)`, and a dead re-assignment at ~497. Neither changes a verdict boolean |
| R16 | **`permutation_null.py` has no `--padj-threshold`.** It hard-codes `dr.responder_mask(padj)` at `:482`, i.e. α = 0.05, while `delta_reproducibility.py` takes `--padj-threshold` (`:1149`). The design fixes the threshold at 0.05, so the two agree today; a change to one would silently not reach the other |

---

## PENDING RUN

Claims whose truth is a number the unfinished run has not produced. None is verdicted, and none
counts as passing. The audit is not complete until each is checked against the run's artifacts,
with those artifacts' checksums recorded here as `docs/audit.md` requires.

| # | Claim | What must be checked once the run lands |
|---|---|---|
| D17 | Three quarters of conditions split one plate against two | Recompute from `rung0_pool_description.csv` (`n_plates_half0`, `n_plates_half1`). The 3/4 figure is inherited from the superseded 32-drug pool and may not survive full extent |
| D30 (count) | Fifty cell-line keys, one of them the literal `NA` | `rung0_pool_description.csv` read with `keep_default_na=False`. If it is 49 or 51, the design's sentence moves — and `summary.ipynb` will print the wrong number regardless until D30's drift is fixed |
| D99 | Three results promoted | `results/rung0-assay-reliability/` after gate 2 |
| D100 | Each provenance record's inputs are the tranche content hash and nothing else | The `.provenance.json` files |
| D101 | The record's arguments carry the inclusion choices and the selection rule | The `.provenance.json` files, against the params sidecars |
| D104 | All-gene correlations will be low | `rung0_reliability.csv`, `all_splithalf_mean_r` |
| D105 | Responder correlations higher than all-gene | `responder_splithalf_mean_r` against `all_splithalf_mean_r` |
| D106 | Noise higher where correlations are lower | Requires a cell that does not yet exist (D102) |
| D107 | Aggregate noise higher in responders | Requires a number the summary does not yet compute (D102) |
| S11 | Both ceilings significantly above the mismatched-pair null | `p_vs_null` and `p_vs_same_drug`, both of which the verification battery does not currently recompute (D75) |

Two further things the run must be read for, which are not claims but would invalidate several
above if they went wrong:

- **Confirm the permutation job scored every drug.** Read the resolved arguments printed by
  `permutation_null.py:455-456`. If the log says a drug count rather than "all drugs (no drug list
  given)", D97's hazard fired and the permutation check is on a different pool from the
  reliability.
- **Confirm the dose column resolved.** `noise_aggregate:237-242` prints a warning and proceeds if
  no dose column is found, in which case D21's claim does not hold for that run.

---

## Fix wave

**Not started. Awaiting the fix wave.**

Twenty-one drift items stand: D6, D8, D22, D29, D30, D35, D39, D44, D56, D60, D63, D70, D74, D75,
D78, D79, D83, D96, D97, D102, D109. Each is fixed, recorded as a dated decision, or ruled not a
defect with its reasoning, and the commit that carried it named here.

---

## Re-audit

**Not started. Awaiting the fix wave.**

A second fresh reader re-checks **only** the twenty-one items verdicted DRIFT above
(`docs/audit.md`, "The cap" — the re-audit does not re-enumerate the design), plus the PENDING RUN
section once the run's artifacts exist, and records its own verdict with the artifacts' checksums.
**The audit is not passed until the re-audit says so.**
