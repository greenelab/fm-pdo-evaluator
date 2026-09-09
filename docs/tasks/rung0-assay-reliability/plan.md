# Rung 0 — implementation plan

> **For agentic workers:** use `superpowers:subagent-driven-development` or
> `superpowers:executing-plans` to work this task by task. Steps use checkbox (`- [ ]`) syntax.

**Goal** Measure two split-half reliabilities of the Tahoe delta — over all genes, and over each
condition's responders — and decompose replicate noise into its between-plate and within-plate
parts, all from one pinned tranche, with the controls and figures the design declares.

**Architecture** The ported `scripts/delta_reproducibility.py` keeps its shape: DuckDB aggregates
the 89 GB table in-engine into a compact per-(line, drug, gene) frame, everything after that is
vectorized numpy over two pivots. Three things are added — a `padj` column in the build so
responders can be selected, a second aggregation grouped by dose so plate variance can be
separated from `lfcSE`, and a per-condition gene mask so one scorer serves both gene sets. Plot
primitives go in a new `src/fmharness/figures.py`; the script decides which figures to write.

**Tech stack** Python 3.13, DuckDB (in-engine aggregation, bounded memory), pandas, numpy,
matplotlib. Alpine `acpu` partition, CPU only.

**Spec** [`design.md`](design.md) · project rules in [`docs/SPEC.md`](../../SPEC.md) · lifecycle in
[`docs/PROCESS.md`](../../PROCESS.md).

## Global constraints

- **No embedded source.** PROCESS §1: this plan carries interfaces, invariants, ordered steps and
  expected test outcomes, never full implementations. Code exists once, in the repository. Five of
  the superseded rung's audit findings were plan-copy versus shipped-copy divergences.
- **Every step declares a positive and a negative control** (SPEC project rule 4), implemented as
  tests that import the real shipped function, never a reimplementation.
- **Every promoted comparison reports its MDE** at α = 0.05, power = 0.80, from the same bootstrap
  as its p-value.
- **No gene panel and no drug list.** Every gene the table carries, every drug with at least two
  distinct plates, all fifty cell-line keys including the literal string `NA`, doses pooled for
  the reliabilities and held fixed for the decomposition.
- **Selection is one-sided.** A condition's responders come from its first plate group alone.
- **Nothing is committed from a run until promotion** (PROCESS §1). Outputs live in
  `docs/tasks/rung0-assay-reliability/` and `figures/` in the working tree, uncommitted, through
  verify, audit and summarise.
- **Local gates before every push**: `uv run pytest -q`, then
  `git ls-files '*.py' '*.ipynb' | xargs uv run ruff check`, the same with `ruff format --check`,
  then `uv run pyright`.
- **Acronyms spelled at first use** in every document this task writes.

---

### Task 1: The build reads every drug, every gene, and `padj`

**Files**
- Modify: `scripts/delta_reproducibility.py` — `build_split_half_frame`, `pool_description`,
  `main`
- Test: `tests/test_rung0_controls.py`

**Interfaces**
- Produces: `build_split_half_frame(paths, target_names, repl_col, tmp, memory_limit) ->
  tuple[pd.DataFrame, str]` where `target_names` may now be `None` or empty, meaning *no drug
  filter*. The returned frame gains one column, `padj0: float` — the **minimum** `padj` over the
  first plate group's (plate, dose) rows for that (line, drug, gene). Existing columns
  (`patient`, `drug`, `gene_name`, `lfc0`, `lfc1`) keep their names and meaning.
- Produces: `pool_description(paths, target_names, repl, tmp)` with the same `None`-means-all
  behaviour, and one added column `n_plates_even: bool` — whether the condition's plate count is
  even, which Task 4 needs for the exact-split Spearman-Brown subset.

**Invariants**
- The minimum is the right aggregate for selection: the rule is "significant in **at least one**
  of the group's rows", so a gene qualifies when its smallest adjusted p-value clears 0.05.
- A gene DESeq2 could not test carries null in `padj` as well as `log2FoldChange`; DuckDB's
  `min` skips nulls, so such genes get `padj0` null and are excluded by the finiteness rule
  rather than by a filter.
- Dropping the drug filter changes the SQL from a `WHERE drug IN (...)` scan to a full scan of
  the same shards. It does not change the shard set read, so the roughly forty-minute build cost
  stands.

- [ ] **Step 1: Write the failing tests**

  In `tests/test_rung0_controls.py`, against the real `build_split_half_frame` over a synthetic
  parquet fixture written to a temp directory (the existing fixture builder in that file is the
  model — extend it rather than writing a second one):
  - `test_build_admits_every_drug_when_no_drug_list_is_given` — a fixture with three drugs, called
    with `target_names=None`, returns conditions for all three. Positive control for **build**.
  - `test_build_carries_the_first_groups_minimum_padj` — a gene planted with `padj` 0.9 on one
    first-group plate and 0.01 on another comes back with `padj0 == 0.01`; a gene planted 0.9 on
    both comes back 0.9. Pins that the aggregate is a minimum, not a mean.
  - `test_build_leaves_untestable_genes_null` — a gene with `baseMean` zero and null statistics
    comes back with null `padj0` and null `lfc0`, and is absent from any scored condition.
  - `test_build_yields_no_conditions_without_plate_replication` — a fixture whose every condition
    has one plate returns an empty scoreable set. Negative control for **build**.

- [ ] **Step 2: Run them and watch them fail**

  `uv run pytest tests/test_rung0_controls.py -k "build" -v`
  Expected: the drug-list and untestable-gene tests fail on the `WHERE drug IN` filter and the
  missing `padj0` column; `TypeError` or `KeyError`, not a wrong number.

- [ ] **Step 3: Make the drug filter optional and add `padj0`**

  In `build_split_half_frame`, build the `WHERE` clause conditionally so an empty or `None`
  `target_names` omits the `drug IN` predicate entirely, and add
  `min(padj) FILTER (WHERE hash(<repl>) % 2 = 0) AS padj0` to the select list. Mirror the same
  conditional predicate in `pool_description`, and add `n_plates % 2 = 0 AS n_plates_even`.
  In `main`, stop requiring a drug source: when neither `--drug-names-file` nor a readable
  `--drugs-cid-file` is given, pass `None` and print `all drugs (no drug list given)`.

- [ ] **Step 4: Run them and watch them pass**

  `uv run pytest tests/test_rung0_controls.py -k "build" -v` — expected PASS, 4 tests.
  Then the whole suite: `uv run pytest -q` — expected 77 passed, 3 failed (the three
  `test_verify_rung0.py` failures the port documented), 4 skipped.

- [ ] **Step 5: Commit**

  `git add scripts/delta_reproducibility.py tests/test_rung0_controls.py`, message stating the
  root cause: the build filtered to a drug list this task does not have, and never read the
  adjusted p-value the selection step needs.

---

### Task 2: Responder selection, one-sided

**Files**
- Modify: `scripts/delta_reproducibility.py` — new `responder_mask`, `score_split_half`
- Test: `tests/test_rung0_controls.py`

**Interfaces**
- Produces: `padj_pivot(de: pd.DataFrame, panel: set[str]) -> pd.DataFrame` — the `padj0` column
  pivoted to the same (condition × gene) shape and column order `score_split_half` produces, so
  the mask lines up with the pivots without a second alignment step.
- Produces: `responder_mask(piv_padj0: pd.DataFrame, alpha: float = 0.05) -> np.ndarray` — a
  boolean array shaped like the pivots, true where that condition's first group called that gene
  differentially expressed. Genes with null `padj0` are false.
- Produces: `score_split_half(de, panel, min_genes=50, *, select=None) -> tuple[np.ndarray,
  pd.DataFrame, pd.DataFrame]` — `select` is an optional boolean mask of the pivots' shape;
  entries false in it are treated as missing for the correlation, exactly as a non-finite value
  is. Passing `select=None` reproduces the current behaviour bit for bit.
- Produces: `masked_rowwise_pearson(a, b, min_genes, *, select=None) -> np.ndarray` with the same
  addition; the mask ANDs into the existing finiteness mask before any moment is computed.

**Invariants**
- **The mask is built from `padj0` only.** Nothing in this path may read `lfc1` or a second-group
  `padj`. This is the one invariant whose violation would silently inflate every responder number.
- Centring happens after masking: the mean subtracted from each row is the mean over the genes
  that row actually scores, never over all genes.
- A condition with fewer than `min_genes` selected genes returns NaN, and NaN conditions are kept
  in the per-condition table honestly rather than dropped.

- [ ] **Step 1: Write the failing tests**

  - `test_selection_recovers_the_planted_responder_set` — a pool with responders planted in a known
    gene subset and `padj0` planted to match returns exactly that subset from `responder_mask`.
    Positive control for **select**.
  - `test_responder_reliability_exceeds_all_gene_reliability_on_a_planted_pool` — on that same
    pool, the mean correlation over the selected genes is higher than over all genes.
  - `test_selection_admits_no_more_than_the_nominal_rate_on_signal_free_data` — a signal-free pool
    with uniform `padj0` selects at most 0.05 of genes plus binomial slack (assert on a seeded
    pool with an explicit tolerance, not on the expectation alone), and the responder correlation
    is not distinguishable from zero. Negative control for **select**.
  - `test_two_sided_selection_inflates_a_signal_free_correlation` — selecting on both halves of the
    same signal-free pool returns a mean correlation materially above the one-sided value. The
    **leakage check**: assert the inflation exists and record the number in the test's message, so
    a refactor that reintroduces two-sided selection turns this test red.
  - `test_select_none_reproduces_the_unselected_scorer_exactly` — `select=None` equals the
    pre-change output element for element.

- [ ] **Step 2: Run them and watch them fail** — `uv run pytest tests/test_rung0_controls.py -k
  "select or responder or leakage" -v`. Expected: `TypeError` on the unknown `select` keyword.

- [ ] **Step 3: Implement**

  Add the `select` keyword to `masked_rowwise_pearson` and thread it through `score_split_half`;
  add `responder_mask`. `score_split_half` gains a fourth returned pivot only if needed — prefer
  returning the `padj0` pivot from a small separate helper so the existing three-tuple signature
  and its callers (including `scripts/permutation_null.py`) keep working unchanged.

- [ ] **Step 4: Run them and watch them pass**, then the full suite as in Task 1.

- [ ] **Step 5: Commit** — message stating that selection reads the first group alone and why the
  two-sided variant is a test rather than an option.

---

### Task 3: Noise decomposition

**Files**
- Modify: `scripts/delta_reproducibility.py` — new `build_noise_frame`, `decompose_noise`
- Test: `tests/test_rung0_controls.py`

**Interfaces**
- Produces: `build_noise_frame(paths, target_names, repl_col, dose_col, tmp, memory_limit) ->
  pd.DataFrame` — one row per (line, drug, dose, gene) with at least two plates, carrying
  `var_lfc` (the sample variance of `log2FoldChange` across that group's plates, `var_samp`),
  `mean_se2` (the mean of `lfcSE * lfcSE` over the same rows), `n_plates`, and `base_mean`.
- Produces: `decompose_noise(noise: pd.DataFrame) -> pd.DataFrame` — the same rows plus
  `sigma2_plate = max(var_lfc - mean_se2, 0)` and
  `between_plate_fraction = sigma2_plate / var_lfc` (null where `var_lfc` is zero).

**Invariants**
- **Dose is a grouping key here, not pooled.** Pooling dose would put a dose effect into
  `var_lfc` and report it as plate noise.
- The estimator is unbiased under plate offsets plus independent sampling error: the sample
  variance across plates has expectation `sigma2_plate + mean(lfcSE^2)` exactly, for any set of
  per-plate standard errors. The floor at zero is what makes it a variance rather than a
  difference, and the negative control checks the floor is reached without being crossed.
- `var_samp`, not `var_pop` — the expectation identity above is for the sample variance.

- [ ] **Step 1: Write the failing tests**

  - `test_decompose_recovers_a_planted_plate_variance` — plates given offsets of known variance on
    top of sampling noise of known `lfcSE` return `sigma2_plate` within tolerance and the planted
    between-plate fraction within tolerance. Positive control for **decompose**. Use enough genes
    that the tolerance is tight, and state it in the assertion.
  - `test_decompose_floors_at_zero_without_going_negative` — plates differing only by the planted
    sampling noise return `sigma2_plate` at zero and never below. Negative control for
    **decompose**.
  - `test_decompose_does_not_charge_a_dose_effect_to_plate_noise` — a pool where each dose has a
    different mean response but plates within a dose are identical returns a plate component at
    zero. This is the test that pins dose as a grouping key.

- [ ] **Step 2: Run them and watch them fail** — `NameError: build_noise_frame`.

- [ ] **Step 3: Implement** the DuckDB aggregation and the arithmetic. The dose column is
  discovered from `DOSE_CANDIDATES` the way `pool_description` already discovers it; if no dose
  column exists, group by (line, drug, gene) and say so in the printed run log, because the
  decomposition then cannot exclude a dose effect and the design's claim would not hold.

- [ ] **Step 4: Run them and watch them pass**, then the full suite.

- [ ] **Step 5: Commit.**

---

### Task 4: Two reliabilities, Spearman-Brown on both

**Files**
- Modify: `scripts/delta_reproducibility.py` — `summarize`, `per_pair_table`,
  `effect_size_terciles`, `main`
- Test: `tests/test_rung0_controls.py`, `tests/test_statistics_known_answers.py`

**Interfaces**
- Produces: `summarize(r, nulls, seed=0, *, label="") -> dict` — unchanged in shape, with every
  key prefixed by `label` when one is given, so one function serves both gene sets and the
  summary row carries `all_*` and `responder_*` families rather than two files.
- Produces: `per_pair_table(piv0, piv1, r, *, r_responder=None, n_responders=None)` — the
  per-condition table gains `r_responder` and `n_responders` columns beside the existing `r`,
  `n_genes_scored` and `mean_abs_delta`.
- Consumes: `spearman_brown(r)` from `src/fmharness/statistics.py`, unchanged.

**Invariants**
- The correction is applied to the **mean over conditions**, not per condition then averaged;
  `2r/(1+r)` is not linear, so the two differ and the design's declared statistic is the mean.
- The even-plate-count subset is computed from `pool_description`'s `n_plates_even`, joined on
  (line, drug) — not re-derived, so one definition of "even" exists.
- Guard `r <= -1`: `spearman_brown` is undefined there and its docstring says callers guard.

- [ ] **Step 1: Write the failing tests**

  - `test_spearman_brown_round_trips_a_planted_full_data_reliability` (known answer) — a pool
    planted at full-data reliability `R` returns a half correlation of `R / (2 - R)` and a
    corrected value back at `R`, for `R` in {0.2, 0.5, 0.8}. Positive control for **score**, and
    the test that checks the correction rather than assuming it.
  - `test_zero_signal_returns_null_and_the_correction_leaves_zero_at_zero` — negative control for
    **score**.
  - `test_summary_carries_both_gene_sets_with_their_own_counts_and_mdes` — the summary dict has
    `all_n_pairs`, `responder_n_pairs`, both `*_splithalf_mean_r`, both `*_spearman_brown_full`,
    both `*_mde_80_vs_diff_drug`, and the responder count is not silently equal to the all-gene
    count.

- [ ] **Step 2: Run them and watch them fail.**

- [ ] **Step 3: Implement** the label prefix, the responder columns, and the even-plate subset
  value `*_spearman_brown_full_even_plates`.

- [ ] **Step 4: Run them and watch them pass**, then the full suite.

- [ ] **Step 5: Commit.**

---

### Task 5: Nulls and the permutation check on both gene sets

**Files**
- Modify: `scripts/delta_reproducibility.py` — `stratified_null_draws`
- Modify: `scripts/permutation_null.py` — `permutation_null`, `stratified_permutation_null`, `main`
- Test: `tests/test_rung0_controls.py`, `tests/test_permutation_null.py`

**Interfaces**
- Produces: `stratified_null_draws(piv0, piv1, n_perm=500, seed=0, min_genes=50, *, select=None)`
  — when `select` is given, a mismatched draw pairing condition *i*'s first group with condition
  *j*'s second group scores over **row *i*'s** selected genes, intersected with the genes finite
  in row *j*'s second group.
- Produces: `stratified_permutation_null(..., select=None)` with the same rule.

**Invariants**
- The selecting row is the one whose first group is used. Using row *j*'s mask, or the union,
  would apply a different selection rule to the null than to the observed value, and the
  comparison would no longer be like for like.
- The permutation stays a permutation with no fixed points; the rename to "permutation" does not
  relax that, and `sample_permutation`'s existing no-fixed-point assertion stays.

- [ ] **Step 1: Write the failing tests**

  - `test_null_ordering_is_recovered_from_planted_structure` — planting separate drug-shared and
    line-specific components recovers matched > same-drug null > different-drug null, on both gene
    sets. Positive control for **null**.
  - `test_signal_free_data_sits_at_its_floors` — negative control for **null**: the observed value
    clears neither floor.
  - `test_a_mismatched_responder_draw_uses_the_first_conditions_mask` — construct two conditions
    with disjoint planted responder sets and assert the null draw's scored gene count equals the
    first condition's responder count, not the second's and not the union's.

- [ ] **Step 2: Run them and watch them fail.**
- [ ] **Step 3: Implement.**
- [ ] **Step 4: Run them and watch them pass**, then the full suite.
- [ ] **Step 5: Commit.**

---

### Task 6: Figures

**Files**
- Create: `src/fmharness/figures.py`
- Modify: `scripts/delta_reproducibility.py` — replaces `write_figure`, `write_per_gene_figure`
- Test: `tests/test_rung0_figures.py`

**Interfaces**

One function per declared figure, each taking the committed table (a DataFrame) and an output
path, returning the path written. Named for the step it serves:

- `fig_build_pool_composition(pool: pd.DataFrame, out: Path)` — plates per condition, conditions
  per line, conditions per drug.
- `fig_build_delta_distributions(real: pd.DataFrame, synthetic: pd.DataFrame, out: Path)` — real
  `log2FoldChange` histograms beside the synthetic control pool's, shared axes.
- `fig_build_untestable_fraction(pool: pd.DataFrame, out: Path)`
- `fig_split_group_sizes(pool: pd.DataFrame, out: Path)` — the one-against-two imbalance.
- `fig_split_shared_gene_counts(per_pair: pd.DataFrame, out: Path)` — threshold at 50 marked.
- `fig_select_padj_and_responder_counts(per_pair, padj_sample, out)` — threshold marked.
- `fig_select_half_overlap(overlap: pd.DataFrame, out: Path)` — caption states it is a diagnostic
  and never an input to selection.
- `fig_select_leakage(leakage: pd.DataFrame, out: Path)` — one-sided beside two-sided on the same
  signal-free pool.
- `fig_score_example_scatters(profiles, index, out)` — first half against second half for the
  example conditions, drawn twice (all genes, responders), each panel's own correlation printed
  and recomputed from the plotted points.
- `fig_score_r_histograms(per_pair, control_per_pair, out)` — both gene sets on real data, the
  positive- and negative-control pools beneath on shared axes, with raw, corrected, and
  even-plate corrected means marked.
- `fig_decompose_fraction(noise: pd.DataFrame, out: Path)`
- `fig_decompose_scatter(noise: pd.DataFrame, out: Path)` — `sigma2_plate` against `mean_se2`,
  identity line drawn.
- `fig_decompose_strata(noise: pd.DataFrame, out: Path)` — by `base_mean` and by response size.
- `fig_null_overlays(per_pair, null_draws, out)` — matched against all three strata, one panel per
  gene set, stratum means marked.
- `fig_null_permutation_vs_bootstrap(perm_means, summary, out)` — the design effect, drawn.
- `fig_control_terciles(terciles: pd.DataFrame, out: Path)` — with confidence intervals.
- `fig_power_mde_curve(per_pair, null_draws, out)` — MDE against condition count, both observed
  counts marked.

**Invariants**
- **Every function takes a table and reads nothing else.** No function recomputes a statistic from
  raw data; if a number appears on a figure it came from a committed table, which is what makes it
  checkable. This is the design's first figure rule, enforced by the signature.
- `matplotlib.use("Agg")` before `pyplot` is imported, as the ported code already does.
- Every axis is labelled with units, and every figure that shows a control says which panel is the
  control in its title or legend.

- [ ] **Step 1: Write the failing tests** in `tests/test_rung0_figures.py`. These are export
  controls, not appearance tests: for each function, build a small table with a known answer, call
  the real function, and assert the file exists, is a non-trivial PNG, and — where the figure
  prints a number — that the number drawn is recoverable. Concretely, assert
  `fig_score_example_scatters` writes a companion `.values.csv` holding the points it plotted and
  the correlation it printed, and that recomputing Pearson from those points reproduces the
  printed value to four decimals. That companion file is what makes the design's "recomputable
  from the points plotted" a checkable claim rather than a promise.

- [ ] **Step 2: Run them and watch them fail** — `ModuleNotFoundError: fmharness.figures`.

- [ ] **Step 3: Implement** `src/fmharness/figures.py`, then replace the two ported figure
  functions in `delta_reproducibility.py` with calls into it and wire every figure into `main`,
  writing to `<out-dir>/figures/`.

- [ ] **Step 4: Run them and watch them pass**, then the full suite.

- [ ] **Step 5: Commit.**

---

### Task 7: Wire the run, and prove it on synthetic data locally

**Files**
- Modify: `scripts/delta_reproducibility.py` — `main`, `_write_params_sidecar`
- Modify: `scripts/alpine/delta_reproducibility.sbatch`, `scripts/alpine/permutation_null.sbatch`
- Test: `tests/test_rung0_controls.py`

**Invariants**
- PROCESS §3: **test on synthetic data before spending cluster time.** `main` runs end to end on
  the synthetic fixture, writing every table and every figure, before anything is submitted.
- The sbatch jobs pass no `--panel-file` and no drug file, and set `--frame-cache`.
- Outputs land in `docs/tasks/rung0-assay-reliability/` with figures under `figures/`.

- [ ] **Step 1: Write the failing test** — `test_main_writes_every_declared_artifact_on_a_synthetic
  _pool`: run `main` against the fixture with `sys.argv` patched, then assert the exact set of
  files the design declares exists — summary, per-condition table, null draws, noise
  decomposition, example profiles and index, pool description, and one file per figure in Task 6.
  A missing artifact fails by name.
- [ ] **Step 2: Run it and watch it fail.**
- [ ] **Step 3: Wire `main`**, extend the params sidecar with the selection rule and its
  threshold, and update both sbatch scripts.
- [ ] **Step 4: Run it and watch it pass**, then the full suite, then the lint, format and type
  gates over tracked files.
- [ ] **Step 5: Commit and push.** Pushing is required before Alpine can pull (PROCESS §2).

---

### Task 8: Run it on Alpine

**Files**
- Create (uncommitted, working tree): `docs/tasks/rung0-assay-reliability/*.csv`,
  `*.params.json`, `figures/*.png`, and the job logs.

- [ ] **Step 1: Confirm the tranche is still on scratch** — `ralpine ls` the DE directory and
  compare the file count against the tranche manifest's 1,026. If scratch has been purged, rebuild
  it with `scripts/alpine/01_pseudobulk_shortcut.sbatch` first and say so in `verification.md`.
- [ ] **Step 2: `ralpine update`** to pull the pushed branch onto the Alpine checkout.
- [ ] **Step 3: `ralpine submit scripts/alpine/delta_reproducibility.sbatch`**, then verify the
  job is queued and note its identifier.
- [ ] **Step 4: When it finishes, submit the permutation job** with
  `--dependency=afterok:<jobid>` and verify the dependency actually attached —
  `ralpine jobinfo <id> | grep Dependency`. A chain that silently drops its dependency runs
  immediately and out of order (PROCESS §2).
- [ ] **Step 5: Pull the outputs and both logs** into the task folder. **Do not commit them**
  (PROCESS §1, "What reaches GitHub, and when").
- [ ] **Step 6: Read the run log before reading the numbers** — the resolved arguments, the
  replicate and dose columns chosen, the condition counts, any warning. A silently degraded run
  invalidates everything downstream.

---

### Task 9: The verification battery, rebuilt

**Files**
- Modify: `scripts/verify_rung0.py`, `tests/test_verify_rung0.py`
- Create: `docs/tasks/rung0-assay-reliability/verify.ipynb`,
  `docs/tasks/rung0-assay-reliability/verification.md`

**Invariants**
- The battery recomputes **every** claim this task will promote, from the artifacts alone,
  printing claim / recomputed / pass-fail, on a laptop in about a minute.
- `verify.ipynb` recomputes **inline and self-contained** — standard-library hashing, direct table
  reads, explicit arithmetic, nothing imported from this project. A notebook that calls the script
  relocates the trust instead of discharging it. The script is its final cross-check cell.
- Committed **without outputs**.
- It carries no figures; figures belong to `summary.ipynb` (PROCESS §3, the two-notebook rule).

- [ ] **Step 0** Record, in `audit_checksums.json` beside the artifacts, the sha256 of every table
  the battery reads. Task 10's audit cites these and Task 12's promotion checks against them.
- [ ] **Step 1** Rewrite the battery's claim list against this task's summary row: both
  reliabilities raw and corrected, the even-plate corrected values, every null mean, both
  p-values, both MDEs, the tercile ordering, the between-plate fraction, the tranche content hash
  recomputed from the manifest, and the per-condition table's own arithmetic.
- [ ] **Step 2** Update `tests/test_verify_rung0.py` — the three tests that fail from the port are
  the ones that pin the claim count and the layer coverage; they turn green here or the battery is
  incomplete. Add the **document** step's negative control:
  `test_a_perturbed_claim_fails_the_battery` — copy the artifacts to a temp directory, alter one
  number in the summary row, run the real battery against them, and require a failure naming that
  claim. A battery that passes on altered evidence is checking nothing.
- [ ] **Step 2b** Add the **promote** step's new refusal to `scripts/promote_result.py` and
  `tests/test_promote_result.py`: promotion takes the checksums the audit recorded and refuses
  when an artifact's current checksum differs. Positive control — matching checksums promote.
  Negative control — a single altered byte is refused by name. This is what closes the window
  the design opened by having the audit read uncommitted artifacts; without it that window is
  unchecked.
- [ ] **Step 3** Write `verify.ipynb` with one cell per claim, recomputing inline.
- [ ] **Step 4** Run the battery and the notebook; paste the commands and their output into
  `verification.md`, with pointers to every table, figure and log.
- [ ] **Step 5** Full suite — expected 3 failed → 0. Commit the code and the documents, not the
  artifacts.

---

### Task 10: The audit

**Files**
- Create: `docs/tasks/rung0-assay-reliability/audit.md`
- Create: `docs/tasks/rung0-assay-reliability/review.md`

- [ ] **Step 1** Code review of the branch's diff (`superpowers:requesting-code-review`), findings
  and dispositions into `review.md`.
- [ ] **Step 2** First audit by a **fresh reader**, following [`docs/audit.md`](../../audit.md):
  numbered claims from `design.md` and SPEC's rung-0 section, both diff directions, three
  verdicts, evidence per claim, and the checksum of every artifact read.
- [ ] **Step 3** Fix wave for every drift item; dispositions recorded.
- [ ] **Step 4** Re-audit by a second fresh reader over the drift items only. The audit is not
  passed until it says so.
- [ ] **Step 5** Commit.

---

### Task 11: The summary notebook — gate 2

**Files**
- Create: `docs/tasks/rung0-assay-reliability/summary.ipynb`

- [ ] **Step 1** Write it in the design's order: the four hypotheses from **Expected result**,
  then build → split → select → score → decompose → null, each with its figures beside the table
  they were drawn from and its controls stated, then the conclusions, then the scripts touched.
- [ ] **Step 2** Say plainly, for each hypothesis, whether it held — including any that did not.
  A hypothesis that failed is the finding, not a defect to soften.
- [ ] **Step 3** Strip outputs, run it clean end to end, confirm every figure renders.
- [ ] **Step 4** Full suite, lint, format, type gates. Commit the notebook only.
- [ ] **Step 5** **Stop. This is gate 2** — Lucas reads the summary before anything is promoted.

---

### Task 12: Promote and open the draft pull request (after gate 2)

- [ ] Promote the three results with `scripts/promote_result.py`; commit them together with the
  figures, logs and provenance records in one change.
- [ ] Confirm every provenance checksum matches the one the audit recorded.
- [ ] `docs/STATE.md` and `README.md` in the same change, with the entry's three links.
- [ ] Open the pull request as a **draft**, description sending the reader to `summary.ipynb`
  first, then `design.md`, then `audit.md`, and separating review surface from generated evidence.
