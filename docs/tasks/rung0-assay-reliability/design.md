# Rung 0 — the reliability of the assay

**Task** `rung0-assay-reliability` · **Status** OPEN · **Branch** `rung0-assay-reliability`
**Steps** build, split, select, score, decompose, null, document, promote.
**Spec** [docs/SPEC.md](../../SPEC.md), rung 0 · **State** [docs/STATE.md](../../STATE.md)
**As of** 2026-09-01. Decision lineage in [`decisions.md`](decisions.md).

## What this rung establishes

How much of a measured drug response is signal rather than assay noise. A model predicting
these responses cannot agree with the measurement better than the measurement agrees with
itself, so this number is the ceiling every later score is read against.

This task supersedes the unmerged branch `rung0-replicate-ceiling`. That branch chose gene- and drug- pools from an archived lineage. 

## Data

One tranche: the Tahoe-100M pseudobulk differential-expression table, registered as
`tahoe100m-pseudobulk-de.v1` with a content hash over all 1,026 downloaded files
(`data/tranches/`; described in [docs/DATA.md](../../DATA.md)). Per cell line, drug, dose, and
plate it carries each gene's log2 fold change against plate-matched solvent controls. 

Each row is one DESeq2 result. This rung reads four of its statistics — `log2FoldChange`, its
standard error `lfcSE`, the Benjamini-Hochberg adjusted p-value `padj`, and `baseMean` — keyed by
`gene_name`, `plate`, `concentration`, `Cell_ID_DepMap` and `drug`. Nine columns of the sixteen;
`stat`, `pvalue`, `n_cells_trt`, `n_cells_ctrl` and the two alternative cell-line identifiers are
not read.

## What is measured

**The quantity.** For each (cell line, drug) condition, the replicate plates are split into two
groups, each group's per-gene log2 fold change is averaged, and the two averaged profiles are
correlated across genes (Pearson and Spearman-Brown corrected, `2r/(1+r)`). Statistics are read against stratified mismatched-condition nulls.

**Two gene sets** That correlation is computed twice from the same split:

- *All-gene reliability* — over every gene the table carries.
- *Responder reliability* — over only the genes the condition's **first** group calls
  differentially expressed: `padj < 0.05` in at least one of that group's (plate, dose) rows.

The hypothesis is that the correlation across all genes
is dominated by the genes that did not respond, sitting at zero plus noise, so it measures how
reproducibly the assay reports a mostly-null profile. Across the responding genes, it measures how
reproducibly the assay reports the response itself.

**Spearman-Brown, on both.** Dorrelations are reported raw and Spearman-Brown corrected, `2r / (1 + r)`.  Three quarters of conditions split one plate against two, so the corrected value is reported again over the conditions
with an even plate count — where the split is exact and the correction is not an approximation —
and the gap between the two is the size of the assumption.

**Noise decomposition** `lfcSE` is the standard error
of one plate's treated-versus-control contrast: cell-sampling error at that row's `n_cells_trt`
and `n_cells_ctrl`. It cannot see plate-to-plate variation — culture day, handling, position. For each (line, drug, dose, gene) with at least two plates, the variance of
`log2FoldChange` across plates has expectation `sigma^2_plate + mean(lfcSE^2)`, so

    sigma^2_plate = var_across_plates(log2FoldChange) - mean(lfcSE^2), floored at zero

estimates the plate component alone. Dose is held fixed here, so a dose effect cannot masquerade
as plate noise. Reported as the fraction of delta variance that is between-plate, aggregated over
genes within a condition and over conditions, stratified by expression (`baseMean`) and by
response size.


**Inclusion rules — defined from the table, not inherited.**

- *All Genes.* Every gene the table carries. A gene contributes to a condition's
  correlation when it has a finite value in both groups.A condition is scored when at least 50
  genes qualify.
- *Responder Genes* in addition to having a finite value in each group, the first group also called it differentially expressed (padj < .05). A condition is scored when at least 50
  genes qualify. 
- *Drugs.* Every drug for which a condition has a plate in **each** hash group, which is what
  makes a split exist. Two distinct plates are necessary but not sufficient — both can hash to
  the same group — so the rule is stated as the engine applies it, not as the count it is often
  mistaken for.
- *Cell lines.* Every key in the table, including the one whose DepMap identifier is missing
  and appears as the literal string `NA` — it is a real line's data and carries a consistent key
  throughout. Fifty keys, counted exactly (`count(DISTINCT ...)`, job 31996456). An earlier
  approximate count read 45 and this clause was briefly corrected to it; the approximation was
  wrong and the original number right. See `decisions.md`, 2026-09-01.
- *Doses.* **Held fixed.** A scoreable unit is a (cell line, drug, dose) triple, and the plate
  split happens inside one. Doses were pooled until the screen was counted: 86.6% of (line,
  drug, dose) combinations sit on a single plate, so splitting a condition's plates while
  pooling dose put **different doses in the two halves for 99.7% of conditions**. That number is
  a dose-to-dose correlation, not a test-retest reliability, and it cannot be the ceiling a
  later rung divides by. Holding dose fixed costs base — 7,641 dose-conditions have two or more
  plates, against 18,350 splittable conditions when dose was pooled — and buys the quantity the
  rung is named for. It also gives the reliabilities and the noise decomposition one inclusion
  rule instead of two.
- *Replicate unit and split.* The plate, within one dose. Plates are assigned to groups by
  `hash(plate) % 2`, one fixed split per dose-condition. Most replicated dose-conditions have
  exactly two plates (7,441 of 7,641), so their groups are one plate against one — equal, which
  is the case Spearman-Brown's correction is exactly right for. Under the earlier dose-pooled
  rule the halves were unequal for 70% of conditions; that assumption now holds for nearly all
  of them.


## Figures, controls and power (project rule 4)

Every step declares three things: a positive control that plants a known answer and requires the
shipped code to recover it, a negative control that feeds signal-free or mismatched data and
requires null, and the figures that show a reviewer what that step did.

Two rules hold for every figure. **A figure is drawn from a committed table**, so a reader can
recompute the number it displays; no figure asserts a value that exists only inside a run. **A
figure that has a control shows it in the same panel or the one beside it**, on shared axes —
real data alone shows what the screen looks like, real data beside a planted answer shows whether
the machinery reads it correctly. Figures are produced by the run, not drawn by hand, and land in
`docs/tasks/rung0-assay-reliability/figures/`; the reviewer meets them in the summary notebook
(`summary.ipynb`) in the order below, which is the order the measurement happens in.

- **build** — the screen as it arrives.
  - *positive*: a synthetic replicate pool with planted structure flows through the real DuckDB
    split-half builder and comes out with the planted shape (expected condition count, both
    groups populated).
  - *negative*: a pool with no plate replication yields no scoreable conditions.
  - *figures*: histogram of plates per condition, and of conditions per cell line and per drug. Histogram of
    `log2FoldChange` for a handful of real conditions, with the synthetic pool's histogram beside
    it. Histogram of the fraction of genes DESeq2 could not test (`baseMean` zero) per
    condition.
- **split** — one condition becomes two half-profiles.
  - *positive*: a planted pool splits into two populated groups.
  - *negative*: a single plate cannot split and yields no scoreable conditions.
  - *figures*: histogram of the two group sizes, which shows the one-plate-against-two imbalance. Histogram of the number of genes
    finite in both halves per condition, with the 50-gene scoring threshold marked.
- **select** — the responder set, chosen from the first half alone.
  - *positive*: with responders planted in a known subset of genes and `padj` planted to match,
    selection from the first group alone recovers that subset, and the responder reliability sits
    above the all-gene reliability on the same pool.
  - *negative*: on a signal-free pool, selection admits genes at no more than the nominal
    false-discovery rate and the responder reliability sits at its null.
  - *leakage check*: selecting on both halves of that same signal-free pool must return a visibly
    inflated correlation.
  - *figures*: histogram of group-1 `padj`, and of responders per condition, with the 50-gene
    threshold marked. Overlap between the first group's responders and the second group's — a
    diagnostic only, never an input to selection, and the figure says so in its caption. A panel
    putting the two-sided leakage value beside the one-sided value on the same signal-free pool,
    which is the clearest statement of why selection is one-sided.
- **score** — the two reliabilities.
  - *positive*: a planted reliability of known closed-form value is recovered by the real scoring
    function within tolerance, on both gene sets, raw and Spearman-Brown corrected — a pool
    planted at full-data reliability `R` must return a half correlation of `R / (2 - R)` and a
    corrected value back at `R`, which tests the correction rather than assuming it.
  - *negative*: planted zero signal returns null, and the correction leaves zero at zero.
  - *figures*: Scatter of first-half against second-half delta for a few example
    conditions spanning the range of reliability, drawn twice — all genes, then that condition's
    responders — with each panel's own correlation printed on it and recomputable from the points
    plotted. Histogram of the per-condition correlation for both gene sets on real data, with the
    same two histograms from the positive-control pool (mass at the planted value) and from the
    negative-control pool (mass at zero) beneath them on shared axes. Raw and Spearman-Brown
    corrected means marked on those histograms, with the even-plate-count subset's corrected mean
    marked alongside, so the equal-halves assumption is read off the figure.
- **decompose** — what kind of noise the ceiling is made of.
  - *positive*: a pool with plate offsets of known variance planted on top of sampling noise of
    known `lfcSE` recovers `sigma^2_plate` within tolerance, and recovers the planted
    between-plate fraction.
  - *negative*: a pool whose plates differ only by the planted sampling noise returns a plate
    component at its floor of zero, and does not go negative.
  - *figures*: histogram of the between-plate fraction of delta variance. Scatter of
    `sigma^2_plate` against `mean(lfcSE^2)` per gene with the identity line drawn, which is where
    a reader sees directly whether plate effects or cell sampling dominates. The same fraction
    stratified by expression (`baseMean`) and by response size. Each with its control pool beside
    it, since a decomposition that cannot recover a planted split is not evidence about the real
    one.
- **null** — what the reliabilities are read against.
  - *positive*: planting separate drug-shared and line-specific components recovers the ordering
    matched > same-drug null > different-drug null.
  - *negative*: signal-free data sits at its floors and the observed clears neither.
  - *figures*: overlaid histograms of the matched correlations against all three null strata, one
    panel per gene set, with each stratum's mean marked — the whole significance claim in one
    picture. The permutation check's null distribution drawn on the same axes as the bootstrap's,
    which is the design effect made visible rather than asserted.
- **document** — the reviewer's path, and the last step before anything is committed as evidence.
  - *positive*: the verification battery recomputes every claim from the run's artifacts alone and
    reports pass.
  - *negative*: a claim perturbed in the summary fails that mechanical check rather than passing
    unnoticed.
  - *figures*: none of its own. This step is where every figure above is placed in front of the
    reviewer, in `summary.ipynb`, each beside the table it was drawn from.
- **promote** — the number becomes citable, after the summary has been read.
  - *positive*: a promoted copy byte-identical to the task-side table passes, and the provenance
    record's checksums recompute from the files and match the ones the audit recorded.
  - *negative*: promotion refuses when the two copies differ, when the record is incomplete, or
    when a checksum has moved since the audit read it.

Three checks span the steps rather than belonging to one:

- **Exports** — each committed evidence table is checked against what it summarises: the
  per-condition table carries each condition's own correlations, responder count and effect size
  (a graded plant a misaligned export would scramble); the decomposition table's between-plate
  fractions recompute from the same rows; the null-draw table preserves each stratum's count and
  mean; each example profile reproduces its own correlation from its exported points, over both
  its full gene set and its marked responders; the build cache returns the frame that was built
  and never a frame built from different inputs.
- **Empirical in-run control** — conditions stratified into thirds by response size; the
  split-half mean must rise across the thirds. An assay that cannot find more reproducibility
  where there is more signal is broken. Its figure is the tercile means with their confidence
  intervals, which shows whether the rise is real or within noise.
- **Power** — every promoted comparison reports its minimum detectable effect (MDE) at α =
  0.05, power = 0.80, from the same null bootstrap as its p-value, and each of the two
  reliabilities reports its own at its own condition count. Its figure is MDE against condition
  count with both observed counts marked, so a reader can see how much of the power is the
  screen's size rather than the effect's.

## Statistics machinery

`src/fmharness/statistics.py`: `bootstrap_aggregate_pvalue` (the observed aggregate against the
bootstrapped null aggregate at the observed sample size), `minimum_detectable_aggregate` (the
same bootstrap read in the other direction), and `spearman_brown` (a half-length reliability
corrected to full length, applied to both correlations). Each has known-answer tests.

## Null strata

Three mismatched-condition nulls, each a set of correlations between one condition's first group
and a different condition's second group: *any pair*; *different drug and line* (the generic-
structure floor); *same drug, different line* (the line-specificity floor — the stricter one,
since two lines given one drug share that drug's generic response). The reported p-values read
each reliability against the second and third. A mismatched draw for the responder reliability
uses the genes the *first* condition's first group selected, intersected with the genes finite in
the second condition's second group: the same selection rule as the matched pair. Because mismatched draws reuse the same half-profiles,
an exact permutation check — 500 permutations of the pairing, once pooled and once within each
stratum — measures the dependence the bootstrap ignores and reports it as a design effect. 

## Run and promotion

One cluster job (`scripts/alpine/delta_reproducibility.sbatch`) with the build cache enabled and
no gene or drug file: `scripts/delta_reproducibility.py` over the tranche on scratch, then
`scripts/permutation_null.py` for the permutation check. Outputs land in the task folder. Three
results are promoted with `scripts/promote_result.py` — the all-gene reliability and the
responder reliability, each raw and Spearman-Brown corrected, with its own condition count, null
p-values and MDE, and the even-plate-count corrected value beside it; and the noise
decomposition, whose between-plate fraction is a finding about the assay in its own right and is
cited as one. Each provenance record's inputs are the tranche content hash and nothing else —
none of the three has a panel to pin — and its arguments record the inclusion choices (all genes;
all splittable drugs; doses pooled) and, for the responder row, the selection rule and its `padj`
threshold.

The reviewer's path through all of it is `summary.ipynb` (PROCESS §1, Summarise): the hypotheses
above, then each step's figures beside the table they were drawn from, then the conclusions. The
claim-by-claim recomputation stays in `verify.ipynb`, which recomputes inline from the committed
artifacts and imports nothing from this project's own code.

## Expected result

- The correlations between plates of differential expression over all genes will be low, since most genes are not affected by the drugs. 
- Correlations between genes that were significantly different (responders) in group 1 will be higher than the correlations over all genes. 
- The noise will be higher in cell line - drug combinations with lower correlations over either all genes or the differentially expressed genes.
- The aggregate noise will be higher in the differentially expressed genes (responders) than over all genes. 

## Ported apparatus

Carried over by path from the superseded, unmerged branch `rung0-replicate-ceiling` on origin. Nothing is re-typed; each file arrives with the tests
that exercise it.

| Path | Role |
|---|---|
| `scripts/delta_reproducibility.py` | the measurement: build, split, select, score, decompose, null, reporting, evidence exports, figures, build cache. The select and decompose steps and the figures are new to this task, along with the supporting tables each figure is drawn from, the engine-side noise aggregation, and the scatter-based pivot that replaced `pivot_table` at this screen's scale — the additions are listed in `decisions.md` rather than claimed to be three |
| `scripts/permutation_null.py` | the permutation check of the dependence assumption. Ported from `scripts/derangement_null.py`; "derangement" is the mathematical name for a permutation that leaves nothing in place, and it is renamed to permutation throughout — file, sbatch script, test, and output names |
| `src/fmharness/statistics.py` | significance, power, Spearman-Brown |
| `scripts/register_tranche.py`, `scripts/promote_result.py`, `src/fmharness/schema/` | provenance machinery |
| `scripts/alpine/ralpine`, `scripts/alpine/*.sbatch` | cluster boundary and jobs |
| `scripts/verify_rung0.py`, `tests/test_verify_rung0.py`, `verify.ipynb` | the executable verification battery and reviewer notebook, rebuilt around this task's outputs |
| `tests/test_statistics_known_answers.py`, `tests/test_rung0_controls.py`, `tests/test_permutation_null.py`, `tests/test_promote_result.py`, `tests/test_register_tranche.py`, `tests/test_ralpine_boundary.py`, `tests/test_download_tahoe.py` | controls and known-answer tests |
| `docs/DATA.md`, `data/tranches/` | dataset registry and the pinned tranche |

## Out of scope

A dose-resolved reliability; sensitivity of either
reliability (all or responders) to the choice of split; an external deposit of the half-profile matrix. 
