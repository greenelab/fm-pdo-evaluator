# fm-pdo-evaluator — project spec

## The question

Foundation models trained on cell-line perturbation atlases are advertised as predicting how a tumour will respond to a drug.
The clinically meaningful test of that claim is not another cell-line benchmark: it is whether the prediction survives the move to **patient-derived tumour organoids** — a different substrate, grown from a real patient, screened on a different assay.

This project measures that move.
It produces two numbers and the distance between them: an **in-silico transfer score**, obtained now on an embargoed organoid screen held out from all model development, and a **prospective score** obtained later against organoids screened after the prediction is made.
The gap between a model's cell-line performance and its organoid performance is the quantity the project exists to report.

## Why a ladder, not a single test

Running a cell-line-trained model directly on organoids yields one number and no diagnosis.
When that number is poor — and it usually is — at least five distinct explanations are consistent with it:

1. the organoid measurement is not reproducible, so nothing could have been predicted;
2. the model cannot predict an unseen cell line even in its own training distribution;
3. the measurement platform changed;
4. the readout changed, from expression response to viability;
5. the substrate genuinely differs — tumour organoids are not cell lines.

Those are five different conclusions demanding five different responses, and one end-to-end number cannot separate them.

So the evaluation is built as a **ladder**: each rung adds exactly one distribution shift to the rung below it, and every rung is scored against the same reproducibility ceiling.
Where the score falls off localises the failure to a specific shift.
A model that clears rungs 0-3 and fails rung 4 has told us something about organoids; a model that fails rung 1 has told us something about the model.

## Three words this project uses precisely

| Word | Means | Lifecycle |
|---|---|---|
| **Rung** | A scientific question — one level of the ladder, adding one distribution shift. The experiment. | Answered, or not yet |
| **Step** | A stage of the evaluation pipeline that every rung passes through: load, build, restrict, split, fit, score, null, promote, release, document. The machinery. | Never "done"; rules attach here |
| **Task** | A unit of work with a spec, an owner and a definition of done. Touches one or more steps, serves one rung or is cross-cutting. | OPEN or CLOSED |

A rung is *what we are asking*; a step is *where in the machinery*; a task is *what someone is doing about it*.
The project rules below bind to steps, which is why they hold for every rung without being restated per rung.

## The ladder

Each rung states the question it settles, the shift it adds, how it is measured, and what a passing result means.
That is the whole design; the detail a rung needs to be run — gene panel, aggregation, its baseline and control lists — is written in that rung's own spec.
**Rungs are built one at a time, and a rung's spec arrives with its implementation**, so that the design of a measurement and the code producing it are reviewed together rather than a year apart.
Nothing above rung 0 is in this repository yet.

---

### Rung 0 — is the target reproducible?


**Question** How much of a measured perturbation response is signal rather than assay noise?
**Adds** Nothing — this is the reference every other rung is read against.
**Measure** Split the replicate pool in half, correlate the two halves' per-(line, drug) expression deltas on the declared gene panel, Spearman-Brown corrected to full-data reliability, against a mismatched-pair null.
**Passing means** A ceiling significantly above that null. A rung above it can then be reported as a fraction of what is achievable rather than as a fraction of a hypothetical 1.0.
**If it fails** Nothing above it is interpretable: an unpredictable target cannot distinguish a bad model from an unmeasurable outcome.

### Rung 1 — can a model predict an unseen cell line, in-distribution?


**Question** Given a cell line the model has never seen, can it predict that line's expression response to a drug?
**Adds** One shift: an unseen cell line. Same platform, same readout, same substrate.
**Measure** Correlation between predicted and measured delta on the held-out line, against a floor that must fail and a planted signal that must be recovered, reported as a fraction of rung 0's ceiling.
**Passing means** The method beats its floor and recovers a stated fraction of the ceiling — the minimum competence claim before any harder shift is worth testing.
**If it fails** The failure is about the model, not about organoids, and no result above this rung can be attributed to the substrate.

### Rung 2 — what does crossing a measurement platform cost?


**Question** How much predictive signal survives when the mapping is fitted on one transcriptomic platform and tested on another?
**Adds** One shift on top of rung 1: the platform (fit on L1000, test on Tahoe).
**Measure** The retained fraction of rung 1's in-platform score, per method, against a scrambled-line control.
**Passing means** Retention clearly separable from the scrambled control — the platform shift costs something quantifiable rather than everything.
**Why it exists** The organoid comparison contains a platform change bundled with the substrate change. This rung prices the platform part on its own, so rung 4's shortfall is not silently attributed to biology.

### Rung 3 — does the representation predict drug response, not just expression?


**Question** Does a representation carry information about how much a drug *kills* a line, beyond what the drug does on average?
**Adds** One shift: the readout changes from expression delta to viability (dose-response AUC).
**Measure** Drug×line interaction, residualising each drug's mean so that only line-specific ranking counts, corrected across all declared variants, reported against the screen-agreement ceiling between independent screens of the same lines.
**Passing means** Interaction significantly above zero after correction, at a stated fraction of the screen-agreement ceiling.
**Why it exists** Viability is the endpoint the organoid rung uses. A representation that predicts expression but not viability would fail rung 4 for reasons that have nothing to do with organoids.

### Rung 4 — does it transfer to patient-derived organoids?


**Question** Does a model trained on cell-line drug response predict drug response in patient-derived tumour organoids?
**Adds** The final shifts: substrate (immortalised line → patient organoid), patients, and assay.
**Measure** The same interaction statistic as rung 3, fitted on cell lines and evaluated on a frozen, embargoed organoid screen, on the drug axis the cell-line training set defines, reported as a fraction of both rung 3 and rung 0.
**Passing means** Either transfer holds at a measurable fraction — the useful positive result — or it does not, and the ladder below localises the loss to the substrate step rather than leaving it unexplained.
**Constraint** The organoid cohort is a frozen holdout under embargo. It is not looked at, tuned on, or unfrozen except under conditions stated in its spec, because a holdout inspected once is no longer a holdout.

---

### Cross-cutting work

Some work belongs to no single rung: the shared evaluation library every rung scores through, the fairness machinery that makes methods comparable, the release gate that keeps patient identifiers out of anything published, and the enforcement that turns the rules below from advice into tests.
Each is specified alongside the rung that first needs it, rather than up front, so that a piece of shared machinery is designed against a real requirement.

## Project rules — what every task must satisfy

Every script that computes or reports a number must satisfy these.
**Amend this section deliberately when a genuine exception is needed — silently violating one is unacceptable.**

Each rule names the harness step it governs and the test that enforces it.
**The step is the unit, not the task**: any task touching that step inherits its tests, whichever rung or dataset the task is about.
A task runs `uv run pytest tests/test_project_rules.py` in full — several are repository- or artifact-wide scans, which catch a violation the task never intended — and adds `-m` for the steps it touched, e.g. `uv run pytest -m "step_split or step_score"`.
The ten step markers are registered in `pyproject.toml`.
Each entry also names the **edge case its test cannot see, and the mechanism that closes it** — a scan over source text is weaker evidence than a behavioural check, and reading the two as equivalent is its own failure mode.
These rules stay generic: which table or driver violates one today is current state, and belongs in `docs/PROJECT_STATE.md` and the task that owns the fix, not here.
A rule a piece of work cannot yet satisfy is exempted in one place, in the test suite, with the reason and the task that will close it — never by quietly not running the check.

1. **One shared way of splitting samples into cross-validation folds, everywhere.**
   `fmharness.deltas.fold_assignment` produces it: sorted, deterministic, and it becomes leave-one-out once there are at least as many folds as cell lines.
   Never write the split by hand in a script (`{ln: i % n_folds ...}`), even when it looks equivalent.
   It stops being equivalent the moment two scripts happen to list the cell lines in a different order, nothing announces that, and two rungs simply stop being comparable.
   **Step** split. **Enforced by** `tests/test_project_rules.py::test_rule_01_fold_split_is_order_free_deterministic_and_degenerates_to_loo`, `::test_rule_01_no_analysis_code_writes_its_own_fold_map` (`-m step_split`).
   **Edge case** a split written in another form (`np.array_split`, a groupby, a hand-rolled loop) is invisible to a text scan.
   **Edge-case test** `::test_rule_01_edge_no_alternative_split_mechanism_in_scripts` — enforced: the named alternative constructors (`KFold`, `ShuffleSplit`, `array_split`, …) fail the scan too, with exemptions listed by file and reason.
   The residual hole, a bespoke loop no pattern matches, closes structurally when the scoring entry point takes folds only from the shared helper.
2. **One correlation and one way of averaging it, applied identically wherever two numbers are compared.**
   If rung A reports a median Pearson and rung B a mean Spearman, then `A/B` or `A - B` measures nothing — it is noise with a sign.
   Fix the correlation and the averaging once, for the comparison, rather than separately inside each script that reports a piece of it.
   **Step** score. **Enforced by** `tests/test_project_rules.py::test_rule_02_transfer_penalty_is_a_difference_of_like_quantities`, `::test_rule_02_ladder_summary_publishes_no_number_without_saying_what_it_is_of` (`-m step_score`).
   **Edge case** nothing yet states each rung's correlation and averaging in a form code can read, so two numbers can be divided by accident and the result still looks like a measurement.
   **Edge-case test** `::test_rule_02_edge_metric_declaration_exists_and_covers_every_rung_table` — the acceptance criterion for that declaration: it fails until one exists and covers every rung that reports a number.
3. **Compare only on the samples, drugs, and genes that every method can cover.**
   Every method in one table must be scored on the same (cell line or organoid, drug) pairs and the same gene panel, never each method against whatever subset it happens to share with the measurements.
   (`restrict_common_support`, `common_gene_panel` + `assert_common_genes`.)
   **Step** build, restrict. **Enforced by** `tests/test_project_rules.py::test_rule_03_every_method_is_scored_on_the_identical_support`, `::test_rule_03_a_method_on_a_different_gene_panel_is_rejected` (`-m "step_build or step_restrict"`).
   **Edge case** the guarantee holds only where a caller remembers to invoke the helpers; a new scoring path that forgets them produces a well-formed table comparing different things.
   **Edge-case test** `::test_rule_03_edge_producing_code_restricts_support_and_panel` — enforced: every promoted comparison table's producing pipeline (the script its own sidecar names, its stage family, and their `fmharness` imports) must invoke both guard families.
   Moving the calls inside the scoring entry point remains the structural close, and turns this scan into a backstop.
4. **Every reported statistic needs a test with a known answer, and that test must call the same function the analysis calls.**
   Plant a signal and require it to be recovered; plant none and require the result to come back null.
   A test that re-derives the calculation inside itself passes even when the function it is supposed to be checking is broken — which has happened here: a control's own test stayed green while the control crashed on the cluster.
   **Step** score. **Enforced by** `tests/test_project_rules.py::test_rule_04_every_reported_statistic_has_a_known_answer_test` (`-m step_score`).
   **Edge case** a test can name a statistic without planting anything, so "a test exists" is weaker than "a signal was planted and recovered."
   **Edge-case test** `::test_rule_04_edge_known_answer_tests_carry_the_marker` — enforced: every reported statistic must be named in a test file declaring `pytest.mark.known_answer`, the author's signed claim that a signal is planted and recovered, reviewed once per file.
5. **The null for an average must be built by resampling that average, not by comparing it against single draws.**
   Use `fmharness.statistics.bootstrap_aggregate_pvalue`.
   A mean or median over n pairs varies roughly √n less than a single pair does, so scoring it against individual null draws inflates p by one to two orders of magnitude.
   The same mistake appeared independently in four scripts in one week before a shared, tested version of it existed.
   **Step** null. **Enforced by** `tests/test_project_rules.py::test_rule_05_aggregate_null_recovers_a_planted_shift_and_stays_null_without_one`, plus the fuller battery in `tests/test_statistics_recover_known_answers.py` (`-m step_null`).
   **Edge case** the defect is a mismatch of kinds, not a text pattern: the same expression is correct when the null holds replicates of the same aggregate and wrong when it holds single draws, so no search distinguishes them.
   **Edge-case test** `::test_rule_05_edge_every_manual_pvalue_site_is_allowlisted` — enforced: any p-value computed outside the shared helper fails unless its site is allowlisted with the verified reason it is correct.
   This retires the open "systematic audit" item: the audit is now the allowlist, enforced instead of pending.
6. **Every comparison table needs a floor, and a positive control wherever one can be built.**
   A floor is something that must fail (shuffled labels, a constant prediction); a positive control is something that must succeed (a planted signal).
   Without both, a table of real methods cannot separate "no method predicts this" from "nothing could have been detected here at all" — and that separation is what gives a null result its meaning.
   **Step** build, null. **Enforced by** `tests/test_project_rules.py::test_rule_06_comparison_tables_carry_a_floor_and_a_positive_control`, one case per discovered comparison table (`-m "step_build or step_null"`).
   **Edge case** a control can be present in the table and still prove nothing — a positive control substituted at scoring time rather than threaded into the fit target scores at chance, and is indistinguishable from a control that genuinely failed.
   **Edge-case test** `::test_rule_06_edge_positive_control_is_recovered_not_just_present` — enforced: in every table carrying both, the planted row must beat every floor row on the table's primary score; a positive control at chance fails.
7. **However much flexibility a model is allowed, every method in the same table gets it the same way.**
   Choose the number of components and the penalty strength by cross-validation, or fix identical values for every method — never one method tuned and another left at a value written into the script, even when that value came first and still looks reasonable.
   **Step** fit. **Enforced by** `tests/test_project_rules.py::test_rule_07_drivers_do_not_pin_capacity_to_a_literal`, one case per submission script (`-m step_fit`).
   **Edge case** a driver can avoid pinning capacity on the command line and still pin it in code.
   **Edge-case test** `::test_rule_07_edge_tables_record_selected_capacity` — the acceptance criterion for that column: it fails until every comparison table records what was selected.
8. **Every promoted result carries its `.provenance.json` record beside it** — commit, job id, the exact arguments used, and checksums of the inputs and the log — written from a checkout with no uncommitted changes, at the commit that produced the result.
   No record, not evidence.
   **Step** promote. **Enforced by** `tests/test_project_rules.py::test_rule_08_every_promoted_result_carries_a_complete_provenance_record` (`-m step_promote`).
   **Edge case** the clean-checkout half cannot be reconstructed after the fact — nothing recorded it.
   **Edge-case test** `::test_rule_08_edge_sidecars_record_clean_tree_and_producing_commit` — the acceptance criterion for recording both at promotion time, since neither can be reconstructed afterwards.
9. **Embargo is checked value by value, on what the cells contain, not on what the columns or the files are named.**
   A column named `line` is not the only place a patient identifier appears, and a table that merely lists another table's column names carries those identifiers as ordinary text.
   `data/release_manifest.yaml`'s three tiers cover every file: anything not explicitly classified stays embargoed.
   **Step** load, release. **Enforced by** `tests/test_project_rules.py::test_rule_09_embargo_gate_sees_an_identifier_inside_a_cell_value` (`-m "step_load or step_release"`).
   **Edge case** an identifier can sit in a cell of a column whose name looks innocuous, which a header-driven scan never opens.
   **Edge-case test** `::test_rule_09_embargo_gate_sees_an_identifier_inside_a_cell_value` — the gate must find an embargoed identifier inside an innocuously-named column's cell, not only in a column named like one.
10. **Reversing an analysis design, a method, or a data source needs a written decision, not just a commit.**
    If work undoes something a spec or plan established — a data source, a scoring rule, a control — write one paragraph in `docs/decisions/` saying what changed and why, however terse and however late.
    The CoderData → custom-loaders reversal (a full revert, `1bfb922`→`93dc76f`, merged as "restore-custom-loaders", 2026-06-16) left *zero* trace in any document here for over two months; the only record was git log and one person's memory.
    See `docs/decisions/2026-06-16-revert-coderdata-loaders.md` for the recovered entry.
    It is the sharpest example of the problem this file exists to prevent, and it must not repeat.
    **Step** document. **Enforced by** `tests/test_project_rules.py::test_rule_10_every_task_and_decision_is_named_in_the_spec_tree` (`-m step_document`).
    **Edge case** a reversal nobody writes down is undetectable mechanically, and no CI heuristic changes that.
    **Edge-case test** none, deliberately: an unwritten reversal is undetectable mechanically, and stays a human read of git log against `docs/decisions/`. The primary test catches only the next failure along — a document that exists but is indexed nowhere.

## Process

How work moves from question to promoted result — the lifecycle, the task folders, the compute boundary, the definition of done — is in [`docs/PROCESS.md`](PROCESS.md).
The short version: every task gets a folder under `docs/tasks/`, a place in the ladder above, and a state entry; every promoted number carries provenance; every rule above is enforced by a named test.
