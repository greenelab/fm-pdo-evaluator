# fm-pm-evaluator — project spec

**As of** 2026-08-31.

## The question

Foundation models trained on cell-line transcriptomic data claim generalization to out-of-distribution tasks.
This repo builds the infrastructure to find the boundary of those claims.
This repo will test whether model prediction survives the move to an unseen cell line, across platforms, in a new modality, and on prospective experimental data.

## The generalization boundary ladder

Running a cell-line-trained model directly on the hardest target yields one number and no diagnosis.
When that number is poor — and it usually is — seven distinct explanations are consistent with it.
Each names a boundary the prediction might not cross, and each is the question one rung settles:

0. the target itself is not reproducible, so nothing could have been predicted;
1. the model cannot predict an unseen cell line, even within its own training distribution;
2. the resolution changed — the model reads single cells and the sample was measured in bulk;
3. the measurement platform changed;
4. the readout changed, from expression response to viability;
5. the substrate changed — patient-derived organoids are not cell lines;
6. the result held only because the outcome was already known when the analytic choices were made.

Seven different conclusions demanding seven different responses, and one end-to-end number cannot separate them.

The evaluation is built as a **ladder**: each rung adds exactly one boundary to the rung below it, and every rung is scored against the reproducibility ceiling.
Where the score falls off localises the failure to a specific boundary.

## Five words this project uses precisely

| Word | Means | Lifecycle |
|---|---|---|
| **Rung** | A scientific question — one level of the ladder, adding one distribution shift. The experiment. | Answered, or not yet |
| **Step** | A stage of the evaluation pipeline that every rung passes through: load, build, restrict, split, fit, score, null, promote, release, document. The machinery. | Never "done"; rules attach here |
| **Task** | A unit of work with a spec, an owner and a definition of done. Touches one or more steps, serves one rung or is cross-cutting. | OPEN or CLOSED |
| **Tranche** | A versioned, content-hashed bundle of data from one source — the material a rung is run on, and what a leakage profile is keyed to. | Ingested once, then immutable |
| **Frame** | The invariants a cross-rung ratio holds fixed: the tranches (by content hash), the declared statistic, the scoring unit and how dose enters it, and the weighting over units. A fraction-of-ceiling is valid only when numerator and denominator share a frame *and* a panel — the denominator being rung 0's ceiling restricted to the numerator's genes, drugs and dose-conditions, under the same weighting. | Versioned by its hashes; a new tranche opens a new frame; a new panel adds a restriction, not a frame; a change to the unit, the dose handling or the weighting is a new statistic, and so a new frame |

A rung is *what we are asking*; a step is *where in the machinery*; a task is *what someone is doing about it*; a tranche is *what it runs on*; a frame is *what a ratio holds fixed*.
The project rules below bind to steps, which is why they hold for every rung without being restated per rung.

Scores from different frames are different quantities.
Bringing in a new tranche — new cell lines, a new screen — does not extend the current frame; it opens a new one, and every cross-rung comparison is recomputed inside it.
That recompute is a re-run of parameterized pipelines over new inputs, not a redesign, and it does not invalidate the old frame: promoted results carry the tranche hashes that say which frame they belong to, and remain citable claims about it.
Panels are not frame elements. Rung 0 measures the ceiling at the assay's full extent — every gene and every splittable drug — and a rung that scores a subset reads against a **restriction** of that ceiling to its own subset, declared before it scores, promoted under the rung that needs it, and listed in `docs/STATE.md` with its panel hashes and motivating rung. A restriction never changes the value it restricts; it is a different question, never a revision.
The one asymmetry is rung 6 — a registered prospective prediction is frozen in its frame by definition, so a frame change after registration is recorded, and new predictions are registered in the new frame going forward.

## Specific Implementation

---

### Rung 0 — is the target reproducible?

**Question** How much of a measured perturbation response is signal rather than assay noise?
**Adds** A reference for every other rung to read against.
**Measure** Split each (line, drug, dose) condition's replicate plates in two and correlate the two groups' per-gene expression deltas — per-condition Pearson, aggregated as the mean over conditions, the one statistic every delta rung reports — Spearman-Brown corrected to full-data reliability, against stratified mismatched-pair nulls and alongside a planted-signal positive control. Dose is part of the condition, never pooled: on this screen dose is very nearly confounded with plate, so a split that pools dose compares one dose against another. The screen replicated its doses unevenly, so the rung reports the mean over every replicated (line, drug, dose) triple with equal weight AND the same statistic at each dose level and with each (line, drug) pair weighted once, promotes the per-triple table those are aggregates of, and declares in its task's decision lineage which aggregate is the ceiling. A later rung reads against a restriction of that aggregate to its own triples, under the same weighting. At the assay's full extent, every drug with plates enough to split, over two gene sets: every gene the table carries, and the genes each condition's first plate group calls differentially expressed — selected from that group alone, never from the group the correlation is scored against. The two are different questions, since across all genes the correlation is dominated by genes that did not respond. Promoted beside them: a decomposition of each delta's variance into its between-plate part and the within-plate standard error the screen's own differential-expression statistics carry, which says whether replicate noise is plate effects or cell sampling — a finding about the assay in its own right, and the thing that decides whether an unreplicated condition can ever carry a ceiling. Rung 0 measures at full extent only; a later rung that scores a subset declares and computes its own restriction of these ceilings, per the frame rules above.
**Passing means** Both ceilings significantly above the mismatched-pair null — the same-drug null holding dose fixed as the condition does. A rung above them can then be reported as a fraction of what is achievable rather than as a fraction of a hypothetical 1.0.
**If it is low** A low ceiling does not stop the rungs above it; it rescales them. A score of 0.2 against a ceiling of 0.2 sits at the limit the assay supports, while the same 0.2 against a ceiling of 0.9 is a large shortfall — reporting a rung without its ceiling makes those two indistinguishable.
**Tasks** [docs/tasks/rung0-assay-reliability/design.md](tasks/rung0-assay-reliability/design.md) — OPEN. Supersedes the unmerged branch `rung0-replicate-ceiling`, whose gene and drug panels could not be justified from anything a reader at rung 0 can open.

### Rung 1 — can a model predict an unseen cell line, in-distribution?


**Question** Given a cell line the model has never seen, can it predict that line's expression response to a drug?
**Adds** One boundary: a cell line the model has not seen — where *unseen* is what the tranche's leakage profile asserts, not what the split assumes. Same platform, readout and substrate.
**Measure** Correlation between predicted and measured delta on the held-out line, scored with rung 0's declared statistic — per-pair Pearson, mean over pairs, so the fraction is a ratio of like quantities — against a floor that must fail and a planted signal that must be recovered, reported as a fraction of rung 0's ceiling restricted to the genes and drugs this rung scores — whichever of rung 0's two ceilings matches the genes it scores — a restriction this rung declares before scoring.
**Passing means** The method beats its floor and recovers a stated fraction of the ceiling — the minimum competence claim, in the easiest setting the ladder offers.
**How it contextualises the rest** A model that fails here fails in-distribution, so a shortfall at any rung above is not evidence about organoids, platforms or readouts — the same weakness is already present with none of those boundaries crossed.

### Rung 2 — can a bulk sample be read by a single-cell model?

**Question** The model reads individual cells; every rung above reads bulk on at least one side. Does a bulk profile, converted into a population the model can read, arrive where the same material's real single cells would?
**Adds** One boundary on top of rung 1: resolution. Same biology, same platform, a population average rather than cell by cell.
**Measure** Agreement between a synthesised population and the real single cells of the same material, calibrated on paired benchmarks that profile the same cultures both ways, then checked on our own lines. The bridge and the benchmarks are named in the rung's spec.
**Passing means** The synthesised population lands near the real cells by whatever representation the rungs above consume, clearing a mismatched-line null.
**How it contextualises the rest** Every rung above carries this conversion whether or not it has been priced. Measured here, the resolution cost separates from rung 3's platform cost; unmeasured, the two are reported together and neither is attributable.

### Rung 3 — what does crossing a measurement platform cost?


**Question** How much predictive signal survives when the mapping is fitted on one transcriptomic platform and tested on another?
**Adds** One boundary on top of rung 2: the platform (fit on L1000, test on Tahoe), with the resolution cost already priced below.
**Measure** The retained fraction of rung 2's score, per method, against a scrambled-line control. Measuring against rung 1 instead would fold the resolution cost into the platform number, which is the confusion rung 2 exists to remove.
**Passing means** Retention clearly separable from the scrambled control — the platform shift costs something quantifiable rather than everything.
**How it contextualises the rest** The organoid comparison contains a platform change bundled with the substrate change. This rung prices the platform part on its own, so rung 5's shortfall is not silently attributed to biology.

### Rung 4 — does the representation predict drug response, not just expression?


**Question** Does a representation carry information about how much a drug *kills* a line, beyond what the drug does on average?
**Adds** One shift: the readout changes from expression delta to viability (dose-response AUC).
**Measure** Drug×line interaction, residualising each drug's mean so that only line-specific ranking counts, corrected across all declared variants, reported against the screen-agreement ceiling between independent screens of the same lines.
**Passing means** Interaction significantly above zero after correction, at a stated fraction of the screen-agreement ceiling.
**How it contextualises the rest** Viability is the endpoint the organoid rung uses. A representation that predicts expression but not viability would fail rung 5 for reasons that have nothing to do with organoids.

### Rung 5 — does it transfer to patient-derived organoids?

**Question** Does a model trained on cell-line drug response predict drug response in patient-derived tumour organoids?
**Adds** Three boundaries at once — substrate, patients, and assay — and this is the one rung that cannot hold to one. No available data separates them: an organoid is patient-derived by definition, and it is screened on a different assay than a cell-line panel. The ladder's discipline breaks here, deliberately and visibly, rather than by an unstated bundling.
**Measure** The same interaction statistic as rung 4, fitted on cell lines and evaluated on the frozen organoid screen, reported as a fraction of both rung 4 and rung 0.
**Passing means** Either transfer holds at a measurable fraction, or it does not — and because the rungs below have priced resolution, platform and readout, what remains here is the substrate, the patients and the assay together, not an unexplained shortfall. Attributing it further needs data that does not exist: cell lines screened on the organoid assay would separate assay from substrate.
**Constraint** A frozen holdout under embargo, and rung 2's bridge is calibrated on cell lines but assumed here. The rung's spec states that assumption, what survives it, and the residual-structure alarm that can break it.

---

### Rung 6 — does it hold when the prediction comes first?

**Question** Does the transfer survive when predictions are filed before the organoids are screened?
**Adds** Time. Nothing about the data or the model changes; the outcome simply does not exist when the prediction is made.
**Measure** The same interaction statistic as rung 5, on organoids screened after a registered prediction, reported alongside rung 5's retrospective score.
**Passing means** The prospective score is consistent with the retrospective one. The gap between them prices what retrospective evaluation was worth.
**Constraint** A prediction counts only if it is registered before the screen is run: written down, committed, and dated in this repository.

---

### Cross-cutting work

Some work belongs to no single rung: the shared evaluation library every rung scores through, the fairness machinery that makes methods comparable, the release gate that keeps patient identifiers out of anything published, and the enforcement that turns the rules below from advice into tests.
Each is specified alongside the rung that first needs it, rather than up front, so that a piece of shared machinery is designed against a real requirement.

## Project rules — what every task must satisfy

A rule, once here, is permanent: **amend this section deliberately when a genuine exception is needed — silently violating one is unacceptable.**

Each rule names the pipeline step it governs and the test that enforces it.
**The step is the unit, not the task**: any task touching that step inherits its tests, whichever rung or dataset the task is about.
A task runs `uv run pytest tests/test_project_rules.py` in full — some of these check the repository or its artifacts rather than one function, and catch a violation the task never intended — and adds `-m` for the steps it touched.
Each step's marker is registered in `pyproject.toml` alongside the test that uses it.
The rules stay generic: which artifact violates one today is current state, and belongs in `docs/STATE.md` and the task that owns the fix, not here.

1. **Every promoted result carries its provenance record beside it in `results/<task-slug>/`.**
   The record is `fmharness.schema.PromotedResult`, so the format is defined once, in code, and validated rather than trusted: the producing commit, seed and determinism flag come from the same `EnvironmentSnapshot` a prediction carries, alongside the arguments, the input and log checksums, the job id, and the artifact's own checksum.
   Three of its fields cannot be recovered later and so must be written at promotion: whether the tree was clean, which commit produced the result, and the result's checksum. That last one is what keeps the artifact from being edited underneath the claim — promotion refuses when the task-folder copy and the promoted copy differ.
   No record, not evidence.
   **Step** promote. **Enforced by** `tests/test_project_rules.py::test_rule_01_every_promoted_result_carries_a_complete_provenance_record` (`-m step_promote`).
   **Edge case** the record can be complete and still untrue: nothing checks that the commit named is the one that ran, or that the tree really was clean. Validation catches an absent field, not a wrong one.
   **Edge-case test** `::test_rule_01_edge_promoted_records_validate_against_the_schema` — every record must parse as `PromotedResult`, so a second provenance format cannot appear alongside the first without failing.

2. **Reversing an analysis design, a method, or a data source is written into the task's own documents, not left to the commit.**
   The old choice and the reason for changing it are dated and appended to that task's `decisions.md`, labeled by the document they amend (project-level decisions go to `docs/decisions.md`), so a task's history sits in one findable place while `design.md` and `plan.md` carry the current position.
   A reversal recorded only in git log is invisible to everyone who reads the documents, which is how a decision gets silently reinstated by the next person to touch it.
   **Step** document. **Enforced by** `tests/test_project_rules.py::test_rule_02_every_task_is_named_in_the_spec_tree` (`-m step_document`).
   **Edge case** the check reads the document, not the work: a reversal carried out in code while the task document is left untouched passes it. Where the reversal moves a promoted number, rule 1's checksum catches it instead; where it changes a method without yet changing a number, nothing does.
   **Edge-case test** `::test_rule_02_edge_non_additive_task_edits_carry_a_dated_entry` — a task document (`design.md`, `plan.md`, or `decisions.md` itself) whose diff against the merge base deletes or rewrites existing lines must also gain a dated entry, in the document itself or in the task's `decisions.md`. Appending needs no entry; rewriting history does.

3. **The README stays in step with the documents it summarises.**
   Most readers open the README and nothing else, so a stale summary misinforms more people than a stale document does.
   It carries the project's question, the ladder, links to the three documents, and the current status — and it is updated in the change that moves any of them, not in a later tidy-up.
   **Step** document. **Enforced by** `tests/test_project_rules.py::test_rule_03_readme_links_to_the_project_documents` (`-m step_document`).
   **Edge case** a README can link correctly and still describe last month's ladder; no test reads prose for accuracy.
   **Edge-case test** `::test_rule_03_edge_readme_is_revisited_when_the_ladder_changes` — a change to the ladder in `docs/SPEC.md` requires `README.md` to have changed too. It proves the summary was revisited, not that it was revised well.

4. **Every measurement step a task touches carries a positive and a negative control, and every promoted comparison reports its detection power.**
   For each pipeline step that produces or transforms a measurement — load, build, restrict, split, fit, score, null; promote, release and document are bookkeeping — the task's `design.md` declares a positive control (plant a known signal, require the real code to recover it) and a negative control (feed signal-free or mismatched data, require null), implemented as known-answer tests that import the shipped functions.
   Controls are placed relative to the measurement's minimum detectable effect, and every promoted comparison reports that MDE at the declared α and power beside its p-value — a null with no MDE cannot be told apart from an underpowered experiment.
   **Step** score, null (declared at design time for every measurement step the task touches). **Enforced by** `tests/test_project_rules.py::test_rule_04_every_task_declares_controls_for_its_measurement_steps` (`-m "step_score or step_null"`).
   **Edge case** the scan reads declarations, not implementations: a Controls section can describe a control no test implements, and a weak plant proves little. The known-answer tests are the implemented half, and they land with the code they validate.
   **Edge-case test** `::test_rule_04_edge_promoted_tasks_have_known_answer_tests` — a promoted result may not exist while the repository carries no known-answer-marked test; skips until the first promotion.

## Process

How work moves from question to promoted result — the lifecycle, the task folders, the compute boundary, the definition of done — is in [`docs/PROCESS.md`](PROCESS.md).
The short version: every task gets a folder under `docs/tasks/`, its own branch cut from the trunk, a place in the ladder above, and a state entry; every task passes a drift audit of its design against what actually landed **before its result is promoted**, to the standard in [`docs/audit.md`](audit.md); every promoted number carries provenance; every rule above is enforced by a named test; and every task carries a plain-language summary notebook of hypothesis, evidence and conclusions, walking the reviewer through each step with its figures, its pull request opening as a draft for one final review round.

## Changes

This document's dated change lineage lives in [`docs/decisions.md`](decisions.md) (moved 2026-08-31).
