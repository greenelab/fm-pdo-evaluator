# fm-pdo-evaluator — project spec

**As of** 2026-08-27.

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

## Three words this project uses precisely

| Word | Means | Lifecycle |
|---|---|---|
| **Rung** | A scientific question — one level of the ladder, adding one distribution shift. The experiment. | Answered, or not yet |
| **Step** | A stage of the evaluation pipeline that every rung passes through: load, build, restrict, split, fit, score, null, promote, release, document. The machinery. | Never "done"; rules attach here |
| **Task** | A unit of work with a spec, an owner and a definition of done. Touches one or more steps, serves one rung or is cross-cutting. | OPEN or CLOSED |
| **Tranche** | A versioned, content-hashed bundle of data from one source — the material a rung is run on, and what a leakage profile is keyed to. | Ingested once, then immutable |

A rung is *what we are asking*; a step is *where in the machinery*; a task is *what someone is doing about it*.
The project rules below bind to steps, which is why they hold for every rung without being restated per rung.

## Specific Implementation

---

### Rung 0 — is the target reproducible?

**Question** How much of a measured perturbation response is signal rather than assay noise?
**Adds** A reference for every other rung to read against.
**Measure** Split the replicate pool in half, correlate the two halves' per-(line, drug) expression deltas on the declared gene panel, Spearman-Brown corrected to full-data reliability, against a mismatched-pair null.
**Passing means** A ceiling significantly above the mismatched-pair null. A rung above it can then be reported as a fraction of what is achievable rather than as a fraction of a hypothetical 1.0.
**If it is low** A low ceiling does not stop the rungs above it; it rescales them. A score of 0.2 against a ceiling of 0.2 sits at the limit the assay supports, while the same 0.2 against a ceiling of 0.9 is a large shortfall — reporting a rung without its ceiling makes those two indistinguishable.

### Rung 1 — can a model predict an unseen cell line, in-distribution?


**Question** Given a cell line the model has never seen, can it predict that line's expression response to a drug?
**Adds** One boundary: a cell line the model has not seen — where *unseen* is what the tranche's leakage profile asserts, not what the split assumes. Same platform, readout and substrate.
**Measure** Correlation between predicted and measured delta on the held-out line, against a floor that must fail and a planted signal that must be recovered, reported as a fraction of rung 0's ceiling.
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
**Measure** The retained fraction of rung 1's in-platform score, per method, against a scrambled-line control.
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
**Adds** The final shifts: substrate, patients, and assay.
**Measure** The same interaction statistic as rung 4, fitted on cell lines and evaluated on the frozen organoid screen, reported as a fraction of both rung 4 and rung 0.
**Passing means** Either transfer holds at a measurable fraction, or it does not and the ladder below localises the loss to the substrate rather than leaving it unexplained.
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
   The old choice and the reason for changing it are dated and appended to the bottom of that task's `design.md` and `plan.md`, so the document a reader arrives at carries its own history.
   A reversal recorded only in git log is invisible to everyone who reads the documents, which is how a decision gets silently reinstated by the next person to touch it.
   **Step** document. **Enforced by** `tests/test_project_rules.py::test_rule_02_every_task_is_named_in_the_spec_tree` (`-m step_document`).
   **Edge case** the check reads the document, not the work: a reversal carried out in code while the task document is left untouched passes it. Where the reversal moves a promoted number, rule 1's checksum catches it instead; where it changes a method without yet changing a number, nothing does.
   **Edge-case test** `::test_rule_02_edge_non_additive_task_edits_carry_a_dated_entry` — a task document whose diff against the merge base deletes or rewrites existing lines must also gain a dated entry at its foot. Appending needs no entry; rewriting history does.

3. **The README stays in step with the documents it summarises.**
   Most readers open the README and nothing else, so a stale summary misinforms more people than a stale document does.
   It carries the project's question, the ladder, links to the three documents, and the current status — and it is updated in the change that moves any of them, not in a later tidy-up.
   **Step** document. **Enforced by** `tests/test_project_rules.py::test_rule_03_readme_links_to_the_project_documents` (`-m step_document`).
   **Edge case** a README can link correctly and still describe last month's ladder; no test reads prose for accuracy.
   **Edge-case test** `::test_rule_03_edge_readme_is_revisited_when_the_ladder_changes` — a change to the ladder in `docs/SPEC.md` requires `README.md` to have changed too. It proves the summary was revisited, not that it was revised well.

## Process

How work moves from question to promoted result — the lifecycle, the task folders, the compute boundary, the definition of done — is in [`docs/PROCESS.md`](PROCESS.md).
The short version: every task gets a folder under `docs/tasks/`, a place in the ladder above, and a state entry; every promoted number carries provenance; every rule above is enforced by a named test.
