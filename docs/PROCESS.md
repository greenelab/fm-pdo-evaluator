# How we operate

**As of** 2026-08-31.

This is the operating manual: the workflow, tooling, and collaboration mechanics for this project.
It answers "what do I do and in what order," not "what must every number respect" (that's `docs/SPEC.md`'s project rules) and not "what's true right now" (that's `docs/STATE.md`).
Read those two after this one, before starting work.

The practice is written down first so those tools get built to it, rather than the practice being reverse-engineered from whatever the tools happened to do.

## 0. Start here, every session

1. Read `docs/SPEC.md`'s spec tree — is there already an OPEN branch covering what you're about to do?
   Extend that task's `docs/tasks/<task-slug>/design.md` or explicitly supersede it; don't open a second task on the same subject.
2. Read `docs/STATE.md` — is the thing you're about to fix already fixed, already known broken with a reason, or already being worked in a session?
3. Check `git log --oneline -20` on the branch you're on.
   Multiple sessions (this machine and the HPC (in this case Alpine), sometimes concurrently) touch this repo; the two documents above are the coordination mechanism, not tribal memory of who's doing what.


## 1. The work lifecycle

**Brainstorm → Design spec → Implementation plan → Execute → Review → Verify → Audit → Summarise → Promote → Close out.**

**The lifecycle has exactly two human gates, and the run between them is autonomous.**

**Gate 1 — the design spec.** Brainstorming ends with the design presented and explicitly approved.
Until that approval nothing is edited and nothing is implemented — a revision request reopens the design, it does not approve it.
This is the gate that carries the weight: the design is where the measurement is decided, and every stage after it is checked against the design rather than against anyone's memory of the conversation.

**Gate 2 — the summary, and it comes before promotion.** Approval of the design opens the run, and plan → execute → review → verify → audit proceed without further asking, stopping at Summarise for human review; promotion happens after that review, not before it.
Each of those stages carries its own mechanical check — tests, controls, the verification battery, the audit's fresh reader — and pausing to ask between them substitutes an unchecked interruption for a checked stage.
Promotion is the one stage held back deliberately: it is what makes a number citable and what puts the evidence into the permanent git record, so it waits until a person has read the finding.
A promoted result that the summary review then changes costs a correction banner and a re-promotion (§6); the same change made before promotion costs nothing.
Three things still stop the run: a change to the design (that reopens gate 1 and needs approval before implementation continues); anything outside §2's standing permissions (force-push, `git reset --hard`, branch deletion, `upstream` or `main`); and a stage that cannot be resolved from the documents, which is reported rather than guessed at.
An audit's "fresh reader" is fresh with respect to the work, not necessarily human — a reader who did not write the code or the design satisfies it.

**Work is organized by task, not by date.**
Every task owns one directory, `docs/tasks/<task-slug>/`, named for what the task does with no date prefix.
The slug is the task's identity: it is what `SPEC.md`'s spec tree keys on and what `STATE.md` links to.
Dates live in git and in each document's own header line.

| Stage | Where it lives | Skill |
|---|---|---|
| Brainstorm | (conversation) | `superpowers:brainstorming` |
| Design spec | `docs/tasks/<task-slug>/design.md` | (write directly, or via brainstorming's output) |
| Implementation plan | `docs/tasks/<task-slug>/plan.md` — interfaces, invariants, ordered steps, and expected test outcomes; never embedded full source. Code exists once, in the repository: a plan that carries the code is a second copy, and a second copy drifts (five of rung 0's audit findings were plan-copy vs shipped-copy divergences) | `superpowers:writing-plans` |
| Execute | code changes | `superpowers:executing-plans` or `superpowers:subagent-driven-development`, with `superpowers:test-driven-development` and `superpowers:systematic-debugging` as needed |
| Review | `docs/tasks/<task-slug>/review.md` — each finding and what was done about it | `superpowers:requesting-code-review` / `receiving-code-review` |
| Verify | `docs/tasks/<task-slug>/verification.md` — the commands run, their output, and pointers to everything the run produced: tables, figures, logs | `superpowers:verification-before-completion` — evidence before "done," always |
| Audit | `docs/tasks/<task-slug>/audit.md` — the design read as numbered, checkable claims, each verified against the landed tree with evidence, every departure classified aligned / recorded / drift; then the fix wave's dispositions; then a **re-audit by a fresh reader** confirming every drift item fixed or recorded. The audit is not passed until the re-audit says so, and it passes **before the summary and before promotion**, so the finding a person reads and the number that becomes citable are both what the documents describe. It reads the run's artifacts in the working tree, before they are committed, and records the checksum of each one it read — the promotion commit's provenance records carry the same checksums, so an artifact that changed between the audit and the commit cannot pass silently. The procedure — both diff directions, the verdict vocabulary, what evidence each verdict needs, and the cap on re-auditing — is [`docs/audit.md`](audit.md), one standard for every task | (clause-by-clause read of `design.md` vs the tree; fresh readers for both passes) |
| Summarise | `docs/tasks/<task-slug>/summary.ipynb` — the finding in plain language, as a notebook the reviewer runs: the hypothesis, then **each measurement step in the order it happens, with its figures beside the table they were drawn from** — what the data looks like, what each control did to data with a known answer, what the statistical tests returned — then the conclusions and the scripts the task touched. Written for the lab member reviewing the pull request, in statistical and biological terms. Committed **without outputs**, like `verify.ipynb`, so the figures a reviewer sees are the ones their own execution produced. The step-by-step figure list is declared in that task's `design.md` before the run, not chosen afterwards from what looked good. **This is gate 2**: it runs against the working tree's artifacts, before those artifacts are committed | (write directly) |
| Promote | `results/<task-slug>/<result>.csv` + a `<result>.provenance.json` record beside it — only the subset the project will cite, copied out of the task folder and made citable. The record carries the result's own checksum, and promotion refuses when the two copies differ. Promotion and the evidence commit are one change: the artifacts the audit and the summary read enter the git record here, once, in their reviewed form | `scripts/promote_result.py` |
| Close out | branch merge/PR. The pull request opens as a **draft**; the task's documents — `verification.md` especially — get one more human review before it is marked ready. Its description opens with the **reading order** below and separates the **review surface** from **generated evidence**, both defined below | `superpowers:finishing-a-development-branch` |

**What a reviewer reads, and in what order.**
The pull request description opens by sending the reader to `summary.ipynb` — run it, and the finding arrives with the figures that support it, step by step, rather than as a claim to be taken on trust.
From there to `design.md`, which says what the measurement was supposed to be, so the reader can judge the finding against its intent rather than against their own expectations.
Then `audit.md`, which is where the two are checked against each other and against what actually landed.
Everything else — `verification.md`, `review.md`, `decisions.md`, the code — is reached from those three when a reader wants to go deeper, and the description says so rather than listing every file.
The order matters because the reverse order is what happens by default: a reviewer opening a pull request sees a diff, forms an impression from the code, and reads the finding last.

**What reaches GitHub, and when.**
Code and documents go up early; evidence goes up once, late.

**Code and documents must be pushed** — that is the only route to the cluster (§2: `git push` → `ralpine update` → `ralpine submit`, never a direct copy).
A task branch therefore accumulates commits of `src/`, `scripts/`, `tests/` and the task's documents from the moment implementation starts, and that is fine: it is the review surface, it is small, and its history is what a reviewer reads to see how the work developed.

**Everything a run produces stays out of the git record until promotion.**
Results land on Alpine scratch, are pulled into the local working tree, and stay uncommitted through Verify, Audit and Summarise.
The verification battery, the audit and `summary.ipynb` all run against those local files.
Only after gate 2 does the evidence enter git, in the promotion change, in the form a person has already read.
This is not tidiness. Git stores every version of a file whole and forever, so a result table committed, corrected and recommitted three times costs three copies of itself in the repository permanently, and no later cleanup removes them; committing once after review costs one.

**The exposure this accepts, stated plainly.**
Between the run and promotion, the evidence exists on Alpine scratch (purgeable) and in one local working tree, and nowhere else.
Losing both costs a re-run, not the result: everything is reproducible from the pinned tranche and the committed code, which is the reproduction chain the tranche records exist to guarantee, and the build cache puts a rerun at about forty minutes rather than a re-derivation from scratch.
What it must never cost is a claim without its artifact — hence the audit's checksum rule above, and hence promotion and the evidence commit being one change rather than two.

**Review surface versus generated evidence.**
The **review surface** is what a person is being asked to judge: the task's documents (`design.md`, `plan.md`, `review.md`, `verification.md`, `audit.md`, `decisions.md`), both notebooks, and the code — `src/`, `scripts/`, `tests/` — plus any change to a project document.
**Generated evidence** is everything a run produced and committed so a reviewer can re-derive the numbers without cluster access: result tables, provenance records, parameter files, figures, cluster logs, tranche manifests.
It is committed to be *checked*, mechanically, by the verification battery — not read line by line, and a reviewer is never asked to eyeball a table of ten thousand correlations.
Generated paths are marked `linguist-generated` in `.gitattributes` so they collapse in the diff view; when a new kind of generated artifact appears, it is added there in the same change that first produces it, or the next pull request buries its review surface under it.
Rung 0's first pull request added roughly 26,000 lines of which roughly 17,500 were generated tables, and the lab review read that as review burden rather than as evidence.

**The lifecycle runs per task and converges; it is not re-entered per run.**
Work added after a stage has passed ripples downstream incrementally — a scoped review of the new commits, an amended `verification.md`, an audit delta over the new surface with the same fresh-reader confirmation, an amended `summary.ipynb` — never a restart of the whole cycle.
Promotion re-triggers only when a promoted claim changes, and that is §6's correction path (re-run, re-promote with fresh provenance, correction banner), not an amendment.
This holds both failure modes off: infinite full-cycle loops on every late addition, and late additions slipping in unreviewed after the gates.

**The ripple is capped.**
Documentary findings batch into a single fix wave; the re-audit re-checks only the items verdicted drift; there is no confirmation pass of a confirmation pass.
Mechanical claims — counts, tallies, whether a cross-reference resolves — are checked by script, never by a fresh reader: agents for judgment, scripts for arithmetic.
Rung 0's audit deltas kept finding documentary defects introduced by the previous fix wave, and one pass was spent reconciling the audit's own clause counting; the cap is what stops prose maintenance from generating its own findings without bound.

Three rules keep the task folder from becoming an orphan:

1. **Create the folder when you write the design spec**, and in that same change add the task's branch to `SPEC.md`'s spec tree, under its rung or under cross-cutting.
   A task folder outside the tree is invisible to the next session, which is the failure mode this structure exists to prevent.
2. **One task, one folder.**
   Extending existing work means editing that task's `design.md`, not starting `<something>-v2`.
   If the scope genuinely changes into different work, open a new task and supersede the old one explicitly (§5).
3. **`STATE.md` links back to the folder.**
   Every status entry names its task's `design.md`, the commits that implemented it, and the outputs it produced (§5) — that link is what makes a number traceable to the intent behind it.

**Decision lineage lives in `decisions.md`, at two levels.**
A dated entry for a change to a task's design or plan goes to that task's `docs/tasks/<task-slug>/decisions.md`, labeled by the document it amends; a dated entry for a change to a project document (this manual, `docs/SPEC.md`) goes to `docs/decisions.md`.
The amended document itself carries only the current position plus a one-line pointer to its lineage — history appended at the foot of a document readers open for current guidance is how the foot outgrows the document.
Every task folder gets a `decisions.md`; SPEC project rule 2 is the binding form.
Decisions that cut across tasks are linked to the `design.md` or other evidence in `docs/tasks/<other-task-slug>/` that drove the change relevant to the current task.


**Promotion is the line between "we ran something" and "it's evidence."**
A result without a `.provenance.json` provenance record is not citable in a write-up, per this project's own standing rule — see project rule 1 in `SPEC.md`.

## 2. Compute: local vs. Alpine

- **Local (this machine)**: interactive development, doc/spec work.
  Run `uv run pytest` here before every push.
- **Alpine**: anything needing real compute.
  Access is exclusively through `./scripts/alpine/ralpine` — never raw `ssh`/`scp`.
  That boundary is enforced in version-controlled script code (`READ_ONLY` allowlist + `reject_metacharacters`), not in a permission string, specifically so it's reviewable.
- **Deployment is git-only**: `git push` → `ralpine update` (git pull --ff-only on the Alpine checkout) → `ralpine submit`.
  Never copy a file to Alpine directly — if Alpine needs a file that isn't code, that's a sign it belongs in `data/` with a download/build script, not a one-off transfer.
- **Chained jobs**: pass `--dependency=afterok:<jobid>` etc. to `ralpine submit`; the wrapper reorders args correctly, but *verify* with `ralpine jobinfo <id> | grep Dependency` after submitting — a chain that silently drops its dependency runs immediately and out of order.
- **Parallelize a scan-bound job as a job array, never as a loop in one process.**
  When the cost of a job is reading a large table (rung 0's four billion rows), find a key that is part of every group the job aggregates over — for rung 0, the gene — and split the work by a hash of that key into slices, so every group lands in exactly one slice and per-slice results concatenate or add back to the whole exactly.
  Run the slices as one Slurm job array (`#SBATCH --array=0-N`), one task per slice, each with the memory a slice needs rather than the memory the whole needs; then one combine job, submitted with `--dependency=afterok:<array id>`, that reads the cached slices and does only what needs the whole.
  Wall clock is one scan instead of N in a row; a task that finds its slice already cached skips, so a failed index is resubmitted alone (`--array=<index>`); and a slice count is a parameter both the array and the combine must agree on, so it lives in one place (`scripts/alpine/rung0_env.sh`).
  Prove the slicing is arithmetic with a test that runs the same fixture in one slice and in several and compares them exactly (`test_partitioned_frame_equals_one_pass`, `test_partitioned_noise_equals_one_pass`), and test the staged path against the one-process path on a fixture before submitting (`test_the_staged_run_produces_the_same_artifacts_as_one_process`).
  More slices than nodes the scheduler will give at once buys nothing: every task still reads the whole table, so past that point extra slices only lower per-task memory.
- **Verify inputs exist before submitting.**
  `ralpine ls`/`ralpine du` the files a job needs first.
- **Pull the job's log into `results/<task-slug>/` when it finishes.**
  A log left on the cluster is unreadable to anyone reviewing the pull request, and is gone when the scratch space is cleared.
  It is the record of what the run actually did — the resolved arguments, the warnings, the wall time — so it belongs beside the result it produced.
- **Standing permissions** : `git commit`/`git push` to `origin` on the current task's branch, effective once that task's design is approved; `ralpine submit`/`cancel`/`update` — all without per-action confirmation.
  `submit` is the working default — an agent stages the run and submits it.
  `cancel` covers this project's own jobs that are stale (superseded by a newer submission, or running code that a later commit has already fixed) or running far past their expected wall time; never another user's job, and never a whole array range, which `ralpine cancel` refuses structurally by taking exactly one numeric id.
  Still confirm before force-push, `git reset --hard`, branch deletion, or anything touching `upstream` or `main`.
  A standing grant is scoped to what it says — if a new branch or a materially different kind of action comes up, ask.

## 3. Verification discipline

- **Test on synthetic data before spending cluster time.**
  Build a small synthetic fixture that exercises the real code path (not a reimplementation of its logic), run it locally, confirm the shape and sign of the result look right, *then* submit to Alpine.
- **Every measurement step carries a positive and a negative control** — project rule 4 in `docs/SPEC.md` is the binding form.
  The positive control plants a known signal and requires the real, shipped code to recover it; the negative control feeds signal-free or mismatched data and requires null.
  Both are declared per step in the task's `design.md` and implemented as known-answer tests that import the real function.
- **State the detection power of every reported comparison.**
  Every promoted comparison carries its minimum detectable effect at the declared α and power, computed from the same null bootstrap as its p-value, and positive-control plants are placed relative to it.
  A null result without its MDE cannot be told apart from an underpowered experiment — the distinction rung 5's small cohort will turn on.
- **Run the full local suite before every push**: `uv run pytest -q`.
- **Run the lint and type gates over tracked files, not over hand-picked directories.**
  `git ls-files '*.py' '*.ipynb' | xargs uv run ruff check`, the same with `ruff format --check`, then `uv run pyright`.
  Continuous integration runs `ruff check .` on a clean checkout, which lints notebook code cells too; a local gate scoped to source directories cannot see a violation in a committed notebook, and the branch went red for a day on exactly that gap.
  Selecting tracked files keeps the untracked working-tree strays out (the reason the scoping existed) without hiding anything continuous integration will see.
  **Run the gates after staging, not before.** `git ls-files` lists what is tracked *now*, so a file added in the same change is invisible to a gate run before `git add` — which is how a notebook reached a commit unformatted while the local check reported every file clean. Continuous integration sees the committed tree and has no such blind spot, so the gate that matters is the one run against what is about to be committed.
- **Run the project-rule tests for what you touched.**
  `uv run pytest tests/test_project_rules.py` every time — they may catch a violation the task never intended — plus `-m` for the steps named in the task's `design.md` header, e.g. `uv run pytest -m "step_split or step_score"`.
  The step-to-rule-to-test mapping sits under each rule in `docs/SPEC.md`, and each step's marker is registered in `pyproject.toml` alongside the test that uses it, so the selection is the same whichever rung or dataset the task is about.
  An exemption in that run is a recorded gap, not a pass: if the task closed one, the test starts passing unexpectedly and the exemption comes out in the same change.
- **`verification-before-completion` applies to claims, not just code**: before saying something is fixed, promoted, or done, show the command and its output.
  "Should work now" is not done.
- **A verification check that is cheap to compute is worth running.**
  When an assumption or a claim can be tested directly for minutes of compute, run the test rather than argue the assumption — the argument convinces, the check settles.
- **Committed evidence is sized to the claim it supports, not to what the run produced.**
  A table is committed when a promoted claim is re-derived from it (the per-condition correlations, the null draws, the screen composition); an illustrative artifact is committed compressed, and carries both its own recomputable statistic and the authoritative full-data one.
  Compress before you subsample: gzip cost a quarter of the plain-text size here and stays exact, while a 2,000-gene subsample of the example profiles put enough sampling error on small correlations to reorder them — cheaper bytes are not worth a figure whose own numbers contradict its caption.
  Git stores every version of a file whole and forever, so a large artifact regenerated each run costs its full size on every rerun; anything genuinely large stays on the cluster, pinned by checksum in the provenance record, and the repository carries what a reviewer needs to check the claim.
- **Every task ships an executable verification entry point.**
  A script that recomputes each promoted claim from the artifacts alone — hashes against the provenance record, statistics from the tables, the summary's evidence table checked against the artifacts mechanically — printing claim / recomputed / pass-fail, runnable on a laptop in about a minute; and a notebook, committed **without outputs**, whose cells recompute the same claims **inline and self-contained** (standard-library hashing, direct table reads, explicit arithmetic — nothing imported from the project's own code), so the reviewer reads exactly what is computed and watches it pass.
  Before promotion it reads those artifacts in the working tree at the paths they will be committed to; after promotion the same script over the same paths is what continuous integration runs on a clean checkout, and the audit's recorded checksums are what tie the two runs to the same bytes.
  A notebook that only calls the script relocates the trust instead of discharging it; the script is the notebook's final cross-check cell and the continuous-integration form (a test runs the same checks), not its body.
  Trust comes from re-derivation the reader performs, not from narrative the writer asserts; a number recomputed at read time cannot drift the way a number transcribed across documents can.
- **Two notebooks, two jobs, no overlap.**
  `verify.ipynb` *proves*: it recomputes each promoted claim inline from the committed artifacts, imports nothing from the project's own code, and carries no figures.
  `summary.ipynb` *explains*: it walks the reviewer through the measurement step by step with the figures that show what the data looks like and what each control did, reading the same committed tables.
  A number that appears in both is computed in `verify.ipynb` and displayed in `summary.ipynb`, never recomputed twice by two methods — two derivations of one number is how a document starts disagreeing with itself.
- **Every figure is drawn from a committed table, and shows its control where one exists.**
  A figure whose values live only inside a run cannot be checked, and a figure of real data alone shows what the screen looks like without showing whether the machinery reads it correctly — so the planted-answer panel sits beside the real one, on shared axes.
  Figures are generated by the run and never drawn by hand (Greene Lab reproducibility standard), and which figures a task produces is declared per step in its `design.md` before the run, so the reviewer sees the evidence the design promised rather than the subset that looked best afterwards.
- **State an assumption with its quantitative exposure, never as a bare hedge.**
  Say what value would have to obtain for the conclusion to change and whether that value is plausible — "the draws are not independent" alone tells a reader nothing; "the dependence would need to inflate the null variance 3,700-fold to lose significance" tells them everything.

## 4. Git and collaboration

- **One task, one branch — and one at a time.**
  Every task gets its own branch, cut from the trunk (`project-docs` until it merges; `main` after that) and named for the task slug, so the branch, the task folder, and the spec-tree entry share one identity.
  Process changes a task motivates are made on the task's branch and land with it, not slipped onto the trunk.
  Finish the work, open the pull request, merge it, then branch again: branching before the current piece lands puts two half-finished versions of the same documents on disk, and every question after that — which number is current, which file is real — has more than one answer.
- **The branch carries every script in its result's provenance chain.**
  Data downloads, input builds, model builds: if a script produced anything the task's result depends on, it is on the task's branch before the result is promoted.
  A result whose producing scripts live elsewhere — an old worktree, an uncommitted file on the cluster — cannot be regenerated from the branch that claims it.
- **Commit messages state the root cause, not just the change**, written so a future reader can tell what happened and why without re-deriving it from the diff.
  Do not put result numbers in a commit message — reference the output and the document that reports it, so there is one place a number can be corrected.
- **Push to `origin`** (the `lagillenwater` fork).
  **Never `upstream`** (`greenelab`) without being explicitly told.
  Greene Lab standard beyond this fork: fork-per-contributor, pull requests only, minimum one lab-member approval, one functional area per pull request.

## 5. Documentation discipline

The rules these serve are in [`docs/SPEC.md`](SPEC.md); this is what they require in practice:

- A new task = a `docs/tasks/<task-slug>/design.md` **and** a spec-tree branch in `SPEC.md` pointing at it, in the same change.
- Reversing an architecture, method, or data source are updated in the task-specific `design.md` and `plan.md` docs, with the old choices and reasons for changing dated and appended to the bottom of the documents. No exceptions.
- A change to what the project asks, how work is done, or where a rung stands = `README.md` updated in the same change.
  The README is the only document most readers open, so it is the one that goes stale first and the one that misleads most when it does. Update the summary, not just the document it summarises.
- A number change = `docs/STATE.md` updated in the same change that produced it, and that entry carries its three links: **spec** (the task's `design.md`), **code** (the commits or files that changed), **outputs** (the promoted artifacts + their provenance record).
  A status claim with no link to the spec that motivated it, the code that produced it, and the artifact it came from is a claim a future session has to re-derive from git log.
- **Landed documents reference only landed work.**
  A rung above the current one is referenced in the sense of `docs/SPEC.md` — what the ladder will ask — never through unlanded code or documents a reviewer of this repository cannot open.
  Work from the archived development lineage is cited only by its archive branch or archive file, labeled as archive.
  A forward reference a reviewer cannot follow reads as missing work, and a worktree reference reads as work that exists when, for the repository, it does not.
  Files on the cluster are cited by repository-relative path plus "on Alpine" (scratch locations in `$USER` form) — an absolute site path embeds a username and a mount layout no reviewer can use.
- **Every acronym is spelled out at its first use, in every report document.**
  A task's design, verification, review, audit, and summary documents are read by people who did not write them; an undefined acronym is a claim the reader cannot check.
  Common file-format and tooling abbreviations are not exempt.

## 6. Definition of done

- **A bug fix**: a real test (imports the actual function) passes locally; if it touches a promoted number, that result is re-run, re-promoted with a fresh provenance record, and the document that reported the old number gets a correction banner — never a silent edit.
- **A new capability**: went through brainstorm → spec → plan in one `docs/tasks/<task-slug>/` folder, indexed in `SPEC.md`, with `STATE.md`'s entry linking spec, code, and outputs. The design was approved at gate 1 and the run proceeded to the summary without further approvals (§1).
- **A claim of "complete"**: backed by command output you actually ran, not by what the code should do.
- **Any task**: the README reflects it. A rung that landed, a rule that arrived, a status that moved — if a reader would learn it from the three documents, the README says it too.
- **A task heading to its pull request**: in this order — the drift audit passed after verification (design read clause by clause against the landed tree, every departure aligned or recorded through a fix wave, a re-audit by a fresh reader confirming it, all in `audit.md` and following [`docs/audit.md`](audit.md)); then the summary notebook written and read by a person (`summary.ipynb`: hypothesis, then every step with its figures and controls, conclusions, scripts touched); then promotion, which is also the change that commits the run's evidence. The pull request opens as a draft after all three, its description sending the reader to the summary first, then the design, then the audit, and separating review surface from generated evidence (§1). It is marked ready only after a final human review of the task documents.

## 7. Failure modes this process exists to prevent

Written down because each has already cost this project real time, and because the lesson outlives the specific confusion.

- **A newer, more authoritative-looking document can itself be the source of drift.**
  Recency is not accuracy: a new spec must be checked against the code and the specs it supersedes, not just written down.
- **Renaming a concept without a forward pointer leaves every prior document silently ambiguous.**
  A code-level rename with a working alias reads as complete while every document using the old name keeps meaning something slightly different.
  A rename is a two-line fix in the documents that used the old name; do it at rename time, not retroactively.
- **A parallel implementation that bypasses an abstraction silently drops that abstraction's guarantees.**
  When a driver is written to call a library directly because the registry path is inconvenient, every guarantee the registry provided — leakage filtering, shared partitions — quietly stops applying, and no output says so.
  If a guarantee matters, it belongs where it cannot be routed around.
- **A reversal with no decision entry is invisible to everyone but git log.**
  Reverting an architecture, a method or a data source and writing nothing down leaves the next reader to rediscover it from commit archaeology, or to reinstate what was deliberately removed.

## Changes

This manual's dated change lineage lives in [`docs/decisions.md`](decisions.md) (moved 2026-08-31).
