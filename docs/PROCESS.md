# How we operate

**As of** 2026-08-27.

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

**Brainstorm → Design spec → Implementation plan → Execute → Review → Verify → Promote.**

**Work is organized by task, not by date.**
Every task owns one directory, `docs/tasks/<task-slug>/`, named for what the task does with no date prefix.
The slug is the task's identity: it is what `SPEC.md`'s spec tree keys on and what `STATE.md` links to.
Dates live in git and in each document's own header line.

| Stage | Where it lives | Skill |
|---|---|---|
| Brainstorm | (conversation) | `superpowers:brainstorming` |
| Design spec | `docs/tasks/<task-slug>/design.md` | (write directly, or via brainstorming's output) |
| Implementation plan | `docs/tasks/<task-slug>/plan.md` | `superpowers:writing-plans` |
| Execute | code changes | `superpowers:executing-plans` or `superpowers:subagent-driven-development`, with `superpowers:test-driven-development` and `superpowers:systematic-debugging` as needed |
| Review | `docs/tasks/<task-slug>/review.md` — each finding and what was done about it | `superpowers:requesting-code-review` / `receiving-code-review` |
| Verify | `docs/tasks/<task-slug>/verification.md` — the commands run, their output, and pointers to everything the run produced: tables, figures, logs | `superpowers:verification-before-completion` — evidence before "done," always |
| Promote | `results/<task-slug>/<result>.csv` + a `<result>.provenance.json` record beside it — only the subset the project will cite, copied out of the task folder and made citable. The record carries the result's own checksum, and promotion refuses when the two copies differ | `scripts/promote_result.py` |
| Close out | branch merge/PR | `superpowers:finishing-a-development-branch` |

Three rules keep the task folder from becoming an orphan:

1. **Create the folder when you write the design spec**, and in that same change add the task's branch to `SPEC.md`'s spec tree, under its rung or under cross-cutting.
   A task folder outside the tree is invisible to the next session, which is the failure mode this structure exists to prevent.
2. **One task, one folder.**
   Extending existing work means editing that task's `design.md`, not starting `<something>-v2`.
   If the scope genuinely changes into different work, open a new task and supersede the old one explicitly (§5).
3. **`STATE.md` links back to the folder.**
   Every status entry names its task's `design.md`, the commits that implemented it, and the outputs it produced (§5) — that link is what makes a number traceable to the intent behind it.

Changes to design decisions and plans are updated in the `design.md` and `plan.md` files, respectively. 
Old approaches are moved to the bottom, where a dated desciption of the change is appended.  
Decisions that cut across tasks are linked to the `design.md` or other evidence in `docs/tasks/<other task slug>/` that drove the change relevant to the current task. 


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
- **Chained jobs**: pass `--dependency=afterok:<jobid>` etc. to `ralpine submit`; the wrapper reorders args correctly, but *verify* with `ralpine run scontrol show job <id> | grep Dependency` after submitting — a chain that silently drops its dependency runs immediately and out of order.
- **Verify inputs exist before submitting.**
  `ralpine ls`/`ralpine run find` the files a job needs first.
- **Pull the job's log into `results/<task-slug>/` when it finishes.**
  A log left on the cluster is unreadable to anyone reviewing the pull request, and is gone when the scratch space is cleared.
  It is the record of what the run actually did — the resolved arguments, the warnings, the wall time — so it belongs beside the result it produced.
- **Standing permissions** : `git commit`/`git push` to `origin`, `ralpine submit`/`cancel`/`update`, all without per-action confirmation.
  `submit` is the working default — an agent stages the run and submits it.
  `cancel` covers this project's own jobs that are stale (superseded by a newer submission, or running code that a later commit has already fixed) or running far past their expected wall time; never another user's job, and never a whole array range, which `ralpine cancel` refuses structurally by taking exactly one numeric id.
  Still confirm before force-push, `git reset --hard`, branch deletion, or anything touching `upstream` or `main`.
  A standing grant is scoped to what it says — if a new branch or a materially different kind of action comes up, ask.

## 3. Verification discipline

- **Test on synthetic data before spending cluster time.**
  Build a small synthetic fixture that exercises the real code path (not a reimplementation of its logic), run it locally, confirm the shape and sign of the result look right, *then* submit to Alpine.
- **Every reported statistic needs a known-answer test.**
  Plant a signal, require it recovered; plant nothing, require null.
  The test imports the real function.
- **Run the full local suite before every push**: `uv run pytest -q`.
- **Run the project-rule tests for what you touched.**
  `uv run pytest tests/test_project_rules.py` every time — they may catch a violation the task never intended — plus `-m` for the steps named in the task's `design.md` header, e.g. `uv run pytest -m "step_split or step_score"`.
  The step-to-rule-to-test mapping sits under each rule in `docs/SPEC.md`, and each step's marker is registered in `pyproject.toml` alongside the test that uses it, so the selection is the same whichever rung or dataset the task is about.
  An exemption in that run is a recorded gap, not a pass: if the task closed one, the test starts passing unexpectedly and the exemption comes out in the same change.
- **`verification-before-completion` applies to claims, not just code**: before saying something is fixed, promoted, or done, show the command and its output.
  "Should work now" is not done.

## 4. Git and collaboration

- **One branch at a time.**
  Finish the work, open the pull request, merge it, then branch again.
  Branching before the current piece lands puts two half-finished versions of the same documents on disk, and every question after that — which number is current, which file is real — has more than one answer.
- **Branch names say what the work is**, not who or when.
  A branch is named for a task and lives as long as that task.
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

## 6. Definition of done

- **A bug fix**: a real test (imports the actual function) passes locally; if it touches a promoted number, that result is re-run, re-promoted with a fresh provenance record, and the document that reported the old number gets a correction banner — never a silent edit.
- **A new capability**: went through brainstorm → spec → plan in one `docs/tasks/<task-slug>/` folder, indexed in `SPEC.md`, with `STATE.md`'s entry linking spec, code, and outputs.
- **A claim of "complete"**: backed by command output you actually ran, not by what the code should do.
- **Any task**: the README reflects it. A rung that landed, a rule that arrived, a status that moved — if a reader would learn it from the three documents, the README says it too.

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
