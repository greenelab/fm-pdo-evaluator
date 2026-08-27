# How we operate

This is the operating manual: the workflow, tooling, and collaboration mechanics for this project.
It answers "what do I do and in what order," not "what must every number respect" (that's `docs/PROJECT_SPEC.md`'s project rules) and not "what's true right now" (that's `docs/PROJECT_STATE.md`).
Read those two after this one, before starting work.

Some of what this manual names does not exist yet — the cluster wrapper, the promotion script, the tests enforcing the project rules.
Each arrives with the rung that first needs it.
The practice is written down first so those tools get built to it, rather than the practice being reverse-engineered from whatever the tools happened to do.

## 0. Start here, every session

1. Read `docs/PROJECT_SPEC.md`'s spec tree — is there already an OPEN branch covering what you're about to do?
   Extend that task's `docs/tasks/<task-slug>/design.md` or explicitly supersede it; don't open a second task on the same subject.
2. Read `docs/PROJECT_STATE.md` — is the thing you're about to fix already fixed, already known broken with a reason, or already being worked in a session you'd be racing?
3. Check `git log --oneline -20` on the branch you're on.
   Multiple sessions (this machine and Alpine, sometimes concurrently) touch this repo; the two documents above are the coordination mechanism, not tribal memory of who's doing what.


## 1. The work lifecycle

**Rung, step, task** are defined in [`docs/PROJECT_SPEC.md`](PROJECT_SPEC.md) and used precisely here: a rung is a scientific question, a step is a stage of the pipeline every rung passes through, a task is a unit of work with a spec and a definition of done.
Tasks open and close; steps never do; rungs are answered.


**Brainstorm → Design spec → Implementation plan → Execute → Review → Verify → Promote.**

**Work is organized by task, not by date.**
Every task owns one directory, `docs/tasks/<task-slug>/`, named for what the task does (`rung0-replicate-ceiling`, `embargo-gate-cell-values`) with no date prefix.
The slug is the task's identity: it is what `PROJECT_SPEC.md`'s spec tree keys on and what `PROJECT_STATE.md` links to.
Dates live in git and in each document's own header line.

| Stage | Where it lives | Skill |
|---|---|---|
| Brainstorm | (conversation) | `superpowers:brainstorming` |
| Design spec | `docs/tasks/<task-slug>/design.md` | (write directly, or via brainstorming's output) |
| Implementation plan | `docs/tasks/<task-slug>/plan.md` | `superpowers:writing-plans` |
| Task-local decisions | `docs/tasks/<task-slug>/decisions.md` | — |
| Execute | code changes | `superpowers:executing-plans` or `superpowers:subagent-driven-development`, with `superpowers:test-driven-development` and `superpowers:systematic-debugging` as needed |
| Review | — | `superpowers:requesting-code-review` / `receiving-code-review` |
| Verify | — | `superpowers:verification-before-completion` — evidence before "done," always |
| Promote | `docs/results/*.csv` + `.provenance.json` sidecar | `scripts/promote_result.py` |
| Cross-task reversal | `docs/decisions/YYYY-MM-DD-<slug>.md` | — |
| Close out | branch merge/PR | `superpowers:finishing-a-development-branch` |

Three rules keep the task folder from becoming another orphan:

1. **Create the folder when you write the design spec**, and in that same change add the task's branch to `PROJECT_SPEC.md`'s spec tree, under its rung or under cross-cutting.
   A task folder outside the tree is invisible to the next session, which is the failure mode this structure exists to prevent.
2. **One task, one folder.**
   Extending existing work means editing that task's `design.md`, not starting `<something>-v2`.
   If the scope genuinely changes into different work, open a new task and supersede the old one explicitly (§5).
3. **`PROJECT_STATE.md` links back to the folder.**
   Every status entry names its task's `design.md`, the commits that implemented it, and the outputs it produced (§5) — that link is what makes a number traceable to the intent behind it.

`docs/decisions/` stays dated, because a decision is an event, not a task: it records what was reversed and when.
Decisions made *inside* a task, that only bind that task, go in its `decisions.md`.

**Promotion is the line between "we ran something" and "it's evidence."**
A result without a `.provenance.json` sidecar is not citable in a write-up, per this project's own standing rule — see project rule 8 in `PROJECT_SPEC.md`.

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
- **Standing permissions** : `git commit`/`git push` to `origin`, `ralpine submit`/`cancel`/`update`, all without per-action confirmation.
  `submit` is the working default — an agent stages the run and submits it.
  `cancel` covers this project's own jobs that are stale (superseded by a newer submission, or running code that a later commit has already fixed) or running far past their expected wall time; never another user's job, and never a whole array range, which `ralpine cancel` refuses structurally by taking exactly one numeric id.
  Still confirm before force-push, `git reset --hard`, branch deletion, or anything touching `upstream` or `main`.
  A standing grant is scoped to what it says — if a new branch or a materially different kind of action comes up, ask.

## 3. Verification discipline

- **Test on synthetic data before spending cluster time.**
  Build a small synthetic fixture that exercises the real code path (not a reimplementation of its logic — see project rule 4), run it locally, confirm the shape and sign of the result look right, *then* submit to Alpine.
- **Every reported statistic needs a known-answer test.**
  Plant a signal, require it recovered; plant nothing, require null.
  The test imports the real function.
- **Run the full local suite before every push**: `uv run pytest -q`.
- **Run the project-rule tests for what you touched.**
  `uv run pytest tests/test_project_rules.py` every time — six of those are repository- or artifact-wide scans, so they catch a violation the task never intended — plus `-m` for the steps named in the task's `design.md` header, e.g. `uv run pytest -m "step_split or step_score"`.
  The step-to-rule-to-test mapping sits under each rule in `docs/PROJECT_SPEC.md`, and the steps are markers registered in `pyproject.toml`, so the selection is the same whichever rung or dataset the task is about.
  An exemption in that run is a recorded gap, not a pass: if the task closed one, the test starts passing unexpectedly and the exemption comes out in the same change.
- **`verification-before-completion` applies to claims, not just code**: before saying something is fixed, promoted, or done, show the command and its output.
  "Should work now" is not done.

## 4. Git and collaboration

- **Branches**: `worktree-<topic>`, checked out under `.claude/worktrees/<topic>/` alongside the main checkout (`git worktree list` shows all active ones — check for collisions before naming a new one).
- **Commit messages state the root cause, not just the change**, written so a future session (or another agent) can tell what happened and why without re-deriving it from the diff — this is the established convention in this repo's own history, not a preference introduced this session. Do not include specific result numbers in the commit message. Reference the outputs and summary documents. 
- **Push to `origin`** (the `lagillenwater` fork).
  **Never `upstream`** (`greenelab`) without being explicitly told.
  Greene Lab standard beyond this fork: fork-per-contributor, PRs only, minimum one lab-member approval, one functional area per PR.
- **Multiple sessions run concurrently** (this branch carries provenance sidecars pointing at *another machine's* `/private/tmp/` scratch path, written the same day by a different session).
  Don't assume you're the only session touching this repo — `docs/PROJECT_STATE.md` and recent git log are how you find out, not memory of who was told to do what.

## 5. Documentation discipline

Full rules in `docs/PROJECT_SPEC.md`'s "Process for new specs" section; the short version:

- A new task = a `docs/tasks/<task-slug>/design.md` **and** a spec-tree branch in `PROJECT_SPEC.md` pointing at it, in the same change.
- Superseding a doc = a one-line dated banner on the *old* doc pointing to the new one, **and** an index-row update in `PROJECT_SPEC.md`, in the same change.
  Not one or the other.
- Reversing an architecture, method, or data source = a decision entry in `docs/decisions/`. No exceptions.
- A number change = `docs/PROJECT_STATE.md` updated in the same change that produced it, and that entry carries its three links: **spec** (the task's `design.md`), **code** (the commits or files that changed), **outputs** (the promoted CSV + sidecar, job id, figure).
  A status claim with no link to the spec that motivated it, the code that produced it, and the artifact it came from is a claim a future session has to re-derive from git log — which is exactly how the CoderData reversal stayed invisible for two months.

## 6. Definition of done

- **A bug fix**: a real test (imports the actual function) passes locally; if it touches a promoted number, that result is re-run, re-promoted with a fresh sidecar, and the document that reported the old number gets a correction banner — never a silent edit.
- **A new capability**: went through brainstorm → spec → plan in one `docs/tasks/<task-slug>/` folder, indexed in `PROJECT_SPEC.md`, with `PROJECT_STATE.md`'s entry linking spec, code, and outputs.
- **A claim of "complete"**: backed by command output you actually ran, not by what the code should do.

## 7. Failure modes this process exists to prevent

Written down because each has already cost this project real time, and because the lesson outlives the specific confusion.

- **A newer, more authoritative-looking document can itself be the source of drift.**
  The most recent design doc is the one everyone cites, which is exactly why an error in it propagates fastest.
  Recency is not accuracy: a new spec must be checked against the code and the specs it supersedes, not just written down.
- **Renaming a concept without a forward pointer leaves every prior document silently ambiguous.**
  A code-level rename with a working alias reads as complete while every document using the old name keeps meaning something slightly different.
  A rename is a two-line fix in the documents that used the old name; do it at rename time, not retroactively.
- **A parallel implementation that bypasses an abstraction silently drops that abstraction's guarantees.**
  When a driver is written to call a library directly because the registry path is inconvenient, every guarantee the registry provided — leakage filtering, shared partitions — quietly stops applying, and no output says so.
  If a guarantee matters, it belongs where it cannot be routed around.
- **A reversal with no decision entry is invisible to everyone but git log.**
  Reverting an architecture, a method or a data source and writing nothing down leaves the next reader to rediscover it from commit archaeology, or to reinstate what was deliberately removed.
  Project rule 10 exists for this; a terse entry written late still beats none.
