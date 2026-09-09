# The drift audit

**As of** 2026-09-01.

The repository-wide standard for the audit stage ([`docs/PROCESS.md`](PROCESS.md) §1). It says
what an audit is, how to run one, and what its record must contain. Each task's own audit record
lives at `docs/tasks/<task-slug>/audit.md`; this document is the procedure those records follow.
Decision lineage in [`docs/decisions.md`](decisions.md).

## What it is

The audit is a **diff check between a task's documents and the tree the task actually landed**.
It reads the design as a list of numbered, checkable claims, checks each one against the code and
artifacts on the branch, and classifies every departure. It is a documentary check, not a code
review: the review stage asks whether the code is good, the verify stage asks whether the numbers
are right, and the audit asks whether the documents describe what is there.

It exists because a design is written before the work and edited during it, and the edits go in
where the writer is looking. What survives is a document that reads as authoritative and is
quietly wrong in places nobody re-read. That is not hypothetical here: the first audit of the
first rung found the ceiling had been measured over 31 of 32 declared compounds, hidden because
the wrong number and the right number happened to coincide.

**The audit passes before the summary and before promotion.** A number promoted ahead of its audit
is citable before anyone has checked that the documents describe how it was produced. The first
rung audited after promoting, found exactly that class of error, and moved the stage.

The audit therefore reads the run's artifacts **in the working tree**, before they are committed
(PROCESS §1, "What reaches GitHub, and when"). It records the checksum of every artifact it read,
and the provenance records written at promotion carry the same checksums — so an artifact that
changed between the audit and the commit is caught rather than assumed identical. An audit that
cites an artifact without its checksum has not checked that artifact.

## The two directions

A diff check that runs one way only is half an audit, and both halves have caught real defects.

**Forward — every claim reaches the tree.** Enumerate the design's claims and verify each against
the branch with named evidence. A claim with no counterpart in the tree is drift, whether the tree
is missing the work or the document is describing work that moved.

**Reverse — every change reaches a claim.** Read `git diff <merge-base>...HEAD --stat`, and for
each changed file ask which claim accounts for it. A file changed under no claim is either
undocumented work or scope that arrived without a decision entry; both are findings. This is the
direction that catches what was added rather than what was dropped, and it is the one an auditor
reading only the design will never run.

## Enumerating claims

Number the claims and keep the numbering stable, prefixed by source (`D` for `design.md`, `S` for
`docs/SPEC.md`, `P` for `plan.md`), so the fix wave and the re-audit can refer to `D14` and mean
one thing across three documents.

A claim is one thing a reader could find false on its own. Split a sentence carrying two
assertions into two claims; keep a genuinely single assertion whole even if it runs long. Do not
merge adjacent claims to shorten the table — merging is how the first rung's audit ended up with
two different clause counts for one pass over the tree, and had to explain the discrepancy in
prose rather than fixing it.

Every claim in the design's scope gets a verdict, including the ones that are obviously fine.
An audit that lists only findings cannot be checked for coverage.

## Verdicts

Exactly three, and each carries its evidence — a command and its output, a recomputed hash, a
recount, a file path with a line number. A verdict with no evidence is an opinion.

| Verdict | Means |
|---|---|
| **ALIGNED** | The tree does what the claim says. Evidence names how that was checked, not that it was. |
| **DEVIATION-RECORDED** | The tree departs from the claim, and a dated entry in the task's `decisions.md` already says so and why. A recorded deviation is a working process, not a defect. |
| **DRIFT** | The tree departs from the claim and nothing records it. Every drift item is fixed or recorded before the audit passes. |

Drift is not a severity. A wrong count and a stale file path are both drift; the fix wave decides
what each one costs.

## Who runs it

A **fresh reader** — an auditor who did not write the code or the design, and who works from the
documents and the tree rather than from the conversation that produced them. Familiarity is the
thing being controlled for: the writer reads the sentence they meant, and an audit is only useful
when someone reads the sentence that is there.

**Mechanical claims are checked by script, never by a reader.** Counts, tallies, checksums,
whether a cross-reference resolves, whether two copies of a table are byte-identical: these go to
a script, and the audit records the command and its output. Readers are for judgment — whether a
paragraph still describes the method — and a reader recounting rows is both slower and less
reliable than `wc -l`. One of the first rung's audit passes was spent reconciling the audit's own
arithmetic, which is the failure this rule prevents.

## The record

`docs/tasks/<task-slug>/audit.md`, in this order:

1. **Header** — the date, the commit audited, and who audited it.
2. **Method** — what was read, what was recomputed, and the command output that supports the pass.
3. **Counts** — total claims and the tally by verdict. One canonical count; if an earlier pass
   counted differently, say which is canonical and why, rather than carrying two.
4. **Clause verdicts** — one table per source document, every claim with its verdict and evidence.
5. **Reverse-direction findings** — changed files no claim accounted for.
6. **Fix wave** — each drift item, its disposition (fixed, recorded as a decision, or ruled not a
   defect with the reasoning), and the commits that carried it.
7. **Re-audit** — a fresh reader confirming every drift item is fixed or recorded, with its own
   verdict. **The audit is not passed until the re-audit says so.**

## The cap

The re-audit re-checks **only the items verdicted drift**. It does not re-enumerate the design, and
there is no confirmation pass of a confirmation pass. Without that cap, each fix wave introduces
small documentary defects the next pass finds, and prose maintenance generates findings without
bound — which is what happened before the rule existed.

Work added after an audit passes gets an **audit delta**: the same procedure over the new surface
only, with the same fresh-reader confirmation, appended to the same `audit.md` under its own dated
heading. The lifecycle converges per task; it is not re-entered per run.
