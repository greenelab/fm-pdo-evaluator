# fm-pm-evaluator

Foundation models trained on cell-line transcriptomic data claim generalization to
out-of-distribution tasks. This repository builds the infrastructure to find the boundary of
those claims: whether a model's prediction survives the move to an unseen cell line, across
platforms, into a new modality, and onto patient-derived tumour organoids.

The evaluation is a **ladder**. Each rung adds exactly one boundary to the rung below it, so a
score that falls off says *which* boundary it fell on rather than that something, somewhere, went
wrong.

| Rung | Question it settles |
|---|---|
| 0 | Is the target reproducible? |
| 1 | Can a model predict an unseen cell line, in-distribution? |
| 2 | Can a bulk sample be read by a single-cell model? |
| 3 | What does crossing a measurement platform cost? |
| 4 | Does the representation predict drug response, not just expression? |
| 5 | Does it transfer to patient-derived organoids? |
| 6 | Does it hold when the prediction comes first? |

Rungs are built one at a time — each task on its own branch, named for the task — and a rung's
spec arrives with its implementation, so the design of a measurement and the code producing it
are reviewed together. Every measurement ships with a positive and a negative control, and every
promoted comparison reports the smallest effect it was powered to detect. Every task passes a
drift audit of its design against what landed before promoting its result — the procedure is
[docs/audit.md](docs/audit.md), one standard for the whole repository — and carries a
plain-language summary notebook that walks a reviewer through the hypothesis, each measurement
step with its figures and controls, and the conclusions, before opening its pull request as a
draft for a final review round. A task has two points where a person signs off — the design spec
before any code is written, and the summary once the result is in — and the stages between them
run on their own checks rather than on approvals. Nothing a run produces enters the git record
until after that second sign-off: results stay on the cluster and in the working tree through
verification, audit and summary, and are committed once, in the form a person has read. Its pull request sends the reader to the summary
first, then the design, then the audit, and separates the files needing human eyes from the
generated tables committed so the numbers can be re-derived without cluster access. Each task also ships an executable verification entry point —
a script plus a notebook that recompute every promoted claim from the committed artifacts, so a
reviewer re-derives the numbers rather than trusting the write-up. The dated history of project
and task decisions lives in [docs/decisions.md](docs/decisions.md) and each task's own
`decisions.md`. All rungs are scored
inside one declared evaluation frame — the same datasets and the same statistic — so their
numbers divide like by like; a new dataset opens a new frame rather than silently shifting the
old one, and a rung scoring a subset of genes or drugs reads against a restriction of the rung
below rather than against a different frame.

## The documents

Three documents carry the project, and they answer different questions. Read them in this order.

- **[docs/SPEC.md](docs/SPEC.md)** — what the project asks and what every result must obey: the
  question, the ladder rung by rung, and the project rules.
- **[docs/PROCESS.md](docs/PROCESS.md)** — how work moves from a question to a promoted result:
  the lifecycle, the task folders, the compute boundary, what "done" requires.
- **[docs/STATE.md](docs/STATE.md)** — where each rung stands right now, and what stands in the
  way.

This README summarises those three. When a task changes what the project asks, how work is done,
or where a rung stands, it updates the README in the same change — see `PROCESS.md` §5.

## Status

Rung 0 — the reliability ceiling every higher rung is read against — is in progress on the
`rung0-assay-reliability` branch ([design](docs/tasks/rung0-assay-reliability/design.md)). It
measures the ceiling at the assay's full extent, every gene and every splittable drug; a higher
rung reads against a restriction of that ceiling to the genes and drugs it scores, declared
before it scores. [`docs/STATE.md`](docs/STATE.md) is authoritative.

## Quickstart

```bash
# Install uv (Python package manager) if you don't already have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies (creates .venv)
uv sync --extra dev

# Run the tests, including the project-rule checks
uv run pytest
```

Dependencies are declared as the code that imports them arrives, so the manifest describes what
the project runs on rather than what it might run on later.

## Datasets

Every dataset the project reads is described in the registry, [docs/DATA.md](docs/DATA.md) —
what it is, the paper behind it, the script and date of its download, and the processing applied
after. Currently registered:

- **Tahoe-100M** (Vevo Therapeutics; Zhang et al. 2025, CC0-1.0) — ~100M single-cell profiles,
  50 cancer cell lines, ~1,100 drug–dose perturbations on 14 plates with plate-matched DMSO
  controls. The delta rungs (0–3) read its pseudobulk differential-expression table, restricted
  to the 32 drugs shared with GDSC2.
- Datasets for the higher rungs — the viability screens and the embargoed organoid cohort the
  ladder names in [docs/SPEC.md](docs/SPEC.md) — enter the registry with their rungs.

## License

BSD-2-Clause Plus Patent License (see [LICENSE](LICENSE)).
