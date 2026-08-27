# Project state

**Reads against** [`docs/PROJECT_SPEC.md`](PROJECT_SPEC.md).
The spec says what each rung must establish and what a passing result means; this document says where each one stands.
It carries no history: a rung's result belongs here, how it came to be belongs in git and in that rung's spec.

**As of** 2026-08-27.

A number not carried here with its provenance record is not evidence.
Promotion means a result in `docs/results/` with a `.provenance.json` beside it recording the commit, job and inputs that produced it — project rule 8.
`docs/results/` does not exist yet, so nothing in this repository is evidence of anything.

---

## Ladder status

| Rung | What the spec requires | Status |
|---|---|---|
| 0 — replicate ceiling | A reproducibility ceiling clearing its null, on the declared panel | Next |
| 1 — held-out line | Prediction beating a floor and recovering a planted signal, as a fraction of rung 0 | Not started |
| 2 — cross-platform | Retention separable from a scrambled-line control | Not started |
| 3 — GDSC2 viability | Interaction above zero after correction, against the screen-agreement ceiling | Not started |
| 4 — organoid viability | The transfer number, on a frozen embargoed holdout | Not started |

A rung closes when its result is promoted with provenance, this table records it, and the project-rule tests for the steps it touches pass.
Rungs are built in order, because each one is read against the one below it: rung 0's ceiling is the denominator for rung 1, and a rung whose denominator is unmeasured reports a ratio to an imaginary 1.0.

## What this repository holds today

| Present | Consequence |
|---|---|
| Schema, determinism and adapter scaffolding, with tests | The apparatus a rung is added to exists; nothing here yet produces a measurement |
| No `docs/results/` | Rung 0 is the first work to promote a number, and the first to be held to the rules in the spec |
| `docs/adapter_contract.md` and `docs/environment.md`, predating this spec | Neither has been reconciled against it. The rung that first depends on either brings it into line rather than a sweep that touches everything at once |

## Where things live

- **Results** `docs/results/*.csv` with a matching `.provenance.json`. No sidecar, no evidence.
- **Figures** `docs/figures/*.png`, generated rather than hand-made.
- **Rung and task specs** `docs/tasks/<slug>/design.md`, one folder per task, arriving with the work it specifies.
- **Decisions** `docs/decisions/YYYY-MM-DD-<slug>.md` for anything reversed or chosen against an alternative.
- **Rules and their tests** `docs/PROJECT_SPEC.md` and `tests/test_project_rules.py`.
