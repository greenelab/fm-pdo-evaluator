# fm-pdo-evaluator

Foundation-model evaluation harness for patient-derived tumor organoid (PDTO) drug-response prediction. Realizing the benefits of foundation models requires careful evaluations that map the boundaries of generalization.

The harness produces a registry-backed report comparing three models against three out-of-distribution split strategies on two PDTO datasets, with negative/positive controls, bootstrap confidence intervals, and a pretraining-leakage exposure profile per result.

See [docs/fm-pdo-evaluator-plan.md](docs/fm-pdo-evaluator-plan.md) for the 3-week plan.

## Quickstart

```bash
# Install uv (Python package manager) if you don't already have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies (creates .venv and uv.lock)
uv sync --extra dev

# Run the tests
uv run pytest
```

## Datasets

- **Soragni 2024** sarcoma PDTOs ([Synapse PDTOSarcoma](https://www.synapse.org/PDTOSarcoma))
- **Yang 2024** primary liver cancer PDOs ([Cancer Cell](https://www.cell.com/cancer-cell/fulltext/S1535-6108(24)00089-8))


## Affiliation

Greene Laboratory, University of Colorado Anschutz Medical Campus.

## License

BSD-2-Clause Plus Patent License (see [LICENSE](LICENSE)).
