# fm-pdo-evaluator

Foundation-model evaluation harness for patient-derived tumor organoid (PDTO) drug-response prediction. Realizing the benefits of foundation models requires careful evaluations that map the boundaries of generalization.

The harness produces a registry-backed report comparing models against three out-of-distribution split strategies on PDTO datasets, with negative/positive controls, bootstrap confidence intervals, and a pretraining-leakage exposure profile per result.


## Quickstart

```bash
# Install uv (Python package manager) if you don't already have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Sync dependencies (creates .venv and uv.lock)
uv sync --extra dev

# Run the tests
uv run pytest
```

## Test Dataset

- **Soragni 2024** sarcoma PDTOs ([Synapse PDTOSarcoma](https://www.synapse.org/PDTOSarcoma))



## License

BSD-2-Clause Plus Patent License (see [LICENSE](LICENSE)).
