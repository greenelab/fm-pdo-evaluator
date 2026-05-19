# Adapter contract

How a model joins the harness. Every model — linear baseline, Tahoe-x1, STATE, and any future addition — implements the same `ModelAdapter` Protocol so the rest of the pipeline (splits, probe, metrics, registry, leakage scan) is model-agnostic.

The contract has three required surfaces (`embed`, `metadata`, `version`) and one optional surface (`predict_native`). Day 7 implements the Protocol and the `linear_baseline` and `MockAdapter` reference implementations; later days add `tahoe_x1` (Day 8) and `state` (Day 11) against the same contract.

## 1. The interface

```python
from typing import Protocol, runtime_checkable
import numpy as np
from anndata import AnnData


@runtime_checkable
class ModelAdapter(Protocol):
    """Wrap a foundation model (or baseline) so the harness can use it."""

    def version(self) -> str:
        """Stable identifier, e.g. `tahoe_x1@v1.0.0`. Embedded in PredictionRecord."""

    def metadata(self) -> "ModelMetadata":
        """Pretraining provenance for the leakage scan."""

    def embed(self, adata: AnnData) -> np.ndarray:
        """Encode samples (rows of `adata`) into a dense embedding matrix.

        Returns an array of shape ``(adata.n_obs, embedding_dim)``. Rows are
        aligned with ``adata.obs_names``. The harness owns input ordering;
        adapters must not reorder rows.
        """

    def predict_native(
        self, adata: AnnData, drug_ids: list[str]
    ) -> np.ndarray | None:
        """Optional native drug-aware head. Return ``None`` if the adapter
        does not expose one. Shape: ``(adata.n_obs, len(drug_ids))``.
        """
```

`ModelMetadata` (a pydantic model declared in `schema/`):

| Field | Required | Purpose |
|---|---|---|
| `pretraining_corpus` | yes | Free-form name (`"tahoe_100m"`, `"none"` for baselines) |
| `pretraining_cutoff_date` | yes | ISO date the corpus was frozen; used to flag dataset-leakage risk |
| `task_signal_in_pretrain` | yes | One of `"none"`, `"adjacent"`, `"direct"` — declares whether the corpus contained drug-response labels |
| `model_weights_hash` | yes for FM models | sha256 of the loaded checkpoint, captured into `EnvironmentSnapshot` |
| `container_digest` | yes for FM models | Apptainer image digest the adapter expects to run inside |

## 2. The probe-based prediction pipeline (default path)

All four matrix rows (linear baseline, Tahoe-x1, STATE, plus the metadata-only control) share an identical probe so the comparison isolates "what the encoder captures":

```
  sample (RNA-seq)  --[ encoder ]-->  embedding  --[ concat drug feat ]-->  [ probe ]  --> P(responder)
                                                                              ^
                                                                              |
                                                                       trained per split-fold
```

- **Encoder** is model-specific. For the linear baseline the encoder is `StandardScaler` (a passthrough; "embedding" == scaled expression). For Tahoe-x1 / STATE it is the pretrained transformer encoder.
- **Drug feature** is a one-hot over the drug crosswalk's canonical IDs at MVP; richer drug descriptors (Morgan fingerprint, ATC class) are a deferred extension.
- **Probe** is a fixed architecture across all models: `StandardScaler → ElasticNetCV` (continuous response) or `LogisticRegressionCV` (binary responder). Declared once in `src/fmharness/probe/linear.py`. The harness — not the adapter — owns the probe.

The adapter's only job is to produce a faithful embedding. Probe training and inference are downstream.

## 3. The native-head path (optional)

Foundation models with a drug-aware head can return a prediction directly:

- Tahoe-x1: trained on perturbation-response prediction; may expose a `predict(baseline_state, drug) → post_state` or scalar head.
- STATE: the ST (state transition) component is exactly this.

When `predict_native` returns a value, the harness records it as a separate row in the registry tagged `prediction_mode="native"`. The probe-based row (`prediction_mode="probe"`) is always produced for fair comparison; the native row is supplementary.

Adapters that do not expose a native head return `None` and only the probe path runs.

## 4. Required behaviors

Every adapter must:

1. **Be deterministic.** Calling `embed` twice on the same input must produce identical output bits-for-bits when `fmharness.utils.determinism.fix_seeds(seed)` has been called. The determinism check on Day 14 will fail loud otherwise.
2. **Not reorder rows.** `adata.obs_names[i]` must correspond to `embedding[i]`.
3. **Declare its container.** GPU adapters run inside the Apptainer image whose digest is returned by `metadata().container_digest`. The harness refuses to record a `PredictionRecord` whose `EnvironmentSnapshot.container_digest` does not match.
4. **Cache by content.** Embedding caches under `data/tranches/{tranche_id}/embeddings/{model_version}/` are keyed by sha256 of the input AnnData bytes plus `model_version`. Adapters should call into `fmharness.utils.cache` rather than rolling their own.
5. **Gracefully refuse mismatched inputs.** If the input gene panel does not match the adapter's expected reference, raise `GenePanelMismatch` rather than silently aligning — gene-panel reconciliation is the loader's responsibility (Day 4), not the adapter's.

## 5. Adding a new model

1. Build / pick the Apptainer image (`containers/<name>.def`). Pin the digest in `containers/digests.json`.
2. Create `src/fmharness/models/wrappers/<name>.py` implementing `ModelAdapter`.
3. Register the adapter in `src/fmharness/models/registry.py` so the CLI's `--model <name>` flag resolves.
4. Add `configs/<name>_{soragni,yang}_{id,lpo,lso}.yaml` (6 files) following the existing pattern.
5. Add a smoke test under `tests/wrappers/test_<name>.py` using the `MockAdapter` round-trip pattern.

## 6. The metadata-only control adapter

The metadata-only baseline is implemented as a `ModelAdapter` too — it lives in `src/fmharness/controls/negative.py` and the same probe is applied. Its `embed()` returns the one-hot concatenation of tissue, subtype, and drug ID; its `metadata()` declares `pretraining_corpus="none"` and `task_signal_in_pretrain="none"`. The leakage scan reports zero overlap for metadata-only rows by construction.

## 7. Reference: the linear baseline

```python
class LinearBaselineAdapter:
    def version(self) -> str:
        return "linear_baseline@v1.0"

    def metadata(self) -> ModelMetadata:
        return ModelMetadata(
            pretraining_corpus="none",
            pretraining_cutoff_date=date(1970, 1, 1),
            task_signal_in_pretrain="none",
        )

    def embed(self, adata: AnnData) -> np.ndarray:
        # The "encoder" is identity; downstream StandardScaler in the probe
        # handles per-feature standardization.
        return np.asarray(adata.X, dtype=np.float32)

    def predict_native(self, adata, drug_ids):
        return None
```

Anything more elaborate than this for the baseline (e.g., PCA, per-gene z-score, gene selection) belongs in the probe pipeline, not the adapter.
