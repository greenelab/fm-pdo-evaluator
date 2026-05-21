# Environment contract

This document describes the harness's environment-control contract: which
parts of the runtime must be reproducible, how each is pinned, and how the
provenance is captured in every prediction record. Without these guarantees,
the determinism check is half-blind to node-to-node
drift and the leakage scan can't be tied to a specific code/data state.

## 1. Containers

Foundation-model inference runs inside Apptainer images pinned by digest in
`containers/digests.json`:

| Image | Definition | Built on | Wraps |
|---|---|---|---|
| `fmharness` | `containers/fmharness.def` | Day 1+2 (skeleton); rebuilt Day 8 | core Python deps |
| `tahoe` | `containers/tahoe.def` | Day 8 | Tahoe-x1 + torch + CUDA |
| `state` | `containers/state.def` (built only if needed) | Day 11 | STATE + torch + CUDA |

STATE reuses the Tahoe container unless torch/CUDA conflicts force a split;
the decision (and the reason) is recorded on Day 11 in this document.

Every `PredictionRecord` carries `EnvironmentSnapshot.container_digest`. A
prediction made outside a pinned container fails determinism check #6.

## 2. Deterministic GPU execution

Every CLI entrypoint calls `fmharness.utils.determinism.fix_seeds(seed)`
before importing torch CUDA functionality. That call:

- seeds `random`, `numpy`, `torch`, `torch.cuda`
- sets `PYTHONHASHSEED`
- sets `CUBLAS_WORKSPACE_CONFIG=:4096:8` (required for cuBLAS determinism)
- calls `torch.use_deterministic_algorithms(True)`

The `CUBLAS_WORKSPACE_CONFIG` env var must be set before CUDA initializes,
so `fix_seeds` is unsafe to call after model load. Callers that load torch
modules at import time must invoke `fix_seeds` before any such import.

`EnvironmentSnapshot.cuda_deterministic` records whether this contract was
active when the prediction ran. Set it to `True` only if `fix_seeds` was
called and `torch.are_deterministic_algorithms_enabled()` returned `True`.

## 3. Secrets

Two secrets the harness needs are kept out of the repo:

- `HUGGINGFACE_TOKEN` — Tahoe-x1 weights download
- `SYNAPSE_PAT` — Soragni 2024 dataset access

Layout:

- Local dev: copy `.env.example` to `.env`, fill in values. `.env` is
  gitignored. Loaded by `pydantic-settings` from the repo root.
- Alpine: `~/.fmharness/secrets` with `chmod 600`. Slurm sbatch headers
  source it explicitly: `source ~/.fmharness/secrets`.

`.env` and any path containing `secrets` are caught by the pre-commit
`detect-private-key` hook and a project-specific token-pattern check
(added Day 13 when secrets handling is fully exercised).

## 4. Static asset versioning

Static assets (the Tahoe-100M drug list, drug crosswalk tables, gene-panel
reconciliation tables, reference FASTA + GTF used by the RNA quantification
pipeline) live under `data/static/` with a `manifest.json` recording sha256
per file. Loaders verify on read and refuse to proceed on mismatch.

The manifest's hash of the reference genome + GENCODE annotation propagates
into `EnvironmentSnapshot.data_commit` via the tranche content hash, so a
mismatch is detectable from any prediction record.

## 5. Pinned Python environment

`uv.lock` (committed) pins every direct and transitive Python dependency.
CI installs from the lock; Alpine inference jobs install from the lock
inside the Apptainer build (`uv sync --frozen --no-dev`).

When upgrading a dependency:

1. Edit `pyproject.toml`.
2. Run `uv lock` to refresh `uv.lock`.
3. Rebuild any Apptainer image whose `%files` includes the lock.
4. Commit both files in the same PR.

## 6. `EnvironmentSnapshot` field reference

| Field | Source | Required for det. check |
|---|---|---|
| `code_commit` | `git rev-parse HEAD` at run time | yes |
| `container_digest` | `apptainer inspect --digest <image>` | yes (GPU runs) |
| `python_version` | `sys.version` | yes |
| `torch_version` | `torch.__version__` if torch available | yes (GPU runs) |
| `cuda_version` | `torch.version.cuda` if CUDA available | yes (GPU runs) |
| `model_weights_hash` | sha256 of the loaded checkpoint bytes | yes |
| `data_commit` | tranche `content_hash` for the inputs | yes |
| `seed` | seed passed into the CLI | yes |
| `cuda_deterministic` | `True` iff `fix_seeds` ran + det algos enabled | yes |
