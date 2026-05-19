"""Provenance models: ``LeakageProfile`` and ``EnvironmentSnapshot``."""

from __future__ import annotations

from datetime import datetime

from pydantic import BaseModel, ConfigDict, Field


class LeakageProfile(BaseModel):
    """Pretraining-exposure summary attached to a model run on a tranche.

    ``drug_overlap_tahoe_100m`` maps drug_id -> bool indicating whether the
    drug appears in the declared Tahoe-100M pretraining corpus.
    ``drug_overlap_fraction`` is the fraction of drugs in the tranche that
    appeared in the corpus. ``declared_corpus_overlap`` carries per-corpus
    overlap fractions for any other corpora declared by the model wrapper.
    ``subtype_prevalence`` records the prevalence of each tranche subtype in
    the pretraining corpus (when knowable).
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    tranche_id: str = Field(min_length=1)
    model_version: str = Field(min_length=1)
    drug_overlap_tahoe_100m: dict[str, bool] = Field(default_factory=dict)
    drug_overlap_fraction: float = Field(ge=0.0, le=1.0)
    declared_corpus_overlap: dict[str, float] | None = None
    subtype_prevalence: dict[str, float] = Field(default_factory=dict)
    generated_at: datetime


class EnvironmentSnapshot(BaseModel):
    """Captured environment for a prediction run.

    Embedded in every ``PredictionRecord``. The combination of
    ``code_commit``, ``container_digest``, ``model_weights_hash``,
    ``data_commit``, and ``seed`` uniquely identifies the inputs that
    produced a prediction; ``cuda_deterministic`` records whether the
    determinism contract was active.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    code_commit: str = Field(min_length=7)
    python_version: str = Field(min_length=1)
    seed: int
    cuda_deterministic: bool
    data_commit: str = Field(min_length=1)
    container_digest: str | None = None
    torch_version: str | None = None
    cuda_version: str | None = None
    model_weights_hash: str | None = None
