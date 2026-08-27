"""Provenance models: ``LeakageProfile``, ``EnvironmentSnapshot`` and ``PromotedResult``."""

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

    model_config = ConfigDict(extra="forbid", frozen=True, protected_namespaces=())

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


class PromotedResult(BaseModel):
    """The record written beside a promoted result, saying how it was produced.

    A prediction carries its provenance inside ``PredictionRecord``; a promoted *artifact* —
    the table a write-up cites — needs the same account at file granularity, which is what this
    model is. It embeds ``EnvironmentSnapshot`` rather than restating it, so a commit is a
    ``code_commit`` and a seed is a ``seed`` wherever provenance is written in this project.

    Three fields exist because they cannot be reconstructed afterwards. ``clean_tree`` records
    whether the working tree carried uncommitted changes at promotion, which no later inspection
    can recover. ``result_sha256`` pins the artifact so an edit made after promotion is
    detectable rather than silent. ``job_id`` is optional because a result computed locally has
    no scheduler to name.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    result: str = Field(min_length=1)
    result_sha256: str = Field(pattern=r"^[0-9a-f]{64}$")
    task: str = Field(min_length=1)
    script: str = Field(min_length=1)
    args: dict[str, str] = Field(default_factory=dict)
    inputs: dict[str, str] = Field(min_length=1)
    log: str | None = None
    log_sha256: str | None = Field(default=None, pattern=r"^[0-9a-f]{64}$")
    job_id: str | None = None
    clean_tree: bool
    environment: EnvironmentSnapshot
    promoted_at: datetime
